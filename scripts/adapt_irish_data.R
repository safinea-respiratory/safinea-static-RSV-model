# ==============================================================
# adapt_irish_data.R
#
# Converts the Irish 2025/26 surveillance extracts into the input
# format this model reads, so the Ireland round can run on the same
# code as the multi-country RespiCompass round rather than on a fork.
#
# INPUT (confidential - never committed; see .git/info/exclude)
#   data/source_ie/RSV_weekly_counts.csv     weekly aggregate admissions
#   data/source_ie/RSV_monthly_prop_age.csv  age split per 4-week period
#   data/source_ie/country_monthly_births.csv  Eurostat DEMO_FMONTH
#
# OUTPUT
#   data/epidemiological/hospitaladmissions_IE.csv
#   data/epidemiological/hospitalburden_agegroups_IE.csv   (confidential)
#   data/population/births_by_month_IE.csv                 (public, Eurostat)
#   data/population/population_estimates_IE.csv            (public, derived)
#
# Two things are worth knowing about the conversion:
#
#   1. The age split is 4-WEEKLY here, where the RespiCompass data is
#      one window per season. That is strictly more information: the
#      sampler constrains the age mix every four weeks instead of once
#      a season, so the draws are tighter. The sampler already iterates
#      over whatever burden windows it is given, so nothing changes in
#      the model.
#
#   2. Age bands are the six the Irish data actually carries. The model
#      takes its bands from the data, so no code cares - but the adult
#      programme can then only target "65+y" as one block.
#
# Run from the project root:  Rscript scripts/adapt_irish_data.R
# ==============================================================

suppressMessages({
  library(dplyr); library(tidyr); library(purrr)
  library(tibble); library(lubridate); library(readr)
})

SRC <- "data/source_ie"
stopifnot(dir.exists(SRC))

COUNTRY <- "Ireland"
ISO2    <- "IE"

# Labels the Irish extract uses -> the labels this model will carry
AGE_MAP <- c(
  "< 3 months"  = "0-2mo",
  "3-5 months"  = "3-5mo",
  "6-11 months" = "6-11mo",
  "1-4 years"   = "1-4y",
  "5-64 years"  = "5-64y",
  "65+ years"   = "65+y"
)

say <- function(...) cat(" ", sprintf(...), "\n", sep = "")


# ---- 1. weekly admissions ------------------------------------------------
# date_wk_floor is the MONDAY of the ISO week; the model keys on the
# week-ENDING date, so add 6 days.
raw_wk <- read_csv(file.path(SRC, "RSV_weekly_counts.csv"),
                   show_col_types = FALSE)

admissions <- raw_wk %>%
  transmute(
    country                     = COUNTRY,
    age_group                   = "total",
    target_end_date             = as.Date(date_wk_floor) + 6,
    year_week                   = sprintf("%d-W%02d",
                                          isoyear(target_end_date),
                                          isoweek(target_end_date)),
    week                        = isoweek(target_end_date),
    year                        = isoyear(target_end_date),
    weekly_rsv_hospitalisations = as.integer(case_counts)
  ) %>%
  arrange(target_end_date)

stopifnot(!anyDuplicated(admissions$target_end_date))

# The source OMITS weeks with no admissions rather than recording a zero
# - its minimum count is 1, and every gap falls in the summer off-season.
# The model needs a contiguous weekly grid (the submission's `horizon` is
# a week index), so those weeks are filled back in as zeros. This changes
# no total; it only restores rows the extract dropped.
#
# Only SHORT gaps are filled. The 2020-2022 hole is missing surveillance,
# not zero admissions, and inventing zeros across it would be a fiction.
MAX_FILL_WEEKS <- 6

fill_weekly_gaps <- function(df) {
  d    <- sort(df$target_end_date)
  gaps <- which(diff(as.numeric(d)) > 7)
  add  <- as.Date(character())

  for (k in gaps) {
    n_missing <- as.numeric(d[k + 1] - d[k]) / 7 - 1
    if (n_missing <= MAX_FILL_WEEKS) {
      add <- c(add, seq(d[k] + 7, d[k + 1] - 7, by = "week"))
    } else {
      say("LEFT UNFILLED: %s .. %s (%d weeks) - too long to be an ",
          format(d[k]), format(d[k + 1]), n_missing)
      say("               off-season zero run; treated as missing data.")
    }
  }
  if (length(add) == 0) return(df)

  say("filled %d omitted week(s) with zero: %s", length(add),
      paste(format(sort(add)), collapse = ", "))

  bind_rows(df, tibble(
    country                     = COUNTRY,
    age_group                   = "total",
    target_end_date             = add,
    year_week                   = sprintf("%d-W%02d", isoyear(add), isoweek(add)),
    week                        = isoweek(add),
    year                        = isoyear(add),
    weekly_rsv_hospitalisations = 0L
  )) %>% arrange(target_end_date)
}

admissions <- fill_weekly_gaps(admissions)

say("weekly admissions: %d weeks, %s .. %s",
    nrow(admissions), min(admissions$target_end_date),
    max(admissions$target_end_date))


# ---- 2. age-stratified burden -------------------------------------------
# The source gives PROPORTIONS over rolling 4-week periods. Multiply by
# the admissions actually observed in each period to recover counts.
#
# date_28days_floor is again a Monday, so the period runs from that
# Monday + 6 (the first week-ending date) for four weeks.
raw_age <- read_csv(file.path(SRC, "RSV_monthly_prop_age.csv"),
                    show_col_types = FALSE)

unknown <- setdiff(unique(raw_age$age_gp_modelling), names(AGE_MAP))
if (length(unknown) > 0) {
  stop("Unmapped age label(s) in RSV_monthly_prop_age.csv: ",
       paste(unknown, collapse = ", "), call. = FALSE)
}

periods <- raw_age %>%
  distinct(date_28days_floor) %>%
  mutate(start_date = as.Date(date_28days_floor) + 6,
         end_date   = start_date + weeks(3))

# Admissions falling inside each 4-week window
period_totals <- periods %>%
  mutate(total = map2_dbl(start_date, end_date, function(s, e) {
    sum(admissions$weekly_rsv_hospitalisations[
      admissions$target_end_date >= s & admissions$target_end_date <= e])
  }))

burden <- raw_age %>%
  mutate(start_date = as.Date(date_28days_floor) + 6) %>%
  left_join(period_totals %>% select(start_date, end_date, total),
            by = "start_date") %>%
  filter(!is.na(total)) %>%
  transmute(country   = COUNTRY,
            age_group = unname(AGE_MAP[age_gp_modelling]),
            start_date,
            end_date,
            # Carried through from the source. Without it the model would
            # derive a season label from each window's own start and end
            # years - fine for one window per season, but here it would
            # cut every seasonal summary into four-week slices.
            season    = season_name,
            total_rsv_hospitalisations = total * proportion) %>%
  arrange(start_date, age_group)

# The proportions must sum to 1 in every period, or the age split does
# not reconcile with the weekly totals it is applied to.
chk <- raw_age %>% group_by(date_28days_floor) %>%
  summarise(p = sum(proportion), .groups = "drop") %>%
  filter(abs(p - 1) > 1e-6)
if (nrow(chk) > 0) {
  stop("Age proportions do not sum to 1 in ", nrow(chk), " period(s), e.g. ",
       chk$date_28days_floor[1], " sums to ", signif(chk$p[1], 6),
       call. = FALSE)
}

# And the reconstructed counts must equal the weekly totals they came from
recon <- burden %>% group_by(start_date) %>%
  summarise(got = sum(total_rsv_hospitalisations), .groups = "drop") %>%
  left_join(period_totals %>% select(start_date, total), by = "start_date") %>%
  filter(abs(got - total) > 1e-6)
if (nrow(recon) > 0) {
  stop("Reconstructed age counts do not match the period totals in ",
       nrow(recon), " period(s).", call. = FALSE)
}

# A window must not straddle two seasons, or its weeks belong to both.
straddle <- burden %>% distinct(start_date, season) %>%
  count(start_date) %>% filter(n > 1)
if (nrow(straddle) > 0) {
  stop("Burden window(s) carry more than one season label: ",
       paste(format(straddle$start_date), collapse = ", "), call. = FALSE)
}

say("burden: %d rows, %d periods x %d bands, %s .. %s",
    nrow(burden), n_distinct(burden$start_date), n_distinct(burden$age_group),
    min(burden$start_date), max(burden$end_date))
say("seasons: %s", paste(sort(unique(burden$season)), collapse = ", "))
say("bands: %s", paste(sort(unique(burden$age_group)), collapse = ", "))


# ---- 3. births -----------------------------------------------------------
# Eurostat DEMO_FMONTH. Real data ends before 2024, so the 2023 monthly
# pattern is repeated forward - the same placeholder the previous Ireland
# round used. This drives dose counts only, not admissions.
raw_b <- read_csv(file.path(SRC, "country_monthly_births.csv"),
                  show_col_types = FALSE)

births_real <- raw_b %>%
  transmute(date   = make_date(as.integer(TIME_PERIOD),
                               match(month, month.name), 1L),
            births = as.numeric(OBS_VALUE)) %>%
  filter(!is.na(date), date < as.Date("2024-01-01"))

PROJECT_TO <- 2024:2027
projected <- map_dfr(PROJECT_TO, function(y) {
  births_real %>% filter(year(date) == 2023) %>%
    mutate(date = date %m+% years(y - 2023))
})

births <- bind_rows(births_real, projected) %>%
  transmute(country = ISO2,
            date,
            month   = month.name[month(date)],
            year    = year(date),
            births) %>%
  arrange(date)

stopifnot(!anyDuplicated(births$date))
say("births: %d months, %s .. %s (projected from %d)",
    nrow(births), min(births$date), max(births$date), 2023)


# ---- 4. population -------------------------------------------------------
# Aggregated from the repo's own RespiCompass population file, collapsing
# its 11 bands into the six the Irish data carries. That file is a
# 2026/27 estimate applied to a 2025/26 model - a small mismatch, but it
# only scales the per-100k impact tables and the adult dose denominator.
pop_src <- read_csv("data/population/population_estimates.csv",
                    show_col_types = FALSE) %>%
  filter(country == ISO2)

COLLAPSE <- c("0-2mo" = "0-2mo", "3-5mo" = "3-5mo", "6-11mo" = "6-11mo",
              "1-4" = "1-4y",
              "5-17" = "5-64y", "18-59" = "5-64y", "60-64" = "5-64y",
              "65-69" = "65+y", "70-74" = "65+y", "75-79" = "65+y",
              "80+" = "65+y")

missing <- setdiff(pop_src$age_group, names(COLLAPSE))
if (length(missing) > 0) {
  stop("population_estimates.csv has band(s) with no collapse rule: ",
       paste(missing, collapse = ", "), call. = FALSE)
}

population <- pop_src %>%
  mutate(age_group = unname(COLLAPSE[age_group])) %>%
  group_by(country, age_group) %>%
  summarise(population = sum(population), .groups = "drop") %>%
  arrange(match(age_group, unname(AGE_MAP)))

stopifnot(setequal(population$age_group, unname(AGE_MAP)))
say("population: %s", paste(sprintf("%s=%s", population$age_group,
                                    format(population$population, big.mark = ",")),
                            collapse = "  "))


# ---- 5. write ------------------------------------------------------------
write_csv(admissions, "data/epidemiological/hospitaladmissions_IE.csv")
write_csv(burden,     "data/epidemiological/hospitalburden_agegroups_IE.csv")
write_csv(births,     "data/population/births_by_month_IE.csv")
write_csv(population, "data/population/population_estimates_IE.csv")

cat("\nWrote:\n")
for (f in c("data/epidemiological/hospitaladmissions_IE.csv",
            "data/epidemiological/hospitalburden_agegroups_IE.csv",
            "data/population/births_by_month_IE.csv",
            "data/population/population_estimates_IE.csv")) {
  cat(sprintf("  %-56s %6d rows\n", f, nrow(read_csv(f, show_col_types = FALSE))))
}
cat("\nDONE\n")
