# ==============================================================
# make_two_seasons.R
#
# Builds a TWO-SEASON version of the epidemiological and births
# inputs by duplicating the single upstream season forward one year.
#
# The RespiCompass target data covers one season (2026/27). To exercise
# the model over a multi-season horizon - in particular to see a single
# adult campaign's protection wane across a following season - season 2
# is an exact copy of season 1 shifted forward.
#
# THE OUTPUT IS DERIVED, NOT REAL. Season 2 is synthetic. The generated
# files carry a `season` column so the two are never confused, and the
# originals are left untouched.
#
# Two different shifts are used, deliberately:
#
#   weekly + burden   +364 days (exactly 52 weeks)
#       Preserves the weekday, so every date stays a week-ending Sunday
#       and season 2's first week lands exactly 7 days after season 1's
#       last. A calendar-year shift would not: 2028 is a leap year, so
#       +1 year is 365 or 366 days and would break the weekly grid.
#
#   births            +1 calendar year
#       Preserves the first-of-month alignment that days_in_month()
#       relies on in build_dose_table(). A 364-day shift would land
#       2026-09-01 on 2027-08-31, which is not a month start.
#
# Run from the project root:  Rscript scripts/make_two_seasons.R
# ==============================================================

suppressMessages({
  library(dplyr); library(tidyr); library(lubridate); library(purrr)
})

SHIFT_DAYS <- 364L   # 52 weeks

adm_in <- "data/epidemiological/hospitaladmissions.csv"
bur_in <- "data/epidemiological/hospitalburden_agegroups.csv"
bir_in <- "data/population/births_by_month.csv"

adm_out <- "data/epidemiological/hospitaladmissions_2seasons.csv"
bur_out <- "data/epidemiological/hospitalburden_agegroups_2seasons.csv"
bir_out <- "data/population/births_by_month_2seasons.csv"

stopifnot(file.exists(adm_in), file.exists(bur_in), file.exists(bir_in))

adm <- read.csv(adm_in) %>% mutate(target_end_date = as.Date(target_end_date))
bur <- read.csv(bur_in) %>% mutate(start_date = as.Date(start_date),
                                   end_date   = as.Date(end_date))
bir <- read.csv(bir_in) %>% mutate(date = as.Date(date))

# ---- season labels, from the burden window ----
s1_start <- min(bur$start_date); s1_end <- max(bur$end_date)
s2_start <- s1_start + SHIFT_DAYS; s2_end <- s1_end + SHIFT_DAYS
lab <- function(a, b) paste0(format(a, "%Y"), "/", format(b, "%Y"))
S1 <- lab(s1_start, s1_end); S2 <- lab(s2_start, s2_end)

cat("season 1:", S1, format(s1_start), "->", format(s1_end), "\n")
cat("season 2:", S2, format(s2_start), "->", format(s2_end), "\n")
cat("shift    :", SHIFT_DAYS, "days =", SHIFT_DAYS / 7, "weeks\n\n")

# ---- weekly admissions ----
# year_week / week / year are recomputed rather than shifted, so the
# ISO labels stay consistent with the new dates. They are unused by the
# model but would otherwise contradict target_end_date.
shift_weekly <- function(d, season) {
  d %>%
    mutate(target_end_date = target_end_date + if (season == S1) 0L else SHIFT_DAYS,
           week      = isoweek(target_end_date),
           year      = isoyear(target_end_date),
           year_week = sprintf("%d-W%02d", year, week),
           season    = season)
}
adm2 <- bind_rows(shift_weekly(adm, S1), shift_weekly(adm, S2)) %>%
  arrange(country, target_end_date) %>%
  select(country, age_group, target_end_date, year_week, week, year,
         weekly_rsv_hospitalisations, season)

# ---- age burden ----
shift_burden <- function(d, season) {
  off <- if (season == S1) 0L else SHIFT_DAYS
  d %>% mutate(start_date = start_date + off,
               end_date   = end_date   + off,
               season     = season)
}
bur2 <- bind_rows(shift_burden(bur, S1), shift_burden(bur, S2)) %>%
  arrange(country, start_date, age_group) %>%
  select(country, age_group, start_date, end_date,
         total_rsv_hospitalisations, season)

# ---- births ----
# No season column: births are a continuous monthly series consumed by
# date, and a month straddling a season boundary has no clean label.
bir2 <- bind_rows(bir, bir %>% mutate(date = date %m+% years(1),
                                      year = year + 1L)) %>%
  arrange(country, date)

# ============================================================ #
# Validate before writing - a silently mis-shifted row would be very
# hard to notice downstream.
# ============================================================ #
fail <- function(...) stop("make_two_seasons: ", ..., call. = FALSE)

# weeks: contiguous 7-day grid, no overlap, seasons join exactly
for (cc in unique(adm2$country)) {
  d <- adm2 %>% filter(country == cc) %>% arrange(target_end_date)
  gaps <- unique(as.numeric(diff(d$target_end_date)))
  if (!identical(gaps, 7)) {
    fail("country ", cc, " has non-7-day gaps: ", paste(gaps, collapse = ", "))
  }
  if (anyDuplicated(d$target_end_date) > 0) fail("duplicate weeks for ", cc)
}
n_wk <- adm2 %>% count(country, season)
if (length(unique(n_wk$n)) != 1) fail("uneven week counts across country/season")

j <- adm2 %>% group_by(season) %>%
  summarise(first = min(target_end_date), last = max(target_end_date), .groups = "drop")
if (j$first[j$season == S2] != j$last[j$season == S1] + 7) {
  fail("season 2 does not start exactly one week after season 1 ends")
}

# burden windows must not overlap
if (s2_start <= s1_end) fail("burden windows overlap at the season boundary")

# season 2 must reproduce season 1's totals exactly
t1 <- adm2 %>% filter(season == S1) %>% group_by(country) %>%
  summarise(v = sum(weekly_rsv_hospitalisations), .groups = "drop")
t2 <- adm2 %>% filter(season == S2) %>% group_by(country) %>%
  summarise(v = sum(weekly_rsv_hospitalisations), .groups = "drop")
if (!isTRUE(all.equal(t1$v, t2$v))) fail("season totals differ between seasons")

# margins must still reconcile within each season
rec <- adm2 %>% group_by(country, season) %>%
  summarise(wk = sum(weekly_rsv_hospitalisations), .groups = "drop") %>%
  inner_join(bur2 %>% group_by(country, season) %>%
               summarise(bu = sum(total_rsv_hospitalisations), .groups = "drop"),
             by = c("country", "season")) %>%
  mutate(d = wk - bu)
if (any(rec$d != 0)) {
  fail(sum(rec$d != 0), " country-seasons where weekly and burden totals disagree")
}

# every week must fall inside its own season's burden window
win <- bur2 %>% distinct(season, start_date, end_date)
chk <- adm2 %>% inner_join(win, by = "season") %>%
  filter(target_end_date < start_date | target_end_date > end_date)
if (nrow(chk) > 0) fail(nrow(chk), " week(s) outside their season's burden window")

# births: contiguous months, doubled, no duplicates
for (cc in unique(bir2$country)) {
  d <- bir2 %>% filter(country == cc) %>% arrange(date)
  if (anyDuplicated(d$date) > 0) fail("duplicate birth months for ", cc)
  if (!all(day(d$date) == 1)) fail("births not on month starts for ", cc)
  if (nrow(d) != 2 * nrow(bir %>% filter(country == cc))) {
    fail("births not doubled for ", cc)
  }
}

# ---- write ----
write.csv(adm2, adm_out, row.names = FALSE)
write.csv(bur2, bur_out, row.names = FALSE)
write.csv(bir2, bir_out, row.names = FALSE)

cat("weekly : ", nrow(adm2), " rows -> ", adm_out, "\n", sep = "")
cat("burden : ", nrow(bur2), " rows -> ", bur_out, "\n", sep = "")
cat("births : ", nrow(bir2), " rows -> ", bir_out, "\n", sep = "")
cat("\nweeks per country per season: ", unique(n_wk$n), "\n", sep = "")
cat("all validation checks passed\n")
