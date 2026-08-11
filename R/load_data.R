# Parse a character date column with an explicit format, strictly.
#
# R's as.Date() is dangerous here: given "01/09/2025" with no format it
# returns 0001-09-20 rather than NA, so a mis-formatted input silently
# produces dates ~2000 years out and every downstream join quietly
# yields zero rows. This wrapper requires the format up front, rejects
# anything that fails to parse, and additionally rejects implausible
# years so that class of error can never pass silently again.
parse_dates_strict <- function(x, format, what, min_year = 1900, max_year = 2100) {

  d <- as.Date(as.character(x), format = format)

  if (any(is.na(d))) {
    bad <- unique(as.character(x)[is.na(d)])
    stop("Could not parse ", what, " with format '", format, "'.",
         "\n  Unparseable value(s): ",
         paste0('"', utils::head(bad, 5), '"', collapse = ", "),
         if (length(bad) > 5) " ..." else "",
         "\n  Fix the data, or correct input_date_formats in the config.",
         call. = FALSE)
  }

  yr <- as.integer(format(d, "%Y"))
  if (any(yr < min_year | yr > max_year)) {
    bad_i <- which(yr < min_year | yr > max_year)
    stop("Parsing ", what, " with format '", format, "' produced ",
         "implausible date(s) - the format is almost certainly wrong.",
         "\n  e.g. \"", as.character(x)[bad_i[1]], "\" -> ", format(d[bad_i[1]]),
         "\n  Expected years between ", min_year, " and ", max_year, ".",
         call. = FALSE)
  }

  d
}


# Rename source age-group labels using the optional alias map.
# Labels with no alias pass through unchanged.
apply_age_aliases <- function(x, aliases) {
  if (is.null(aliases) || length(aliases) == 0) return(x)
  m   <- unlist(aliases)
  out <- unname(m[x])
  ifelse(is.na(out), x, out)
}


# Load and parse the YAML configuration file.
#
# Returns the raw cfg list augmented with pre-parsed dates and two
# normalised sub-lists, cfg$infant and cfg$adult, so that callers never
# reach into the nested YAML structure or re-parse dates.
#
# The infant and adult programmes are kept strictly separate: they share
# no parameters, and each sub-list is self-contained.
load_config <- function(path = "config/static_model.yaml") {
  cfg <- yaml::read_yaml(path)

  cfg$anchor        <- as.Date(cfg$submission_horizon_anchor)
  cfg$births_cutoff <- ymd(cfg$births_data_cutoff)

  # ---- Infant programme (maternal / birth-dose) ----
  inf <- cfg$infant_vaccination

  # Waning indexed by AGE BAND, because for this product age equals
  # time since vaccination.
  waning_df <- data.frame(
    age_group = names(inf$waning_by_band),
    waning    = unlist(inf$waning_by_band, use.names = FALSE),
    row.names = NULL
  )

  # Month bounds per band, [min inclusive, max exclusive).
  age_bounds <- data.frame(
    age_group = names(inf$age_bounds),
    min_mo    = vapply(inf$age_bounds, function(b) as.numeric(b[[1]]), numeric(1)),
    max_mo    = vapply(inf$age_bounds, function(b) as.numeric(b[[2]]), numeric(1)),
    row.names = NULL
  )

  cfg$infant <- list(
    vacc_IE_mean    = inf$vacc_IE$mean,
    vacc_IE_sd      = inf$vacc_IE$sd,
    waning_df       = waning_df,
    age_bounds      = age_bounds,
    vacc_start      = ymd(unlist(inf$windows$start)),
    vacc_end        = ymd(unlist(inf$windows$end)),
    baseline_uptake = inf$baseline_uptake,
    scenarios       = inf$scenarios
  )

  # ---- Adult programme (calendar campaign) ----
  # No VE level parameter here by design: the waning ensemble carries
  # the uncertainty, one whole curve per Monte-Carlo sample.
  adu <- cfg$adult_vaccination

  # Campaigns: normalise dates and fill in equal shares when omitted.
  campaigns <- lapply(adu$campaigns, function(cp) {
    list(start   = as.Date(cp$start),
         end     = as.Date(cp$end),
         share   = if (is.null(cp$share)) NA_real_ else as.numeric(cp$share),
         profile = if (is.null(cp$profile)) "uniform" else cp$profile)
  })
  if (length(campaigns) > 0 && all(vapply(campaigns, function(c) is.na(c$share), logical(1)))) {
    eq <- 1 / length(campaigns)
    campaigns <- lapply(campaigns, function(c) { c$share <- eq; c })
  }

  cfg$adult <- list(
    eligible_age_groups = unlist(adu$eligible_age_groups),
    waning_curves_path  = adu$waning_curves,
    ve_target           = adu$ve_target,
    ve_beyond_curve     = adu$ve_beyond_curve,
    campaigns           = campaigns,
    baseline_coverage   = adu$baseline_coverage,
    scenarios           = adu$scenarios
  )

  cfg
}


# Load observed RSV weekly admissions and the monthly age-split proportions.
#
# Returns a named list:
#   $admissions — weekly aggregate counts
#                 (target_end_date, season_name, weekly_rsv_hospitalisations)
#   $burden     — per-age 4-week counts
#                 (burden_start_date, burden_end_date, age_group,
#                  total_rsv_hospitalisations)
#
# The burden table is derived by multiplying the raw monthly age
# proportions by the weekly counts aggregated over the matching 4-week
# window. Age-group labels come from the data itself; the only
# transformation applied is the optional alias map from the config.
load_epidemiological_data <- function(
    cfg,
    weekly_file  = "data/epidemiological/RSV_weekly_counts.csv",
    monthly_file = "data/epidemiological/RSV_monthly_prop_age.csv") {

  fmt_weekly  <- cfg$input_date_formats$weekly_counts
  fmt_monthly <- cfg$input_date_formats$monthly_age

  admissions <- read.csv(weekly_file) %>%
    mutate(target_end_date             = parse_dates_strict(
                                           date_wk_floor, fmt_weekly,
                                           "RSV_weekly_counts.csv$date_wk_floor") + 6,
           weekly_rsv_hospitalisations = case_counts) %>%
    select(target_end_date, season_name, weekly_rsv_hospitalisations) %>%
    setDT()

  raw_age <- read.csv(monthly_file, fileEncoding = "UTF-8-BOM") %>%
    mutate(.period_date = parse_dates_strict(
                            date_28days_floor, fmt_monthly,
                            "RSV_monthly_prop_age.csv$date_28days_floor"))

  # 4-week period boundaries derived from the date_28days_floor column
  periods <- raw_age %>%
    distinct(.period_date) %>%
    mutate(period_start = .period_date + 6,
           period_end   = .period_date + 6 + weeks(3))

  # Aggregate weekly counts within each 4-week window so proportions can
  # be applied to recover per-age-band counts
  weekly_4wk <- admissions %>%
    crossing(periods) %>%
    filter(target_end_date >= period_start & target_end_date <= period_end) %>%
    group_by(period_start) %>%
    summarise(total_4wk = sum(weekly_rsv_hospitalisations, na.rm = TRUE),
              .groups = "drop") %>%
    rename(date = period_start)

  burden <- raw_age %>%
    mutate(date      = .period_date + 6,
           age_group = apply_age_aliases(age_gp_modelling, cfg$age_group_aliases),
           value     = NA) %>%
    select(date, age_group, value, proportion) %>%
    left_join(weekly_4wk, by = "date") %>%
    mutate(value                      = total_4wk * proportion,
           burden_start_date          = date,
           burden_end_date            = date + weeks(3),
           total_rsv_hospitalisations = value) %>%
    filter(!is.na(value)) %>%
    select(burden_start_date, burden_end_date, age_group, total_rsv_hospitalisations) %>%
    setDT()

  list(admissions = admissions, burden = burden)
}


# Load monthly births and project the 2023 seasonal pattern forward
# for years where real data is unavailable (births_project_years).
#
# Births data is only available pre-2024 at time of writing.
# The 2023 monthly pattern is repeated as a placeholder that assumes
# year-on-year stability in the seasonal birth distribution.
# Replace with real projected births when they become available.
load_births_data <- function(
    cfg,
    births_file = "data/population/country_monthly_births.csv") {

  births_raw <- read.csv(births_file, fileEncoding = "UTF-8-BOM") %>%
    mutate(date    = make_date(year  = as.integer(TIME_PERIOD),
                               month = match(month, month.name),
                               day   = 1),
           country = geo,
           births  = OBS_VALUE) %>%
    select(country, date, births) %>%
    filter(!is.na(date), date < cfg$births_cutoff) %>%
    setDT()

  projected <- map_dfr(
    cfg$births_project_years,
    ~ births_raw %>%
        filter(year(date) == 2023) %>%
        mutate(date = date + years(.x - 2023))
  )

  bind_rows(births_raw, projected) %>% arrange(country, date)
}


# Load the ADULT vaccine waning ensemble.
#
# The source file holds one complete VE-over-time curve per `rep`. Each
# Monte-Carlo sample takes one entire curve, so within-sample
# correlation across months is preserved by construction and no
# distributional assumption is required. This is deliberately unlike the
# infant programme, whose uncertainty is a parametric draw on vacc_IE.
#
# Reps are independent realisations, so the first n_draws are used.
# Errors if n_draws exceeds the ensemble size.
#
# Returns a data.frame(sample, month, ve) where `sample` is 1..n_draws
# and `ve` is the column named by cfg$adult$ve_target.
load_waning_curves <- function(cfg, n_draws) {

  path      <- cfg$adult$waning_curves_path
  ve_target <- cfg$adult$ve_target

  if (!file.exists(path)) {
    stop("Adult waning curve file not found: ", path, call. = FALSE)
  }

  curves <- read.csv(path)

  required <- c("rep", "month", ve_target)
  absent   <- setdiff(required, names(curves))
  if (length(absent) > 0) {
    stop("Waning curve file is missing required column(s): ",
         paste(absent, collapse = ", "),
         "\n  file: ", path,
         "\n  columns found: ", paste(names(curves), collapse = ", "),
         call. = FALSE)
  }

  reps   <- sort(unique(curves$rep))
  n_reps <- length(reps)

  if (n_draws > n_reps) {
    stop("n_draws (", n_draws, ") exceeds the number of reps in the adult ",
         "waning ensemble (", n_reps, ").",
         "\n  Reduce n_draws, or supply an ensemble with more reps.",
         "\n  file: ", path,
         call. = FALSE)
  }

  # Every rep must supply a complete, gap-free curve starting at month 0
  per_rep <- curves %>%
    group_by(rep) %>%
    summarise(n_months   = n(),
              min_month  = min(month),
              n_distinct = n_distinct(month),
              .groups    = "drop")

  expected_months <- max(per_rep$n_months)
  bad <- per_rep %>%
    filter(min_month != 0 | n_months != n_distinct | n_months != expected_months)

  if (nrow(bad) > 0) {
    stop("Waning curves are malformed for rep(s): ",
         paste(utils::head(bad$rep, 10), collapse = ", "),
         if (nrow(bad) > 10) " ..." else "",
         "\n  Each rep must start at month 0 and carry the same number of ",
         "unique, non-duplicated months.",
         "\n  file: ", path,
         call. = FALSE)
  }

  ve_vals <- curves[[ve_target]]
  if (any(is.na(ve_vals)) || any(ve_vals < 0 | ve_vals > 1)) {
    stop("Column '", ve_target, "' contains NA or values outside [0, 1].",
         "\n  file: ", path, call. = FALSE)
  }

  keep <- reps[seq_len(n_draws)]

  curves %>%
    filter(rep %in% keep) %>%
    mutate(sample = match(rep, keep),
           ve     = .data[[ve_target]]) %>%
    select(sample, month, ve) %>%
    arrange(sample, month)
}
