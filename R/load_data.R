# Load and parse the YAML configuration file.
# Returns the raw cfg list augmented with pre-parsed date objects and
# the waning data frame so callers never need to re-parse them.
load_config <- function(path = "config/static_model.yaml") {
  cfg <- yaml::read_yaml(path)

  cfg$vacc_start      <- ymd(cfg$vaccination_start)
  cfg$vacc_end        <- ymd(cfg$vaccination_end)
  cfg$anchor          <- as.Date(cfg$submission_horizon_anchor)
  cfg$births_cutoff   <- ymd(cfg$births_data_cutoff)

  # Waning data frame consumed by apply_scenario() via left_join on age_group
  cfg$waning_function <- data.frame(
    age_group = names(cfg$waning_by_band),
    waning    = unlist(cfg$waning_by_band),
    row.names = NULL
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
# The burden table is derived by multiplying the raw monthly age proportions
# by the weekly counts aggregated over the matching 4-week window, then
# translating RespiCompass age labels to the project's band convention.
load_epidemiological_data <- function(
    cfg,
    weekly_file  = "data/epidemiological/RSV_weekly_counts.csv",
    monthly_file = "data/epidemiological/RSV_monthly_prop_age.csv") {

  admissions <- read.csv(weekly_file) %>%
    mutate(target_end_date             = as.Date(date_wk_floor) + 6,
           weekly_rsv_hospitalisations = case_counts) %>%
    select(target_end_date, season_name, weekly_rsv_hospitalisations) %>%
    setDT()

  raw_age <- read.csv(monthly_file, fileEncoding = "UTF-8-BOM")

  # 4-week period boundaries derived from the date_28days_floor column
  periods <- raw_age %>%
    distinct(date_28days_floor) %>%
    mutate(period_start = as.Date(date_28days_floor) + 6,
           period_end   = as.Date(date_28days_floor) + 6 + weeks(3))

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
    mutate(date      = as.Date(date_28days_floor) + 6,
           age_group = age_gp_modelling,
           value     = NA) %>%
    select(date, age_group, value, proportion) %>%
    left_join(weekly_4wk, by = "date") %>%
    mutate(value                      = total_4wk * proportion,
           burden_start_date          = date,
           burden_end_date            = date + weeks(3),
           total_rsv_hospitalisations = value) %>%
    filter(!is.na(value)) %>%
    # Translate RespiCompass labels to the project's age-band convention
    mutate(age_group = case_when(
      age_group == "< 3 months"  ~ "0-2mo",
      age_group == "3-5 months"  ~ "3-5mo",
      age_group == "6-11 months" ~ "6-11mo",
      age_group == "1-4 years"   ~ "1-4y",
      age_group == "5-64 years"  ~ "5-64y",
      age_group == "65+ years"   ~ "65+y",
      TRUE ~ age_group)) %>%
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
