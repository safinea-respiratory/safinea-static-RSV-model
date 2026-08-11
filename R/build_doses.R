# Build the INFANT administered-doses table for all scenarios.
#
# Weekly dose counts are derived from monthly births aligned to
# ISO-week-ending Sundays, zeroed outside the vaccination windows, then
# scaled by each scenario's uptake fraction.
#
# The adult programme will need its own dose track, driven by campaign
# coverage against an adult population denominator rather than births.
#
# The resulting table is crossed with the submission's output_type_id grid
# so the dose rows share the same sample-index structure as the
# hospitalisation rows.
#
# Parameters:
#   baseline_df     – output of simulate_weekly_age_fixed_margins, used to
#                     determine the target_end_date grid for the dose table
#   births_df       – output of load_births_data()
#   vacc_start/end  – Date vectors (one element per season)
#   cfg             – parsed config list (needs baseline_uptake, scenarios,
#                     round_id)
#   anchor          – submission horizon anchor date
#   output_type_ids – character vector of output_type_id values from the
#                     hospitalisation submission table
build_dose_table <- function(baseline_df, births_df,
                             vacc_start, vacc_end,
                             cfg, anchor, output_type_ids) {

  dates_df <- baseline_df %>%
    mutate(target_end_date = as.Date(target_end_date)) %>%
    distinct(target_end_date) %>%
    arrange(target_end_date)

  # Align monthly birth counts to ISO week-ending Sunday, then zero out
  # any weeks that fall outside the vaccination windows
  base_doses <- dates_df %>%
    left_join(
      births_df %>%
        mutate(target_end_date = floor_date(as.Date(date), "week", week_start = 1) + 6) %>%
        select(country, target_end_date, births),
      by = "target_end_date"
    ) %>%
    mutate(
      target    = "administered_doses",
      pop_group = "undefined",
      value     = ifelse(is.na(births), 0, births),
      in_window = map_lgl(target_end_date, ~ any(.x >= vacc_start & .x <= vacc_end)),
      value     = ifelse(in_window, value, 0)
    ) %>%
    select(target_end_date, value, target, pop_group)

  scale_doses <- function(uptake, scenario_id) {
    base_doses %>% mutate(value = value * uptake, scenario_id = scenario_id)
  }

  crossing(
    bind_rows(
      scale_doses(cfg$infant$baseline_uptake,     "baseline"),
      scale_doses(cfg$infant$scenarios$no_vacc,   "no_vacc"),
      scale_doses(cfg$infant$scenarios$high_vacc, "high_vacc")
    ),
    round_id       = cfg$round_id,
    output_type    = "sample",
    output_type_id = output_type_ids
  ) %>%
    mutate(horizon = as.integer((target_end_date - anchor) / 7))
}
