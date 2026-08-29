# Add age-aggregated totals per immunisation status.
# RespiCompass expects a total_<imm_status> pop_group per
# (round, scenario, target, horizon, output_type_id) in addition to
# the per-age rows, so the hub can compute all-ages summaries.
add_age_totals <- function(df) {
  totals <- df %>%
    group_by(round_id, scenario_id, target, horizon, target_end_date,
             output_type, output_type_id, imm_status) %>%
    summarise(value = sum(value), .groups = "drop") %>%
    mutate(pop_group = paste0("total_", imm_status))
  bind_rows(df, totals)
}


# Reshape a scenario result from apply_scenario() into the column
# layout needed by assemble_submission().
scenario_to_submission_shape <- function(df, scenario_id) {
  df %>%
    select(target_end_date, age_group, sample, value, immunisation) %>%
    mutate(scenario = scenario_id)
}


# Combine scenario results into the RespiCompass submission format.
#
# scenario_results: named list of apply_scenario() outputs, keyed by the
#   scenario_id to emit. Any number of scenarios is supported; the set is
#   driven entirely by the `scenarios` block in the config.
#
# A zero-coverage scenario needs no special handling: apply_scenario()
# already gives it yes = 0 and no = total = observed.
#
# Final columns:
#   round_id, scenario_id, target, pop_group, horizon,
#   target_end_date, output_type, output_type_id, value
assemble_submission <- function(scenario_results, round_id, anchor) {

  imap_dfr(scenario_results, ~ scenario_to_submission_shape(.x, .y)) %>%
    mutate(target         = "rsv_hospitalisations",
           horizon        = as.integer((target_end_date - anchor) / 7),
           round_id       = round_id,
           output_type_id = as.character(sample),
           output_type    = "sample") %>%
    rename(scenario_id = scenario) %>%
    select(round_id, scenario_id, target, age_group, horizon, target_end_date,
           output_type, output_type_id, value, immunisation) %>%
    mutate(imm_status = case_when(
             immunisation == "yes"   ~ "immYes",
             immunisation == "no"    ~ "immNo",
             immunisation == "total" ~ "immTotal"),
           pop_group = paste0(age_group, "_", imm_status)) %>%
    select(-immunisation, -age_group) %>%
    add_age_totals() %>%
    select(-imm_status)
}
