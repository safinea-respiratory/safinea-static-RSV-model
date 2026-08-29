# Build the INFANT administered-doses table for all scenarios.
#
# Monthly births are spread evenly across their days, masked to the
# vaccination windows, then re-aggregated to ISO weeks and scaled by each
# scenario's uptake.
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
#   births_df       – output of load_births_data(); one row per month, with
#                     `date` the first of the month
#   vacc_start/end  – Date vectors (one element per season)
#   cfg             – parsed config; uses cfg$infant$baseline_uptake and
#                     cfg$infant$scenarios
#   round_id        – submission round identifier
#   anchor          – submission horizon anchor date
#   output_type_ids – character vector of output_type_id values from the
#                     hospitalisation submission table
build_dose_table <- function(baseline_df, births_df,
                             vacc_start, vacc_end,
                             cfg, round_id, anchor, output_type_ids) {

  dates_df <- baseline_df %>%
    mutate(target_end_date = as.Date(target_end_date)) %>%
    distinct(target_end_date) %>%
    arrange(target_end_date)

  # Spread each month's births evenly across its days.
  #
  # Assigning a whole month to the single ISO week containing its 1st -
  # as this did previously - leaves three weeks in four empty and silently
  # discards any month whose week falls outside the modelled grid. On the
  # Ireland 2026/27 data that lost 8.5 % of births and left only 5 of 51
  # weeks non-zero.
  daily <- map_dfr(seq_len(nrow(births_df)), function(i) {
    m  <- births_df[i, ]
    nd <- as.integer(days_in_month(m$date))
    tibble(day        = seq(m$date, by = "day", length.out = nd),
           births_day = m$births / nd)
  })

  # Mask per DAY, so a week straddling a window boundary is credited only
  # for the days actually inside the window.
  weekly <- daily %>%
    mutate(in_window       = map_lgl(day, ~ any(.x >= vacc_start & .x <= vacc_end)),
           births_day      = ifelse(in_window, births_day, 0),
           target_end_date = floor_date(day, "week", week_start = 1) + 6) %>%
    group_by(target_end_date) %>%
    summarise(births = sum(births_day), .groups = "drop")

  # Births in a vaccination window but outside the modelled weeks would
  # otherwise disappear silently in the join below.
  off_grid <- weekly %>% filter(!target_end_date %in% dates_df$target_end_date)
  if (sum(off_grid$births) > 1e-9) {
    warning("build_dose_table: ", format(round(sum(off_grid$births))),
            " in-window births fall outside the modelled weeks and are ",
            "excluded from the dose table.", call. = FALSE)
  }

  base_doses <- dates_df %>%
    left_join(weekly, by = "target_end_date") %>%
    mutate(target    = "administered_doses",
           pop_group = "undefined",
           value     = coalesce(births, 0)) %>%
    select(target_end_date, value, target, pop_group)

  # One dose track per configured scenario, scaled by that scenario's
  # infant uptake. Adult doses are not included yet - they need the
  # population denominator (step 5).
  sc <- cfg$scenarios_df

  crossing(
    bind_rows(lapply(seq_len(nrow(sc)), function(i) {
      base_doses %>% mutate(value       = value * sc$infant_uptake[i],
                            scenario_id = sc$id[i])
    })),
    round_id       = round_id,
    output_type    = "sample",
    output_type_id = output_type_ids
  ) %>%
    mutate(horizon = as.integer((target_end_date - anchor) / 7))
}
