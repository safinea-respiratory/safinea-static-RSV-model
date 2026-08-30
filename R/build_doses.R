# Build the administered-doses table for all scenarios.
#
# Doses come from BOTH programmes, each with its own denominator:
#
#   infant – monthly births, spread evenly across their days, masked to
#            the vaccination windows per day, re-aggregated to ISO weeks
#            and scaled by the scenario's infant uptake.
#   adult  – the eligible age bands' population multiplied by that week's
#            NEW coverage increment from the campaign schedule. Because
#            coverage is one-off and cumulative, the weekly increment is
#            exactly the number of people newly vaccinated, which is what
#            a dose count means.
#
# The two are summed per week. Doses are deterministic - neither
# programme's dose count depends on the Monte-Carlo draw - so the result
# is crossed with the output_type_id grid to match the hospitalisation
# rows' shape.
#
# Parameters:
#   baseline_df     – output of simulate_weekly_age_fixed_margins, used to
#                     determine the target_end_date grid for the dose table
#   births_df       – output of load_births_data(); one row per month, with
#                     `date` the first of the month
#   population_df   – output of load_population_data(); population by age band
#   vacc_start/end  – Date vectors (one element per season), infant windows
#   cfg             – parsed config; uses cfg$scenarios_df, cfg$adult
#   round_id        – submission round identifier
#   anchor          – submission horizon anchor date
#   output_type_ids – character vector of output_type_id values from the
#                     hospitalisation submission table
build_dose_table <- function(baseline_df, births_df, population_df,
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

  infant_base <- dates_df %>%
    left_join(weekly, by = "target_end_date") %>%
    mutate(births = coalesce(births, 0)) %>%
    select(target_end_date, births)

  weeks_vec <- dates_df$target_end_date

  # ---- adult denominator ----
  # Coverage is applied uniformly across a scenario's targeted bands, so
  # the denominator is their combined population. It is computed PER
  # SCENARIO because scenarios may target different bands - retargeting
  # from 65+ to 75+ shrinks the denominator as well as the effect.
  adult_doses_for <- function(coverage, bands) {
    adult_pop <- population_df %>%
      filter(age_group %in% bands) %>%
      summarise(p = sum(population)) %>%
      pull(p)

    sched <- build_campaign_schedule(cfg$adult$campaigns, coverage, weeks_vec)
    tibble(target_end_date = weeks_vec) %>%
      left_join(sched %>% rename(target_end_date = week), by = "target_end_date") %>%
      mutate(adult = coalesce(delta, 0) * adult_pop) %>%
      select(target_end_date, adult)
  }

  # One dose track per configured scenario: infant births x infant uptake,
  # plus adult population x that week's coverage increment.
  sc <- cfg$scenarios_df

  crossing(
    bind_rows(lapply(seq_len(nrow(sc)), function(i) {
      infant_base %>%
        left_join(adult_doses_for(sc$adult_coverage[i], sc$adult_age_groups[[i]]),
                  by = "target_end_date") %>%
        mutate(value       = births * sc$infant_uptake[i] + adult,
               target      = "administered_doses",
               pop_group   = "undefined",
               scenario_id = sc$id[i]) %>%
        select(target_end_date, value, target, pop_group, scenario_id)
    })),
    round_id       = round_id,
    output_type    = "sample",
    output_type_id = output_type_ids
  ) %>%
    mutate(horizon = as.integer((target_end_date - anchor) / 7))
}
