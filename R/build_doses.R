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

  # ---- adult doses, per age band ----
  # Coverage is applied uniformly across a scenario's targeted bands, so
  # each band gets that week's coverage increment times its OWN
  # population. Summing the bands first would throw away exactly the
  # breakdown we want.
  adult_doses_for <- function(coverage, bands) {
    pops <- population_df %>%
      filter(age_group %in% bands) %>%
      select(age_group, population)
    if (nrow(pops) == 0) {
      return(tibble(target_end_date = as.Date(character()),
                    age_group = character(), doses = numeric()))
    }
    sched <- build_campaign_schedule(cfg$adult$campaigns, coverage, weeks_vec)
    crossing(tibble(target_end_date = weeks_vec), pops) %>%
      left_join(sched %>% rename(target_end_date = week), by = "target_end_date") %>%
      mutate(doses = coalesce(delta, 0) * population) %>%
      select(target_end_date, age_group, doses)
  }

  sc <- cfg$scenarios_df

  # Every band any scenario vaccinates. The dose grid must be the same
  # for every scenario - a scenario that does not target a band still
  # emits a zero row for it - or the submission grid would be ragged.
  dose_bands <- sort(unique(c(
    unlist(lapply(seq_len(nrow(sc)), function(i)
      if (sc$adult_coverage[i] > 0) sc$adult_age_groups[[i]] else character(0))),
    if (any(sc$infant_uptake > 0)) cfg$infant$dose_age_group else character(0)
  )))

  if (length(dose_bands) == 0) {
    stop("No scenario vaccinates anyone, so there are no dose age groups ",
         "to report.", call. = FALSE)
  }

  # One dose track per scenario: infant births x infant uptake attributed
  # to the birth cohort, plus adult population x coverage increment for
  # each targeted band.
  per_band <- bind_rows(lapply(seq_len(nrow(sc)), function(i) {

    inf <- if (sc$infant_uptake[i] > 0) {
      infant_base %>%
        transmute(target_end_date,
                  age_group = cfg$infant$dose_age_group,
                  doses     = births * sc$infant_uptake[i])
    } else tibble()

    adu <- adult_doses_for(sc$adult_coverage[i], sc$adult_age_groups[[i]])

    # Fill the full band grid so every scenario carries every band.
    crossing(target_end_date = weeks_vec, age_group = dose_bands) %>%
      left_join(bind_rows(inf, adu) %>%
                  group_by(target_end_date, age_group) %>%
                  summarise(doses = sum(doses), .groups = "drop"),
                by = c("target_end_date", "age_group")) %>%
      mutate(value       = coalesce(doses, 0),
             scenario_id = sc$id[i]) %>%
      select(target_end_date, age_group, value, scenario_id)
  }))

  # An all-ages row alongside the per-band ones, mirroring the
  # total_imm* rows on the hospitalisation side. Anything summing doses
  # must use one or the other, never both.
  totals <- per_band %>%
    group_by(target_end_date, scenario_id) %>%
    summarise(value = sum(value), .groups = "drop") %>%
    mutate(age_group = "total")

  crossing(
    bind_rows(per_band, totals) %>%
      transmute(target_end_date, value,
                target    = "administered_doses",
                pop_group = age_group,
                scenario_id),
    round_id       = round_id,
    output_type    = "sample",
    output_type_id = output_type_ids
  ) %>%
    mutate(horizon = as.integer((target_end_date - anchor) / 7))
}
