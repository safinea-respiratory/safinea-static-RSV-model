# Run the full pipeline for a single country.
#
# Everything downstream of data loading is country-agnostic, so this is
# a thin wrapper: load that country's data, sample, apply scenarios,
# format, and stamp the `location` column.
#
# Parameters:
#   cfg          – parsed config
#   country_name – full name, used to select rows from target-data
#   country_iso2 – two-letter code, used for auxiliary-data and emitted
#                  as the submission's `location`
#   adult_waning – the waning ensemble, loaded once and shared
#   quiet        – suppress the informational validation messages
#                  (they are identical for every country)
#
# Returns a list:
#   $submission      – hospitalisations + doses, with `location`
#   $baseline_df     – the raw Monte-Carlo samples
#   $adult_protection – named list of coverage/residual-VE per scenario
run_country <- function(cfg, country_name, country_iso2,
                        adult_waning, raw, quiet = FALSE) {

  epi    <- load_epidemiological_data(cfg, country_name, raw)
  validate_age_groups(cfg, epi, quiet = quiet)
  births <- load_births_data(cfg, country_iso2, raw)

  # ---- Baseline Monte-Carlo ----
  baseline_df <- simulate_weekly_age_fixed_margins(
    weekly_df = epi$admissions %>% filter(target_end_date   >= ymd(cfg$data_start)),
    age_df    = epi$burden     %>% filter(burden_start_date >= ymd(cfg$data_start)),
    n_draws   = cfg$n_draws,
    seed      = cfg$mc_seed
  )

  if (nrow(baseline_df) == 0) {
    stop("No baseline samples produced for ", country_name,
         ". Check that data_start (", cfg$data_start, ") precedes the ",
         "available data.", call. = FALSE)
  }

  # ---- Infant vaccination scenarios ----
  run_scenario <- function(uptake) {
    apply_scenario(
      df                   = baseline_df,
      IE_mean              = cfg$infant$vacc_IE_mean,
      IE_sd                = cfg$infant$vacc_IE_sd,
      vacc_uptake          = uptake,
      vacc_start           = cfg$infant$vacc_start,
      vacc_end             = cfg$infant$vacc_end,
      vacc_uptake_baseline = cfg$infant$baseline_uptake,
      waning_df            = cfg$infant$waning_df,
      age_bounds           = cfg$infant$age_bounds
    )
  }

  submission_pre <- assemble_submission(
    baseline_df   = baseline_df,
    scenario_A_df = run_scenario(cfg$infant$scenarios$no_vacc),
    scenario_B_df = run_scenario(cfg$infant$scenarios$high_vacc),
    round_id      = cfg$round_id,
    anchor        = cfg$anchor
  )

  # ---- Administered doses ----
  doses_df <- build_dose_table(
    baseline_df     = baseline_df,
    births_df       = births,
    vacc_start      = cfg$infant$vacc_start,
    vacc_end        = cfg$infant$vacc_end,
    cfg             = cfg,
    round_id        = cfg$round_id,
    anchor          = cfg$anchor,
    output_type_ids = unique(submission_pre$output_type_id)
  )

  # ---- Adult coverage & protection (computed, not yet applied) ----
  adult_prot <- function(total_coverage) {
    build_adult_protection(
      target_weeks        = unique(baseline_df$target_end_date),
      waning_curves       = adult_waning,
      campaigns           = cfg$adult$campaigns,
      total_coverage      = total_coverage,
      ve_beyond_curve     = cfg$adult$ve_beyond_curve,
      eligible_age_groups = cfg$adult$eligible_age_groups
    )
  }

  submission <- bind_rows(submission_pre, doses_df) %>%
    mutate(location = country_iso2)

  list(
    submission  = submission,
    baseline_df = baseline_df,
    adult_protection = list(
      baseline  = adult_prot(cfg$adult$baseline_coverage),
      no_vacc   = adult_prot(cfg$adult$scenarios$no_vacc),
      high_vacc = adult_prot(cfg$adult$scenarios$high_vacc)
    )
  )
}


# Run every country listed in the config and bind the results.
#
# Countries are independent: each gets its own data, its own sampler
# draws and its own submission rows, distinguished by `location`. The
# same mc_seed is used throughout, so the vaccine-effectiveness draws
# are shared across countries - defensible, since it is the same vaccine.
run_all_countries <- function(cfg, adult_waning, raw) {

  cs <- cfg$countries_df

  results <- map(seq_len(nrow(cs)), function(i) {
    message("Running ", cs$name[i], " (", cs$iso2[i], ") ",
            "[", i, "/", nrow(cs), "]")
    run_country(cfg, cs$name[i], cs$iso2[i], adult_waning, raw,
                quiet = (i > 1))
  })
  names(results) <- cs$iso2

  list(
    submission = bind_rows(map(results, "submission")) %>% as.data.table(),
    by_country = results
  )
}
