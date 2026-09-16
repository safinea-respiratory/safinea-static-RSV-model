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
  popn   <- load_population_data(cfg, country_iso2, raw)

  # ---- Baseline Monte-Carlo ----
  baseline_df <- simulate_weekly_age_fixed_margins(
    weekly_df = epi$admissions %>% filter(target_end_date   >= ymd(cfg$data_start)),
    age_df    = epi$burden     %>% filter(burden_start_date >= ymd(cfg$data_start)),
    n_draws   = cfg$n_draws,
    seed      = cfg$mc_seed
  )

  if (nrow(baseline_df) == 0) {
    n_wk  <- sum(epi$admissions$target_end_date   >= ymd(cfg$data_start))
    n_bur <- sum(epi$burden$burden_start_date >= ymd(cfg$data_start))
    stop("No baseline samples produced for ", country_name, ".",
         "\n  After filtering on data_start (", cfg$data_start, "): ",
         n_wk, " weekly row(s), ", n_bur, " burden row(s).",
         "\n  Unfiltered the data spans:",
         "\n    burden window:     ", format(min(epi$burden$burden_start_date)),
         " .. ", format(max(epi$burden$burden_end_date)),
         "\n    weekly admissions: ", format(min(epi$admissions$target_end_date)),
         " .. ", format(max(epi$admissions$target_end_date)),
         "\n  The sampler needs both sides non-empty.",
         call. = FALSE)
  }

  # ---- Protection tables ----
  # Both programmes emit the same (coverage, residual_ve) contract, and
  # the config forces them onto disjoint age bands, so the two tables can
  # simply be stacked. apply_scenario() then needs to know nothing about
  # birth cohorts or campaigns.
  grid  <- baseline_df %>% distinct(age_group, target_end_date, sample)
  weeks <- unique(baseline_df$target_end_date)

  # One VE draw per sample, shared across scenarios so a given sample
  # index means the same VE world in every arm.
  ve_draws <- tibble(
    sample  = sort(unique(baseline_df$sample)),
    vacc_IE = rnorm(length(unique(baseline_df$sample)),
                    mean = cfg$infant$vacc_IE_mean,
                    sd   = cfg$infant$vacc_IE_sd)
  )

  # The infant table's expensive part - the rowwise birth-window overlap -
  # does not depend on uptake, so it is built ONCE here and each scenario
  # just scales it. With 16 scenarios this is the difference between one
  # rowwise pass per country and seventeen.
  infant_base <- build_infant_base(grid,
                                   cfg$infant$vacc_start, cfg$infant$vacc_end,
                                   cfg$infant$age_bounds, cfg$infant$waning_df,
                                   ve_draws)

  # adult_bands varies per scenario: a scenario may retarget the adult
  # programme (e.g. 75+ instead of 65+) at the same coverage.
  protection_for <- function(infant_uptake, adult_coverage, adult_bands) {
    bind_rows(
      scale_infant_protection(infant_base, infant_uptake),
      build_adult_protection(weeks, adult_waning, cfg$adult$campaigns,
                             adult_coverage, cfg$adult$ve_beyond_curve,
                             adult_bands)
    )
  }

  # Coverage already embedded in the observed data - drives the
  # back-calculation for every scenario. This uses the programme-wide
  # bands, not any scenario's override: it describes what actually
  # happened, not a hypothetical targeting.
  protection_baseline <- protection_for(cfg$infant$baseline_uptake,
                                        cfg$adult$baseline_coverage,
                                        cfg$adult$eligible_age_groups)

  # ---- Scenarios ----
  sc <- cfg$scenarios_df
  scenario_results <- setNames(
    lapply(seq_len(nrow(sc)), function(i) {
      apply_scenario(baseline_df,
                     protection_for(sc$infant_uptake[i], sc$adult_coverage[i],
                                    sc$adult_age_groups[[i]]),
                     protection_baseline)
    }),
    sc$id
  )

  submission_pre <- assemble_submission(
    scenario_results = scenario_results,
    round_id         = cfg$round_id,
    anchor           = cfg$anchor
  )

  # ---- Administered doses ----
  doses_df <- build_dose_table(
    baseline_df     = baseline_df,
    births_df       = births,
    population_df   = popn,
    vacc_start      = cfg$infant$vacc_start,
    vacc_end        = cfg$infant$vacc_end,
    cfg             = cfg,
    round_id        = cfg$round_id,
    anchor          = cfg$anchor,
    output_type_ids = unique(submission_pre$output_type_id)
  )

  submission <- bind_rows(submission_pre, doses_df) %>%
    mutate(location = country_iso2)

  # Adult coverage / residual VE per scenario, kept for inspection via
  # plot_adult_protection(). These now feed the admissions above rather
  # than sitting unused.
  adult_protection <- setNames(
    lapply(seq_len(nrow(sc)), function(i) {
      build_adult_protection(weeks, adult_waning, cfg$adult$campaigns,
                             sc$adult_coverage[i], cfg$adult$ve_beyond_curve,
                             sc$adult_age_groups[[i]])
    }),
    sc$id
  )

  list(
    submission       = submission,
    baseline_df      = baseline_df,
    adult_protection = adult_protection
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
