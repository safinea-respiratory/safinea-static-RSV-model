# ==============================================================
# launch.R
# Static RSV Model — RespiCompass 2025/2026
#
# Sampling-based (non-ODE) model that:
#   1. Loads observed weekly RSV admissions and monthly age splits.
#   2. Generates Monte-Carlo (week × age_group) admission tables that
#      preserve both marginals exactly (stats::r2dtable).
#   3. Applies per-scenario vaccination corrections (uptake × VE × waning).
#   4. Formats results into the RespiCompass submission schema and
#      builds a parallel administered-doses table.
#
# Configuration: config/static_model.yaml
# Run from the project root: Rscript launch.R
# ==============================================================

library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(lubridate)
library(data.table)
library(yaml)
library(ggplot2)
library(nanoparquet)

source("R/utils.R")
source("R/simulate_margins.R")
source("R/apply_scenario.R")
source("R/load_data.R")
source("R/format_submission.R")
source("R/build_doses.R")
source("R/plots.R")


# ---- Config ----------------------------------------------------
cfg <- load_config("config/static_model.yaml")


# ---- Data ------------------------------------------------------
epi    <- load_epidemiological_data(cfg)
births <- load_births_data(cfg)


# ---- Baseline Monte-Carlo --------------------------------------
# Draw n_draws contingency tables per 4-week period. Each draw
# respects the observed weekly totals (row margins) and the
# observed per-age-band 4-week totals (column margins) exactly.
baseline_df <- simulate_weekly_age_fixed_margins(
  weekly_df = epi$admissions %>% filter(target_end_date   >= ymd(cfg$data_start)),
  age_df    = epi$burden     %>% filter(burden_start_date >= ymd(cfg$data_start)),
  n_draws   = cfg$n_draws,
  seed      = cfg$mc_seed
)


# ---- Vaccination scenarios -------------------------------------
# run_scenario() back-calculates the no-vaccination counterfactual
# from the observed baseline (using cfg$baseline_uptake) then
# re-applies the requested uptake to give the scenario admissions.
run_scenario <- function(uptake) {
  apply_scenario(
    df                   = baseline_df,
    IE_mean              = cfg$vacc_IE$mean,
    IE_sd                = cfg$vacc_IE$sd,
    vacc_uptake          = uptake,
    vacc_start           = cfg$vacc_start,
    vacc_end             = cfg$vacc_end,
    vacc_uptake_baseline = cfg$baseline_uptake,
    waning_df            = cfg$waning_function
  )
}

# no_vacc   — counterfactual: zero vaccination uptake
# high_vacc — counterfactual: 95 % uptake
scenario_A_df <- run_scenario(cfg$scenarios$no_vacc)
scenario_B_df <- run_scenario(cfg$scenarios$high_vacc)


# ---- RespiCompass submission format ----------------------------
submission_pre <- assemble_submission(
  baseline_df   = baseline_df,
  scenario_A_df = scenario_A_df,
  scenario_B_df = scenario_B_df,
  round_id      = cfg$round_id,
  anchor        = cfg$anchor
)


# ---- Administered doses ----------------------------------------
doses_df <- build_dose_table(
  baseline_df     = baseline_df,
  births_df       = births,
  vacc_start      = cfg$vacc_start,
  vacc_end        = cfg$vacc_end,
  cfg             = cfg,
  anchor          = cfg$anchor,
  output_type_ids = unique(submission_pre$output_type_id)
)


# ---- Final output ----------------------------------------------
submission <- bind_rows(submission_pre, doses_df) %>% as.data.table()


# ---- Diagnostic plots (uncomment to view) ----------------------
# plot_baseline_samples(baseline_df)
# plot_scenario_comparison(submission_pre)
# plot_age_breakdown(submission_pre, scenario = "baseline")
# plot_dose_schedule(doses_df)


# ---- Persist (uncomment to write) ------------------------------
# write_parquet(
#   submission,
#   "output/respiCompass_2025_2026_results_staticModel.parquet",
#   compression = "gzip"
# )
