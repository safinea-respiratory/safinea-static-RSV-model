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
# The infant and adult vaccination programmes are INDEPENDENT: different
# products, eligibility rules, waning data and uncertainty models. Only
# the infant programme is applied at present; the adult waning ensemble
# is loaded and validated but not yet wired in.
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
source("R/validate.R")
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

# Fail fast: every age band named in the config must exist in the data,
# and the two programmes must cover disjoint bands. Runs before any
# modelling because the joins it protects fail silently.
validate_config(cfg, epi)

# Adult waning ensemble: one whole VE-over-time curve per sample.
# Loaded and validated now; not yet consumed by the model.
adult_waning <- load_waning_curves(cfg, n_draws = cfg$n_draws)


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


# ---- Infant vaccination scenarios ------------------------------
# run_scenario() back-calculates the no-vaccination counterfactual
# from the observed baseline (using cfg$infant$baseline_uptake) then
# re-applies the requested uptake to give the scenario admissions.
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

# no_vacc   — counterfactual: zero vaccination uptake
# high_vacc — counterfactual: 95 % uptake
scenario_A_df <- run_scenario(cfg$infant$scenarios$no_vacc)
scenario_B_df <- run_scenario(cfg$infant$scenarios$high_vacc)


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
  vacc_start      = cfg$infant$vacc_start,
  vacc_end        = cfg$infant$vacc_end,
  cfg             = cfg,
  anchor          = cfg$anchor,
  output_type_ids = unique(submission_pre$output_type_id)
)


# ---- Final output ----------------------------------------------
submission <- bind_rows(submission_pre, doses_df) %>% as.data.table()


# ---- Diagnostic plots (uncomment to view) ----------------------
# plot_baseline_samples(baseline_df)
# plot_scenario_comparison(submission_pre)
# plot_age_breakdown(submission_pre, scenario = "baseline",
#                    age_order = unlist(cfg$age_group_order))
# plot_dose_schedule(doses_df)


# ---- Persist (uncomment to write) ------------------------------
# write_parquet(
#   submission,
#   "output/respiCompass_2025_2026_results_staticModel.parquet",
#   compression = "gzip"
# )
