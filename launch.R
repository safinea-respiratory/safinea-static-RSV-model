# ==============================================================
# launch.R
# Static RSV Model — RespiCompass 2026/2027
#
# Sampling-based (non-ODE) model that:
#   1. Loads observed weekly RSV admissions and the age-stratified
#      seasonal burden, per country.
#   2. Generates Monte-Carlo (week × age_group) admission tables that
#      preserve both marginals exactly (stats::r2dtable).
#   3. Applies per-scenario vaccination corrections (uptake × VE × waning).
#   4. Formats results into the RespiCompass submission schema and
#      builds a parallel administered-doses table.
#
# Runs once per country listed in the config and binds the results,
# distinguished by the `location` column.
#
# The infant and adult vaccination programmes are INDEPENDENT: different
# products, eligibility rules, waning data and uncertainty models. Only
# the infant programme is applied at present; adult coverage and
# residual VE are computed but not yet applied to admissions.
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
source("R/adult_protection.R")
source("R/load_data.R")
source("R/format_submission.R")
source("R/build_doses.R")
source("R/run_model.R")
source("R/submission_output.R")
source("R/plots.R")


# ---- Config ----------------------------------------------------
cfg <- load_config("config/static_model.yaml")

# Every input CSV is read once here and filtered per country downstream.
raw <- load_raw_inputs(cfg)

# Fail fast: every configured country must resolve in every input file,
# and the adult settings must be coherent, before any modelling starts.
validate_global_config(cfg, raw)

# Adult waning ensemble: one whole VE-over-time curve per sample.
# Loaded once and shared across countries.
adult_waning <- load_waning_curves(cfg, n_draws = cfg$n_draws)


# ---- Run every configured country ------------------------------
results    <- run_all_countries(cfg, adult_waning, raw)
submission <- results$submission


# ---- Check the output before anything consumes it ---------------
# Internal consistency: identifiers match the config, the grid is
# complete, no duplicates, immYes + immNo == immTotal, and the all-ages
# totals equal the sum over bands. Errors list every problem at once.
validate_submission(submission, cfg)


# ---- Diagnostic plots (uncomment to view) ----------------------
# Per-country objects live in results$by_country[["IE"]]
# ie <- results$by_country[["IE"]]
#
# plot_baseline_samples(ie$baseline_df)
# plot_scenario_comparison(submission %>% filter(location == "IE"))
# plot_age_breakdown(submission %>% filter(location == "IE"),
#                    scenario  = "no_vacc",
#                    age_order = unlist(cfg$age_group_order))
# plot_dose_schedule(submission %>% filter(location == "IE",
#                                          target == "administered_doses"))
# plot_adult_protection(ie$adult_protection$adult_70)


# ---- Persist ---------------------------------------------------
# Writes output/<round_id>_staticModel.parquet, creating the directory if
# needed. Pass `path` to override. output/ is gitignored.
write_submission(submission, cfg)
