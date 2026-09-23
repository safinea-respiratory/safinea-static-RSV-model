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
# products, eligibility rules, waning data and uncertainty models. Both
# feed admissions and doses. For the 2026/27 round infant uptake is zero
# in every scenario, so only the adult programme has an effect.
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
library(arrow)
library(nanoparquet)

source("R/utils.R")
source("R/validate.R")
source("R/simulate_margins.R")
source("R/apply_scenario.R")
source("R/adult_protection.R")
source("R/scenario_engine.R")
source("R/load_data.R")
source("R/format_submission.R")
source("R/build_doses.R")
source("R/run_model.R")
source("R/submission_output.R")
source("R/scenario_impact.R")
source("R/impact_plots.R")
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
# Streaming: each country validates its own rows, reduces them to its
# impact tables, writes its own parquet part, and returns only that
# summary. Nothing large is ever held.
#
# The whole submission is ~14 GB at 28 countries, 16 scenarios and 100
# draws. Building it in memory needs more RAM than most machines have,
# and under a cluster it was also the bottleneck - serialising 500 MB per
# country back to the parent cost about what the parallelism saved.
#
# workers = 1 is serial. It caps itself against free memory either way.
PARTS <- file.path("output", "parts")

results <- run_all_countries(cfg, adult_waning, raw,
                             workers = 6,
                             out_dir = PARTS)

# Every configured country must have produced a part. Everything else -
# grid completeness, duplicates, immYes + immNo == immTotal, the all-ages
# totals - was already checked per country as it was built.
validate_submission_set(results$summary, cfg)


# ---- Diagnostic plots (uncomment to view) ----------------------
# These need the per-country objects, which streaming does not keep. To
# use them, re-run one country on its own:
#
#   ie <- run_country(cfg, "Ireland", "IE", adult_waning, raw)
#
# plot_baseline_samples(ie$baseline_df)
# plot_scenario_comparison(submission %>% filter(location == "IE"))
# plot_age_breakdown(submission %>% filter(location == "IE"),
#                    scenario  = "no_vacc",
#                    age_order = unlist(cfg$age_group_order))
# plot_dose_schedule(submission %>% filter(location == "IE",
#                                          target == "administered_doses"))
# plot_adult_protection(ie$adult_protection$adult_70)


# ---- Scenario impact summaries ---------------------------------
# Reduces the sample-level submission to the standard impact tables:
# admissions averted and doses, absolute and per 100k of total or
# eligible population. Every interval is over PAIRED samples - scenario
# minus baseline within a draw - so the pairing is not thrown away.
# Already computed per country and bound here - every group_by and join
# in compute_scenario_impact() carries `location`, so splitting by
# country and binding gives the identical tables.
impact <- results$impact
write_scenario_impact(impact, "output/3_results")

# Reference band for the relative-change plots, one per scenario AND
# season: -100 x coverage x mean VE, spanning the range of dose ages that
# season contains. With a single autumn campaign that is roughly months
# 0-11 in the first season and 12-23 in the second.
all_weeks <- sort(unique(parse_dates_strict(
  raw$admissions$target_end_date, cfg$input_date_formats$weekly_counts,
  "hospitaladmissions.csv$target_end_date")))
all_weeks <- all_weeks[all_weeks >= ymd(cfg$data_start)]
expected  <- expected_reduction(cfg, week_season_map(cfg, raw, all_weeks))

impact_figures <- build_impact_plots(impact, expected,
                                     dir = "output/3_results/figures")


# ---- Persist ---------------------------------------------------
# Concatenate the per-country parts into the single dense file the hub
# expects. Arrow writes it row group by row group, so peak memory is one
# country rather than the whole 14 GB.
write_submission_parts(results$parts, cfg)

