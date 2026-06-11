# Static RSV Model — RespiCompass 2025/2026

A sampling-based model that produces RSV hospitalisation projections under
user-specified vaccination scenarios. It was developed as part of the 2025/2026
[RespiCompass](https://github.com/european-modelling-hubs/RespiCompass) modelling round on RSV interventions.

---

## How it works

1. **Load** observed weekly RSV admissions and monthly age-distribution proportions.
2. **Sample** — for each 4-week burden period, draw `n_draws` contingency tables
   of (week × age group) admissions that preserve both marginals exactly
   (`stats::r2dtable`). This captures the uncertainty in how weekly totals
   distribute across age bands.
3. **Apply scenarios** — each scenario back-calculates the no-vaccination
   counterfactual from the observed data using the baseline uptake, then
   re-applies the scenario uptake combined with a per-draw vaccine effectiveness
   (sampled from a normal distribution) and age-specific waning.
4. **Format** results into the RespiCompass submission schema and build a
   parallel administered-doses table driven by projected births.

---

## Repository structure

```
launch.R                          # entry point — run this
config/
  static_model.yaml               # all tunable parameters
data/
  epidemiological/
    RSV_weekly_counts.csv         # observed weekly RSV admissions
    RSV_monthly_prop_age.csv      # monthly age-split proportions
  population/
    country_monthly_births.csv    # monthly births
R/
  utils.R                         # round_preserve_sum
  simulate_margins.R              # fixed-margin Monte-Carlo sampler
  apply_scenario.R                # vaccination scenario logic + birth-window helper
  load_data.R                     # config and data loading functions
  format_submission.R             # RespiCompass submission formatting
  build_doses.R                   # administered-doses table
  plots.R                         # diagnostic plots
```

---

## Quick start

```r
# From the project root in R or via Rscript:
source("launch.R")
```

Required packages: `dplyr`, `tidyr`, `purrr`, `tibble`, `lubridate`,
`data.table`, `yaml`, `ggplot2`, `nanoparquet`.

---

## Configuration

All parameters live in [`config/static_model.yaml`](config/static_model.yaml).

| Parameter | Description |
|-----------|-------------|
| `vacc_IE` | Vaccine effectiveness distribution (`mean`, `sd`) |
| `waning_by_band` | Residual susceptibility after waning, per age band (1 = no protection) |
| `vaccination_start/end` | Vaccination window dates (one entry per season) |
| `baseline_uptake` | RSV vaccine uptake already reflected in the observed data |
| `round_id` | RespiCompass round identifier |
| `submission_horizon_anchor` | Anchor date for computing the `horizon` column |
| `n_draws` / `mc_seed` | Monte-Carlo sample count and RNG seed |
| `data_start` | Earliest date included in the baseline sampler |

---

## Scenarios

| Scenario | Uptake | Purpose |
|----------|--------|---------|
| `baseline` | observed | Observed data expressed in submission format |
| `no_vacc` | 0 % | Counterfactual — no RSV vaccination |
| `high_vacc` | 95 % | Counterfactual — near-universal uptake |

---

## Output

`submission` is a `data.table` with columns:

`round_id`, `scenario_id`, `target`, `pop_group`, `horizon`,
`target_end_date`, `output_type`, `output_type_id`, `value`

`target` is either `rsv_hospitalisations` or `administered_doses`.  
`pop_group` follows the pattern `<age_band>_<imm_status>` (e.g. `0-2mo_immTotal`)
plus all-ages aggregates (`total_immNo`, `total_immYes`, `total_immTotal`).

Uncomment the `write_parquet` call at the bottom of `launch.R` to persist the
output to disk.

---

## Diagnostic plots

Four plot functions are defined in [`R/plots.R`](R/plots.R) and can be called
after running `launch.R`:

| Function | Shows |
|----------|-------|
| `plot_baseline_samples(baseline_df)` | MC uncertainty bands on weekly total admissions |
| `plot_scenario_comparison(submission_pre)` | Median + 95 % CI per scenario |
| `plot_age_breakdown(submission_pre, scenario)` | Age-group breakdown for one scenario |
| `plot_dose_schedule(doses_df)` | Weekly administered doses by scenario |
