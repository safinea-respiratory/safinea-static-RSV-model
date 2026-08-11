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

## Two independent vaccination programmes

The infant and adult programmes are different products and are modelled
separately. They share no parameters, and the config requires them to cover
disjoint age bands.

| | Infant (maternal / birth-dose) | Adult |
|---|---|---|
| Eligibility | Birth cohort ∩ vaccination window | Calendar campaign coverage |
| Waning indexed by | **Age band** — age equals time since dose | **Months since dose** (0–36) |
| Waning source | `waning_by_band` in the config | `data/vaccine/waning_curves.csv` |
| Uncertainty | Parametric draw on `vacc_IE` | Empirical — one whole curve per sample |
| Coverage accrual | Per season, birth-driven | One-off and cumulative across seasons |
| Status | Applied to admissions | **Coverage computed, not yet applied** |

The adult programme's coverage and residual VE are produced by convolving the
campaign uptake curve with the VE ensemble:

```
coverage(t)    = Σ  δ(w)              for vaccination weeks w ≤ t
protection(t)  = Σ  δ(w) · VE(t − w)
residual_ve(t) = protection(t) / coverage(t)
```

Because expected admissions are *linear* in VE, this coverage-weighted mean is
exact rather than an approximation. `protection(t)` is precisely the
`coverage × residual_VE` term the scenario arithmetic needs, which is why
wiring it in (step 4) leaves the surrounding algebra untouched.

---

## Age bands come from the data

The model does not define age bands. Whatever labels appear in
`RSV_monthly_prop_age.csv` become the bands, and flow through to the
submission's `pop_group`. The config only *references* them — any band named
in the config but absent from the data is a hard error at load, as is any
overlap between the two programmes. This matters because the joins involved
otherwise fail silently, yielding `NA` values rather than an error.

---

## Repository structure

```
launch.R                          # entry point — run this
config/
  static_model.yaml               # all tunable parameters
data/
  epidemiological/
    RSV_weekly_counts.csv         # observed weekly RSV admissions
    RSV_monthly_prop_age.csv      # monthly age-split proportions (defines age bands)
  population/
    country_monthly_births.csv    # monthly births
  vaccine/
    waning_curves.csv             # adult VE ensemble, 500 curves (from RespiCompass)
R/
  utils.R                         # round_preserve_sum
  validate.R                      # fail-fast config/data consistency checks
  simulate_margins.R              # fixed-margin Monte-Carlo sampler
  apply_scenario.R                # INFANT programme: birth-window + scenario logic
  adult_protection.R              # ADULT programme: campaign x waning convolution
  load_data.R                     # config, data and waning-curve loaders
  format_submission.R             # RespiCompass submission formatting
  build_doses.R                   # administered-doses table
  plots.R                         # diagnostic plots
```

`data/vaccine/waning_curves.csv` is vendored from RespiCompass — see
[data/vaccine/README.md](data/vaccine/README.md) for provenance. It is the only
external dataset in the repo; nothing is fetched over the network at runtime.

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

**Top level**

| Parameter | Description |
|-----------|-------------|
| `input_date_formats` | strptime format for each input CSV's date column. Parsed strictly — see note below |
| `age_group_aliases` | Optional rename of source age labels. Leave empty if the data already uses the desired labels |
| `age_group_order` | Display order for plots. Optional; the fallback is alphabetical, which orders age bands wrongly |
| `round_id` | RespiCompass round identifier |
| `submission_horizon_anchor` | Anchor date for computing the `horizon` column |
| `n_draws` / `mc_seed` | Monte-Carlo sample count and RNG seed |
| `data_start` | Earliest date included in the baseline sampler |

**`infant_vaccination`**

| Parameter | Description |
|-----------|-------------|
| `vacc_IE` | Vaccine effectiveness distribution (`mean`, `sd`) |
| `waning_by_band` | Residual protection per age band (1 = full initial protection, 0 = none left). Bands not listed default to 0 |
| `age_bounds` | Age span of each band in months. List every band the programme has reached, **including fully-waned ones** — these bounds also drive the vaccinated/unvaccinated split |
| `windows` | Vaccination window dates (one entry per season) |
| `baseline_uptake` | Uptake already reflected in the observed data |
| `scenarios` | Per-scenario uptake |

**`adult_vaccination`**

| Parameter | Description |
|-----------|-------------|
| `eligible_age_groups` | Bands the adult programme covers. Change this to retarget the programme |
| `waning_curves` | Path to the VE ensemble |
| `ve_target` | Which VE column to use (`VE_sev`) |
| `ve_beyond_curve` | Protection past month 36 — `zero` or `hold_last` |
| `campaigns` | Campaign windows with `share` (portion of total coverage, must sum to 1) and `profile` (`uniform`) |
| `baseline_coverage` | Cumulative coverage already reflected in the observed data |
| `scenarios` | Total cumulative coverage of the eligible population, per scenario |

> **Dates are declared, never guessed.** `as.Date("01/09/2025")` does not
> return `NA` in R — it silently returns `0001-09-20`. Declaring the format in
> `input_date_formats`, combined with a plausible-year check, turns that class
> of silent corruption into a load-time error.

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
| `plot_adult_protection(adult_prot_high_vacc)` | Adult coverage, residual VE and effective protection over time |
