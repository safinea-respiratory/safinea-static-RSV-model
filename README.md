# Static RSV Model — RespiCompass 2026/2027

A sampling-based model that produces RSV hospitalisation projections under
user-specified vaccination scenarios. It was developed as part of the 2026/2027
[RespiCompass](https://github.com/european-modelling-hubs/RespiCompass) modelling round on RSV interventions.

---

## How it works

1. **Load** observed weekly RSV admissions and the age-stratified seasonal burden
   (both from RespiCompass target-data, 28 EU/EEA countries).
2. **Sample** — for each burden window, draw `n_draws` contingency tables
   of (week × age group) admissions that preserve both marginals exactly
   (`stats::r2dtable`). This captures the uncertainty in how weekly totals
   distribute across age bands.
3. **Apply scenarios** — each scenario back-calculates the no-vaccination
   counterfactual from the observed data using the baseline uptake, then
   re-applies the scenario uptake combined with a per-draw vaccine effectiveness
   (sampled from a normal distribution) and age-specific waning.
4. **Format** results into the RespiCompass submission schema and build a
   parallel administered-doses table covering both programmes.

### Administered doses

Both programmes contribute, each against its own denominator:

| | Denominator | Weekly dose count |
|---|---|---|
| **Infant** | Monthly births | births × scenario `infant_uptake`, masked to the vaccination windows |
| **Adult** | Population of `eligible_age_groups` | population × that week's **new** coverage increment |

The adult figure works because coverage is one-off and cumulative — the weekly
increment *is* the number of people newly vaccinated, which is what a dose
count means. Over a whole campaign the doses sum to exactly
`population × adult_coverage`.

The infant track spreads each month's births evenly across its days, masks them
to the vaccination windows *per day* (so a week straddling a boundary is
credited only for the days inside it), then re-aggregates to ISO weeks. Births
falling in a window but outside the modelled weeks are reported with a warning
rather than silently dropped.

For Ireland 2026/27, `adult_70` gives 597,489 doses over the 13 campaign weeks —
70 % of the 853,555 people aged 65+ — averting 701 admissions, or roughly **853
doses per admission averted**.

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
| Status | Wired, dormant this round (uptake 0) | **Applied to admissions** |

The adult programme's coverage and residual VE are produced by convolving the
campaign uptake curve with the VE ensemble:

```
coverage(t)    = Σ  δ(w)              for vaccination weeks w ≤ t
protection(t)  = Σ  δ(w) · VE(t − w)
residual_ve(t) = protection(t) / coverage(t)
```

Because expected admissions are *linear* in VE, this coverage-weighted mean is
exact rather than an approximation.

Both programmes emit the **same two-column contract** — `coverage` and
`residual_ve` per (age band, week, sample) — and because the config forces them
onto disjoint bands, the two tables are simply stacked. `apply_scenario()` then
knows nothing about birth cohorts or campaigns; it just runs:

```
no_vax = observed / (1 − coverage_baseline × residual_ve_baseline)
yes    = (1 − residual_ve) × coverage       × no_vax
no     =                     (1 − coverage) × no_vax
```

Bands covered by neither programme are absent from both tables, default to zero
coverage, and pass through at their observed values.

---

## Multiple countries

The model runs once per country listed in the config and binds the results,
distinguished by the `location` column:

```yaml
countries:
  - {name: "Ireland",  iso2: "IE"}
  - {name: "Austria",  iso2: "AT"}
```

Each country is fully independent — its own data, its own sampler draws, its own
submission rows. Roughly 19 s per country.

**Two identifiers are needed** because RespiCompass is not internally consistent
about country naming:

| Identifier | Used by |
|---|---|
| `name` (`Ireland`) | target-data — `hospitaladmissions.csv`, `hospitalburden_agegroups.csv` |
| `iso2` (`IE`) | auxiliary-data — births, population — **and** the submission's `location` column |

Every configured country is checked against every input file *before* any
modelling starts, so a typo fails immediately rather than partway through a long
run. Unknown values error with the available list.

`mc_seed` is shared across countries, so the vaccine-effectiveness draws are
common to all of them — defensible, since it is the same vaccine.

---

## Age bands come from the data

The model does not define age bands. Whatever labels appear in
`hospitalburden_agegroups.csv` become the bands, and flow through to the
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
    hospitaladmissions.csv        # weekly RSV admissions, all ages
    hospitalburden_agegroups.csv  # seasonal age-stratified burden (defines age bands)
  population/
    births_by_month.csv           # monthly births, mapped to the scenario period
    population_estimates.csv      # population by age band — adult dose denominator
  vaccine/
    waning_curves.csv             # adult VE ensemble, 500 curves
R/
  utils.R                         # round_preserve_sum
  validate.R                      # fail-fast config/data consistency checks
  simulate_margins.R              # fixed-margin Monte-Carlo sampler
  apply_scenario.R                # INFANT programme: birth-window + scenario logic
  adult_protection.R              # ADULT programme: campaign x waning convolution
  load_data.R                     # config, data and waning-curve loaders
  format_submission.R             # RespiCompass submission formatting
  build_doses.R                   # administered-doses table
  run_model.R                     # per-country pipeline + multi-country driver
  plots.R                         # diagnostic plots
```

The epidemiological and vaccine datasets are all vendored from RespiCompass —
see [data/epidemiological/README.md](data/epidemiological/README.md) and
[data/vaccine/README.md](data/vaccine/README.md) for provenance and commit
hashes. Nothing is fetched over the network at runtime.

> **Time resolution caveat.** `hospitalburden_agegroups.csv` gives one total per
> age band for the *whole season*, so the sampler constrains the age mix across
> the season rather than within it. Weekly totals stay exact in every draw, but
> the per-week age split carries a median coefficient of variation of ~0.38 on
> the Ireland 2026/27 data. See the data README for detail.

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
| `countries` | List of `{name, iso2}` to model. Unknown values error with the available list |
| `input_date_formats` | strptime format for each input CSV's date column. Parsed strictly — see note below |
| `births_file` / `population_file` | Paths to the auxiliary demographic data |
| `age_group_aliases` | Optional rename of source age labels. Leave empty if the data already uses the desired labels |
| `age_group_order` | Display order for plots. Optional; the fallback is alphabetical, which orders age bands wrongly |
| `scenarios` | List of `{id, infant_uptake, adult_coverage}`. Drives the whole scenario set — `id` becomes `scenario_id` in the submission |
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
| `baseline_uptake` | Uptake already reflected in the observed data (0 this round) |

**`adult_vaccination`**

| Parameter | Description |
|-----------|-------------|
| `eligible_age_groups` | Bands the adult programme covers. Change this to retarget the programme |
| `waning_curves` | Path to the VE ensemble |
| `ve_target` | Which VE column to use (`VE_sev`) |
| `ve_beyond_curve` | Protection past month 36 — `zero` or `hold_last` |
| `campaigns` | Campaign windows with `share` (portion of total coverage, must sum to 1) and `profile` (`uniform`) |
| `baseline_coverage` | Adult coverage already reflected in the observed data (0 this round) |

> **Dates are declared, never guessed.** `as.Date("01/09/2025")` does not
> return `NA` in R — it silently returns `0001-09-20`. Declaring the format in
> `input_date_formats`, combined with a plausible-year check, turns that class
> of silent corruption into a load-time error.

---

## Scenarios

Scenarios are defined entirely by the `scenarios` block in the config — any
number, with each programme's uptake set independently. The 2026/27 round
varies adult coverage of 65+ against a no-vaccination reference:

| `scenario_id` | Infant uptake | Adult coverage | Purpose |
|----------|---------------|----------------|---------|
| `no_vacc` | 0 % | 0 % | Reference — no vaccination. Reproduces the observed data exactly |
| `adult_20` | 0 % | 20 % | Low adult uptake |
| `adult_70` | 0 % | 70 % | High adult uptake |

Infant uptake is zero throughout: the infant programme is out of scope this
round, so infant age bands pass through at their observed values. Its machinery
is retained — give a scenario a non-zero `infant_uptake` and set
`infant_vaccination.baseline_uptake` to the real coverage to re-enable it.

**Illustrative effect** (Ireland, draw 1). The four adult bands carry 1,380 of
5,116 admissions:

| Scenario | Adult-band admissions | Averted |
|---|---|---|
| `no_vacc` | 1,380 | — |
| `adult_20` | 1,180 | 200 (14.5 %) |
| `adult_70` | 679 | 701 (50.8 %) |

The effect is exactly linear in coverage — 0.7/0.2 = 3.5, and 701/200 = 3.5 —
because expected admissions are linear in VE, so the coverage-weighted mean is
exact rather than approximate.

---

## Output

`submission` is a `data.table` with columns:

`round_id`, `scenario_id`, `target`, `location`, `pop_group`, `horizon`,
`target_end_date`, `output_type`, `output_type_id`, `value`

`location` is the ISO2 country code. Per-country intermediates (baseline draws,
adult protection) are kept in `results$by_country[["IE"]]`.

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
| `plot_adult_protection(ie$adult_protection$adult_70)` | Adult coverage, residual VE and effective protection over time |
