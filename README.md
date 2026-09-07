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
   parallel administered-doses table covering all three programmes.

### Administered doses

All three programmes contribute, each against its own denominator:

| | Denominator | Weekly dose count |
|---|---|---|
| **Infant** | Monthly births | births × scenario `infant_uptake`, masked to the vaccination windows |
| **Catch-up** | Births in the cohort's birth window | cohort size × that week's coverage increment, split across the bands the cohort occupies **when dosed** |
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

## Three vaccination programmes

| | Infant (birth dose) | Adult | Infant catch-up |
|---|---|---|---|
| Vaccine | infant product | **different product** | infant product |
| Eligibility | birth cohort ∩ vaccination window | **fixed age bands** | **age range on the day of the dose** |
| Waning indexed by | **age band** — age equals time since dose | months since dose | months since dose |
| Waning source | `waning_by_band` in the config | `data/vaccine/waning_curves.csv` | `waning_by_band`, re-indexed onto months since dose |
| Uncertainty | parametric draw on `vacc_IE` | empirical — one whole curve per sample | parametric draw on `vacc_IE` |
| Coverage accrual | per season, birth-driven | one-off and cumulative | one-off, single campaign |
| Status | wired, dormant this round (uptake 0) | **applied** | **applied** |

The infant and adult programmes are different products and must cover disjoint
age bands — the config enforces it. The catch-up is the *same product* as the
birth dose given on a calendar campaign, so it shares that programme's `vacc_IE`
and `waning_by_band` and deliberately shares its age bands too.

### Why the catch-up needs its own eligibility rule

A static band list cannot express it. The adult programme assumes the vaccinated
group *is* an age band and stays it — true for a 67-year-old, who is still in
`65-69` six months later. A catch-up defines the group once, on the day of the
dose, and that group then **ages across bands** while the band refills with
children who were never eligible.

For an under-6-months campaign, the share of each band that is actually
vaccinated moves like this:

| Weeks after campaign | `0-2mo` | `3-5mo` | `6-11mo` |
|---|---|---|---|
| 0 | 100 % | 100 % | 0 % |
| 13 | 0 % | 100 % | 51 % |
| 26 | 0 % | 0 % | 100 % |

So `6-11mo` — a band you would never have listed — ends up holding all of the
protection, and `0-2mo` holds none. `eligibility_birth_cohort()` in
[`R/campaign_protection.R`](R/campaign_protection.R) computes those fractions
from the overlap of two birth-date windows.

Note that doses and protection are counted in *different* bands: doses where the
children were when injected, protection where they are when it acts.

### The shared convolution

Both campaign programmes produce coverage and residual VE by convolving the
campaign uptake curve with their VE curve, weighted by eligibility:

```
coverage(b,t)    = Σ  δ(w) · w_elig(b,t,w)     for vaccination weeks w ≤ t
protection(b,t)  = Σ  δ(w) · w_elig(b,t,w) · VE(t − w)
residual_ve(b,t) = protection(b,t) / coverage(b,t)
```

`w_elig` is the eligibility weight: identically 1 for the adult programme, so
both sums collapse to the plain campaign convolution and its numbers are
unchanged (verified bit-for-bit against the previous implementation).

Because expected admissions are *linear* in VE, this coverage-weighted mean is
exact rather than an approximation — and for the same reason the eligibility
weight can be folded into the same sum.

All three programmes emit the **same two-column contract** — `coverage` and
`residual_ve` per (age band, week, sample). `combine_protection()` merges them,
summing within a band as programmes reaching disjoint people, and errors if the
combined coverage exceeds 1. `apply_scenario()` then knows nothing about birth
cohorts or campaigns; it just runs:

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
  campaign_protection.R           # ADULT + CATCH-UP: campaign x waning convolution
  load_data.R                     # config, data and waning-curve loaders
  format_submission.R             # RespiCompass submission formatting
  build_doses.R                   # administered-doses table
  run_model.R                     # per-country pipeline + multi-country driver
  submission_output.R             # output validation + parquet writing
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
| `scenarios` | List of `{id, infant_uptake, adult_coverage, catchup_coverage}`. Drives the whole scenario set — `id` becomes `scenario_id` in the submission |
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
| `windows` | Vaccination window dates (one entry per season). Must not overlap — checked at load, because uptake is set per window |
| `baseline_uptake` | Uptake already reflected in the observed data (0 this round). Either one value for all windows, or **one value per window** |

Infant uptake — both `baseline_uptake` and each scenario's `infant_uptake` — may
be given per vaccination window, so a programme that ramped up over successive
seasons is expressible as it actually happened:

```yaml
baseline_uptake: [0.45, 0.83]     # first season, second season
```

Coverage is then

```
coverage(band, t) = Σ  uptake[s] · overlap[s]  /  width of the band's birth window
```

summed over seasons `s`, rather than one rate times the total overlap. A cohort
straddling a programme change therefore carries a **blend** of both years,
weighted by how much of it falls in each — a `1-4` band spanning 180 days of a
20 % season and 89 days of a 90 % season comes out at 43.2 %, not at either rate.

**`adult_vaccination`**

| Parameter | Description |
|-----------|-------------|
| `eligible_age_groups` | Bands the adult programme covers. Change this to retarget the programme |
| `waning_curves` | Path to the VE ensemble |
| `ve_target` | Which VE column to use (`VE_sev`) |
| `ve_beyond_curve` | Protection past month 36 — `zero` or `hold_last` |
| `campaigns` | Campaign windows with `share` (portion of total coverage, must sum to 1) and `profile` (`uniform`) |
| `baseline_coverage` | Adult coverage already reflected in the observed data (0 this round) |

**`catchup_vaccination`**

No VE or waning settings: it is the same product as the birth dose and takes
`vacc_IE` and `waning_by_band` from `infant_vaccination`.

| Parameter | Description |
|-----------|-------------|
| `age_at_campaign_months` | `[min, max)` age in months eligible **on the day of the dose**. `[0, 6]` is "every child under six months old that day" |
| `age_bounds` | Every band the cohort can occupy over the horizon — not just the ones it starts in. Checked at load: bounds that run out before the last modelled week are a hard error, because the cohort would silently lose its protection part-way through |
| `campaigns` | Campaign windows, same shape as the adult programme's |
| `baseline_coverage` | Catch-up coverage already reflected in the observed data (0 — no such campaign has run) |

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

### Validation and writing

`launch.R` validates the submission before writing it. The checks are all
internal — derived from the config and the data, with no external schema:

| Check | Catches |
|---|---|
| Required columns present | Structural breakage |
| No `NA`, `Inf` or negative values | Arithmetic gone wrong upstream |
| `round_id`, `scenario_id`, `location` match the config | Renamed or stray identifiers |
| `horizon == (target_end_date − anchor)/7` | Anchor or date drift |
| No duplicate rows on the identifying key | Join fan-out — the bug class that inflated adult admissions 4× |
| `immYes + immNo == immTotal` | Broken stratification |
| `total_*` equals the sum over age bands | Aggregation errors |
| Grid completeness per target | Dropped rows |
| `administered_doses` uses only `pop_group = "undefined"` | Target/pop_group mismatch |

All failures are collected and reported **together**, so one run tells you
everything that is wrong rather than stopping at the first problem.

The file is then written to `output/<round_id>_staticModel.parquet` (~0.9 MB
for one country), creating the directory if needed. `output/` is gitignored.
Pass `path` to `write_submission()` to override.

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
