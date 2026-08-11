# Vaccine data

## waning_curves.csv

Vaccine-effectiveness-over-time ensemble for the **adult** RSV vaccination
programme.

### Provenance

| | |
|---|---|
| Source | [RespiCompass `auxiliary-data/waning-immunity`](https://github.com/european-modelling-hubs/RespiCompass/tree/main/auxiliary-data/waning-immunity) |
| Commit | `142d08c97df01fc42110757c5e1bbf716a37bea4` |
| Committed | 2026-07-28 |
| Retrieved | 2026-08-11 |

Vendored into this repository deliberately: the model reads no data over the
network, so runs stay reproducible and offline. Re-download and update the
commit hash above if RespiCompass revises the curves.

### Structure

500 reps × 37 months = 18,500 rows.

| Column | Type | Description |
|--------|------|-------------|
| `rep` | integer 1–500 | Uncertainty sample. One `rep` is one complete, plausible VE-over-time curve |
| `month` | integer 0–36 | Months elapsed since vaccination |
| `VE_inf` | float [0, 1] | Protection against infection |
| `VE_sev` | float [0, 1] | Protection against severe disease / hospitalisation |

### How the model uses it

**One whole curve per Monte-Carlo sample.** Sample *d* takes `rep = d`.

This matters. Drawing VE independently at each month would produce
non-monotonic curves and would artificially narrow the aggregate uncertainty.
Taking an entire realisation preserves the correlation across months by
construction, and requires no distributional assumption at all.

Reps are independent realisations, so the model uses the first `n_draws` of
them. `n_draws` must not exceed the ensemble size (500) — exceeding it is an
error rather than a silent recycling of curves.

`VE_sev` is the column used (`adult_vaccination.ve_target`). It is
unconditional protection against RSV hospitalisation, which is exactly the
multiplier this model needs for its hospitalisation target. `VE_inf` is loaded
but unused: in a static model with no transmission term it has no effect, and
in particular the indirect benefit of adult vaccination in reducing
transmission to infants **cannot** be represented here.

Behaviour past month 36 is set by `adult_vaccination.ve_beyond_curve`
(`zero` or `hold_last`).

### Not the infant programme

The infant (maternal / birth-dose) programme is completely independent — a
different product with its own eligibility rules, its own uncertainty model
(a parametric draw on `vacc_IE`), and its own waning representation
(`infant_vaccination.waning_by_band`, indexed by *age band* rather than by
months since dose, because for that product age and time-since-vaccination
are the same thing). The two share no parameters, and the config enforces
that they cover disjoint age bands.
