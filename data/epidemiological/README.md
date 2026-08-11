# Epidemiological data

Two files are required, spanning multiple RSV seasons.

**Age bands are defined by these files, not by the code.** Whatever labels
appear in `age_gp_modelling` become the model's age bands and flow through to
the submission's `pop_group` column. The config only *references* them — any
band named in `config/static_model.yaml` that is missing from the data is a
hard error at load.

**Date formats are declared, not guessed.** Each file's date format is set in
`input_date_formats` in the config and parsed strictly. This matters: R's
`as.Date("01/09/2025")` does not return `NA`, it silently returns `0001-09-20`,
which previously caused every downstream join to yield zero rows without any
error. Declared formats plus a plausible-year check make that impossible.

---

## RSV_weekly_counts.csv

Weekly aggregate RSV hospitalisation counts. One row per week.
Date format: `input_date_formats.weekly_counts`.

| Column | Type | Description |
|--------|------|-------------|
| `date_wk_floor` | date | Monday of the ISO week (the model adds 6 days to get the week-ending Sunday) |
| `season_name` | string | RSV season label, e.g. `2025/2026` |
| `case_counts` | integer | Total RSV hospitalisations for that week |

**Example**

```
date_wk_floor,season_name,case_counts
2025-09-01,2025/2026,5
2025-09-08,2025/2026,7
```

---

## RSV_monthly_prop_age.csv

Age-distribution proportions over rolling 4-week periods. One row per
(period × age group). Save with **UTF-8 BOM** encoding
(`fileEncoding = "UTF-8-BOM"` in R) — the age labels may contain
non-ASCII characters.
Date format: `input_date_formats.monthly_age`.

| Column | Type | Description |
|--------|------|-------------|
| `date_28days_floor` | date | Monday of the 4-week period (the model adds 6 days to align with the weekly data) |
| `age_gp_modelling` | string | Age band label — defines the model's bands (see below) |
| `proportion` | float [0, 1] | Fraction of RSV admissions in this age band during the period |
| `season_name` | string | RSV season label, e.g. `2025/2026` |

Proportions for a given `date_28days_floor` must sum to 1.0 across all
age bands.

### Age band labels

The labels here are authoritative. If they don't match the labels you want
in the submission, either change them in this file, or map them with the
optional `age_group_aliases` block in the config:

```yaml
age_group_aliases:
  "< 3 months": "0-2mo"
```

Leave `age_group_aliases` empty when the file already uses the desired
labels. RespiCompass currently uses:

`0-2mo, 3-5mo, 6-11mo, 1-4, 5-17, 18-59, 60-64, 65-69, 70-74, 75-79, 80+`

Finer adult bands matter for the adult vaccination programme: a 75+ campaign
can only be represented correctly if `75-79` and `80+` exist as separate
bands. Lumping them into a single `65+` band and scaling by a population
fraction understates the effect, because RSV hospitalisation risk rises
steeply with age.

**Example**

```
date_28days_floor,age_gp_modelling,proportion,season_name
01/09/2025,< 3 months,0.1,2025/2026
01/09/2025,3-5 months,0.3,2025/2026
01/09/2025,6-11 months,0.2,2025/2026
01/09/2025,1-4 years,0.166667,2025/2026
01/09/2025,5-64 years,0.133333,2025/2026
01/09/2025,65+ years,0.1,2025/2026
```
