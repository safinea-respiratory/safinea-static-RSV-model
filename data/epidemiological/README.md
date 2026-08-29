# Epidemiological data

Both files come from RespiCompass target-data and cover **28 EU/EEA countries**.
`country` in the config selects which one is modelled; an unknown country is a
load-time error listing the available options.

### Provenance

| | |
|---|---|
| Source | [RespiCompass `target-data`](https://github.com/european-modelling-hubs/RespiCompass/tree/main/target-data) |
| Commit | `3b4f627cac6b87e5533e4dfb6d42ac2fe30c3f35` |
| Committed | 2026-07-03 |
| Retrieved | 2026-08-11 |
| Season | 2026/27 |

Vendored deliberately — the model reads nothing over the network, so runs stay
reproducible and offline.

**Age bands are defined by these files, not by the code.** Whatever labels
appear in `age_group` become the model's bands and flow through to the
submission's `pop_group`. The config only *references* them; any band named in
the config but missing from the data is a hard error at load.

**Date formats are declared, not guessed** (`input_date_formats` in the config).
R's `as.Date("01/09/2025")` does not return `NA` — it silently returns
`0001-09-20`. Declared formats plus a plausible-year check make that class of
silent corruption impossible.

---

## hospitaladmissions.csv

Weekly aggregate RSV hospitalisations, all ages. One row per (country × week).

| Column | Type | Description |
|--------|------|-------------|
| `country` | string | Country name — filtered by `country` in the config |
| `age_group` | string | Always `total` — unused by the model |
| `target_end_date` | date | Week-ending date (Sunday). Used directly, no offset applied |
| `year_week` | string | ISO year-week, e.g. `2026-W36` — unused |
| `week` / `year` | integer | ISO week and year — unused |
| `weekly_rsv_hospitalisations` | integer | Total RSV hospitalisations that week |

**⚠ ISO week 53 is absent.** The series runs 2026-09-06 → 2027-08-29, which spans
52 week-endings, but the file supplies **51** — the week ending **2027-01-03** is
missing for all 28 countries. Consecutive rows are therefore 7 days apart except
for one 14-day gap.

This is consistent across countries, so it looks like a deliberate 52-week season
convention rather than corrupt data. Two consequences: `horizon` skips one integer,
and any births falling in that week cannot be credited to the dose table.
`build_dose_table()` warns with the exact number rather than dropping them
silently (1,007 births for Ireland).

**Example**

```
country,age_group,target_end_date,year_week,week,year,weekly_rsv_hospitalisations
Ireland,total,2026-09-06,2026-W36,36,2026,4
Ireland,total,2026-09-13,2026-W37,37,2026,6
```

---

## hospitalburden_agegroups.csv

Age-stratified RSV hospitalisation totals over the season. One row per
(country × age band). These are **absolute counts, not proportions** — the model
uses them directly.

| Column | Type | Description |
|--------|------|-------------|
| `country` | string | Country name — filtered by `country` in the config |
| `age_group` | string | Age band label — defines the model's bands |
| `start_date` | date | Start of the burden window |
| `end_date` | date | End of the burden window |
| `total_rsv_hospitalisations` | integer | Total RSV hospitalisations in that band over the window |

### Age bands

`0-2mo, 3-5mo, 6-11mo, 1-4, 5-17, 18-59, 60-64, 65-69, 70-74, 75-79, 80+`

The finer adult bands are what make an age-targeted adult campaign expressible.
A 75+ programme is `["75-79", "80+"]` in
`adult_vaccination.eligible_age_groups` — no population-fraction approximation,
which would otherwise bias the estimate because RSV hospitalisation risk rises
steeply with age.

### ⚠ Time resolution

This file gives **one total per band for the entire season**. The fixed-margin
sampler therefore draws a single (weeks × ages) table per season rather than one
per 4-week period, so the age mix is constrained only *across* the season, not
within it.

That is an honest reflection of what this data constrains — but the consequence
is real and measurable. On the Ireland 2026/27 data the per-week age split has a
median coefficient of variation of **0.38** across draws, with the 90th
percentile at **0.80**.

Weekly *totals* remain exact in every draw (they are a fixed margin), so
all-ages outputs are unaffected. It is the age-specific rows that carry this
uncertainty. Period-level age data, if available nationally, would tighten it
considerably.

**Example**

```
country,age_group,start_date,end_date,total_rsv_hospitalisations
Ireland,0-2mo,2026-08-31,2027-08-29,1006
Ireland,3-5mo,2026-08-31,2027-08-29,617
Ireland,80+,2026-08-31,2027-08-29,717
```

Marginals must reconcile: for Ireland 2026/27 both files sum to 5,116. When they
disagree the sampler scales one side (`reconcile` in
`simulate_weekly_age_fixed_margins`).
