# Epidemiological data

Two files are required, spanning multiple RSV seasons.

---

## RSV_weekly_counts.csv

Weekly aggregate RSV hospitalisation counts. One row per week.

| Column | Type | Description |
|--------|------|-------------|
| `date_wk_floor` | `YYYY-MM-DD` | Monday of the ISO week  |
| `season_name` | string | RSV season label, e.g. `2018/2019` |
| `case_counts` | integer | Total RSV hospitalisations for that week |

**Example**

```
date_wk_floor,season_name,case_counts
2018-05-21,2018/2019,5
2018-05-28,2018/2019,3
```

---

## RSV_monthly_prop_age.csv

Age-distribution proportions over rolling 4-week periods. One row per
(period × age group). The file must be saved with **UTF-8 BOM** encoding
(`fileEncoding = "UTF-8-BOM"` in R) because the age labels contain
non-ASCII characters.

| Column | Type | Description |
|--------|------|-------------|
| `date_28days_floor` | `YYYY-MM-DD` | Monday of the 4-week period (the model adds 6 days to align with the weekly data) |
| `age_gp_modelling` | string | Age band — must be one of the six labels below |
| `proportion` | float [0, 1] | Fraction of RSV admissions in this age band during the period |
| `season_name` | string | RSV season label, e.g. `2018/2019` |

**Expected age band labels**

| `age_gp_modelling` | Mapped to |
|--------------------|-----------|
| `< 3 months` | `0-2mo` |
| `3-5 months` | `3-5mo` |
| `6-11 months` | `6-11mo` |
| `1-4 years` | `1-4y` |
| `5-64 years` | `5-64y` |
| `65+ years` | `65+y` |

Proportions for a given `date_28days_floor` must sum to 1.0 across all six
age bands.

**Example**

```
date_28days_floor,age_gp_modelling,proportion,season_name
2018-05-21,< 3 months,0.461538462,2018/2019
2018-05-21,3-5 months,0.153846154,2018/2019
2018-05-21,6-11 months,0.153846154,2018/2019
2018-05-21,1-4 years,0.076923077,2018/2019
2018-05-21,5-64 years,0.076923077,2018/2019
2018-05-21,65+ years,0.076923077,2018/2019
```
