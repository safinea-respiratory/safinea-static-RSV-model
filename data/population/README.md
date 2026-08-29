# Population data

Both files come from RespiCompass auxiliary-data and are keyed by
**ISO2 country code** (`IE`, `AT`, …) — *not* the full country names used by
the epidemiological target-data. This is why each entry under `countries` in
the config carries both a `name` and an `iso2`.

### Provenance

| | |
|---|---|
| Source | [RespiCompass `auxiliary-data`](https://github.com/european-modelling-hubs/RespiCompass/tree/main/auxiliary-data) |
| Commit | `142d08c97df01fc42110757c5e1bbf716a37bea4` |
| Committed | 2026-07-28 |
| Retrieved | 2026-08-11 |

---

## births_by_month.csv

Monthly live births, already mapped onto the **2026-09 → 2027-08** scenario
period so every model uses a consistent birth cohort. 30 countries × 12 months.

| Column | Type | Description |
|--------|------|-------------|
| `country` | string | **ISO2** code |
| `date` | date | First day of the month |
| `month` / `year` | string / integer | Convenience columns — unused by the model |
| `births` | numeric | Live births that month |

Derived from Eurostat `DEMO_FMONTH`, using the most recent available year
(2025, falling back to 2024/2023).

Because RespiCompass supplies this already aligned to the scenario period, the
model no longer projects a historical monthly pattern forward — that logic and
its `births_data_cutoff` / `births_project_years` settings have been removed.

**Example**

```
country,date,month,year,births
IE,2026-09-01,September,2026,4555.0
IE,2026-10-01,October,2026,4698.0
```

---

## population_estimates.csv

Population by country × age band, using the **same 11 bands** as the burden
data. 30 countries × 11 bands.

| Column | Type | Description |
|--------|------|-------------|
| `country` | string | **ISO2** code |
| `age_group` | string | Age band — matches the epidemiological bands |
| `population` | integer | Number of individuals |

**Not yet consumed by the model.** It is loaded and validated, and is the
denominator the adult administered-doses track will need: adult doses are
`eligible population × weekly coverage increment`, whereas infant doses come
from births. For Ireland the eligible 65+ population totals 853,555 across
`65-69`, `70-74`, `75-79` and `80+`.

**Example**

```
country,age_group,population
IE,0-2mo,13961
IE,65-69,255030
IE,80+,206177
```
