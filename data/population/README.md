# Population data

## country_monthly_births.csv

Monthly live birth counts for the modelled country, sourced from
[Eurostat DEMO_FMONTH](https://ec.europa.eu/eurostat/databrowser/view/DEMO_FMONTH).
One row per (year × month). The file must be saved with **UTF-8 BOM** encoding.

The model uses this file to estimate the weekly size of the infant birth cohort
eligible for RSV immunisation. Real data is currently available up to the date
set in `births_data_cutoff` in `config/static_model.yaml`; later years are
projected by repeating the 2023 monthly pattern.

The model reads four columns and ignores the rest:

| Column | Type | Description |
|--------|------|-------------|
| `month` | string | Full English month name (e.g. `January`, `February`) |
| `geo` | string | Country name — carried through as the `country` identifier |
| `TIME_PERIOD` | integer | Calendar year (e.g. `2023`) |
| `OBS_VALUE` | integer | Number of live births in that month and year |

**Example** (Eurostat export format — extra columns are present but unused)

```
DATAFLOW,LAST UPDATE,freq,unit,month,geo,TIME_PERIOD,OBS_VALUE,OBS_FLAG,CONF_STATUS
ESTAT:DEMO_FMONTH(1.0),29/01/2026 23:00,Annual,Number,January,Ireland,2023,4551,,
ESTAT:DEMO_FMONTH(1.0),29/01/2026 23:00,Annual,Number,February,Ireland,2023,4102,,
```
