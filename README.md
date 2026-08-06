# Global Warming Attribution Project

A reproducible analysis pipeline that builds a global surface temperature
record from raw land + ocean data, and statistically tests three candidate
explanations for the observed warming: rising atmospheric CO2, volcanic
activity, and solar irradiance.

## 1. What's in this project

```
project/
├── README.md                  <- you are here
├── requirements.txt
├── data/raw/                  <- put the 6 source files here (see below)
├── scripts/
│   ├── utils.py                          shared stats/plot helpers
│   ├── 01_process_temperature.py         land + ocean -> global temp record
│   ├── 02_co2_temperature_correlation.py CO2 vs temperature statistics
│   ├── 03_volcanic_vs_co2.py             volcanic activity vs CO2
│   ├── 04_solar_vs_co2_temperature.py    solar irradiance vs CO2/temp
│   ├── 05_multivariate_regression.py     all three drivers, one model
│   └── run_all.py                        runs the five steps in order
└── outputs/
    ├── data/       <- CSVs and .txt statistics summaries
    └── figures/    <- all PNG plots (8 total)
```

## 2. Source data (place these in `data/raw/` before running)

| File | Contents |
|---|---|
| `Land_Temperature.asc` | Global land-average annual temperature anomaly, 1850-2025 |
| `Global_Temperature.asc` | Independent pre-blended land+ocean reference product (used only for validation) |
| `Sea_Temperature.nc` | ERSSTv5 gridded monthly sea-surface temperature, 89x180 grid, 1854-2026 |
| `Solar_Irradiance.nc` | NOAA/LASP yearly Total Solar Irradiance reconstruction, 1610-2025 |
| `Global_C02_Annual.csv` | NOAA GML global annual mean atmospheric CO2, 1979-2025 |
| `Volcanic_Activity.tsv` | Smithsonian Global Volcanism Program eruption catalog, 1850-2025 |

These are the six files you uploaded; the scripts expect them unmodified,
under `data/raw/`, with exactly these filenames.

## 3. How to run it

```bash
# 1. (recommended) create a virtual environment
python3 -m venv venv
source venv/bin/activate        # Windows: venv\Scripts\activate

# 2. install dependencies
pip install -r requirements.txt

# 3. make sure the 6 raw files are in data/raw/, then run the pipeline
cd scripts
python3 run_all.py
```

This takes roughly 1-2 minutes, dominated by regridding the ~150 MB
sea-surface-temperature NetCDF file in step 1. You can also run any single
step on its own (e.g. `python3 03_volcanic_vs_co2.py`) as long as the
earlier steps that produce its inputs have been run at least once.

All results are written to `outputs/data/` (CSVs + plain-text statistics
summaries) and `outputs/figures/` (PNG plots, 8 total).

**The single most important output** is
`outputs/data/global_temperature_drift.csv` — the from-scratch global
land+ocean surface temperature series with propagated uncertainty
(standard-deviation error bars), built entirely from the raw land and
gridded-ocean inputs.

## 4. Method summary

**Step 1 — building the global temperature record.**
`Land_Temperature.asc` is already a pre-computed global land average, so it
is used as-is (with a rolling-residual proxy for its year-to-year
uncertainty, since the file ships no native uncertainty column).
`Sea_Temperature.nc` is genuinely gridded: 89 x 180 = up to 16,020
individual ocean "locations", each with its own monthly time series. For
every grid cell we compute a monthly anomaly relative to that cell's own
1971-2000 climatology, average to annual, and then combine **all** ocean
locations into one global ocean series using latitude (cos-weighted) area
weighting — the spread across cells in a given year becomes that year's
ocean error bar. Land and ocean are then combined using the true Earth
land/ocean area fractions (29.1% / 70.9%), propagating uncertainty in
quadrature. The result is cross-checked against the independent
`Global_Temperature.asc` reference product (Pearson r = 0.998 — see
`fig2_validation_vs_reference.png`).

**Step 2 — CO2 vs temperature.**
Because both series trend upward together, a naive correlation would look
almost perfect even for two unrelated trending series. So the script
reports three things: (a) the level correlation with a p-value corrected
for serial autocorrelation (effective-N method, Bretherton et al. 1999),
(b) a first-difference (year-over-year change) correlation, which removes
the shared trend entirely, and (c) a -5..+5 year lag-correlation scan.

**Step 3 — volcanic activity vs CO2.**
Builds two independent annual indices from the eruption catalog (raw
eruption count, and a VEI-weighted magnitude index ~10^VEI per eruption),
tests whether either has a genuine secular trend over 1850-2025, and
correlates both against CO2 over the years the two datasets overlap.

**Step 4 — solar irradiance vs CO2 and temperature.**
Tests the TSI trend separately over the full industrial period and the
modern satellite era, then correlates TSI against both CO2 and
temperature.

**Step 5 — multivariate regression.**
All three standardized predictors (CO2, TSI, volcanic magnitude) are
entered into one OLS model of temperature, with Newey-West (HAC)
standard errors to account for autocorrelation, so their relative
contributions can be compared directly and fairly.

## 5. Outputs

| File | Description |
|---|---|
| `outputs/data/global_temperature_drift.csv` | **The global surface temperature series** (year, land/ocean/global anomaly, error bars) |
| `outputs/data/land_annual.csv`, `ocean_annual.csv` | Intermediate land-only / ocean-only series |
| `outputs/data/co2_temperature_merged.csv` | Merged CO2 + temperature table |
| `outputs/data/co2_temperature_lag_correlation.csv` | Lag-correlation scan results |
| `outputs/data/volcanic_annual_index.csv` | Annual eruption count + magnitude index |
| `outputs/data/solar_tsi_annual.csv` | Cleaned yearly TSI series |
| `outputs/data/*_stats.txt` | Plain-text statistical summaries for each step |
| `outputs/figures/fig1..fig8*.png` | All plots (see report for descriptions) |

See **`CASE_REPORT.md`** (in the project root) for the full written-up
findings with all figures embedded.
