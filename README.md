# 🌍 Global Warming Attribution Project

A reproducible, statistics-first analysis of whether the observed rise in
global surface temperature since 1850 is better explained by rising
atmospheric CO2, volcanic activity, or solar irradiance — built entirely
from raw land, ocean, atmospheric, volcanic, and solar datasets.

!\[Python](https://img.shields.io/badge/python-3.10%2B-blue)
!\[License](https://img.shields.io/badge/license-MIT-green)
!\[Status](https://img.shields.io/badge/status-active-success)

\---

## Table of Contents

* [Project Overview](#project-overview)
* [Objectives](#objectives)
* [Methodology](#methodology)
* [Getting Started](#getting-started)
* [Challenges](#challenges)
* [Results](#results)
* [Conclusion](#conclusion)
* [Future Work](#future-work)
* [Project Structure](#project-structure)
* [Support](#support)
* [Maintainers](#maintainers)
* [Acknowledgements](#acknowledgements)

\---

## Project Overview

This project builds a global surface temperature record from scratch —
combining a land-average time series with a gridded sea-surface
temperature dataset (\~16,000 ocean "locations") — and then statistically
tests three competing explanations for the warming trend it finds:

1. **Rising atmospheric CO2** (the anthropogenic hypothesis)
2. **Increased volcanic activity** (a natural CO2-source alternative)
3. **Increased solar irradiance** (a natural energy-source alternative)

Rather than asserting a conclusion, the pipeline runs rigorous,
autocorrelation-aware statistics (OLS trends with confidence intervals,
Mann-Kendall non-parametric trend tests, effective-sample-size-corrected
correlations, and a multivariate regression with HAC standard errors) and
lets the data speak. Every number in the final report is reproducible by
re-running the pipeline against the raw source files.

## Objectives

* \[x] Build a global-surface temperature anomaly record with proper
spatial coverage of the ocean (not just a handful of stations)
* \[x] Attach honest, statistically-derived error bars / uncertainty to
every time series produced
* \[x] Quantify the correlation between atmospheric CO2 and temperature
using methods that correct for the fact that both series are
trending and autocorrelated
* \[x] Test whether volcanic activity can account for the CO2 increase
* \[x] Test whether solar irradiance can account for the temperature
increase
* \[x] Combine all three candidate drivers into a single model to compare
their relative explanatory power

## Methodology

|Step|Script|What it does|
|-|-|-|
|1|`01\\\_process\\\_temperature.py`|Converts gridded monthly SST into per-cell annual anomalies (vs. 1971–2000 climatology), area-weights all ocean cells into a global ocean series, combines with the land series (29.1% / 70.9% area weighting), propagates uncertainty in quadrature, and validates the result against an independent reference product|
|2|`02\\\_co2\\\_temperature\\\_correlation.py`|Correlates CO2 against temperature three ways: raw levels (autocorrelation-corrected), year-over-year first differences, and a −5..+5 year lag scan|
|3|`03\\\_volcanic\\\_vs\\\_co2.py`|Builds annual eruption-count and VEI-weighted magnitude indices from the Smithsonian Global Volcanism Program catalog, tests their own secular trend, and correlates them against CO2|
|4|`04\\\_solar\\\_vs\\\_co2\\\_temperature.py`|Tests the Total Solar Irradiance trend over the full industrial era vs. the modern satellite era, and correlates TSI against CO2 and temperature|
|5|`05\\\_multivariate\\\_regression.py`|Puts standardized CO2, TSI, and volcanic magnitude into one OLS regression against temperature, using Newey-West (HAC) standard errors to account for autocorrelation|

**Key statistical safeguards used throughout:**

* Ordinary least squares trends always reported with 95% confidence
intervals, not just a point estimate
* Mann-Kendall non-parametric trend test as a cross-check on OLS (robust
to non-normal residuals / outliers)
* Effective-sample-size correction (Bretherton et al., 1999) applied to
every correlation between two autocorrelated/trending series, since a
naive Pearson p-value on trending climate series is misleadingly small
* Area-weighted (not simple-average) spatial statistics for the gridded
ocean data
* HAC/Newey-West standard errors in the multivariate regression

## Getting Started

### Prerequisites

* Python 3.10+
* \~200 MB free disk space (mainly for the sea-surface-temperature NetCDF file)

### Installation

```bash
# clone the repository
git clone https://github.com/<your-username>/global-warming-attribution.git
cd global-warming-attribution

# create and activate a virtual environment
python3 -m venv venv
source venv/bin/activate        # Windows: venv\\\\Scripts\\\\activate

# install dependencies
pip install -r requirements.txt
```

### Data setup

Place the six raw source files in `data/raw/` (see
[Project Structure](#project-structure) for the exact expected filenames).
`Sea\\\_Temperature.nc` (\~150 MB) is not tracked in this repository — download
it separately and drop it into `data/raw/` before running.

### Running the pipeline

```bash
cd scripts
python3 run\\\_all.py
```

This runs all five steps in order (\~1–2 minutes, dominated by regridding
the sea-surface NetCDF file) and writes every CSV, statistics summary, and
figure to `outputs/`. Each step can also be run individually, e.g.:

```bash
python3 03\\\_volcanic\\\_vs\\\_co2.py
```

## Challenges

* **Autocorrelation inflating apparent significance.** Two smoothly
trending time series (like CO2 and temperature) will show a very high
Pearson correlation almost regardless of whether they're causally
related. This was addressed by computing an effective sample size from
each series' lag-1 autocorrelation before testing significance, and by
cross-checking the level correlation against a first-difference
correlation that removes the shared trend entirely.
* **Reporting bias in historical volcanic records.** Raw eruption *counts*
show an artificial upward trend simply because global volcanic
monitoring has improved since 1850 — small eruptions in remote regions
are far more likely to be recorded today. This was mitigated by also
computing a VEI-magnitude-weighted index, which is dominated by large,
historically well-documented eruptions and is far less sensitive to
this bias.
* **Missing native uncertainty in the land temperature file.** The
provided land-temperature series ships no per-year uncertainty column.
A rolling-residual proxy (spread around a 15-year smoothed mean) was
used as a reasonable stand-in error bar.
* **Large, sparsely-documented NetCDF files.** The gridded ocean file
needed a per-cell climatology and month-by-month anomaly calculation
before any spatial averaging was valid — averaging absolute
temperatures directly (without converting to anomalies first) would
have produced a meaningless series dominated by geography rather than
by warming.
* **Short overlap between datasets.** The CO2 record only starts in 1979,
which caps every CO2-involving statistical test at 47 years of data —
a real limitation, discussed openly in the results rather than hidden.

## Results

* Global surface temperature has risen **+0.059 °C/decade** (95% CI:
+0.052 to +0.065) since 1854, confirmed by both OLS (R² = 0.66,
p ≈ 2×10⁻⁴¹) and a Mann-Kendall trend test (p < 10⁻³⁰). The
reconstruction validates at **r = 0.998** against an independent
reference product.
* Atmospheric CO2 has risen **+18.9 ppm/decade** since 1979 (R² = 0.99)
and correlates strongly with temperature (r = 0.95), remaining
significant (p = 0.012) even after correcting for autocorrelation.
* Volcanic activity shows **no significant secular trend** once
historical under-reporting is accounted for (magnitude index:
R² = 0.0003, p = 0.77), and **does not correlate with CO2**
(corrected p = 0.33–0.45).
* Solar irradiance has been **statistically flat since 1980**
(p = 0.95) and does not correlate significantly with CO2 or with
recent temperature.
* In a combined regression, CO2's standardized effect on temperature
(+0.268) is roughly **10x larger** than solar's (+0.026, only
marginally significant), and volcanic activity's effect is
statistically indistinguishable from zero.

Full write-up, with every supporting figure, is in
[`CASE\\\_REPORT.md`](CASE_REPORT.md).

## Conclusion

Every statistical test in this pipeline points the same direction: the
observed warming trend tracks atmospheric CO2 closely, even after
correcting for shared autocorrelation, while the two leading natural
alternatives — volcanic activity and solar irradiance — show neither a
comparable long-term trend nor a significant correlation with either CO2
or temperature over the periods tested.

**A necessary caveat:** correlation and regression, however carefully
corrected, cannot establish causation on their own. The causal case for
CO2-driven warming ultimately rests on the independently well-established
physics of CO2's infrared absorption spectrum (the greenhouse effect),
which is outside the scope of a purely statistical analysis. What this
project shows is that the statistical pattern in the data is fully
consistent with the CO2-driven warming hypothesis and is not consistent
with the volcanic or solar alternatives, for the specific datasets
analyzed here.

## Future Work

* Extend the CO2 record further back using ice-core proxy data to test
correlations over the full 1850–2025 industrial period rather than
just 1979–2025
* Add ocean heat content (0–2000m) as a fourth independent line of
evidence, since it is far less affected by short-term atmospheric
noise than surface temperature
* Incorporate volcanic aerosol optical depth (not just eruption
count/VEI) to test short-term volcanic *cooling* effects on
temperature, distinct from the CO2-source question tested here
* Add ENSO (El Niño/La Niña) indices as a control variable in the
multivariate regression to explain more of the year-to-year
variability not accounted for by CO2, solar, or volcanic activity
* Regional/hemispheric breakdowns of the ocean grid, rather than a single
global ocean average
* Package the pipeline as an installable CLI (`pip install .`) with a
single `attribution-pipeline run` entry point

## Project Structure

```
project/
├── README.md                  <- this file
├── CASE\\\_REPORT.md              <- full written results with figures
├── requirements.txt
├── data/
│   └── raw/
│       ├── Land\\\_Temperature.asc
│       ├── Global\\\_Temperature.asc
│       ├── Sea\\\_Temperature.nc        (not tracked in git — see Getting Started)
│       ├── Solar\\\_Irradiance.nc
│       ├── Global\\\_C02\\\_Annual.csv
│       └── Volcanic\\\_Activity.tsv
├── scripts/
│   ├── utils.py                          # shared stats/plotting helpers
│   ├── 01\\\_process\\\_temperature.py
│   ├── 02\\\_co2\\\_temperature\\\_correlation.py
│   ├── 03\\\_volcanic\\\_vs\\\_co2.py
│   ├── 04\\\_solar\\\_vs\\\_co2\\\_temperature.py
│   ├── 05\\\_multivariate\\\_regression.py
│   └── run\\\_all.py
└── outputs/
    ├── data/       <- CSVs + plain-text statistics summaries
    └── figures/    <- all generated plots (PNG)
```

## Support

If you run into an issue or have a question:

* Open a [GitHub Issue](../../issues) with a clear description and, if
applicable, the console output or traceback
* For questions about the underlying datasets, see the source citations
in [Acknowledgements](#acknowledgements)

## Maintainers

* **\[Your Name]** — [@your-github-handle](https://github.com/your-github-handle)

*(Update this section with your own name/handle before publishing.)*

## Acknowledgements

This project builds on publicly available data from:

* **NOAA Global Monitoring Laboratory** — global annual mean atmospheric
CO2 (gml.noaa.gov/ccgg/trends/)
* **NOAA NCEI — ERSSTv5** — Extended Reconstructed Sea Surface
Temperature v5, gridded ocean data
* **NOAA/LASP** — Total Solar Irradiance composite reconstruction
(Coddington \& Lean)
* **Smithsonian Institution — Global Volcanism Program** — historical
eruption catalog with VEI classifications
* Land and blended land+ocean temperature anomaly products used for
validation

Thanks to the maintainers of `numpy`, `pandas`, `scipy`, `xarray`,
`matplotlib`, and `statsmodels`, on which this entire analysis depends.

