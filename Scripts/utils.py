"""
utils.py
--------
Shared constants, I/O helpers and statistical routines used by every
script in this project.

Statistical routines implemented here (kept dependency-light, numpy/scipy
only, except the multivariate regression step which uses statsmodels):

- ols_trend_with_ci(): linear trend + slope confidence interval
- mann_kendall_test(): non-parametric monotonic-trend test (robust to
  non-normal residuals and a handful of outliers)
- weighted_mean_std(): area-weighted mean & standard deviation across a
  spatial grid for a single time slice
- effective_n_ar1(): effective sample size correction for autocorrelated
  time series (Bretherton et al. 1999 / Bartlett), used to avoid
  overstating the significance of correlations between two trending
  series
- pearson_autocorr_corrected(): Pearson r, its naive p-value, and a
  p-value corrected for serial correlation using the effective N above
"""

import numpy as np
from scipy import stats

# ----------------------------------------------------------------------
# Project-wide constants
# ----------------------------------------------------------------------
RAW_DIR = "../data/raw"
FIG_DIR = "../outputs/figures"
DATA_OUT_DIR = "../outputs/data"

# Fraction of Earth's surface that is land vs ocean (standard geographic
# constants; used to area-weight the land-only and ocean-only series into
# a single global surface average).
LAND_FRACTION = 0.291
OCEAN_FRACTION = 1.0 - LAND_FRACTION

# Baseline period used to convert absolute SST to anomalies (matches the
# climatology period documented in the ERSSTv5 file metadata).
SST_BASELINE_YEARS = (1971, 2000)

PLOT_STYLE = {
    "figure.figsize": (10, 6),
    "figure.dpi": 120,
    "axes.grid": True,
    "grid.alpha": 0.3,
    "font.size": 11,
    "axes.titlesize": 13,
    "axes.titleweight": "bold",
}


def apply_style():
    import matplotlib.pyplot as plt
    plt.rcParams.update(PLOT_STYLE)


# ----------------------------------------------------------------------
# Trend statistics
# ----------------------------------------------------------------------
def ols_trend_with_ci(x, y, alpha=0.05):
    """
    Ordinary-least-squares linear trend of y on x, with a (1-alpha)
    confidence interval on the slope.

    Returns dict with slope, intercept, slope_se, ci_low, ci_high,
    r_value, p_value, per_decade (slope*10) and per_decade CI.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    x, y = x[mask], y[mask]
    n = len(x)

    res = stats.linregress(x, y)
    slope, intercept, r, p, se = res.slope, res.intercept, res.rvalue, res.pvalue, res.stderr
    tcrit = stats.t.ppf(1 - alpha / 2, df=n - 2)
    ci_low, ci_high = slope - tcrit * se, slope + tcrit * se

    return {
        "n": n,
        "slope_per_year": slope,
        "intercept": intercept,
        "slope_se": se,
        "ci_low_per_year": ci_low,
        "ci_high_per_year": ci_high,
        "per_decade": slope * 10,
        "per_decade_ci_low": ci_low * 10,
        "per_decade_ci_high": ci_high * 10,
        "r_value": r,
        "r_squared": r ** 2,
        "p_value": p,
    }


def mann_kendall_test(y):
    """
    Non-parametric Mann-Kendall trend test. Returns (S, Z, p_value, trend)
    where trend is one of 'increasing', 'decreasing', 'no trend' at
    alpha=0.05. Robust alternative to OLS that does not assume normally
    distributed residuals or homoscedasticity -- a useful cross-check on
    the OLS trend above.
    """
    y = np.asarray(y, dtype=float)
    y = y[np.isfinite(y)]
    n = len(y)
    s = 0
    for k in range(n - 1):
        s += np.sum(np.sign(y[k + 1:] - y[k]))

    # variance, with a simple tie correction
    unique, counts = np.unique(y, return_counts=True)
    tie_term = np.sum(counts * (counts - 1) * (2 * counts + 5))
    var_s = (n * (n - 1) * (2 * n + 5) - tie_term) / 18.0

    if s > 0:
        z = (s - 1) / np.sqrt(var_s)
    elif s < 0:
        z = (s + 1) / np.sqrt(var_s)
    else:
        z = 0.0

    p = 2 * (1 - stats.norm.cdf(abs(z)))
    if p < 0.05 and s > 0:
        trend = "increasing"
    elif p < 0.05 and s < 0:
        trend = "decreasing"
    else:
        trend = "no significant trend"

    return {"S": s, "Z": z, "p_value": p, "trend": trend}


# ----------------------------------------------------------------------
# Spatial (grid) statistics
# ----------------------------------------------------------------------
def weighted_mean_std(values, weights):
    """
    Weighted mean and weighted standard deviation across a spatial grid
    (e.g. all ocean grid cells contributing to a global mean in one year).
    values, weights: 1-D arrays of equal length (already flattened,
    NaNs removed by the caller).
    """
    values = np.asarray(values, dtype=float)
    weights = np.asarray(weights, dtype=float)
    mask = np.isfinite(values) & np.isfinite(weights) & (weights > 0)
    values, weights = values[mask], weights[mask]
    if len(values) == 0:
        return np.nan, np.nan, 0
    wmean = np.average(values, weights=weights)
    wvar = np.average((values - wmean) ** 2, weights=weights)
    return wmean, np.sqrt(wvar), len(values)


# ----------------------------------------------------------------------
# Correlation statistics that account for serial correlation
# ----------------------------------------------------------------------
def effective_n_ar1(x, y):
    """
    Effective sample size for the correlation between two autocorrelated
    time series, following Bretherton et al. (1999):

        N_eff = N * (1 - r1x*r1y) / (1 + r1x*r1y)

    where r1x, r1y are the lag-1 autocorrelations of each series. This is
    important here because both temperature and CO2 (and, to a lesser
    extent, TSI) are strongly autocorrelated / trending series -- a naive
    Pearson correlation on the raw levels will report an artificially
    tiny p-value. N_eff gives a more honest number of independent
    observations to test significance against.
    """
    def lag1_autocorr(v):
        v = np.asarray(v, dtype=float)
        v = v - np.mean(v)
        return np.sum(v[:-1] * v[1:]) / np.sum(v ** 2)

    n = len(x)
    r1x, r1y = lag1_autocorr(x), lag1_autocorr(y)
    denom = 1 + r1x * r1y
    if denom <= 0:
        return n  # degenerate case, fall back to naive N
    n_eff = n * (1 - r1x * r1y) / denom
    return max(3, min(n, n_eff))


def pearson_autocorr_corrected(x, y):
    """
    Pearson correlation coefficient together with:
      - the naive two-sided p-value (assumes N independent points)
      - an autocorrelation-corrected p-value using N_eff (see above)

    Returns a dict. Use the corrected p-value as the primary result when
    both series are trending/autocorrelated (as climate time series are).
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    x, y = x[mask], y[mask]
    n = len(x)

    r, p_naive = stats.pearsonr(x, y)
    n_eff = effective_n_ar1(x, y)

    # t-statistic recomputed with the effective sample size
    if n_eff > 2 and abs(r) < 1:
        t_eff = r * np.sqrt((n_eff - 2) / (1 - r ** 2))
        p_corrected = 2 * (1 - stats.t.cdf(abs(t_eff), df=n_eff - 2))
    else:
        p_corrected = np.nan

    return {
        "n": n,
        "n_eff": n_eff,
        "r": r,
        "p_naive": p_naive,
        "p_autocorr_corrected": p_corrected,
    }


def first_difference(y):
    """Year-over-year change; removes shared long-term trend so that the
    residual correlation reflects shared short-term co-variability rather
    than the fact that both series simply go up over time."""
    y = np.asarray(y, dtype=float)
    return np.diff(y)
