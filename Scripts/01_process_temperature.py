"""
01_process_temperature.py
--------------------------
Builds the global surface temperature record from the raw inputs:

  1. LAND: Land_Temperature.asc is already a global-land-average annual
     anomaly series (one number per year). We keep it as-is but derive a
     year-by-year variability estimate (a proxy "error bar") from the
     residual of the series about a 15-year LOWESS-like rolling smooth,
     since the file does not ship its own per-year uncertainty (columns
     3-6 are all -999 = missing).

  2. OCEAN: Sea_Temperature.nc is a gridded product (ERSSTv5): monthly
     sea-surface temperature on an 89 x 180 lat/lon grid, 1854-2026. This
     is genuinely "temperature vs time for each location" -- each grid
     cell is one location. For each grid cell we:
        a. convert to a monthly anomaly relative to its own 1971-2000
           climatology (per calendar month, per cell)
        b. average the 12 months of each year -> one annual anomaly per
           cell per year
     We then combine ALL ocean locations into a single global ocean
     average per year using latitude (cos-weighted) area weighting, and
     report the area-weighted standard deviation across cells as the
     spatial spread (error bar) for that year.

  3. GLOBAL: land and ocean annual series are combined into a single
     global-surface average using the true land/ocean area fractions
     (29.1% / 70.9%), with errors propagated in quadrature.

  4. We also load the provided Global_Temperature.asc (an independent,
     pre-blended land+ocean product) purely as a cross-check / sanity
     comparison against our own from-scratch combination.

Outputs
-------
  outputs/data/land_annual.csv
  outputs/data/ocean_annual.csv
  outputs/data/global_temperature_drift.csv   <-- the requested
                                                    "different file"
  outputs/figures/fig1_land_ocean_global_temperature.png
  outputs/figures/fig2_validation_vs_reference.png
"""
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt

from utils import (
    RAW_DIR, FIG_DIR, DATA_OUT_DIR, LAND_FRACTION, OCEAN_FRACTION,
    SST_BASELINE_YEARS, apply_style, ols_trend_with_ci, mann_kendall_test,
    weighted_mean_std,
)

apply_style()


# ----------------------------------------------------------------------
# 1. LAND
# ----------------------------------------------------------------------
def load_land():
    arr = np.loadtxt(f"{RAW_DIR}/Land_Temperature.asc")
    years, anom = arr[:, 0].astype(int), arr[:, 1]
    df = pd.DataFrame({"year": years, "land_anomaly_C": anom})

    # Proxy year-to-year uncertainty: residual std about a rolling
    # 15-year centered mean (captures interannual variability not
    # explained by the slow trend/smooth).
    roll = df["land_anomaly_C"].rolling(15, center=True, min_periods=5).mean()
    resid_std = (df["land_anomaly_C"] - roll).std()
    df["land_anomaly_err_C"] = resid_std
    return df


# ----------------------------------------------------------------------
# 2. OCEAN  (gridded, per-location -> combined)
# ----------------------------------------------------------------------
def load_ocean():
    ds = xr.open_dataset(f"{RAW_DIR}/Sea_Temperature.nc")
    sst = ds["sst"]  # (time, lat, lon), degrees C, monthly

    # restrict to full calendar years only
    years_all = sst["time"].dt.year.values
    full_years = [y for y in np.unique(years_all)
                  if np.sum(years_all == y) == 12]
    sst = sst.sel(time=sst["time"].dt.year.isin(full_years))

    # climatology per grid cell per calendar month over the baseline period
    y0, y1 = SST_BASELINE_YEARS
    baseline = sst.sel(time=slice(f"{y0}-01-01", f"{y1}-12-31"))
    clim = baseline.groupby(baseline["time"].dt.month).mean("time")

    anom = sst.groupby(sst["time"].dt.month) - clim  # monthly anomaly, per cell
    annual = anom.groupby(anom["time"].dt.year).mean("time")  # (year, lat, lon)

    # cos(latitude) area weights, broadcast over lon
    lat = annual["lat"].values
    lon = annual["lon"].values
    latw = np.cos(np.deg2rad(lat))
    weight_grid = np.outer(latw, np.ones_like(lon))  # (lat, lon)

    years = annual["year"].values
    means, stds, ncells = [], [], []
    for i, yr in enumerate(years):
        slice_vals = annual.isel(year=i).values  # (lat, lon)
        wmean, wstd, n = weighted_mean_std(slice_vals.ravel(), weight_grid.ravel())
        means.append(wmean)
        stds.append(wstd)
        ncells.append(n)

    df = pd.DataFrame({
        "year": years.astype(int),
        "ocean_anomaly_C": means,
        "ocean_anomaly_spatial_std_C": stds,
        "ocean_n_locations": ncells,
    })
    return df


# ----------------------------------------------------------------------
# 3. COMBINE -> GLOBAL SURFACE AVERAGE
# ----------------------------------------------------------------------
def combine(land_df, ocean_df):
    df = pd.merge(land_df, ocean_df, on="year", how="inner")
    df["global_anomaly_C"] = (
        LAND_FRACTION * df["land_anomaly_C"] + OCEAN_FRACTION * df["ocean_anomaly_C"]
    )
    # propagate error in quadrature (land error is a constant proxy value,
    # ocean error is the spatial std of that year's grid-cell values)
    df["global_anomaly_err_C"] = np.sqrt(
        (LAND_FRACTION * df["land_anomaly_err_C"]) ** 2
        + (OCEAN_FRACTION * df["ocean_anomaly_spatial_std_C"]) ** 2
    )
    return df


def load_reference():
    arr = np.loadtxt(f"{RAW_DIR}/Global_Temperature.asc")
    return pd.DataFrame({"year": arr[:, 0].astype(int), "reference_anomaly_C": arr[:, 1]})


def main():
    print("Loading land data...")
    land = load_land()
    print("Loading + gridding ocean data (this can take ~30-60s)...")
    ocean = load_ocean()
    print("Combining into global surface average...")
    glob = combine(land, ocean)

    land.to_csv(f"{DATA_OUT_DIR}/land_annual.csv", index=False)
    ocean.to_csv(f"{DATA_OUT_DIR}/ocean_annual.csv", index=False)
    glob.to_csv(f"{DATA_OUT_DIR}/global_temperature_drift.csv", index=False)
    print(f"Saved: {DATA_OUT_DIR}/global_temperature_drift.csv  "
          f"({len(glob)} years, {glob.year.min()}-{glob.year.max()})")

    # --- trend statistics ---
    trend = ols_trend_with_ci(glob["year"], glob["global_anomaly_C"])
    mk = mann_kendall_test(glob["global_anomaly_C"])
    print("\nGLOBAL SURFACE TEMPERATURE TREND")
    print(f"  OLS trend: {trend['per_decade']:+.4f} C/decade "
          f"(95% CI [{trend['per_decade_ci_low']:+.4f}, {trend['per_decade_ci_high']:+.4f}]), "
          f"R^2={trend['r_squared']:.3f}, p={trend['p_value']:.2e}")
    print(f"  Mann-Kendall: {mk['trend']} (Z={mk['Z']:.2f}, p={mk['p_value']:.2e})")

    with open(f"{DATA_OUT_DIR}/temperature_trend_stats.txt", "w") as f:
        f.write("GLOBAL SURFACE TEMPERATURE TREND STATISTICS\n")
        f.write("=" * 50 + "\n")
        f.write(f"Period: {glob.year.min()}-{glob.year.max()} (n={trend['n']} years)\n\n")
        f.write("Ordinary Least Squares:\n")
        f.write(f"  slope           = {trend['per_decade']:+.4f} C/decade\n")
        f.write(f"  95% CI          = [{trend['per_decade_ci_low']:+.4f}, "
                 f"{trend['per_decade_ci_high']:+.4f}] C/decade\n")
        f.write(f"  R^2             = {trend['r_squared']:.4f}\n")
        f.write(f"  p-value         = {trend['p_value']:.3e}\n\n")
        f.write("Mann-Kendall non-parametric trend test:\n")
        f.write(f"  result          = {mk['trend']}\n")
        f.write(f"  Z statistic     = {mk['Z']:.3f}\n")
        f.write(f"  p-value         = {mk['p_value']:.3e}\n")

    # --- Figure 1: land / ocean / global with error bars ---
    fig, axes = plt.subplots(3, 1, figsize=(11, 12), sharex=True)

    axes[0].errorbar(land["year"], land["land_anomaly_C"],
                      yerr=land["land_anomaly_err_C"], fmt="o-", ms=3,
                      lw=1, capsize=2, color="#c0392b", ecolor="#e6a19a",
                      elinewidth=0.8, label="Land annual anomaly")
    axes[0].set_title("Land surface temperature anomaly (global-land average)")
    axes[0].set_ylabel("Anomaly (C)")
    axes[0].legend(loc="upper left")

    axes[1].errorbar(ocean["year"], ocean["ocean_anomaly_C"],
                      yerr=ocean["ocean_anomaly_spatial_std_C"], fmt="o-", ms=3,
                      lw=1, capsize=2, color="#2471a3", ecolor="#a9cce3",
                      elinewidth=0.8, label="Ocean annual anomaly\n(area-weighted over grid cells)")
    axes[1].set_title("Sea-surface temperature anomaly\n"
                       "(combined from all ocean grid-cell locations)")
    axes[1].set_ylabel("Anomaly (C)")
    axes[1].legend(loc="upper left")

    axes[2].errorbar(glob["year"], glob["global_anomaly_C"],
                      yerr=glob["global_anomaly_err_C"], fmt="o-", ms=3,
                      lw=1, capsize=2, color="#1e8449", ecolor="#a9dfbf",
                      elinewidth=0.8, label="Global surface anomaly\n(29.1% land + 70.9% ocean)")
    trend_line = trend["intercept"] + trend["slope_per_year"] * glob["year"]
    axes[2].plot(glob["year"], trend_line, "--", color="black", lw=1.5,
                 label=f"OLS trend: {trend['per_decade']:+.3f} C/decade")
    axes[2].set_title("Global surface temperature anomaly (land + ocean combined)")
    axes[2].set_ylabel("Anomaly (C)")
    axes[2].set_xlabel("Year")
    axes[2].legend(loc="upper left")

    plt.tight_layout()
    plt.savefig(f"{FIG_DIR}/fig1_land_ocean_global_temperature.png")
    plt.close()

    # --- Figure 2: validation against provided reference product ---
    ref = load_reference()
    comp = pd.merge(glob[["year", "global_anomaly_C"]], ref, on="year", how="inner")
    r = np.corrcoef(comp["global_anomaly_C"], comp["reference_anomaly_C"])[0, 1]

    fig, ax = plt.subplots(figsize=(11, 5))
    ax.plot(comp["year"], comp["global_anomaly_C"], label="This project's global average\n(built from land + gridded ocean)", lw=1.5)
    ax.plot(comp["year"], comp["reference_anomaly_C"], label="Provided reference product\n(Global_Temperature.asc)", lw=1.5, ls="--")
    ax.set_title(f"Validation: reconstructed global average vs. reference product\n(Pearson r = {r:.3f})")
    ax.set_xlabel("Year")
    ax.set_ylabel("Anomaly (C)")
    ax.legend(loc="upper left")
    plt.tight_layout()
    plt.savefig(f"{FIG_DIR}/fig2_validation_vs_reference.png")
    plt.close()

    print(f"\nValidation correlation vs reference product: r = {r:.3f}")
    print("Figures saved to outputs/figures/")


if __name__ == "__main__":
    main()
