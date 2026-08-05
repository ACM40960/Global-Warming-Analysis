"""
04_solar_vs_co2_temperature.py
--------------------------------
Tests the second alternative hypothesis: "could the warming be driven by
increasing solar output (Total Solar Irradiance, TSI) rather than by
rising CO2?"

Solar_Irradiance.nc supplies a yearly reconstructed TSI series,
1610-2025 (NOAA/LASP composite, Coddington & Lean). We:

  1. Compute the secular trend in TSI over the full industrial period
     (1850-2025) AND separately over the modern instrumental-record
     period (1980-2025, satellite era) -- the two can behave very
     differently, which matters a lot for this question.
  2. Correlate TSI against global temperature and against CO2 over their
     respective overlap periods, with autocorrelation correction as
     before.
  3. Plot all three normalized series together (TSI, CO2, temperature)
     to visualize the "divergence" (or lack of it) between solar activity
     and warming in recent decades.
"""
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt

from utils import (
    RAW_DIR, FIG_DIR, DATA_OUT_DIR, apply_style, ols_trend_with_ci,
    pearson_autocorr_corrected,
)

apply_style()


def load_tsi():
    ds = xr.open_dataset(f"{RAW_DIR}/Solar_Irradiance.nc", use_cftime=True)
    years = np.array([t.year for t in ds["time"].values])
    tsi = ds["TSI"].values
    unc = ds["TSI_UNC"].values if "TSI_UNC" in ds else np.full_like(tsi, np.nan)
    return pd.DataFrame({"year": years, "tsi_wm2": tsi, "tsi_unc_wm2": unc})


def load_global_temp():
    return pd.read_csv(f"{DATA_OUT_DIR}/global_temperature_drift.csv")


def load_co2():
    df = pd.read_csv(f"{RAW_DIR}/Global_C02_Annual.csv", comment="#")
    df = df.rename(columns={"mean": "co2_ppm", "unc": "co2_unc_ppm"})
    df["year"] = df["year"].astype(int)
    return df[["year", "co2_ppm", "co2_unc_ppm"]].dropna(subset=["co2_ppm"])


def main():
    tsi = load_tsi()
    temp = load_global_temp()
    co2 = load_co2()
    tsi.to_csv(f"{DATA_OUT_DIR}/solar_tsi_annual.csv", index=False)

    # --- 1. TSI trends over two windows ---
    tsi_industrial = tsi[(tsi.year >= 1850) & (tsi.year <= 2025)]
    tsi_modern = tsi[(tsi.year >= 1980) & (tsi.year <= 2025)]
    trend_industrial = ols_trend_with_ci(tsi_industrial["year"], tsi_industrial["tsi_wm2"])
    trend_modern = ols_trend_with_ci(tsi_modern["year"], tsi_modern["tsi_wm2"])

    print("SOLAR IRRADIANCE (TSI) TREND")
    print(f"  1850-2025: {trend_industrial['per_decade']:+.4f} W/m^2 per decade "
          f"(95% CI [{trend_industrial['per_decade_ci_low']:+.4f}, "
          f"{trend_industrial['per_decade_ci_high']:+.4f}]), "
          f"R^2={trend_industrial['r_squared']:.3f}, p={trend_industrial['p_value']:.2e}")
    print(f"  1980-2025 (satellite era): {trend_modern['per_decade']:+.4f} W/m^2 per decade "
          f"(95% CI [{trend_modern['per_decade_ci_low']:+.4f}, "
          f"{trend_modern['per_decade_ci_high']:+.4f}]), "
          f"R^2={trend_modern['r_squared']:.3f}, p={trend_modern['p_value']:.2e}")

    # --- 2. correlations ---
    merged_temp = pd.merge(tsi, temp, on="year", how="inner")
    corr_tsi_temp = pearson_autocorr_corrected(merged_temp["tsi_wm2"], merged_temp["global_anomaly_C"])

    merged_co2 = pd.merge(tsi, co2, on="year", how="inner")
    corr_tsi_co2 = pearson_autocorr_corrected(merged_co2["tsi_wm2"], merged_co2["co2_ppm"])

    # modern-era (satellite, 1980-2025) correlation of TSI vs temperature -- this is
    # the period where the "solar-driven warming" hypothesis is most often tested,
    # since it's the best-observed TSI record
    merged_modern = merged_temp[(merged_temp.year >= 1980)]
    corr_modern = pearson_autocorr_corrected(merged_modern["tsi_wm2"], merged_modern["global_anomaly_C"])

    print(f"\nTSI vs global temperature, {merged_temp.year.min()}-{merged_temp.year.max()} "
          f"(n={len(merged_temp)}): r={corr_tsi_temp['r']:.3f}, "
          f"corrected p={corr_tsi_temp['p_autocorr_corrected']:.3f} (N_eff={corr_tsi_temp['n_eff']:.1f})")
    print(f"TSI vs global temperature, satellite era 1980-{merged_modern.year.max()} "
          f"(n={len(merged_modern)}): r={corr_modern['r']:.3f}, "
          f"corrected p={corr_modern['p_autocorr_corrected']:.3f} (N_eff={corr_modern['n_eff']:.1f})")
    print(f"TSI vs CO2, {merged_co2.year.min()}-{merged_co2.year.max()} "
          f"(n={len(merged_co2)}): r={corr_tsi_co2['r']:.3f}, "
          f"corrected p={corr_tsi_co2['p_autocorr_corrected']:.3f} (N_eff={corr_tsi_co2['n_eff']:.1f})")

    with open(f"{DATA_OUT_DIR}/solar_stats.txt", "w") as f:
        f.write("SOLAR IRRADIANCE (TSI) ANALYSIS\n" + "=" * 50 + "\n")
        f.write("TSI secular trend:\n")
        f.write(f"  1850-2025: {trend_industrial['per_decade']:+.4f} W/m^2/decade "
                f"(R^2={trend_industrial['r_squared']:.3f}, p={trend_industrial['p_value']:.2e})\n")
        f.write(f"  1980-2025: {trend_modern['per_decade']:+.4f} W/m^2/decade "
                f"(R^2={trend_modern['r_squared']:.3f}, p={trend_modern['p_value']:.2e})\n\n")
        f.write(f"TSI vs global temperature ({merged_temp.year.min()}-{merged_temp.year.max()}, "
                f"n={len(merged_temp)}): r={corr_tsi_temp['r']:.4f}, "
                f"corrected p={corr_tsi_temp['p_autocorr_corrected']:.4f} "
                f"(N_eff={corr_tsi_temp['n_eff']:.1f})\n")
        f.write(f"TSI vs global temperature, satellite era only (1980-{merged_modern.year.max()}, "
                f"n={len(merged_modern)}): r={corr_modern['r']:.4f}, "
                f"corrected p={corr_modern['p_autocorr_corrected']:.4f} "
                f"(N_eff={corr_modern['n_eff']:.1f})\n")
        f.write(f"TSI vs CO2 ({merged_co2.year.min()}-{merged_co2.year.max()}, n={len(merged_co2)}): "
                f"r={corr_tsi_co2['r']:.4f}, corrected p={corr_tsi_co2['p_autocorr_corrected']:.4f} "
                f"(N_eff={corr_tsi_co2['n_eff']:.1f})\n")

    # --- Figure 6: TSI full record with industrial period highlighted ---
    fig, ax = plt.subplots(figsize=(11, 5))
    ax.plot(tsi["year"], tsi["tsi_wm2"], color="#e67e22", lw=1)
    ax.fill_between(tsi["year"], tsi["tsi_wm2"] - tsi["tsi_unc_wm2"],
                     tsi["tsi_wm2"] + tsi["tsi_unc_wm2"], color="#f5cba7", alpha=0.5,
                     label="uncertainty band")
    ax.axvspan(1850, 2025, color="grey", alpha=0.08, label="industrial period")
    ax.set_title("Total Solar Irradiance, 1610-2025 (NOAA/LASP reconstruction)")
    ax.set_xlabel("Year")
    ax.set_ylabel("TSI (W/m^2)")
    ax.legend(loc="upper left")
    plt.tight_layout()
    plt.savefig(f"{FIG_DIR}/fig6_solar_tsi_full_record.png")
    plt.close()

    # --- Figure 7: normalized comparison, TSI vs CO2 vs temperature ---
    def zscore(s):
        return (s - s.mean()) / s.std()

    common = merged_co2.merge(temp[["year", "global_anomaly_C"]], on="year", how="inner")
    fig, ax = plt.subplots(figsize=(11, 5))
    ax.plot(common["year"], zscore(common["tsi_wm2"]), label="Solar irradiance (TSI, z-score)",
            color="#e67e22", lw=1.5)
    ax.plot(common["year"], zscore(common["co2_ppm"]), label="CO2 (z-score)",
            color="#7d3c98", lw=1.5)
    ax.plot(common["year"], zscore(common["global_anomaly_C"]), label="Global temperature (z-score)",
            color="#1e8449", lw=1.5)
    ax.set_title(f"Normalized comparison, {common.year.min()}-{common.year.max()}\n"
                 "(each series scaled to zero mean, unit std, for shape comparison only)")
    ax.set_xlabel("Year")
    ax.set_ylabel("Standardized anomaly (z-score)")
    ax.legend(loc="upper left")
    plt.tight_layout()
    plt.savefig(f"{FIG_DIR}/fig7_normalized_tsi_co2_temperature.png")
    plt.close()

    print("\nFigures saved to outputs/figures/")


if __name__ == "__main__":
    main()
