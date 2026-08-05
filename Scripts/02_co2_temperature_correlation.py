"""
02_co2_temperature_correlation.py
----------------------------------
Loads the NOAA global annual mean CO2 record and correlates it against
the global surface temperature series built in step 1.

Because both CO2 and temperature are strongly autocorrelated / trending
series, a naive Pearson correlation on their raw levels will look
almost perfectly significant even if the *shared driver* were something
else entirely (two things that both rise monotonically will always
correlate well). To guard against this we compute, and report,
THREE separate statistics:

  1. Level correlation, with an autocorrelation-corrected p-value
     (effective-N following Bretherton et al. 1999).
  2. First-difference correlation: correlates YEAR-OVER-YEAR CHANGES in
     CO2 against year-over-year changes in temperature. This removes the
     shared long-term trend and asks a much sharper question: do the
     wiggles/accelerations in CO2 track the wiggles in temperature on a
     year-to-year basis?
  3. Lag correlation: correlates temperature against CO2 shifted by
     -5..+5 years, to check whether the relationship is contemporaneous
     or lagged.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats

from utils import (
    RAW_DIR, FIG_DIR, DATA_OUT_DIR, apply_style, ols_trend_with_ci,
    pearson_autocorr_corrected, first_difference,
)

apply_style()


def load_co2():
    df = pd.read_csv(f"{RAW_DIR}/Global_C02_Annual.csv", comment="#")
    df = df.rename(columns={"year": "year", "mean": "co2_ppm", "unc": "co2_unc_ppm"})
    df = df[["year", "co2_ppm", "co2_unc_ppm"]].dropna(subset=["co2_ppm"])
    df["year"] = df["year"].astype(int)
    return df


def load_global_temp():
    return pd.read_csv(f"{DATA_OUT_DIR}/global_temperature_drift.csv")


def main():
    co2 = load_co2()
    temp = load_global_temp()
    df = pd.merge(co2, temp, on="year", how="inner").sort_values("year").reset_index(drop=True)
    df.to_csv(f"{DATA_OUT_DIR}/co2_temperature_merged.csv", index=False)
    print(f"Overlap period: {df.year.min()}-{df.year.max()} (n={len(df)} years)")

    # --- CO2 trend ---
    co2_trend = ols_trend_with_ci(df["year"], df["co2_ppm"])
    print(f"\nCO2 trend: {co2_trend['per_decade']:+.3f} ppm/decade "
          f"(R^2={co2_trend['r_squared']:.3f}, p={co2_trend['p_value']:.2e})")

    # --- 1. Level correlation (autocorrelation-corrected) ---
    level_corr = pearson_autocorr_corrected(df["co2_ppm"], df["global_anomaly_C"])
    print("\n[1] LEVEL correlation (CO2 ppm vs global temperature anomaly):")
    print(f"    r = {level_corr['r']:.3f}")
    print(f"    naive p-value                = {level_corr['p_naive']:.2e}  (N={level_corr['n']})")
    print(f"    autocorr-corrected p-value    = {level_corr['p_autocorr_corrected']:.2e}  "
          f"(N_eff={level_corr['n_eff']:.1f})")

    # --- 2. First-difference correlation ---
    d_co2 = first_difference(df["co2_ppm"].values)
    d_temp = first_difference(df["global_anomaly_C"].values)
    diff_corr = pearson_autocorr_corrected(d_co2, d_temp)
    print("\n[2] FIRST-DIFFERENCE correlation (year-over-year changes):")
    print(f"    r = {diff_corr['r']:.3f}")
    print(f"    naive p-value                = {diff_corr['p_naive']:.2e}  (N={diff_corr['n']})")
    print(f"    autocorr-corrected p-value    = {diff_corr['p_autocorr_corrected']:.2e}  "
          f"(N_eff={diff_corr['n_eff']:.1f})")

    # --- 3. Lag correlation ---
    lags = range(-5, 6)
    lag_results = []
    for lag in lags:
        # positive lag: CO2 leads temperature by `lag` years
        if lag >= 0:
            x = df["co2_ppm"].values[:len(df) - lag] if lag > 0 else df["co2_ppm"].values
            y = df["global_anomaly_C"].values[lag:]
        else:
            x = df["co2_ppm"].values[-lag:]
            y = df["global_anomaly_C"].values[:len(df) + lag]
        r, p = stats.pearsonr(x, y)
        lag_results.append({"lag_years": lag, "r": r, "p_naive": p})
    lag_df = pd.DataFrame(lag_results)
    lag_df.to_csv(f"{DATA_OUT_DIR}/co2_temperature_lag_correlation.csv", index=False)
    best_lag = lag_df.loc[lag_df["r"].idxmax()]
    print(f"\n[3] LAG correlation: strongest r={best_lag.r:.3f} at lag={int(best_lag.lag_years)} years "
          "(positive lag = CO2 change precedes temperature change)")

    # --- save summary ---
    with open(f"{DATA_OUT_DIR}/co2_temperature_correlation_stats.txt", "w") as f:
        f.write("CO2 vs GLOBAL TEMPERATURE CORRELATION\n" + "=" * 50 + "\n")
        f.write(f"Overlap period: {df.year.min()}-{df.year.max()} (n={len(df)})\n\n")
        f.write(f"CO2 trend: {co2_trend['per_decade']:+.3f} ppm/decade "
                f"(R^2={co2_trend['r_squared']:.3f}, p={co2_trend['p_value']:.2e})\n\n")
        f.write("[1] Level correlation:\n")
        f.write(f"    r = {level_corr['r']:.4f}\n")
        f.write(f"    naive p = {level_corr['p_naive']:.3e} (N={level_corr['n']})\n")
        f.write(f"    autocorrelation-corrected p = {level_corr['p_autocorr_corrected']:.3e} "
                f"(N_eff={level_corr['n_eff']:.1f})\n\n")
        f.write("[2] First-difference (year-over-year change) correlation:\n")
        f.write(f"    r = {diff_corr['r']:.4f}\n")
        f.write(f"    naive p = {diff_corr['p_naive']:.3e} (N={diff_corr['n']})\n")
        f.write(f"    autocorrelation-corrected p = {diff_corr['p_autocorr_corrected']:.3e} "
                f"(N_eff={diff_corr['n_eff']:.1f})\n\n")
        f.write("[3] Lag correlation (CO2 leading temperature by N years):\n")
        f.write(lag_df.to_string(index=False))
        f.write("\n")

    # --- Figure 3: dual-axis time series ---
    fig, ax1 = plt.subplots(figsize=(11, 5))
    ax2 = ax1.twinx()
    ax1.errorbar(df["year"], df["global_anomaly_C"], yerr=df["global_anomaly_err_C"],
                 fmt="o-", ms=3, lw=1.2, capsize=2, color="#1e8449", ecolor="#a9dfbf",
                 elinewidth=0.8, label="Global temperature anomaly")
    ax2.errorbar(df["year"], df["co2_ppm"], yerr=df["co2_unc_ppm"], fmt="s--", ms=3, lw=1.2,
                 capsize=2, color="#7d3c98", ecolor="#d2b4de", elinewidth=0.8,
                 label="Atmospheric CO2")
    ax1.set_xlabel("Year")
    ax1.set_ylabel("Temperature anomaly (C)", color="#1e8449")
    ax2.set_ylabel("CO2 (ppm)", color="#7d3c98")
    ax1.tick_params(axis="y", colors="#1e8449")
    ax2.tick_params(axis="y", colors="#7d3c98")
    ax1.set_title(f"Global temperature and atmospheric CO2, {df.year.min()}-{df.year.max()}\n"
                  f"(level correlation r={level_corr['r']:.3f}, "
                  f"autocorr-corrected p={level_corr['p_autocorr_corrected']:.1e})")
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc="upper left")
    plt.tight_layout()
    plt.savefig(f"{FIG_DIR}/fig3_co2_vs_temperature_timeseries.png")
    plt.close()

    # --- Figure 4: scatter with regression + first-difference scatter ---
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    axes[0].scatter(df["co2_ppm"], df["global_anomaly_C"], c=df["year"], cmap="viridis", s=25)
    slope, intercept, r, p, se = stats.linregress(df["co2_ppm"], df["global_anomaly_C"])
    xs = np.linspace(df["co2_ppm"].min(), df["co2_ppm"].max(), 50)
    axes[0].plot(xs, intercept + slope * xs, "k--", lw=1.5)
    cbar = plt.colorbar(axes[0].collections[0], ax=axes[0])
    cbar.set_label("Year")
    axes[0].set_xlabel("CO2 (ppm)")
    axes[0].set_ylabel("Temperature anomaly (C)")
    axes[0].set_title(f"Levels: r={level_corr['r']:.3f}\n"
                       f"(corrected p={level_corr['p_autocorr_corrected']:.1e})")

    axes[1].scatter(d_co2, d_temp, color="#c0392b", s=25)
    slope2, intercept2, r2, p2, se2 = stats.linregress(d_co2, d_temp)
    xs2 = np.linspace(d_co2.min(), d_co2.max(), 50)
    axes[1].plot(xs2, intercept2 + slope2 * xs2, "k--", lw=1.5)
    axes[1].set_xlabel("Year-over-year change in CO2 (ppm/yr)")
    axes[1].set_ylabel("Year-over-year change in temperature (C/yr)")
    axes[1].set_title(f"First differences: r={diff_corr['r']:.3f}\n"
                       f"(corrected p={diff_corr['p_autocorr_corrected']:.1e})")
    plt.tight_layout()
    plt.savefig(f"{FIG_DIR}/fig4_co2_temperature_scatter.png")
    plt.close()

    print("\nFigures saved to outputs/figures/")


if __name__ == "__main__":
    main()
