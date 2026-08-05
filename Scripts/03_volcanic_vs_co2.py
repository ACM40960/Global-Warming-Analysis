"""
03_volcanic_vs_co2.py
----------------------
Tests the alternative hypothesis: "could the rise in CO2 (and hence
temperature) be explained by an increase in volcanic activity over the
industrial era, rather than by human emissions?"

We build two independent annual volcanic-activity indices from the
Smithsonian Global Volcanism Program eruption catalog (1850-2025):
  - eruption COUNT per year (all confirmed eruptions)
  - a VEI-weighted MAGNITUDE INDEX per year: sum(10**VEI) across that
    year's eruptions, which approximates the (highly nonlinear) relation
    between the Volcanic Explosivity Index and erupted mass/gas volume.

We then:
  1. Test whether either volcanic index itself has a significant secular
     (multi-century) trend -- if volcanic activity isn't trending up,
     it cannot explain a monotonic rise in CO2.
  2. Correlate both volcanic indices against atmospheric CO2 over the
     period where the two datasets overlap (1979-2025), again using the
     autocorrelation-corrected significance test.
  3. Plot both series for visual comparison.

Context (not derived from these files, standard geochemical accounting
widely reported by USGS/IPCC): human fossil-fuel and cement emissions
are estimated at roughly 100x the CO2 output of all of Earth's
volcanoes combined, averaged over a year. That figure isn't computed
here -- it's mentioned in the report as independent corroborating
context for the correlation result below.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from utils import (
    RAW_DIR, FIG_DIR, DATA_OUT_DIR, apply_style, ols_trend_with_ci,
    mann_kendall_test, pearson_autocorr_corrected,
)

apply_style()


def load_volcanic():
    df = pd.read_csv(f"{RAW_DIR}/Volcanic_Activity.tsv", sep="\t", header=0, skiprows=[1])
    df = df[["Year", "VEI"]].dropna(subset=["Year"])
    df["Year"] = df["Year"].astype(int)
    # missing VEI entries: treat conservatively as VEI=1 (small, mostly
    # historically under-documented eruptions) rather than dropping them
    df["VEI_filled"] = df["VEI"].fillna(1)

    years = np.arange(1850, 2026)
    counts = df.groupby("Year").size().reindex(years, fill_value=0)
    magnitude = df.groupby("Year")["VEI_filled"].apply(lambda v: np.sum(10.0 ** v)).reindex(years, fill_value=0)

    out = pd.DataFrame({
        "year": years,
        "eruption_count": counts.values,
        "volcanic_magnitude_index": magnitude.values,
    })
    return out


def load_co2():
    df = pd.read_csv(f"{RAW_DIR}/Global_C02_Annual.csv", comment="#")
    df = df.rename(columns={"mean": "co2_ppm", "unc": "co2_unc_ppm"})
    df["year"] = df["year"].astype(int)
    return df[["year", "co2_ppm", "co2_unc_ppm"]].dropna(subset=["co2_ppm"])


def main():
    volc = load_volcanic()
    co2 = load_co2()
    volc.to_csv(f"{DATA_OUT_DIR}/volcanic_annual_index.csv", index=False)

    # --- 1. secular trend in volcanic activity itself ---
    trend_count = ols_trend_with_ci(volc["year"], volc["eruption_count"])
    mk_count = mann_kendall_test(volc["eruption_count"])
    trend_mag = ols_trend_with_ci(volc["year"], volc["volcanic_magnitude_index"])
    mk_mag = mann_kendall_test(volc["volcanic_magnitude_index"])

    print("VOLCANIC ACTIVITY TREND, 1850-2025")
    print(f"  Eruption count/decade: {trend_count['per_decade']:+.3f} "
          f"(R^2={trend_count['r_squared']:.3f}, p={trend_count['p_value']:.2e}); "
          f"Mann-Kendall: {mk_count['trend']} (p={mk_count['p_value']:.2e})")
    print(f"  Magnitude index/decade: {trend_mag['per_decade']:+.2e} "
          f"(R^2={trend_mag['r_squared']:.3f}, p={trend_mag['p_value']:.2e}); "
          f"Mann-Kendall: {mk_mag['trend']} (p={mk_mag['p_value']:.2e})")

    # --- 2. correlation with CO2 over overlap period ---
    merged = pd.merge(volc, co2, on="year", how="inner")
    corr_count = pearson_autocorr_corrected(merged["eruption_count"], merged["co2_ppm"])
    corr_mag = pearson_autocorr_corrected(merged["volcanic_magnitude_index"], merged["co2_ppm"])

    print(f"\nCorrelation with CO2, {merged.year.min()}-{merged.year.max()} (n={len(merged)}):")
    print(f"  eruption count vs CO2:      r={corr_count['r']:.3f}, "
          f"corrected p={corr_count['p_autocorr_corrected']:.3f} (N_eff={corr_count['n_eff']:.1f})")
    print(f"  magnitude index vs CO2:     r={corr_mag['r']:.3f}, "
          f"corrected p={corr_mag['p_autocorr_corrected']:.3f} (N_eff={corr_mag['n_eff']:.1f})")

    with open(f"{DATA_OUT_DIR}/volcanic_co2_stats.txt", "w") as f:
        f.write("VOLCANIC ACTIVITY vs CO2\n" + "=" * 50 + "\n")
        f.write("Secular trend in volcanic activity, 1850-2025:\n")
        f.write(f"  eruption count:   {trend_count['per_decade']:+.3f}/decade, "
                f"R^2={trend_count['r_squared']:.3f}, p={trend_count['p_value']:.2e} "
                f"-> {mk_count['trend']}\n")
        f.write(f"  magnitude index:  {trend_mag['per_decade']:+.2e}/decade, "
                f"R^2={trend_mag['r_squared']:.3f}, p={trend_mag['p_value']:.2e} "
                f"-> {mk_mag['trend']}\n\n")
        f.write(f"Correlation with CO2 ({merged.year.min()}-{merged.year.max()}, n={len(merged)}):\n")
        f.write(f"  eruption count vs CO2:   r={corr_count['r']:.4f}, "
                f"corrected p={corr_count['p_autocorr_corrected']:.4f} "
                f"(N_eff={corr_count['n_eff']:.1f})\n")
        f.write(f"  magnitude idx vs CO2:    r={corr_mag['r']:.4f}, "
                f"corrected p={corr_mag['p_autocorr_corrected']:.4f} "
                f"(N_eff={corr_mag['n_eff']:.1f})\n")

    # --- Figure 5 ---
    fig, axes = plt.subplots(3, 1, figsize=(11, 10), sharex=True)
    axes[0].bar(volc["year"], volc["eruption_count"], color="#7f8c8d", width=1)
    axes[0].set_title("Annual volcanic eruption count, 1850-2025 (Smithsonian GVP catalog)")
    axes[0].set_ylabel("Eruptions/year")

    axes[1].bar(volc["year"], volc["volcanic_magnitude_index"], color="#b03a2e", width=1)
    axes[1].set_title("VEI-weighted volcanic magnitude index (sum of 10^VEI per year)")
    axes[1].set_ylabel("Magnitude index")
    axes[1].set_yscale("log")

    ax2 = axes[2]
    ax2b = ax2.twinx()
    ax2.plot(merged["year"], merged["volcanic_magnitude_index"], color="#b03a2e", label="Volcanic magnitude index")
    ax2b.errorbar(merged["year"], merged["co2_ppm"], yerr=merged["co2_unc_ppm"], color="#7d3c98",
                  fmt="s--", ms=3, capsize=2, label="CO2 (ppm)")
    ax2.set_ylabel("Volcanic magnitude idx", color="#b03a2e")
    ax2b.set_ylabel("CO2 (ppm)", color="#7d3c98")
    ax2.set_xlabel("Year")
    ax2.set_title(f"Overlap period comparison: r={corr_mag['r']:.3f} "
                  f"(corrected p={corr_mag['p_autocorr_corrected']:.3f}, not significant)"
                  if corr_mag['p_autocorr_corrected'] > 0.05 else
                  f"Overlap period comparison: r={corr_mag['r']:.3f} "
                  f"(corrected p={corr_mag['p_autocorr_corrected']:.3f})")
    lines1, labels1 = ax2.get_legend_handles_labels()
    lines2, labels2 = ax2b.get_legend_handles_labels()
    ax2.legend(lines1 + lines2, labels1 + labels2, loc="upper left")

    plt.tight_layout()
    plt.savefig(f"{FIG_DIR}/fig5_volcanic_vs_co2.png")
    plt.close()
    print("\nFigure saved to outputs/figures/fig5_volcanic_vs_co2.png")


if __name__ == "__main__":
    main()
