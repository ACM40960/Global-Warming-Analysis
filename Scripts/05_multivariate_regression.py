"""
05_multivariate_regression.py
-------------------------------
Puts all three candidate drivers in a single multiple regression against
global temperature, over the period where all four series overlap
(1979-2025, limited by the CO2 record):

    temperature ~ CO2 + TSI + volcanic_magnitude_index

We standardize every predictor (z-score) first, so the fitted
coefficients are directly comparable "effect sizes" regardless of each
variable's native units. We use Newey-West (HAC) standard errors instead
of ordinary OLS standard errors, because climate time series are
strongly autocorrelated and OLS would otherwise understate the
uncertainty on each coefficient.

This does not, by itself, "prove" causality (no regression can) -- but
it directly answers the question the project asked: of the three
candidate explanations for the observed warming, which one actually
carries statistical weight once the others are controlled for?
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import statsmodels.api as sm

from utils import RAW_DIR, FIG_DIR, DATA_OUT_DIR, apply_style

apply_style()


def main():
    temp = pd.read_csv(f"{DATA_OUT_DIR}/global_temperature_drift.csv")
    co2 = pd.read_csv(f"{RAW_DIR}/Global_C02_Annual.csv", comment="#").rename(
        columns={"mean": "co2_ppm"})[["year", "co2_ppm"]]
    co2["year"] = co2["year"].astype(int)
    tsi = pd.read_csv(f"{DATA_OUT_DIR}/solar_tsi_annual.csv")
    volc = pd.read_csv(f"{DATA_OUT_DIR}/volcanic_annual_index.csv")

    df = temp.merge(co2, on="year").merge(tsi, on="year").merge(volc, on="year")
    df = df.dropna(subset=["global_anomaly_C", "co2_ppm", "tsi_wm2", "volcanic_magnitude_index"])
    print(f"Multivariate regression sample: {df.year.min()}-{df.year.max()} (n={len(df)})")

    def z(s):
        return (s - s.mean()) / s.std()

    X = pd.DataFrame({
        "CO2": z(df["co2_ppm"]),
        "Solar_TSI": z(df["tsi_wm2"]),
        "Volcanic_magnitude": z(np.log1p(df["volcanic_magnitude_index"])),
    })
    X = sm.add_constant(X)
    y = df["global_anomaly_C"].values

    model = sm.OLS(y, X).fit(cov_type="HAC", cov_kwds={"maxlags": 3})
    print(model.summary())

    with open(f"{DATA_OUT_DIR}/multivariate_regression_summary.txt", "w") as f:
        f.write("MULTIVARIATE REGRESSION: temperature ~ CO2 + Solar TSI + Volcanic magnitude\n")
        f.write(f"(standardized predictors, Newey-West HAC standard errors, "
                f"period {df.year.min()}-{df.year.max()}, n={len(df)})\n")
        f.write("=" * 70 + "\n")
        f.write(str(model.summary()))

    # --- Figure 8: standardized coefficients with 95% CI ---
    params = model.params.drop("const")
    conf = model.conf_int().drop("const")
    conf.columns = ["ci_low", "ci_high"]
    err = np.vstack([params - conf["ci_low"], conf["ci_high"] - params])

    fig, ax = plt.subplots(figsize=(8, 5))
    colors = ["#7d3c98" if p == "CO2" else "#e67e22" if p == "Solar_TSI" else "#b03a2e"
              for p in params.index]
    ax.bar(params.index, params.values, yerr=err, capsize=6, color=colors, alpha=0.85)
    ax.axhline(0, color="black", lw=0.8)
    ax.set_ylabel("Standardized regression coefficient\n(effect on temperature anomaly, per 1 SD change)")
    ax.set_title(f"Relative contribution of each candidate driver to warming\n"
                 f"({df.year.min()}-{df.year.max()}, multivariate OLS, HAC standard errors)")
    for i, p in enumerate(params.index):
        pval = model.pvalues[p]
        stars = "***" if pval < 0.001 else "**" if pval < 0.01 else "*" if pval < 0.05 else "n.s."
        ax.text(i, params.values[i] + np.sign(params.values[i]) * (err[1 if params.values[i] >= 0 else 0][i] + 0.03),
                stars, ha="center", fontsize=12, fontweight="bold")
    plt.tight_layout()
    plt.savefig(f"{FIG_DIR}/fig8_multivariate_regression_coefficients.png")
    plt.close()

    print(f"\nR^2 = {model.rsquared:.3f}, Adj R^2 = {model.rsquared_adj:.3f}")
    print("Standardized coefficients (effect size per 1 SD change in predictor):")
    for p in params.index:
        print(f"  {p:20s}: {params[p]:+.4f}  [{conf.loc[p,'ci_low']:+.4f}, {conf.loc[p,'ci_high']:+.4f}]  "
              f"p={model.pvalues[p]:.3e}")

    print("\nFigure saved to outputs/figures/fig8_multivariate_regression_coefficients.png")


if __name__ == "__main__":
    main()
