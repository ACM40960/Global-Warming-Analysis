"""
run_all.py
----------
Runs the full analysis pipeline end to end, in the correct order
(each step reads outputs written by the previous ones).

Usage:
    cd scripts
    python3 run_all.py
"""
import subprocess
import sys

STEPS = [
    "01_process_temperature.py",
    "02_co2_temperature_correlation.py",
    "03_volcanic_vs_co2.py",
    "04_solar_vs_co2_temperature.py",
    "05_multivariate_regression.py",
]

if __name__ == "__main__":
    for step in STEPS:
        print("\n" + "=" * 70)
        print(f"RUNNING {step}")
        print("=" * 70)
        result = subprocess.run([sys.executable, step])
        if result.returncode != 0:
            print(f"\n!! {step} failed, stopping pipeline.")
            sys.exit(result.returncode)
    print("\n" + "=" * 70)
    print("PIPELINE COMPLETE. See ../outputs/data/ and ../outputs/figures/")
    print("=" * 70)
