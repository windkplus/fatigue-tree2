"""
02_monte_carlo_baseline.py
--------------------------
First Monte Carlo run for paper #2.

For each (decay_ratio, wind_distribution) combination, draw N=10,000 eccentric
circular cavity geometries from the Phase B priors and compute the per-sample
fatigue life ratio rho = T_eccentric / T_concentric.

Decisions confirmed in Phase B:
  - Wind direction:  isotropic + Seoul/Jeju empirical rose (3 scenarios)
  - Cavity offset:   isotropic (Beta(2, 2) on e_norm; Uniform on angle)
  - Joint:           independent (cavity orientation independent of wind)
  - S-N exponent:    1/b = 8 (lower-bound, b = 0.125) and 1/b = 10 (mean, b = 0.10)

Outputs
-------
  data/generated/mc_baseline_lower.parquet   (b = 0.125, lower-bound S-N)
  data/generated/mc_baseline_mean.parquet    (b = 0.10,  mean S-N)
  output/tables/mc_baseline_summary.csv      (percentiles per scenario)
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from eccentric_decay import (  # noqa: E402
    GeometryPriors,
    reference_distributions,
    sample_eccentric_circular,
    life_ratios,
    summarize,
)


# ─────────────────────────────────────────────────────────────────────
# Configuration
# ─────────────────────────────────────────────────────────────────────
SEED = 42
N_SAMPLES = 10_000
N_THETA = 360
DECAY_RATIOS = [0.3, 0.5, 0.6, 0.7, 0.8]
SN_EXPONENTS = {           # name -> b
    "lower": 0.125,        # lower-bound S-N (paper #1)
    "mean":  0.10,         # mean S-N (paper #1)
}
PRIORS = GeometryPriors(ecc_norm_alpha=2.0, ecc_norm_beta=2.0)

# Output paths
DATA_DIR = ROOT / "data" / "generated"
TBL_DIR = ROOT / "output" / "tables"
DATA_DIR.mkdir(parents=True, exist_ok=True)
TBL_DIR.mkdir(parents=True, exist_ok=True)


# ─────────────────────────────────────────────────────────────────────
# Build wind direction distributions ONCE (shared across all decay levels)
# ─────────────────────────────────────────────────────────────────────
wind_dists = reference_distributions(n_theta=N_THETA)
for name, (theta, w) in wind_dists.items():
    print(f"  {name:>10}  weights sum = {w.sum():.6f}  "
          f"mode at theta = {np.degrees(theta[np.argmax(w)]):.0f} deg")


# ─────────────────────────────────────────────────────────────────────
# Run Monte Carlo
# ─────────────────────────────────────────────────────────────────────
all_records = []                  # for the long-format summary CSV
all_arrays = {}                   # (b_name, wind_name, dr) -> rho array

for b_name, b in SN_EXPONENTS.items():
    print(f"\n=== S-N exponent: 1/b = {1.0/b:.0f}  (b = {b}, '{b_name}') ===")
    for dr in DECAY_RATIOS:
        sections, geom = sample_eccentric_circular(
            decay_ratio=dr, n_samples=N_SAMPLES, R=1.0,
            priors=PRIORS, seed=SEED + int(dr * 100))
        for w_name, (theta, weights) in wind_dists.items():
            rhos = life_ratios(sections, dr, theta, weights, b=b)
            stats = summarize(rhos)
            all_arrays[(b_name, w_name, dr)] = rhos
            row = dict(b_name=b_name, b=b, wind=w_name, decay_ratio=dr, **stats)
            all_records.append(row)
            print(f"  decay={dr:.2f}  wind={w_name:>10}  "
                  f"median rho = {stats['p50']:.3f}   "
                  f"5th = {stats['p5']:.3f}   95th = {stats['p95']:.3f}")


# ─────────────────────────────────────────────────────────────────────
# Persist results
# ─────────────────────────────────────────────────────────────────────
summary_df = pd.DataFrame(all_records)
summary_csv = TBL_DIR / "mc_baseline_summary.csv"
summary_df.to_csv(summary_csv, index=False, float_format="%.6f")
print(f"\nSummary CSV saved: {summary_csv}")

# Per-sample arrays as long-format parquet (one file per S-N exponent)
for b_name in SN_EXPONENTS:
    rows = []
    for (bb, w_name, dr), rhos in all_arrays.items():
        if bb != b_name:
            continue
        for r in rhos:
            rows.append(dict(wind=w_name, decay_ratio=dr, rho=r))
    df = pd.DataFrame(rows)
    out = DATA_DIR / f"mc_baseline_{b_name}.parquet"
    try:
        df.to_parquet(out, index=False)
        print(f"Per-sample data saved: {out}")
    except Exception as exc:
        # parquet engines may be missing - fall back to compressed CSV
        out_csv = out.with_suffix(".csv.gz")
        df.to_csv(out_csv, index=False, float_format="%.6f", compression="gzip")
        print(f"Per-sample data saved (CSV fallback): {out_csv}  ({exc})")


# ─────────────────────────────────────────────────────────────────────
# Compact summary table for the paper (lower-bound only, isotropic + roses)
# ─────────────────────────────────────────────────────────────────────
print("\n--- Lower-bound S-N (b = 0.125) life ratio rho = T_ecc / T_conc ---")
sub = summary_df[summary_df["b_name"] == "lower"].copy()
piv = sub.pivot_table(index="decay_ratio", columns="wind",
                      values=["p5", "p50", "p95"])
print(piv.round(3).to_string())
