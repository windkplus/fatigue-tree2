"""
06_policy_matrix.py
-------------------
Table 4 of paper #2: risk-based inspection policy matrix.

For each (region, decay class, S-N curve) cell, derive a recommended
inspection interval

    delta_t  =  T_p5  / safety_factor

where T_p5 is the 5th-percentile of the eccentric-decay fatigue life
distribution, and safety_factor = 2 (i.e. inspect at half the lower-bound).

Output is a long-format CSV plus a heat-mapped figure ready for the paper.

Anchor values
-------------
Paper #1 Table 2 provides T_concentric for Seoul only. For Incheon and Jeju,
absolute lives are scaled from Seoul using the regional damage ratio.

  D(region) / D(Seoul) = ∫p_region(U) U^(2/b) dU  /  ∫p_Seoul(U) U^(2/b) dU

This integral was computed in `weibull_fatigue.py` of paper #1; the output
file `inherited/parameters/weibull_fatigue_results.csv` lists D_annual for
each region (sound tree). We therefore have

  T(region, r_d/R) = T(Seoul, r_d/R) * (D_annual(Seoul) / D_annual(region))
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]

# ─────────────────────────────────────────────────────────────────────
# Anchor values (paper #1 Table 2, Seoul, ginkgo, H=8m, DBH=20cm)
# ─────────────────────────────────────────────────────────────────────
T_CONC_SEOUL = {
    # r_d/R : (T_mean_yr, T_lower_yr)
    0.30: (1569.8, 53.2),
    0.50: ( 336.4, 25.5),
    0.60: (  62.0,  9.5),
    0.70: (   6.9,  1.7),
    0.80: (   0.4,  0.1),
}

# Regional damage ratios (D_annual_region / D_annual_Seoul)
# Read from paper #1 weibull_fatigue_results.csv
weibull = pd.read_csv(ROOT / "data" / "inherited" / "parameters" /
                      "weibull_fatigue_results.csv")
print("Paper #1 regional Weibull results:")
print(weibull.to_string(index=False))

D_SEOUL = float(weibull[weibull["site"].str.startswith("서울")]["D_annual"].iloc[0])
REGIONAL_DAMAGE = {
    "Seoul":   float(weibull[weibull["site"].str.startswith("서울")]["D_annual"].iloc[0]),
    "Incheon": float(weibull[weibull["site"].str.startswith("인천")]["D_annual"].iloc[0]),
    "Jeju":    float(weibull[weibull["site"].str.startswith("제주")]["D_annual"].iloc[0]),
}
REGION_RATIOS = {r: D_SEOUL / d for r, d in REGIONAL_DAMAGE.items()}
# T(region) = T(Seoul) * D_Seoul / D_region   (longer life when damage is smaller)
print("\nRegional life multipliers vs Seoul:")
for r, ratio in REGION_RATIOS.items():
    print(f"  {r:>10}  T_region / T_Seoul = {ratio:.4f}")


# ─────────────────────────────────────────────────────────────────────
# Load eccentric MC rho samples (Phase B baseline, isotropic)
# ─────────────────────────────────────────────────────────────────────
df_low  = pd.read_csv(ROOT / "data" / "generated" / "mc_baseline_lower.csv.gz")
df_low  = df_low[df_low["wind"] == "isotropic"].copy()
df_mean = pd.read_csv(ROOT / "data" / "generated" / "mc_baseline_mean.csv.gz")
df_mean = df_mean[df_mean["wind"] == "isotropic"].copy()


# ─────────────────────────────────────────────────────────────────────
# Build the policy matrix
# ─────────────────────────────────────────────────────────────────────
SAFETY_FACTOR = 2.0
rows = []

for region, life_mult in REGION_RATIOS.items():
    for sn_name, df_rho in [("lower", df_low), ("mean", df_mean)]:
        for dr in T_CONC_SEOUL.keys():
            T_seoul = T_CONC_SEOUL[dr][0 if sn_name == "mean" else 1]
            T_concentric = T_seoul * life_mult
            rho = df_rho[df_rho["decay_ratio"] == dr]["rho"].values
            if len(rho) == 0:
                continue
            T_arr = T_concentric * rho
            T_p5  = float(np.percentile(T_arr,  5))
            T_p50 = float(np.percentile(T_arr, 50))
            T_p95 = float(np.percentile(T_arr, 95))
            rec_interval = T_p5 / SAFETY_FACTOR

            # human-readable category
            if T_p5 < 0.1:
                category = "REMOVE (immediate)"
            elif T_p5 < 1.0:
                category = "REMOVE (urgent)"
            elif T_p5 < 5.0:
                category = "ANNUAL inspection"
            elif T_p5 < 30.0:
                category = "BI-ANNUAL"
            else:
                category = "STANDARD (5-yr)"

            rows.append(dict(
                region=region, sn=sn_name, decay_ratio=dr,
                T_concentric_yr=T_concentric,
                T_p5_yr=T_p5, T_p50_yr=T_p50, T_p95_yr=T_p95,
                inspection_interval_yr=max(rec_interval, 0.05),
                category=category,
            ))

policy = pd.DataFrame(rows)
out_csv = ROOT / "output" / "tables" / "policy_matrix.csv"
policy.to_csv(out_csv, index=False, float_format="%.4f")
print(f"\nSaved: {out_csv}")


# ─────────────────────────────────────────────────────────────────────
# Pretty-print policy table (lower-bound only — most conservative)
# ─────────────────────────────────────────────────────────────────────
print("\n" + "=" * 78)
print("  TABLE 4: Risk-Based Inspection Policy Matrix")
print("  (Lower-bound S-N, ginkgo, H=8m DBH=20cm, eccentric MC 5th-percentile)")
print("=" * 78)
print(f"{'Region':<10}{'r_d/R':>8}{'T_conc':>10}{'T_p5':>9}"
      f"{'recom Δt':>11}  category")
print("-" * 78)
for _, r in policy[policy["sn"] == "lower"].iterrows():
    print(f"{r['region']:<10}{r['decay_ratio']:>8.2f}"
          f"{r['T_concentric_yr']:>10.2f}{r['T_p5_yr']:>9.3f}"
          f"{r['inspection_interval_yr']:>11.2f}  {r['category']}")


# ─────────────────────────────────────────────────────────────────────
# Visualisation: heatmap of recommended intervals + categorical map
# ─────────────────────────────────────────────────────────────────────
plt.rcParams.update({
    "font.family": "Times New Roman",
    "axes.titlesize": 11, "axes.labelsize": 10,
    "legend.fontsize": 9, "xtick.labelsize": 9, "ytick.labelsize": 9,
})

REGIONS = ["Seoul", "Incheon", "Jeju"]
DECAYS = sorted(T_CONC_SEOUL.keys())

fig, axes = plt.subplots(1, 2, figsize=(13, 5))

# Panel 1: numeric heatmap of recommended interval (lower-bound)
ax = axes[0]
mat_int = np.zeros((len(REGIONS), len(DECAYS)))
for i, region in enumerate(REGIONS):
    for j, dr in enumerate(DECAYS):
        sub = policy[(policy["region"] == region) &
                     (policy["sn"] == "lower") &
                     (policy["decay_ratio"] == dr)]
        if len(sub):
            mat_int[i, j] = sub["inspection_interval_yr"].values[0]

im = ax.imshow(mat_int, aspect="auto", cmap="RdYlGn",
               norm=plt.matplotlib.colors.LogNorm(vmin=0.05, vmax=30.0))
ax.set_xticks(range(len(DECAYS)))
ax.set_xticklabels([f"{d:.1f}" for d in DECAYS])
ax.set_yticks(range(len(REGIONS)))
ax.set_yticklabels(REGIONS)
ax.set_xlabel(r"decay ratio  $r_d / R$")
ax.set_ylabel("region")
ax.set_title(r"Recommended inspection interval $\Delta t$ [yr]"
             "\n(Lower-bound S-N, ginkgo H=8m DBH=20cm)")
for i in range(len(REGIONS)):
    for j in range(len(DECAYS)):
        v = mat_int[i, j]
        if v < 0.1:
            txt = "<0.1"
        elif v < 1:
            txt = f"{v:.2f}"
        elif v < 10:
            txt = f"{v:.1f}"
        else:
            txt = f"{v:.0f}"
        ax.text(j, i, txt, ha="center", va="center", fontsize=9,
                color="white" if v < 0.5 else "black")
plt.colorbar(im, ax=ax, label=r"$\Delta t$ [yr]")

# Panel 2: categorical map
ax = axes[1]
cat_levels = ["REMOVE (immediate)", "REMOVE (urgent)",
              "ANNUAL inspection", "BI-ANNUAL", "STANDARD (5-yr)"]
cat_to_int = {c: i for i, c in enumerate(cat_levels)}
mat_cat = np.zeros((len(REGIONS), len(DECAYS)), dtype=int)
for i, region in enumerate(REGIONS):
    for j, dr in enumerate(DECAYS):
        sub = policy[(policy["region"] == region) &
                     (policy["sn"] == "lower") &
                     (policy["decay_ratio"] == dr)]
        if len(sub):
            mat_cat[i, j] = cat_to_int[sub["category"].values[0]]

cmap = plt.matplotlib.colors.ListedColormap(
    ["#7a0000", "#c43a3a", "#f4a445", "#ffe699", "#9bd99b"])
im2 = ax.imshow(mat_cat, aspect="auto", cmap=cmap, vmin=0, vmax=4)
ax.set_xticks(range(len(DECAYS)))
ax.set_xticklabels([f"{d:.1f}" for d in DECAYS])
ax.set_yticks(range(len(REGIONS)))
ax.set_yticklabels(REGIONS)
ax.set_xlabel(r"decay ratio  $r_d / R$")
ax.set_title("Risk category\n(Lower-bound S-N)")
for i in range(len(REGIONS)):
    for j in range(len(DECAYS)):
        c_int = mat_cat[i, j]
        ax.text(j, i, cat_levels[c_int].split()[0], ha="center", va="center",
                fontsize=9, color="white" if c_int < 2 else "black",
                fontweight="bold")
# legend bar
cbar = plt.colorbar(im2, ax=ax, ticks=range(5))
cbar.ax.set_yticklabels(cat_levels, fontsize=8)

fig.tight_layout()
out_fig = ROOT / "output" / "figures" / "fig_policy_matrix.png"
fig.savefig(out_fig, dpi=200, bbox_inches="tight")
print(f"\nSaved: {out_fig}")
plt.close(fig)
