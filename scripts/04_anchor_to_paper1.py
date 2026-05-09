"""
04_anchor_to_paper1.py
----------------------
Combine the geometric life-ratio MC samples (rho) with paper #1's calibrated
concentric fatigue lives (T_conc) to produce absolute-life distributions.

Scheme
------
For each (r_d/R, S-N curve), paper #1 provides T_conc as a single number
(from rainflow + Miner + Weibull integration of the full stress history).

Our MC produced rho_i = T_ecc_i / T_conc, which is geometry-only and is exact
in the multiplicative-stress-scaling regime regardless of the wind/Weibull
integration details. Hence

    T_ecc_i = T_conc * rho_i

is the absolute-life distribution under the same Weibull regime as paper #1.

Source of T_conc
----------------
  data/inherited/parameters/final_results_table.csv  Table 2
  (H = 8 m, DBH = 20 cm, ginkgo, Seoul Weibull k=1.6, c=4.0)

  decay r_d/R :  0.00    0.30    0.50    0.60    0.70    0.75    0.80
  Mean   S-N : 1942.7  1569.8  336.4   62.0    6.9     1.7     0.4
  Lower  S-N :   59.7    53.2   25.5    9.5    1.7     0.5     0.1

Outputs
-------
  data/generated/mc_absolute_lives.csv.gz
  output/tables/absolute_life_summary.csv
  output/figures/fig_absolute_distribution_rd07.png
  output/figures/fig_absolute_percentile.png
  output/figures/fig_inspection_threshold.png
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "data" / "generated"
TBL_DIR  = ROOT / "output" / "tables"
FIG_DIR  = ROOT / "output" / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)
TBL_DIR.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({
    "font.family":      "Times New Roman",
    "axes.titlesize":   11,
    "axes.labelsize":   10,
    "legend.fontsize":  9,
    "xtick.labelsize":  9,
    "ytick.labelsize":  9,
    "axes.spines.top":   False,
    "axes.spines.right": False,
})


# ─────────────────────────────────────────────────────────────────────
# Paper #1 anchor (Seoul, H=8m, DBH=20cm, ginkgo)
# ─────────────────────────────────────────────────────────────────────
T_CONC_PAPER1 = {
    # r_d/R : (T_mean_yr, T_lower_yr)
    0.00: (1942.7, 59.7),
    0.30: (1569.8, 53.2),
    0.50: ( 336.4, 25.5),
    0.60: (  62.0,  9.5),
    0.70: (   6.9,  1.7),
    0.75: (   1.7,  0.5),
    0.80: (   0.4,  0.1),
}


# ─────────────────────────────────────────────────────────────────────
# Load MC rho samples (we have decay levels 0.3, 0.5, 0.6, 0.7, 0.8)
# ─────────────────────────────────────────────────────────────────────
df_lower = pd.read_csv(DATA_DIR / "mc_baseline_lower.csv.gz")
df_mean  = pd.read_csv(DATA_DIR / "mc_baseline_mean.csv.gz")

# Restrict to isotropic wind (Phase B showed wind dist collapses out)
df_lower = df_lower[df_lower["wind"] == "isotropic"].copy()
df_mean  = df_mean[df_mean["wind"] == "isotropic"].copy()


# ─────────────────────────────────────────────────────────────────────
# Build absolute-life arrays
# ─────────────────────────────────────────────────────────────────────
rows = []
abs_arrays = {}                          # (sn, dr) -> ndarray of years

for sn_name, df in [("lower", df_lower), ("mean", df_mean)]:
    for dr in sorted(df["decay_ratio"].unique()):
        if dr not in T_CONC_PAPER1:
            continue
        rho = df[df["decay_ratio"] == dr]["rho"].values
        T_mean_anchor, T_lower_anchor = T_CONC_PAPER1[dr]
        T_anchor = T_lower_anchor if sn_name == "lower" else T_mean_anchor
        T_arr = T_anchor * rho
        abs_arrays[(sn_name, dr)] = T_arr
        rows.append(dict(
            sn=sn_name, decay_ratio=dr,
            T_concentric_anchor=T_anchor,
            T_p5  = float(np.percentile(T_arr, 5)),
            T_p25 = float(np.percentile(T_arr, 25)),
            T_p50 = float(np.percentile(T_arr, 50)),
            T_p75 = float(np.percentile(T_arr, 75)),
            T_p95 = float(np.percentile(T_arr, 95)),
            frac_below_5yr  = float((T_arr <  5.0).mean()),
            frac_below_10yr = float((T_arr < 10.0).mean()),
            frac_below_30yr = float((T_arr < 30.0).mean()),
            frac_below_100yr= float((T_arr <100.0).mean()),
        ))

summary = pd.DataFrame(rows)
out_summary = TBL_DIR / "absolute_life_summary.csv"
summary.to_csv(out_summary, index=False, float_format="%.4f")
print(f"Summary saved: {out_summary}")

# Long-format per-sample CSV
long_rows = []
for (sn, dr), arr in abs_arrays.items():
    for v in arr:
        long_rows.append(dict(sn=sn, decay_ratio=dr, T_yr=v))
df_long = pd.DataFrame(long_rows)
out_long = DATA_DIR / "mc_absolute_lives.csv.gz"
df_long.to_csv(out_long, index=False, float_format="%.6f", compression="gzip")
print(f"Per-sample saved: {out_long}")


# ─────────────────────────────────────────────────────────────────────
# Headline table for paper
# ─────────────────────────────────────────────────────────────────────
print("\n--- Absolute fatigue life (years), Seoul, H=8m DBH=20cm ginkgo ---")
print("\nLower-bound S-N:")
sub = summary[summary["sn"] == "lower"][[
    "decay_ratio", "T_concentric_anchor",
    "T_p5", "T_p25", "T_p50", "T_p75", "T_p95",
    "frac_below_5yr", "frac_below_30yr"]]
print(sub.to_string(index=False, float_format=lambda x: f"{x:8.3f}"))

print("\nMean S-N:")
sub = summary[summary["sn"] == "mean"][[
    "decay_ratio", "T_concentric_anchor",
    "T_p5", "T_p25", "T_p50", "T_p75", "T_p95",
    "frac_below_5yr", "frac_below_30yr"]]
print(sub.to_string(index=False, float_format=lambda x: f"{x:9.3f}"))


# ─────────────────────────────────────────────────────────────────────
# Figure 1 — Absolute life distribution at the Mattheck threshold
# ─────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
DR_FOCUS = 0.7
for ax, sn_name, title in [
    (axes[0], "lower", "Lower-bound S-N"),
    (axes[1], "mean",  "Mean S-N"),
]:
    arr = abs_arrays[(sn_name, DR_FOCUS)]
    T_anchor = T_CONC_PAPER1[DR_FOCUS][0 if sn_name == "mean" else 1]
    ax.hist(arr, bins=80, color="tab:blue", alpha=0.6,
            density=True, edgecolor="white", linewidth=0.4,
            label="eccentric MC (n=10,000)")
    ax.axvline(T_anchor, color="k", lw=2,
               label=f"paper #1 concentric = {T_anchor:.2f} yr")
    p50 = np.percentile(arr, 50)
    p5  = np.percentile(arr,  5)
    ax.axvline(p50, color="tab:red", lw=1.5, ls="--",
               label=f"median = {p50:.3f} yr")
    ax.axvline(p5,  color="tab:orange", lw=1.5, ls=":",
               label=f"5th = {p5:.3f} yr")
    ax.set_xlabel("absolute fatigue life [yr]")
    ax.set_ylabel("density")
    ax.set_title(f"{title}  --  $r_d/R = {DR_FOCUS}$")
    ax.set_xlim(0, T_anchor * 1.15)
    ax.legend(fontsize=8, loc="upper right")
    ax.grid(alpha=0.3)
fig.suptitle(rf"Absolute fatigue life distribution at $r_d/R = {DR_FOCUS}$"
             "\n(Seoul, H=8 m, DBH=20 cm, ginkgo)", fontsize=12)
fig.tight_layout()
out1 = FIG_DIR / "fig_absolute_distribution_rd07.png"
fig.savefig(out1, dpi=200, bbox_inches="tight")
print(f"\nSaved: {out1}")
plt.close(fig)


# ─────────────────────────────────────────────────────────────────────
# Figure 2 — Percentile bands of absolute years vs r_d/R
# ─────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
for ax, sn_name, title in [
    (axes[0], "lower", "Lower-bound S-N"),
    (axes[1], "mean",  "Mean S-N"),
]:
    drs = sorted([dr for sn, dr in abs_arrays if sn == sn_name])
    p5  = np.array([np.percentile(abs_arrays[(sn_name, dr)],  5) for dr in drs])
    p25 = np.array([np.percentile(abs_arrays[(sn_name, dr)], 25) for dr in drs])
    p50 = np.array([np.percentile(abs_arrays[(sn_name, dr)], 50) for dr in drs])
    p75 = np.array([np.percentile(abs_arrays[(sn_name, dr)], 75) for dr in drs])
    p95 = np.array([np.percentile(abs_arrays[(sn_name, dr)], 95) for dr in drs])
    T_anchor = np.array([T_CONC_PAPER1[dr][0 if sn_name == "mean" else 1]
                         for dr in drs])

    ax.fill_between(drs, p5, p95, alpha=0.20, color="tab:blue", label="5-95%")
    ax.fill_between(drs, p25, p75, alpha=0.35, color="tab:blue", label="25-75%")
    ax.plot(drs, p50, "o-", color="tab:blue", lw=2, ms=6, label="median (eccentric)")
    ax.plot(drs, T_anchor, "s--", color="k", lw=1.5, ms=5,
            label="concentric (paper #1)")
    ax.axvline(0.7, color="red", lw=1, ls="--", alpha=0.5,
               label="Mattheck threshold")
    # inspection threshold lines
    ax.axhline(5,  color="orange",  lw=1, ls=":", alpha=0.7)
    ax.text(0.31, 5*1.1,  "5-yr inspection cycle",  fontsize=8, color="orange")
    ax.axhline(30, color="green",   lw=1, ls=":", alpha=0.7)
    ax.text(0.31, 30*1.1, "30-yr (urban tree)",    fontsize=8, color="green")
    ax.set_yscale("log")
    ax.set_xlabel(r"decay ratio  $r_d / R$")
    ax.set_ylabel("absolute fatigue life [yr]")
    ax.set_title(title)
    ax.grid(which="both", alpha=0.3)
    ax.legend(loc="upper right", fontsize=8)
fig.suptitle("Absolute fatigue life: eccentric MC vs concentric (paper #1)\n"
             "(Seoul, H=8 m, DBH=20 cm, ginkgo)", fontsize=12)
fig.tight_layout()
out2 = FIG_DIR / "fig_absolute_percentile.png"
fig.savefig(out2, dpi=200, bbox_inches="tight")
print(f"Saved: {out2}")
plt.close(fig)


# ─────────────────────────────────────────────────────────────────────
# Figure 3 — Inspection threshold (fraction below T_threshold)
# ─────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(8, 5))
T_THRESHOLDS = [1, 5, 10, 30, 100]
COLORS = plt.cm.viridis(np.linspace(0.1, 0.9, len(T_THRESHOLDS)))
for sn_name, ls in [("lower", "-"), ("mean", "--")]:
    drs = sorted([dr for sn, dr in abs_arrays if sn == sn_name])
    for T_th, color in zip(T_THRESHOLDS, COLORS):
        frac = np.array([(abs_arrays[(sn_name, dr)] < T_th).mean()
                         for dr in drs])
        ax.plot(drs, frac, ls=ls, color=color, lw=2,
                marker="o" if ls == "-" else "s", ms=5,
                label=f"{sn_name}, T<{T_th}yr")
ax.axvline(0.7, color="red", lw=1, ls="--", alpha=0.5)
ax.text(0.71, 0.05, "Mattheck", color="red", fontsize=9)
ax.set_xlabel(r"decay ratio  $r_d / R$")
ax.set_ylabel("fraction of trees with life $<$ threshold")
ax.set_title("Probability that an eccentric tree falls below an inspection-cycle threshold")
ax.legend(fontsize=8, ncol=2, loc="upper left")
ax.grid(alpha=0.3)
ax.set_ylim(0, 1.05)
fig.tight_layout()
out3 = FIG_DIR / "fig_inspection_threshold.png"
fig.savefig(out3, dpi=200, bbox_inches="tight")
print(f"Saved: {out3}")
plt.close(fig)
