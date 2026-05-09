"""
03_visualize_mc.py
------------------
Publication-quality figures for paper #2 Phase B Monte Carlo results.

Inputs
------
  data/generated/mc_baseline_lower.csv.gz   (b = 0.125, 30 scenarios x 10,000 rho)
  data/generated/mc_baseline_mean.csv.gz    (b = 0.10)
  output/tables/mc_baseline_summary.csv

Outputs
-------
  output/figures/fig_mc_cdf.png          : life-ratio CDF for each (r_d/R, wind)
  output/figures/fig_mc_percentile.png   : percentile bands vs r_d/R
  output/figures/fig_mc_heatmap.png      : life-shortening factor heatmap
  output/figures/fig_mc_concentric_vs_ecc.png : key narrative figure
  output/tables/mc_policy_table.csv      : concise table for the paper
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
# Load
# ─────────────────────────────────────────────────────────────────────
df_lower = pd.read_csv(DATA_DIR / "mc_baseline_lower.csv.gz")
df_mean  = pd.read_csv(DATA_DIR / "mc_baseline_mean.csv.gz")
summary  = pd.read_csv(TBL_DIR / "mc_baseline_summary.csv")

print(f"Loaded lower-bound : {len(df_lower):,d} rows")
print(f"Loaded mean        : {len(df_mean):,d} rows")
print(f"Summary scenarios  : {len(summary)} rows\n")


# ─────────────────────────────────────────────────────────────────────
# Figure 1 — CDF panels (one per r_d/R) for the LOWER-BOUND S-N
# ─────────────────────────────────────────────────────────────────────
DECAYS = sorted(df_lower["decay_ratio"].unique())
WINDS = ["isotropic", "seoul", "jeju"]
COLORS = {"isotropic": "tab:blue", "seoul": "tab:orange", "jeju": "tab:red"}

fig, axes = plt.subplots(1, len(DECAYS), figsize=(4 * len(DECAYS), 4),
                         sharey=True)
for ax, dr in zip(axes, DECAYS):
    for w in WINDS:
        rhos = np.sort(df_lower[(df_lower["decay_ratio"] == dr) &
                                (df_lower["wind"] == w)]["rho"].values)
        cdf = np.arange(1, len(rhos) + 1) / len(rhos)
        ax.plot(rhos, cdf, lw=1.5, color=COLORS[w], label=w)
    ax.axvline(1.0, color="k", lw=0.8, ls=":", label="concentric (rho=1)")
    ax.set_xlabel(r"life ratio $\rho = T_{\rm ecc}\,/\,T_{\rm conc}$")
    ax.set_title(rf"$r_d/R = {dr:.1f}$")
    ax.set_xlim(0, 1.1)
    ax.grid(alpha=0.3)
    if ax is axes[0]:
        ax.set_ylabel("CDF")
        ax.legend(loc="lower right")
fig.suptitle("Fatigue life ratio CDF (lower-bound S-N, b=0.125)", fontsize=12)
fig.tight_layout()
out1 = FIG_DIR / "fig_mc_cdf.png"
fig.savefig(out1, dpi=200, bbox_inches="tight")
print(f"Saved: {out1}")
plt.close(fig)


# ─────────────────────────────────────────────────────────────────────
# Figure 2 — Percentile bands vs r_d/R (lower vs mean side-by-side)
# ─────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), sharey=True)

for ax, df, title in [
    (axes[0], df_lower, "Lower-bound S-N (b=0.125)"),
    (axes[1], df_mean,  "Mean S-N (b=0.10)"),
]:
    # use isotropic only - wind dists are virtually identical
    sub = df[df["wind"] == "isotropic"]
    p5  = sub.groupby("decay_ratio")["rho"].quantile(0.05).values
    p25 = sub.groupby("decay_ratio")["rho"].quantile(0.25).values
    p50 = sub.groupby("decay_ratio")["rho"].quantile(0.50).values
    p75 = sub.groupby("decay_ratio")["rho"].quantile(0.75).values
    p95 = sub.groupby("decay_ratio")["rho"].quantile(0.95).values
    drs = sorted(sub["decay_ratio"].unique())

    ax.fill_between(drs, p5,  p95, alpha=0.20, color="tab:blue", label="5-95%")
    ax.fill_between(drs, p25, p75, alpha=0.35, color="tab:blue", label="25-75%")
    ax.plot(drs, p50, "o-", color="tab:blue", lw=2, ms=6, label="median")
    ax.axhline(1.0, color="k", lw=0.8, ls=":", label="concentric")
    ax.axvline(0.7, color="red", lw=1, ls="--", alpha=0.5,
               label="Mattheck threshold")
    ax.set_xlabel(r"decay ratio  $r_d / R$")
    ax.set_title(title)
    ax.set_yscale("log")
    ax.grid(which="both", alpha=0.3)
    ax.legend(loc="lower left", fontsize=8)
axes[0].set_ylabel(r"life ratio  $\rho$")
fig.suptitle("Eccentric-decay fatigue life relative to concentric (Mattheck) baseline",
             fontsize=12)
fig.tight_layout()
out2 = FIG_DIR / "fig_mc_percentile.png"
fig.savefig(out2, dpi=200, bbox_inches="tight")
print(f"Saved: {out2}")
plt.close(fig)


# ─────────────────────────────────────────────────────────────────────
# Figure 3 — Heatmap: 1/rho (life-shortening factor) by (r_d/R, percentile)
# ─────────────────────────────────────────────────────────────────────
percentile_names = ["p5", "p25", "p50", "p75", "p95"]
percentile_labels = ["5%", "25%", "50% (median)", "75%", "95%"]

fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
for ax, b_name, title in [
    (axes[0], "lower", "Lower-bound S-N (b=0.125)"),
    (axes[1], "mean",  "Mean S-N (b=0.10)"),
]:
    sub = summary[(summary["b_name"] == b_name) & (summary["wind"] == "isotropic")]
    mat = np.zeros((len(percentile_names), len(DECAYS)))
    for i, p in enumerate(percentile_names):
        for j, dr in enumerate(DECAYS):
            row = sub[sub["decay_ratio"] == dr]
            mat[i, j] = 1.0 / row[p].values[0]    # life-shortening factor
    im = ax.imshow(mat, aspect="auto", cmap="YlOrRd",
                   norm=plt.matplotlib.colors.LogNorm(vmin=1.0, vmax=mat.max()))
    ax.set_xticks(range(len(DECAYS)))
    ax.set_xticklabels([f"{d:.1f}" for d in DECAYS])
    ax.set_yticks(range(len(percentile_names)))
    ax.set_yticklabels(percentile_labels)
    ax.set_xlabel(r"decay ratio  $r_d / R$")
    ax.set_ylabel("sample percentile")
    ax.set_title(title)
    for i in range(len(percentile_names)):
        for j in range(len(DECAYS)):
            v = mat[i, j]
            txt = f"{v:.0f}x" if v >= 10 else f"{v:.1f}x"
            ax.text(j, i, txt, ha="center", va="center", fontsize=8,
                    color="white" if v > 10 else "black")
    plt.colorbar(im, ax=ax, label=r"$1/\rho$  (life-shortening factor)")
fig.suptitle("Fatigue life shortening relative to concentric assumption",
             fontsize=12)
fig.tight_layout()
out3 = FIG_DIR / "fig_mc_heatmap.png"
fig.savefig(out3, dpi=200, bbox_inches="tight")
print(f"Saved: {out3}")
plt.close(fig)


# ─────────────────────────────────────────────────────────────────────
# Figure 4 — Narrative figure: concentric vs eccentric distribution
# ─────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
DR_FOCUS = 0.7
for ax, df, title in [
    (axes[0], df_lower, "Lower-bound S-N (b=0.125)"),
    (axes[1], df_mean,  "Mean S-N (b=0.10)"),
]:
    sub = df[(df["decay_ratio"] == DR_FOCUS) & (df["wind"] == "isotropic")]
    ax.hist(sub["rho"], bins=80, color="tab:blue", alpha=0.6,
            density=True, edgecolor="white", linewidth=0.4,
            label="eccentric MC (n=10,000)")
    ax.axvline(1.0, color="k", lw=2, label="concentric (Mattheck)")
    p5  = sub["rho"].quantile(0.05)
    p50 = sub["rho"].quantile(0.50)
    ax.axvline(p50, color="tab:red", lw=1.5, ls="--",
               label=f"median = {p50:.3f}")
    ax.axvline(p5, color="tab:orange", lw=1.5, ls=":",
               label=f"5th = {p5:.3f}")
    ax.set_xlabel(r"life ratio  $\rho = T_{\rm ecc} / T_{\rm conc}$")
    ax.set_ylabel("density")
    ax.set_title(f"{title}  --  $r_d/R = {DR_FOCUS}$")
    ax.set_xlim(0, 1.1)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
fig.suptitle("Distribution of fatigue life ratio at the Mattheck threshold "
             rf"($r_d/R = {DR_FOCUS}$)",
             fontsize=12)
fig.tight_layout()
out4 = FIG_DIR / "fig_mc_concentric_vs_ecc.png"
fig.savefig(out4, dpi=200, bbox_inches="tight")
print(f"Saved: {out4}")
plt.close(fig)


# ─────────────────────────────────────────────────────────────────────
# Compact policy table for the paper
# ─────────────────────────────────────────────────────────────────────
rows = []
for b_name in ("lower", "mean"):
    for dr in DECAYS:
        for w in WINDS:
            r = summary[(summary["b_name"] == b_name) &
                        (summary["wind"] == w) &
                        (summary["decay_ratio"] == dr)].iloc[0]
            rows.append(dict(
                S_N=b_name, wind=w, decay_ratio=dr,
                rho_p5=r["p5"], rho_p50=r["p50"], rho_p95=r["p95"],
                shorten_p5=1.0 / r["p5"],
                shorten_p50=1.0 / r["p50"],
                shorten_p95=1.0 / r["p95"],
            ))
policy = pd.DataFrame(rows)
out_policy = TBL_DIR / "mc_policy_table.csv"
policy.to_csv(out_policy, index=False, float_format="%.4f")
print(f"\nPolicy table saved: {out_policy}")

# Print the headline table
print("\n--- Headline: lower-bound S-N, isotropic wind ---")
sub = policy[(policy["S_N"] == "lower") & (policy["wind"] == "isotropic")]
print(sub[["decay_ratio", "rho_p5", "rho_p50", "rho_p95",
           "shorten_p5", "shorten_p50", "shorten_p95"]]
      .to_string(index=False, float_format=lambda x: f"{x:.3f}"))
