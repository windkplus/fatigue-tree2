"""
05_sensitivity.py
-----------------
Phase C2 sensitivity analysis.

Three sub-analyses:

  (1) Beta-prior sensitivity on cavity normalised eccentricity.
      Tests whether the headline conclusions of Phase B are robust to the
      choice of geometric prior. Compares Beta(1,1) [uniform], Beta(2,2)
      [baseline], Beta(2,5) [low-ecc skew], Beta(5,2) [high-ecc skew].

  (2) Elliptical cavity aspect-ratio sweep.
      Adds a non-circular cavity dimension at fixed cavity area. Tests
      whether elongated decay pockets are more or less severe than the
      circular baseline at the same area-equivalent r_d/R.

  (3) Variance decomposition.
      Conditional-variance Sobol-style indices estimate how much of the
      total variance in rho is explained by ecc_norm vs. ecc_angle vs.
      aspect (when included).

Outputs
-------
  output/figures/fig_sens_prior.png
  output/figures/fig_sens_aspect.png
  output/figures/fig_sens_variance.png
  output/tables/sensitivity_summary.csv
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from eccentric_decay import (  # noqa: E402
    DecaySection,
    GeometryPriors,
    isotropic_pdf,
    saf,
    sample_eccentric_circular,
    sample_elliptical,
    life_ratios,
    summarize,
)


# ─────────────────────────────────────────────────────────────────────
# Configuration
# ─────────────────────────────────────────────────────────────────────
N_SAMPLES   = 5000
N_THETA     = 360
B_LOWER     = 0.125
DECAY_FOCUS = [0.5, 0.7]

PRIOR_VARIANTS = {
    "Beta(1,1) uniform":    (1.0, 1.0),
    "Beta(2,2) baseline":   (2.0, 2.0),
    "Beta(2,5) low-ecc":    (2.0, 5.0),
    "Beta(5,2) high-ecc":   (5.0, 2.0),
}

ASPECT_VARIANTS = {        # b / a (i.e., minor / major). Aspect_ratio == 1 -> circular.
    "circular":     1.0,
    "ellipse 1:2": 0.5,
    "ellipse 1:3": 1.0 / 3.0,
}

OUT_TBL = ROOT / "output" / "tables"
OUT_FIG = ROOT / "output" / "figures"
OUT_TBL.mkdir(parents=True, exist_ok=True)
OUT_FIG.mkdir(parents=True, exist_ok=True)

theta_iso, w_iso = isotropic_pdf(N_THETA)


# ────────────────────────────────────────────────────────────────────
# (1) Beta prior sensitivity
# ────────────────────────────────────────────────────────────────────
print("=" * 72)
print("  C2-1  Beta prior sensitivity")
print("=" * 72)
records_prior = []
arrays_prior = {}        # (variant, dr) -> rho array

for variant, (a, b) in PRIOR_VARIANTS.items():
    pri = GeometryPriors(ecc_norm_alpha=a, ecc_norm_beta=b)
    for dr in DECAY_FOCUS:
        secs, _ = sample_eccentric_circular(
            decay_ratio=dr, n_samples=N_SAMPLES, R=1.0,
            priors=pri, seed=100 + int(dr * 100))
        rhos = life_ratios(secs, dr, theta_iso, w_iso, b=B_LOWER)
        arrays_prior[(variant, dr)] = rhos
        s = summarize(rhos)
        rec = dict(prior=variant, alpha=a, beta=b, decay_ratio=dr, **s)
        records_prior.append(rec)
        print(f"  {variant:<24} r_d/R={dr}  median={s['p50']:.3f}  "
              f"5th={s['p5']:.3f}  95th={s['p95']:.3f}")


# ────────────────────────────────────────────────────────────────────
# (2) Elliptical aspect ratio sweep
# ────────────────────────────────────────────────────────────────────
print()
print("=" * 72)
print("  C2-2  Elliptical aspect ratio (area-equivalent)")
print("=" * 72)
records_aspect = []
arrays_aspect = {}        # (variant, dr) -> rho array

for variant, ratio_b_over_a in ASPECT_VARIANTS.items():
    for dr in DECAY_FOCUS:
        if variant == "circular":
            secs, _ = sample_eccentric_circular(
                decay_ratio=dr, n_samples=N_SAMPLES, R=1.0,
                priors=GeometryPriors(ecc_norm_alpha=2.0, ecc_norm_beta=2.0),
                seed=200 + int(dr * 100))
        else:
            # Build a custom GeometryPriors with deterministic aspect
            pri = GeometryPriors(
                ecc_norm_alpha=2.0, ecc_norm_beta=2.0,
                aspect_log_mu=np.log(1.0 / ratio_b_over_a),  # a/b = 1/ratio
                aspect_log_sigma=0.0,                         # deterministic shape
            )
            secs, _ = sample_elliptical(
                decay_ratio_equiv=dr, n_samples=N_SAMPLES, R=1.0,
                priors=pri, seed=300 + int(dr * 100))
        if not secs:
            print(f"  {variant:<14} r_d/R={dr}  (no valid samples - cavity too large)")
            continue
        rhos = life_ratios(secs, dr, theta_iso, w_iso, b=B_LOWER)
        arrays_aspect[(variant, dr)] = rhos
        s = summarize(rhos)
        rec = dict(aspect=variant, decay_ratio=dr, **s)
        records_aspect.append(rec)
        print(f"  {variant:<14} r_d/R={dr}  median={s['p50']:.3f}  "
              f"5th={s['p5']:.3f}  95th={s['p95']:.3f}  n={s['n']}")


# ────────────────────────────────────────────────────────────────────
# (3) Variance decomposition (conditional-variance Sobol)
# ────────────────────────────────────────────────────────────────────
print()
print("=" * 72)
print("  C2-3  Variance decomposition (which variable drives rho?)")
print("=" * 72)

# Use the baseline circular case at r_d/R = 0.7 with explicit (ecc_norm, ecc_angle) record
def conditional_variance_indices(rhos: np.ndarray,
                                  variables: dict[str, np.ndarray],
                                  n_bins: int = 20) -> dict[str, float]:
    """
    First-order Sobol index estimator via binned conditional means:
        S_i = Var[ E[Y | X_i] ] / Var[Y]
    Returns one S_i per input variable.
    """
    var_y = float(np.var(rhos, ddof=0))
    out = {}
    for name, x in variables.items():
        bins = np.linspace(x.min(), x.max(), n_bins + 1)
        idx = np.digitize(x, bins[1:-1])
        cond_mean = np.array([rhos[idx == k].mean() if (idx == k).any()
                              else 0.0 for k in range(n_bins)])
        cond_count = np.array([(idx == k).sum() for k in range(n_bins)])
        # weighted variance of conditional means
        w = cond_count / cond_count.sum()
        var_cm = float(np.sum(w * (cond_mean - rhos.mean()) ** 2))
        out[name] = var_cm / var_y if var_y > 0 else 0.0
    return out


# Re-sample with explicit recording of the variables
rng = np.random.default_rng(123)
N_VD = 10_000

# circular eccentric case
ecc_norm_arr  = rng.beta(2.0, 2.0, N_VD)
ecc_angle_arr = rng.uniform(0.0, 2.0 * np.pi, N_VD)

rho_vd = np.zeros(N_VD)
for i, (en, ea) in enumerate(zip(ecc_norm_arr, ecc_angle_arr)):
    sec = DecaySection.eccentric_circular(R=1.0, decay_ratio=0.7,
                                          ecc_norm=en, ecc_angle=ea)
    props = sec.section_properties()
    s_arr = np.array([saf(1.0, props, t) for t in theta_iso])
    d = float(np.sum(w_iso * s_arr ** (1.0 / B_LOWER)))
    rho_vd[i] = (1.0 / (1.0 - 0.7 ** 4)) ** (1.0 / B_LOWER) / d

S_circ = conditional_variance_indices(
    rho_vd,
    {"ecc_norm": ecc_norm_arr, "ecc_angle": ecc_angle_arr})
print(f"  Circular eccentric (r_d/R=0.7, n={N_VD}):")
for k, v in S_circ.items():
    print(f"    S_{k} = {v:.3f}")
print(f"    sum first-order = {sum(S_circ.values()):.3f}  "
      f"(remainder is interaction)")

# elliptical case (with aspect added)
ecc_norm_e  = rng.beta(2.0, 2.0, N_VD)
ecc_angle_e = rng.uniform(0.0, 2.0 * np.pi, N_VD)
phi_e       = rng.uniform(0.0, 2.0 * np.pi, N_VD)
aspect_e    = np.exp(rng.normal(0.0, 0.4, N_VD))

rho_e = np.zeros(N_VD)
keep = np.ones(N_VD, dtype=bool)
r_eq = 0.7
for i in range(N_VD):
    a = np.sqrt(aspect_e[i]) * r_eq
    b = r_eq / np.sqrt(aspect_e[i])
    if max(a, b) >= 0.999:
        keep[i] = False
        continue
    sec = DecaySection.elliptical(R=1.0, a_ratio=a, b_ratio=b,
                                  ecc_norm=ecc_norm_e[i],
                                  ecc_angle=ecc_angle_e[i],
                                  phi=phi_e[i])
    props = sec.section_properties()
    s_arr = np.array([saf(1.0, props, t) for t in theta_iso])
    d = float(np.sum(w_iso * s_arr ** (1.0 / B_LOWER)))
    rho_e[i] = (1.0 / (1.0 - 0.7 ** 4)) ** (1.0 / B_LOWER) / d

S_ell = conditional_variance_indices(
    rho_e[keep],
    {"ecc_norm": ecc_norm_e[keep],
     "ecc_angle": ecc_angle_e[keep],
     "aspect": aspect_e[keep],
     "phi": phi_e[keep]})
print(f"  Elliptical (r_d/R=0.7 area-equiv, n={keep.sum()}):")
for k, v in S_ell.items():
    print(f"    S_{k} = {v:.3f}")
print(f"    sum first-order = {sum(S_ell.values()):.3f}  "
      f"(remainder is interaction)")


# ────────────────────────────────────────────────────────────────────
# Save tables
# ────────────────────────────────────────────────────────────────────
df_prior  = pd.DataFrame(records_prior)
df_aspect = pd.DataFrame(records_aspect)
df_var = pd.DataFrame([
    dict(case="circular_rd07", **{f"S_{k}": v for k, v in S_circ.items()},
         total=sum(S_circ.values())),
    dict(case="elliptical_rd07", **{f"S_{k}": v for k, v in S_ell.items()},
         total=sum(S_ell.values())),
])

with pd.ExcelWriter(OUT_TBL / "sensitivity_summary.xlsx") if False else \
     open(OUT_TBL / "sensitivity_summary.csv", "w", encoding="utf-8") as f:
    f.write("# Phase C2-1 prior sensitivity\n")
df_prior.to_csv(OUT_TBL / "sensitivity_prior.csv", index=False, float_format="%.6f")
df_aspect.to_csv(OUT_TBL / "sensitivity_aspect.csv", index=False, float_format="%.6f")
df_var.to_csv(OUT_TBL / "sensitivity_variance.csv", index=False, float_format="%.4f")
print()
print(f"Saved: {OUT_TBL / 'sensitivity_prior.csv'}")
print(f"Saved: {OUT_TBL / 'sensitivity_aspect.csv'}")
print(f"Saved: {OUT_TBL / 'sensitivity_variance.csv'}")


# ────────────────────────────────────────────────────────────────────
# Visualisations
# ────────────────────────────────────────────────────────────────────
plt.rcParams.update({
    "font.family": "Times New Roman",
    "axes.titlesize": 11, "axes.labelsize": 10,
    "legend.fontsize": 8, "xtick.labelsize": 9, "ytick.labelsize": 9,
    "axes.spines.top": False, "axes.spines.right": False,
})

# (1) Prior sensitivity figure
fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), sharey=True)
COLORS_PRIOR = ["tab:blue", "tab:orange", "tab:green", "tab:red"]
for ax, dr in zip(axes, DECAY_FOCUS):
    for (variant, _), color in zip(PRIOR_VARIANTS.items(), COLORS_PRIOR):
        rhos = np.sort(arrays_prior[(variant, dr)])
        cdf = np.arange(1, len(rhos) + 1) / len(rhos)
        ax.plot(rhos, cdf, lw=1.8, color=color, label=variant)
    ax.axvline(1.0, color="k", lw=0.8, ls=":", label="concentric")
    ax.set_xlabel(r"life ratio  $\rho$")
    ax.set_xlim(0, 1.1)
    ax.set_title(rf"$r_d/R = {dr}$")
    ax.grid(alpha=0.3)
    if ax is axes[0]:
        ax.set_ylabel("CDF")
        ax.legend(loc="lower right")
fig.suptitle("Sensitivity to Beta prior on normalised eccentricity (lower-bound S-N)",
             fontsize=12)
fig.tight_layout()
fig.savefig(OUT_FIG / "fig_sens_prior.png", dpi=200, bbox_inches="tight")
print(f"Saved: {OUT_FIG / 'fig_sens_prior.png'}")
plt.close(fig)

# (2) Aspect ratio figure
fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), sharey=True)
COLORS_ASP = ["tab:blue", "tab:orange", "tab:red"]
for ax, dr in zip(axes, DECAY_FOCUS):
    for (variant, _), color in zip(ASPECT_VARIANTS.items(), COLORS_ASP):
        if (variant, dr) not in arrays_aspect:
            continue
        rhos = np.sort(arrays_aspect[(variant, dr)])
        cdf = np.arange(1, len(rhos) + 1) / len(rhos)
        ax.plot(rhos, cdf, lw=1.8, color=color, label=variant)
    ax.axvline(1.0, color="k", lw=0.8, ls=":", label="concentric")
    ax.set_xlabel(r"life ratio  $\rho$")
    ax.set_xlim(0, 1.1)
    ax.set_title(rf"$r_d/R$ area-equiv = ${dr}$")
    ax.grid(alpha=0.3)
    if ax is axes[0]:
        ax.set_ylabel("CDF")
        ax.legend(loc="lower right")
fig.suptitle("Sensitivity to cavity aspect ratio (lower-bound S-N, baseline prior)",
             fontsize=12)
fig.tight_layout()
fig.savefig(OUT_FIG / "fig_sens_aspect.png", dpi=200, bbox_inches="tight")
print(f"Saved: {OUT_FIG / 'fig_sens_aspect.png'}")
plt.close(fig)

# (3) Variance decomposition bar chart
fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
for ax, (case, S) in zip(axes,
                          [("Circular eccentric (rd/R=0.7)", S_circ),
                           ("Elliptical (rd/R=0.7 area-equiv)", S_ell)]):
    names = list(S.keys())
    vals  = list(S.values())
    colors = ["tab:blue", "tab:orange", "tab:green", "tab:red"][:len(names)]
    ax.bar(names, vals, color=colors)
    rem = max(0.0, 1.0 - sum(vals))
    ax.bar(["interactions"], [rem], color="lightgray", hatch="//")
    ax.set_ylim(0, 1.0)
    ax.set_ylabel("first-order Sobol index  $S_i$")
    ax.set_title(case)
    for x, v in enumerate(vals + [rem]):
        ax.text(x, v + 0.02, f"{v:.2f}", ha="center", fontsize=9)
    ax.grid(axis="y", alpha=0.3)
fig.suptitle("Which variable drives the life ratio?  (variance decomposition)",
             fontsize=12)
fig.tight_layout()
fig.savefig(OUT_FIG / "fig_sens_variance.png", dpi=200, bbox_inches="tight")
print(f"Saved: {OUT_FIG / 'fig_sens_variance.png'}")
plt.close(fig)
