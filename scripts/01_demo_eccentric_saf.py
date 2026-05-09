"""
01_demo_eccentric_saf.py
------------------------
Demonstrate how the Stress Amplification Factor (SAF) depends on:
  (a) cavity eccentricity at fixed decay ratio
  (b) wind direction (polar plot)
  (c) decay ratio sweep, with concentric vs eccentric vs elliptical comparison

Outputs
-------
  output/figures/eccentric_demo.png
  output/tables/eccentric_demo_summary.csv

This is a sanity-check / first-look figure -- not yet the Monte Carlo result.
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
    saf_directional,
    saf_worst,
)


OUT_FIG = ROOT / "output" / "figures" / "eccentric_demo.png"
OUT_CSV = ROOT / "output" / "tables" / "eccentric_demo_summary.csv"
OUT_FIG.parent.mkdir(parents=True, exist_ok=True)
OUT_CSV.parent.mkdir(parents=True, exist_ok=True)

R = 1.0    # unit-radius stem


# ────────────────────────────────────────────────────────────────────
# Figure layout: 2 rows x 3 cols
# ────────────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(14, 9))

# ── (a) Cross-section schematic for 3 representative cases
ax1 = fig.add_subplot(2, 3, 1, aspect="equal")
ax1.set_title("Cross-section schematic\n(r_d/R = 0.5)")

theta_circ = np.linspace(0, 2 * np.pi, 200)
ax1.plot(R * np.cos(theta_circ), R * np.sin(theta_circ),
         "k-", lw=2, label="outer (intact)")

cases = [
    ("concentric",      DecaySection.concentric(R, 0.5),                "tab:blue"),
    ("eccentric e/R=0.5", DecaySection.eccentric_circular(R, 0.5, 0.5),  "tab:orange"),
    ("elliptical 2:1",  DecaySection.elliptical(R, a_ratio=0.6, b_ratio=0.3,
                                                ecc_norm=0.4), "tab:red"),
]
for label, sec, color in cases:
    cav = sec.cavity
    if hasattr(cav, "r_d"):
        xs = cav.e_x + cav.r_d * np.cos(theta_circ)
        ys = cav.e_y + cav.r_d * np.sin(theta_circ)
    else:
        c, s = np.cos(cav.phi), np.sin(cav.phi)
        local_x = cav.a * np.cos(theta_circ)
        local_y = cav.b * np.sin(theta_circ)
        xs = cav.e_x + c * local_x - s * local_y
        ys = cav.e_y + s * local_x + c * local_y
    ax1.fill(xs, ys, alpha=0.4, color=color, label=label)
    props = sec.section_properties()
    ax1.plot(props["x_c"], props["y_c"], "x", color=color, ms=10, mew=2)

ax1.set_xlim(-1.2, 1.2)
ax1.set_ylim(-1.2, 1.2)
ax1.legend(fontsize=8, loc="upper right")
ax1.grid(alpha=0.3)
ax1.set_xlabel("x / R")
ax1.set_ylabel("y / R")
ax1.text(0.02, 0.98, "x = net centroid", transform=ax1.transAxes,
         fontsize=8, va="top", style="italic")


# ── (b) SAF vs wind direction (polar)
ax2 = fig.add_subplot(2, 3, 2, projection="polar")
ax2.set_title("SAF vs wind direction\n(r_d/R = 0.5)", pad=18)
for label, sec, color in cases:
    props = sec.section_properties()
    theta, saf_arr = saf_directional(R, props, n_theta=720)
    ax2.plot(theta, saf_arr, color=color, lw=2, label=label)
ax2.set_theta_zero_location("E")
ax2.set_theta_direction(1)
ax2.legend(fontsize=7, loc="upper right", bbox_to_anchor=(1.4, 1.1))


# ── (c) Worst-case SAF vs eccentricity at several decay ratios
ax3 = fig.add_subplot(2, 3, 3)
ecc_arr = np.linspace(0, 0.95, 40)
ratios_to_show = [0.3, 0.5, 0.7, 0.8]
records = []
for ratio in ratios_to_show:
    saf_vals = []
    for e in ecc_arr:
        sec = DecaySection.eccentric_circular(R, ratio, ecc_norm=e)
        s_worst, _ = saf_worst(R, sec.section_properties(), n_theta=180)
        saf_vals.append(s_worst)
        records.append(dict(decay_ratio=ratio, ecc_norm=e,
                            saf_worst=s_worst,
                            mattheck=1.0 / (1.0 - ratio ** 4)))
    saf_vals = np.array(saf_vals)
    ax3.plot(ecc_arr, saf_vals, lw=2, label=f"r_d/R = {ratio:.1f}")
    # mark concentric (Mattheck) reference as dashed
    ax3.axhline(1.0 / (1.0 - ratio ** 4), color=ax3.lines[-1].get_color(),
                ls="--", lw=1, alpha=0.5)
ax3.set_xlabel("normalized eccentricity  e / (R - r_d)")
ax3.set_ylabel("worst-case SAF")
ax3.set_title("Eccentricity amplifies stress\n(dashed = concentric Mattheck)")
ax3.legend(fontsize=8)
ax3.grid(alpha=0.3)


# ── (d) SAF vs decay ratio: concentric vs eccentric (e/R-r_d=0.5) vs elliptical (2:1)
ax4 = fig.add_subplot(2, 3, 4)
decay_arr = np.linspace(0.05, 0.85, 60)
saf_conc, saf_ecc, saf_ell = [], [], []
for ratio in decay_arr:
    sc = DecaySection.concentric(R, ratio)
    se = DecaySection.eccentric_circular(R, ratio, ecc_norm=0.5)
    # elliptical with same area as the circular cavity (area-equivalent)
    # area = pi r_d^2 = pi a b -> if aspect ratio 2:1, then a = sqrt(2)*r_d, b = r_d / sqrt(2)
    r_d = ratio * R
    a = np.sqrt(2.0) * r_d
    b = r_d / np.sqrt(2.0)
    sl = DecaySection.elliptical(R, a_ratio=a / R, b_ratio=b / R,
                                 ecc_norm=0.5, phi=0.0)
    saf_conc.append(saf_worst(R, sc.section_properties(), n_theta=180)[0])
    saf_ecc.append(saf_worst(R, se.section_properties(), n_theta=180)[0])
    if a < R and b < R:
        saf_ell.append(saf_worst(R, sl.section_properties(), n_theta=180)[0])
    else:
        saf_ell.append(np.nan)
ax4.plot(decay_arr, saf_conc, lw=2, color="tab:blue", label="concentric")
ax4.plot(decay_arr, saf_ecc, lw=2, color="tab:orange",
         label="eccentric (e_norm=0.5)")
ax4.plot(decay_arr, saf_ell, lw=2, color="tab:red",
         label="elliptical 2:1, area-equiv (e_norm=0.5)")
ax4.axvline(0.7, color="gray", ls=":", lw=1)
ax4.text(0.71, 1.5, "Mattheck\nthreshold", fontsize=8, color="gray")
ax4.set_xlabel("decay ratio  r_d / R  (or area-equivalent)")
ax4.set_ylabel("worst-case SAF")
ax4.set_title("SAF vs decay ratio")
ax4.set_yscale("log")
ax4.legend(fontsize=8)
ax4.grid(which="both", alpha=0.3)


# ── (e) Centroid shift (x_c) vs eccentricity
ax5 = fig.add_subplot(2, 3, 5)
for ratio in ratios_to_show:
    xc_arr = []
    for e in ecc_arr:
        sec = DecaySection.eccentric_circular(R, ratio, ecc_norm=e)
        xc_arr.append(sec.section_properties()["x_c"])
    ax5.plot(ecc_arr, np.array(xc_arr) / R, lw=2, label=f"r_d/R = {ratio:.1f}")
ax5.set_xlabel("normalized eccentricity  e / (R - r_d)")
ax5.set_ylabel("net centroid shift  x_c / R")
ax5.set_title("Net centroid shifts\nopposite to cavity")
ax5.legend(fontsize=8)
ax5.grid(alpha=0.3)
ax5.axhline(0, color="k", lw=0.6)


# ── (f) Worst-direction angle (orientation of failure plane)
ax6 = fig.add_subplot(2, 3, 6)
for ratio in ratios_to_show:
    theta_worst = []
    for e in ecc_arr:
        sec = DecaySection.eccentric_circular(R, ratio, ecc_norm=e)
        _, t = saf_worst(R, sec.section_properties(), n_theta=720)
        # report angle in degrees, mod 180 (bending failure has 2-fold symmetry)
        theta_worst.append(np.degrees(t) % 180.0)
    ax6.plot(ecc_arr, theta_worst, lw=2, label=f"r_d/R = {ratio:.1f}")
ax6.set_xlabel("normalized eccentricity  e / (R - r_d)")
ax6.set_ylabel("worst-direction theta_w  [deg, mod 180]")
ax6.set_title("Cavity offset selects\na single failure direction")
ax6.legend(fontsize=8)
ax6.grid(alpha=0.3)
ax6.set_ylim(-5, 185)


plt.tight_layout()
plt.savefig(OUT_FIG, dpi=180, bbox_inches="tight")
print(f"Figure saved: {OUT_FIG}")


# ────────────────────────────────────────────────────────────────────
# Summary table
# ────────────────────────────────────────────────────────────────────
df = pd.DataFrame.from_records(records)
df["uplift_pct"] = (df["saf_worst"] / df["mattheck"] - 1.0) * 100.0
df.to_csv(OUT_CSV, index=False, float_format="%.6f")
print(f"Summary CSV saved: {OUT_CSV}")

print("\n--- Worst-case SAF uplift over Mattheck (concentric) baseline [%] ---")
ecc_targets = [0.0, 0.25, 0.5, 0.75, 0.95]
summary_rows = []
for ratio in ratios_to_show:
    row = {"r_d/R": ratio}
    for et in ecc_targets:
        sub = df[(df["decay_ratio"] == ratio)]
        idx = (sub["ecc_norm"] - et).abs().idxmin()
        row[f"e_norm={et:.2f}"] = sub.loc[idx, "uplift_pct"]
    summary_rows.append(row)
print(pd.DataFrame(summary_rows).round(2).to_string(index=False))
