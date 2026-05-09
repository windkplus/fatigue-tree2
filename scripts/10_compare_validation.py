"""
10_compare_validation.py
------------------------
Compare ABAQUS-computed cross-section properties with our closed-form
predictions for the six paper-#2 Appendix-A validation cases.

Inputs
------
  data/abaqus_validation/validation_results.csv   (from ABAQUS office PC)

Outputs
-------
  output/tables/abaqus_validation_comparison.csv
  output/figures/fig_abaqus_validation.png
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from eccentric_decay import DecaySection  # noqa: E402

R = 1.0

CASES = [
    ("val_rd50_e000", 0.5, 0.00),
    ("val_rd50_e050", 0.5, 0.50),
    ("val_rd50_e095", 0.5, 0.95),
    ("val_rd70_e000", 0.7, 0.00),
    ("val_rd70_e050", 0.7, 0.50),
    ("val_rd70_e095", 0.7, 0.95),
]

# Load ABAQUS results
abq_csv = ROOT / "data" / "abaqus_validation" / "validation_results.csv"
abq = pd.read_csv(abq_csv).set_index("job")
print("ABAQUS validation rows loaded:")
print(abq.to_string())
print()


# ─────────────────────────────────────────────────────────────────────
# Build comparison table
# ─────────────────────────────────────────────────────────────────────
records = []
QUANTS = ["area", "x_c", "y_c", "I_xx", "I_yy", "I_xy"]

# Mapping of ABAQUS column -> our property dict key
ABQ_TO_OURS = {
    "area": "A",
    "x_c":  "x_c",
    "y_c":  "y_c",
    "I_xx": "I_xx",
    "I_yy": "I_yy",
    "I_xy": "I_xy",
}

for name, rdR, e_norm in CASES:
    row = abq.loc[name]
    sec = DecaySection.eccentric_circular(R, rdR, ecc_norm=e_norm)
    props = sec.section_properties()
    rec = dict(case=name, r_dR=rdR, e_norm=e_norm)
    for col, key in ABQ_TO_OURS.items():
        a_val = float(row[col])
        c_val = float(props[key])
        # absolute error (good for near-zero values)
        abs_err = abs(a_val - c_val)
        # relative error (only when denominator non-trivial)
        denom = max(abs(c_val), 1e-9)
        rel_err = abs_err / denom if abs(c_val) > 1e-6 else float("nan")
        rec[f"{col}_abq"] = a_val
        rec[f"{col}_cf"]  = c_val
        rec[f"{col}_abs"] = abs_err
        rec[f"{col}_rel"] = rel_err
    records.append(rec)

cmp = pd.DataFrame(records)
out_csv = ROOT / "output" / "tables" / "abaqus_validation_comparison.csv"
cmp.to_csv(out_csv, index=False, float_format="%.8f")
print(f"\nComparison saved: {out_csv}")


# ─────────────────────────────────────────────────────────────────────
# Pretty-print headline table
# ─────────────────────────────────────────────────────────────────────
def fmt_pair(a, c, abs_err):
    return f"{a:+10.5f} / {c:+10.5f}  (abs {abs_err:.2e})"


print("\n" + "=" * 90)
print("  ABAQUS (CPS4R, ~15k elems) vs closed-form")
print("=" * 90)
header = f"{'case':<16}{'quant':>7}  {'ABAQUS / closed-form (abs err)':<46}{'rel err':>12}"
print(header)
print("-" * 90)

max_abs = {q: 0.0 for q in QUANTS}
max_rel = 0.0
for r in records:
    for col in QUANTS:
        a = r[f"{col}_abq"]
        c = r[f"{col}_cf"]
        abs_e = r[f"{col}_abs"]
        rel_e = r[f"{col}_rel"]
        rel_str = f"{rel_e:.2e}" if not np.isnan(rel_e) else "  (zero)"
        print(f"{r['case']:<16}{col:>7}  "
              f"{a:+10.5f} / {c:+10.5f}  (abs {abs_e:.2e})  {rel_str:>12}")
        max_abs[col] = max(max_abs[col], abs_e)
        if not np.isnan(rel_e):
            max_rel = max(max_rel, rel_e)
    print()
print("-" * 90)
print(f"  Max absolute error per quantity:")
for q, v in max_abs.items():
    print(f"    {q:>5}: {v:.2e}")
print(f"  Maximum relative error overall: {max_rel:.2e}")


# ─────────────────────────────────────────────────────────────────────
# Visualisation
# ─────────────────────────────────────────────────────────────────────
plt.rcParams.update({
    "font.family": "Times New Roman",
    "axes.titlesize": 11, "axes.labelsize": 10,
    "legend.fontsize": 9, "xtick.labelsize": 9, "ytick.labelsize": 9,
    "axes.spines.top": False, "axes.spines.right": False,
})

fig, axes = plt.subplots(1, 2, figsize=(12, 5.0))

# (a) Bar chart of relative errors per quantity per case
ax = axes[0]
case_labels = [
    f"$r_d/R={r['r_dR']}$\n$e={r['e_norm']:.2f}$" for r in records]
quants_to_plot = ["area", "x_c", "I_xx", "I_yy"]
COLORS = ["tab:blue", "tab:orange", "tab:green", "tab:red"]
n_q = len(quants_to_plot)
n_c = len(records)
bar_w = 0.8 / n_q
xs = np.arange(n_c)
for i, q in enumerate(quants_to_plot):
    rels = []
    for r in records:
        rv = r[f"{q}_rel"]
        # for cases where x_c=0 (concentric), rel err is undefined -> use abs
        if np.isnan(rv):
            rv = r[f"{q}_abs"]
        rels.append(max(rv, 1e-10))           # log-safe
    ax.bar(xs + i * bar_w - 0.4 + bar_w / 2,
           rels, width=bar_w, color=COLORS[i], label=q)
ax.set_yscale("log")
ax.set_xticks(xs)
ax.set_xticklabels(case_labels, fontsize=8)
ax.axhline(1e-3, color="gray", lw=1, ls=":", alpha=0.7)
ax.text(n_c - 0.5, 1.2e-3, "0.1 %", color="gray", ha="right", fontsize=8)
ax.set_ylabel("relative error (or absolute, for zero quantities)")
ax.set_title("(a) Per-quantity error  ABAQUS  vs  closed-form")
ax.legend(loc="upper left", ncol=2)
ax.grid(axis="y", which="both", alpha=0.3)

# (b) ABAQUS vs closed-form scatter
ax = axes[1]
all_cf  = []
all_abq = []
for r in records:
    for q in ["area", "x_c", "y_c", "I_xx", "I_yy", "I_xy"]:
        all_cf.append(r[f"{q}_cf"])
        all_abq.append(r[f"{q}_abq"])
all_cf = np.array(all_cf)
all_abq = np.array(all_abq)
ax.plot([all_cf.min() - 0.1, all_cf.max() + 0.1],
        [all_cf.min() - 0.1, all_cf.max() + 0.1],
        "k--", lw=1, label="y = x")
ax.scatter(all_cf, all_abq, c="tab:red", s=40, alpha=0.7,
           label=f"all 36 quantities (6 cases x 6 props)")
ax.set_xlabel("closed-form")
ax.set_ylabel("ABAQUS")
ax.set_title("(b) ABAQUS  vs  closed-form (parity plot)")
ax.legend()
ax.grid(alpha=0.3)
ax.set_aspect("equal")

fig.suptitle(f"Figure A1.  ABAQUS (CPS4R, ~15k elements) vs closed-form section properties.\n"
             f"Maximum relative error across all 6 cases x 6 quantities: "
             f"{max_rel:.2e}.",
             fontsize=11)
fig.tight_layout()
out_fig = ROOT / "output" / "figures" / "fig_abaqus_validation.png"
fig.savefig(out_fig, dpi=200, bbox_inches="tight")
print(f"\nFigure saved: {out_fig}")
plt.close(fig)
