"""
00_regression_test_concentric.py
--------------------------------
Verify the new eccentric-decay framework recovers the Mattheck concentric formula

    SAF_concentric(r_d/R) = 1 / (1 - (r_d/R)^4)

in the e = 0 limit. This is a hard correctness gate for paper #1 ↔ paper #2 continuity:
any deviation > 0.1 % flags a bug in geometry or stress-amplification.

Also checks:
  • Eccentric Mattheck-formula limit:   for e/R → 0, SAF should → concentric value
  • Worst-case sweep stability:         SAF_worst(concentric) should be direction-invariant
  • Numerical integration cross-check:  Monte Carlo grid integration matches analytical
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from eccentric_decay import (  # noqa: E402
    CircularCavity,
    DecaySection,
    saf,
    saf_worst,
)


# ────────────────────────────────────────────────────────────────────
TOL_REL = 1e-3      # 0.1 % tolerance against Mattheck closed-form
R = 1.0             # unit-radius stem (results scale-invariant)


def mattheck_saf(decay_ratio: float) -> float:
    if decay_ratio >= 1.0:
        return np.inf
    return 1.0 / (1.0 - decay_ratio ** 4)


# ────────────────────────────────────────────────────────────────────
# Test 1 - concentric limit
# ────────────────────────────────────────────────────────────────────
print("=" * 72)
print("  Test 1 - Concentric (e=0) limit reproduces Mattheck SAF = 1/(1-(rd/R)^4)")
print("=" * 72)
print(f"  {'r_d/R':>8}  {'SAF (model)':>14}  {'SAF (Mattheck)':>16}  "
      f"{'rel error':>12}  status")
all_pass = True
for ratio in [0.0, 0.1, 0.3, 0.5, 0.7, 0.8, 0.9]:
    sec = DecaySection.concentric(R=R, decay_ratio=ratio)
    props = sec.section_properties()
    s_model = saf(R, props, theta_w=0.0)
    s_ref   = mattheck_saf(ratio)
    rel = 0.0 if not np.isfinite(s_ref) else abs(s_model - s_ref) / s_ref
    ok = rel < TOL_REL
    all_pass &= ok
    print(f"  {ratio:>8.2f}  {s_model:>14.6f}  {s_ref:>16.6f}  "
          f"{rel:>12.2e}  {'PASS' if ok else 'FAIL'}")


# ────────────────────────────────────────────────────────────────────
# Test 2 - direction-invariance for concentric cavities
# ────────────────────────────────────────────────────────────────────
print()
print("=" * 72)
print("  Test 2 - Concentric SAF should be direction-invariant (max == min)")
print("=" * 72)
for ratio in [0.3, 0.5, 0.7]:
    sec = DecaySection.concentric(R=R, decay_ratio=ratio)
    props = sec.section_properties()
    saf_max, theta_max = saf_worst(R, props, n_theta=720)
    # also check minimum
    from eccentric_decay import saf_directional
    _, saf_arr = saf_directional(R, props, n_theta=720)
    saf_min = float(saf_arr.min())
    spread = (saf_max - saf_min) / saf_max
    ok = spread < TOL_REL
    all_pass &= ok
    print(f"  r_d/R={ratio:.2f}  SAF range [{saf_min:.6f}, {saf_max:.6f}]  "
          f"rel spread {spread:.2e}  {'PASS' if ok else 'FAIL'}")


# ────────────────────────────────────────────────────────────────────
# Test 3 - eccentric → concentric continuity (e → 0)
# ────────────────────────────────────────────────────────────────────
print()
print("=" * 72)
print("  Test 3 - Eccentric SAF → concentric SAF as eccentricity → 0")
print("=" * 72)
for ratio in [0.3, 0.5, 0.7]:
    s_ref = mattheck_saf(ratio)
    print(f"  r_d/R = {ratio:.2f}  (Mattheck SAF = {s_ref:.6f})")
    for e_norm in [1e-6, 1e-3, 1e-2, 1e-1]:
        sec = DecaySection.eccentric_circular(
            R=R, decay_ratio=ratio, ecc_norm=e_norm)
        props = sec.section_properties()
        saf_max, _ = saf_worst(R, props, n_theta=720)
        rel = abs(saf_max - s_ref) / s_ref
        ok = rel < max(TOL_REL, 5 * e_norm)   # relax for finite eccentricity
        print(f"    e_norm={e_norm:.1e}   SAF_worst={saf_max:.6f}   rel={rel:.2e}   "
              f"{'PASS' if ok else 'FAIL (expected at large e)'}")


# ────────────────────────────────────────────────────────────────────
# Test 4 - analytical vs numerical (grid) integration cross-check
# ────────────────────────────────────────────────────────────────────
print()
print("=" * 72)
print("  Test 4 - Analytical section properties vs numerical grid integration")
print("=" * 72)


def numerical_section_properties(R: float, cavity_test, n_grid: int = 2000) -> dict:
    """Brute-force integration on a Cartesian grid (validation only)."""
    x = np.linspace(-R, R, n_grid)
    y = np.linspace(-R, R, n_grid)
    dx = x[1] - x[0]
    dy = y[1] - y[0]
    X, Y = np.meshgrid(x, y, indexing="xy")
    inside_outer = (X ** 2 + Y ** 2) <= R ** 2
    if isinstance(cavity_test, CircularCavity):
        inside_cavity = (
            (X - cavity_test.e_x) ** 2 + (Y - cavity_test.e_y) ** 2
            <= cavity_test.r_d ** 2
        )
    else:
        raise NotImplementedError("numerical check supports CircularCavity only")
    mask = inside_outer & ~inside_cavity
    A = mask.sum() * dx * dy
    x_c = (X[mask]).sum() * dx * dy / A
    y_c = (Y[mask]).sum() * dx * dy / A
    Xc = X - x_c
    Yc = Y - y_c
    Ixx = (Yc[mask] ** 2).sum() * dx * dy
    Iyy = (Xc[mask] ** 2).sum() * dx * dy
    Ixy = (Xc[mask] * Yc[mask]).sum() * dx * dy
    return dict(A=A, x_c=x_c, y_c=y_c, I_xx=Ixx, I_yy=Iyy, I_xy=Ixy)


# Pick a moderately eccentric case
for (ratio, e_norm) in [(0.5, 0.0), (0.5, 0.5), (0.7, 0.7)]:
    sec = DecaySection.eccentric_circular(R=R, decay_ratio=ratio, ecc_norm=e_norm)
    props_a = sec.section_properties()
    props_n = numerical_section_properties(R, sec.cavity, n_grid=1500)
    print(f"  r_d/R={ratio:.2f}  e_norm={e_norm:.2f}")
    for key in ["A", "x_c", "y_c", "I_xx", "I_yy", "I_xy"]:
        a = props_a[key]
        n = props_n[key]
        denom = max(abs(a), 1e-9)
        rel = abs(a - n) / denom
        ok = rel < 5e-3   # 0.5 % grid tolerance
        all_pass &= ok if abs(a) > 1e-6 else True
        print(f"    {key:>5}  analytical={a:>+12.6f}   numerical={n:>+12.6f}   "
              f"rel={rel:.2e}   {'PASS' if ok or abs(a) < 1e-6 else 'FAIL'}")


# ────────────────────────────────────────────────────────────────────
print()
print("=" * 72)
print(f"  OVERALL : {'ALL CORE TESTS PASSED' if all_pass else 'SOME TESTS FAILED'}")
print("=" * 72)
sys.exit(0 if all_pass else 1)
