"""
stress_amplification.py — Section modulus and Stress Amplification Factor (SAF)
for eccentric / asymmetric tree cross-sections under wind bending.

Bending convention
------------------
A unit horizontal wind force at the canopy applies a base moment M about an axis
perpendicular to the load direction. Let θ_w be the load direction (CCW from x-axis).

For a point (x, y) in centroidal coordinates, the perpendicular distance to the
bending axis (a line through the centroid at angle θ_w + π/2) — measured ALONG
the load direction — is

    d(x, y; θ_w) = x cos θ_w + y sin θ_w.

Hence the moment of inertia about the bending axis is

    I(θ_w) = ∫ d² dA
           = sin² θ_w · I_xx + cos² θ_w · I_yy + 2 sin θ_w cos θ_w · I_xy,

with I_xx = ∫ y² dA, I_yy = ∫ x² dA, I_xy = ∫ x y dA (engineering convention).

The maximum |d| over the outer circle (radius R, centered at the un-shifted origin
which lies at (-x_c, -y_c) in centroidal coordinates) is

    |d|_max = R + |x_c cos θ_w + y_c sin θ_w|,

so the section modulus is

    S(θ_w) = I(θ_w) / |d|_max,

and the stress amplification factor relative to the sound section S₀ = π R³ / 4 is

    SAF(θ_w) = S₀ / S(θ_w).

In the concentric limit (x_c = y_c = 0, I_xx = I_yy = π (R⁴ - r_d⁴)/4, I_xy = 0)
this reduces to SAF = 1 / (1 - (r_d/R)⁴), matching Mattheck (1994).
"""

from __future__ import annotations

from typing import Tuple

import numpy as np


def section_modulus(R: float, props: dict, theta_w: float) -> float:
    """Net-section modulus for a wind load applied in direction θ_w (rad)."""
    Ixx = props["I_xx"]
    Iyy = props["I_yy"]
    Ixy = props["I_xy"]
    x_c = props["x_c"]
    y_c = props["y_c"]

    cw, sw = np.cos(theta_w), np.sin(theta_w)
    I_bend = sw * sw * Ixx + cw * cw * Iyy + 2.0 * sw * cw * Ixy

    offset_along_load = x_c * cw + y_c * sw
    c_max = R + abs(offset_along_load)

    return I_bend / c_max


def saf(R: float, props: dict, theta_w: float) -> float:
    """Stress amplification factor at a single wind direction."""
    S_sound = np.pi * R ** 3 / 4.0
    S_dec   = section_modulus(R, props, theta_w)
    return S_sound / S_dec


def saf_directional(R: float, props: dict, n_theta: int = 360) -> np.ndarray:
    """
    Sweep SAF over the full circle of wind directions.

    Returns
    -------
    theta : (n_theta,) array of directions in [0, 2π)
    saf_arr : (n_theta,) array of SAF values
    """
    theta = np.linspace(0.0, 2.0 * np.pi, n_theta, endpoint=False)
    saf_arr = np.array([saf(R, props, t) for t in theta])
    return theta, saf_arr


def saf_worst(R: float, props: dict, n_theta: int = 360) -> Tuple[float, float]:
    """Return (worst-case SAF, corresponding wind direction in radians)."""
    theta, saf_arr = saf_directional(R, props, n_theta=n_theta)
    i = int(np.argmax(saf_arr))
    return float(saf_arr[i]), float(theta[i])
