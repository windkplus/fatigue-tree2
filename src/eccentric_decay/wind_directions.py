"""
wind_directions.py - Wind direction distributions for fatigue Monte Carlo.

Three families are provided:

  isotropic_pdf()      : uniform over [0, 2π); the most defensible / conservative
                         baseline when no site-specific wind rose is available.

  vonmises_pdf()       : circular Gaussian (von Mises) with concentration κ.
                         κ = 0 reproduces isotropic; large κ models a single
                         dominant direction (e.g., typhoon corridor).

  empirical_wind_rose(): user-supplied frequency dictionary keyed by direction
                         (degrees, meteorological convention -- direction the
                         wind comes FROM).

The Korean rose values used here (SEOUL_WIND_ROSE, JEJU_WIND_ROSE) are
APPROXIMATE placeholders synthesized from KMA general climatology; they should
be replaced with site-specific KMA AWS frequency data weighted by the
high-wind tail (e.g. ≥ 4 m/s) before publication, since fatigue is dominated
by the upper-velocity bins.

A "damage factor" helper combines a directional SAF curve with a wind-direction
distribution and the S-N exponent 1/b:

    D_factor = E_θ [ SAF(θ)^(1/b) ]                                  (1)

By the multiplicative scaling of cyclic stress, the Miner damage of a stress
history scaled by SAF(θ) is the baseline damage times SAF(θ)^(1/b). Equation
(1) therefore gives the *expected* damage multiplier under the prescribed wind
direction distribution -- and the corresponding fatigue life is

    T_eff = T_baseline / D_factor.

Reference
---------
  Korea Meteorological Administration (KMA) climatological normals.
  Mardia, K.V. & Jupp, P.E. (2000). Directional Statistics. Wiley.
"""

from __future__ import annotations

from typing import Dict, Tuple

import numpy as np


# ─────────────────────────────────────────────────────────────────────
# Approximate 8-bin wind roses (placeholder - replace with KMA AWS data)
#
# Keys are angles in degrees (meteorological "FROM" direction).
# Values are relative frequencies (will be normalized).
# ─────────────────────────────────────────────────────────────────────

SEOUL_WIND_ROSE: Dict[float, float] = {
    # NW dominant in winter (typhoon-relevant high-wind tail biased here)
    0.0:    8.0,   # N
    45.0:   10.0,  # NE
    90.0:   12.0,  # E
    135.0:  18.0,  # SE
    180.0:  10.0,  # S
    225.0:  12.0,  # SW
    270.0:  15.0,  # W
    315.0:  15.0,  # NW
}

JEJU_WIND_ROSE: Dict[float, float] = {
    # N/NE dominant year-round; also notable SW track during typhoon season
    0.0:    18.0,  # N
    45.0:   20.0,  # NE
    90.0:   12.0,  # E
    135.0:  10.0,  # SE
    180.0:   8.0,  # S
    225.0:  12.0,  # SW
    270.0:  10.0,  # W
    315.0:  10.0,  # NW
}


# ─────────────────────────────────────────────────────────────────────
# Distribution generators
# ─────────────────────────────────────────────────────────────────────

def isotropic_pdf(n_theta: int = 360) -> Tuple[np.ndarray, np.ndarray]:
    """Uniform circular distribution: theta_i = 2π i / n, w_i = 1/n."""
    theta = np.linspace(0.0, 2.0 * np.pi, n_theta, endpoint=False)
    weights = np.full(n_theta, 1.0 / n_theta)
    return theta, weights


def vonmises_pdf(kappa: float, mu: float = 0.0,
                 n_theta: int = 360) -> Tuple[np.ndarray, np.ndarray]:
    """
    Discretised von Mises (circular Gaussian) distribution.

    κ = 0    -> isotropic
    κ = 0.5  -> mild concentration
    κ = 2    -> strong concentration (~75% within 90 deg sector around mu)
    """
    theta = np.linspace(0.0, 2.0 * np.pi, n_theta, endpoint=False)
    if kappa < 1e-9:
        weights = np.full(n_theta, 1.0 / n_theta)
        return theta, weights
    # unnormalised von Mises pdf: exp(kappa * cos(theta - mu)) / (2*pi*I0(kappa))
    # The I0 normalising constant cancels when we re-normalise discretely below.
    log_pdf = kappa * np.cos(theta - mu)
    log_pdf -= log_pdf.max()                  # stabilise
    pdf = np.exp(log_pdf)
    weights = pdf / pdf.sum()
    return theta, weights


def empirical_wind_rose(rose: Dict[float, float],
                        n_theta: int = 360) -> Tuple[np.ndarray, np.ndarray]:
    """
    Discretise an N-bin wind rose onto a uniform `n_theta`-point grid.

    The bin boundaries are the midpoints between consecutive rose-key angles;
    each grid point is assigned to the closest rose key (cyclic distance).
    """
    theta = np.linspace(0.0, 2.0 * np.pi, n_theta, endpoint=False)
    keys_deg = np.array(sorted(rose.keys()))
    keys_rad = np.deg2rad(keys_deg)
    freqs = np.array([rose[k] for k in keys_deg], dtype=float)
    # find the closest rose key for every grid point (cyclic min-distance)
    diff = (theta[:, None] - keys_rad[None, :] + np.pi) % (2.0 * np.pi) - np.pi
    closest = np.argmin(np.abs(diff), axis=1)
    weights = freqs[closest]
    weights = weights / weights.sum()
    return theta, weights


# ─────────────────────────────────────────────────────────────────────
# Damage factor helper
# ─────────────────────────────────────────────────────────────────────

def damage_factor(saf_arr: np.ndarray, weights: np.ndarray, b: float) -> float:
    """
    E_θ [ SAF(θ)^(1/b) ] given a precomputed SAF directional sweep.

    Parameters
    ----------
    saf_arr : (n_theta,) array of SAF values, ALREADY ALIGNED with `weights`.
    weights : (n_theta,) wind direction probabilities (sum to 1).
    b       : S-N curve exponent (lower-bound 0.125, mean 0.10).

    Returns
    -------
    D : scalar damage multiplier.
    """
    return float(np.sum(weights * saf_arr ** (1.0 / b)))


# ─────────────────────────────────────────────────────────────────────
# Convenience: build all three reference distributions
# ─────────────────────────────────────────────────────────────────────

def reference_distributions(n_theta: int = 360) -> Dict[str, Tuple[np.ndarray, np.ndarray]]:
    """Return {name: (theta, weights)} for isotropic, Seoul, Jeju."""
    return {
        "isotropic": isotropic_pdf(n_theta),
        "seoul":     empirical_wind_rose(SEOUL_WIND_ROSE, n_theta),
        "jeju":      empirical_wind_rose(JEJU_WIND_ROSE, n_theta),
    }
