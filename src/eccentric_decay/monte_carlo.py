"""
monte_carlo.py - Monte Carlo sampling of eccentric decay geometries and the
resulting fatigue-life distribution.

Sampling priors (defaults; user-overridable)
-------------------------------------------
  ecc_norm  ~ Beta(alpha=2, beta=2)        (symmetric about e_norm = 0.5)
  ecc_angle ~ Uniform(0, 2π)               (cavity offset direction is random)

  optional (elliptical extension):
  aspect    ~ LogNormal(0, 0.4)            (ratio a/b of the cavity)
  phi       ~ Uniform(0, 2π)               (ellipse orientation)

Damage-multiplier formulation
-----------------------------
For a fatigue framework that scales linearly with stress amplitude (Miner +
power-law S-N with exponent 1/b), the damage of a stress history multiplied
by SAF(theta) is

    D(theta) = SAF(theta)^(1/b) * D_baseline.

Hence under a wind direction distribution p(theta),

    D_eff = E_theta [ D(theta) ] = D_factor * D_baseline,
    D_factor = E_theta [ SAF(theta)^(1/b) ],

and the corresponding fatigue life is T_eff = T_baseline / D_factor.

We report results as the *life ratio* with respect to the concentric
(Mattheck) baseline at the same r_d/R, i.e.

    rho = T_eccentric / T_concentric
        = SAF_concentric^(1/b) / E_theta [ SAF(theta)^(1/b) ].

This normalisation is independent of the absolute T_baseline and lets the
results be cross-cited with paper #1 directly.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Tuple

import numpy as np

from .geometry import DecaySection
from .stress_amplification import saf
from .wind_directions import damage_factor


# ─────────────────────────────────────────────────────────────────────
# Default priors
# ─────────────────────────────────────────────────────────────────────

@dataclass
class GeometryPriors:
    """Probabilistic priors on cavity geometry."""

    # Eccentric circular cavity (always sampled)
    ecc_norm_alpha: float = 2.0          # Beta(α, β) on normalized eccentricity
    ecc_norm_beta: float = 2.0
    ecc_angle: str = "uniform"           # only "uniform" supported (per Phase B decision)

    # Optional elliptical extension (turn on with elliptical=True in sampler)
    aspect_log_mu: float = 0.0           # LogNormal(μ, σ) on a/b ratio
    aspect_log_sigma: float = 0.4
    phi_dist: str = "uniform"


# ─────────────────────────────────────────────────────────────────────
# Samplers
# ─────────────────────────────────────────────────────────────────────

def sample_eccentric_circular(
    decay_ratio: float,
    n_samples: int,
    R: float = 1.0,
    priors: GeometryPriors = None,
    seed: int = 42,
) -> Tuple[List[DecaySection], Dict[str, np.ndarray]]:
    """Draw `n_samples` eccentric-circular decay geometries at fixed `decay_ratio`."""
    rng = np.random.default_rng(seed)
    pr = priors or GeometryPriors()
    ecc_norm = rng.beta(pr.ecc_norm_alpha, pr.ecc_norm_beta, n_samples)
    ecc_angle = rng.uniform(0.0, 2.0 * np.pi, n_samples)
    sections = [
        DecaySection.eccentric_circular(R, decay_ratio,
                                        ecc_norm=float(e), ecc_angle=float(a))
        for e, a in zip(ecc_norm, ecc_angle)
    ]
    return sections, dict(ecc_norm=ecc_norm, ecc_angle=ecc_angle)


def sample_elliptical(
    decay_ratio_equiv: float,
    n_samples: int,
    R: float = 1.0,
    priors: GeometryPriors = None,
    seed: int = 42,
) -> Tuple[List[DecaySection], Dict[str, np.ndarray]]:
    """
    Draw elliptical decay geometries with the SAME area as a circular cavity
    of radius decay_ratio_equiv * R.

    For aspect = a / b, a*b = (decay_ratio_equiv*R)^2, hence
        a = sqrt(aspect) * (decay_ratio_equiv * R)
        b = (decay_ratio_equiv * R) / sqrt(aspect)
    """
    rng = np.random.default_rng(seed)
    pr = priors or GeometryPriors()
    ecc_norm = rng.beta(pr.ecc_norm_alpha, pr.ecc_norm_beta, n_samples)
    ecc_angle = rng.uniform(0.0, 2.0 * np.pi, n_samples)
    aspect = np.exp(rng.normal(pr.aspect_log_mu, pr.aspect_log_sigma, n_samples))
    phi = rng.uniform(0.0, 2.0 * np.pi, n_samples)

    sections: List[DecaySection] = []
    r_eq = decay_ratio_equiv * R
    accepted = []
    for e, a_off, asp, p in zip(ecc_norm, ecc_angle, aspect, phi):
        a = np.sqrt(asp) * r_eq
        b = r_eq / np.sqrt(asp)
        if max(a, b) >= R - 1e-3:
            # cavity would breach the outer wall - reject sample
            continue
        sections.append(DecaySection.elliptical(
            R, a_ratio=a / R, b_ratio=b / R,
            ecc_norm=float(e), ecc_angle=float(a_off), phi=float(p)))
        accepted.append(True)
    accepted = np.array(accepted)
    return sections, dict(ecc_norm=ecc_norm[accepted] if len(accepted) else ecc_norm,
                          ecc_angle=ecc_angle[accepted] if len(accepted) else ecc_angle,
                          aspect=aspect[accepted] if len(accepted) else aspect,
                          phi=phi[accepted] if len(accepted) else phi)


# ─────────────────────────────────────────────────────────────────────
# Damage / life-ratio computation
# ─────────────────────────────────────────────────────────────────────

def saf_curve(R: float, props: dict, theta_arr: np.ndarray) -> np.ndarray:
    """SAF(theta) for an array of wind directions (vectorisable in future)."""
    return np.array([saf(R, props, t) for t in theta_arr])


def life_ratios(
    sections: List[DecaySection],
    decay_ratio: float,
    theta_arr: np.ndarray,
    weights: np.ndarray,
    b: float,
) -> np.ndarray:
    """
    Compute per-sample fatigue life ratio rho = T_eccentric / T_concentric.

    rho < 1  -> eccentricity SHORTENS life vs. concentric assumption
    rho = 1  -> matches concentric
    rho > 1  -> eccentric sample is benign (rare; only when offset reduces
                 stress along the prevailing wind direction)
    """
    saf_conc = 1.0 / (1.0 - decay_ratio ** 4)
    d_conc = saf_conc ** (1.0 / b)
    R = sections[0].R
    rhos = np.zeros(len(sections))
    for i, sec in enumerate(sections):
        props = sec.section_properties()
        s_arr = saf_curve(R, props, theta_arr)
        d_ecc = damage_factor(s_arr, weights, b)
        rhos[i] = d_conc / d_ecc
    return rhos


def damage_factors(
    sections: List[DecaySection],
    theta_arr: np.ndarray,
    weights: np.ndarray,
    b: float,
) -> np.ndarray:
    """E_theta[SAF^(1/b)] per sample (raw damage multiplier vs. sound)."""
    R = sections[0].R
    out = np.zeros(len(sections))
    for i, sec in enumerate(sections):
        s_arr = saf_curve(R, sec.section_properties(), theta_arr)
        out[i] = damage_factor(s_arr, weights, b)
    return out


# ─────────────────────────────────────────────────────────────────────
# Summaries
# ─────────────────────────────────────────────────────────────────────

def summarize(arr: np.ndarray,
              percentiles=(2.5, 5, 25, 50, 75, 95, 97.5)) -> Dict[str, float]:
    """Return mean/std and the requested percentiles."""
    out = dict(mean=float(np.mean(arr)),
               std=float(np.std(arr)),
               n=int(len(arr)))
    for p in percentiles:
        out[f"p{p}"] = float(np.percentile(arr, p))
    return out
