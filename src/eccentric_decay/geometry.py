"""
geometry.py — Eccentric decay cross-section models for tree stems.

The outer boundary is assumed to remain a circle of radius R (intact bark/wood).
The decay cavity may be:
  - concentric circular  (e = 0; recovers Mattheck 1994 model in the e→0 limit)
  - eccentric circular   (offset center)
  - elliptical           (eccentric + non-isotropic)

All section properties are returned with respect to the NET-section centroid.

Conventions
-----------
  Coordinates : (x, y) lie in the cross-sectional plane perpendicular to the stem axis.
                The outer disk is centered at the origin in pre-shift coordinates;
                the cavity may be offset.
  Areas       : positive scalar for cavities; subtracted from the outer disk.
  Inertias    : I_xx = ∫ y² dA, I_yy = ∫ x² dA, I_xy = ∫ x y dA  (engineering convention).
  Returned    : I_xx, I_yy, I_xy about the NET centroid (post-shift).

Reference
---------
  Mattheck, C. & Breloer, H. (1994). The Body Language of Trees. HMSO, London.
  Boresi, A.P. & Schmidt, R.J. (2003). Advanced Mechanics of Materials, 6th ed.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Tuple

import numpy as np


# ─────────────────────────────────────────────────────────────────────
# Cavity primitives
# ─────────────────────────────────────────────────────────────────────

@dataclass
class CircularCavity:
    """Circular decay cavity, possibly offset from the stem center."""
    r_d: float                 # cavity radius [m]
    e_x: float = 0.0           # x-offset of cavity center [m]
    e_y: float = 0.0           # y-offset of cavity center [m]

    def area(self) -> float:
        return np.pi * self.r_d ** 2

    def centroid(self) -> Tuple[float, float]:
        return (self.e_x, self.e_y)

    def moments_about_origin(self) -> Tuple[float, float, float]:
        """(I_xx, I_yy, I_xy) of the cavity area about the stem origin."""
        I_c = np.pi * self.r_d ** 4 / 4              # centroidal (both axes equal)
        A   = self.area()
        I_xx = I_c + A * self.e_y ** 2               # parallel axis
        I_yy = I_c + A * self.e_x ** 2
        I_xy = A * self.e_x * self.e_y               # I_xy_centroidal = 0 for a circle
        return I_xx, I_yy, I_xy


@dataclass
class EllipticalCavity:
    """Elliptical cavity with semi-axes (a, b), center (e_x, e_y), orientation phi."""
    a: float                   # semi-axis 1 [m]
    b: float                   # semi-axis 2 [m]
    e_x: float = 0.0
    e_y: float = 0.0
    phi: float = 0.0           # CCW rotation of the a-axis from x-axis [rad]

    def area(self) -> float:
        return np.pi * self.a * self.b

    def centroid(self) -> Tuple[float, float]:
        return (self.e_x, self.e_y)

    def moments_about_origin(self) -> Tuple[float, float, float]:
        # Centroidal moments in the local (a, b) frame
        I_aa = np.pi * self.a * self.b ** 3 / 4     # about the a-axis (depends on b)
        I_bb = np.pi * self.a ** 3 * self.b / 4     # about the b-axis (depends on a)
        # Rotate to (x, y) frame (Boresi & Schmidt eq. 8-7)
        c, s = np.cos(self.phi), np.sin(self.phi)
        I_xx_c = I_aa * c ** 2 + I_bb * s ** 2
        I_yy_c = I_aa * s ** 2 + I_bb * c ** 2
        I_xy_c = (I_aa - I_bb) * s * c
        # Parallel axis to stem origin
        A = self.area()
        I_xx = I_xx_c + A * self.e_y ** 2
        I_yy = I_yy_c + A * self.e_x ** 2
        I_xy = I_xy_c + A * self.e_x * self.e_y
        return I_xx, I_yy, I_xy


# ─────────────────────────────────────────────────────────────────────
# Composite cross-section
# ─────────────────────────────────────────────────────────────────────

@dataclass
class DecaySection:
    """
    Tree cross-section: intact circular outer of radius R minus a decay cavity.
    """
    R: float                                     # outer radius [m]
    cavity: Optional[object] = None              # CircularCavity, EllipticalCavity, or None

    def section_properties(self) -> dict:
        """
        Net-section properties about the NET centroid:
            A, x_c, y_c, I_xx, I_yy, I_xy

        Notes
        -----
        Outer disk contributes A_out = π R², I_xx_out = I_yy_out = π R⁴/4 (about origin),
        and zero first moment. Cavity is subtracted via the negative-area trick.
        """
        A_out = np.pi * self.R ** 2
        I_out = np.pi * self.R ** 4 / 4

        if self.cavity is None:
            return dict(A=A_out, x_c=0.0, y_c=0.0,
                        I_xx=I_out, I_yy=I_out, I_xy=0.0)

        A_c = self.cavity.area()
        e_x, e_y = self.cavity.centroid()
        I_cxx, I_cyy, I_cxy = self.cavity.moments_about_origin()

        if A_c >= A_out:
            raise ValueError(f"Cavity area ({A_c:.4f}) ≥ outer area ({A_out:.4f})")

        A_net = A_out - A_c

        # Centroid of net section (using composite-area rule)
        x_c = -A_c * e_x / A_net
        y_c = -A_c * e_y / A_net

        # Second moments about the stem origin (outer minus cavity)
        I_xx_o = I_out - I_cxx
        I_yy_o = I_out - I_cyy
        I_xy_o = 0.0   - I_cxy

        # Shift to the NET centroid (parallel axis: I_centroid = I_origin - A·d²)
        I_xx = I_xx_o - A_net * y_c ** 2
        I_yy = I_yy_o - A_net * x_c ** 2
        I_xy = I_xy_o - A_net * x_c * y_c

        return dict(A=A_net, x_c=x_c, y_c=y_c,
                    I_xx=I_xx, I_yy=I_yy, I_xy=I_xy)

    # convenience constructors -------------------------------------------------
    @classmethod
    def concentric(cls, R: float, decay_ratio: float) -> "DecaySection":
        """Mattheck-style concentric circular cavity, r_d = decay_ratio · R."""
        return cls(R=R, cavity=CircularCavity(r_d=decay_ratio * R))

    @classmethod
    def eccentric_circular(cls, R: float, decay_ratio: float,
                           ecc_norm: float, ecc_angle: float = 0.0) -> "DecaySection":
        """
        Eccentric circular cavity.

        Parameters
        ----------
        R           : outer radius [m]
        decay_ratio : r_d / R (cavity radius normalized)
        ecc_norm    : geometric eccentricity normalized so that ecc_norm = 1
                      places the cavity tangent to the outer wall (max physical eccentricity).
                      i.e. e = ecc_norm · (R - r_d).
        ecc_angle   : angle of the offset direction from x-axis [rad].
        """
        r_d = decay_ratio * R
        e   = ecc_norm * (R - r_d)
        return cls(R=R, cavity=CircularCavity(
            r_d=r_d, e_x=e * np.cos(ecc_angle), e_y=e * np.sin(ecc_angle)))

    @classmethod
    def elliptical(cls, R: float, a_ratio: float, b_ratio: float,
                   ecc_norm: float = 0.0, ecc_angle: float = 0.0,
                   phi: float = 0.0) -> "DecaySection":
        """
        Elliptical cavity with semi-axes (a_ratio · R, b_ratio · R).

        Eccentricity is interpreted as in `eccentric_circular`, where the equivalent
        "radius" used to bound the offset is max(a, b), so ecc_norm = 1 places the
        ellipse tangent to the outer wall along its major axis.
        """
        a = a_ratio * R
        b = b_ratio * R
        bound = max(a, b)
        e = ecc_norm * (R - bound)
        return cls(R=R, cavity=EllipticalCavity(
            a=a, b=b,
            e_x=e * np.cos(ecc_angle), e_y=e * np.sin(ecc_angle),
            phi=phi))
