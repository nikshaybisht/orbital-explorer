"""Exact hydrogenic atomic orbitals -- the analytic solutions of the Schrodinger
equation for a one-electron atom.

This is the foundation of the whole tool. Every shape in the Orbitron gallery
(1s, 2p_z, 3d_xy, 4f ...) is one of these functions:

    psi(r, theta, phi) = R_{n,l}(r) * Y_{l,m}(theta, phi)

We use the *real* spherical harmonics chemists actually draw and name
(p_x/p_y/p_z, d_xy, d_x2-y2, ...) rather than the complex Y_l^m, because the
real combinations are the ones with the familiar lobed shapes.

Lengths are in Bohr radii (a0 = 1) and the nuclear charge Z is adjustable, so a
hydrogen 1s (Z=1) and a He+ 1s (Z=2) come out correctly scaled. Everything is
normalised so that the integral of |psi|^2 over all space is 1.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy.special import factorial, genlaguerre

# --------------------------------------------------------------------------- #
# Real spherical harmonics, written directly as polynomials in the components
# of the unit direction vector (ux, uy, uz) = (x/r, y/r, z/r).
#
# Each entry maps a chemist's orbital label to (l, function). The functions are
# normalised so that the integral of Y^2 over the unit sphere equals 1, which
# means |psi|^2 integrates to 1 over all space (the radial part is normalised
# separately below). Coefficients are the standard real-harmonic constants,
# e.g. p uses sqrt(3/4pi), so they are exact, not fitted.
# --------------------------------------------------------------------------- #

_PI = math.pi


def _c(value: float) -> float:
    return math.sqrt(value)


ANGULAR: dict[str, tuple[int, "callable"]] = {
    # l = 0
    "s": (0, lambda ux, uy, uz: np.full_like(ux, _c(1 / (4 * _PI)))),
    # l = 1
    "px": (1, lambda ux, uy, uz: _c(3 / (4 * _PI)) * ux),
    "py": (1, lambda ux, uy, uz: _c(3 / (4 * _PI)) * uy),
    "pz": (1, lambda ux, uy, uz: _c(3 / (4 * _PI)) * uz),
    # l = 2
    "dz2": (2, lambda ux, uy, uz: _c(5 / (16 * _PI)) * (3 * uz**2 - 1)),
    "dxz": (2, lambda ux, uy, uz: _c(15 / (4 * _PI)) * ux * uz),
    "dyz": (2, lambda ux, uy, uz: _c(15 / (4 * _PI)) * uy * uz),
    "dxy": (2, lambda ux, uy, uz: _c(15 / (4 * _PI)) * ux * uy),
    "dx2-y2": (2, lambda ux, uy, uz: _c(15 / (16 * _PI)) * (ux**2 - uy**2)),
    # l = 3
    "fz3": (3, lambda ux, uy, uz: _c(7 / (16 * _PI)) * uz * (5 * uz**2 - 3)),
    "fxz2": (3, lambda ux, uy, uz: _c(21 / (32 * _PI)) * ux * (5 * uz**2 - 1)),
    "fyz2": (3, lambda ux, uy, uz: _c(21 / (32 * _PI)) * uy * (5 * uz**2 - 1)),
    "fxyz": (3, lambda ux, uy, uz: _c(105 / (4 * _PI)) * ux * uy * uz),
    "fz(x2-y2)": (3, lambda ux, uy, uz: _c(105 / (16 * _PI)) * uz * (ux**2 - uy**2)),
    "fx(x2-3y2)": (3, lambda ux, uy, uz: _c(35 / (32 * _PI)) * ux * (ux**2 - 3 * uy**2)),
    "fy(3x2-y2)": (3, lambda ux, uy, uz: _c(35 / (32 * _PI)) * uy * (3 * ux**2 - uy**2)),
}

# Orbital labels grouped by subshell, in the order chemists list them.
_BY_L: dict[int, list[str]] = {
    0: ["s"],
    1: ["px", "py", "pz"],
    2: ["dz2", "dxz", "dyz", "dxy", "dx2-y2"],
    3: ["fz3", "fxz2", "fyz2", "fxyz", "fz(x2-y2)", "fx(x2-3y2)", "fy(3x2-y2)"],
}

_SUBSHELL = {0: "s", 1: "p", 2: "d", 3: "f"}


@dataclass(frozen=True)
class Orbital:
    """A single hydrogenic orbital, e.g. Orbital(n=3, label='dxy', Z=1)."""

    n: int
    label: str
    Z: float = 1.0

    def __post_init__(self) -> None:
        if self.label not in ANGULAR:
            raise ValueError(f"unknown orbital label {self.label!r}")
        if self.n < 1:
            raise ValueError("principal quantum number n must be >= 1")
        if self.l >= self.n:
            raise ValueError(
                f"invalid orbital: l={self.l} requires n>{self.l}, got n={self.n}"
            )

    @property
    def l(self) -> int:  # noqa: E743 - matches physics notation
        return ANGULAR[self.label][0]

    @property
    def name(self) -> str:
        """Human label like '3dxy' or '2pz'."""
        return f"{self.n}{self.label}"

    @property
    def radial_nodes(self) -> int:
        return self.n - self.l - 1

    @property
    def angular_nodes(self) -> int:
        return self.l


def available_orbitals(n: int) -> list[str]:
    """Orbital labels allowed for principal quantum number n (l = 0..n-1)."""
    labels: list[str] = []
    for l in range(min(n, 4)):  # f-orbitals (l=3) are the highest we draw
        labels.extend(_BY_L[l])
    return labels


def radial(n: int, l: int, r: np.ndarray, Z: float = 1.0) -> np.ndarray:
    """Radial wavefunction R_{n,l}(r), normalised so int R^2 r^2 dr = 1.

    Built from the associated Laguerre polynomial, exactly as in the analytic
    hydrogen solution; no approximation.
    """
    r = np.asarray(r, dtype=float)
    rho = 2.0 * Z * r / n
    norm = math.sqrt((2.0 * Z / n) ** 3 * factorial(n - l - 1) / (2.0 * n * factorial(n + l)))
    laguerre = genlaguerre(n - l - 1, 2 * l + 1)(rho)
    return norm * np.exp(-rho / 2.0) * rho**l * laguerre


def _unit_vectors(x, y, z, r):
    """x/r, y/r, z/r with the r=0 singularity set to 0 (safe: |psi|=0 there for l>=1)."""
    with np.errstate(divide="ignore", invalid="ignore"):
        ux = np.where(r > 0, x / r, 0.0)
        uy = np.where(r > 0, y / r, 0.0)
        uz = np.where(r > 0, z / r, 0.0)
    return ux, uy, uz


def wavefunction(orbital: Orbital, x, y, z) -> np.ndarray:
    """Evaluate the real wavefunction psi on Cartesian coordinate arrays.

    Sign is meaningful: positive and negative lobes are the orbital's phase,
    which is what gets coloured differently in the visualisation.
    """
    x, y, z = np.asarray(x, float), np.asarray(y, float), np.asarray(z, float)
    r = np.sqrt(x**2 + y**2 + z**2)
    ux, uy, uz = _unit_vectors(x, y, z, r)
    l, angular_fn = ANGULAR[orbital.label]
    return radial(orbital.n, l, r, orbital.Z) * angular_fn(ux, uy, uz)


def probability_density(orbital: Orbital, x, y, z) -> np.ndarray:
    """|psi|^2 -- the probability density of finding the electron."""
    psi = wavefunction(orbital, x, y, z)
    return psi**2


def suggested_extent(orbital: Orbital) -> float:
    """A reasonable half-width (in Bohr radii) for a cubic sampling box.

    Sized to enclose ~99% of the density while staying tight enough that the
    orbital fills the view (the old n^2 rule left high-n orbitals tiny in a huge
    empty box). Scales like n(n+1.5)/Z, matching how <r> grows with n.
    """
    return 2.2 * orbital.n * (orbital.n + 1.5) / orbital.Z + 3.0
