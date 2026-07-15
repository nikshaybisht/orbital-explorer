"""Tests that the hydrogenic orbitals obey the physics they must obey.

We don't compare against a reference image; we check the defining mathematical
properties: orbitals are normalised, mutually orthogonal, and have the right
number of nodes. If those hold, the shapes are necessarily correct.
"""

import numpy as np
import pytest

from orbital_explorer.atomic import (
    Orbital,
    available_orbitals,
    radial,
    suggested_extent,
    wavefunction,
)


def _grid(orb, n=130, shrink=1.0):
    L = suggested_extent(orb) * shrink
    axis = np.linspace(-L, L, n)
    dx = axis[1] - axis[0]
    X, Y, Z = np.meshgrid(axis, axis, axis, indexing="ij")
    return X, Y, Z, dx


@pytest.mark.parametrize("n,label", [(1, "s"), (2, "s"), (2, "pz"), (3, "dz2"), (3, "dxy")])
def test_normalised(n, label):
    orb = Orbital(n=n, label=label)
    X, Y, Z, dx = _grid(orb)
    integral = np.sum(wavefunction(orb, X, Y, Z) ** 2) * dx**3
    assert integral == pytest.approx(1.0, abs=0.03)


@pytest.mark.parametrize(
    "a,b",
    [(("2", "pz"), ("2", "s")), (("2", "px"), ("2", "py")), (("3", "dxy"), ("3", "dz2"))],
)
def test_orthogonal(a, b):
    oa = Orbital(n=int(a[0]), label=a[1])
    ob = Orbital(n=int(b[0]), label=b[1])
    X, Y, Z, dx = _grid(oa)
    overlap = np.sum(wavefunction(oa, X, Y, Z) * wavefunction(ob, X, Y, Z)) * dx**3
    assert overlap == pytest.approx(0.0, abs=0.02)


@pytest.mark.parametrize(
    "n,l,expected", [(1, 0, 0), (2, 0, 1), (3, 0, 2), (4, 0, 3), (2, 1, 0), (3, 1, 1), (3, 2, 0)]
)
def test_radial_node_count(n, l, expected):
    r = np.linspace(1e-4, 80, 6000)
    R = radial(n, l, r)
    sign_changes = int(np.sum(np.diff(np.sign(R)) != 0))
    assert sign_changes == expected


def test_node_properties():
    orb = Orbital(n=4, label="dxy")
    assert orb.l == 2
    assert orb.radial_nodes == 4 - 2 - 1
    assert orb.angular_nodes == 2
    assert orb.name == "4dxy"


def test_invalid_orbital_rejected():
    with pytest.raises(ValueError):
        Orbital(n=1, label="pz")  # l=1 needs n>=2
    with pytest.raises(ValueError):
        Orbital(n=2, label="not-an-orbital")


def test_available_orbitals():
    assert available_orbitals(1) == ["s"]
    assert available_orbitals(2) == ["s", "px", "py", "pz"]
    assert "dxy" in available_orbitals(3)
    assert "fxyz" in available_orbitals(4)
