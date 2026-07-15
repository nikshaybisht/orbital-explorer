"""Huckel results checked against the analytic textbook values.

Path graph P_N eigenvalues are 2cos(k*pi/(N+1)); cycle C_N are 2cos(2*pi*k/N).
If our matrix and solver are right, the eigenvalues, electron counts and
HOMO/LUMO filling must reproduce these exactly.
"""

import numpy as np
import pytest

from orbital_explorer.huckel import solve
from orbital_explorer.molecule import resolve


def _xs(smiles):
    info = resolve(smiles)
    res = solve(info.mol, info.pi_atom_indices)
    return res, sorted((o.x for o in res.orbitals), reverse=True)


def test_ethylene():
    res, xs = _xs("C=C")
    assert xs == pytest.approx([1.0, -1.0])
    assert res.n_pi_electrons == 2
    assert [o.occupation for o in res.orbitals] == [2, 0]
    assert res.gap_beta == pytest.approx(2.0)


def test_butadiene():
    res, xs = _xs("C=CC=C")
    assert xs == pytest.approx([1.618034, 0.618034, -0.618034, -1.618034], abs=1e-5)
    assert res.n_pi_electrons == 4
    assert [o.occupation for o in res.orbitals] == [2, 2, 0, 0]


def test_benzene_aromatic_sextet():
    res, xs = _xs("c1ccccc1")
    assert xs == pytest.approx([2, 1, 1, -1, -1, -2], abs=1e-6)
    assert res.n_pi_electrons == 6
    assert [o.occupation for o in res.orbitals] == [2, 2, 2, 0, 0, 0]
    # Degenerate HOMO pair and degenerate LUMO pair.
    assert res.orbitals[res.homo].x == pytest.approx(1.0)
    assert res.orbitals[res.lumo].x == pytest.approx(-1.0)


@pytest.mark.parametrize(
    "smiles,n_pi,occ_homo",
    [("[CH2]C=C", 3, 1), ("[CH2+]C=C", 2, 2), ("[CH2-]C=C", 4, 2)],
)
def test_allyl_charge_states(smiles, n_pi, occ_homo):
    res, _ = _xs(smiles)
    assert res.n_pi_electrons == n_pi
    assert res.orbitals[res.homo].occupation == occ_homo


def test_cyclobutadiene_is_triplet_diradical():
    res, xs = _xs("C1=CC=C1")
    assert xs == pytest.approx([2, 0, 0, -2], abs=1e-6)
    # Two singly-occupied nonbonding orbitals (Hund) -> antiaromatic diradical.
    singly = [o for o in res.orbitals if o.occupation == 1]
    assert len(singly) == 2
    assert all(abs(o.x) < 1e-6 for o in singly)


def test_cyclopentadienyl_anion_aromatic():
    res, _ = _xs("[CH-]1C=CC=C1")
    assert res.n_pi_electrons == 6
    assert res.charge == -1


def test_coefficients_normalised():
    # Each MO eigenvector is normalised; the set is orthonormal.
    info = resolve("c1ccccc1")
    res = solve(info.mol, info.pi_atom_indices)
    C = np.column_stack([o.coefficients for o in res.orbitals])
    assert np.allclose(C.T @ C, np.eye(C.shape[1]), atol=1e-6)
