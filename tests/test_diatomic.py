"""Diatomic MO results checked against textbook chemistry.

The headline cases an intro course tests on: bond orders, the N2 vs O2 ordering
switch, and O2's paramagnetism (two unpaired electrons).
"""

import pytest
from rdkit import Chem

from orbital_explorer.diatomic import solve_diatomic


def _result(smiles):
    mol = Chem.MolFromSmiles(smiles) or Chem.MolFromSmiles(smiles, sanitize=False)
    return solve_diatomic(mol, Chem.GetFormalCharge(mol))


def _order(res):
    return [mo.name for mo in res.orbitals]


@pytest.mark.parametrize(
    "smiles,n_val,bond_order,unpaired,magnetism",
    [
        ("[H][H]", 2, 1.0, 0, "diamagnetic"),    # H2
        ("[He][He]", 4, 0.0, 0, "diamagnetic"),  # He2 (no bond)
        ("N#N", 10, 3.0, 0, "diamagnetic"),      # N2
        ("O=O", 12, 2.0, 2, "paramagnetic"),     # O2 -- the famous one
        ("FF", 14, 1.0, 0, "diamagnetic"),       # F2
        ("[C-]#[O+]", 10, 3.0, 0, "diamagnetic"),  # CO
        ("[O+]=[O]", 11, 2.5, 1, "paramagnetic"),  # O2+
        ("F", 8, 1.0, 0, "diamagnetic"),         # HF ("F" -> H-F)
    ],
)
def test_diatomic_properties(smiles, n_val, bond_order, unpaired, magnetism):
    res = _result(smiles)
    assert res.n_valence_electrons == n_val
    assert res.bond_order == bond_order
    assert res.n_unpaired == unpaired
    assert res.magnetism == magnetism


def test_o2_unpaired_electrons_are_in_pi_star():
    """O2's two unpaired electrons must sit in the degenerate pi*(2p)."""
    res = _result("O=O")
    pistar = next(mo for mo in res.orbitals if mo.name == "π*2p")
    assert pistar.occupations == [1, 1]


def test_n2_vs_o2_ordering_switch():
    """The s-p mixing switch: N2 has pi(2p) below sigma(2p); O2 the reverse."""
    n2, o2 = _order(_result("N#N")), _order(_result("O=O"))
    assert n2.index("π2p") < n2.index("σ2p")   # mixed order
    assert o2.index("σ2p") < o2.index("π2p")   # unmixed order


def test_homo_identity():
    assert _result("N#N").orbitals[_result("N#N").homo].name == "σ2p"
    o2 = _result("O=O")
    assert o2.orbitals[o2.homo].name == "π*2p"
    hf = _result("F")
    assert hf.orbitals[hf.homo].name == "π(lp)"  # nonbonding lone pair


def test_unsupported_atom_is_graceful():
    res = _result("[Fe]")  # not a diatomic input, but solve must not raise
    assert res.note  # some explanatory note, no exception
