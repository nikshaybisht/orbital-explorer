"""Qualitative molecular-orbital diagrams for diatomic molecules.

Hueckel theory only covers conjugated pi-systems, so the canonical MO-theory
examples -- H2, N2, O2, CO, HF -- need a different model. This module builds the
standard *qualitative* diatomic MO picture taught in every intro course:

* a fixed ladder of sigma/pi MOs from the atoms' valence s and p orbitals,
* the all-important s-p mixing switch (for Li2..N2 the pi(2p) pair sits BELOW
  sigma(2p); for O2, F2, Ne2 the sigma(2p) drops below the pi pair),
* Aufbau + Hund filling (so O2 comes out with two unpaired electrons -> it is
  correctly predicted paramagnetic),
* bond order, HOMO/LUMO and magnetism.

Electrons are counted from atomic numbers and the overall charge, not from RDKit
bond orders -- that is both the chemically correct way to populate an MO diagram
and far more robust for radicals (NO), ions (O2+) and noble gases (He2).

Scope: period 1-2 diatomics (homo- and heteronuclear) and H-X molecules where X
is in the table (HF, HCl, LiH). Anything else returns a result with a clear note
instead of raising.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from rdkit import Chem

# Valence atomic-orbital energies in eV (negative = bound). Standard
# orbital-energy / VSIP values; only the relative ordering matters qualitatively,
# but real numbers keep the side columns of the diagram defensible.
_AO_ENERGY: dict[str, tuple[int, dict[str, float]]] = {
    "H": (1, {"1s": -13.6}),
    "He": (2, {"1s": -24.6}),
    "Li": (3, {"2s": -5.4, "2p": -3.5}),
    "Be": (4, {"2s": -9.3, "2p": -6.0}),
    "B": (5, {"2s": -14.0, "2p": -8.3}),
    "C": (6, {"2s": -19.4, "2p": -10.7}),
    "N": (7, {"2s": -25.6, "2p": -13.2}),
    "O": (8, {"2s": -32.4, "2p": -15.9}),
    "F": (9, {"2s": -40.2, "2p": -18.7}),
    "Ne": (10, {"2s": -48.5, "2p": -21.6}),
    "Na": (11, {"3s": -5.1, "3p": -3.0}),
    "Cl": (17, {"3s": -25.3, "3p": -13.7}),
}
# Core (non-valence) electrons that we do not draw as MOs.
_CORE = {"H": 0, "He": 0, "Li": 2, "Be": 2, "B": 2, "C": 2, "N": 2, "O": 2,
         "F": 2, "Ne": 2, "Na": 10, "Cl": 10}
_PERIOD1 = {"H", "He"}
_PERIOD2 = {"Li", "Be", "B", "C", "N", "O", "F", "Ne"}

# MO ladders: (name, character, degeneracy), most stable first.
_TEMPLATE_1S = [("σ1s", "bond", 1), ("σ*1s", "anti", 1)]
_TEMPLATE_MIXED = [  # Li2..N2 (and CO, NO): pi(2p) BELOW sigma(2p)
    ("σ2s", "bond", 1), ("σ*2s", "anti", 1),
    ("π2p", "bond", 2), ("σ2p", "bond", 1),
    ("π*2p", "anti", 2), ("σ*2p", "anti", 1),
]
_TEMPLATE_UNMIXED = [  # O2, F2, Ne2: sigma(2p) BELOW pi(2p)
    ("σ2s", "bond", 1), ("σ*2s", "anti", 1),
    ("σ2p", "bond", 1), ("π2p", "bond", 2),
    ("π*2p", "anti", 2), ("σ*2p", "anti", 1),
]


@dataclass
class DiatomicMO:
    index: int            # 0 = most stable
    name: str             # "σ2s", "π2p", ...
    character: str        # "bond" | "anti" | "nonbond"
    degeneracy: int       # 1 or 2
    occupations: list[int] = field(default_factory=list)  # per component, 0/1/2

    @property
    def n_electrons(self) -> int:
        return sum(self.occupations)

    @property
    def capacity(self) -> int:
        return 2 * self.degeneracy


@dataclass
class DiatomicResult:
    symbols: tuple[str, str]
    is_homonuclear: bool
    charge: int
    n_valence_electrons: int
    orbitals: list[DiatomicMO]
    bond_order: float
    n_unpaired: int
    magnetism: str                       # "paramagnetic" | "diamagnetic"
    homo: int | None
    lumo: int | None
    mixing: bool                         # s-p mixing on (pi below sigma)?
    ao_levels: dict[str, list[tuple[str, float]]]  # symbol -> [(subshell, eV)]
    note: str = ""

    @property
    def supported(self) -> bool:
        return bool(self.orbitals)


def _sp_mixing(sym_a: str, sym_b: str) -> bool:
    """True => strong 2s-2p mixing => pi(2p) sits below sigma(2p).

    Standard qualitative rule: mixing is significant for the early period-2
    elements and turns off once both atoms are O or later. So Li2..N2, and the
    heteronuclear CO/NO, use the mixed order; O2/F2/Ne2 use the unmixed order.
    """
    late = {"O", "F", "Ne"}
    return not (sym_a in late and sym_b in late)


def _fill(orbitals: list[DiatomicMO], n_electrons: int) -> None:
    """Aufbau + Hund filling, in place. One electron into each degenerate
    component before pairing, so half-filled degenerate shells come out as
    unpaired electrons (the origin of O2's paramagnetism)."""
    remaining = n_electrons
    for mo in orbitals:
        if remaining <= 0:
            mo.occupations = [0] * mo.degeneracy
            continue
        if remaining >= mo.capacity:
            mo.occupations = [2] * mo.degeneracy
            remaining -= mo.capacity
            continue
        occ = [0] * mo.degeneracy
        for _ in range(2):  # pass 1: singles (Hund); pass 2: pair up
            for i in range(mo.degeneracy):
                if remaining <= 0:
                    break
                occ[i] += 1
                remaining -= 1
        mo.occupations = occ


def _hx_template(x: str) -> list[tuple[str, str, int]]:
    """MO ladder for an H-X molecule (H bonds to heavy atom X).

    p-block X (has p valence electrons): an X ns lone pair, the H-X sigma bond,
    a degenerate pi lone pair (X's other p orbitals), and sigma*. s-block X
    (Li/Be/Na, no p electrons): just the sigma bond and sigma*.
    """
    valence = _AO_ENERGY[x][0] - _CORE[x]
    if valence >= 3:  # p-block: B, C, N, O, F, ... (has occupied np)
        return [("ns(lp)", "nonbond", 1), ("σ", "bond", 1),
                ("π(lp)", "nonbond", 2), ("σ*", "anti", 1)]
    return [("σ", "bond", 1), ("σ*", "anti", 1)]


def _template(sym_a: str, sym_b: str) -> tuple[list[tuple[str, str, int]], bool] | None:
    """Choose the MO ladder. Returns (template, mixing) or None if unsupported."""
    pair = {sym_a, sym_b}
    if pair <= _PERIOD1:
        return _TEMPLATE_1S, _sp_mixing(sym_a, sym_b)
    if "H" in pair:  # H-X
        x = sym_b if sym_a == "H" else sym_a
        return _hx_template(x), False
    if {sym_a, sym_b} <= _PERIOD2:
        mixing = _sp_mixing(sym_a, sym_b)
        return (_TEMPLATE_MIXED if mixing else _TEMPLATE_UNMIXED), mixing
    return None  # period-3 heavy-heavy etc. -- not covered


def solve_diatomic(mol: Chem.Mol, charge: int = 0) -> DiatomicResult:
    """Build the qualitative MO diagram data for a diatomic molecule."""
    syms = [a.GetSymbol() for a in Chem.AddHs(mol).GetAtoms()]
    a, b = (syms + ["?", "?"])[:2]
    homonuclear = a == b
    blank = dict(symbols=(a, b), is_homonuclear=homonuclear, charge=charge,
                 n_valence_electrons=0, orbitals=[], bond_order=0.0, n_unpaired=0,
                 magnetism="diamagnetic", homo=None, lumo=None, mixing=False,
                 ao_levels={})

    if a not in _AO_ENERGY or b not in _AO_ENERGY:
        missing = a if a not in _AO_ENERGY else b
        return DiatomicResult(**blank, note=f"MO diagram covers H-Ne diatomics (plus HCl); {missing} is not supported yet.")

    chosen = _template(a, b)
    if chosen is None:
        return DiatomicResult(**blank, note="Diatomic diagrams currently cover period 1-2 atoms and H-X molecules.")
    template, mixing = chosen

    n_val = (_AO_ENERGY[a][0] - _CORE[a]) + (_AO_ENERGY[b][0] - _CORE[b]) - charge
    orbitals = [DiatomicMO(index=i, name=n, character=c, degeneracy=d)
                for i, (n, c, d) in enumerate(template)]
    _fill(orbitals, n_val)

    bonding = sum(mo.n_electrons for mo in orbitals if mo.character == "bond")
    antibonding = sum(mo.n_electrons for mo in orbitals if mo.character == "anti")
    bond_order = (bonding - antibonding) / 2
    n_unpaired = sum(1 for mo in orbitals for occ in mo.occupations if occ == 1)
    homo = max((mo.index for mo in orbitals if mo.n_electrons > 0), default=None)
    lumo = next((mo.index for mo in orbitals if mo.n_electrons == 0), None)

    ao_levels = {a: list(_AO_ENERGY[a][1].items())}
    if not homonuclear:
        ao_levels[b] = list(_AO_ENERGY[b][1].items())

    return DiatomicResult(
        symbols=(a, b), is_homonuclear=homonuclear, charge=charge,
        n_valence_electrons=n_val, orbitals=orbitals, bond_order=bond_order,
        n_unpaired=n_unpaired,
        magnetism="paramagnetic" if n_unpaired else "diamagnetic",
        homo=homo, lumo=lumo, mixing=mixing, ao_levels=ao_levels,
    )
