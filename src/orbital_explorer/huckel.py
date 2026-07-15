"""Huckel molecular-orbital theory for conjugated pi-systems.

Huckel theory is the simplest model that still captures the essential physics of
a conjugated pi-system: it treats each p-orbital as a basis function, puts the
atom's energy (alpha) on the diagonal and the neighbour interaction (beta) on
bonds, then diagonalises. The eigenvalues are MO energies, written E = alpha + x*beta
(beta is negative, so the *largest* x is the most stable orbital). The
eigenvectors are the LCAO coefficients -- literally how much each atom's
p-orbital contributes to each molecular orbital, i.e. the "orbital sharing".

It's exact for the textbook cases (ethylene, allyl, butadiene, benzene, ...) and
qualitatively right for heteroatom systems via the standard Streitwieser
parameters. We deliberately keep it light: no SCF, no basis sets, runs instantly,
and the numbers match what a student computes by hand.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from rdkit import Chem

# Streitwieser-style heteroatom parameters, in units of beta.
# alpha_X = alpha + h * beta ;  beta_CX = k * beta.   (carbon: h=0, k=1)
_HETERO = {
    # symbol, "type": (h, pi_electrons_donated)
    ("N", "pyridine"): (0.5, 1),   # =N- , two connections
    ("N", "pyrrole"): (1.5, 2),    # -NH-/-NR-, three connections (donates lone pair)
    ("O", "furan"): (2.0, 2),      # -O- in a ring (donates lone pair)
    ("O", "carbonyl"): (1.0, 1),   # =O
    ("S", "thiophene"): (1.0, 2),
}
_K_CX = {"N": 0.8, "O": 0.8, "S": 0.7}

_DEGEN_TOL = 1e-6


@dataclass
class MolecularOrbital:
    index: int            # 0 = most stable
    x: float              # energy as E = alpha + x*beta
    energy_label: str     # e.g. "alpha + 1.618 beta"
    coefficients: np.ndarray  # LCAO coefficient per pi atom
    occupation: int       # 0, 1 or 2 electrons


@dataclass
class HuckelResult:
    pi_atom_indices: list[int]      # RDKit atom indices, in matrix order
    symbols: list[str]
    orbitals: list[MolecularOrbital]  # ordered most stable -> least stable
    n_pi_electrons: int
    homo: int | None                # index into orbitals
    lumo: int | None
    charge: int

    @property
    def gap_beta(self) -> float | None:
        """HOMO-LUMO gap in units of |beta| (positive)."""
        if self.homo is None or self.lumo is None:
            return None
        return self.orbitals[self.homo].x - self.orbitals[self.lumo].x


def _atom_parameters(atom: Chem.Atom) -> tuple[float, int]:
    """Return (h, pi_electrons_donated) for an atom in the pi framework."""
    sym = atom.GetSymbol()
    if sym == "C":
        return 0.0, 1
    deg = atom.GetTotalDegree()
    if sym == "N":
        kind = "pyrrole" if deg >= 3 else "pyridine"
        return _HETERO[("N", kind)]
    if sym == "O":
        kind = "carbonyl" if deg == 1 else "furan"
        return _HETERO[("O", kind)]
    if sym == "S":
        return _HETERO[("S", "thiophene")]
    return 0.0, 1


def build_matrix(mol: Chem.Mol, pi_atoms: list[int]) -> tuple[np.ndarray, int]:
    """Build the Huckel matrix (in beta units) and the pi-electron count.

    Returns (M, n_pi_electrons). M is symmetric: diagonal = h, off-diagonal = k
    for bonded pi atoms.
    """
    n = len(pi_atoms)
    idx_of = {a: i for i, a in enumerate(pi_atoms)}
    M = np.zeros((n, n))

    donated = 0
    for a in pi_atoms:
        atom = mol.GetAtomWithIdx(a)
        h, e = _atom_parameters(atom)
        M[idx_of[a], idx_of[a]] = h
        donated += e

    for bond in mol.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if i in idx_of and j in idx_of:
            si = mol.GetAtomWithIdx(i).GetSymbol()
            sj = mol.GetAtomWithIdx(j).GetSymbol()
            hetero = {si, sj} - {"C"}
            k = _K_CX.get(next(iter(hetero)), 1.0) if hetero else 1.0
            M[idx_of[i], idx_of[j]] = M[idx_of[j], idx_of[i]] = k

    n_pi = donated - Chem.GetFormalCharge(mol)
    return M, n_pi


def _fill(energies_x: np.ndarray, n_electrons: int) -> list[int]:
    """Aufbau + Hund occupation across (possibly degenerate) levels.

    energies_x is ordered most stable first (descending x). Returns occupation
    (0/1/2) per orbital.
    """
    n = len(energies_x)
    occ = [0] * n
    remaining = n_electrons
    i = 0
    while remaining > 0 and i < n:
        # Gather a degenerate group.
        j = i
        while j + 1 < n and abs(energies_x[j + 1] - energies_x[i]) < _DEGEN_TOL:
            j += 1
        group = list(range(i, j + 1))
        capacity = 2 * len(group)
        if remaining >= capacity:
            for g in group:
                occ[g] = 2
            remaining -= capacity
        else:
            # One electron each first (Hund), then pair up.
            for g in group:
                if remaining <= 0:
                    break
                occ[g] += 1
                remaining -= 1
            for g in group:
                if remaining <= 0:
                    break
                occ[g] += 1
                remaining -= 1
        i = j + 1
    return occ


def solve(mol: Chem.Mol, pi_atoms: list[int]) -> HuckelResult:
    """Diagonalise the Huckel matrix and assign occupations / HOMO / LUMO."""
    M, n_pi = build_matrix(mol, pi_atoms)
    eigvals, eigvecs = np.linalg.eigh(M)  # ascending eigenvalues

    # Most stable orbital = largest x (because beta < 0). Reorder descending.
    order = np.argsort(eigvals)[::-1]
    xs = eigvals[order]
    vecs = eigvecs[:, order]

    occ = _fill(xs, n_pi)

    orbitals: list[MolecularOrbital] = []
    for k, x in enumerate(xs):
        orbitals.append(
            MolecularOrbital(
                index=k,
                x=float(x),
                energy_label=_energy_label(x),
                coefficients=vecs[:, k],
                occupation=occ[k],
            )
        )

    homo = max((k for k in range(len(orbitals)) if occ[k] > 0), default=None)
    lumo = min((k for k in range(len(orbitals)) if occ[k] == 0), default=None)

    symbols = [mol.GetAtomWithIdx(a).GetSymbol() for a in pi_atoms]
    return HuckelResult(
        pi_atom_indices=list(pi_atoms),
        symbols=symbols,
        orbitals=orbitals,
        n_pi_electrons=n_pi,
        homo=homo,
        lumo=lumo,
        charge=Chem.GetFormalCharge(mol),
    )


def _energy_label(x: float) -> str:
    if abs(x) < _DEGEN_TOL:
        return "α"  # alpha (nonbonding)
    sign = "+" if x > 0 else "−"
    return f"α {sign} {abs(x):.3f} β"
