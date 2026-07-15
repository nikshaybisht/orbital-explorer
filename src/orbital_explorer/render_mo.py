"""3D view of a molecular orbital -- the "orbital sharing" picture.

For a chosen MO we place a p-orbital on every conjugated atom, weighted by that
atom's LCAO coefficient (size = |coefficient|, colour = sign), and add them up:

    psi_MO(r) = sum_i  c_i * p_z,i(r)

Rendering the isosurface of that sum shows exactly how the atomic orbitals share:
where neighbouring lobes have the same phase they merge into a bonding region;
where they have opposite phase a node appears between the atoms (antibonding).
The molecular skeleton (atoms + bonds) is drawn underneath for orientation.

The p-orbital exponent is chosen for visual clarity (lobes large enough to show
overlap), so this is a faithful *qualitative* LCAO picture, not a quantitative
electron density.
"""

from __future__ import annotations

import numpy as np
import plotly.graph_objects as go
from rdkit import Chem, rdBase
from rdkit.Chem import AllChem

from .huckel import HuckelResult
from .isosurface import phase_meshes

# CPK-ish colours and covalent radii (Angstrom), matching the convention used
# elsewhere in the user's tooling.
_CPK = {
    "H": "#e8e8e8", "C": "#404040", "N": "#3050f8", "O": "#ff2010",
    "S": "#e6c000", "P": "#ff8000", "F": "#50e050", "Cl": "#20d020",
}
_RADIUS = {"H": 0.31, "C": 0.76, "N": 0.71, "O": 0.66, "S": 1.05, "P": 1.07, "F": 0.57}


def embed(mol: Chem.Mol) -> Chem.Mol:
    """Add hydrogens and generate an optimised 3D conformer (indices preserved)."""
    with rdBase.BlockLogs():
        m = Chem.AddHs(Chem.Mol(mol))
        params = AllChem.ETKDGv3()
        params.randomSeed = 0xF00D
        if AllChem.EmbedMolecule(m, params) != 0:
            AllChem.EmbedMolecule(m, useRandomCoords=True, randomSeed=0xF00D)
        if m.GetNumConformers() == 0:
            raise ValueError("RDKit could not generate 3D coordinates for this molecule.")
        try:
            AllChem.MMFFOptimizeMolecule(m)
        except Exception:
            pass
    return m


def _coords(m: Chem.Mol) -> np.ndarray:
    conf = m.GetConformer()
    return np.array(
        [[conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y, conf.GetAtomPosition(i).z]
         for i in range(m.GetNumAtoms())]
    )


def _plane_fit(points: np.ndarray) -> tuple[np.ndarray, float]:
    """Best-fit plane through the conjugated atoms.

    Returns (unit normal, non-planarity). vh[-1] from the SVD is the direction
    of least spread (already unit length). Non-planarity is the smallest singular
    value relative to the largest: ~0 for a flat system, larger when atoms pucker
    out of plane, in which case a single p-orbital axis is only approximate.
    """
    centred = points - points.mean(axis=0)
    _, s, vh = np.linalg.svd(centred)
    non_planarity = float(s[-1] / (s[0] + 1e-9))
    return vh[-1], non_planarity


def _mo_field(grid: np.ndarray, centers: np.ndarray, coeffs, normal, zeta: float):
    """psi = sum_i c_i * (p_z lobe at atom i). grid is (..., 3)."""
    psi = np.zeros(grid.shape[:-1])
    for c, R in zip(coeffs, centers):
        d = grid - R
        along = d @ normal               # signed distance along the p-axis -> the two lobes
        rad = np.linalg.norm(d, axis=-1)
        psi += c * along * np.exp(-zeta * rad)
    return psi


def _skeleton_traces(m: Chem.Mol, coords: np.ndarray) -> list[go.Scatter3d]:
    bx, by, bz = [], [], []
    for bond in m.GetBonds():
        a, b = coords[bond.GetBeginAtomIdx()], coords[bond.GetEndAtomIdx()]
        bx += [a[0], b[0], None]
        by += [a[1], b[1], None]
        bz += [a[2], b[2], None]
    bonds = go.Scatter3d(x=bx, y=by, z=bz, mode="lines",
                         line=dict(color="#9aa6bd", width=4), hoverinfo="skip", showlegend=False)

    syms = [a.GetSymbol() for a in m.GetAtoms()]
    atoms = go.Scatter3d(
        x=coords[:, 0], y=coords[:, 1], z=coords[:, 2], mode="markers",
        marker=dict(
            size=[max(6, _RADIUS.get(s, 0.7) * 16) for s in syms],
            color=[_CPK.get(s, "#ff60c0") for s in syms],
            line=dict(color="#000", width=0.5),
        ),
        text=syms, hoverinfo="text", showlegend=False,
    )
    return [bonds, atoms]


def mo_figure(
    mol: Chem.Mol,
    result: HuckelResult,
    mo_index: int,
    n_grid: int = 76,
    zeta: float = 1.1,
    iso_fraction: float = 0.12,
    margin: float = 3.0,
) -> go.Figure:
    """Interactive 3D figure of molecular orbital `mo_index` for this molecule."""
    if len(result.pi_atom_indices) < 2:
        raise ValueError("Need at least two conjugated atoms to draw a molecular orbital.")
    mo_index = max(0, min(int(mo_index), len(result.orbitals) - 1))

    m = embed(mol)
    coords = _coords(m)
    centers = coords[result.pi_atom_indices]
    normal, non_planarity = _plane_fit(centers)
    orbital = result.orbitals[mo_index]

    lo, hi = coords.min(axis=0) - margin, coords.max(axis=0) + margin
    gx = np.linspace(lo[0], hi[0], n_grid)
    gy = np.linspace(lo[1], hi[1], n_grid)
    gz = np.linspace(lo[2], hi[2], n_grid)
    X, Y, Z = np.meshgrid(gx, gy, gz, indexing="ij")
    grid = np.stack([X, Y, Z], axis=-1)

    psi = _mo_field(grid, centers, orbital.coefficients, normal, zeta)
    peak = float(np.abs(psi).max()) or 1.0
    lobes = phase_meshes(psi, (gx, gy, gz), iso=iso_fraction * peak, opacity=0.3)

    fig = go.Figure(data=[*lobes, *_skeleton_traces(m, coords)])
    _style(fig, result, orbital, non_planarity)
    return fig


def _character(orbital, result: HuckelResult) -> str:
    if result.homo is not None and orbital.index == result.homo:
        return "HOMO"
    if result.lumo is not None and orbital.index == result.lumo:
        return "LUMO"
    if orbital.x > 1e-6:
        return "bonding"
    if orbital.x < -1e-6:
        return "antibonding"
    return "nonbonding"


def _style(fig: go.Figure, result: HuckelResult, orbital, non_planarity: float = 0.0) -> None:
    occ = {0: "empty", 1: "singly occupied", 2: "filled"}[orbital.occupation]
    title = (
        f"MO ψ{orbital.index + 1} of {len(result.orbitals)}: "
        f"{orbital.energy_label}  ({_character(orbital, result)}, {occ})"
    )
    axis_kw = dict(showbackground=False, showticklabels=False, title="",
                   gridcolor="#1c2330", zerolinecolor="#1c2330")
    annotations = []
    if non_planarity > 0.15:  # atoms deviate noticeably from one plane
        annotations.append(dict(
            text="non-planar: p-orbital axis is approximate",
            showarrow=False, xref="paper", yref="paper", x=0.5, y=0,
            font=dict(color="#e0a85a", size=12),
        ))
    fig.update_layout(
        title=dict(text=title, x=0.5, font=dict(color="#e8edf6", size=16)),
        paper_bgcolor="#0b0e17",
        scene=dict(xaxis=axis_kw, yaxis=axis_kw, zaxis=axis_kw, aspectmode="data"),
        margin=dict(l=0, r=0, t=44, b=0),
        annotations=annotations,
    )


def save_html(mol: Chem.Mol, result: HuckelResult, mo_index: int, path: str, **kwargs) -> str:
    fig = mo_figure(mol, result, mo_index, **kwargs)
    fig.write_html(path, include_plotlyjs="cdn")
    return path
