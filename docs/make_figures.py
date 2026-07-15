"""Regenerate the static images embedded in README.md.

The app itself is interactive (Plotly in the browser); these PNGs are just
previews, rendered offline with matplotlib from the same wavefunction / Huckel
code so they always match the tool.

    python docs/make_figures.py   ->   docs/images/*.png
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from skimage import measure

SRC = Path(__file__).resolve().parent.parent / "src"
sys.path.insert(0, str(SRC))

from orbital_explorer import render_mo as RMO  # noqa: E402
from orbital_explorer.atomic import Orbital, suggested_extent, wavefunction  # noqa: E402
from orbital_explorer.diatomic import solve_diatomic  # noqa: E402
from orbital_explorer.diatomic_diagram import diatomic_diagram  # noqa: E402
from orbital_explorer.huckel import solve  # noqa: E402
from orbital_explorer.mo_diagram import mo_diagram  # noqa: E402
from orbital_explorer.molecule import resolve  # noqa: E402
from rdkit import Chem  # noqa: E402

OUT = Path(__file__).resolve().parent / "images"
OUT.mkdir(parents=True, exist_ok=True)
BG, FG = "#0b0e17", "#e8edf6"
BLUE, GOLD = "#2f7bf6", "#f4c430"


def _add_mesh(ax, values, axes_origin, spacing, iso):
    for level, color in ((iso, BLUE), (-iso, GOLD)):
        try:
            verts, faces, _, _ = measure.marching_cubes(values, level=level, spacing=spacing)
        except (ValueError, RuntimeError):
            continue
        verts = verts + axes_origin
        pc = Poly3DCollection(verts[faces], alpha=0.45, linewidths=0)
        pc.set_facecolor(color)
        ax.add_collection3d(pc)


def _equal_3d(ax, half, center=(0, 0, 0)):
    cx, cy, cz = center
    ax.set_xlim(cx - half, cx + half)
    ax.set_ylim(cy - half, cy + half)
    ax.set_zlim(cz - half, cz + half)
    ax.set_box_aspect((1, 1, 1))
    ax.set_axis_off()


def atomic_montage():
    specs = [(2, "pz"), (3, "dz2"), (3, "dxy"), (4, "fxyz")]
    fig = plt.figure(figsize=(13, 3.7))
    fig.patch.set_facecolor(BG)
    for i, (n, label) in enumerate(specs, 1):
        orb = Orbital(n=n, label=label)
        ext = suggested_extent(orb)
        axis = np.linspace(-ext, ext, 72)
        X, Y, Z = np.meshgrid(axis, axis, axis, indexing="ij")
        psi = wavefunction(orb, X, Y, Z)
        ax = fig.add_subplot(1, 4, i, projection="3d")
        ax.set_facecolor(BG)
        _add_mesh(ax, psi, np.array([-ext, -ext, -ext]), (axis[1] - axis[0],) * 3,
                  0.10 * np.abs(psi).max())
        _equal_3d(ax, ext * 0.72)
        ax.set_title(orb.name, color=FG, fontsize=14)
    fig.tight_layout()
    fig.savefig(OUT / "atomic_orbitals.png", dpi=130, facecolor=BG)
    plt.close(fig)


def mo_diagram_img():
    info = resolve("benzene")
    res = solve(info.mol, info.pi_atom_indices)
    fig = mo_diagram(res, title="benzene  (C6H6)", is_monocycle=True)
    fig.savefig(OUT / "mo_diagram_benzene.png", dpi=130, facecolor=fig.get_facecolor())
    plt.close(fig)


def mo_3d_img():
    info = resolve("benzene")
    res = solve(info.mol, info.pi_atom_indices)
    m = RMO.embed(info.mol)
    coords = RMO._coords(m)
    centers = coords[res.pi_atom_indices]
    normal, _ = RMO._plane_fit(centers)
    orb = res.orbitals[res.homo]

    lo, hi = coords.min(axis=0) - 3, coords.max(axis=0) + 3
    gx = np.linspace(lo[0], hi[0], 72)
    gy = np.linspace(lo[1], hi[1], 72)
    gz = np.linspace(lo[2], hi[2], 72)
    X, Y, Z = np.meshgrid(gx, gy, gz, indexing="ij")
    psi = RMO._mo_field(np.stack([X, Y, Z], axis=-1), centers, orb.coefficients, normal, 1.1)

    fig = plt.figure(figsize=(6.2, 5.6))
    fig.patch.set_facecolor(BG)
    ax = fig.add_subplot(111, projection="3d")
    ax.set_facecolor(BG)
    _add_mesh(ax, psi, np.array([gx[0], gy[0], gz[0]]),
              (gx[1] - gx[0], gy[1] - gy[0], gz[1] - gz[0]), 0.12 * np.abs(psi).max())
    for bond in m.GetBonds():
        a, b = coords[bond.GetBeginAtomIdx()], coords[bond.GetEndAtomIdx()]
        ax.plot([a[0], b[0]], [a[1], b[1]], [a[2], b[2]], color="#9aa6bd", lw=2)
    span = (coords.max(0) - coords.min(0)).max() / 2 + 2.5
    c = coords.mean(0)
    _equal_3d(ax, span, center=(c[0], c[1], c[2]))
    ax.set_title("benzene HOMO (pi orbital sharing)", color=FG, fontsize=13)
    fig.savefig(OUT / "mo_3d_benzene.png", dpi=130, facecolor=BG)
    plt.close(fig)


def diatomic_img():
    mol = Chem.MolFromSmiles("O=O")
    res = solve_diatomic(mol, Chem.GetFormalCharge(mol))
    fig = diatomic_diagram(res, title="O₂")
    fig.savefig(OUT / "diatomic_o2.png", dpi=130, facecolor=fig.get_facecolor())
    plt.close(fig)


if __name__ == "__main__":
    atomic_montage()
    mo_diagram_img()
    mo_3d_img()
    diatomic_img()
    print("wrote images to", OUT)
