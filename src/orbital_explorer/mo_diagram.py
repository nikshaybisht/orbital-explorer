"""Molecular-orbital energy-level diagram and the numbers that go with it.

Draws the classic vertical MO ladder for a Huckel result: bonding orbitals below
the alpha line, antibonding above, electrons shown as up/down arrows, and the
HOMO and LUMO called out with the gap between them. `mo_summary` returns the
quantitative facts a student needs: pi-electron count, total pi-energy,
delocalisation (resonance) energy, HOMO/LUMO energies and the gap.
"""

from __future__ import annotations

from dataclasses import dataclass

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from .huckel import HuckelResult

_BG = "#0b0e17"
_FG = "#e8edf6"
_BOND = "#5aa0ff"
_ANTI = "#ff7676"
_NONBOND = "#c9cf6a"


@dataclass
class MOSummary:
    n_pi_electrons: int
    total_pi_energy_beta: float       # coefficient of beta in E_pi = n*alpha + (this)*beta
    delocalization_energy_beta: float | None
    homo_x: float | None
    lumo_x: float | None
    gap_beta: float | None
    aromaticity: str | None           # "aromatic (4n+2)", "antiaromatic (4n)", or None


def mo_summary(result: HuckelResult, is_monocycle: bool = False) -> MOSummary:
    """Compute the headline quantities for a Huckel result."""
    total = sum(o.occupation * o.x for o in result.orbitals)

    # Delocalisation energy vs the same number of isolated ethylene pi bonds
    # (each filled ethylene pi contributes 2*beta). Only meaningful closed-shell.
    open_shell = any(o.occupation == 1 for o in result.orbitals)
    deloc = None
    if not open_shell:
        reference = result.n_pi_electrons  # (n/2 bonds) * 2 electrons * x=1
        deloc = total - reference

    aromaticity = None
    if is_monocycle and result.n_pi_electrons > 0:
        m = result.n_pi_electrons
        if m % 4 == 2:
            aromaticity = f"aromatic (4n+2, n={(m - 2) // 4})"
        elif m % 4 == 0:
            aromaticity = f"antiaromatic (4n, n={m // 4})"

    return MOSummary(
        n_pi_electrons=result.n_pi_electrons,
        total_pi_energy_beta=total,
        delocalization_energy_beta=deloc,
        homo_x=result.orbitals[result.homo].x if result.homo is not None else None,
        lumo_x=result.orbitals[result.lumo].x if result.lumo is not None else None,
        gap_beta=result.gap_beta,
        aromaticity=aromaticity,
    )


def _group_levels(result: HuckelResult, tol: float = 1e-6):
    """Group orbital indices by (near-)equal energy, most stable first."""
    groups: list[list[int]] = []
    for o in result.orbitals:
        if groups and abs(result.orbitals[groups[-1][0]].x - o.x) < tol:
            groups[-1].append(o.index)
        else:
            groups.append([o.index])
    return groups


def _draw_electrons(ax, xc: float, y: float, occ: int) -> None:
    dy = 0.16
    arrow = dict(arrowstyle="-|>", lw=1.6, mutation_scale=11, color=_FG)
    if occ == 1:
        ax.annotate("", xy=(xc, y + dy), xytext=(xc, y - dy), arrowprops=arrow)
    elif occ == 2:
        ax.annotate("", xy=(xc - 0.07, y + dy), xytext=(xc - 0.07, y - dy), arrowprops=arrow)
        ax.annotate("", xy=(xc + 0.07, y - dy), xytext=(xc + 0.07, y + dy), arrowprops=arrow)


def mo_diagram(result: HuckelResult, title: str = "", is_monocycle: bool = False):
    """Build a matplotlib Figure of the MO energy-level diagram."""
    fig, ax = plt.subplots(figsize=(6.4, 7.2))
    fig.patch.set_facecolor(_BG)
    ax.set_facecolor(_BG)

    groups = _group_levels(result)
    seg_w = 0.5
    spacing = 0.7

    # alpha (nonbonding) reference line.
    ax.axhline(0, color="#39507a", lw=1, ls="--")
    ax.text(2.7, 0, "α", color="#7a8aa8", va="center", fontsize=12)

    for group in groups:
        x_energy = result.orbitals[group[0]].x
        y = -x_energy  # bonding (x>0) sits below alpha
        color = _BOND if x_energy > 1e-6 else _ANTI if x_energy < -1e-6 else _NONBOND
        n = len(group)
        centers = [(i - (n - 1) / 2) * spacing for i in range(n)]
        for idx, xc in zip(group, centers):
            o = result.orbitals[idx]
            ax.plot([xc - seg_w / 2, xc + seg_w / 2], [y, y], color=color, lw=3, solid_capstyle="round")
            _draw_electrons(ax, xc, y, o.occupation)
            if idx == result.homo:
                ax.text(xc, y - 0.32, "HOMO", color=_BOND, ha="center", fontsize=9, weight="bold")
            if idx == result.lumo:
                ax.text(xc, y + 0.30, "LUMO", color=_ANTI, ha="center", fontsize=9, weight="bold")
        ax.text(2.7, y, result.orbitals[group[0]].energy_label, color=_FG, va="center", fontsize=10)

    # HOMO-LUMO gap arrow (sits in the empty space left of the levels).
    if result.homo is not None and result.lumo is not None:
        yh, yl = -result.orbitals[result.homo].x, -result.orbitals[result.lumo].x
        ax.annotate(
            "", xy=(-1.4, yl), xytext=(-1.4, yh),
            arrowprops=dict(arrowstyle="<->", color="#9ad08a", lw=1.4),
        )
        ax.text(-1.55, (yh + yl) / 2, f"gap\n{result.gap_beta:.3f}|β|",
                color="#9ad08a", ha="right", va="center", fontsize=9)

    ax.set_xlim(-3.0, 3.6)
    ax.set_ylim(-2.6, 2.6)
    # Numeric y-ticks would read as -x (misleading), so we show a direction arrow
    # instead and let each level's "α + x β" label carry the quantitative value.
    ax.set_ylabel("Energy", color=_FG)
    ax.annotate("", xy=(-2.85, 2.3), xytext=(-2.85, -2.3),
                arrowprops=dict(arrowstyle="-|>", color="#7a8aa8", lw=1.2))
    ax.text(-2.75, 2.15, "less stable", color="#7a8aa8", fontsize=8, rotation=90, va="top")
    ax.text(-2.75, -2.15, "more stable", color="#7a8aa8", fontsize=8, rotation=90, va="bottom")
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)

    s = mo_summary(result, is_monocycle)
    sub = f"{result.n_pi_electrons} π e⁻"
    if s.delocalization_energy_beta is not None:
        sub += f"   •   E_deloc = {s.delocalization_energy_beta:.3f} β"
    if s.aromaticity:
        sub += f"   •   {s.aromaticity}"
    ax.set_title(f"{title}\n{sub}" if title else sub, color=_FG, fontsize=13)

    fig.tight_layout()
    return fig
