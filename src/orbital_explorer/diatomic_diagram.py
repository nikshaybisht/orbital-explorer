"""The classic three-column molecular-orbital diagram for a diatomic.

Atomic orbitals of the two atoms sit in the left and right columns (the more
electronegative atom's levels lower); the molecular orbitals they form sit in the
centre, filled with electron arrows. Bonding MOs are blue, antibonding red,
nonbonding/lone-pair olive -- the same palette as the Huckel diagram.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from .diatomic import DiatomicResult
from .mo_diagram import _ANTI, _BG, _BOND, _FG, _NONBOND, _draw_electrons

_MO_X = 0.0
_AO_X = 2.7
_SEG = 0.16
_DEGEN_DX = 0.33


def _mo_color(character: str) -> str:
    return {"bond": _BOND, "anti": _ANTI}.get(character, _NONBOND)


def _mo_shell(name: str) -> str | None:
    if any(t in name for t in ("1s", "2s", "3s")) or name == "ns(lp)":
        return "s"
    if "π(lp)" in name or any(t in name for t in ("2p", "3p")):
        return "p"
    return None  # H-X sigma forms from both s and p -> tie to everything


def diatomic_diagram(result: DiatomicResult, title: str = ""):
    """Build a matplotlib Figure of the diatomic MO diagram."""
    fig, ax = plt.subplots(figsize=(7.2, 7.6))
    fig.patch.set_facecolor(_BG)
    ax.set_facecolor(_BG)

    n = len(result.orbitals)
    y_top = n - 1

    # --- atomic-orbital columns, positioned by energy (lower energy = lower) ---
    all_e = [e for levels in result.ao_levels.values() for _, e in levels]
    emin, emax = min(all_e), max(all_e)
    span = (emax - emin) or 1.0

    def ao_y(energy: float) -> float:
        return 0.3 + (energy - emin) / span * (y_top - 0.6)

    sym_a, sym_b = result.symbols
    columns = [(-_AO_X, sym_a, result.ao_levels[sym_a])]
    levels_b = result.ao_levels.get(sym_b, result.ao_levels[sym_a])  # homonuclear reuses A
    columns.append((_AO_X, sym_b, levels_b))

    ao_points: list[tuple[float, float, str]] = []  # (x_inner, y, shell)
    for x, sym, levels in columns:
        ax.text(x, y_top + 1.0, sym, color=_FG, ha="center", fontsize=15, weight="bold")
        for sub, energy in levels:
            y = ao_y(energy)
            ax.plot([x - 0.3, x + 0.3], [y, y], color="#8595b5", lw=2.5, solid_capstyle="round")
            label_x = x - 0.5 if x < 0 else x + 0.5
            ax.text(label_x, y, sub, color="#9aa6bd", va="center",
                    ha="right" if x < 0 else "left", fontsize=9)
            ao_points.append((x + (0.3 if x < 0 else -0.3), y, sub[-1]))

    # --- tie-lines: AO -> the MOs of the same shell (faint guides) ---
    for mo in result.orbitals:
        shell = _mo_shell(mo.name)
        for x_inner, y_ao, ao_shell in ao_points:
            if shell is None or shell == ao_shell:
                edge = -_SEG if x_inner < 0 else _SEG
                ax.plot([x_inner, _MO_X + edge], [y_ao, mo.index],
                        color="#39507a", lw=0.6, ls=":", alpha=0.5, zorder=0)

    # --- central MO column with electrons ---
    for mo in result.orbitals:
        y = mo.index
        color = _mo_color(mo.character)
        centers = [-_DEGEN_DX, _DEGEN_DX] if mo.degeneracy == 2 else [0.0]
        for xc, occ in zip(centers, mo.occupations):
            ax.plot([xc - _SEG, xc + _SEG], [y, y], color=color, lw=3.5, solid_capstyle="round")
            _draw_electrons(ax, xc, y, occ)
        ax.text(0.62, y, mo.name, color=_FG, va="center", fontsize=9)
        if mo.index == result.homo:
            ax.text(-0.62, y, "HOMO", color=_BOND, ha="right", va="center", fontsize=8, weight="bold")
        if mo.index == result.lumo:
            ax.text(-0.62, y, "LUMO", color=_ANTI, ha="right", va="center", fontsize=8, weight="bold")

    # --- frame & caption ---
    ax.set_xlim(-3.7, 3.7)
    ax.set_ylim(-0.9, y_top + 1.6)
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    ax.annotate("", xy=(-3.5, y_top + 0.5), xytext=(-3.5, -0.5),
                arrowprops=dict(arrowstyle="-|>", color="#7a8aa8", lw=1.2))
    ax.text(-3.62, y_top / 2, "Energy", color="#7a8aa8", rotation=90, va="center", fontsize=10)

    bo = result.bond_order
    sub = f"bond order {bo:g}   •   {result.magnetism}"
    if result.n_unpaired:
        sub += f" ({result.n_unpaired} unpaired e⁻)"
    if any(mo.name == "σ2p" for mo in result.orbitals):
        sub += "\n" + ("s-p mixing on (pi2p below sigma2p)" if result.mixing else "no s-p mixing (sigma2p below pi2p)")
    ax.set_title(f"{title}\n{sub}" if title else sub, color=_FG, fontsize=13)

    fig.tight_layout()
    return fig
