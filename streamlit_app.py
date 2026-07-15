"""Streamlit front-end for free hosting (Streamlit Community Cloud).

Deploy: https://share.streamlit.io → New app → this repo → main → streamlit_app.py

Local::

    pip install -e ".[cloud]"
    streamlit run streamlit_app.py
"""

from __future__ import annotations

import io
import sys
from pathlib import Path

# Editable / source-tree import without install
_SRC = Path(__file__).resolve().parent / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import streamlit as st

from orbital_explorer.atomic import Orbital, available_orbitals
from orbital_explorer.diatomic import solve_diatomic
from orbital_explorer.diatomic_diagram import diatomic_diagram
from orbital_explorer.huckel import solve
from orbital_explorer.mo_diagram import mo_diagram, mo_summary
from orbital_explorer.molecule import is_monocyclic_pi, resolve
from orbital_explorer.render_atomic import orbital_figure
from orbital_explorer.render_mo import mo_figure

st.set_page_config(
    page_title="Orbital Explorer",
    page_icon="⚛️",
    layout="wide",
    initial_sidebar_state="expanded",
    menu_items={
        "Get help": "https://github.com/nikshaybisht/orbital-explorer#install",
        "Report a bug": "https://github.com/nikshaybisht/orbital-explorer/issues",
        "About": (
            "Orbital Explorer v0.1.1 — teaching tool for undergrad MO theory. "
            "Hydrogenic AOs, Huckel pi MOs, diatomic diagrams. MIT License."
        ),
    },
)

_CSS = """
<style>
  .block-container { padding-top: 1.2rem; }
  h1 { font-weight: 700; letter-spacing: -0.02em; }
  .oe-muted { color: #9aa6bd; font-size: 0.95rem; }
</style>
"""
st.markdown(_CSS, unsafe_allow_html=True)


def _fig_png(fig) -> bytes:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", facecolor=fig.get_facecolor(), dpi=120, bbox_inches="tight")
    plt.close(fig)
    buf.seek(0)
    return buf.getvalue()


@st.cache_data(show_spinner=False)
def _atomic_fig(n: int, label: str, mode: str, z: float):
    orb = Orbital(n=n, label=label, Z=z)
    return orbital_figure(orb, mode=mode, n_grid=72)


@st.cache_data(show_spinner="Solving molecule…")
def _resolve_mol(text: str):
    return resolve(text)


def main() -> None:
    st.title("Orbital Explorer")
    st.markdown(
        '<p class="oe-muted">Atomic and molecular orbitals for undergrad MO revision: '
        "hydrogenic s/p/d/f, Huckel pi MOs, and diatomic diagrams. "
        "Teaching model, not a full quantum chemistry package.</p>",
        unsafe_allow_html=True,
    )

    tab_atomic, tab_mo, tab_diagram = st.tabs(
        ["Atomic orbitals", "Molecular orbitals (sharing)", "MO diagram & info"]
    )

    with tab_atomic:
        c1, c2, c3, c4 = st.columns(4)
        n = c1.selectbox("n", list(range(1, 7)), index=1)
        labels = available_orbitals(n)
        label = c2.selectbox("orbital", labels, index=min(2, len(labels) - 1) if labels else 0)
        z = c3.number_input("Z (nuclear charge)", min_value=1.0, max_value=20.0, value=1.0, step=0.5)
        mode = c4.selectbox("view", ["phase", "density"], index=0)
        try:
            fig = _atomic_fig(int(n), label, mode, float(z))
            st.plotly_chart(fig, use_container_width=True)
            st.caption(f"Exact hydrogenic ψ for **{Orbital(n=int(n), label=label, Z=float(z)).name}**.")
        except Exception as exc:  # show error in UI
            st.error(f"Could not render orbital: {exc}")

    with tab_mo:
        mol_text = st.text_input(
            "Molecule (name or SMILES)",
            value="benzene",
            help="Examples: benzene, pyridine, C=CC=C, O2, allyl cation",
        )
        if st.button("Solve", type="primary", key="solve_mo") or mol_text:
            try:
                info = _resolve_mol(mol_text.strip())
                st.info(info.note or info.formula or mol_text)
                if info.category == "diatomic" or (
                    info.mol is not None
                    and info.n_heavy_atoms == 2
                    and not info.pi_atom_indices
                ):
                    st.warning(
                        "This looks like a **diatomic**. Open the **MO diagram** tab "
                        "for the three-column diagram (bond order, magnetism)."
                    )
                elif not info.pi_atom_indices or info.mol is None:
                    st.warning("No conjugated π system detected for 3D orbital-sharing view.")
                else:
                    res = solve(info.mol, info.pi_atom_indices)
                    n_mo = len(res.orbitals)
                    default_mo = res.homo if res.homo is not None else 0
                    mo_idx = st.slider(
                        "MO index (0 = lowest energy)",
                        0,
                        max(0, n_mo - 1),
                        value=int(default_mo),
                    )
                    fig = mo_figure(info.mol, res, mo_idx, n_grid=48)
                    st.plotly_chart(fig, use_container_width=True)
                    orb = res.orbitals[mo_idx]
                    tag = "HOMO" if mo_idx == res.homo else ("LUMO" if mo_idx == res.lumo else "")
                    st.caption(
                        f"MO {mo_idx} {tag} · E = {orb.energy_label} · "
                        "Blue/gold = LCAO phase; merged lobes = bonding, nodes = antibonding."
                    )
            except Exception as exc:  # noqa: BLE001
                st.error(f"Could not solve molecule: {exc}")

    with tab_diagram:
        mol_text = st.text_input(
            "Molecule (name or SMILES)",
            value="O2",
            key="diagram_mol",
            help="Try O2, N2, CO, benzene, butadiene…",
        )
        if st.button("Build diagram", type="primary", key="solve_diagram") or mol_text:
            try:
                info = _resolve_mol(mol_text.strip())
                if info.mol is None:
                    st.error(info.note or "Could not resolve that molecule.")
                    return
                # Prefer diatomic engine when applicable
                dia = solve_diatomic(info.mol, info.charge)
                if dia.supported:
                    fig = diatomic_diagram(dia, title=info.formula or mol_text)
                    st.image(_fig_png(fig), use_container_width=True)
                    c1, c2, c3, c4 = st.columns(4)
                    c1.metric("Bond order", f"{dia.bond_order:g}")
                    c2.metric("Unpaired e⁻", dia.n_unpaired)
                    c3.metric("Magnetism", dia.magnetism)
                    homo = dia.orbitals[dia.homo].name if dia.homo is not None else "n/a"
                    lumo = dia.orbitals[dia.lumo].name if dia.lumo is not None else "n/a"
                    c4.metric("HOMO / LUMO", f"{homo} / {lumo}")
                    if dia.note:
                        st.caption(dia.note)
                elif info.pi_atom_indices:
                    res = solve(info.mol, info.pi_atom_indices)
                    cyclic = is_monocyclic_pi(info)
                    fig = mo_diagram(res, title=info.formula or "", is_monocycle=cyclic)
                    st.image(_fig_png(fig), use_container_width=True)
                    s = mo_summary(res, cyclic)
                    c1, c2, c3 = st.columns(3)
                    c1.metric("π electrons", s.n_pi_electrons)
                    c2.metric(
                        "HOMO-LUMO gap",
                        f"{s.gap_beta:.3f} |β|" if s.gap_beta is not None else "n/a",
                    )
                    c3.metric("Aromaticity", s.aromaticity or "n/a")
                    if info.note:
                        st.caption(info.note)
                else:
                    st.warning(dia.note or "No supported MO diagram for this input.")
            except Exception as exc:  # noqa: BLE001
                st.error(f"Could not build diagram: {exc}")

    st.divider()
    st.caption(
        "Source: [github.com/nikshaybisht/orbital-explorer](https://github.com/nikshaybisht/orbital-explorer) · "
        "MIT License · v0.1.1 · "
        "[Static demo](https://nikshaybisht.github.io/orbital-explorer/)"
    )


if __name__ == "__main__":
    main()
