"""Dash UI: atomic orbitals, Huckel MO sharing view, and MO diagrams.

    python -m orbital_explorer.app
"""

from __future__ import annotations

import base64
import io

import matplotlib
import plotly.graph_objects as go

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from dash import Dash, Input, Output, State, dcc, html

from .atomic import Orbital, available_orbitals
from .diatomic import solve_diatomic
from .diatomic_diagram import diatomic_diagram
from .huckel import solve
from .molecule import is_monocyclic_pi, resolve
from .mo_diagram import mo_diagram, mo_summary
from .render_atomic import orbital_figure
from .render_mo import mo_figure

_BG = "#0b0e17"
_PANEL = "#121726"
_FG = "#e8edf6"
_ACCENT = "#5aa0ff"

_LABEL = {"color": _FG, "fontSize": "0.85rem", "marginRight": "0.4rem"}
_CTRL = {"backgroundColor": _PANEL, "color": _FG, "border": "1px solid #2a3447",
         "borderRadius": "6px", "padding": "0.3rem"}


def _empty(message: str) -> go.Figure:
    fig = go.Figure()
    fig.add_annotation(text=message, showarrow=False, font=dict(color=_FG, size=16))
    fig.update_layout(paper_bgcolor=_BG, plot_bgcolor=_BG,
                      xaxis=dict(visible=False), yaxis=dict(visible=False),
                      margin=dict(l=0, r=0, t=0, b=0))
    return fig


def _fig_to_img(fig) -> str:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", facecolor=fig.get_facecolor(), dpi=110)
    plt.close(fig)
    buf.seek(0)
    return "data:image/png;base64," + base64.b64encode(buf.read()).decode()


# --------------------------------------------------------------------------- #
# Layout
# --------------------------------------------------------------------------- #
def _atomic_tab() -> dcc.Tab:
    return dcc.Tab(
        label="Atomic orbitals", value="atomic",
        style={"backgroundColor": _PANEL, "color": _FG, "border": "none"},
        selected_style={"backgroundColor": _BG, "color": _ACCENT, "borderTop": f"2px solid {_ACCENT}"},
        children=[html.Div(style={"padding": "1rem"}, children=[
            html.Div(style={"display": "flex", "gap": "1.2rem", "alignItems": "center",
                            "flexWrap": "wrap", "marginBottom": "0.6rem"}, children=[
                html.Div([html.Span("n", style=_LABEL),
                          dcc.Dropdown(id="orb-n", options=[{"label": i, "value": i} for i in range(1, 7)],
                                       value=2, clearable=False, style={"width": "80px", **_CTRL})]),
                html.Div([html.Span("orbital", style=_LABEL),
                          dcc.Dropdown(id="orb-label", value="pz", clearable=False,
                                       style={"width": "150px", **_CTRL})]),
                html.Div([html.Span("Z (nuclear charge)", style=_LABEL),
                          dcc.Input(id="orb-z", type="number", value=1, min=1, max=20, step=1,
                                    style={"width": "70px", **_CTRL})]),
                html.Div([html.Span("show", style=_LABEL),
                          dcc.RadioItems(id="orb-mode",
                                         options=[{"label": " phase (±)", "value": "phase"},
                                                  {"label": " density |ψ|²", "value": "density"}],
                                         value="phase", inline=True,
                                         labelStyle={"color": _FG, "marginRight": "0.9rem"},
                                         style={"display": "inline-block"})]),
            ]),
            dcc.Loading(dcc.Graph(id="atomic-graph", style={"height": "70vh"}), color=_ACCENT),
        ])],
    )


def _mo_tab() -> dcc.Tab:
    return dcc.Tab(
        label="Molecular orbitals", value="mo",
        style={"backgroundColor": _PANEL, "color": _FG, "border": "none"},
        selected_style={"backgroundColor": _BG, "color": _ACCENT, "borderTop": f"2px solid {_ACCENT}"},
        children=[html.Div(style={"padding": "1rem"}, children=[
            html.Div(children=[
                html.Span("molecular orbital", style=_LABEL),
                dcc.Dropdown(id="mo-select", style={"width": "360px", "display": "inline-block",
                                                    "verticalAlign": "middle", **_CTRL}),
            ]),
            dcc.Loading(dcc.Graph(id="mo-graph", style={"height": "64vh"}), color=_ACCENT),
            html.P("Blue and gold lobes are opposite phases. Where neighbouring lobes match, "
                   "the orbital is bonding; a node between atoms means antibonding.",
                   style={"color": "#9aa6bd", "fontSize": "0.8rem"}),
        ])],
    )


def _diagram_tab() -> dcc.Tab:
    return dcc.Tab(
        label="MO diagram & info", value="diagram",
        style={"backgroundColor": _PANEL, "color": _FG, "border": "none"},
        selected_style={"backgroundColor": _BG, "color": _ACCENT, "borderTop": f"2px solid {_ACCENT}"},
        children=[html.Div(style={"padding": "1rem", "display": "flex", "gap": "1.5rem",
                                  "flexWrap": "wrap"}, children=[
            html.Img(id="diagram-img", style={"maxHeight": "78vh", "borderRadius": "8px"}),
            html.Div(id="diagram-info", style={"color": _FG, "minWidth": "260px"}),
        ])],
    )


_DARK_CSS = """
<style>
  .Select-control, .Select-menu-outer, .Select-menu, .Select-option,
  .VirtualizedSelectOption { background-color:#121726 !important; color:#e8edf6 !important; }
  .Select-value-label, .Select-placeholder, .Select-input > input,
  .Select--single > .Select-control .Select-value { color:#e8edf6 !important; }
  .Select-option.is-focused, .VirtualizedSelectFocusedOption { background-color:#1c2742 !important; }
  .Select-arrow { border-color:#9aa6bd transparent transparent; }
</style>
"""


def build_app() -> Dash:
    app = Dash(__name__, title="Orbital Explorer", suppress_callback_exceptions=True)
    app.index_string = (
        "<!DOCTYPE html><html><head>{%metas%}<title>{%title%}</title>"
        "{%favicon%}{%css%}" + _DARK_CSS + "</head><body>{%app_entry%}"
        "<footer>{%config%}{%scripts%}{%renderer%}</footer></body></html>"
    )
    app.layout = html.Div(style={"backgroundColor": _BG, "minHeight": "100vh",
                                 "fontFamily": "system-ui, sans-serif"}, children=[
        html.Div(style={"padding": "0.8rem 1rem"}, children=[
            html.H2("Orbital Explorer", style={"color": _FG, "margin": 0}),
            html.Span("atomic orbitals, Huckel pi MOs, and diatomic diagrams for MO revision",
                      style={"color": "#9aa6bd", "fontSize": "0.85rem"}),
            html.Div(style={"display": "flex", "gap": "0.6rem", "alignItems": "center",
                            "flexWrap": "wrap", "marginTop": "0.6rem"}, children=[
                html.Span("molecule", style=_LABEL),
                dcc.Input(id="mol-input", type="text", value="benzene", debounce=True,
                          autoComplete="off",
                          placeholder="name or SMILES, e.g. benzene or C=CC=C",
                          style={"width": "320px", **_CTRL}),
                html.Button("Solve", id="mol-submit", n_clicks=0,
                            style={"backgroundColor": _ACCENT, "color": "#04101f", "border": "none",
                                   "borderRadius": "6px", "padding": "0.4rem 1rem", "cursor": "pointer"}),
                html.Span(id="mol-note", style={"color": "#9aa6bd", "marginLeft": "0.6rem",
                                                "fontSize": "0.85rem"}),
            ]),
        ]),
        dcc.Tabs(id="tabs", value="atomic", children=[_atomic_tab(), _mo_tab(), _diagram_tab()]),
        dcc.Store(id="mol-store"),
    ])
    _register(app)
    return app


# --------------------------------------------------------------------------- #
# Callbacks
# --------------------------------------------------------------------------- #
def _register(app: Dash) -> None:

    @app.callback(Output("orb-label", "options"), Output("orb-label", "value"),
                  Input("orb-n", "value"), State("orb-label", "value"))
    def _labels(n, current):
        opts = [{"label": lab, "value": lab} for lab in available_orbitals(n)]
        value = current if any(o["value"] == current for o in opts) else opts[-1]["value"]
        return opts, value

    @app.callback(Output("atomic-graph", "figure"),
                  Input("orb-n", "value"), Input("orb-label", "value"),
                  Input("orb-z", "value"), Input("orb-mode", "value"))
    def _atomic(n, label, z, mode):
        try:
            orb = Orbital(n=int(n), label=label, Z=float(z or 1))
        except (ValueError, TypeError):
            return _empty("pick a valid orbital (l must be < n)")
        return orbital_figure(orb, mode=mode)

    @app.callback(Output("mol-store", "data"), Output("mol-note", "children"),
                  Input("mol-submit", "n_clicks"), Input("mol-input", "value"))
    def _solve(_clicks, text):
        info = resolve(text or "")
        note = f"{info.formula or '?'} · {info.category} · {info.note}"
        return {"smiles": info.smiles, "category": info.category}, note

    @app.callback(Output("mo-select", "options"), Output("mo-select", "value"),
                  Input("mol-store", "data"))
    def _selector(data):
        if not data or data.get("category") != "pi-system":
            return [], None
        info = resolve(data["smiles"])
        res = solve(info.mol, info.pi_atom_indices)
        options = [
            {"label": f"ψ{o.index + 1}: {o.energy_label}"
                      f"  [{'●●' if o.occupation == 2 else '●' if o.occupation == 1 else '○'}]",
             "value": o.index}
            for o in res.orbitals
        ]
        return options, (res.homo if res.homo is not None else 0)

    @app.callback(Output("mo-graph", "figure"),
                  Input("mol-store", "data"), Input("mo-select", "value"))
    def _mo(data, mo_index):
        if data and data.get("category") == "diatomic":
            return _empty("Diatomics: see the MO diagram tab (the 3D σ/π view isn't drawn here yet).")
        if not data or data.get("category") != "pi-system" or mo_index is None:
            return _empty("type a conjugated molecule (benzene, butadiene, pyridine …) and press Solve")
        try:
            info = resolve(data["smiles"])
            res = solve(info.mol, info.pi_atom_indices)
            return mo_figure(info.mol, res, int(mo_index))
        except Exception as exc:  # never leave the graph silently blank
            return _empty(f"Couldn't render this orbital: {exc}")

    @app.callback(Output("diagram-img", "src"), Output("diagram-info", "children"),
                  Input("mol-store", "data"))
    def _diagram(data):
        if not data:
            return None, _hint("Enter a molecule above to see its MO diagram.")
        category = data.get("category")
        try:
            if category == "pi-system":
                return _pi_diagram_view(resolve(data["smiles"]))
            if category == "diatomic":
                return _diatomic_diagram_view(resolve(data["smiles"]))
        except Exception as exc:  # show the problem instead of a blank panel
            return None, _hint(f"Couldn't build the diagram: {exc}", warn=True)
        return None, _hint("No MO diagram for this one. Try a conjugated molecule "
                           "(benzene, butadiene…) or a diatomic (O₂, N₂, CO).")


def _hint(message: str, warn: bool = False):
    return html.Span(message, style={"color": "#e0a85a" if warn else "#9aa6bd"})


def _table(rows):
    return html.Table([html.Tbody([
        html.Tr([html.Td(k, style={"padding": "0.25rem 0.8rem 0.25rem 0", "color": "#9aa6bd"}),
                 html.Td(str(v), style={"padding": "0.25rem 0", "fontWeight": "600"})])
        for k, v in rows
    ])], style={"borderCollapse": "collapse"})


def _pi_diagram_view(info):
    res = solve(info.mol, info.pi_atom_indices)
    cyclic = is_monocyclic_pi(info)
    fig = mo_diagram(res, title=info.formula or "", is_monocycle=cyclic)
    s = mo_summary(res, cyclic)
    homo_lbl = res.orbitals[res.homo].energy_label if res.homo is not None else "n/a"
    lumo_lbl = res.orbitals[res.lumo].energy_label if res.lumo is not None else "n/a"
    total = s.total_pi_energy_beta
    rows = [
        ("π electrons", s.n_pi_electrons),
        ("HOMO", homo_lbl),
        ("LUMO", lumo_lbl),
        ("HOMO-LUMO gap", f"{s.gap_beta:.3f} |β|" if s.gap_beta is not None else "n/a"),
        ("total π energy", f"{res.n_pi_electrons}α {'+' if total >= 0 else '−'} {abs(total):.3f} β"),
        ("delocalisation E", f"{s.delocalization_energy_beta:.3f} β"
                             if s.delocalization_energy_beta is not None else "open shell"),
        ("aromaticity", s.aromaticity or "n/a"),
    ]
    return _fig_to_img(fig), [html.H4("Key numbers", style={"marginTop": 0}), _table(rows),
                              _hint(info.note)]


def _diatomic_diagram_view(info):
    res = solve_diatomic(info.mol, info.charge)
    if not res.supported:
        return None, _hint(res.note)
    fig = diatomic_diagram(res, title=info.formula or "")
    homo = res.orbitals[res.homo].name if res.homo is not None else "n/a"
    lumo = res.orbitals[res.lumo].name if res.lumo is not None else "n/a"
    rows = [
        ("valence electrons", res.n_valence_electrons),
        ("bond order", f"{res.bond_order:g}"),
        ("HOMO", homo),
        ("LUMO", lumo),
        ("unpaired e⁻", res.n_unpaired),
        ("magnetism", res.magnetism),
    ]
    if any(o.name == "σ2p" for o in res.orbitals):
        rows.append(("s-p mixing", "on (π2p below σ2p)" if res.mixing else "off (σ2p below π2p)"))
    return _fig_to_img(fig), [html.H4("Key numbers", style={"marginTop": 0}), _table(rows)]


def main() -> None:
    """Run the Dash app. Uses PORT and HOST from the environment if set."""
    import os

    port = int(os.environ.get("PORT", "8050"))
    # Local default stays loopback; production (Docker/Render) should set HOST=0.0.0.0
    host = os.environ.get("HOST", "127.0.0.1")
    debug = os.environ.get("DEBUG", "").lower() in {"1", "true", "yes"}
    build_app().run(debug=debug, host=host, port=port)


if __name__ == "__main__":
    main()
