"""Interactive 3D rendering of atomic orbitals with Plotly.

Two viewing modes:

* "phase"   -- the signed wavefunction psi. Positive and negative lobes get
               different colours (blue / gold), exactly like the Orbitron
               gallery. This is the view that makes bonding intuition click,
               because the lobe signs are what add or cancel when atoms bond.
* "density" -- the probability cloud |psi|^2, a single-colour isosurface
               answering "where is the electron likely to be?".

Rendering uses marching-cubes isosurfaces (plotly.graph_objects.Isosurface),
which run entirely in the browser via WebGL -- rotate, zoom and pan are free.
"""

from __future__ import annotations

import numpy as np
import plotly.graph_objects as go

from .atomic import Orbital, probability_density, suggested_extent, wavefunction
from .isosurface import density_meshes, phase_meshes


def orbital_figure(
    orbital: Orbital,
    mode: str = "phase",
    n_grid: int = 96,
    iso_fraction: float = 0.10,
    extent: float | None = None,
) -> go.Figure:
    """Build an interactive Plotly figure for a single atomic orbital.

    iso_fraction sets the isosurface level as a fraction of the peak value, so
    the same setting gives a comparable "size" of cloud across orbitals. The
    surface is extracted as a smooth triangle mesh (see isosurface.py).
    """
    if extent is None:
        extent = suggested_extent(orbital)
    axis = np.linspace(-extent, extent, n_grid)
    X, Y, Z = np.meshgrid(axis, axis, axis, indexing="ij")

    if mode == "density":
        values = probability_density(orbital, X, Y, Z)
        peak = float(values.max())
        traces = density_meshes(values, (axis, axis, axis),
                                levels=[(0.5 * peak, 0.4), (0.2 * peak, 0.2), (0.05 * peak, 0.1)])
    else:  # phase
        values = wavefunction(orbital, X, Y, Z)
        peak = float(np.abs(values).max())
        traces = phase_meshes(values, (axis, axis, axis), iso=iso_fraction * peak)

    fig = go.Figure(data=traces)
    _style(fig, orbital, extent)
    return fig


def _style(fig: go.Figure, orbital: Orbital, extent: float) -> None:
    title = (
        f"{orbital.name}  "
        f"(n={orbital.n}, l={orbital.l}, Z={orbital.Z:g})  "
        f"({orbital.radial_nodes} radial / {orbital.angular_nodes} angular nodes)"
    )
    axis_kw = dict(
        showbackground=True, backgroundcolor="#0b0e17",
        gridcolor="#222a3a", zerolinecolor="#39507a",
        range=[-extent, extent], title="",
    )
    fig.update_layout(
        title=dict(text=title, x=0.5, font=dict(color="#e8edf6", size=18)),
        paper_bgcolor="#0b0e17",
        scene=dict(
            xaxis=axis_kw, yaxis=axis_kw, zaxis=axis_kw,
            aspectmode="cube",
            annotations=[],
        ),
        margin=dict(l=0, r=0, t=48, b=0),
    )


def show(orbital: Orbital, **kwargs) -> go.Figure:
    """Render the orbital and open it in the default browser."""
    fig = orbital_figure(orbital, **kwargs)
    fig.show()
    return fig


def save_html(orbital: Orbital, path: str, **kwargs) -> str:
    """Render the orbital to a self-contained interactive HTML file."""
    fig = orbital_figure(orbital, **kwargs)
    fig.write_html(path, include_plotlyjs="cdn")
    return path
