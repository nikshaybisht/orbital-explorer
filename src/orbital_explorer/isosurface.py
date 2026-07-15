"""Smooth isosurface meshes from a scalar field.

Instead of handing the whole voxel grid to Plotly (blocky, and tens of MB per
figure), we extract just the surface with marching cubes and render it as a
triangle mesh with *interpolated* normals (flatshading=False) plus a lighting
model. The result is a smooth, glossy lobe and a far smaller figure.

Used by both the atomic-orbital and molecular-orbital views so they look and
behave the same.
"""

from __future__ import annotations

import numpy as np
import plotly.graph_objects as go
from skimage import measure

# Phase colours: positive lobe blue, negative lobe gold (Orbitron convention).
POS_COLOR = "#1f77ff"
NEG_COLOR = "#f4c430"
DENSITY_COLOR = "#2ec4b6"

# Flat, low-specular lighting (close to Plotly's default) so translucent surfaces
# read as airy and see-through rather than as a glossy solid. Smoothness comes
# from the mesh geometry, not from shiny highlights.
_LIGHTING = dict(ambient=0.82, diffuse=0.6, specular=0.05, roughness=0.6, fresnel=0.15)
_LIGHTPOS = dict(x=120, y=200, z=180)


def _mesh_at_level(values, level, spacing, origin, color, opacity):
    """One Mesh3d for a single isolevel, or None if the level isn't crossed."""
    try:
        verts, faces, _normals, _ = measure.marching_cubes(values, level=level, spacing=spacing)
    except (ValueError, RuntimeError):
        return None  # level outside the data range (e.g. 1s has no negative lobe)
    verts = verts + origin
    return go.Mesh3d(
        x=verts[:, 0], y=verts[:, 1], z=verts[:, 2],
        i=faces[:, 0], j=faces[:, 1], k=faces[:, 2],
        color=color, opacity=opacity, flatshading=False,
        lighting=_LIGHTING, lightposition=_LIGHTPOS,
        hoverinfo="skip", showscale=False,
    )


def _geometry(axes):
    xs, ys, zs = axes
    spacing = (float(xs[1] - xs[0]), float(ys[1] - ys[0]), float(zs[1] - zs[0]))
    origin = np.array([xs[0], ys[0], zs[0]])
    return spacing, origin


def phase_meshes(values: np.ndarray, axes, iso: float, opacity: float = 0.32) -> list[go.Mesh3d]:
    """Two smooth lobes: the +iso surface (blue) and the -iso surface (gold).

    Rendered translucent (opacity < 1) so the far lobe and any overlap show
    through -- the see-through look that reads better than solid surfaces.
    """
    spacing, origin = _geometry(axes)
    traces = []
    for level, color in ((iso, POS_COLOR), (-iso, NEG_COLOR)):
        mesh = _mesh_at_level(values, level, spacing, origin, color, opacity)
        if mesh is not None:
            traces.append(mesh)
    return traces


def density_meshes(values: np.ndarray, axes, levels) -> list[go.Mesh3d]:
    """Nested translucent shells of |psi|^2. `levels` is a list of (iso, opacity)."""
    spacing, origin = _geometry(axes)
    traces = []
    for iso, opacity in levels:
        mesh = _mesh_at_level(values, iso, spacing, origin, DENSITY_COLOR, opacity)
        if mesh is not None:
            traces.append(mesh)
    return traces
