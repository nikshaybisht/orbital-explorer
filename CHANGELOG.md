# Changelog

Format loosely follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).
Versioning follows [SemVer](https://semver.org/spec/v2.0.0.html).

## [0.1.1] (2026-07-16)

### Demo / docs

- Expanded GitHub Pages gallery: multiple interactive atomic orbitals (2pz, 2px, 3dz2, 3dxy, 4fxyz)
- Interactive benzene HOMO HTML when RDKit is available in the build
- Dark chrome headers on Plotly demo pages with back-link to home
- Open Graph / Twitter meta tags on the demo landing page
- Streamlit Cloud hints: `packages.txt`, `runtime.txt` (Python 3.12)
- Issue templates for bugs and feature ideas

### Fixed / polished

- Pages deploy copies all `docs/demo/*.html` files
- Version bump to 0.1.1 across package metadata and Streamlit footer

## [0.1.0] (2026-07-15)

First public freeze.

### Packaging / project

- MIT license
- `pyproject.toml` with console script `orbital-explorer`
- Version ranges in requirements; pytest and streamlit included
- `wsgi.py`, `Dockerfile`, `render.yaml`
- `streamlit_app.py` for Streamlit Community Cloud
- CI: pytest on Ubuntu and Windows (3.11, 3.12), Docker build
- `HOST` / `PORT` env support
- Demo site + tour GIF on GitHub Pages
- CONTRIBUTING, CV notes, launch notes

### Science already in the tree

- Exact hydrogenic atomic orbitals (s/p/d/f) as 3D isosurfaces
- Huckel pi-MO solver with 3D LCAO view and energy diagrams
- Qualitative diatomic MO diagrams (H-Ne, H-X) with s-p mixing and correct O2 paramagnetism
- Name / SMILES resolution (including PubChem when online)

[0.1.1]: https://github.com/nikshaybisht/orbital-explorer/releases/tag/v0.1.1
[0.1.0]: https://github.com/nikshaybisht/orbital-explorer/releases/tag/v0.1.0
