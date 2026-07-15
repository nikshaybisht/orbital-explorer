# Launch notes

Public demo: https://nikshaybisht.github.io/orbital-explorer/

| Day | Task | Status |
|-----|------|--------|
| 1 | MIT LICENSE | done |
| 1 | pytest + version pins | done |
| 2 | Tag v0.1.0 | done |
| 3-4 | Public URL | done (GitHub Pages) |
| 5 | Demo GIF (benzene HOMO + O2) | done (`docs/images/demo.gif`) |
| 6 | Post | drafts below (post on X yourself) |
| 7 | CV + freeze | done (`docs/CV.md`, issue #1) |

---

## Extra deploy options

### Streamlit Cloud

1. https://share.streamlit.io with GitHub login
2. New app → `nikshaybisht/orbital-explorer` → `main` → `streamlit_app.py`
3. Put the `*.streamlit.app` URL in the README if you want a full hosted UI

### Render

Blueprint from `render.yaml` on https://dashboard.render.com

### Docker

```bash
docker build -t orbital-explorer .
docker run --rm -p 8050:8050 -e PORT=8050 orbital-explorer
```

---

## Screen recording (~75 s), if you want video later

| Time | Do this | Rough line |
|------|---------|------------|
| 0-10s | Title / demo site | "Orbital Explorer. Atomic and molecular orbitals in the browser." |
| 10-25s | Atomic: n=2 pz, then n=3 dz2 | "Exact hydrogenic wavefunctions. Phase lobes, not just a density blob." |
| 25-50s | MO tab, benzene, HOMO | "Huckel: p orbitals sharing on the ring. Same phase bonds, opposite phase gives a node." |
| 50-70s | Diagram tab, O2 | "Standard diatomic diagram. Bond order 2, two unpaired electrons, paramagnetic." |
| 70-80s | Repo + demo link | "MIT, free. Links in the description." |

---

## X / LinkedIn draft (copy-paste)

```
I built Orbital Explorer while revising undergrad MO theory.

- exact hydrogenic s/p/d/f in 3D
- Huckel pi MOs (try benzene HOMO)
- diatomic diagrams with O2 paramagnetism done properly

MIT, open source
https://github.com/nikshaybisht/orbital-explorer
demo: https://nikshaybisht.github.io/orbital-explorer/

#chemistry #compchem #python
```

Shorter X version:

```
Orbital Explorer: free undergrad MO visualiser.

Hydrogenic AOs, Huckel pi MOs (benzene HOMO), diatomic ladders (O2 paramagnetic).

https://github.com/nikshaybisht/orbital-explorer
https://nikshaybisht.github.io/orbital-explorer/
```

### Reddit (r/chemistry)

**Title:** Free undergrad MO visualiser (Huckel + diatomics + hydrogenic orbitals)

**Body:**

I put together a small teaching app for MO concepts. Type a molecule (name or SMILES) and you get a Huckel energy ladder plus a 3D "orbital sharing" view, or a classic diatomic diagram (O2 comes out paramagnetic with two unpaired electrons). Atomic tab is exact hydrogenic psi.

Not a replacement for ORCA/Gaussian. It is deliberately textbook-level so it installs with pip and runs in a browser.

Repo: https://github.com/nikshaybisht/orbital-explorer  
Demo: https://nikshaybisht.github.io/orbital-explorer/  
License: MIT. Happy to hear about chemistry edge cases.

---

## Freeze (until 2026-07-29)

Tracked in issue #1.

- OK: bugfixes, docs, dep pins, deploy
- Not OK: new features
- Hotfixes as `v0.1.x`
