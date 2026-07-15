# Orbital Explorer: Project Export

A complete write-up of the project: what it is, the science behind it, how the
code is organised, the design decisions and their rationale, what is correct vs.
approximate, how it is tested, and how it was built.

> **One line:** pip-installable MO visualiser. Hydrogenic AOs, Huckel pi systems,
> and diatomic diagrams, drawn in 3D. Aimed at undergrad revision, not production QC.

Repository: <https://github.com/nikshaybisht/orbital-explorer>

---

## 1. Motivation

The project was inspired by the look of the [qmsolve](https://github.com/quantum-visualizations/qmsolve)
gallery and the classic Orbitron atomic-orbital images, but with a different
goal: instead of a general Schrödinger-equation solver, build a **chemistry-input
driven** tool where you type an atom or molecule and immediately see

* the orbital shapes,
* how atomic orbitals combine ("orbital sharing") into molecular orbitals,
* the molecular-orbital (MO) energy diagram with HOMO/LUMO and the key numbers.

A hard constraint shaped every technical choice: it must **install and run on
Windows with plain `pip`, no conda, no WSL**. That ruled out the heavy
quantum-chemistry backends (PySCF needs WSL; Psi4 needs conda) and pushed the
design toward methods that are *exact or qualitatively correct* yet lightweight:
analytic hydrogenic wavefunctions, Hückel MO theory, and a curated diatomic MO
model.

---

## 2. What it does

The app (`python run_app.py` → <http://127.0.0.1:8050>) is a three-tab Dash
application.

### Tab 1: Atomic orbitals
Pick `n`, the orbital (`s, p, d, f`: the chemist labels `dz2`, `dxy`, `fxyz`…)
and the nuclear charge `Z`. See the exact hydrogenic orbital in 3D, as either the
**phase** view (blue/gold ± lobes) or the **probability density** `|ψ|²` cloud.

### Tab 2: Molecular orbitals ("orbital sharing")
Type a molecule (name or SMILES). It is auto-classified, solved with **Hückel
theory**, and any MO is shown in 3D as the p-orbitals combining across the
framework: matching lobes = bonding, a node between atoms = antibonding.

### Tab 3: MO diagram & info
* **Conjugated molecule:** the Hückel energy ladder with electron arrows,
  HOMO/LUMO, gap, π-electron count, total π-energy, delocalisation (resonance)
  energy and aromaticity (4n+2 / 4n).
* **Diatomic** (`O2`, `N2`, `CO`, `HF`, `He2`, `B2`…): the classic three-column
  diagram with bond order, HOMO/LUMO and magnetism.

---

## 3. The science and methods

### 3.1 Atomic orbitals: *exact*
The hydrogenic wavefunction is the analytic solution of the one-electron
Schrödinger equation:

```
ψ(r,θ,φ) = R_{n,l}(r) · Y_{l,m}(θ,φ)
```

* **Radial part** `R_{n,l}(r)` is built from the associated Laguerre polynomial
  (`scipy.special.genlaguerre`), normalised so `∫ R² r² dr = 1`.
* **Angular part** uses the **real** spherical harmonics chemists actually draw
  (p_x/p_y/p_z, d_xy, d_x²−y², …), written directly as normalised polynomials in
  the unit-vector components: so the lobed shapes and their signs are exact, not
  fitted. Lengths are in Bohr radii; `Z` is adjustable (H 1s vs He⁺ 1s).

This part is **exact**: tests confirm `∫|ψ|² dV = 1.000`, the radial node count
equals `n−l−1`, and angular nodes equal `l`.

### 3.2 Molecular orbitals (π): *qualitative Hückel*
For conjugated π-systems the tool uses **Hückel MO theory**: each p-orbital is a
basis function, the atom energy `α` sits on the diagonal and the neighbour
interaction `β` on bonds; diagonalising gives MO energies `E = α + x·β` (β < 0,
so the *largest* x is most stable) and the eigenvectors are the LCAO
coefficients: literally the "orbital sharing." Heteroatoms use standard
Streitwieser parameters. This is **exact for the model** and reproduces the
textbook results (ethylene, allyl, butadiene, benzene, cyclobutadiene,
cyclopentadienyl, [n]annulenes) but is *qualitative*: energies are in units of β,
not kJ/mol, and it covers conjugated π-systems only.

The 3D "orbital sharing" view places a p-orbital on each atom weighted by its
LCAO coefficient (size = |c|, colour = sign) and renders the isosurface of the
sum `ψ_MO = Σ c_i p_{z,i}`: so bonding regions merge and antibonding nodes
appear, with the molecular skeleton drawn underneath.

### 3.3 Diatomic MOs: *qualitative but correct ordering*
Hückel does not cover the canonical MO examples (H₂, N₂, O₂, CO, HF), so a
dedicated module builds the standard qualitative diatomic picture from a curated
table of valence-AO energies. It gets the teaching-critical details right:

* the **s-p mixing switch**: π2p sits *below* σ2p for Li₂-N₂ (and CO, NO), and
  σ2p drops below π2p for O₂, F₂, Ne₂;
* **Aufbau + Hund filling**, so half-filled degenerate shells produce unpaired
  electrons;
* **bond order**, **HOMO/LUMO** and **magnetism**.

Electrons are counted from atomic numbers + overall charge (not RDKit bond
orders): both the chemically correct way to populate an MO diagram and far more
robust for radicals (NO), ions (O₂⁺) and noble gases (He₂). Verified results:
O₂ → bond order 2, **paramagnetic** (2 unpaired e⁻ in π*2p); N₂ → 3, diamagnetic;
CO → 3 with HOMO σ2p / LUMO π*2p; **B₂ paramagnetic, C₂ diamagnetic** (the cases
that only come out right because of the s-p mixing); He₂/Be₂/Ne₂ → bond order 0.

---

## 4. Architecture

A modular, single-purpose pipeline (`src/orbital_explorer/`):

| Module | Responsibility |
|---|---|
| `atomic.py` | Exact hydrogenic wavefunctions: real spherical harmonics, radial part, `Orbital` dataclass, node properties. |
| `isosurface.py` | Smooth surface meshes via `skimage.measure.marching_cubes` → `plotly.Mesh3d` with interpolated normals (the shared 3D rendering primitive). |
| `render_atomic.py` | Builds the interactive atomic-orbital figure (phase / density). |
| `molecule.py` | Input → RDKit molecule: SMILES, offline name table, OPSIN, PubChem fallback, `[n]annulene` generator; auto-classification (atom / diatomic / π-system / σ-only) with honest conjugation labelling. |
| `huckel.py` | Hückel solver: builds the matrix, diagonalises, assigns occupations (Aufbau+Hund), HOMO/LUMO, coefficients. |
| `mo_diagram.py` | π MO energy-level diagram (matplotlib) + summary numbers (delocalisation energy, aromaticity). |
| `render_mo.py` | 3D LCAO "orbital sharing" view (embed geometry, plane fit, summed p-orbital field). |
| `diatomic.py` | Curated diatomic MO engine: AO data, s-p mixing rule, templates, filling, `DiatomicResult`. |
| `diatomic_diagram.py` | Three-column diatomic MO diagram (matplotlib), reusing the π-diagram theme. |
| `app.py` | The Dash multi-tab UI, callbacks, dark theme. |

`tests/` holds the pytest suite; `docs/make_figures.py` regenerates the README
images offline from the same physics code.

---

## 5. Key design decisions

* **Methods over machinery.** Analytic atomic orbitals + Hückel + a curated
  diatomic model give correct-or-qualitative chemistry with zero heavy
  dependencies: the Windows/pip constraint made this the right call, not a
  compromise.
* **Mesh rendering, not voxels.** Early versions embedded the full volume in the
  Plotly figure (~10-12 MB each, blocky). Switching to marching-cubes meshes made
  surfaces smooth *and* shrank figures ~10-48×.
* **Translucency.** Plotly `Mesh3d` transparency only reads as see-through at low
  opacity (≈0.3) with flat (low-specular) lighting; smoothness comes from the
  mesh geometry, not glossy shading. This keeps the molecular skeleton visible
  *through* the orbital lobes.
* **Robust name resolution.** SMILES → offline name table → OPSIN → PubChem REST,
  so most chemical names work online while SMILES and common names work offline.
  Special handling: bare `annulene` → [18]annulene (the prototypical aromatic
  one), `[n]annulene` generates the correct ring; `CO`/`NO` are read as the
  diatomics (not the methanol/hydroxylamine SMILES).
* **Honest chemistry.** Classification distinguishes a truly conjugated system
  from a single isolated π bond (acetone, ethylene) or several separate π units
  (1,4-pentadiene), and flags orthogonal-π (allene/alkyne) and non-planar systems
  as approximate rather than silently showing a misleading picture.
* **Fail loud, not blank.** Heavy callbacks are wrapped so an unembeddable
  molecule or bad input shows a clear message instead of a silent blank graph.

---

## 6. Scope and honest limitations

**Covered:** all hydrogenic atomic orbitals (s/p/d/f, any n, any Z); Hückel for
conjugated π-systems including heteroatoms, ions and radicals; qualitative MO
diagrams for H-Ne diatomics (homo- and heteronuclear) plus H-X molecules.

**Not covered / approximate:**
* Quantitative (ab-initio / kJ·mol⁻¹) energies: everything is qualitative.
* σ-framework MOs of polyatomic molecules.
* Period-3+ d-orbital diatomics; transition-metal complexes.
* Orthogonal π systems (allenes, alkynes): only the in-plane set is shown.
* Non-planar conjugated systems: a single p-orbital axis is used and flagged as
  approximate.

These limits are surfaced *in the app* (the classification note) rather than
hidden, because a teaching tool that quietly shows wrong chemistry is worse than
one that says "this case is approximate."

---

## 7. Testing

`python -m pytest`: **39 tests**, all checking physics directly rather than
comparing to a saved image:

* **Atomic:** normalisation `∫|ψ|²=1`, orthogonality, radial/angular node counts.
* **Hückel:** eigenvalues vs analytic textbook values (benzene `α±β, α±2β`;
  butadiene `±0.618, ±1.618`), electron counts for allyl cation/radical/anion,
  cyclobutadiene as a triplet diradical, aromatic cyclopentadienyl anion,
  orthonormal coefficient sets.
* **Diatomic:** bond orders and magnetism (O₂ paramagnetic with 2 unpaired e⁻,
  He₂ bond order 0), the N₂ vs O₂ ordering switch, CO bond order 3, HF, O₂⁺.

---

## 8. Install & run

Requires Python 3.11+.

```bash
pip install -r requirements.txt   # numpy, scipy, scikit-image, rdkit, plotly, dash, matplotlib, requests
python run_app.py                 # then open http://127.0.0.1:8050
python -m pytest                  # run the test suite
```

Things to try: `3dz2`, `4fxyz` (atomic); `benzene`, `naphthalene`, `pyridine`,
`C=CC=C`, `[18]annulene`, `azulene` (π-systems); `O2`, `N2`, `CO`, `NO`, `HF`,
`B2`, `He2` (diatomics).

**Tech stack:** NumPy · SciPy · scikit-image · RDKit · Plotly · Dash · matplotlib · requests.

---

## 9. How it was built (development history)

1. **Atomic engine**: exact hydrogenic wavefunctions; verified normalisation and
   node counts.
2. **Atomic 3D viewer**: Plotly isosurfaces, phase/density modes.
3. **Molecule input + classification**: RDKit, offline names.
4. **Hückel solver**: validated against analytic eigenvalues.
5. **MO energy diagram**: matplotlib ladder with HOMO/LUMO, gap, delocalisation,
   aromaticity.
6. **3D LCAO "orbital sharing"** view.
7. **Dash multi-tab app** unifying everything.
8. **Restructure** into `orgo-sim/ChemResolve` + `orgo-sim/orbital-explorer`;
   published to GitHub.
9. **Smoothness**: marching-cubes meshes (smaller + smoother figures).
10. **Translucency**: tuned opacity/lighting so lobes are see-through.
11. **Robust names**: PubChem fallback; `[n]annulene` family; `annulene` →
    [18]annulene fix.
12. **Honest chemistry, crash-proofing, and diatomic MO diagrams**: the latest
    round: correct conjugation labels, graceful error handling, and the full
    diatomic MO feature.

---

## 10. Future work

* Diatomic **3D σ/π orbital** rendering (currently only the energy diagram).
* Better heteroatom Hückel parameters; explicit handling of charged heteroatoms
  (e.g. pyridinium).
* Performance for large π-systems and high-n orbitals (grid + marching-cubes
  cost).
* Optional export of figures/data; colour-blind-friendly phase palette.
