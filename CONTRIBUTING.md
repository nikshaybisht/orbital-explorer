# Contributing

Thanks for helping with a small MO theory visualiser.

## Setup

Python 3.11 to 3.13.

```bash
git clone https://github.com/nikshaybisht/orbital-explorer.git
cd orbital-explorer
python -m venv .venv

# Windows
.venv\Scripts\activate
# macOS / Linux
source .venv/bin/activate

pip install -U pip
pip install -e ".[dev]"
```

### Windows + RDKit

PyPI has Windows wheels. Prefer:

```bash
pip install -e ".[dev]"
```

If `pip install rdkit` fails (usually a bad/mismatched Python install):

1. Use 64-bit CPython 3.11 or 3.12 from [python.org](https://www.python.org/downloads/).
2. `python -m pip install -U pip setuptools wheel`
3. `pip install "rdkit>=2024.3.1"`
4. Fallback: [conda-forge](https://anaconda.org/conda-forge/rdkit)  
   `conda install -c conda-forge rdkit`, then `pip install -e ".[dev]" --no-deps` and install the other packages separately.

## Run

```bash
python run_app.py
# or
orbital-explorer
```

http://127.0.0.1:8050

## Tests

```bash
python -m pytest
```

Keep the physics tests green: orbital normalisation, Huckel textbook eigenvalues, diatomic bond orders / magnetism (O2 paramagnetic, etc.).

## Pull requests

1. One fix or feature per PR.
2. Add or update tests when behaviour changes.
3. Prefer a correct chemical label over a prettier wrong one.
4. Run `pytest` before opening the PR.

## Conduct

Be decent. Correct chemistry beats clever code.
