"""Launch the Orbital Explorer app.

    python run_app.py

Then open http://127.0.0.1:8050 in your browser.

If the package is installed (``pip install -e .``), you can also run::

    orbital-explorer
    python -m orbital_explorer.app

For production, prefer gunicorn via ``wsgi.py`` (see Dockerfile / README).
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

# Make src/ importable without an editable install (dev convenience).
sys.path.insert(0, str(Path(__file__).resolve().parent / "src"))

from orbital_explorer.app import main  # noqa: E402

if __name__ == "__main__":
    # Local runner defaults to loopback; override with HOST/PORT if needed.
    os.environ.setdefault("HOST", "127.0.0.1")
    os.environ.setdefault("PORT", "8050")
    main()
