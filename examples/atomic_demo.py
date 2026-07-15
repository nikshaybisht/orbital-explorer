"""Generate a gallery of interactive atomic-orbital pages.

Run it with no arguments to build a small representative gallery, or pass an
orbital to render just that one:

    python examples/atomic_demo.py            # gallery -> examples/output/index.html
    python examples/atomic_demo.py 3 dz2      # render 3dz2 only
    python examples/atomic_demo.py 4 fxyz density

Each page is a self-contained HTML file you can open in any browser and rotate.
"""

from __future__ import annotations

import sys
import webbrowser
from pathlib import Path

# Make the package importable without installing it.
SRC = Path(__file__).resolve().parent.parent / "src"
sys.path.insert(0, str(SRC))

from orbital_explorer.atomic import Orbital  # noqa: E402
from orbital_explorer.render_atomic import orbital_figure  # noqa: E402

OUTPUT = Path(__file__).resolve().parent / "output"

# A tour that shows off shells, phase lobes and radial nodes.
GALLERY = [
    (1, "s"), (2, "s"), (2, "pz"),
    (3, "dz2"), (3, "dxy"), (3, "dx2-y2"),
    (4, "fz3"), (4, "fxyz"),
]


def render_one(n: int, label: str, mode: str = "phase") -> Path:
    OUTPUT.mkdir(exist_ok=True)
    orb = Orbital(n=n, label=label)
    fig = orbital_figure(orb, mode=mode)
    path = OUTPUT / f"{orb.name}_{mode}.html"
    fig.write_html(path, include_plotlyjs="cdn")
    return path


def build_gallery() -> Path:
    OUTPUT.mkdir(exist_ok=True)
    links = []
    for n, label in GALLERY:
        path = render_one(n, label)
        links.append((Orbital(n=n, label=label).name, path.name))
    index = OUTPUT / "index.html"
    items = "\n".join(
        f'<li><a href="{fname}">{name}</a></li>' for name, fname in links
    )
    index.write_text(
        "<!doctype html><meta charset='utf-8'>"
        "<title>Orbital Explorer - atomic gallery</title>"
        "<style>body{background:#0b0e17;color:#e8edf6;font-family:system-ui;"
        "padding:2rem}a{color:#5aa0ff;font-size:1.2rem}li{margin:.4rem 0}</style>"
        "<h1>Atomic orbital gallery</h1><ul>" + items + "</ul>",
        encoding="utf-8",
    )
    return index


def main() -> None:
    args = sys.argv[1:]
    if len(args) >= 2:
        mode = args[2] if len(args) >= 3 else "phase"
        path = render_one(int(args[0]), args[1], mode)
        print("wrote", path)
        webbrowser.open(path.as_uri())
    else:
        index = build_gallery()
        print("wrote gallery ->", index)
        webbrowser.open(index.as_uri())


if __name__ == "__main__":
    main()
