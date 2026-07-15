"""Build launch demo media: GIF + interactive HTML gallery for GitHub Pages.

    python docs/build_demo_media.py

Outputs:
  docs/images/demo.gif
  docs/demo/index.html
  docs/demo/atomic_*.html
  docs/demo/mo_benzene_homo.html   (if RDKit works)
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image, ImageDraw, ImageFont

ROOT = Path(__file__).resolve().parent.parent
SRC = ROOT / "src"
sys.path.insert(0, str(SRC))

from orbital_explorer.atomic import Orbital  # noqa: E402
from orbital_explorer.render_atomic import orbital_figure  # noqa: E402

IMAGES = ROOT / "docs" / "images"
DEMO = ROOT / "docs" / "demo"
DEMO.mkdir(parents=True, exist_ok=True)
IMAGES.mkdir(parents=True, exist_ok=True)

BG = (11, 14, 23)
FG = (232, 237, 246)
ACCENT = (90, 160, 255)
GOLD = (244, 196, 48)

# Pre-baked interactive demos (no server needed on GitHub Pages)
ATOMIC_DEMOS: list[tuple[str, int, str, str]] = [
    # filename stem, n, label, short blurb
    ("atomic_2pz", 2, "pz", "Classic dumbbell. Phase (blue/gold) is what bonds or cancels."),
    ("atomic_2px", 2, "px", "Same shape as 2pz, rotated onto x."),
    ("atomic_3dz2", 3, "dz2", "d-orbital with a doughnut; two nodes around the axis."),
    ("atomic_3dxy", 3, "dxy", "Four lobes in the xy plane."),
    ("atomic_4fxyz", 4, "fxyz", "f-orbital example for higher n."),
]


def _font(size: int):
    for name in (
        "C:/Windows/Fonts/segoeui.ttf",
        "C:/Windows/Fonts/arial.ttf",
        "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
        "/System/Library/Fonts/Helvetica.ttc",
    ):
        p = Path(name)
        if p.exists():
            try:
                return ImageFont.truetype(str(p), size=size)
            except OSError:
                continue
    return ImageFont.load_default()


def _slide(base: Image.Image, title: str, subtitle: str, w: int = 1280, h: int = 720) -> Image.Image:
    """Letterbox a PNG onto a dark 1280x720 slide with captions."""
    canvas = Image.new("RGB", (w, h), BG)
    draw = ImageDraw.Draw(canvas)
    title_font = _font(36)
    sub_font = _font(22)
    draw.text((40, 28), title, fill=ACCENT, font=title_font)
    draw.text((40, 78), subtitle, fill=FG, font=sub_font)

    area = (40, 120, w - 40, h - 40)
    aw, ah = area[2] - area[0], area[3] - area[1]
    img = base.convert("RGB")
    img.thumbnail((aw, ah), Image.Resampling.LANCZOS)
    x = area[0] + (aw - img.width) // 2
    y = area[1] + (ah - img.height) // 2
    canvas.paste(img, (x, y))
    draw.text(
        (40, h - 36),
        "Orbital Explorer | MIT | github.com/nikshaybisht/orbital-explorer",
        fill=(154, 166, 189),
        font=_font(16),
    )
    return canvas


def _title_card(w: int = 1280, h: int = 720) -> Image.Image:
    canvas = Image.new("RGB", (w, h), BG)
    draw = ImageDraw.Draw(canvas)
    draw.text((80, 220), "Orbital Explorer", fill=FG, font=_font(64))
    draw.text((80, 310), "MO theory you can actually look at", fill=ACCENT, font=_font(28))
    draw.text(
        (80, 380),
        "Hydrogenic AOs, Huckel pi MOs, diatomic diagrams",
        fill=(154, 166, 189),
        font=_font(22),
    )
    draw.text((80, 480), "Free, open source, MIT", fill=GOLD, font=_font(24))
    return canvas


def _end_card(w: int = 1280, h: int = 720) -> Image.Image:
    canvas = Image.new("RGB", (w, h), BG)
    draw = ImageDraw.Draw(canvas)
    draw.text((80, 240), "Try it", fill=ACCENT, font=_font(48))
    draw.text((80, 320), "github.com/nikshaybisht/orbital-explorer", fill=FG, font=_font(28))
    draw.text(
        (80, 380),
        "pip install -r requirements.txt && python run_app.py",
        fill=(154, 166, 189),
        font=_font(22),
    )
    draw.text((80, 440), "streamlit run streamlit_app.py", fill=(154, 166, 189), font=_font(22))
    draw.text(
        (80, 520),
        "nikshaybisht.github.io/orbital-explorer",
        fill=GOLD,
        font=_font(22),
    )
    return canvas


def build_gif() -> Path:
    frames: list[Image.Image] = []
    frames.append(_title_card())

    mapping = [
        ("atomic_orbitals.png", "1 · Atomic orbitals", "Exact hydrogenic psi: s / p / d / f phase lobes"),
        ("mo_3d_benzene.png", "2 · Benzene HOMO", "Huckel LCAO: p-orbitals sharing across the ring"),
        ("mo_diagram_benzene.png", "3 · Benzene MO ladder", "HOMO/LUMO, gap, aromaticity (4n+2)"),
        ("diatomic_o2.png", "4 · O2 MO diagram", "Bond order 2 · two unpaired e- · paramagnetic"),
    ]
    for name, title, sub in mapping:
        path = IMAGES / name
        if not path.exists():
            print(f"skip missing {path}")
            continue
        frames.append(_slide(Image.open(path), title, sub))

    frames.append(_end_card())

    out = IMAGES / "demo.gif"
    sequenced: list[Image.Image] = []
    durations: list[int] = []
    for i, fr in enumerate(frames):
        hold = 2500 if i in (0, len(frames) - 1) else 3500
        sequenced.append(fr)
        durations.append(hold)

    sequenced[0].save(
        out,
        save_all=True,
        append_images=sequenced[1:],
        duration=durations,
        loop=0,
        optimize=False,
    )
    print(f"wrote {out} ({out.stat().st_size // 1024} KB)")
    return out


def _write_plotly_html(fig, path: Path, title: str) -> None:
    """Wrap Plotly HTML with a dark chrome header + back link."""
    body = fig.to_html(include_plotlyjs="cdn", full_html=False, config={"displayModeBar": True})
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1" />
  <title>{title} · Orbital Explorer</title>
  <meta name="description" content="{title} — interactive 3D orbital from Orbital Explorer." />
  <meta property="og:title" content="{title} · Orbital Explorer" />
  <meta property="og:description" content="Interactive 3D orbital visualisation for undergrad MO theory." />
  <style>
    :root {{ color-scheme: dark; }}
    body {{ margin: 0; font-family: system-ui, sans-serif; background: #0b0e17; color: #e8edf6; }}
    header {{ padding: 0.75rem 1rem; border-bottom: 1px solid #2a3447;
             display: flex; flex-wrap: wrap; gap: 0.75rem; align-items: center; }}
    a {{ color: #5aa0ff; text-decoration: none; }}
    a:hover {{ text-decoration: underline; }}
    .title {{ font-weight: 600; }}
    .muted {{ color: #9aa6bd; font-size: 0.9rem; }}
    .plot {{ height: calc(100vh - 56px); }}
  </style>
</head>
<body>
  <header>
    <a href="index.html">← Demo home</a>
    <span class="title">{title}</span>
    <span class="muted">drag to rotate · scroll to zoom</span>
    <a class="muted" href="https://github.com/nikshaybisht/orbital-explorer" style="margin-left:auto">GitHub</a>
  </header>
  <div class="plot">{body}</div>
</body>
</html>
"""
    path.write_text(html, encoding="utf-8")
    print(f"wrote {path}")


def build_atomic_gallery() -> list[Path]:
    paths: list[Path] = []
    for stem, n, label, _blurb in ATOMIC_DEMOS:
        fig = orbital_figure(Orbital(n=n, label=label), mode="phase", n_grid=56)
        out = DEMO / f"{stem}.html"
        title = f"{n}{label}" if not label[0].isdigit() else label
        # Prefer chemist-facing name from Orbital
        orb = Orbital(n=n, label=label)
        _write_plotly_html(fig, out, orb.name)
        paths.append(out)
    # Keep legacy filename used by older links / README
    legacy = DEMO / "atomic_2pz.html"
    if (DEMO / "atomic_2pz.html").exists() or True:
        # atomic_2pz already written as stem atomic_2pz
        pass
    return paths


def build_benzene_mo_html() -> Path | None:
    """Interactive benzene HOMO if RDKit is available."""
    try:
        from orbital_explorer.huckel import solve
        from orbital_explorer.molecule import resolve
        from orbital_explorer.render_mo import mo_figure
    except Exception as exc:  # noqa: BLE001
        print(f"skip benzene MO HTML (import): {exc}")
        return None

    try:
        info = resolve("benzene")
        if info.mol is None or not info.pi_atom_indices:
            print("skip benzene MO HTML: resolve failed")
            return None
        res = solve(info.mol, info.pi_atom_indices)
        mo_idx = res.homo if res.homo is not None else 0
        fig = mo_figure(info.mol, res, mo_idx, n_grid=44)
        out = DEMO / "mo_benzene_homo.html"
        _write_plotly_html(fig, out, "Benzene HOMO (Huckel LCAO)")
        return out
    except Exception as exc:  # noqa: BLE001
        print(f"skip benzene MO HTML (runtime): {exc}")
        return None


def build_index(has_benzene: bool) -> Path:
    """Landing page for GitHub Pages (image paths use ../images for repo; CI rewrites)."""
    atomic_cards = []
    for stem, n, label, blurb in ATOMIC_DEMOS:
        orb = Orbital(n=n, label=label)
        atomic_cards.append(
            f"""      <a class="tile" href="{stem}.html">
        <strong>{orb.name}</strong>
        <span class="muted">{blurb}</span>
      </a>"""
        )
    atomic_grid = "\n".join(atomic_cards)

    benzene_tile = ""
    if has_benzene:
        benzene_tile = """
      <a class="tile highlight" href="mo_benzene_homo.html">
        <strong>Benzene HOMO (3D)</strong>
        <span class="muted">Interactive Huckel LCAO. Rotate the shared p-lobes.</span>
      </a>"""

    out = DEMO / "index.html"
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1" />
  <title>Orbital Explorer · demo</title>
  <meta name="description" content="Interactive 3D atomic and molecular orbitals for undergrad MO theory. Free open-source teaching tool." />
  <meta property="og:title" content="Orbital Explorer" />
  <meta property="og:description" content="3D AOs, Huckel pi MOs, diatomic diagrams for MO theory revision." />
  <meta property="og:type" content="website" />
  <meta property="og:url" content="https://nikshaybisht.github.io/orbital-explorer/" />
  <meta property="og:image" content="https://nikshaybisht.github.io/orbital-explorer/images/atomic_orbitals.png" />
  <meta name="twitter:card" content="summary_large_image" />
  <link rel="canonical" href="https://nikshaybisht.github.io/orbital-explorer/" />
  <style>
    :root {{ color-scheme: dark; }}
    body {{ margin: 0; font-family: system-ui, sans-serif; background: #0b0e17; color: #e8edf6; line-height: 1.55; }}
    a {{ color: #5aa0ff; }}
    header {{ padding: 2rem 1.25rem 1rem; max-width: 980px; margin: 0 auto; }}
    h1 {{ margin: 0 0 0.4rem; font-size: 2rem; letter-spacing: -0.02em; }}
    .muted {{ color: #9aa6bd; }}
    main {{ max-width: 980px; margin: 0 auto; padding: 0 1.25rem 3rem; }}
    .card {{ background: #121726; border: 1px solid #2a3447; border-radius: 12px; padding: 1rem 1.1rem; margin: 1rem 0; }}
    img, video {{ max-width: 100%; border-radius: 8px; display: block; }}
    .row {{ display: grid; gap: 1rem; grid-template-columns: repeat(auto-fit, minmax(260px, 1fr)); }}
    .tiles {{ display: grid; gap: 0.75rem; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); }}
    .tile {{ display: flex; flex-direction: column; gap: 0.35rem; padding: 0.9rem 1rem;
            background: #0b0e17; border: 1px solid #2a3447; border-radius: 10px;
            text-decoration: none; color: inherit; transition: border-color 0.15s, transform 0.15s; }}
    .tile:hover {{ border-color: #5aa0ff; transform: translateY(-1px); }}
    .tile.highlight {{ border-color: #5aa0ff55; }}
    .tile strong {{ color: #e8edf6; }}
    .btn {{ display: inline-block; background: #5aa0ff; color: #04101f; text-decoration: none;
           font-weight: 600; padding: 0.55rem 1rem; border-radius: 8px; margin: 0.3rem 0.4rem 0.3rem 0; }}
    .btn.secondary {{ background: transparent; color: #5aa0ff; border: 1px solid #5aa0ff; }}
    footer {{ max-width: 980px; margin: 0 auto; padding: 1rem 1.25rem 2rem; color: #9aa6bd; font-size: 0.9rem; }}
    code, pre {{ font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace; font-size: 0.88rem; }}
    pre {{ background: #0b0e17; padding: 1rem; border-radius: 8px; overflow: auto; border: 1px solid #2a3447; }}
    h2 {{ margin: 0.2rem 0 0.75rem; font-size: 1.2rem; }}
  </style>
</head>
<body>
  <header>
    <h1>Orbital Explorer</h1>
    <p class="muted">Undergrad MO theory, visual: hydrogenic AOs, Huckel pi systems, diatomic diagrams.</p>
    <p>
      <a class="btn" href="https://github.com/nikshaybisht/orbital-explorer">GitHub</a>
      <a class="btn secondary" href="https://github.com/nikshaybisht/orbital-explorer/releases/tag/v0.1.1">v0.1.1</a>
      <a class="btn secondary" href="https://github.com/nikshaybisht/orbital-explorer#install">Install</a>
    </p>
  </header>
  <main>
    <div class="card">
      <h2>Interactive gallery (no install)</h2>
      <p class="muted">Rotate in the browser. Full molecule solver still needs the local or Streamlit app.</p>
      <div class="tiles">
{atomic_grid}
{benzene_tile}
      </div>
    </div>

    <div class="card">
      <h2>Quick tour</h2>
      <p class="muted">Atomic orbitals → benzene HOMO → O2 (paramagnetic).</p>
      <img src="../images/demo.gif" alt="Demo tour: atomic orbitals, benzene HOMO, O2 MO diagram" width="960" height="540" />
    </div>

    <div class="row">
      <div class="card">
        <h3>Benzene HOMO</h3>
        <img src="../images/mo_3d_benzene.png" alt="Benzene HOMO 3D orbital sharing" />
        <p class="muted">Huckel LCAO: p orbitals sharing on the ring.</p>
      </div>
      <div class="card">
        <h3>O2 MO diagram</h3>
        <img src="../images/diatomic_o2.png" alt="O2 diatomic molecular orbital diagram" />
        <p class="muted">Bond order 2, two unpaired electrons, paramagnetic.</p>
      </div>
    </div>

    <div class="card">
      <h2>Run the full app</h2>
      <pre>git clone https://github.com/nikshaybisht/orbital-explorer.git
cd orbital-explorer
pip install -e ".[all]"
python run_app.py
# or: streamlit run streamlit_app.py</pre>
      <p class="muted">
        Dash on port 8050 · Streamlit for free cloud hosts ·
        <a href="https://github.com/nikshaybisht/orbital-explorer#docker">Docker</a> ·
        <a href="https://render.com/deploy?repo=https://github.com/nikshaybisht/orbital-explorer">Render</a>
      </p>
    </div>

    <div class="card">
      <h2>What this is</h2>
      <p class="muted">
        Teaching tool: exact hydrogenic wavefunctions, qualitative Huckel and diatomic MO diagrams.
        Not Gaussian / ORCA. Honest scope, 39 chemistry unit tests, CI on Ubuntu and Windows.
      </p>
    </div>
  </main>
  <footer>
    MIT · Nikshay Bisht ·
    <a href="https://github.com/nikshaybisht/orbital-explorer">source</a>
  </footer>
</body>
</html>
"""
    out.write_text(html, encoding="utf-8")
    print(f"wrote {out}")
    return out


def main() -> None:
    build_gif()
    build_atomic_gallery()
    has_benzene = build_benzene_mo_html() is not None
    build_index(has_benzene=has_benzene)
    print("demo media ready")


if __name__ == "__main__":
    main()
