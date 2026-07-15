"""Build launch demo media: GIF + interactive atomic HTML (no RDKit required).

    python docs/build_demo_media.py

Outputs:
  docs/images/demo.gif         : benzene HOMO + O2 diagram storyboard
  docs/demo/index.html         : public GitHub Pages landing page
  docs/demo/atomic_2pz.html    : interactive Plotly atomic orbital
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

    # Fit image into remaining area
    area = (40, 120, w - 40, h - 40)
    aw, ah = area[2] - area[0], area[3] - area[1]
    img = base.convert("RGB")
    img.thumbnail((aw, ah), Image.Resampling.LANCZOS)
    x = area[0] + (aw - img.width) // 2
    y = area[1] + (ah - img.height) // 2
    canvas.paste(img, (x, y))
    # footer
    draw.text((40, h - 36), "Orbital Explorer | MIT | github.com/nikshaybisht/orbital-explorer",
              fill=(154, 166, 189), font=_font(16))
    return canvas


def _title_card(w: int = 1280, h: int = 720) -> Image.Image:
    canvas = Image.new("RGB", (w, h), BG)
    draw = ImageDraw.Draw(canvas)
    draw.text((80, 220), "Orbital Explorer", fill=FG, font=_font(64))
    draw.text((80, 310), "MO theory you can actually look at", fill=ACCENT, font=_font(28))
    draw.text((80, 380), "Hydrogenic AOs, Huckel pi MOs, diatomic diagrams", fill=(154, 166, 189), font=_font(22))
    draw.text((80, 480), "Free, open source, MIT", fill=GOLD, font=_font(24))
    return canvas


def _end_card(w: int = 1280, h: int = 720) -> Image.Image:
    canvas = Image.new("RGB", (w, h), BG)
    draw = ImageDraw.Draw(canvas)
    draw.text((80, 240), "Try it", fill=ACCENT, font=_font(48))
    draw.text((80, 320), "github.com/nikshaybisht/orbital-explorer", fill=FG, font=_font(28))
    draw.text((80, 380), "pip install -r requirements.txt && python run_app.py", fill=(154, 166, 189), font=_font(22))
    draw.text((80, 440), "streamlit run streamlit_app.py", fill=(154, 166, 189), font=_font(22))
    return canvas


def build_gif() -> Path:
    frames: list[Image.Image] = []
    frames.append(_title_card())

    mapping = [
        ("atomic_orbitals.png", "1 · Atomic orbitals", "Exact hydrogenic psi: s / p / d / f phase lobes"),
        ("mo_3d_benzene.png", "2 · Benzene HOMO", "Hückel LCAO: p-orbitals sharing across the ring"),
        ("mo_diagram_benzene.png", "3 · Benzene MO ladder", "HOMO/LUMO, gap, aromaticity (4n+2)"),
        ("diatomic_o2.png", "4 · O₂ MO diagram", "Bond order 2 · two unpaired e⁻ · paramagnetic"),
    ]
    for name, title, sub in mapping:
        path = IMAGES / name
        if not path.exists():
            print(f"skip missing {path}")
            continue
        frames.append(_slide(Image.open(path), title, sub))

    frames.append(_end_card())

    out = IMAGES / "demo.gif"
    # ~75s at ~4s/frame for 6 frames is short; use 18-20 frames by duplication for pacing
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


def build_atomic_html() -> Path:
    fig = orbital_figure(Orbital(n=2, label="pz"), mode="phase", n_grid=64)
    out = DEMO / "atomic_2pz.html"
    fig.write_html(str(out), include_plotlyjs="cdn", full_html=True, config={"displayModeBar": True})
    print(f"wrote {out}")
    return out


def build_index() -> Path:
    out = DEMO / "index.html"
    # Relative paths for GitHub Pages project site
    html = """<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1" />
  <title>Orbital Explorer demo</title>
  <meta name="description" content="3D atomic and molecular orbitals for undergrad MO theory." />
  <style>
    :root { color-scheme: dark; }
    body { margin: 0; font-family: system-ui, sans-serif; background: #0b0e17; color: #e8edf6; line-height: 1.5; }
    a { color: #5aa0ff; }
    header { padding: 2rem 1.25rem 1rem; max-width: 960px; margin: 0 auto; }
    h1 { margin: 0 0 0.4rem; font-size: 2rem; }
    .muted { color: #9aa6bd; }
    main { max-width: 960px; margin: 0 auto; padding: 0 1.25rem 3rem; }
    .card { background: #121726; border: 1px solid #2a3447; border-radius: 12px; padding: 1rem; margin: 1rem 0; }
    img, video { max-width: 100%; border-radius: 8px; display: block; }
    .row { display: grid; gap: 1rem; grid-template-columns: repeat(auto-fit, minmax(260px, 1fr)); }
    .btn { display: inline-block; background: #5aa0ff; color: #04101f; text-decoration: none;
           font-weight: 600; padding: 0.55rem 1rem; border-radius: 8px; margin: 0.3rem 0.4rem 0.3rem 0; }
    .btn.secondary { background: transparent; color: #5aa0ff; border: 1px solid #5aa0ff; }
    footer { max-width: 960px; margin: 0 auto; padding: 1rem 1.25rem 2rem; color: #9aa6bd; font-size: 0.9rem; }
  </style>
</head>
<body>
  <header>
    <h1>Orbital Explorer</h1>
    <p class="muted">Undergrad MO theory, visual: AOs, Huckel, diatomics.</p>
    <p>
      <a class="btn" href="https://github.com/nikshaybisht/orbital-explorer">GitHub</a>
      <a class="btn secondary" href="https://github.com/nikshaybisht/orbital-explorer/releases/tag/v0.1.0">v0.1.0</a>
      <a class="btn secondary" href="atomic_2pz.html">Interactive 2p<sub>z</sub></a>
    </p>
  </header>
  <main>
    <div class="card">
      <h2>Quick tour</h2>
      <p class="muted">Atomic orbitals, then benzene HOMO, then O2 (paramagnetic).</p>
      <img src="../images/demo.gif" alt="Demo tour: atomic orbitals, benzene HOMO, O2 MO diagram" />
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
      <pre style="background:#0b0e17;padding:1rem;border-radius:8px;overflow:auto">pip install -r requirements.txt
python run_app.py
# or
streamlit run streamlit_app.py</pre>
      <p class="muted">Dash on port 8050, or Streamlit if you prefer that UI.</p>
    </div>
  </main>
  <footer>
    MIT · Nikshay Bisht · Huckel / hydrogenic teaching models, not production QC.
  </footer>
</body>
</html>
"""
    out.write_text(html, encoding="utf-8")
    print(f"wrote {out}")
    return out


def main() -> None:
    build_gif()
    build_atomic_html()
    build_index()
    print("demo media ready")


if __name__ == "__main__":
    main()
