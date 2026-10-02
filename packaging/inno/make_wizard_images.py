# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Generate the branded Inno Setup wizard images for pycsamt-desktop.

Replaces Inno's stock "cardboard box + CD" artwork with pyCSAMT's own
identity (the navy/amber palette of ``pycsamt-v2-splash.png``, the circular
wave symbol, and the subsurface-resistivity hero strip
``hero-carousel.webp``):

* ``wizard_large_<pct>.png`` -- the tall side panel shown on the Welcome and
  Finished pages (164x314 px at 100 % DPI);
* ``wizard_small_<pct>.png`` -- the corner badge shown on every inner page
  (55x55 px at 100 % DPI, transparent background so it sits cleanly on both
  the light and dark ``windows11`` wizard styles).

One file per DPI step (100/125/150/175/200/250 %) is written so Inno never
has to upscale a bitmap -- ``pycsamt_desktop.iss`` lists them all and Setup
picks the closest match at run time.

Usage (needs Pillow + PySide6, both already in any desktop build env)::

    python packaging/inno/make_wizard_images.py

Outputs land in ``packaging/inno/assets/`` and are committed, so building
the installer itself never depends on this script having been run.
"""

from __future__ import annotations

import sys
from pathlib import Path

from PIL import Image, ImageDraw, ImageFilter, ImageFont

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parent.parent
ICONS = REPO_ROOT / "pycsamt" / "app" / "desktop" / "resources" / "icons"
HERO = ICONS / "hero-carousel.webp"
SYMBOL_SVG = ICONS / "pycsamt-v2-symbol.svg"
OUT_DIR = HERE / "assets"

SCALES = (100, 125, 150, 175, 200, 250)
LARGE_BASE = (164, 314)
SMALL_BASE = 55

# Palette sampled from pycsamt-v2-splash.png / the symbol SVG.
NAVY_TOP = (9, 30, 82)
NAVY_BOTTOM = (22, 64, 150)
AMBER = (247, 165, 49)
SKY = (160, 196, 255)
WHITE = (255, 255, 255)

_FONT_DIRS = (Path("C:/Windows/Fonts"), Path("/usr/share/fonts/truetype/dejavu"))


def _font(size: int, bold: bool = False) -> ImageFont.FreeTypeFont:
    names = (
        ("segoeuib.ttf", "DejaVuSans-Bold.ttf")
        if bold
        else ("segoeui.ttf", "DejaVuSans.ttf")
    )
    for d in _FONT_DIRS:
        for n in names:
            if (d / n).exists():
                return ImageFont.truetype(str(d / n), size)
    return ImageFont.load_default()


def _render_symbol(px: int) -> Image.Image:
    """Rasterize the circular wave symbol at ``px`` x ``px`` via Qt's SVG
    renderer (Pillow has no SVG support)."""
    from PySide6.QtCore import QRectF, Qt
    from PySide6.QtGui import QGuiApplication, QImage, QPainter
    from PySide6.QtSvg import QSvgRenderer

    _app = QGuiApplication.instance() or QGuiApplication(sys.argv[:1])  # noqa: F841
    renderer = QSvgRenderer(str(SYMBOL_SVG))
    vb = renderer.viewBoxF()
    # Oversample x4 then downscale with LANCZOS for clean small badges.
    big = px * 4
    img = QImage(big, big, QImage.Format_ARGB32)
    img.fill(Qt.transparent)
    painter = QPainter(img)
    painter.setRenderHint(QPainter.Antialiasing)
    s = big / max(vb.width(), vb.height())
    w, h = vb.width() * s, vb.height() * s
    renderer.render(painter, QRectF((big - w) / 2, (big - h) / 2, w, h))
    painter.end()
    buf = img.constBits().tobytes()
    pil = Image.frombuffer("RGBA", (big, big), buf, "raw", "BGRA", 0, 1)
    # Crop to the *coloured* content: the viewBox has transparent margins
    # and a stray white stroked element that would off-centre the disc.
    ink = Image.eval(pil.convert("RGB").convert("L"), lambda v: 255 if v < 230 else 0)
    ink = Image.composite(ink, Image.new("L", pil.size, 0), pil.getchannel("A"))
    pil = pil.crop(ink.getbbox())
    # Drop anything outside the disc (that stray stroke) with a round mask.
    mask = Image.new("L", pil.size, 0)
    ImageDraw.Draw(mask).ellipse((0, 0, pil.width - 1, pil.height - 1), fill=255)
    pil.putalpha(Image.composite(pil.getchannel("A"), mask, mask))
    side = max(pil.size)
    sq = Image.new("RGBA", (side, side), (0, 0, 0, 0))
    sq.paste(pil, ((side - pil.width) // 2, (side - pil.height) // 2))
    return sq.resize((px, px), Image.LANCZOS)


def _vgradient(size, top, bottom) -> Image.Image:
    w, h = size
    grad = Image.new("RGB", (1, h))
    for y in range(h):
        t = y / max(h - 1, 1)
        grad.putpixel((0, y), tuple(int(a + (b - a) * t) for a, b in zip(top, bottom)))
    return grad.resize((w, h))


def make_large(scale: int) -> Image.Image:
    w, h = (round(v * scale / 100) for v in LARGE_BASE)
    k = scale / 100
    canvas = _vgradient((w, h), NAVY_TOP, NAVY_BOTTOM).convert("RGBA")

    # Bottom ~48 %: the hero strip's layered-resistivity section, cropped
    # to the panel's aspect so its amber layer sits mid-panel.
    hero = Image.open(HERO).convert("RGB")
    band_h = int(h * 0.47)
    crop_w = int(hero.height * w / band_h)
    x0 = min(int(hero.width * 0.70), hero.width - crop_w)
    strip = hero.crop((x0, 0, x0 + crop_w, hero.height)).resize(
        (w, band_h), Image.LANCZOS
    )
    # The hero art's top edge is near-white sky: fade it in from the navy
    # (smoothstep) so the amber layer emerges instead of a pasted rectangle.
    fade = Image.new("L", (w, band_h))
    for y in range(band_h):
        t = min(1.0, y / (band_h * 0.42))
        fade.paste(int(255 * t * t * (3 - 2 * t)), (0, y, w, y + 1))
    canvas.paste(strip.convert("RGBA"), (0, h - band_h), fade)

    # Soft glow behind the symbol, echoing the splash screen.
    sym_px = int(84 * k)
    cx, cy = w // 2, int(78 * k)
    glow_layer = Image.new("RGBA", (w, h), (0, 0, 0, 0))
    gd = ImageDraw.Draw(glow_layer)
    r = int(sym_px * 0.62)
    gd.ellipse((cx - r, cy - r, cx + r, cy + r), fill=(120, 170, 255, 150))
    glow_layer = glow_layer.filter(ImageFilter.GaussianBlur(14 * k))
    canvas = Image.alpha_composite(canvas, glow_layer)

    # White disc behind the symbol so its white waves read as waves.
    disc = Image.new("RGBA", (w, h), (0, 0, 0, 0))
    ImageDraw.Draw(disc).ellipse(
        (cx - sym_px // 2 - int(3 * k), cy - sym_px // 2 - int(3 * k),
         cx + sym_px // 2 + int(3 * k), cy + sym_px // 2 + int(3 * k)),
        fill=WHITE + (255,),
    )
    canvas = Image.alpha_composite(canvas, disc)
    sym = _render_symbol(sym_px)
    canvas.alpha_composite(sym, (cx - sym_px // 2, cy - sym_px // 2))

    # Wordmark: "py" sky-blue, "CSAMT" amber -- same split as the splash.
    d = ImageDraw.Draw(canvas)
    f_word = _font(int(27 * k), bold=True)
    py_w = d.textlength("py", font=f_word)
    cs_w = d.textlength("CSAMT", font=f_word)
    tx = (w - (py_w + cs_w)) / 2
    ty = int(128 * k)
    d.text((tx, ty), "py", font=f_word, fill=SKY)
    d.text((tx + py_w, ty), "CSAMT", font=f_word, fill=AMBER)

    f_sub = _font(int(9.5 * k))
    for i, line in enumerate(("MT · AMT · CSAMT · CSEM", "Desktop Suite")):
        lw = d.textlength(line, font=f_sub)
        d.text(((w - lw) / 2, int((166 + 13 * i) * k)), line, font=f_sub,
               fill=(214, 226, 250) if i == 0 else AMBER)
    return canvas.convert("RGB")


def make_small(scale: int) -> Image.Image:
    px = round(SMALL_BASE * scale / 100)
    pad = max(2, round(px * 0.06))
    inner = px - 2 * pad
    out = Image.new("RGBA", (px, px), (0, 0, 0, 0))
    disc = Image.new("RGBA", (px * 4, px * 4), (0, 0, 0, 0))
    ImageDraw.Draw(disc).ellipse((0, 0, px * 4 - 1, px * 4 - 1), fill=WHITE + (255,))
    disc = disc.resize((px, px), Image.LANCZOS)
    out.alpha_composite(disc)
    out.alpha_composite(_render_symbol(inner), (pad, pad))
    return out


def main() -> int:
    OUT_DIR.mkdir(exist_ok=True)
    for s in SCALES:
        make_large(s).save(OUT_DIR / f"wizard_large_{s}.png", optimize=True)
        make_small(s).save(OUT_DIR / f"wizard_small_{s}.png", optimize=True)
    print(f"Wrote {2 * len(SCALES)} wizard images to {OUT_DIR}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
