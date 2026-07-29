"""Build the pyCSAMT 2.0 social launch video, teaser GIF, and cover."""

from __future__ import annotations

from pathlib import Path

import imageio_ffmpeg
from PIL import Image, ImageDraw, ImageEnhance, ImageFilter, ImageFont

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "social_media" / "pycsamt_v2_launch"
SIZE = 1080
FPS = 24

NAVY = (5, 18, 32)
CYAN = (35, 201, 220)
ORANGE = (246, 145, 52)
WHITE = (247, 250, 252)
MUTED = (184, 205, 218)

ASSETS = {
    "logo": ROOT / "docs/source/_static/logo/pycsamt-v2.png",
    "hero": ROOT / "docs/source/_static/images/home/hero-carousel.webp",
    "data": ROOT
    / "docs/source/_static/applications/web/home-dashboard-loaded.png",
    "qc": ROOT / "docs/source/_static/applications/web/quality-control.png",
    "correction": ROOT / "docs/source/_static/applications/web/correction.png",
    "ai3d": ROOT
    / "docs/source/images/tutorials/essential_3d_ai_inversion/l18_ai3d_topography_block.png",
    "agents": ROOT / "docs/source/_static/applications/web/agents-chat.png",
    "pipeline": ROOT
    / "docs/source/images/tutorials/run_pipeline_from_config/pipeline_configured_chain.png",
    "map3d": ROOT
    / "docs/source/_static/applications/mapview/inversion-3d-conductive.png",
    "hydro": ROOT
    / "docs/source/images/user_guide/interpretation/hydro_aquifer_characterization.png",
    "webgif": ROOT
    / "docs/source/_static/applications/web/web-walkthrough.gif",
    "mapgif": ROOT
    / "docs/source/_static/applications/mapview/mapview-walkthrough.gif",
}


def font(size: int, bold: bool = False) -> ImageFont.FreeTypeFont:
    names = [
        "C:/Windows/Fonts/arialbd.ttf"
        if bold
        else "C:/Windows/Fonts/arial.ttf",
        "C:/Windows/Fonts/segoeuib.ttf"
        if bold
        else "C:/Windows/Fonts/segoeui.ttf",
    ]
    for name in names:
        if Path(name).exists():
            return ImageFont.truetype(name, size)
    return ImageFont.load_default()


def load(path: Path) -> Image.Image:
    return Image.open(path).convert("RGB")


def fit_cover(
    im: Image.Image, size: int = SIZE, zoom: float = 1.0, pan: float = 0.5
) -> Image.Image:
    w, h = im.size
    scale = max(size / w, size / h) * zoom
    resized = im.resize(
        (round(w * scale), round(h * scale)), Image.Resampling.LANCZOS
    )
    left_max = max(0, resized.width - size)
    top_max = max(0, resized.height - size)
    left = round(left_max * pan)
    top = round(top_max * (1.0 - pan))
    return resized.crop((left, top, left + size, top + size))


def tint_background(im: Image.Image, darkness: float = 0.44) -> Image.Image:
    im = ImageEnhance.Color(im).enhance(0.84)
    overlay = Image.new("RGB", im.size, NAVY)
    return Image.blend(im, overlay, darkness)


def gradient_panel(base: Image.Image, top: int = 610) -> Image.Image:
    rgba = base.convert("RGBA")
    panel = Image.new("RGBA", rgba.size, (0, 0, 0, 0))
    px = panel.load()
    for y in range(top, SIZE):
        alpha = int(235 * ((y - top) / max(1, SIZE - top)) ** 0.62)
        for x in range(SIZE):
            px[x, y] = (*NAVY, alpha)
    return Image.alpha_composite(rgba, panel).convert("RGB")


def centered(
    draw: ImageDraw.ImageDraw, text: str, y: int, fnt, fill=WHITE, spacing=8
):
    box = draw.multiline_textbbox(
        (0, 0), text, font=fnt, spacing=spacing, align="center"
    )
    draw.multiline_text(
        ((SIZE - (box[2] - box[0])) / 2, y),
        text,
        font=fnt,
        fill=fill,
        spacing=spacing,
        align="center",
    )


def add_brand(draw: ImageDraw.ImageDraw):
    draw.rounded_rectangle(
        (54, 52, 280, 104),
        18,
        fill=(*NAVY, 225),
        outline=(*CYAN, 220),
        width=2,
    )
    draw.text((76, 64), "pyCSAMT 2.0", font=font(27, True), fill=WHITE)


def scene_frame(
    asset: Path, title: str, subtitle: str, progress: float
) -> Image.Image:
    src = load(asset)
    zoom = 1.02 + 0.05 * progress
    bg = fit_cover(src, zoom=zoom, pan=0.25 + 0.5 * progress)
    bg = gradient_panel(bg, 570)
    draw = ImageDraw.Draw(bg, "RGBA")
    add_brand(draw)
    draw.rounded_rectangle(
        (58, 665, 1022, 982),
        28,
        fill=(*NAVY, 222),
        outline=(*CYAN, 165),
        width=2,
    )
    draw.rectangle((92, 708, 200, 716), fill=ORANGE)
    draw.text((92, 742), title, font=font(52, True), fill=WHITE)
    draw.multiline_text(
        (92, 820), subtitle, font=font(31), fill=MUTED, spacing=13
    )
    return bg


def title_frame(progress: float) -> Image.Image:
    bg = tint_background(
        fit_cover(load(ASSETS["hero"]), zoom=1.03 + 0.04 * progress), 0.56
    )
    bg = bg.filter(ImageFilter.GaussianBlur(1.1))
    logo = Image.open(ASSETS["logo"]).convert("RGBA")
    logo.thumbnail((650, 180), Image.Resampling.LANCZOS)
    canvas = bg.convert("RGBA")
    canvas.alpha_composite(logo, ((SIZE - logo.width) // 2, 150))
    draw = ImageDraw.Draw(canvas, "RGBA")
    centered(draw, "AU-DELÀ DES CHATBOTS", 420, font(61, True), CYAN)
    centered(
        draw,
        "L’intelligence artificielle au service\nde l’exploration du sous-sol",
        520,
        font(46, True),
    )
    draw.rounded_rectangle((265, 745, 815, 825), 30, fill=(*ORANGE, 235))
    centered(
        draw,
        "LA GÉOPHYSIQUE ENTRE DANS UNE NOUVELLE ÈRE",
        771,
        font(20, True),
        NAVY,
    )
    return canvas.convert("RGB")


def applications_frame(progress: float) -> Image.Image:
    bg = tint_background(
        fit_cover(load(ASSETS["hydro"]), zoom=1.02 + 0.04 * progress), 0.52
    )
    draw = ImageDraw.Draw(bg, "RGBA")
    add_brand(draw)
    centered(draw, "DES APPLICATIONS CONCRÈTES", 165, font(49, True), CYAN)
    items = [
        "EAUX SOUTERRAINES",
        "EXPLORATION MINIÈRE",
        "GÉOTHERMIE",
        "ENVIRONNEMENT & GÉOTECHNIQUE",
    ]
    y = 335
    for i, item in enumerate(items):
        color = CYAN if i % 2 == 0 else ORANGE
        draw.rounded_rectangle(
            (155, y, 925, y + 90),
            24,
            fill=(*NAVY, 220),
            outline=(*color, 210),
            width=3,
        )
        centered(draw, item, y + 27, font(29, True), WHITE)
        y += 120
    return bg


def cta_frame(progress: float) -> Image.Image:
    bg = Image.new("RGB", (SIZE, SIZE), NAVY)
    draw = ImageDraw.Draw(bg, "RGBA")
    for radius in range(520, 60, -40):
        alpha = max(3, int(32 * (1 - radius / 520)))
        draw.ellipse(
            (
                SIZE // 2 - radius,
                SIZE // 2 - radius,
                SIZE // 2 + radius,
                SIZE // 2 + radius,
            ),
            outline=(*CYAN, alpha),
            width=3,
        )
    logo = Image.open(ASSETS["logo"]).convert("RGBA")
    logo.thumbnail((720, 190), Image.Resampling.LANCZOS)
    canvas = bg.convert("RGBA")
    canvas.alpha_composite(logo, ((SIZE - logo.width) // 2, 190))
    draw = ImageDraw.Draw(canvas, "RGBA")
    centered(
        draw, "PYTHON • OPEN SOURCE • SCIENTIFIC AI", 445, font(28, True), CYAN
    )
    centered(
        draw,
        "De la donnée de terrain\nà la décision scientifique",
        545,
        font(49, True),
        WHITE,
    )
    draw.rounded_rectangle((250, 745, 830, 850), 32, fill=ORANGE)
    centered(draw, "DÉCOUVREZ pycsamt.org", 778, font(31, True), NAVY)
    centered(
        draw,
        "Construisons nos propres outils scientifiques.",
        920,
        font(25),
        MUTED,
    )
    return canvas.convert("RGB")


SCENES = [
    (4.5, lambda p: title_frame(p)),
    (
        5.0,
        lambda p: scene_frame(
            ASSETS["data"],
            "DES DONNÉES AU SOUS-SOL",
            "CSAMT • AMT • MT • TDEM\nUne plateforme scientifique intégrée",
            p,
        ),
    ),
    (
        5.5,
        lambda p: scene_frame(
            ASSETS["qc"],
            "CONTRÔLER & CORRIGER",
            "Qualité • Static shift • Strike\nTenseur de phase • Traçabilité",
            p,
        ),
    ),
    (
        6.0,
        lambda p: scene_frame(
            ASSETS["ai3d"],
            "INVERSER EN 1D, 2D & 3D",
            "Méthodes classiques • Deep learning\nRéseaux informés par la physique",
            p,
        ),
    ),
    (
        7.0,
        lambda p: scene_frame(
            ASSETS["agents"],
            "COLLABORER AVEC DES AGENTS",
            "Décrivez votre objectif en langage naturel\nDonnées → QC → Inversion → Rapport",
            p,
        ),
    ),
    (
        4.5,
        lambda p: scene_frame(
            ASSETS["map3d"],
            "VOIR, COMPRENDRE, DÉCIDER",
            "Cartes • Profils • Sections\nInterprétation et visualisation 3D",
            p,
        ),
    ),
    (4.5, lambda p: applications_frame(p)),
    (5.0, lambda p: cta_frame(p)),
]


def iter_video_frames():
    for duration, renderer in SCENES:
        count = round(duration * FPS)
        for i in range(count):
            p = i / max(1, count - 1)
            frame = renderer(p)
            fade = min(1.0, i / 10, (count - i - 1) / 10)
            if fade < 1:
                frame = Image.blend(
                    Image.new("RGB", frame.size, NAVY), frame, max(0, fade)
                )
            yield frame


def write_mp4(path: Path):
    writer = imageio_ffmpeg.write_frames(
        str(path),
        (SIZE, SIZE),
        fps=FPS,
        codec="libx264",
        pix_fmt_in="rgb24",
        pix_fmt_out="yuv420p",
        output_params=[
            "-crf",
            "20",
            "-preset",
            "medium",
            "-movflags",
            "+faststart",
        ],
    )
    writer.send(None)
    try:
        for frame in iter_video_frames():
            writer.send(frame.tobytes())
    finally:
        writer.close()


def write_cover(path: Path):
    frame = title_frame(0.35)
    draw = ImageDraw.Draw(frame, "RGBA")
    draw.rounded_rectangle(
        (330, 890, 750, 958),
        22,
        fill=(*NAVY, 225),
        outline=(*CYAN, 210),
        width=2,
    )
    centered(draw, "VERSION 2.0 DISPONIBLE", 910, font(24, True), WHITE)
    frame.save(path, quality=94, optimize=True)


def write_teaser(path: Path):
    frames = []
    teaser = [
        (title_frame, 18),
        (
            lambda p: scene_frame(
                ASSETS["agents"],
                "DES AGENTS SCIENTIFIQUES",
                "Du contrôle qualité à l’interprétation",
                p,
            ),
            18,
        ),
        (
            lambda p: scene_frame(
                ASSETS["ai3d"],
                "INVERSIONS 1D • 2D • 3D",
                "Géophysique classique & intelligence artificielle",
                p,
            ),
            18,
        ),
        (cta_frame, 22),
    ]
    for renderer, count in teaser:
        for i in range(count):
            im = renderer(i / max(1, count - 1)).resize(
                (720, 720), Image.Resampling.LANCZOS
            )
            frames.append(
                im.convert("P", palette=Image.Palette.ADAPTIVE, colors=160)
            )
    frames[0].save(
        path,
        save_all=True,
        append_images=frames[1:],
        duration=150,
        loop=0,
        optimize=True,
        disposal=2,
    )


def main():
    missing = [str(p) for p in ASSETS.values() if not p.exists()]
    if missing:
        raise FileNotFoundError("Missing assets:\n" + "\n".join(missing))
    OUT.mkdir(parents=True, exist_ok=True)
    write_cover(OUT / "pycsamt_v2_cover.jpg")
    write_teaser(OUT / "pycsamt_v2_teaser.gif")
    write_mp4(OUT / "pycsamt_v2_demo_fr.mp4")
    print(f"Created launch assets in {OUT}")


if __name__ == "__main__":
    main()
