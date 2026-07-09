"""Sphinx extension: compact "hub" landing page for the examples gallery.

sphinx-gallery renders ``examples/index`` as one long page holding every
thumbnail of every section (70+ cards).  This extension swaps that body
for a grid of one *animated* card per section — each card carousels
through its section's plot thumbnails with a typewriter caption — while
keeping the generated ``toctree`` tail verbatim so the sidebar and
navigation hierarchy are untouched.  The per-section pages that
sphinx-gallery already generates (``examples/<section>/index``) are the
click-through targets.

Everything is derived automatically, no hand-maintained lists:

* ``docs/examples/<section>/README.txt``   -> card title + blurb
* ``docs/examples/<section>/plot_*.py``    -> one slide per script,
  caption = first docstring line (same title sphinx-gallery uses)
* ``sphinx_gallery_conf["subsection_order"]`` -> card order
* thumbnails are copied to ``<outdir>/_gallery_hub/<section>/`` at the
  end of the build, so slide URLs never depend on how sphinx
  deduplicates names inside ``_images/``.

The manifest is embedded inline as JSON and rendered client-side by
``_static/js/gallery-hub.js`` (styled by ``_static/css/gallery-hub.css``).
"""

from __future__ import annotations

import json
import re
import shutil
from pathlib import Path

HUB_DOCNAME = "examples/index"
THUMB_OUTDIR = "_gallery_hub"

# Landing-page intro; the toctree tail of the generated page is appended
# below this (see _swap_source).  Mirrors docs/examples/README.txt, whose
# prose no longer reaches the rendered landing page.
HUB_RST = """\
:orphan:

.. _general_examples:

Examples
========

Runnable, fully rendered examples for the pyCSAMT tool families,
organised by section. Every script is executed while the documentation
is built, so each thumbnail and figure stays in sync with the code.

Each card below previews the plots of one section — open it to browse
that section's full example gallery.

For guided, end-to-end workflow walkthroughs see the
:doc:`tutorials </tutorials/index>`; for the complete API reference see
the :doc:`API docs </api/index>`.
"""

# --------------------------------------------------------------------------
# RST -> plain text for card blurbs / slide captions
# --------------------------------------------------------------------------

_ROLE_RE = re.compile(r":[\w.+:-]+:`([^`]+)`")


def _clean_rst(text: str) -> str:
    """Strip inline RST markup (roles, literals, emphasis) for display."""

    def _role(match: re.Match) -> str:
        inner = match.group(1)
        if "<" in inner:  # :ref:`Label <target>` -> Label
            inner = inner.split("<", 1)[0].strip()
        if inner.startswith("~"):  # ~pkg.mod.name -> name
            inner = inner[1:].rsplit(".", 1)[-1]
        return inner

    text = _ROLE_RE.sub(_role, text)
    text = re.sub(r"``([^`]+)``", r"\1", text)
    text = re.sub(r"\*\*([^*]+)\*\*", r"\1", text)
    text = re.sub(r"\*([^*]+)\*", r"\1", text)
    return re.sub(r"\s+", " ", text).strip()


_UNDERLINE_RE = re.compile(r"^([=\-~^\"'`#*+_])\1{2,}\s*$")


def _readme_meta(readme: Path) -> tuple[str, str]:
    """Return (title, first paragraph) from a section README.txt."""
    lines = readme.read_text(encoding="utf-8").splitlines()
    title, blurb_lines, i = "", [], 0
    while i < len(lines) - 1:
        line, nxt = lines[i].strip(), lines[i + 1].strip()
        if line and not line.startswith("..") and _UNDERLINE_RE.match(nxt):
            title = line
            i += 2
            break
        i += 1
    while i < len(lines) and not lines[i].strip():
        i += 1
    while i < len(lines) and lines[i].strip():
        blurb_lines.append(lines[i].strip())
        i += 1
    return _clean_rst(title), _clean_rst(" ".join(blurb_lines))


_DOCSTRING_RE = re.compile(r'^\s*[rRuUbB]{0,2}("""|\'\'\')', re.MULTILINE)


def _script_title(script: Path) -> str:
    """First docstring line — the title sphinx-gallery shows for the card."""
    text = script.read_text(encoding="utf-8", errors="replace")
    match = _DOCSTRING_RE.search(text)
    if not match:
        return script.stem
    for line in text[match.end():].splitlines():
        line = line.strip()
        if line and not _UNDERLINE_RE.match(line):
            return _clean_rst(line)
    return script.stem


# --------------------------------------------------------------------------
# Manifest
# --------------------------------------------------------------------------

def _examples_src(app) -> Path:
    """The gallery *source* tree (docs/examples)."""
    return (Path(app.srcdir) / app.config.sphinx_gallery_conf["examples_dirs"][0]).resolve()


def _section_names(app) -> list[str]:
    """Section directory names in the configured display order."""
    src = _examples_src(app)
    present = sorted(
        d.name for d in src.iterdir()
        if d.is_dir() and (d / "README.txt").exists() and list(d.glob("plot_*.py"))
    )
    order = app.config.sphinx_gallery_conf.get("subsection_order")
    ordered = [
        Path(entry).name
        for entry in getattr(order, "ordered_list", [])
        if Path(entry).name not in ("*", "..")
    ]
    names = [n for n in ordered if n in present]
    names += [n for n in present if n not in names]
    return names


def _manifest(app) -> list[dict]:
    """One entry per section; URLs are relative to examples/index.html."""
    src = _examples_src(app)
    sections = []
    for name in _section_names(app):
        sec_dir = src / name
        title, blurb = _readme_meta(sec_dir / "README.txt")
        thumbs = Path(app.srcdir) / "examples" / name / "images" / "thumb"
        slides = []
        for script in sorted(sec_dir.glob("plot_*.py")):
            thumb = thumbs / f"sphx_glr_{script.stem}_thumb.png"
            if not thumb.exists():
                continue
            slides.append({
                "title": _script_title(script),
                "url": f"{name}/{script.stem}.html",
                "thumb": f"../{THUMB_OUTDIR}/{name}/{thumb.name}",
            })
        sections.append({
            "id": name,
            "title": title or name,
            "blurb": blurb,
            "url": f"{name}/index.html",
            "count": len(list(sec_dir.glob("plot_*.py"))),
            "slides": slides,
        })
    return sections


# --------------------------------------------------------------------------
# Sphinx hooks
# --------------------------------------------------------------------------

def _swap_source(app, docname: str, source: list) -> None:
    """Replace the generated gallery index body, keep its toctree tail."""
    if docname != HUB_DOCNAME:
        return
    match = re.search(r"^\.\. toctree::", source[0], flags=re.MULTILINE)
    tail = source[0][match.start():] if match else ""

    sections = _manifest(app)
    # `<` must not appear inside the inline <script> payload.
    payload = json.dumps(sections, ensure_ascii=False).replace("<", "\\u003c")
    fallback = "".join(
        f'<li><a href="{s["url"]}">{s["title"]} ({s["count"]} examples)</a></li>'
        for s in sections
    )
    hub_html = (
        '   <div id="pycsamt-gallery-hub" class="pgh-grid">\n'
        f"   <noscript><ul>{fallback}</ul></noscript>\n"
        "   </div>\n"
        f"   <script>window.__PYCSAMT_GALLERY_HUB__ = {payload};</script>\n"
    )
    source[0] = f"{HUB_RST}\n.. raw:: html\n\n{hub_html}\n\n{tail}"


def _force_reread(app, env, docnames: list) -> None:
    """Re-read the hub page every build so _swap_source always runs."""
    if HUB_DOCNAME in env.found_docs and HUB_DOCNAME not in docnames:
        docnames.append(HUB_DOCNAME)


def _copy_thumbs(app, exception) -> None:
    """Publish gallery thumbnails under a stable, collision-free path."""
    if exception is not None or app.builder.name != "html":
        return
    for name in _section_names(app):
        thumbs = Path(app.srcdir) / "examples" / name / "images" / "thumb"
        if not thumbs.is_dir():
            continue
        dest = Path(app.outdir) / THUMB_OUTDIR / name
        dest.mkdir(parents=True, exist_ok=True)
        for png in thumbs.glob("*.png"):
            shutil.copy2(png, dest / png.name)


def setup(app):
    app.connect("env-before-read-docs", _force_reread)
    app.connect("source-read", _swap_source)
    app.connect("build-finished", _copy_thumbs)
    return {"parallel_read_safe": True, "parallel_write_safe": True}
