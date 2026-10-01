# -*- mode: python ; coding: utf-8 -*-
"""
PyInstaller spec for the full pycsamt-desktop suite
(:mod:`pycsamt.app.desktop`).

Build with (from anywhere; paths below are all resolved relative to
this file's own directory)::

    pyinstaller packaging/pyinstaller/pycsamt_desktop.spec

or via the wrapper scripts in this directory
(``build_desktop.ps1`` / ``build_desktop.sh``), which also clean old
build artifacts first. Output lands in ``dist/pycsamt-desktop/``
(onedir build -- see ``packaging/pyinstaller/README.md`` for why onedir
over onefile).

Unlike ``pycsamt_converter.spec``, this build is **not** a small,
narrowly-scoped slice of the library. `pycsamt-desktop` reaches:

- ``pycsamt.format`` / ``pycsamt.emtf`` / ``pycsamt.geology`` /
  ``pycsamt.models.{occam2d,modem,mare2dem}`` / ``pycsamt.metadata`` /
  ``pycsamt.gis`` / ``pycsamt.seg`` (same base surface the converter needs),
- ``pycsamt.ai`` (AI Inversion + the 7-tool AI-processing suite) and
  ``pycsamt.agents`` (every registered agent, incl. the 2-D/3-D AI-inversion
  agents and the imputer/uncertainty/distortion/ts-denoise agents added in
  Phase 3) -- these pull in **torch**, **tensorflow**, and **scikit-learn**,
  all via deep, function-local imports (confirmed by grep: e.g.
  ``pycsamt/ai/processing/anomaly.py`` does
  ``from tensorflow.keras import Model, layers`` and
  ``from torch.utils.data import DataLoader, TensorDataset`` *inside*
  methods, never at module scope). Unlike the converter, **this build
  cannot exclude torch/tensorflow** -- doing so would silently break every
  QC/Denoising/Anomaly/Imputer/Distortion/Uncertainty agent and both AI
  Inversion dimensions the moment a user actually clicks one.
- ``pycsamt.airborne`` / ``pycsamt.interp`` / ``pycsamt.map`` /
  ``pycsamt.forward`` (Airborne EM, Interpretation, the PCSF 3-D Studio's
  Plotly figures, and the Maxwell-adapter physics backends),
- ``pycsamt.app.converter`` (embedded, launched in-process by
  ``MainWindow._open_converter_app()`` -- see Phase 9) and
  ``pycsamt.app.agent_master`` (launched as a **subprocess**, not
  in-process -- see ``entry_desktop.py`` and
  ``pycsamt/app/desktop/agent_master_bridge.py`` for why it still needs to
  be part of *this* same frozen build rather than excluded).

``pycsamt.app.mapview`` is included since 2026-09-25: the 3-D window's
"Open in Map View" re-executes this binary with ``--mapview-server``
(``pycsamt/app/desktop/mapview_bridge.py``), like Agent Master.

Confirmed genuinely unreachable from `pycsamt-desktop` and safe to
exclude: ``pycsamt.app.web`` (the Dash web app).

Expect this build to be substantially larger than the converter's ~600 MB
-- torch alone typically adds several hundred MB even CPU-only.
"""

from pathlib import Path

from PyInstaller.utils.hooks import collect_submodules

REPO_ROOT = Path(SPECPATH).resolve().parent.parent
DESKTOP_RESOURCES_DIR = REPO_ROOT / "pycsamt" / "app" / "desktop" / "resources"
DESKTOP_ICONS_DIR = DESKTOP_RESOURCES_DIR / "icons"
CONVERTER_RESOURCES_DIR = REPO_ROOT / "pycsamt" / "app" / "converter" / "resources"
CONVERTER_ICONS_DIR = CONVERTER_RESOURCES_DIR / "icons"
AGENT_MASTER_ASSETS_DIR = REPO_ROOT / "pycsamt" / "app" / "agent_master" / "assets"
MAPVIEW_ASSETS_DIR = REPO_ROOT / "pycsamt" / "app" / "mapview" / "assets"

block_cipher = None
PKG_DIR = REPO_ROOT / "pycsamt"

# Compile output left in the vendored solver sources by local builds
# (objects, modules, DLLs, executables) -- never bundled; the Solver Builder
# rebuilds from the sources themselves.
_BUILD_OUTPUT = {".o", ".mod", ".dll", ".exe", ".a", ".so", ".obj", ".lib",
                 ".pyc"}


def _package_data(rel, patterns=("*",), recursive=False):
    """(file, dest-dir) pairs for pycsamt/<rel> -- the runtime data files
    ``pyproject.toml``'s ``package-data`` ships in the wheel.  A frozen
    build does not read package-data, so these were missing: the EPSG
    table (``gis/epsg.npy``: coordinate conversions failed), the logging
    config, metadata schemas, borehole templates and the Solver Builder's
    scripts and solver sources."""
    root = PKG_DIR / rel
    if root.is_file():
        return [(str(root), str(Path("pycsamt") / Path(rel).parent))]
    out = []
    for pattern in patterns:
        for f in (root.rglob(pattern) if recursive else root.glob(pattern)):
            if (not f.is_file() or f.suffix.lower() in _BUILD_OUTPUT
                    or "__pycache__" in f.parts or "tests" in f.parts):
                continue
            dest = Path("pycsamt") / f.parent.relative_to(PKG_DIR)
            out.append((str(f), str(dest)))
    return out


PACKAGE_DATAS = (
    _package_data("log/p.configlog.yml")
    + _package_data("metadata", ("*.json", "*.yml"))
    + _package_data("format/borehole", ("*.md", "*.json", "*.csv"),
                    recursive=True)
    + _package_data("gis", ("*.npy",))
    + _package_data("models/_solver_build", recursive=True)
    + _package_data("models/occam2d/_source", recursive=True)
    + _package_data("models/modem/_source", recursive=True)
)


def _collect(*packages):
    """collect_submodules over several packages, dropping their test suites
    (never needed at runtime, and several import pytest -- excluded below,
    which would otherwise fail PyInstaller's per-module analysis)."""
    modules = []
    for package in packages:
        modules.extend(collect_submodules(package))
    return [m for m in modules if ".tests" not in m]


hidden_imports = _collect(
    "pycsamt.app.mapview",  # "Open in Map View" (3-D window) subprocess
    "pycsamt.format",
    "pycsamt.emtf",
    "pycsamt.geology",
    "pycsamt.models.occam2d",
    "pycsamt.models.modem",
    "pycsamt.models.mare2dem",
    "pycsamt.metadata",
    "pycsamt.gis",
    "pycsamt.seg",
    "pycsamt.ai",
    "pycsamt.agents",
    "pycsamt.airborne",
    "pycsamt.interp",
    "pycsamt.map",
    "pycsamt.forward",
    "pycsamt.app.desktop",
    "pycsamt.app.converter",
    "pycsamt.app.agent_master",
) + [
    # Reached only via function-local imports throughout pycsamt.ai (see
    # module docstring above) -- PyInstaller's bytecode scan generally
    # catches these too, but the whole point of this spec's approach
    # elsewhere is not depending on that; listed explicitly so the top
    # packages are unambiguously pulled in and their own PyInstaller
    # hooks (bundled with each package / pyinstaller-hooks-contrib) do
    # the rest of the submodule/binary discovery.
    "torch",
    "tensorflow",
    "sklearn",
]

# Confirmed genuinely unreachable (see module docstring) -- excluding them
# keeps the frozen build meaningfully smaller without removing anything a
# desktop-app user can actually click on.
excludes = [
    "pycsamt.app.web",
    "pytest",
    "IPython",
    "jupyter",
    "notebook",
]

a = Analysis(
    [str(REPO_ROOT / "packaging" / "pyinstaller" / "entry_desktop.py")],
    pathex=[str(REPO_ROOT)],
    binaries=[],
    datas=[
        (str(DESKTOP_ICONS_DIR), "pycsamt/app/desktop/resources/icons"),
        (str(DESKTOP_RESOURCES_DIR / "light_theme.qss"), "pycsamt/app/desktop/resources"),
        (str(DESKTOP_RESOURCES_DIR / "dark_theme.qss"), "pycsamt/app/desktop/resources"),
        (str(DESKTOP_RESOURCES_DIR / "pycsamt-v2-splash.png"), "pycsamt/app/desktop/resources"),
        (str(CONVERTER_ICONS_DIR), "pycsamt/app/converter/resources/icons"),
        (str(CONVERTER_RESOURCES_DIR / "light_theme.qss"), "pycsamt/app/converter/resources"),
        (str(CONVERTER_RESOURCES_DIR / "dark_theme.qss"), "pycsamt/app/converter/resources"),
        (str(AGENT_MASTER_ASSETS_DIR), "pycsamt/app/agent_master/assets"),
        (str(MAPVIEW_ASSETS_DIR), "pycsamt/app/mapview/assets"),
    ] + PACKAGE_DATAS,
    hiddenimports=hidden_imports,
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=excludes,
    noarchive=False,
    cipher=block_cipher,
)

# ---------------------------------------------------------------------------
# Work around PyInstaller + PySide6 shared-runtime-DLL version collisions.
# ---------------------------------------------------------------------------
# Identical fixup to pycsamt_converter.spec's own -- see that file's
# comment for the full root-cause explanation (NumPy/SciPy/scikit-learn
# and PySide6/shiboken6 each vendor their own copies of a handful of
# shared MSVC runtime DLLs; whichever copy loads first wins that shared
# DLL name for the whole process, and a non-Qt package's older/different
# copy winning breaks shiboken6/Qt6Core.dll with "the specified procedure
# could not be found" the instant QtWidgets is imported). This build adds
# torch/tensorflow to the mix, which vendor their own copies of the same
# handful of DLL basenames too -- the same generic "force every copy to
# match whichever one PySide6/shiboken6 vendors" fix covers them as well,
# no torch/tensorflow-specific handling needed.
import os as _os
from collections import defaultdict as _defaultdict


def _dest_basename(dest):
    return _os.path.basename(dest.replace("/", "\\")).lower()


def _dest_top_dir(dest):
    parts = dest.replace("/", "\\").split("\\")
    return parts[0].lower() if len(parts) > 1 else ""


_by_basename = _defaultdict(list)
for _dest, _src, _kind in a.binaries:
    if _dest.lower().endswith(".dll"):
        _by_basename[_dest_basename(_dest)].append((_dest, _src, _kind))

_canonical_src = {}
for _basename, _entries in _by_basename.items():
    if len(_entries) < 2:
        continue
    _distinct_srcs = {_src for _, _src, _ in _entries}
    if len(_distinct_srcs) < 2:
        continue  # same bytes everywhere already, nothing to fix
    _preferred = next(
        (
            (_dest, _src)
            for _dest, _src, _ in _entries
            if _dest_top_dir(_dest) in ("pyside6", "shiboken6")
        ),
        None,
    )
    if _preferred is not None:
        _canonical_src[_basename] = _preferred[1]

if _canonical_src:
    a.binaries = [
        (
            (dest, _canonical_src[_dest_basename(dest)], kind)
            if _dest_basename(dest) in _canonical_src
            else (dest, src, kind)
        )
        for dest, src, kind in a.binaries
    ]

# ---------------------------------------------------------------------------
# Drop a stray, ancient ICU pulled in from an unrelated Anaconda *base* env.
# ---------------------------------------------------------------------------
# Same fixup as pycsamt_converter.spec -- see that file's comment. Rebuild
# on a clean (non-Anaconda-base-polluted) machine may make this a no-op;
# harmless either way since it only drops entries that match these exact
# stray basenames.
_STRAY_ICU_DLLS = {"icuuc.dll", "icudt58.dll", "icuin58.dll", "icuuc58.dll"}
a.binaries = [
    (dest, src, kind)
    for dest, src, kind in a.binaries
    if _dest_basename(dest) not in _STRAY_ICU_DLLS
]

pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name="pycsamt-desktop",
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=False,
    console=False,
    # pycsamt.logo.ico, not pycsamt.ico -- matches the running app's own
    # window/taskbar icon (main_window.py, __main__.py via
    # branding.LOGO_ICO). Note: this file currently embeds only one
    # 32x31 resolution, unlike pycsamt.ico's full 16/32/48/64/128/256
    # set, so the exe's icon will look softer in Explorer's "large
    # icons" view than pycsamt.ico did.
    icon=str(DESKTOP_ICONS_DIR / "pycsamt.logo.ico"),
)

coll = COLLECT(
    exe,
    a.binaries,
    a.zipfiles,
    a.datas,
    strip=False,
    upx=False,
    name="pycsamt-desktop",
)
