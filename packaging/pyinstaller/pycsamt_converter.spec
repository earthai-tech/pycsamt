# -*- mode: python ; coding: utf-8 -*-
"""
PyInstaller spec for pyCSAMT Format Studio (the standalone format
converter app, :mod:`pycsamt.app.converter`).

Build with (from anywhere; paths below are all resolved relative to
this file's own directory)::

    pyinstaller packaging/pyinstaller/pycsamt_converter.spec

or via the wrapper scripts in this directory
(``build_converter.ps1`` / ``build_converter.sh``), which also clean
old build artifacts first. Output lands in ``dist/pycsamt-converter/``
(onedir build -- see the module docstring in
``pycsamt/app/converter/__init__.py`` for why this app is a separate,
freezable build from the full ``pycsamt-desktop`` suite).
"""

from pathlib import Path

from PyInstaller.utils.hooks import collect_submodules

REPO_ROOT = Path(SPECPATH).resolve().parent.parent
RESOURCES_DIR = REPO_ROOT / "pycsamt" / "app" / "converter" / "resources"
ICONS_DIR = RESOURCES_DIR / "icons"

block_cipher = None

# The converter app only reaches pycsamt.format / pycsamt.emtf / a few
# pycsamt.models.*.results modules (see pycsamt/app/converter/jobs.py) --
# but several of those do lazy, function-local imports, and
# `pycsamt.format`'s own __init__ eagerly pulls in its whole namespace
# (adapters, borehole, geology, structure, pointset, ...). Rather than
# hand-list every leaf module, collect entire submodule trees for the
# packages the app actually touches -- cheap (these are the packages
# already required for a correct build) and safe against a single
# missed function-local import breaking the frozen app at runtime.
def _collect(*packages):
    """collect_submodules over several packages, dropping their test suites
    (never needed at runtime, and a couple import pytest -- excluded below,
    which would otherwise fail PyInstaller's per-module analysis)."""
    modules = []
    for package in packages:
        modules.extend(collect_submodules(package))
    return [m for m in modules if ".tests" not in m]


hidden_imports = _collect(
    "pycsamt.format",
    "pycsamt.emtf",
    "pycsamt.geology",
    "pycsamt.models.occam2d",
    "pycsamt.models.modem",
    "pycsamt.models.mare2dem",
    "pycsamt.metadata",
    "pycsamt.gis",
    "pycsamt.seg",
)

# Packages the converter app's import graph never reaches (verified by
# the app's headless test/click-through suite importing only
# pycsamt.format / pycsamt.emtf / pycsamt.models.*.results): the AI/DL
# stack, the other pyCSAMT desktop/web/mapview/agent apps, and the test
# runner. Excluding them keeps the frozen build a few hundred MB
# smaller. If a future converter feature needs one of these, drop it
# from this list rather than fighting a "module not found" at runtime.
excludes = [
    "torch",
    "tensorflow",
    "pycsamt.ai",
    "pycsamt.app.desktop",
    "pycsamt.app.web",
    "pycsamt.app.mapview",
    "pycsamt.app.agent_master",
    "pytest",
    "IPython",
    "jupyter",
    "notebook",
]

a = Analysis(
    [str(REPO_ROOT / "packaging" / "pyinstaller" / "entry_converter.py")],
    pathex=[str(REPO_ROOT)],
    binaries=[],
    datas=[
        (str(ICONS_DIR), "pycsamt/app/converter/resources/icons"),
        (str(RESOURCES_DIR / "light_theme.qss"), "pycsamt/app/converter/resources"),
        (str(RESOURCES_DIR / "dark_theme.qss"), "pycsamt/app/converter/resources"),
    ],
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
# NumPy/SciPy/scikit-learn and PySide6/shiboken6 each vendor their own
# copies of a handful of shared MSVC runtime DLLs (MSVCP140.dll,
# VCRUNTIME140.dll, VCRUNTIME140_1.dll, VCOMP140.dll, ...). PyInstaller
# collects every package's copy into the onedir build side by side
# (e.g. root-level, PySide6/, shiboken6/, sklearn/.libs/); whichever
# gets loaded into the process first wins that shared DLL *name* for
# the whole process. If a non-Qt package's (older/differently built)
# copy wins, shiboken6's own native module -- which every PySide6
# submodule depends on -- or Qt6Core.dll itself silently resolves its
# runtime import to the wrong version and fails with "the specified
# procedure could not be found" the instant QtWidgets is imported
# (reproduced directly with ctypes.WinDLL against a frozen build while
# writing this spec: MSVCP140/VCRUNTIME140/VCRUNTIME140_1 differed
# across 4 locations, VCOMP140 across 2).
#
# Rather than hand-list the exact filenames (fragile -- the next
# PySide6/NumPy/scikit-learn release could vendor a new one), find
# every DLL basename that appears more than once in the bundle *and*
# has a copy under PySide6/ or shiboken6/, and force every other copy
# of that basename to be byte-identical to that one -- these are
# strict-superset-compatible MSVC redistributable runtimes, so a newer
# Qt-vendored copy is safe for NumPy/SciPy/scikit-learn too.
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
# PyInstaller's binary-dependency walk found Qt6Core.dll's `icuuc.dll`
# import by searching PATH rather than PySide6's own package directory
# (PySide6 6.11 doesn't vendor one -- it ships `icudtl.dat` instead, and
# runs fine on a machine with no icuuc.dll at all), and picked up
# %CONDA_BASE%/Library/bin/icuuc.dll + icudt58.dll: ICU 58, an
# ~decade-old build from an entirely different, unrelated environment.
# Its symbol layout doesn't match what this Qt6Core.dll expects, so
# bundling it turns a normally-harmless situation (import present, but
# never actually resolved against a wrong version at runtime) into a
# hard "the specified procedure could not be found" failure the instant
# QtWidgets is imported. Exclude it; Qt6Core.dll loads fine without it.
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
    name="pycsamt-converter",
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=False,
    console=False,
    icon=str(ICONS_DIR / "pycsamt.ico"),
)

coll = COLLECT(
    exe,
    a.binaries,
    a.zipfiles,
    a.datas,
    strip=False,
    upx=False,
    name="pycsamt-converter",
)
