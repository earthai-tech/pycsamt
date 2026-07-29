"""Shared bundled-sample-data loader for the pyCSAMT example galleries.

Not an example itself: sphinx-gallery is configured (see
``sphinx_gallery_conf["ignore_pattern"]`` in ``docs/source/conf.py``) to
skip files whose name starts with ``_``. sphinx-gallery executes every
``plot_*.py`` script from its own source directory, so any script can do
``from _datasets import load_survey`` / ``load_sites``.

.. important::

   sphinx-gallery runs *all* gallery sections in one Python process and
   caches imported modules in ``sys.modules`` by name. Because more than
   one section ships a ``_datasets.py`` helper, whichever loads first wins
   for the whole build. This file is therefore a **superset** kept
   byte-identical with ``docs/examples/survey_diagnostics/_datasets.py`` —
   it exposes every name any section needs (``load_survey``, ``load_sites``,
   ``dataset_path``, ``line_groups``) over the union of all bundled
   datasets, so the collision is harmless. Edit both copies together.
"""

from __future__ import annotations

import os
from pathlib import Path

from pycsamt.api import read_edis
from pycsamt.emtools._core import (
    _iter_items,
    _name,
    ensure_sites,
)

# name -> path segments under <repo_root>/data (union of every section)
_DATASETS = {
    "amt_willy": ("AMT", "WILLY_DATA"),  # all 5 lines, 128 stations
    "amt_l18plt": ("AMT", "WILLY_DATA", "L18PLT"),
    "amt_l22plt": ("AMT", "WILLY_DATA", "L22PLT"),
    "amt_l26plt": ("AMT", "WILLY_DATA", "L26PLT"),
    "amt_l30plt": ("AMT", "WILLY_DATA", "L30PLT"),
    "amt_l34plt": ("AMT", "WILLY_DATA", "L34PLT"),
    "mt_kap03": ("MT", "kap03lmt_edis"),  # SAMTEX line, real tipper
    "mt_spectra": ("MT", "SPECTRA"),
}


def repo_root() -> Path:
    """Repo root, set by ``docs/source/conf.py`` during a docs build."""
    root = os.environ.get("PYCSAMT_DOCS_REPO_ROOT")
    return Path(root) if root else Path.cwd()


def dataset_path(name: str) -> Path:
    """Resolve one of the bundled dataset names to its directory."""
    try:
        parts = _DATASETS[name]
    except KeyError:
        raise KeyError(
            f"Unknown bundled dataset {name!r}; choose from "
            f"{sorted(_DATASETS)}"
        ) from None
    return repo_root() / "data" / Path(*parts)


def load_survey(name: str, **kwargs):
    """Read a bundled dataset into an ``APISurvey`` via :func:`read_edis`."""
    kwargs.setdefault("recursive", False)
    kwargs.setdefault("progress", False)
    return read_edis(dataset_path(name), **kwargs)


def load_sites(name: str, recursive: bool = False, **kwargs):
    """Read a bundled dataset into a ready-to-plot ``Sites`` object."""
    return ensure_sites(str(dataset_path(name)), recursive=recursive, **kwargs)


def line_groups(sites):
    """Map each WILLY_DATA profile label to its sorted station names."""
    names = [_name(ed, i) for i, ed in enumerate(_iter_items(sites))]
    return {
        "L18": sorted(n for n in names if n.startswith("18-")),
        "L22": sorted(n for n in names if n.startswith("22-")),
        "L26": sorted(n for n in names if n.startswith("26-")),
        "L30": sorted(n for n in names if n.startswith("30-")),
        "L34": sorted(n for n in names if n.startswith("34-")),
    }
