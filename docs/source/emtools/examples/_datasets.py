"""Shared bundled-sample-data loader for ``pycsamt.emtools`` gallery scripts.

Not an example itself: sphinx-gallery is configured (see
``sphinx_gallery_conf["ignore_pattern"]`` in ``docs/source/conf.py``) to
skip this file. sphinx-gallery executes every ``plot_*.py`` script from
its own source directory, so any other script in this folder can just
``from _datasets import load_survey``.
"""

from __future__ import annotations

import os
from pathlib import Path

from pycsamt.api import read_edis

# name -> path segments under <repo_root>/data
_DATASETS = {
    "amt_l18plt": ("AMT", "WILLY_DATA", "L18PLT"),
    "amt_l22plt": ("AMT", "WILLY_DATA", "L22PLT"),
    "mt_kap03": ("MT", "kap03lmt_edis"),
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
    """Read a bundled sample line into an ``APISurvey``.

    Parameters
    ----------
    name : {"amt_l18plt", "amt_l22plt", "mt_kap03"}
        See the module docstring / ``data/AMT/WILLY_DATA/README.md`` and
        ``data/MT/README.md`` for what each dataset is.
    **kwargs
        Forwarded to :func:`pycsamt.api.read_edis` (``recursive=False``
        and ``progress=False`` are applied by default).
    """
    kwargs.setdefault("recursive", False)
    kwargs.setdefault("progress", False)
    return read_edis(dataset_path(name), **kwargs)
