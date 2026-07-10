"""Shared helpers for the ``pycsamt.pipeline`` example gallery.

Not an example itself (sphinx-gallery skips ``_``-prefixed files). Provides
the sample survey line every pipeline example runs on, a scratch output
directory, and a context manager that silences the benign ``rho/phi``
recompute ERROR logs the frequency steps emit. The helper name is unique to
avoid the shared-``sys.modules`` collision that same-named gallery helpers
would otherwise hit.
"""

from __future__ import annotations

import contextlib
import logging
import os
import tempfile
from pathlib import Path

from pycsamt.emtools._core import ensure_sites


def _repo_root() -> Path:
    root = os.environ.get("PYCSAMT_DOCS_REPO_ROOT")
    return Path(root) if root else Path.cwd()


def demo_sites(n: int = 8):
    """Load the first *n* stations of the WILLY_DATA L22 line as ``Sites``."""
    data_dir = _repo_root() / "data" / "AMT" / "WILLY_DATA" / "L22PLT"
    paths = sorted(data_dir.glob("*.edi"))[:n]
    return ensure_sites([str(p) for p in paths])


def scratch_dir(prefix: str = "pycsamt_pipe_") -> Path:
    """A throwaway output directory (never inside the repo)."""
    return Path(tempfile.mkdtemp(prefix=prefix))


def station_rho(sites, component: str = "xy"):
    """Return ``{station: (period_s, rho_a)}`` for one impedance component.

    Used by the advanced examples to draw before/after and multi-strategy
    apparent-resistivity comparisons directly from a ``Sites`` object.
    """
    import numpy as np

    from pycsamt.emtools._core import _iter_items, _name

    ij = {"xx": (0, 0), "xy": (0, 1), "yx": (1, 0), "yy": (1, 1)}[component]
    out = {}
    for i, ed in enumerate(_iter_items(sites)):
        freq = np.asarray(ed.freq, dtype=float)
        rho = np.asarray(ed.rho, dtype=float)[:, ij[0], ij[1]]
        out[_name(ed, i)] = (1.0 / freq, rho)
    return out


@contextlib.contextmanager
def quiet_logs():
    """Silence the ``Z`` logger's benign 'Failed to compute rho/phi' ERRORs.

    The frequency steps re-set Z with a differently-shaped error tensor,
    which the impedance layer logs (harmlessly) at ERROR. This keeps the
    gallery's captured output focused on the pipeline itself.
    """
    logging.disable(logging.ERROR)
    try:
        yield
    finally:
        logging.disable(logging.NOTSET)
