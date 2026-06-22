# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""MARE2DEM ``.emdata`` / ``.EMResp`` file access.

This module re-exports the full parser from
:mod:`pycsamt.models.mare2dem.iotools.emdata` and adds a thin
:class:`EMData` compatibility wrapper that matches the original stub
API (``path``, ``header``, ``data``, ``write``).
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np

from .iotools.emdata import (
    CSEMConfig,
    DCConfig,
    EMDataFile,
    MTConfig,
    UTMOrigin,
    read_emdata,
    write_emdata,
)

__all__ = [
    # full objects
    "EMDataFile",
    "UTMOrigin",
    "CSEMConfig",
    "MTConfig",
    "DCConfig",
    # convenience functions
    "read_emdata",
    "write_emdata",
    # compatibility wrapper
    "EMData",
]


class EMData:
    """Thin wrapper around :class:`EMDataFile` for backwards compatibility.

    Provides the same ``path``, ``header``, ``data``, and ``write()``
    interface as the original stub while delegating to the full parser.

    Parameters
    ----------
    path : path-like, optional
        Source file. Read immediately when provided.
    """

    def __init__(self, path: str | Path | None = None, **kwargs):
        self._em: EMDataFile | None = None
        self.path: Path | None = None
        if path is not None:
            self._em = read_emdata(path)
            self.path = self._em.path

    # ------------------------------------------------------------------
    # Legacy attribute access
    # ------------------------------------------------------------------

    @property
    def header(self) -> dict[str, Any]:
        if self._em is None:
            return {}
        return {
            "Format": self._em.format,
            "Comment": self._em.comment,
            "UTM": self._em.utm.__dict__ if self._em.utm else {},
        }

    @property
    def data(self) -> np.ndarray | None:
        return self._em.data if self._em is not None else None

    @property
    def n_data(self) -> int:
        return self._em.n_data if self._em is not None else 0

    def write(self, path: str | Path) -> Path:
        """Write the stored ``.emdata`` to *path*."""
        if self._em is None:
            raise RuntimeError("No EMData loaded — nothing to write.")
        return write_emdata(self._em, path)

    def __repr__(self) -> str:
        return f"EMData(path={self.path}, n_data={self.n_data})"
