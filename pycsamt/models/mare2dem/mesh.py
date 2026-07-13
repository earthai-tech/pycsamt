# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""MARE2DEM resistivity model and Triangle ``.poly`` mesh access.

Re-exports the full parsers from :mod:`iotools.resistivity` and
:mod:`iotools.poly` and provides thin compatibility wrappers that
preserve the original stub API (``path``, ``header``, ``write``).
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np

from .iotools.poly import PolyFile, read_poly, write_poly
from .iotools.resistivity import (
    ResistivityFile,
    read_resistivity,
    write_resistivity,
)

__all__ = [
    # full objects
    "ResistivityFile",
    "PolyFile",
    # convenience functions
    "read_resistivity",
    "write_resistivity",
    "read_poly",
    "write_poly",
    # compatibility wrappers
    "ResistivityModel",
    "PolyMesh",
]


class ResistivityModel:
    """Thin wrapper around :class:`ResistivityFile` for backwards compatibility.

    Parameters
    ----------
    path : path-like, optional
        Source ``.resistivity`` file.
    """

    def __init__(self, path: str | Path | None = None, **kwargs):
        self._rf: ResistivityFile | None = None
        self.path: Path | None = None
        if path is not None:
            self._rf = read_resistivity(path)
            self.path = Path(path).resolve()

    @property
    def header(self) -> dict[str, Any]:
        if self._rf is None:
            return {}
        return {
            "version": self._rf.version,
            "anisotropy": self._rf.anisotropy,
            "target_misfit": self._rf.target_misfit,
            "num_regions": self._rf.num_regions,
        }

    @property
    def n_elements(self) -> int:
        return self._rf.num_regions if self._rf is not None else 0

    @property
    def n_nodes(self) -> int:
        return 0  # not stored in the .resistivity file itself

    def write(self, path: str | Path) -> Path:
        """Write the model to *path*."""
        if self._rf is None:
            # write a minimal stub
            rf = ResistivityFile(resistivity_file=str(path))
            return write_resistivity(rf, path)
        self._rf.resistivity_file = str(path)
        return write_resistivity(self._rf, path)

    @classmethod
    def halfspace(
        cls,
        log10_rho: float = 0.0,
        *,
        n_nodes: int = 0,
    ) -> ResistivityModel:
        """Return a homogeneous half-space resistivity model stub."""
        m = cls()
        m._rf = ResistivityFile()
        m._rf.resistivity = np.array([[10.0**log10_rho]])
        m._rf.free_parameter = np.zeros((1, 1))
        m._rf.bounds = np.array([[-2.0, 5.0]])
        m._rf.prejudice = np.zeros((1, 2))
        return m

    def __repr__(self) -> str:
        return (
            f"ResistivityModel(path={self.path}, "
            f"n_elements={self.n_elements})"
        )


class PolyMesh:
    """Thin wrapper around :class:`PolyFile` for backwards compatibility.

    Parameters
    ----------
    path : path-like, optional
        Source ``.poly`` file.
    """

    def __init__(self, path: str | Path | None = None, **kwargs):
        self._pf: PolyFile | None = None
        self.path: Path | None = None
        if path is not None:
            self._pf = read_poly(path)
            self.path = Path(path).resolve()

    @property
    def vertices(self) -> np.ndarray | None:
        return self._pf.nodes if self._pf is not None else None

    @property
    def segments(self) -> np.ndarray | None:
        return self._pf.segments if self._pf is not None else None

    @property
    def holes(self) -> np.ndarray | None:
        return self._pf.holes if self._pf is not None else None

    def write(self, path: str | Path) -> Path:
        if self._pf is None:
            pf = PolyFile()
            return write_poly(pf, path)
        return write_poly(self._pf, path)

    def __repr__(self) -> str:
        nv = len(self._pf.nodes) if self._pf is not None else 0
        return f"PolyMesh(path={self.path}, n_vertices={nv})"
