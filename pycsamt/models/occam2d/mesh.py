# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""OccamMesh — build and parse the Occam2DMesh file.

The mesh defines the 2-D finite-element grid used by the forward solver.
It is built from station positions and depth parameterisation specified in
``OccamConfig``.

Entry points
------------
``OccamMesh.from_data(occam_data, config)``
    Construct a mesh consistent with station positions in a data file.
``OccamMesh.read(path)``
    Parse an existing ``Occam2DMesh`` file.
``OccamMesh.write(path)``
    Serialise to disk in the native PW2D mesh format.

Mesh file format (PW2D)
-----------------------
Two header integers: n_xcells  n_zcells
Followed by ``n_zcells`` rows, each listing cell widths / thicknesses.
Column lines encode: [layer_flag  n_cols  col_widths...]
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional, Union

import numpy as np

from .base   import OccamBase
from .config import OccamConfig

PathLike = Union[str, Path]

__all__ = ["OccamMesh"]


class OccamMesh(OccamBase):
    """Occam2D mesh container.

    Attributes
    ----------
    x_nodes : np.ndarray, shape (n_xcells+1,)
        Horizontal node positions in metres (profile direction).
    z_nodes : np.ndarray, shape (n_zcells+1,)
        Depth node positions in metres (positive downward).
    cell_widths : np.ndarray, shape (n_zcells, n_xcells)
        Width of each cell per row (PW2D representation).
    n_airlayers : int
        Number of rows treated as air layers.
    """

    def __init__(
        self,
        config: Optional[OccamConfig] = None,
        **kwargs,
    ):
        super().__init__(**kwargs)
        self.config      = config or OccamConfig()
        self.x_nodes:      np.ndarray = np.array([])
        self.z_nodes:      np.ndarray = np.array([])
        self.cell_widths:  np.ndarray = np.empty((0, 0))
        self.n_airlayers:  int        = self.config.n_airlayers

    # ------------------------------------------------------------------
    # Construction
    # ------------------------------------------------------------------
    @classmethod
    def from_data(
        cls,
        data: "OccamData",  # noqa: F821
        config: Optional[OccamConfig] = None,
        **kwargs,
    ) -> "OccamMesh":
        """Build a mesh from station positions encoded in *data*.

        Parameters
        ----------
        data : OccamData
            Populated data object (station offsets required).
        config : OccamConfig, optional
            Mesh parameters (n_layers, cell sizes, padding, …).
        """
        # TODO: implement
        # 1. Derive x-cell widths from station offsets + padding
        # 2. Build z-cell thicknesses with geometric growth (depth_scale)
        # 3. Prepend n_airlayers air rows
        raise NotImplementedError("OccamMesh.from_data — not yet implemented")

    # ------------------------------------------------------------------
    # I/O
    # ------------------------------------------------------------------
    @classmethod
    def read(cls, path: PathLike, **kwargs) -> "OccamMesh":
        """Parse an existing ``Occam2DMesh`` file."""
        # TODO: implement PW2D format parser
        raise NotImplementedError("OccamMesh.read — not yet implemented")

    def write(self, path: PathLike) -> Path:
        """Write mesh to *path* in PW2D format."""
        # TODO: implement
        raise NotImplementedError("OccamMesh.write — not yet implemented")

    # ------------------------------------------------------------------
    # Convenience
    # ------------------------------------------------------------------
    @property
    def n_xcells(self) -> int:
        return max(0, len(self.x_nodes) - 1)

    @property
    def n_zcells(self) -> int:
        return max(0, len(self.z_nodes) - 1)

    @property
    def n_params(self) -> int:
        """Number of free model parameters (non-air, non-padding cells)."""
        # TODO: compute after mesh is built
        return 0
