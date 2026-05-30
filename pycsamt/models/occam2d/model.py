# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""OccamModel — read and write the Occam2DModel file.

The model file links the mesh to the inversion parameterisation.  It
specifies which mesh cells share the same free parameter (the
``Num Layers`` / column-mapping block) and global model metadata.

Entry points
------------
``OccamModel.from_mesh(mesh, config)``
    Build a model definition consistent with a given mesh.
``OccamModel.read(path)``
    Parse an existing ``Occam2DModel`` file.
``OccamModel.write(path)``
    Serialise to disk in OCCAM2MTMOD_1.0 format.

File format (OCCAM2MTMOD_1.0)
------------------------------
Format:         OCCAM2MTMOD_1.0
Model Name:     <name>
Description:    <description>
Mesh File:      <mesh_filename>
Mesh Type:      PW2D
Statics File:   none
Prejudice File: none
Binding Offset: <float>
Num Layers:     <n>
  <n_params_row>  <n_cols_row>
  <param_codes…>
  ...
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional, Union

import numpy as np

from .base   import OccamBase
from .config import OccamConfig

PathLike = Union[str, Path]

__all__ = ["OccamModel"]


class OccamModel(OccamBase):
    """Occam2D model-definition container.

    Attributes
    ----------
    name : str
    description : str
    mesh_file : str
        Relative path to the mesh file (referenced inside the model file).
    binding_offset : float
        Horizontal offset of the binding column.
    n_layers : int
        Number of free-parameter layers.
    param_map : np.ndarray
        Integer array encoding which parameter index each cell belongs to.
    """

    def __init__(
        self,
        name: str = "MODEL MADE BY PYCSAMT",
        description: str = "SMOOTH INVERSION",
        config: Optional[OccamConfig] = None,
        **kwargs,
    ):
        super().__init__(**kwargs)
        self.name            = name
        self.description     = description
        self.config          = config or OccamConfig()
        self.mesh_file       = self.config.mesh_file
        self.binding_offset  = 0.0
        self.n_layers        = 0
        self.param_map:  np.ndarray = np.empty((0,), dtype=int)

    # ------------------------------------------------------------------
    # Construction
    # ------------------------------------------------------------------
    @classmethod
    def from_mesh(
        cls,
        mesh: "OccamMesh",  # noqa: F821
        config: Optional[OccamConfig] = None,
        **kwargs,
    ) -> "OccamModel":
        """Build a model definition from a populated mesh.

        Each non-air column of cells is assigned a unique free parameter
        per depth layer, giving a standard smooth-model parameterisation.
        """
        # TODO: implement
        # 1. Iterate over mesh layers (excluding air layers)
        # 2. Assign contiguous parameter indices
        # 3. Set n_layers, param_map
        raise NotImplementedError("OccamModel.from_mesh — not yet implemented")

    # ------------------------------------------------------------------
    # I/O
    # ------------------------------------------------------------------
    @classmethod
    def read(cls, path: PathLike, **kwargs) -> "OccamModel":
        """Parse an existing ``Occam2DModel`` file."""
        # TODO: implement
        raise NotImplementedError("OccamModel.read — not yet implemented")

    def write(self, path: PathLike) -> Path:
        """Write model to *path* in OCCAM2MTMOD_1.0 format."""
        # TODO: implement
        raise NotImplementedError("OccamModel.write — not yet implemented")

    @property
    def n_params(self) -> int:
        """Total number of free inversion parameters."""
        return int(self.param_map.max()) if self.param_map.size else 0
