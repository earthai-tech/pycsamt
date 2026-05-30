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
::

    Format:           OCCAM2MTMOD_1.0
    Model Name:       <name>
    Description:      <description>
    Mesh File:        <mesh_filename>
    Mesh Type:        PW2D
    Statics File:     none
    Prejudice File:   none
    Binding Offset:   <float>
    Num Layers:       <N>
      <n_merge>  <n_cols>
      <n_cols integers — parameter codes for each model column>
      ... (repeated N times)
    NO. EXCEPTIONS:   <K>

Parameter codes
---------------
Each integer in a layer's data row encodes the cell type:

* **7** — boundary column (tied to the binding value, not a free parameter)
* **even ≥ 2** — active column; the exact value encodes how many mesh columns
  are merged into one model column (2 = 1:1, 4 = 2:1, …)

The total cell count across all layers equals the ``Param Count`` recorded in
the Startup / .iter file.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional, Union

import numpy as np

from .base   import OccamBase
from .config import OccamConfig

PathLike = Union[str, Path]

__all__ = ["OccamModel"]

_FORMAT_TAG = "OCCAM2MTMOD_1.0"

# Header keyword → attribute name (upper-cased key for matching)
_HEADER_MAP: Dict[str, str] = {
    "FORMAT":          "format_str",
    "MODEL NAME":      "name",
    "DESCRIPTION":     "description",
    "MESH FILE":       "mesh_file",
    "MESH TYPE":       "mesh_type",
    "STATICS FILE":    "statics_file",
    "PREJUDICE FILE":  "prejudice_file",
    "BINDING OFFSET":  "binding_offset",
    "NUM LAYERS":      "n_layers",
    "NO. EXCEPTIONS":  "n_exceptions",
}

_FLOAT_FIELDS = {"binding_offset"}
_INT_FIELDS   = {"n_layers", "n_exceptions"}


# -----------------------------------------------------------------------
# Low-level parser
# -----------------------------------------------------------------------

def _parse_model(path: Path) -> dict:
    """Parse an OCCAM2MTMOD_1.0 file.

    Returns
    -------
    dict
        Keys: all ``_HEADER_MAP`` values plus ``"layers"``
        (``list[dict]`` with keys ``n_merge``, ``n_cols``, ``params``).

    Raises
    ------
    FileNotFoundError
    ValueError
        If the Format tag is absent or wrong.
    """
    if not path.exists():
        raise FileNotFoundError(f"Occam2DModel file not found: {path}")

    result: dict = {v: None for v in _HEADER_MAP.values()}
    result["layers"] = []

    with path.open("r", errors="replace") as fh:
        lines = [l.rstrip("\n") for l in fh]

    i   = 0
    N   = len(lines)

    while i < N:
        raw = lines[i].strip()
        i  += 1

        if not raw:
            continue

        # Header line?
        if ":" in raw:
            raw_key, _, raw_val = raw.partition(":")
            key  = raw_key.strip().upper()
            val  = raw_val.strip()
            attr = _HEADER_MAP.get(key)
            if attr is not None:
                if attr in _INT_FIELDS:
                    try:
                        result[attr] = int(float(val))
                    except ValueError:
                        result[attr] = 0
                elif attr in _FLOAT_FIELDS:
                    try:
                        result[attr] = float(val)
                    except ValueError:
                        result[attr] = 0.0
                elif attr == "format_str":
                    result[attr] = val
                    if val.upper() != _FORMAT_TAG:
                        raise ValueError(
                            f"Expected format '{_FORMAT_TAG}', got '{val}' in {path}"
                        )
                else:
                    result[attr] = val

                # After Num Layers we switch to reading layer blocks
                if attr == "n_layers" and result["n_layers"]:
                    n_layers = result["n_layers"]
                    for _ in range(n_layers):
                        # Control line: n_merge  n_cols
                        while i < N and not lines[i].strip():
                            i += 1
                        ctrl = lines[i].strip().split()
                        i += 1
                        n_merge = int(ctrl[0])
                        n_cols  = int(ctrl[1])
                        # Data line: n_cols integers
                        while i < N and not lines[i].strip():
                            i += 1
                        params = np.array(
                            [int(x) for x in lines[i].strip().split()],
                            dtype=np.int32,
                        )
                        i += 1
                        result["layers"].append(
                            {"n_merge": n_merge, "n_cols": n_cols, "params": params}
                        )

    if result["format_str"] is None:
        raise ValueError(
            f"File does not contain a valid OCCAM2MTMOD_1.0 header: {path}"
        )

    return result


# -----------------------------------------------------------------------
# OccamModel
# -----------------------------------------------------------------------

class OccamModel(OccamBase):
    """Occam2D model-definition container.

    Attributes
    ----------
    format_str : str
    name : str
    description : str
    mesh_file : str
        Filename of the associated mesh (referenced inside the model file).
    mesh_type : str
    statics_file : str
    prejudice_file : str
    binding_offset : float
        Horizontal offset of the binding column.
    n_layers : int
        Number of model layers (depth intervals).
    layers : list[dict]
        Per-layer parameter specification.  Each entry has keys:

        ``n_merge`` : int
            Number of mesh z-rows merged into this model layer.
        ``n_cols`` : int
            Number of model columns in this layer.
        ``params`` : np.ndarray(int), shape (n_cols,)
            Parameter codes for each column (7 = boundary, even ≥ 2 = active).

    n_exceptions : int
    """

    def __init__(
        self,
        name: str = "MODEL MADE BY PYCSAMT",
        description: str = "SMOOTH INVERSION",
        config: Optional[OccamConfig] = None,
        **kwargs,
    ):
        super().__init__(**kwargs)
        self.format_str:     str               = _FORMAT_TAG
        self.name:           str               = name
        self.description:    str               = description
        self.config:         OccamConfig       = config or OccamConfig()
        self.mesh_file:      str               = self.config.mesh_file
        self.mesh_type:      str               = "PW2D"
        self.statics_file:   str               = "none"
        self.prejudice_file: str               = "none"
        self.binding_offset: float             = 0.0
        self.n_layers:       int               = 0
        self.layers:         List[dict]        = []
        self.n_exceptions:   int               = 0

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

        Parameterisation
        ----------------
        Each active (non-air) mesh z-row becomes one model layer
        (``n_merge = 1``).  Horizontal grouping:

        * Left boundary column: code ``7`` (spans 7 mesh x-cells).
        * Interior columns: code ``2`` (one model column per 2 mesh x-cells).
        * Right boundary column: code ``7``.

        This requires ``n_xcells`` to be even so that
        ``(n_xcells − 14) / 2`` is an integer.  The mesh builder
        :meth:`OccamMesh.from_data` guarantees this.

        Parameters
        ----------
        mesh : OccamMesh
            Populated mesh object.
        config : OccamConfig, optional
            Run configuration.

        Returns
        -------
        OccamModel

        Raises
        ------
        ValueError
            If ``n_xcells`` < 14 or the model would have no active layers.
        """
        cfg      = config or OccamConfig()
        obj      = cls(config=cfg, **kwargs)
        n_xcells = mesh.n_xcells
        n_air    = mesh.n_airlayers
        n_active = mesh.n_zcells - n_air

        if n_active <= 0:
            raise ValueError(
                f"OccamModel.from_mesh: no active layers "
                f"(n_zcells={mesh.n_zcells}, n_airlayers={n_air})"
            )
        if n_xcells < 14:
            raise ValueError(
                f"OccamModel.from_mesh: n_xcells={n_xcells} < 14; "
                "mesh too narrow for boundary columns"
            )

        # Build per-layer code array: [7, 2, 2, ..., 2, 7]
        n_interior = (n_xcells - 14) // 2
        remainder  = (n_xcells - 14) - 2 * n_interior
        interior   = [2] * n_interior
        if remainder and n_interior > 0:
            # Absorb remainder into the last interior column (keep it even)
            interior[-1] += remainder
        elif remainder:
            interior = [2 + remainder]
        codes = np.array([7] + interior + [7], dtype=np.int32)

        assert int(codes.sum()) == n_xcells, (
            f"code sum {codes.sum()} != n_xcells {n_xcells}"
        )
        n_cols = len(codes)

        layers: List[dict] = []
        for _ in range(n_active):
            layers.append({
                "n_merge": 1,
                "n_cols":  n_cols,
                "params":  codes.copy(),
            })

        obj.mesh_file    = cfg.mesh_file
        obj.n_layers     = n_active
        obj.layers       = layers
        obj.n_exceptions = 0

        if obj.verbose:
            obj.logger.info(
                "OccamModel.from_mesh: %d layers, %d cols/layer, %d total params",
                n_active, n_cols, obj.n_params,
            )
        return obj

    # ------------------------------------------------------------------
    # I/O
    # ------------------------------------------------------------------
    @classmethod
    def read(cls, path: PathLike, **kwargs) -> "OccamModel":
        """Parse an existing ``Occam2DModel`` file.

        Parameters
        ----------
        path : path-like
            Path to the ``Occam2DModel`` file.

        Returns
        -------
        OccamModel

        Raises
        ------
        FileNotFoundError
            If *path* does not exist.
        ValueError
            If the Format tag is absent or wrong.
        """
        p      = Path(path)
        parsed = _parse_model(p)
        obj    = cls(**kwargs)

        for attr in _HEADER_MAP.values():
            val = parsed.get(attr)
            if val is not None and hasattr(obj, attr):
                setattr(obj, attr, val)

        obj.layers = parsed["layers"]

        if obj.verbose:
            obj.logger.info(
                "OccamModel.read: %d layers, %d total cells from %s",
                obj.n_layers, obj.n_params, p,
            )
        return obj

    def write(self, path: PathLike) -> Path:
        """Write model to *path* in OCCAM2MTMOD_1.0 format.

        Returns the resolved path.
        """
        p = Path(path)
        p.parent.mkdir(parents=True, exist_ok=True)

        _W = 18   # width for the keyword field including colon

        def _kv(kw: str, val) -> str:
            key = f"{kw}:"
            return f"{key:<{_W}}{val}\n"

        lines: list[str] = [
            _kv("Format",          self.format_str),
            _kv("Model Name",      self.name),
            _kv("Description",     self.description),
            _kv("Mesh File",       self.mesh_file),
            _kv("Mesh Type",       self.mesh_type),
            _kv("Statics File",    self.statics_file),
            _kv("Prejudice File",  self.prejudice_file),
            _kv("Binding Offset",  f"{self.binding_offset:.1f}"),
            _kv("Num Layers",      self.n_layers),
        ]

        for layer in self.layers:
            n_merge  = layer["n_merge"]
            n_cols   = layer["n_cols"]
            params   = layer["params"]
            lines.append(f"     {n_merge}   {n_cols}\n")
            lines.append("    " + "    ".join(str(v) for v in params) + "\n")

        lines.append(f"NO. EXCEPTIONS:   {self.n_exceptions}\n")

        with p.open("w") as fh:
            fh.writelines(lines)
        return p

    # ------------------------------------------------------------------
    # Derived
    # ------------------------------------------------------------------
    @property
    def n_params(self) -> int:
        """Total number of model cells (equals ``Param Count`` in the iter file)."""
        return sum(layer["n_cols"] for layer in self.layers)

    @property
    def n_free_params(self) -> int:
        """Number of free (non-boundary) model cells."""
        return sum(
            int((layer["params"] != 7).sum()) for layer in self.layers
        )
