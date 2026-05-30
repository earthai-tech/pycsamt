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
::

    <comment line>
    0  num_h_nodes  num_v_nodes  n_airlayers  0  mesh_type
    <x-cell widths: num_h_nodes-1 floats, wrapped across any number of lines>
    <blank line — optional>
    <z-layer thicknesses: num_v_nodes-1 floats, wrapped across any number of lines>
    0
    <cell-type character rows: num_v_nodes rows of num_h_nodes-1 characters each>
        '?' = free model parameter, other = fixed / air / boundary
"""

from __future__ import annotations

from pathlib import Path
from typing import List, Optional, Union

import numpy as np

from .base   import OccamBase
from .config import OccamConfig

PathLike = Union[str, Path]

__all__ = ["OccamMesh"]


# -----------------------------------------------------------------------
# Low-level parser
# -----------------------------------------------------------------------

def _parse_mesh(path: Path) -> dict:
    """Parse an Occam2DMesh PW2D file.

    Returns
    -------
    dict with keys:
        comment, n_airlayers, x_widths, z_widths,
        x_nodes, z_nodes, cell_rows (list[str])

    Raises
    ------
    FileNotFoundError
        If *path* does not exist.
    ValueError
        If the control line cannot be parsed.
    """
    if not path.exists():
        raise FileNotFoundError(f"Occam2DMesh file not found: {path}")

    with path.open("r", errors="replace") as fh:
        lines = fh.readlines()

    if len(lines) < 2:
        raise ValueError(f"Mesh file too short: {path}")

    comment = lines[0].strip()

    # Line 1: control integers — "0  num_h_nodes  num_v_nodes  n_airlayers  0  mesh_type"
    ctrl = lines[1].strip().split()
    if len(ctrl) < 3:
        raise ValueError(f"Cannot parse mesh control line: {lines[1]!r}")
    num_h_nodes  = int(ctrl[1])
    num_v_nodes  = int(ctrl[2])
    n_airlayers  = int(ctrl[3]) if len(ctrl) > 3 else 0
    n_xcells     = num_h_nodes - 1
    n_zcells     = num_v_nodes - 1

    # Collect float tokens until the binding sentinel ('0' alone) then
    # collect char-matrix rows for the rest of the file.
    float_tokens: list[float] = []
    cell_rows:    list[str]   = []
    in_char = False

    for raw in lines[2:]:
        stripped = raw.strip()
        if not stripped:
            continue

        if in_char:
            cell_rows.append(stripped)
            continue

        # Binding sentinel: a line that contains only '0' and we have
        # already read all expected float tokens.
        if stripped == "0" and len(float_tokens) >= n_xcells + n_zcells:
            in_char = True
            continue

        for tok in stripped.split():
            try:
                float_tokens.append(float(tok))
            except ValueError:
                pass

    x_widths = np.array(float_tokens[:n_xcells],            dtype=float)
    z_widths = np.array(float_tokens[n_xcells: n_xcells + n_zcells], dtype=float)

    x_nodes  = np.concatenate([[0.0], np.cumsum(x_widths)])
    z_nodes  = np.concatenate([[0.0], np.cumsum(z_widths)])

    return {
        "comment":     comment,
        "n_airlayers": n_airlayers,
        "x_widths":    x_widths,
        "z_widths":    z_widths,
        "x_nodes":     x_nodes,
        "z_nodes":     z_nodes,
        "cell_rows":   cell_rows,
    }


# -----------------------------------------------------------------------
# OccamMesh
# -----------------------------------------------------------------------

class OccamMesh(OccamBase):
    """Occam2D mesh container.

    Attributes
    ----------
    comment : str
        First line of the mesh file (provenance comment).
    x_widths : np.ndarray, shape (n_xcells,)
        Horizontal cell widths in metres.
    z_widths : np.ndarray, shape (n_zcells,)
        Vertical layer thicknesses in metres.
    x_nodes : np.ndarray, shape (n_xcells+1,)
        Cumulative horizontal node positions (leading 0).
    z_nodes : np.ndarray, shape (n_zcells+1,)
        Cumulative depth node positions (leading 0, positive down).
    cell_rows : list[str]
        Raw character rows from the mesh file; each character encodes the
        cell type at that column position (``'?'`` = free parameter).
    n_airlayers : int
        Number of rows treated as air layers.
    """

    def __init__(
        self,
        config: Optional[OccamConfig] = None,
        **kwargs,
    ):
        super().__init__(**kwargs)
        self.config:       OccamConfig = config or OccamConfig()
        self.comment:      str         = ""
        self.x_widths:     np.ndarray  = np.array([])
        self.z_widths:     np.ndarray  = np.array([])
        self.x_nodes:      np.ndarray  = np.array([])
        self.z_nodes:      np.ndarray  = np.array([])
        self.cell_rows:    List[str]   = []
        self.n_airlayers:  int         = self.config.n_airlayers

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

        The mesh uses the standard Occam2D PW2D layout:

        * 7 horizontal padding cells on each side (geometric expansion),
          matching the boundary-column code ``7`` used by
          :meth:`OccamModel.from_mesh`.
        * Station-zone cells with width ``≈ config.cell_size_horizontal``.
        * ``config.n_airlayers`` thin air rows at the top.
        * ``config.n_layers`` active rows with geometric depth growth.

        Parameters
        ----------
        data : OccamData
            Populated data object (``data.offsets`` required).
        config : OccamConfig, optional
            Mesh parameters (n_layers, cell sizes, padding, …).

        Returns
        -------
        OccamMesh
        """
        cfg = config or OccamConfig()
        obj = cls(config=cfg, **kwargs)

        offsets = np.sort(np.asarray(data.offsets, dtype=float))
        if offsets.size == 0:
            raise ValueError("OccamMesh.from_data: data has no station offsets")

        cell_h  = cfg.cell_size_horizontal
        n_pad   = 7   # fixed: matches boundary code 7 in OccamModel.from_mesh

        # ---- station-zone x-cells ----------------------------------------
        station_widths: list[float] = []
        if offsets.size == 1:
            station_widths = [cell_h]
        else:
            for i in range(len(offsets) - 1):
                gap    = float(offsets[i + 1] - offsets[i])
                n_cell = max(1, round(gap / cell_h))
                w      = gap / n_cell
                station_widths.extend([w] * n_cell)

        # Ensure even count so that (n_xcells - 14) is even for code-2 interior cols
        if len(station_widths) % 2 != 0:
            station_widths.append(station_widths[-1] if station_widths else cell_h)

        # ---- horizontal padding (7 cells each side, geometrically expanding)
        pad = [cell_h * float(2 ** (k + 1)) for k in range(n_pad)]
        left_pad  = list(reversed(pad))   # widest cell at outer edge
        right_pad = pad

        x_widths = np.array(left_pad + station_widths + right_pad, dtype=float)

        # ---- vertical cells ----------------------------------------------
        n_air    = cfg.n_airlayers
        n_active = cfg.n_layers
        air_h    = cfg.cell_size_vertical_top

        air_w = [air_h] * n_air
        z_w: list[float] = []
        thick = float(cfg.cell_size_vertical_top)
        for _ in range(n_active):
            z_w.append(thick)
            thick *= cfg.depth_scale

        z_widths = np.array(air_w + z_w, dtype=float)

        # ---- char matrix (4 rows per z-cell) ----------------------------
        n_xcells = len(x_widths)
        n_zcells = len(z_widths)
        cell_rows: list[str] = []
        for iz in range(n_zcells):
            row = ("0" if iz < n_air else "?") * n_xcells
            for _ in range(4):
                cell_rows.append(row)

        # ---- populate object --------------------------------------------
        obj.comment     = "MESH FILE Created by pycsamt.models.occam2d"
        obj.x_widths    = x_widths
        obj.z_widths    = z_widths
        obj.x_nodes     = np.concatenate([[0.0], np.cumsum(x_widths)])
        obj.z_nodes     = np.concatenate([[0.0], np.cumsum(z_widths)])
        obj.cell_rows   = cell_rows
        obj.n_airlayers = n_air

        if obj.verbose:
            obj.logger.info(
                "OccamMesh.from_data: %d×%d cells (%d airlayers) from %d stations",
                n_xcells, n_zcells, n_air, offsets.size,
            )
        return obj

    # ------------------------------------------------------------------
    # I/O
    # ------------------------------------------------------------------
    @classmethod
    def read(cls, path: PathLike, **kwargs) -> "OccamMesh":
        """Parse an existing ``Occam2DMesh`` PW2D file.

        Parameters
        ----------
        path : path-like
            Path to the ``Occam2DMesh`` file.

        Returns
        -------
        OccamMesh

        Raises
        ------
        FileNotFoundError
            If *path* does not exist.
        ValueError
            If the control line cannot be parsed.
        """
        p   = Path(path)
        d   = _parse_mesh(p)
        obj = cls(**kwargs)

        obj.comment     = d["comment"]
        obj.n_airlayers = d["n_airlayers"]
        obj.x_widths    = d["x_widths"]
        obj.z_widths    = d["z_widths"]
        obj.x_nodes     = d["x_nodes"]
        obj.z_nodes     = d["z_nodes"]
        obj.cell_rows   = d["cell_rows"]

        if obj.verbose:
            obj.logger.info(
                "OccamMesh.read: %d×%d cells (%d char rows) from %s",
                obj.n_xcells, obj.n_zcells, len(obj.cell_rows), p,
            )
        return obj

    def write(self, path: PathLike) -> Path:
        """Write mesh to *path* in PW2D format.

        Returns the resolved path.
        """
        p = Path(path)
        p.parent.mkdir(parents=True, exist_ok=True)

        n_h = self.n_xcells + 1   # num_h_nodes
        n_v = self.n_zcells + 1   # num_v_nodes
        comment = self.comment or "MESH FILE Created by pycsamt"

        lines: list[str] = []
        lines.append(f"{comment}\n")
        lines.append(f"   0  {n_h}  {n_v}  {self.n_airlayers}  0  2\n")

        def _write_floats(vals: np.ndarray, per_row: int = 8) -> None:
            for i in range(0, len(vals), per_row):
                chunk = vals[i: i + per_row]
                lines.append("  " + "  ".join(f"{v:10.4f}" for v in chunk) + "\n")

        _write_floats(self.x_widths)
        lines.append("\n")
        _write_floats(self.z_widths)
        lines.append("    0\n")
        for row in self.cell_rows:
            lines.append(row + "\n")

        with p.open("w") as fh:
            fh.writelines(lines)
        return p

    # ------------------------------------------------------------------
    # Convenience
    # ------------------------------------------------------------------
    @property
    def n_xcells(self) -> int:
        """Number of horizontal cells."""
        return max(0, len(self.x_nodes) - 1)

    @property
    def n_zcells(self) -> int:
        """Number of vertical layers."""
        return max(0, len(self.z_nodes) - 1)

    @property
    def n_params(self) -> int:
        """Number of free model parameters (cells coded as ``'?'``)."""
        if not self.cell_rows:
            return 0
        return sum(row.count("?") for row in self.cell_rows)
