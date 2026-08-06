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
from typing import Union

import numpy as np

from .base import OccamBase
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
    num_h_nodes = int(ctrl[1])
    num_v_nodes = int(ctrl[2])
    n_airlayers = int(ctrl[3]) if len(ctrl) > 3 else 0
    n_xcells = num_h_nodes - 1
    n_zcells = num_v_nodes - 1

    # Collect float tokens until the binding sentinel ('0' alone) then
    # collect char-matrix rows for the rest of the file.
    float_tokens: list[float] = []
    cell_rows: list[str] = []
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

    x_widths = np.array(float_tokens[:n_xcells], dtype=float)
    z_widths = np.array(
        float_tokens[n_xcells : n_xcells + n_zcells], dtype=float
    )

    x_nodes = np.concatenate([[0.0], np.cumsum(x_widths)])
    z_nodes = np.concatenate([[0.0], np.cumsum(z_widths)])

    return {
        "comment": comment,
        "n_airlayers": n_airlayers,
        "x_widths": x_widths,
        "z_widths": z_widths,
        "x_nodes": x_nodes,
        "z_nodes": z_nodes,
        "cell_rows": cell_rows,
    }


# -----------------------------------------------------------------------
# OccamMesh
# -----------------------------------------------------------------------


class OccamMesh(OccamBase):
    r"""Represent the Occam2D PW2D finite-element mesh.

    ``OccamMesh`` stores the two-dimensional grid consumed by
    the Occam2D forward solver. Horizontal cell widths define
    the profile direction. Vertical widths define air and
    earth layers, and character rows mark whether cells are
    fixed, air, or boundary cells in the PW2D mesh format.

    Node coordinates are cumulative sums of cell widths:

    .. math::

        x_j = \sum_{i=0}^{j-1} \Delta x_i,
        \qquad
        z_k = \sum_{i=0}^{k-1} \Delta z_i.

    Depth :math:`z` is positive downward. Mesh construction
    uses station offsets from :class:`OccamData`, horizontal
    padding on both profile ends, optional air layers, and a
    geometrically expanding earth-layer thickness sequence.

    Parameters
    ----------
    config : OccamConfig, optional
        Configuration object controlling the number of active
        layers, number of air layers, near-surface cell sizes,
        depth scaling, and horizontal padding. If omitted, a
        default :class:`OccamConfig` is created.
    verbose : int or bool, default 0
        Verbosity level inherited from :class:`OccamBase`.
        Positive values enable progress messages through the
        instance logger.
    logger : logging.Logger, optional
        Logger used for progress and diagnostic messages. If
        omitted, a class-specific PyCSAMT logger is created.

    Attributes
    ----------
    comment : str
        First line of the mesh file, usually a provenance
        comment beginning with ``"MESH FILE"``.
    x_widths : numpy.ndarray of float, shape (n_xcells,)
        Horizontal cell widths in metres.
    z_widths : numpy.ndarray of float, shape (n_zcells,)
        Vertical layer thicknesses in metres.
    x_nodes : numpy.ndarray of float, shape (n_xcells + 1,)
        Cumulative horizontal node positions in metres.
    z_nodes : numpy.ndarray of float, shape (n_zcells + 1,)
        Cumulative depth node positions in metres, positive
        downward.
    cell_rows : list[str]
        Raw PW2D cell-type rows. Each character encodes the
        cell type at one horizontal position. The ``"?"``
        character marks cells that may contribute to free
        inversion parameters.
    n_airlayers : int
        Number of rows treated as air layers.

    Notes
    -----
    The mesh file stores widths rather than absolute node
    coordinates. :attr:`x_nodes` and :attr:`z_nodes` are
    reconstructed by cumulative summation when reading or
    building a mesh. The generated mesh uses seven padding
    cells on each side to match the boundary-column code used
    by :meth:`OccamModel.from_mesh`.

    See Also
    --------
    OccamData
        Provides station offsets used to build the mesh.
    OccamModel.from_mesh
        Converts mesh cells into inversion-parameter columns.
    InputBuilder
        Builds data, mesh, model, and startup files together.

    Examples
    --------
    Build a mesh from an Occam data file:

    >>> from pycsamt.models.occam2d import OccamData
    >>> from pycsamt.models.occam2d import OccamMesh
    >>> data = OccamData.read("occam_run/OccamDataFile.dat")
    >>> mesh = OccamMesh.from_data(data)
    >>> mesh.write("occam_run/Occam2DMesh")

    Read an existing PW2D mesh:

    >>> from pycsamt.models.occam2d import OccamMesh
    >>> mesh = OccamMesh.read("occam_run/Occam2DMesh")
    >>> mesh.n_xcells, mesh.n_zcells

    References
    ----------
    .. [OccamMesh-1] deGroot-Hedlin, C., and Constable, S.,
       "Occam's inversion to generate smooth, two-dimensional
       models from magnetotelluric data", Geophysics, 55(12),
       1613-1624, 1990.
    .. [OccamMesh-2] Constable, S. C., Parker, R. L., and Constable,
       C. G., "Occam's inversion: A practical algorithm for
       generating smooth models from electromagnetic sounding
       data", Geophysics, 52(3), 289-300, 1987.
    """

    def __init__(
        self,
        config: OccamConfig | None = None,
        **kwargs,
    ):
        super().__init__(**kwargs)
        self.config: OccamConfig = config or OccamConfig()
        self.comment: str = ""
        self.x_widths: np.ndarray = np.array([])
        self.z_widths: np.ndarray = np.array([])
        self.x_nodes: np.ndarray = np.array([])
        self.z_nodes: np.ndarray = np.array([])
        self.cell_rows: list[str] = []
        self.n_airlayers: int = self.config.n_airlayers

    # ------------------------------------------------------------------
    # Construction
    # ------------------------------------------------------------------
    @classmethod
    def from_data(
        cls,
        data: OccamData,  # noqa: F821
        config: OccamConfig | None = None,
        **kwargs,
    ) -> OccamMesh:
        """Build a PW2D mesh from Occam data offsets.

        The method creates a finite-element mesh spanning the
        profile described by ``data.offsets``. It uses seven
        padding cells on each side, station-zone cells near
        the configured horizontal cell size, optional air
        layers, and geometrically expanding earth layers.

        The horizontal padding is chosen to match the boundary
        columns used by :meth:`OccamModel.from_mesh`. Interior
        station-zone cells are adjusted so the number of cells
        remains compatible with the model parameter grouping.

        Parameters
        ----------
        data : OccamData
            Populated data object. Its ``offsets`` array must
            contain station chainages in metres. Offsets are
            sorted before mesh construction.
        config : OccamConfig, optional
            Configuration object controlling mesh geometry.
            The builder uses ``cell_size_horizontal``,
            ``n_airlayers``, ``n_layers``,
            ``cell_size_vertical_top``, and ``depth_scale``.
            If omitted, a default :class:`OccamConfig` is
            created.
        **kwargs
            Additional keyword arguments forwarded to the
            ``OccamMesh`` constructor. Use this for
            ``verbose`` or ``logger``.

        Returns
        -------
        OccamMesh
            Mesh object ready to be written as ``Occam2DMesh``
            or passed to :meth:`OccamModel.from_mesh`.

        Raises
        ------
        ValueError
            Raised when the data object contains no station
            offsets.

        See Also
        --------
        OccamData.from_edi
            Creates the offsets used by this method.
        OccamModel.from_mesh
            Builds the inversion-parameter mapping.

        Examples
        --------
        >>> from pycsamt.models.occam2d import OccamData
        >>> from pycsamt.models.occam2d import OccamMesh
        >>> data = OccamData.read("OccamDataFile.dat")
        >>> mesh = OccamMesh.from_data(data)
        >>> mesh.n_airlayers
        """
        cfg = config or OccamConfig()
        obj = cls(config=cfg, **kwargs)

        offsets = np.sort(np.asarray(data.offsets, dtype=float))
        if offsets.size == 0:
            raise ValueError(
                "OccamMesh.from_data: data has no station offsets"
            )

        cell_h = cfg.cell_size_horizontal
        n_pad = 7  # fixed: matches boundary code 7 in OccamModel.from_mesh

        # ---- station-zone x-cells ----------------------------------------
        station_widths: list[float] = []
        if offsets.size == 1:
            station_widths = [cell_h]
        else:
            for i in range(len(offsets) - 1):
                gap = float(offsets[i + 1] - offsets[i])
                n_cell = max(1, round(gap / cell_h))
                w = gap / n_cell
                station_widths.extend([w] * n_cell)

        # Ensure even count so that (n_xcells - 14) is even for code-2 interior cols
        if len(station_widths) % 2 != 0:
            station_widths.append(
                station_widths[-1] if station_widths else cell_h
            )

        # ---- horizontal padding (7 cells each side, geometrically expanding)
        pad = [cell_h * float(2 ** (k + 1)) for k in range(n_pad)]
        left_pad = list(reversed(pad))  # widest cell at outer edge
        right_pad = pad

        x_widths = np.array(left_pad + station_widths + right_pad, dtype=float)

        # ---- vertical cells ----------------------------------------------
        n_air = cfg.n_airlayers
        n_active = cfg.n_layers
        air_h = cfg.cell_size_vertical_top

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
        obj.comment = "MESH FILE Created by pycsamt.models.occam2d"
        obj.x_widths = x_widths
        obj.z_widths = z_widths
        obj.x_nodes = np.concatenate([[0.0], np.cumsum(x_widths)])
        obj.z_nodes = np.concatenate([[0.0], np.cumsum(z_widths)])
        obj.cell_rows = cell_rows
        obj.n_airlayers = n_air

        if obj.verbose:
            obj.logger.info(
                "OccamMesh.from_data: %d×%d cells (%d airlayers) from %d stations",
                n_xcells,
                n_zcells,
                n_air,
                offsets.size,
            )
        return obj

    # ------------------------------------------------------------------
    # I/O
    # ------------------------------------------------------------------
    @classmethod
    def read(cls, path: PathLike, **kwargs) -> OccamMesh:
        """Read an existing ``Occam2DMesh`` PW2D file.

        The reader parses the comment line, control line,
        horizontal widths, vertical widths, air-layer count,
        and cell-type character rows. Node arrays are rebuilt
        from cumulative sums of widths.

        Parameters
        ----------
        path : path-like
            Path to the mesh file. The value may be a string,
            :class:`pathlib.Path`, or any object accepted by
            :class:`pathlib.Path`.
        **kwargs
            Additional keyword arguments forwarded to the
            ``OccamMesh`` constructor before parsed values are
            attached. Use this for ``config``, ``verbose``, or
            ``logger``.

        Returns
        -------
        OccamMesh
            Parsed mesh container with widths, nodes, and cell
            rows populated.

        Raises
        ------
        FileNotFoundError
            Raised when ``path`` does not exist.
        ValueError
            Raised when the file is too short or the control
            line cannot be parsed.

        Examples
        --------
        >>> from pycsamt.models.occam2d import OccamMesh
        >>> mesh = OccamMesh.read("occam_run/Occam2DMesh")
        >>> mesh.x_nodes.shape
        """
        p = Path(path)
        d = _parse_mesh(p)
        obj = cls(**kwargs)

        obj.comment = d["comment"]
        obj.n_airlayers = d["n_airlayers"]
        obj.x_widths = d["x_widths"]
        obj.z_widths = d["z_widths"]
        obj.x_nodes = d["x_nodes"]
        obj.z_nodes = d["z_nodes"]
        obj.cell_rows = d["cell_rows"]

        if obj.verbose:
            obj.logger.info(
                "OccamMesh.read: %d×%d cells (%d char rows) from %s",
                obj.n_xcells,
                obj.n_zcells,
                len(obj.cell_rows),
                p,
            )
        return obj

    def write(self, path: PathLike) -> Path:
        """Write this mesh in PW2D format.

        The writer serializes the current comment, control
        values, horizontal widths, vertical widths, and
        cell-type rows to the native Occam2D mesh format.
        Parent directories are created before writing.

        Parameters
        ----------
        path : path-like
            Destination path for the mesh file. The value may
            be a string, :class:`pathlib.Path`, or any object
            accepted by :class:`pathlib.Path`.

        Returns
        -------
        pathlib.Path
            Path to the file that was written.

        See Also
        --------
        OccamMesh.read
            Parses mesh files written by this method.
        InputBuilder.build
            Calls this method during input-file generation.

        Examples
        --------
        >>> from pycsamt.models.occam2d import OccamMesh
        >>> mesh = OccamMesh.read("source/Occam2DMesh")
        >>> written = mesh.write("copy/Occam2DMesh")
        """
        p = Path(path)
        p.parent.mkdir(parents=True, exist_ok=True)

        n_h = self.n_xcells + 1  # num_h_nodes
        n_v = self.n_zcells + 1  # num_v_nodes
        comment = self.comment or "MESH FILE Created by pycsamt"

        lines: list[str] = []
        lines.append(f"{comment}\n")
        lines.append(f"   0  {n_h}  {n_v}  {self.n_airlayers}  0  2\n")

        def _write_floats(vals: np.ndarray, per_row: int = 8) -> None:
            for i in range(0, len(vals), per_row):
                chunk = vals[i : i + per_row]
                lines.append(
                    "  " + "  ".join(f"{v:10.4f}" for v in chunk) + "\n"
                )

        _write_floats(self.x_widths)
        lines.append("\n")
        _write_floats(self.z_widths)
        lines.append("    0\n")
        for row in self.cell_rows:
            lines.append(row + "\n")

        with p.open("w") as fh:
            fh.writelines(lines)
        self.path = p
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
