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
    0  num_h_nodes  num_v_nodes  n_fixed_resistivities  0  mesh_type
    <x-cell widths: num_h_nodes-1 floats, wrapped across any number of lines>
    <blank line — optional>
    <z-layer thicknesses: num_v_nodes-1 floats, wrapped across any number of lines>
    0
    <cell-type character rows: 4*(num_v_nodes-1) rows of num_h_nodes-1 characters
     each, 4 identical rows per z-layer>
        '?' = free model parameter, other = fixed / air / boundary

The 4th control-line field is the number of *fixed resistivities*
(``nrfix`` in the reference Fortran reader), not the air-layer count --
pyCSAMT does not write fixed-resistivity blocks, so this is always ``0``.
Air layers are not stored as a header count at all; they are inferred from
the cell-type rows themselves (leading z-layers whose 4 rows are entirely
``'0'``).
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

    # Line 1: control integers —
    # "0  num_h_nodes  num_v_nodes  n_fixed_resistivities  0  mesh_type"
    # The 4th field is nrfix (fixed resistivities), not an air-layer count;
    # see the module docstring. Air layers are recovered below from the
    # cell-type rows themselves.
    ctrl = lines[1].strip().split()
    if len(ctrl) < 3:
        raise ValueError(f"Cannot parse mesh control line: {lines[1]!r}")
    num_h_nodes = int(ctrl[1])
    num_v_nodes = int(ctrl[2])
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

    # Air layers are written as 4 identical all-'0' rows per z-layer
    # (see from_data()); recover the count by scanning leading layers.
    n_airlayers = 0
    for layer_start in range(0, len(cell_rows) - 3, 4):
        layer_rows = cell_rows[layer_start : layer_start + 4]
        if all(set(row) <= {"0"} for row in layer_rows):
            n_airlayers += 1
        else:
            break

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

    #: Fixed horizontal padding-cell count on each profile end. Matches
    #: the boundary code 7 used by :meth:`OccamModel.from_mesh` and the
    #: hardcoded ``n_pad`` in :meth:`from_data` -- an architectural
    #: invariant, not something inferred per-mesh, so it also applies to
    #: meshes rebuilt via :meth:`read`.
    N_PAD: int = 7

    def cell_centers_survey_x(self) -> np.ndarray:
        """Horizontal cell-center coordinates in survey (offset) space.

        :attr:`x_widths`/:attr:`x_nodes` are zero-based at the *outer*
        edge of the left padding, not at the survey's own zero offset
        (see :class:`OccamData`'s ``offsets``, which starts at 0 after
        normalization). Comparing or resampling a solved model against
        anything expressed in survey/offset coordinates -- true models,
        AI display grids, station chainage -- must shift by the total
        left-padding width, or the comparison silently lands entirely
        inside the padding zone.

        Returns
        -------
        numpy.ndarray of float, shape (n_xcells,)
            Cell-center x-coordinates, in metres, in the same
            zero-based convention as :class:`OccamData`'s ``offsets``.
        """
        centers = np.cumsum(self.x_widths) - 0.5 * self.x_widths
        left_pad_width = float(np.sum(self.x_widths[: self.N_PAD]))
        return centers - left_pad_width

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
            ``n_airlayers``, ``n_layers``, ``max_depth``,
            ``cell_size_vertical_top``, and ``depth_scale``.
            Earth layers expand geometrically from
            ``cell_size_vertical_top`` by ``depth_scale`` and
            stop -- truncating the last layer if needed -- once
            cumulative depth reaches ``max_depth`` (default
            1500 m) or ``n_layers`` layers have been added,
            whichever comes first. If omitted, a default
            :class:`OccamConfig` is created.
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
        n_pad = cls.N_PAD  # fixed: matches boundary code 7 in OccamModel.from_mesh

        # Two stations closer together than this are for practical purposes
        # coincident (e.g. a repeat/QC measurement at the same physical
        # site under a different station suffix). Left unguarded, the gap
        # still produces one degenerate sub-metre-wide station-zone cell
        # (max(1, round(gap/cell_h)) never returns 0), which builds and
        # writes without error here but breaks the external Occam2D
        # Fortran solver's own free-parameter count deep inside its own
        # I/O ("Model file error: # params <> # free bricks") with no
        # indication of which station caused it. Fail fast here instead,
        # with an actionable message, rather than downstream in a
        # compiled binary's stderr.
        min_gap = 0.01 * cell_h
        if offsets.size > 1:
            gaps = np.diff(offsets)
            too_close = np.nonzero(gaps < min_gap)[0]
            if too_close.size:
                pairs = ", ".join(
                    f"offsets[{i}]={offsets[i]:.3f} m, "
                    f"offsets[{i + 1}]={offsets[i + 1]:.3f} m "
                    f"(gap={gaps[i]:.3f} m)"
                    for i in too_close
                )
                raise ValueError(
                    "OccamMesh.from_data: station offsets closer than "
                    f"{min_gap:.3f} m (1% of cell_size_horizontal="
                    f"{cell_h:g} m) found: {pairs}. These are almost "
                    "certainly repeat/QC measurements at the same "
                    "physical site rather than distinct along-line "
                    "positions; drop one station from each such pair "
                    "before building the mesh."
                )

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

        # ---- vertical cells ------------------------------------------------
        # Occam2D's own Fortran solver adds its own air layer(s) internally
        # (the `addair` subroutine in MT2D.f90) purely for the forward
        # boundary condition; a properly formed Occam2DMesh/Occam2DModel
        # file pair does not represent air in the finite-element region
        # grid at all -- confirmed against a real, working mtpy-generated
        # reference (data/occam2D/Tonkeng), where every mesh z-row is a
        # free ("?") earth cell and sum(irz) across the declared model
        # layers equals the mesh's total z-cell count exactly, with no
        # offset for air. Writing `cfg.n_airlayers` extra fixed ("0") rows
        # here, as an earlier version of this method did, desynchronizes
        # that count: `OccamModel.from_mesh` declares only the earth
        # layers (`mesh.n_zcells - mesh.n_airlayers`), but the free-brick
        # counting loop inside Occam2D's own MT2D.f90 walks mesh z-rows
        # starting from row 1 for exactly `nlay` (the *declared* earth
        # layer count) steps -- so with extra unaccounted air rows
        # prepended, it silently reads across the wrong z-range and comes
        # up short by exactly `n_airlayers` layers' worth of free
        # parameters, aborting with "Model file error: # params <> # free
        # bricks" deep inside the compiled binary. `cfg.n_airlayers` /
        # `cfg.cell_size_vertical_top` are therefore not applicable to the
        # Occam2D mesh file and are intentionally not used below; they
        # remain on `OccamConfig` only because the same config class is
        # shared with the 2-D startup/model builders' docstrings and CLI.
        n_active = cfg.n_layers
        max_depth = float(getattr(cfg, "max_depth", 0.0) or 0.0)

        z_w: list[float] = []
        thick = float(cfg.cell_size_vertical_top)
        cum_depth = 0.0
        for _ in range(n_active):
            if max_depth > 0.0 and cum_depth >= max_depth:
                break
            if max_depth > 0.0:
                thick = min(thick, max_depth - cum_depth)
            if thick <= 0.0:
                break
            z_w.append(thick)
            cum_depth += thick
            thick *= cfg.depth_scale

        z_widths = np.array(z_w, dtype=float)

        # ---- char matrix (4 rows per z-cell, all free earth cells) --------
        n_xcells = len(x_widths)
        n_zcells = len(z_widths)
        cell_rows: list[str] = []
        for _iz in range(n_zcells):
            row = "?" * n_xcells
            for _ in range(4):
                cell_rows.append(row)

        # ---- populate object --------------------------------------------
        obj.comment = "MESH FILE Created by pycsamt.models.occam2d"
        obj.x_widths = x_widths
        obj.z_widths = z_widths
        obj.x_nodes = np.concatenate([[0.0], np.cumsum(x_widths)])
        obj.z_nodes = np.concatenate([[0.0], np.cumsum(z_widths)])
        obj.cell_rows = cell_rows
        obj.n_airlayers = 0

        if obj.verbose:
            obj.logger.info(
                "OccamMesh.from_data: %d×%d cells (no air rows written; "
                "Occam2D adds its own internally) from %d stations",
                n_xcells,
                n_zcells,
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
        # 4th field is nrfix (fixed resistivities), always 0 here — see the
        # module docstring. It is NOT the air-layer count.
        lines.append(f"   0  {n_h}  {n_v}  0  0  2\n")

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


def resample_rho_to_grid(
    rho_2d: np.ndarray,
    mesh: OccamMesh,
    x: np.ndarray,
    z: np.ndarray,
) -> np.ndarray:
    """Resample a solved Occam2D resistivity model onto a regular grid.

    ``rho_2d`` is defined on the mesh's own irregular, padding-inclusive
    cell grid. ``x``/``z`` are typically a regular display or comparison
    grid expressed in survey coordinates (station-offset-relative, e.g.
    an AI benchmark's true-model grid or a fixed display grid) -- *not*
    the mesh's own zero-at-outer-padding coordinate system. Resampling
    must therefore go through :meth:`OccamMesh.cell_centers_survey_x`,
    not the mesh's raw node coordinates, or every query lands inside the
    (typically tens of kilometres of) horizontal padding and is clamped
    to a single, flat, constant-per-row value regardless of the model's
    real structure.

    Parameters
    ----------
    rho_2d : numpy.ndarray, shape (n_zcells, n_xcells)
        Resistivity (or log-resistivity) on the mesh's own grid,
        including any air-layer rows at the top (stripped internally
        using ``mesh.n_airlayers``).
    mesh : OccamMesh
        Mesh the model was solved on.
    x, z : numpy.ndarray
        Target cell-center coordinates, in metres, in survey (offset)
        and depth (positive-down) coordinates respectively.

    Returns
    -------
    numpy.ndarray, shape (z.size, x.size)
        ``rho_2d`` resampled onto the ``x``/``z`` grid.
    """
    n_air = int(mesh.n_airlayers)
    earth = np.asarray(rho_2d, dtype=float)[n_air:, :]
    z_widths = np.asarray(mesh.z_widths, dtype=float)[n_air:]
    source_x = mesh.cell_centers_survey_x()
    source_z = np.cumsum(z_widths) - 0.5 * z_widths
    horizontal = np.vstack([np.interp(x, source_x, row) for row in earth])
    return np.vstack(
        [np.interp(z, source_z, horizontal[:, c]) for c in range(x.size)]
    ).T
