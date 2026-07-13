# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Read, write, and build ModEM 3-D model files.

The module defines :class:`ModEmModel3D`, the container for
three-dimensional ModEM starting and iteration models.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Union

import numpy as np

from .base import ModEmBase
from .config import ModEmConfig
from .doc import _modem_param_docs as _params

if TYPE_CHECKING:
    from .data import ModEmData

PathLike = Union[str, Path]

__all__ = ["ModEmModel3D"]


# ----------------------------------------------------------------------
# Low-level parser
# ----------------------------------------------------------------------


def _parse_model3d(path: Path) -> dict:
    """Parse a ModEM 3-D model file into plain Python values."""
    with path.open("r", errors="replace") as fh:
        lines = [ln.rstrip("\n") for ln in fh]

    N = len(lines)
    i = 0
    # skip blank / comment lines before the control line
    while i < N and (
        not lines[i].strip() or lines[i].strip().startswith("#")
    ):
        i += 1

    ctrl = lines[i].split()
    i += 1
    nx = int(ctrl[0])
    ny = int(ctrl[1])
    nz = int(ctrl[2])
    if len(ctrl) > 3 and ctrl[3].lstrip("-").isdigit():
        n_air = int(ctrl[3])
    else:
        n_air = 0
    # log_type token is the last alphabetic token in ctrl
    log_type = "LOGE"
    for tok in ctrl:
        if tok.isalpha() or tok.upper() in ("LOGE", "LOG10", "LINEAR"):
            log_type = tok.upper()

    # collect all float tokens until we have nx + ny + nz of them
    float_tokens: list[float] = []
    target = nx + ny + nz
    while i < N and len(float_tokens) < target:
        ln = lines[i].strip()
        i += 1
        if not ln or ln.startswith("#"):
            continue
        for tok in ln.split():
            try:
                float_tokens.append(float(tok))
            except ValueError:
                pass

    x_widths = np.array(float_tokens[:nx], dtype=float)
    y_widths = np.array(float_tokens[nx : nx + ny], dtype=float)
    z_widths = np.array(float_tokens[nx + ny : target], dtype=float)

    # read rho grid: Nz * Ny rows of Nx values
    total = nx * ny * nz
    rho_flat: list[float] = []
    while i < N and len(rho_flat) < total:
        ln = lines[i].strip()
        i += 1
        if not ln or ln.startswith("#"):
            continue
        # Stop at possible centre coordinates or rotation lines.
        parts = ln.split()
        if len(parts) <= 3 and all(_is_simple_numeric(p) for p in parts):
            # Could be trailing centre/rotation. Stop only if we have
            # most of the data we need
            if len(rho_flat) >= total - nx:
                break
        for tok in parts:
            try:
                rho_flat.append(float(tok))
            except ValueError:
                pass

    rho_arr = np.array(
        rho_flat[:total],
        dtype=float,
    ).reshape(nz, ny, nx)

    if log_type == "LOG10":
        rho_loge = rho_arr * np.log(10.0)
    elif log_type == "LINEAR":
        with np.errstate(divide="ignore", invalid="ignore"):
            rho_loge = np.log(np.where(rho_arr > 0, rho_arr, np.nan))
    else:
        rho_loge = rho_arr  # already LOGE

    return {
        "nx": nx,
        "ny": ny,
        "nz": nz,
        "n_air": n_air,
        "log_type": log_type,
        "x_widths": x_widths,
        "y_widths": y_widths,
        "z_widths": z_widths,
        "rho_loge": rho_loge,
    }


def _is_simple_numeric(s: str) -> bool:
    try:
        float(s)
        return True
    except ValueError:
        return False


# ----------------------------------------------------------------------
# ModEmModel3D
# ----------------------------------------------------------------------


class ModEmModel3D(ModEmBase):
    def __init__(self, config: ModEmConfig | None = None, **kwargs):
        """Initialize an empty 3-D ModEM model container."""
        super().__init__(**kwargs)
        self.config: ModEmConfig = config or ModEmConfig()
        self.x_widths: np.ndarray = np.array([])
        self.y_widths: np.ndarray = np.array([])
        self.z_widths: np.ndarray = np.array([])
        self.rho_loge: np.ndarray = np.empty((0, 0, 0))
        self.n_air: int = 0
        self.log_type: str = "LOGE"

    # ------------------------------------------------------------------
    # Derived
    # ------------------------------------------------------------------
    @property
    def nx(self) -> int:
        """Number of model cells along the x direction."""
        return len(self.x_widths)

    @property
    def ny(self) -> int:
        """Number of model cells along the y direction."""
        return len(self.y_widths)

    @property
    def nz(self) -> int:
        """Number of vertical layers, including air layers."""
        return len(self.z_widths)

    @property
    def x_nodes(self) -> np.ndarray:
        """Cumulative x-node coordinates in metres."""
        return np.concatenate([[0.0], np.cumsum(self.x_widths)])

    @property
    def y_nodes(self) -> np.ndarray:
        """Cumulative y-node coordinates in metres."""
        return np.concatenate([[0.0], np.cumsum(self.y_widths)])

    @property
    def z_nodes(self) -> np.ndarray:
        """Cumulative vertical node depths in metres."""
        return np.concatenate([[0.0], np.cumsum(self.z_widths)])

    @property
    def rho_linear(self) -> np.ndarray:
        """Return resistivity in linear ohm metres."""
        return np.exp(self.rho_loge)

    @property
    def shape(self) -> tuple[int, int, int]:
        """Return ``(nz, ny, nx)`` for the resistivity grid."""
        return (self.nz, self.ny, self.nx)

    # ------------------------------------------------------------------
    # Construction
    # ------------------------------------------------------------------
    @classmethod
    def halfspace(
        cls,
        data: ModEmData,
        config: ModEmConfig | None = None,
        **kwargs,
    ) -> ModEmModel3D:
        """Build a uniform half-space starting model.

        The method derives the horizontal x and y grids from
        station coordinates in a :class:`ModEmData` object. It
        fills earth cells with ``config.initial_rho`` and assigns
        air layers a high fixed resistivity. Resistivity is
        stored internally as natural logarithm values.

        Parameters
        ----------
        data : ModEmData
            Populated data object. The builder uses
            ``data.x_coords`` and ``data.y_coords`` to define
            the station-zone grid in the horizontal plane.
        config : ModEmConfig, optional
            Configuration supplying 3-D grid dimensions,
            padding, vertical growth, air-layer count, and
            starting half-space resistivity. If omitted, a
            default :class:`ModEmConfig` is used.
        **kwargs : dict
            Additional keyword arguments forwarded to
            :class:`ModEmModel3D`.

        Returns
        -------
        ModEmModel3D
            Populated 3-D model object ready to be written as a
            ModEM ``.ws`` model file.

        Examples
        --------
        >>> from pycsamt.models.modem.model3d import ModEmModel3D
        >>> model = ModEmModel3D.halfspace(data, config=cfg)
        >>> model.shape == model.rho_loge.shape
        True

        Notes
        -----
        The station-zone widths are derived separately for x
        and y from sorted unique station coordinates. If only
        one coordinate exists in a direction, one station-zone
        cell is created using ``config.cell_size_h``. Padding
        cells grow by powers of two away from the station zone.
        """
        cfg = config or ModEmConfig()
        obj = cls(config=cfg, **kwargs)

        # -- horizontal grids ------------------------------------------
        if hasattr(data, "x_coords"):
            xs = np.sort(np.unique(data.x_coords))
        else:
            xs = np.array([0.0])
        if hasattr(data, "y_coords"):
            ys = np.sort(np.unique(data.y_coords))
        else:
            ys = np.array([0.0])

        cell_h = cfg.cell_size_h
        n_pad = cfg.n_padding_xy

        def _station_cells(coords: np.ndarray) -> list:
            if coords.size <= 1:
                return [cell_h]
            widths = []
            for i in range(len(coords) - 1):
                gap = float(coords[i + 1] - coords[i])
                n_cell = max(1, round(gap / cell_h))
                widths.extend([gap / n_cell] * n_cell)
            return widths

        pad = [cell_h * float(2 ** (k + 1)) for k in range(n_pad)]

        x_station = _station_cells(xs)
        y_station = _station_cells(ys)

        x_w = np.array(
            list(reversed(pad)) + x_station + pad,
            dtype=float,
        )
        y_w = np.array(
            list(reversed(pad)) + y_station + pad,
            dtype=float,
        )

        # -- vertical grid ---------------------------------------------
        n_air = cfg.n_airlayers
        n_active = cfg.nz
        air_h = cfg.cell_size_v_top
        air_z = [air_h] * n_air
        earth_z: list[float] = []
        thick = float(cfg.cell_size_v_top)
        for _ in range(n_active):
            earth_z.append(thick)
            thick *= cfg.depth_scale
        z_w = np.array(air_z + earth_z, dtype=float)

        nz_tot = len(z_w)
        ny_tot = len(y_w)
        nx_tot = len(x_w)

        rho_val = float(np.log(cfg.initial_rho))
        rho_grid = np.full(
            (nz_tot, ny_tot, nx_tot),
            rho_val,
            dtype=float,
        )
        if n_air > 0:
            rho_grid[:n_air, :, :] = np.log(1e12)

        obj.x_widths = x_w
        obj.y_widths = y_w
        obj.z_widths = z_w
        obj.n_air = n_air
        obj.rho_loge = rho_grid

        if obj.verbose:
            obj.logger.info(
                "ModEmModel3D.halfspace: %d x %d x %d grid, rho=%.1f ohm m",
                nx_tot,
                ny_tot,
                nz_tot,
                cfg.initial_rho,
            )
        return obj

    # ------------------------------------------------------------------
    # I/O
    # ------------------------------------------------------------------
    @classmethod
    def read(cls, path: PathLike, **kwargs) -> ModEmModel3D:
        """Parse an existing ModEM 3-D model file.

        Parameters
        ----------
        path : path-like
            Path to a ModEM 3-D ``.ws`` model file. The file may
            store resistivity as ``LOGE``, ``LOG10``, or
            ``LINEAR``; values are converted to natural-log
            resistivity internally.
        **kwargs : dict
            Additional keyword arguments forwarded to
            :class:`ModEmModel3D`.

        Returns
        -------
        ModEmModel3D
            Parsed model with x, y, and z cell widths,
            natural-log resistivity, air-layer count, and
            log-type metadata.

        Raises
        ------
        FileNotFoundError
            If ``path`` does not exist.

        Examples
        --------
        >>> from pycsamt.models.modem.model3d import ModEmModel3D
        >>> model = ModEmModel3D.read("m0.ws")
        >>> model.rho_loge.shape == model.shape
        True
        """
        p = Path(path)
        if not p.exists():
            msg = f"ModEM 3D model file not found: {p}"
            raise FileNotFoundError(msg)
        d = _parse_model3d(p)
        obj = cls(**kwargs)
        obj.x_widths = d["x_widths"]
        obj.y_widths = d["y_widths"]
        obj.z_widths = d["z_widths"]
        obj.rho_loge = d["rho_loge"]
        obj.n_air = d["n_air"]
        obj.log_type = d["log_type"]
        if obj.verbose:
            obj.logger.info(
                "ModEmModel3D.read: %d x %d x %d from %s",
                obj.nx,
                obj.ny,
                obj.nz,
                p,
            )
        return obj

    def write(self, path: PathLike) -> Path:
        """Write the model to ``path`` in ModEM 3-D format.

        Parameters
        ----------
        path : path-like
            Destination model file. Parent directories are
            created before writing. Existing files are
            overwritten.

        Returns
        -------
        pathlib.Path
            Path passed to the writer, converted to
            :class:`pathlib.Path`.

        Notes
        -----
        The writer emits the WinGLink/ModEM ``.ws`` style grid:
        header, x widths, y widths, z widths, then
        ``nz * ny`` rows of ``nx`` resistivity values. It writes
        the values stored in ``rho_loge`` and the encoding label
        stored in ``log_type``.

        Examples
        --------
        >>> from pycsamt.models.modem.model3d import ModEmModel3D
        >>> model = ModEmModel3D.halfspace(data, config=cfg)
        >>> path = model.write("m0.ws")
        >>> path.name
        'm0.ws'
        """
        p = Path(path)
        p.parent.mkdir(parents=True, exist_ok=True)

        def _floats(arr: np.ndarray, per_row: int = 12) -> str:
            rows = []
            for i in range(0, len(arr), per_row):
                rows.append(
                    "  "
                    + "  ".join(f"{v:>12.4f}" for v in arr[i : i + per_row])
                )
            return "\n".join(rows) + "\n"

        lines: list[str] = [
            (
                f"  {self.nx}  {self.ny}  {self.nz}  "
                f"{self.n_air}  {self.log_type}\n"
            ),
        ]
        lines.append(_floats(self.x_widths))
        lines.append(_floats(self.y_widths))
        lines.append(_floats(self.z_widths))
        lines.append("\n")

        for iz in range(self.nz):
            for iy in range(self.ny):
                lines.append(
                    "  "
                    + "  ".join(
                        f"{v:>12.5E}" for v in self.rho_loge[iz, iy, :]
                    )
                    + "\n"
                )

        with p.open("w") as fh:
            fh.writelines(lines)
        return p


ModEmModel3D.__doc__ = rf"""
Represent a ModEM three-dimensional resistivity model.

``ModEmModel3D`` stores the grid and resistivity values used by
the ModEM 3-D executable. The model contains x, y, and z cell
widths, a count of air layers at the top of the mesh, and a
resistivity volume with shape ``(nz, ny, nx)``. Resistivity is
stored internally as natural logarithms, so a linear
resistivity :math:`\rho` is represented as

.. math::

   m = \ln(\rho).

The ASCII ``.ws`` format contains a header with ``nx``, ``ny``,
``nz``, optional air-layer count, and log encoding. It then
stores cell widths followed by the resistivity volume written
layer by layer.

Parameters
----------
{_params.common.config}
{_params.common.verbose}
{_params.common.logger}

Attributes
----------
config : ModEmConfig
    Configuration object used when constructing half-space
    models.
x_widths : numpy.ndarray, shape (nx,)
    Cell widths in metres along the ModEM x direction.
y_widths : numpy.ndarray, shape (ny,)
    Cell widths in metres along the ModEM y direction.
z_widths : numpy.ndarray, shape (nz,)
    Vertical layer thicknesses in metres. The array includes
    air layers followed by active earth layers.
rho_loge : numpy.ndarray, shape (nz, ny, nx)
    Natural-log resistivity values. Earth cells in a default
    half-space are initialized as
    :math:`\ln(\rho_0)`, where :math:`\rho_0` is
    ``config.initial_rho``.
n_air : int
    Number of air layers included at the top of the model.
log_type : str, default "LOGE"
    Encoding label written to the ModEM file header. ``"LOGE"``
    is the standard internal representation used by this
    class. Readers also accept ``"LOG10"`` and ``"LINEAR"`` and
    convert them to natural logarithms.

Derived Properties
------------------
nx, ny, nz : int
    Number of cells in the x, y, and z directions.
x_nodes, y_nodes, z_nodes : numpy.ndarray
    Cumulative node coordinates in metres for each direction.
rho_linear : numpy.ndarray, shape (nz, ny, nx)
    Resistivity values in linear ohm metres, obtained as
    :math:`\exp(\mathtt{{rho\_loge}})`.
shape : tuple of int
    Convenience tuple ``(nz, ny, nx)`` matching
    ``rho_loge.shape``.

Notes
-----
The :meth:`halfspace` constructor creates a uniform starting
model. Air layers are assigned a very high resistivity, while
earth cells are assigned ``config.initial_rho``. Horizontal
padding grows away from the station footprint in both x and y
directions so artificial boundaries are moved away from the
survey.

See Also
--------
ModEmConfig
    Supplies 3-D grid and initial-resistivity settings.
ModEmData
    Provides station coordinates used by :meth:`halfspace`.
ModEmCovariance
    Uses this model geometry to build 3-D covariance masks.
InputBuilder
    Builds and writes a 3-D starting model automatically.
ModEmRunner
    Passes the written model file to the ModEM executable.

Examples
--------
Build a 3-D half-space model:

>>> from pycsamt.models.modem.config import ModEmConfig
>>> from pycsamt.models.modem.model3d import ModEmModel3D
>>> cfg = ModEmConfig(mode="3d", initial_rho=100.0)
>>> model = ModEmModel3D.halfspace(data, config=cfg)
>>> model.shape == model.rho_loge.shape
True

Read and write a ModEM ``.ws`` model file:

>>> model = ModEmModel3D.read("m0.ws")
>>> path = model.write("m0_copy.ws")
>>> path.name
'm0_copy.ws'

Inspect linear resistivity:

>>> rho = model.rho_linear
>>> rho.shape == model.shape
True

References
----------
.. [1] Egbert, G. D., and Kelbert, A., "Computational
   recipes for electromagnetic inverse problems", Geophysical
   Journal International, 189(1), 251-267, 2012,
   doi:10.1111/j.1365-246X.2011.05347.x.
.. [2] Kelbert, A., Meqbel, N., Egbert, G. D., and Tandon,
   K., "ModEM: A modular system for inversion of
   electromagnetic geophysical data", Computers and
   Geosciences, 66, 40-53, 2014,
   doi:10.1016/j.cageo.2014.01.010.
"""
