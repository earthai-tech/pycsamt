# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Read, write, and build ModEM 2-D model files.

The module defines :class:`ModEmModel2D`, the container for
two-dimensional ModEM starting or iteration models.
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

__all__ = ["ModEmModel2D"]


# ----------------------------------------------------------------------
# Low-level parser
# ----------------------------------------------------------------------


def _parse_model2d(path: Path) -> dict:
    """Parse a ModEM 2-D model file into plain Python values."""
    with path.open("r", errors="replace") as fh:
        lines = [line.rstrip("\n") for line in fh]

    # Skip blank / comment lines to find control line
    i = 0
    N = len(lines)
    while i < N and (not lines[i].strip() or lines[i].strip().startswith("#")):
        i += 1

    ctrl = lines[i].split()
    i += 1
    nx = int(ctrl[0])
    nz = int(ctrl[1])
    log_type = ctrl[2].upper() if len(ctrl) > 2 else "LOGE"

    # Collect floats until the model block-count sentinel.
    float_tokens: list[float] = []
    while i < N:
        ln = lines[i].strip()
        i += 1
        if not ln or ln.startswith("#"):
            continue
        # Block-count sentinel: a lone integer
        parts = ln.split()
        if len(parts) == 1:
            try:
                int(parts[0])
                break
            except ValueError:
                pass
        for tok in parts:
            try:
                float_tokens.append(float(tok))
            except ValueError:
                pass

    x_widths = np.array(float_tokens[:nx], dtype=float)
    z_widths = np.array(float_tokens[nx : nx + nz], dtype=float)

    # Read rho grid: nz rows each of nx values
    rho_flat: list[float] = []
    while i < N and len(rho_flat) < nx * nz:
        ln = lines[i].strip()
        i += 1
        if not ln or ln.startswith("#"):
            continue
        for tok in ln.split():
            try:
                rho_flat.append(float(tok))
            except ValueError:
                pass

    rho_log = np.array(rho_flat[: nx * nz], dtype=float).reshape(nz, nx)

    # Convert to loge if needed
    if log_type == "LOG10":
        rho_loge = rho_log * np.log(10.0)
    elif log_type == "LINEAR":
        with np.errstate(divide="ignore", invalid="ignore"):
            rho_loge = np.log(np.where(rho_log > 0, rho_log, np.nan))
    else:
        rho_loge = rho_log  # already LOGE

    return {
        "nx": nx,
        "nz": nz,
        "log_type": log_type,
        "x_widths": x_widths,
        "z_widths": z_widths,
        "rho_loge": rho_loge,
    }


# ----------------------------------------------------------------------
# ModEmModel2D
# ----------------------------------------------------------------------


class ModEmModel2D(ModEmBase):
    def __init__(self, config: ModEmConfig | None = None, **kwargs):
        """Initialize an empty 2-D ModEM model container."""
        super().__init__(**kwargs)
        self.config: ModEmConfig = config or ModEmConfig()
        self.x_widths: np.ndarray = np.array([])
        self.z_widths: np.ndarray = np.array([])
        self.rho_loge: np.ndarray = np.empty((0, 0))
        self.log_type: str = "LOGE"

    # ------------------------------------------------------------------
    # Derived
    # ------------------------------------------------------------------
    @property
    def nx(self) -> int:
        """Number of horizontal model cells."""
        return len(self.x_widths)

    @property
    def nz(self) -> int:
        """Number of vertical layers, including air layers."""
        return len(self.z_widths)

    @property
    def x_nodes(self) -> np.ndarray:
        """Cumulative horizontal node coordinates in metres."""
        return np.concatenate([[0.0], np.cumsum(self.x_widths)])

    @property
    def z_nodes(self) -> np.ndarray:
        """Cumulative vertical node depths in metres."""
        return np.concatenate([[0.0], np.cumsum(self.z_widths)])

    @property
    def rho_linear(self) -> np.ndarray:
        """Return resistivity in linear ohm metres."""
        return np.exp(self.rho_loge)

    # ------------------------------------------------------------------
    # Construction
    # ------------------------------------------------------------------
    @classmethod
    def halfspace(
        cls,
        data: ModEmData,
        config: ModEmConfig | None = None,
        **kwargs,
    ) -> ModEmModel2D:
        """Build a uniform half-space starting model.

        The method derives the horizontal grid from station
        offsets in a :class:`ModEmData` object. It fills earth
        cells with ``config.initial_rho`` and assigns air layers
        a high fixed resistivity. Resistivity is stored
        internally as natural logarithm values.

        Parameters
        ----------
        data : ModEmData
            Populated data object. The builder uses
            ``data.offsets`` to determine the station-zone
            width and to place lateral padding cells.
        config : ModEmConfig, optional
            Configuration supplying 2-D grid dimensions,
            padding, layer growth, air-layer count, and
            starting half-space resistivity. If omitted, a
            default :class:`ModEmConfig` is used.
        **kwargs : dict
            Additional keyword arguments forwarded to
            :class:`ModEmModel2D`.

        Returns
        -------
        ModEmModel2D
            Populated 2-D model object ready to be written as a
            ModEM model file.

        Examples
        --------
        >>> from pycsamt.models.modem.model2d import ModEmModel2D
        >>> model = ModEmModel2D.halfspace(data, config=cfg)
        >>> model.rho_loge.shape == (model.nz, model.nx)
        True

        Notes
        -----
        The station-zone cell widths are derived from gaps
        between sorted station offsets. If only one station is
        present, one station-zone cell is created using
        ``config.cell_size_h_2d``. Padding cells grow by powers
        of two away from the station zone.
        """
        cfg = config or ModEmConfig()
        obj = cls(config=cfg, **kwargs)

        offsets = np.sort(data.offsets)
        n_pad = cfg.n_padding_x_2d

        # station-zone cells
        cell_h = cfg.cell_size_h_2d
        station_widths: list[float] = []
        if offsets.size <= 1:
            station_widths = [cell_h]
        else:
            for i in range(len(offsets) - 1):
                gap = float(offsets[i + 1] - offsets[i])
                n_cell = max(1, round(gap / cell_h))
                station_widths.extend([gap / n_cell] * n_cell)
        if len(station_widths) % 2 != 0:
            station_widths.append(
                station_widths[-1] if station_widths else cell_h
            )

        # padding
        pad = [cell_h * float(2 ** (k + 1)) for k in range(n_pad)]
        x_w = np.array(
            list(reversed(pad)) + station_widths + pad,
            dtype=float,
        )

        # vertical
        n_air = cfg.n_airlayers_2d
        n_active = cfg.nz_2d
        air_h = cfg.cell_size_v_top_2d
        air_z = [air_h] * n_air
        z_w: list[float] = []
        thick = float(cfg.cell_size_v_top_2d)
        for _ in range(n_active):
            z_w.append(thick)
            thick *= cfg.depth_scale_2d
        z_w_arr = np.array(air_z + z_w, dtype=float)

        nz_total = len(z_w_arr)
        nx_total = len(x_w)
        rho_val = float(np.log(cfg.initial_rho))
        rho_grid = np.full((nz_total, nx_total), rho_val, dtype=float)
        # Air layers set to very high resistivity
        if n_air > 0:
            rho_grid[:n_air, :] = np.log(1e12)

        obj.x_widths = x_w
        obj.z_widths = z_w_arr
        obj.rho_loge = rho_grid

        if obj.verbose:
            obj.logger.info(
                "ModEmModel2D.halfspace: %d x %d grid, rho=%.1f ohm m",
                nx_total,
                nz_total,
                cfg.initial_rho,
            )
        return obj

    # ------------------------------------------------------------------
    # I/O
    # ------------------------------------------------------------------
    @classmethod
    def read(cls, path: PathLike, **kwargs) -> ModEmModel2D:
        """Parse an existing ModEM 2-D model file.

        Parameters
        ----------
        path : path-like
            Path to a ModEM 2-D model file. The file may store
            resistivity as ``LOGE``, ``LOG10``, or ``LINEAR``;
            values are converted to natural-log resistivity
            internally.
        **kwargs : dict
            Additional keyword arguments forwarded to
            :class:`ModEmModel2D`.

        Returns
        -------
        ModEmModel2D
            Parsed model with cell widths, layer thicknesses,
            natural-log resistivity, and log-type metadata.

        Raises
        ------
        FileNotFoundError
            If ``path`` does not exist.

        Examples
        --------
        >>> from pycsamt.models.modem.model2d import ModEmModel2D
        >>> model = ModEmModel2D.read("m0.rho")
        >>> model.rho_loge.shape == (model.nz, model.nx)
        True
        """
        p = Path(path)
        if not p.exists():
            msg = f"ModEM 2D model file not found: {p}"
            raise FileNotFoundError(msg)
        d = _parse_model2d(p)
        obj = cls(**kwargs)
        obj.x_widths = d["x_widths"]
        obj.z_widths = d["z_widths"]
        obj.rho_loge = d["rho_loge"]
        obj.log_type = d["log_type"]
        if obj.verbose:
            obj.logger.info(
                "ModEmModel2D.read: %d x %d from %s",
                obj.nx,
                obj.nz,
                p,
            )
        return obj

    def write(self, path: PathLike) -> Path:
        """Write the model to ``path`` in ModEM 2-D format.

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
        The writer emits one parameter block and writes the
        values stored in ``rho_loge``. The ``log_type`` header
        is written from the object, so callers should keep it
        consistent with the numerical encoding of ``rho_loge``.

        Examples
        --------
        >>> from pycsamt.models.modem.model2d import ModEmModel2D
        >>> model = ModEmModel2D.halfspace(data, config=cfg)
        >>> path = model.write("m0.rho")
        >>> path.name
        'm0.rho'
        """
        p = Path(path)
        p.parent.mkdir(parents=True, exist_ok=True)

        def _floats(arr: np.ndarray, per_row: int = 12) -> str:
            rows = []
            for i in range(0, len(arr), per_row):
                rows.append(
                    "  "
                    + "  ".join(f"{v:>10.4f}" for v in arr[i : i + per_row])
                )
            return "\n".join(rows) + "\n"

        lines: list[str] = [
            f"  {self.nx}  {self.nz}  {self.log_type}\n",
        ]
        lines.append(_floats(self.x_widths))
        lines.append(_floats(self.z_widths))
        lines.append("  1\n")  # block count

        for iz in range(self.nz):
            lines.append(
                "  "
                + "  ".join(f"{v:>12.5E}" for v in self.rho_loge[iz, :])
                + "\n"
            )

        with p.open("w") as fh:
            fh.writelines(lines)
        return p


ModEmModel2D.__doc__ = rf"""
Represent a ModEM two-dimensional resistivity model.

``ModEmModel2D`` stores the model grid and resistivity values
used by the ModEM 2-D executable. The horizontal axis follows
the survey profile, while the vertical axis contains air layers
and active earth layers. Resistivity values are stored
internally as natural logarithms, so a linear resistivity
:math:`\rho` is represented as

.. math::

   m = \ln(\rho).

The ASCII file format contains a header with ``nx``, ``nz``,
and log encoding, followed by horizontal cell widths, vertical
layer thicknesses, a block count, and the resistivity grid.

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
    Horizontal cell widths in metres. The array includes
    lateral padding and station-zone cells.
z_widths : numpy.ndarray, shape (nz,)
    Vertical layer thicknesses in metres. The array includes
    air layers followed by active earth layers.
rho_loge : numpy.ndarray, shape (nz, nx)
    Natural-log resistivity values. Earth cells in a default
    half-space are initialized as
    :math:`\ln(\rho_0)`, where :math:`\rho_0` is
    ``config.initial_rho``.
log_type : str, default "LOGE"
    Encoding label written to the ModEM file header. ``"LOGE"``
    is the standard internal representation used by this
    class. Readers also accept ``"LOG10"`` and ``"LINEAR"`` and
    convert them to natural logarithms.

Derived Properties
------------------
nx : int
    Number of horizontal cells.
nz : int
    Number of vertical layers, including air layers.
x_nodes : numpy.ndarray, shape (nx + 1,)
    Cumulative horizontal node coordinates in metres.
z_nodes : numpy.ndarray, shape (nz + 1,)
    Cumulative vertical node depths in metres.
rho_linear : numpy.ndarray, shape (nz, nx)
    Resistivity values in linear ohm metres, obtained as
    :math:`\exp(\mathtt{{rho\_loge}})`.

Notes
-----
The :meth:`halfspace` constructor creates a conservative
starting model. Air layers are assigned a very high
resistivity, while earth cells are assigned
``config.initial_rho``. Horizontal padding grows away from the
station zone so artificial boundaries are moved away from the
profile.

See Also
--------
ModEmConfig
    Supplies 2-D grid and initial-resistivity settings.
ModEmData
    Provides station offsets used by :meth:`halfspace`.
InputBuilder
    Builds and writes a 2-D starting model automatically.
ModEmRunner
    Passes the written model file to the ModEM executable.

Examples
--------
Build a 2-D half-space model:

>>> from pycsamt.models.modem.config import ModEmConfig
>>> from pycsamt.models.modem.model2d import ModEmModel2D
>>> cfg = ModEmConfig(mode="2d", initial_rho=100.0)
>>> model = ModEmModel2D.halfspace(data, config=cfg)
>>> model.rho_loge.shape == (model.nz, model.nx)
True

Read and write a ModEM model file:

>>> model = ModEmModel2D.read("m0.rho")
>>> path = model.write("m0_copy.rho")
>>> path.name
'm0_copy.rho'

Inspect linear resistivity:

>>> rho = model.rho_linear
>>> rho.shape == model.rho_loge.shape
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
