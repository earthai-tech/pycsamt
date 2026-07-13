# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Read, write, and build ModEM covariance files.

The module defines :class:`ModEmCovariance`, the container for
3-D ModEM smoothing coefficients, smoothing exceptions, and
integer region masks.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Union

import numpy as np

from .base import ModEmBase
from .config import ModEmConfig
from .doc import _modem_param_docs as _params

if TYPE_CHECKING:
    from .model3d import ModEmModel3D

PathLike = Union[str, Path]

__all__ = ["ModEmCovariance"]

_HEADER_LINES = 16

# ----------------------------------------------------------------------
# Low-level parser
# ----------------------------------------------------------------------


def _parse_covariance(path: Path) -> dict:
    """Parse a ModEM covariance file into plain Python values."""
    with path.open("r", errors="replace") as fh:
        lines = [ln.rstrip("\n") for ln in fh]

    N = len(lines)
    i = 0

    # -- 16-line comment header ----------------------------------------
    header_comment: list[str] = []
    while i < N and len(header_comment) < _HEADER_LINES:
        header_comment.append(lines[i])
        i += 1

    # -- skip blanks, then read "Nx Ny NzEarth" ------------------------
    while i < N and not lines[i].strip():
        i += 1
    dims = lines[i].split()
    nx_earth = int(dims[0])
    ny_earth = int(dims[1])
    nz_earth = int(dims[2])
    i += 1

    # -- collect smoothing floats: 2*NzEarth + 1 values ----------------
    n_smooth_vals = 2 * nz_earth + 1
    smooth_vals: list[float] = []
    while i < N and len(smooth_vals) < n_smooth_vals:
        ln = lines[i].strip()
        i += 1
        if not ln:
            continue
        for tok in ln.split():
            try:
                smooth_vals.append(float(tok))
                if len(smooth_vals) == n_smooth_vals:
                    break
            except ValueError:
                pass

    smooth_x = np.array(smooth_vals[:nz_earth], dtype=float)
    smooth_y = np.array(
        smooth_vals[nz_earth : 2 * nz_earth],
        dtype=float,
    )
    if len(smooth_vals) > 2 * nz_earth:
        smooth_z = float(smooth_vals[2 * nz_earth])
    else:
        smooth_z = 0.1

    # -- n_smooth_iter (first non-blank integer) -----------------------
    while i < N and not lines[i].strip():
        i += 1
    n_smooth_iter = int(lines[i].strip())
    i += 1

    # -- n_exceptions --------------------------------------------------
    while i < N and not lines[i].strip():
        i += 1
    n_exceptions = int(lines[i].strip())
    i += 1

    exceptions: list[tuple] = []
    for _ in range(n_exceptions):
        while i < N and not lines[i].strip():
            i += 1
        parts = lines[i].split()
        i += 1
        exceptions.append((int(parts[0]), int(parts[1]), float(parts[2])))

    # -- mask blocks: (layer_start, layer_end, Ny by Nx grid) ----------
    mask_blocks: list[dict] = []
    while i < N:
        # find a line that looks like two integers (layer range)
        while i < N and not lines[i].strip():
            i += 1
        if i >= N:
            break
        parts = lines[i].strip().split()
        if len(parts) >= 2:
            try:
                l_start = int(parts[0])
                l_end = int(parts[1])
                i += 1
            except ValueError:
                i += 1
                continue
        else:
            i += 1
            continue

        # read Nx rows of Ny integers  (mask indexed as [ix, iy])
        mask_rows: list[list[int]] = []
        while i < N and len(mask_rows) < nx_earth:
            ln = lines[i].strip()
            i += 1
            if not ln:
                continue
            row_parts = ln.split()
            # stop if this looks like a new block header (2 ints only)
            if len(row_parts) == 2:
                try:
                    int(row_parts[0])
                    int(row_parts[1])
                    i -= 1  # put back
                    break
                except ValueError:
                    pass
            try:
                mask_rows.append([int(v) for v in row_parts[:ny_earth]])
            except ValueError:
                pass

        if mask_rows:
            mask_blocks.append(
                {
                    "layer_start": l_start,
                    "layer_end": l_end,
                    "mask": np.array(mask_rows, dtype=np.int32),
                }
            )

    return {
        "nx_earth": nx_earth,
        "ny_earth": ny_earth,
        "nz_earth": nz_earth,
        "smooth_x": smooth_x,
        "smooth_y": smooth_y,
        "smooth_z": smooth_z,
        "n_smooth_iter": n_smooth_iter,
        "exceptions": exceptions,
        "mask_blocks": mask_blocks,
    }


# ----------------------------------------------------------------------
# ModEmCovariance
# ----------------------------------------------------------------------

_HEADER_TEMPLATE = """\
+---------------------------------------------------------------+
| ModEM model covariance for recursive autoregression smoothing.|
| Integer masks divide the model into smoothing regions.        |
| Mask 0 is air; mask 9 is ocean. Smoothing involving air or    |
| ocean is turned off automatically by ModEM.                   |
| Exceptions can override smoothing between two mask regions.   |
| Set an exception value to zero to turn smoothing off.         |
| 1. Grid dimensions excluding air layers: Nx Ny NzEarth.       |
| 2. Smoothing in X: NzEarth real values.                       |
| 3. Smoothing in Y: NzEarth real values.                       |
| 4. Vertical smoothing: one real value.                        |
| 5. Number of smoothing applications: one integer >= 0.        |
| 6. Number of exceptions: one integer >= 0.                    |
| 7. Exceptions: mask_a mask_b value.                           |
| 8. Layer range and Nx by Ny mask block, repeated as needed.   |
+---------------------------------------------------------------+"""


class ModEmCovariance(ModEmBase):
    def __init__(self, config: ModEmConfig | None = None, **kwargs):
        """Initialize an empty covariance object."""
        super().__init__(**kwargs)
        cfg = config or ModEmConfig()
        self.config: ModEmConfig = cfg
        self.nx_earth: int = 0
        self.ny_earth: int = 0
        self.nz_earth: int = 0
        self.smooth_x: np.ndarray = np.array([])
        self.smooth_y: np.ndarray = np.array([])
        self.smooth_z: float = cfg.smooth_z
        self.n_smooth_iter: int = cfg.n_smooth_iter
        self.exceptions: list[tuple[int, int, float]] = []
        self.mask_blocks: list[dict] = []

    # ------------------------------------------------------------------
    # Construction
    # ------------------------------------------------------------------
    @classmethod
    def from_model(
        cls,
        model: ModEmModel3D,
        config: ModEmConfig | None = None,
        **kwargs,
    ) -> ModEmCovariance:
        """Build a uniform covariance object from a 3-D model.

        The generated object assigns one active earth region
        with mask value ``1`` over all earth layers and uses
        the smoothing values stored in ``config``. Air layers
        are excluded from the covariance grid because the ModEM
        covariance dimensions are earth-only.

        Parameters
        ----------
        model : ModEmModel3D
            Populated 3-D model used to determine earth-grid
            dimensions. The model must expose ``nx``, ``ny``,
            ``nz``, and ``n_air``. The number of covariance
            layers is computed as ``model.nz - model.n_air``.
        config : ModEmConfig, optional
            Configuration supplying ``smooth_x``, ``smooth_y``,
            ``smooth_z``, and ``n_smooth_iter``. If omitted, a
            default :class:`ModEmConfig` is used.
        **kwargs : dict
            Additional keyword arguments forwarded to
            :class:`ModEmCovariance`, commonly ``verbose`` or
            ``logger``.

        Returns
        -------
        ModEmCovariance
            Covariance object with one mask block spanning all
            earth layers.

        Examples
        --------
        >>> from pycsamt.models.modem.covariance import ModEmCovariance
        >>> cov = ModEmCovariance.from_model(model, config=cfg)
        >>> cov.nz_earth == model.nz - model.n_air
        True
        >>> len(cov.mask_blocks)
        1

        Notes
        -----
        ``from_model`` creates a simple uniform regularization
        region. Users who need geological domains, ocean masks,
        or smoothing exceptions can edit ``mask_blocks`` and
        ``exceptions`` before calling :meth:`write`.
        """
        cfg = config or ModEmConfig()
        obj = cls(config=cfg, **kwargs)

        n_air = model.n_air
        nz_earth = model.nz - n_air
        nx_earth = model.nx
        ny_earth = model.ny

        obj.nx_earth = nx_earth
        obj.ny_earth = ny_earth
        obj.nz_earth = nz_earth
        obj.smooth_x = np.full(nz_earth, cfg.smooth_x)
        obj.smooth_y = np.full(nz_earth, cfg.smooth_y)
        obj.smooth_z = cfg.smooth_z
        obj.n_smooth_iter = cfg.n_smooth_iter
        obj.exceptions = []
        obj.mask_blocks = [
            {
                "layer_start": 1,
                "layer_end": nz_earth,
                "mask": np.ones((nx_earth, ny_earth), dtype=np.int32),
            }
        ]
        return obj

    # ------------------------------------------------------------------
    # I/O
    # ------------------------------------------------------------------
    @classmethod
    def read(cls, path: PathLike, **kwargs) -> ModEmCovariance:
        """Parse an existing ModEM covariance file.

        Parameters
        ----------
        path : path-like
            Path to a ModEM covariance file. The file must
            contain the standard header, earth-grid dimensions,
            smoothing rows, exception count, and one or more
            mask blocks.
        **kwargs : dict
            Additional keyword arguments forwarded to
            :class:`ModEmCovariance`.

        Returns
        -------
        ModEmCovariance
            Parsed covariance object with smoothing arrays,
            exception tuples, and mask blocks populated.

        Raises
        ------
        FileNotFoundError
            If ``path`` does not exist.

        Examples
        --------
        >>> from pycsamt.models.modem.covariance import ModEmCovariance
        >>> cov = ModEmCovariance.read("ModEM.cov")
        >>> cov.smooth_x.shape == (cov.nz_earth,)
        True
        """
        p = Path(path)
        if not p.exists():
            msg = f"ModEM covariance file not found: {p}"
            raise FileNotFoundError(msg)

        parsed = _parse_covariance(p)
        obj = cls(**kwargs)
        obj.nx_earth = parsed["nx_earth"]
        obj.ny_earth = parsed["ny_earth"]
        obj.nz_earth = parsed["nz_earth"]
        obj.smooth_x = parsed["smooth_x"]
        obj.smooth_y = parsed["smooth_y"]
        obj.smooth_z = parsed["smooth_z"]
        obj.n_smooth_iter = parsed["n_smooth_iter"]
        obj.exceptions = parsed["exceptions"]
        obj.mask_blocks = parsed["mask_blocks"]

        if obj.verbose:
            obj.logger.info("ModEmCovariance.read: loaded from %s", p)
        return obj

    def write(self, path: PathLike) -> Path:
        """Write the covariance object to a ModEM file.

        Parameters
        ----------
        path : path-like
            Destination covariance file. Parent directories are
            created before writing. Existing files are
            overwritten.

        Returns
        -------
        pathlib.Path
            Path passed to the writer, converted to
            :class:`pathlib.Path`.

        Raises
        ------
        KeyError
            If a mask block lacks ``"layer_start"``,
            ``"layer_end"``, or ``"mask"``.

        Examples
        --------
        >>> from pycsamt.models.modem.covariance import ModEmCovariance
        >>> cov = ModEmCovariance.from_model(model, config=cfg)
        >>> path = cov.write("covariance.cov")
        >>> path.name
        'covariance.cov'

        Notes
        -----
        The writer preserves the mask-block orientation used by
        this module: each mask array is shaped
        ``(nx_earth, ny_earth)`` and is written row by row.
        """
        p = Path(path)
        p.parent.mkdir(parents=True, exist_ok=True)

        def _float_row(arr: np.ndarray, per_row: int = 20) -> str:
            parts = []
            for i in range(0, len(arr), per_row):
                parts.append(" ".join(f"{v}" for v in arr[i : i + per_row]))
            return "\n".join(parts)

        lines: list[str] = [_HEADER_TEMPLATE + "\n"]
        lines.append("\n")
        lines.append(f"{self.nx_earth} {self.ny_earth} {self.nz_earth}\n")
        lines.append("\n")
        lines.append(_float_row(self.smooth_x) + "\n")
        lines.append(_float_row(self.smooth_y) + "\n")
        lines.append(f"{self.smooth_z}\n")
        lines.append("\n")
        lines.append(f"{self.n_smooth_iter}\n")
        lines.append("\n")
        lines.append(f"{len(self.exceptions)}\n")
        for exc in self.exceptions:
            lines.append(f"{exc[0]} {exc[1]} {exc[2]}\n")
        lines.append("\n")

        for blk in self.mask_blocks:
            lines.append(f"{blk['layer_start']} {blk['layer_end']}\n")
            mask = blk["mask"]
            for row in mask:
                lines.append(" ".join(str(v) for v in row) + "\n")
            lines.append("\n")

        with p.open("w") as fh:
            fh.writelines(lines)
        return p


ModEmCovariance.__doc__ = rf"""
Represent a ModEM 3-D covariance and smoothing file.

``ModEmCovariance`` stores the regularization information used
by ModEM 3-D inversions. The covariance file defines earth-only
grid dimensions, per-layer horizontal smoothing coefficients,
one vertical smoothing coefficient, optional smoothing
exceptions between mask regions, and integer mask blocks that
divide the model into smoothing domains.

In ModEM's regularized objective, covariance controls the model
regularization term:

.. math::

   \Phi_m(m) =
   \| W_m (m - m_0) \|_2^2 ,

where :math:`W_m` is shaped by smoothing weights and region
masks from this file. Larger smoothing values generally favour
models that vary more gradually in the corresponding direction.
Mask exceptions can reduce or disable smoothing across
geological boundaries, air, or ocean regions [1]_, [2]_.

Parameters
----------
{_params.common.config}
{_params.common.verbose}
{_params.common.logger}

Attributes
----------
config : ModEmConfig
    Configuration object used to initialize smoothing defaults.
nx_earth : int
    Number of earth cells along the ModEM x direction. Air
    layers are excluded from this dimension.
ny_earth : int
    Number of earth cells along the ModEM y direction. Air
    layers are excluded from this dimension.
nz_earth : int
    Number of earth layers represented by the covariance file.
    This is normally ``model.nz - model.n_air``.
smooth_x : numpy.ndarray, shape (nz_earth,)
    Per-layer smoothing coefficients in the x direction.
smooth_y : numpy.ndarray, shape (nz_earth,)
    Per-layer smoothing coefficients in the y direction.
{_params.covariance.smooth_z}
{_params.covariance.n_smooth_iter}
exceptions : list of tuple of int, int, float
    Smoothing overrides between two mask identifiers. A tuple
    ``(a, b, value)`` sets the smoothing across the boundary
    between mask regions ``a`` and ``b``. Values of ``0`` turn
    off smoothing across that boundary.
mask_blocks : list of dict
    Layer-group masks. Each block contains ``"layer_start"``,
    ``"layer_end"``, and ``"mask"``. Layer indices are
    one-based and inclusive, matching the ModEM text format.
    Each mask array has shape ``(nx_earth, ny_earth)``.

Notes
-----
The covariance file begins with a 16-line explanatory header,
then stores the earth-grid dimensions, smoothing arrays,
vertical smoothing value, number of smoothing iterations,
exception rows, and one or more mask blocks.

Mask values follow the ModEM convention:

* ``0`` is reserved for air;
* ``9`` is reserved for ocean;
* ``1`` through ``8`` are user-defined model regions.

Smoothing involving air and ocean is disabled by ModEM
automatically. Additional boundaries can be controlled through
``exceptions``.

See Also
--------
ModEmConfig
    Supplies default smoothing and iteration values.
ModEmModel3D
    Provides the grid dimensions used by :meth:`from_model`.
InputBuilder
    Creates a covariance file automatically for 3-D builds.
ModEmRunner
    Passes the covariance file to the ModEM 3-D executable.

Examples
--------
Create a uniform covariance object from a 3-D model:

>>> from pycsamt.models.modem.covariance import ModEmCovariance
>>> cov = ModEmCovariance.from_model(model, config=cfg)
>>> cov.nz_earth == model.nz - model.n_air
True
>>> cov.mask_blocks[0]["layer_start"]
1

Disable smoothing between two user-defined regions:

>>> cov.exceptions.append((1, 2, 0.0))
>>> cov.write("ModEM.cov")
PosixPath('ModEM.cov')

Read an existing covariance file:

>>> loaded = ModEmCovariance.read("ModEM.cov")
>>> loaded.n_smooth_iter >= 0
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
