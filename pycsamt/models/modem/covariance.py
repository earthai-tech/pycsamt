# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""ModEmCovariance — read and write the ModEM covariance file.

The covariance file controls smoothing in the ModEM inversion.

File format
-----------
::

    +-- 16-line comment header --+

    Nx  Ny  NzEarth           ← grid dims *excluding* air layers

    smooth_x_1 ... smooth_x_NzEarth    ← per-layer smoothing in X
    smooth_y_1 ... smooth_y_NzEarth    ← per-layer smoothing in Y
    smooth_z                           ← scalar vertical smoothing

    n_smooth_iter

    n_exceptions
    mask_a  mask_b  value    ← one line per exception

    layer_start  layer_end   ← mask block header (1-based layer indices)
    <Ny rows of Nx integer mask values>
    [additional blocks …]

Mask values: 0 = air (smoothing disabled automatically),
9 = ocean, 1–8 user-defined regions.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import List, Optional, Tuple, Union

import numpy as np

from .base   import ModEmBase
from .config import ModEmConfig

PathLike = Union[str, Path]

__all__ = ["ModEmCovariance"]

_HEADER_LINES = 16

# ---------------------------------------------------------------------------
# Low-level parser
# ---------------------------------------------------------------------------

def _parse_covariance(path: Path) -> dict:
    with path.open("r", errors="replace") as fh:
        lines = [ln.rstrip("\n") for ln in fh]

    N = len(lines)
    i = 0

    # -- 16-line comment header ------------------------------------------
    header_comment: list[str] = []
    while i < N and len(header_comment) < _HEADER_LINES:
        header_comment.append(lines[i])
        i += 1

    # -- skip blanks, then read "Nx Ny NzEarth" --------------------------
    while i < N and not lines[i].strip():
        i += 1
    dims = lines[i].split()
    nx_earth = int(dims[0])
    ny_earth = int(dims[1])
    nz_earth = int(dims[2])
    i += 1

    # -- collect smoothing floats: 2*NzEarth + 1 values ------------------
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

    smooth_x = np.array(smooth_vals[:nz_earth],                dtype=float)
    smooth_y = np.array(smooth_vals[nz_earth:2 * nz_earth],    dtype=float)
    smooth_z = float(smooth_vals[2 * nz_earth]) if len(smooth_vals) > 2 * nz_earth else 0.1

    # -- n_smooth_iter (first non-blank integer) -------------------------
    while i < N and not lines[i].strip():
        i += 1
    n_smooth_iter = int(lines[i].strip())
    i += 1

    # -- n_exceptions ----------------------------------------------------
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

    # -- mask blocks: (layer_start, layer_end, Ny×Nx grid) ---------------
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
                l_end   = int(parts[1])
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
                    int(row_parts[0]); int(row_parts[1])
                    i -= 1   # put back
                    break
                except ValueError:
                    pass
            try:
                mask_rows.append([int(v) for v in row_parts[:ny_earth]])
            except ValueError:
                pass

        if mask_rows:
            mask_blocks.append({
                "layer_start": l_start,
                "layer_end":   l_end,
                "mask":        np.array(mask_rows, dtype=np.int32),
            })

    return {
        "nx_earth":     nx_earth,
        "ny_earth":     ny_earth,
        "nz_earth":     nz_earth,
        "smooth_x":     smooth_x,
        "smooth_y":     smooth_y,
        "smooth_z":     smooth_z,
        "n_smooth_iter": n_smooth_iter,
        "exceptions":   exceptions,
        "mask_blocks":  mask_blocks,
    }


# ---------------------------------------------------------------------------
# ModEmCovariance
# ---------------------------------------------------------------------------

_HEADER_TEMPLATE = """\
+-----------------------------------------------------------------------------+
| This file defines model covariance for a recursive autoregression scheme.   |
| The model space may be divided into distinct areas using integer masks.     |
| Mask 0 is reserved for air; mask 9 is reserved for ocean. Smoothing between |
| air, ocean and the rest of the model is turned off automatically. You can   |
| also define exceptions to override smoothing between any two model areas.   |
| To turn off smoothing set it to zero. This header is 16 lines long.         |
| 1. Grid dimensions excluding air layers (Nx, Ny, NzEarth)                   |
| 2. Smoothing in the X direction (NzEarth real values)                       |
| 3. Smoothing in the Y direction (NzEarth real values)                       |
| 4. Vertical smoothing (1 real value)                                        |
| 5. Number of times the smoothing should be applied (1 integer >= 0)         |
| 6. Number of exceptions (1 integer >= 0)                                    |
| 7. Exceptions in the form e.g. 2 3 0. (to turn off smoothing between 3 & 4) |
| 8. Two integer layer indices and Nx x Ny block of masks, repeated as needed.|
+-----------------------------------------------------------------------------+"""


class ModEmCovariance(ModEmBase):
    """ModEM covariance / smoothing file container.

    Attributes
    ----------
    nx_earth, ny_earth, nz_earth : int
        Grid dimensions (excluding air layers).
    smooth_x : np.ndarray, shape (nz_earth,)
        Per-layer smoothing coefficient in X.
    smooth_y : np.ndarray, shape (nz_earth,)
        Per-layer smoothing coefficient in Y.
    smooth_z : float
        Vertical smoothing coefficient.
    n_smooth_iter : int
        Number of smoothing iterations.
    exceptions : list of (int, int, float)
        Pairs of mask IDs and the smoothing value between them.
    mask_blocks : list of dict
        Each entry has ``'layer_start'``, ``'layer_end'``, and
        ``'mask'`` (np.ndarray of shape (nx_earth, ny_earth)).
    """

    def __init__(self, config: Optional[ModEmConfig] = None, **kwargs):
        super().__init__(**kwargs)
        cfg = config or ModEmConfig()
        self.config:        ModEmConfig     = cfg
        self.nx_earth:      int             = 0
        self.ny_earth:      int             = 0
        self.nz_earth:      int             = 0
        self.smooth_x:      np.ndarray      = np.array([])
        self.smooth_y:      np.ndarray      = np.array([])
        self.smooth_z:      float           = cfg.smooth_z
        self.n_smooth_iter: int             = cfg.n_smooth_iter
        self.exceptions:    List[Tuple]     = []
        self.mask_blocks:   List[dict]      = []

    # ------------------------------------------------------------------
    # Construction
    # ------------------------------------------------------------------
    @classmethod
    def from_model(
        cls,
        model: "ModEmModel3D",    # noqa: F821
        config: Optional[ModEmConfig] = None,
        **kwargs,
    ) -> "ModEmCovariance":
        """Build a uniform covariance object from a ``ModEmModel3D``.

        All earth cells get mask value 1 (one region, uniform smoothing).
        Air layers (mask 0) are handled automatically by ModEM.
        """
        cfg = config or ModEmConfig()
        obj = cls(config=cfg, **kwargs)

        n_air     = model.n_air
        nz_earth  = model.nz - n_air
        nx_earth  = model.nx
        ny_earth  = model.ny

        obj.nx_earth      = nx_earth
        obj.ny_earth      = ny_earth
        obj.nz_earth      = nz_earth
        obj.smooth_x      = np.full(nz_earth, cfg.smooth_x)
        obj.smooth_y      = np.full(nz_earth, cfg.smooth_y)
        obj.smooth_z      = cfg.smooth_z
        obj.n_smooth_iter = cfg.n_smooth_iter
        obj.exceptions    = []
        obj.mask_blocks   = [{
            "layer_start": 1,
            "layer_end":   nz_earth,
            "mask":        np.ones((nx_earth, ny_earth), dtype=np.int32),
        }]
        return obj

    # ------------------------------------------------------------------
    # I/O
    # ------------------------------------------------------------------
    @classmethod
    def read(cls, path: PathLike, **kwargs) -> "ModEmCovariance":
        """Parse an existing ModEM covariance file."""
        p = Path(path)
        if not p.exists():
            raise FileNotFoundError(f"ModEM covariance file not found: {p}")

        parsed = _parse_covariance(p)
        obj    = cls(**kwargs)
        obj.nx_earth      = parsed["nx_earth"]
        obj.ny_earth      = parsed["ny_earth"]
        obj.nz_earth      = parsed["nz_earth"]
        obj.smooth_x      = parsed["smooth_x"]
        obj.smooth_y      = parsed["smooth_y"]
        obj.smooth_z      = parsed["smooth_z"]
        obj.n_smooth_iter = parsed["n_smooth_iter"]
        obj.exceptions    = parsed["exceptions"]
        obj.mask_blocks   = parsed["mask_blocks"]

        if obj.verbose:
            obj.logger.info("ModEmCovariance.read: loaded from %s", p)
        return obj

    def write(self, path: PathLike) -> Path:
        """Write covariance file to *path*.

        Returns the resolved path.
        """
        p = Path(path)
        p.parent.mkdir(parents=True, exist_ok=True)

        def _float_row(arr: np.ndarray, per_row: int = 20) -> str:
            parts = []
            for i in range(0, len(arr), per_row):
                parts.append(" ".join(f"{v}" for v in arr[i:i+per_row]))
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
