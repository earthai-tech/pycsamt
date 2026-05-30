# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""ModEmModel2D — read, write, and build the ModEM 2-D model file.

File format
-----------
::

    N_x  N_z  LOGE          ← header: n horizontal cells, n vertical cells, log type
    <x-cell widths: N_x floats, wrapped across lines>
    <z-layer thicknesses: N_z floats, wrapped across lines>
    1                        ← number of model parameter blocks (always 1)
    <N_z rows of N_x loge(ρ) values>

The resistivity values are stored as ``ln(ρ)`` (natural log, LOGE) or
optionally ``log10(ρ)`` (LOG10) or linear (LINEAR).  ``LOGE`` is the
standard ModEM format.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional, Union

import numpy as np

from .base   import ModEmBase
from .config import ModEmConfig

PathLike = Union[str, Path]

__all__ = ["ModEmModel2D"]


# ---------------------------------------------------------------------------
# Low-level parser
# ---------------------------------------------------------------------------

def _parse_model2d(path: Path) -> dict:
    with path.open("r", errors="replace") as fh:
        lines = [l.rstrip("\n") for l in fh]

    # Skip blank / comment lines to find control line
    i = 0
    N = len(lines)
    while i < N and (not lines[i].strip() or lines[i].strip().startswith("#")):
        i += 1

    ctrl   = lines[i].split()
    i     += 1
    nx     = int(ctrl[0])
    nz     = int(ctrl[1])
    log_type = ctrl[2].upper() if len(ctrl) > 2 else "LOGE"

    # Collect all floats until we see a lone integer that is "1" (block count)
    float_tokens: list[float] = []
    n_blocks = 1
    while i < N:
        ln = lines[i].strip()
        i += 1
        if not ln or ln.startswith("#"):
            continue
        # Block-count sentinel: a lone integer
        parts = ln.split()
        if len(parts) == 1:
            try:
                n_blocks = int(parts[0])
                break
            except ValueError:
                pass
        for tok in parts:
            try:
                float_tokens.append(float(tok))
            except ValueError:
                pass

    x_widths = np.array(float_tokens[:nx],       dtype=float)
    z_widths = np.array(float_tokens[nx:nx + nz], dtype=float)

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

    rho_log = np.array(rho_flat[:nx * nz], dtype=float).reshape(nz, nx)

    # Convert to loge if needed
    if log_type == "LOG10":
        rho_loge = rho_log * np.log(10.0)
    elif log_type == "LINEAR":
        with np.errstate(divide="ignore", invalid="ignore"):
            rho_loge = np.log(np.where(rho_log > 0, rho_log, np.nan))
    else:
        rho_loge = rho_log   # already LOGE

    return {
        "nx": nx, "nz": nz, "log_type": log_type,
        "x_widths": x_widths,
        "z_widths": z_widths,
        "rho_loge": rho_loge,
    }


# ---------------------------------------------------------------------------
# ModEmModel2D
# ---------------------------------------------------------------------------

class ModEmModel2D(ModEmBase):
    """ModEM 2-D model container.

    Attributes
    ----------
    x_widths : np.ndarray, shape (nx,)
        Horizontal cell widths (m).
    z_widths : np.ndarray, shape (nz,)
        Vertical layer thicknesses (m).
    rho_loge : np.ndarray, shape (nz, nx)
        Natural-log resistivity values ln(Ω·m).
    log_type : str
        ``'LOGE'`` (default).
    """

    def __init__(self, config: Optional[ModEmConfig] = None, **kwargs):
        super().__init__(**kwargs)
        self.config:    ModEmConfig = config or ModEmConfig()
        self.x_widths:  np.ndarray  = np.array([])
        self.z_widths:  np.ndarray  = np.array([])
        self.rho_loge:  np.ndarray  = np.empty((0, 0))
        self.log_type:  str         = "LOGE"

    # ------------------------------------------------------------------
    # Derived
    # ------------------------------------------------------------------
    @property
    def nx(self) -> int:
        return len(self.x_widths)

    @property
    def nz(self) -> int:
        return len(self.z_widths)

    @property
    def x_nodes(self) -> np.ndarray:
        return np.concatenate([[0.0], np.cumsum(self.x_widths)])

    @property
    def z_nodes(self) -> np.ndarray:
        return np.concatenate([[0.0], np.cumsum(self.z_widths)])

    @property
    def rho_linear(self) -> np.ndarray:
        """Resistivity in Ω·m (linear)."""
        return np.exp(self.rho_loge)

    # ------------------------------------------------------------------
    # Construction
    # ------------------------------------------------------------------
    @classmethod
    def halfspace(
        cls,
        data: "ModEmData",       # noqa: F821
        config: Optional[ModEmConfig] = None,
        **kwargs,
    ) -> "ModEmModel2D":
        """Build a uniform half-space starting model.

        Parameters
        ----------
        data : ModEmData
            Populated data object (``data.offsets`` required).
        config : ModEmConfig, optional

        Returns
        -------
        ModEmModel2D
        """
        cfg = config or ModEmConfig()
        obj = cls(config=cfg, **kwargs)

        offsets = np.sort(data.offsets)
        n_pad   = cfg.n_padding_x_2d

        # station-zone cells
        cell_h      = cfg.cell_size_h_2d
        station_widths: list[float] = []
        if offsets.size <= 1:
            station_widths = [cell_h]
        else:
            for i in range(len(offsets) - 1):
                gap    = float(offsets[i + 1] - offsets[i])
                n_cell = max(1, round(gap / cell_h))
                station_widths.extend([gap / n_cell] * n_cell)
        if len(station_widths) % 2 != 0:
            station_widths.append(station_widths[-1] if station_widths else cell_h)

        # padding
        pad   = [cell_h * float(2 ** (k + 1)) for k in range(n_pad)]
        x_w   = np.array(list(reversed(pad)) + station_widths + pad, dtype=float)

        # vertical
        n_air    = cfg.n_airlayers_2d
        n_active = cfg.nz_2d
        air_h    = cfg.cell_size_v_top_2d
        air_z    = [air_h] * n_air
        z_w: list[float] = []
        thick = float(cfg.cell_size_v_top_2d)
        for _ in range(n_active):
            z_w.append(thick)
            thick *= cfg.depth_scale_2d
        z_w_arr = np.array(air_z + z_w, dtype=float)

        nz_total = len(z_w_arr)
        nx_total = len(x_w)
        rho_val  = float(np.log(cfg.initial_rho))
        rho_grid = np.full((nz_total, nx_total), rho_val, dtype=float)
        # Air layers set to very high resistivity
        if n_air > 0:
            rho_grid[:n_air, :] = np.log(1e12)

        obj.x_widths = x_w
        obj.z_widths = z_w_arr
        obj.rho_loge = rho_grid

        if obj.verbose:
            obj.logger.info(
                "ModEmModel2D.halfspace: %d×%d grid, rho=%.1f Ω·m",
                nx_total, nz_total, cfg.initial_rho,
            )
        return obj

    # ------------------------------------------------------------------
    # I/O
    # ------------------------------------------------------------------
    @classmethod
    def read(cls, path: PathLike, **kwargs) -> "ModEmModel2D":
        """Parse an existing ModEM 2D model file."""
        p = Path(path)
        if not p.exists():
            raise FileNotFoundError(f"ModEM 2D model file not found: {p}")
        d   = _parse_model2d(p)
        obj = cls(**kwargs)
        obj.x_widths = d["x_widths"]
        obj.z_widths = d["z_widths"]
        obj.rho_loge = d["rho_loge"]
        obj.log_type = d["log_type"]
        if obj.verbose:
            obj.logger.info(
                "ModEmModel2D.read: %d×%d from %s", obj.nx, obj.nz, p,
            )
        return obj

    def write(self, path: PathLike) -> Path:
        """Write model to *path* in ModEM 2D format.

        Returns the resolved path.
        """
        p = Path(path)
        p.parent.mkdir(parents=True, exist_ok=True)

        def _floats(arr: np.ndarray, per_row: int = 12) -> str:
            rows = []
            for i in range(0, len(arr), per_row):
                rows.append("  " + "  ".join(f"{v:>10.4f}" for v in arr[i:i+per_row]))
            return "\n".join(rows) + "\n"

        lines: list[str] = [
            f"  {self.nx}  {self.nz}  {self.log_type}\n",
        ]
        lines.append(_floats(self.x_widths))
        lines.append(_floats(self.z_widths))
        lines.append("  1\n")   # block count

        for iz in range(self.nz):
            lines.append("  " + "  ".join(
                f"{v:>12.5E}" for v in self.rho_loge[iz, :]
            ) + "\n")

        with p.open("w") as fh:
            fh.writelines(lines)
        return p
