# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""ModEmModel3D — read, write, and build the ModEM 3-D model file.

File format (WinGLink / ModEM WS format)
-----------------------------------------
::

    # optional comment
    Nx  Ny  Nz  [n_air]  LOGE       ← header
    <Nx x-widths, wrapped>
    <Ny y-widths, wrapped>
    <Nz z-widths, wrapped>
    <blank line>
    <Nz layers × Ny rows of Nx loge(ρ) values>
    [cx  cy  cz]                     ← optional grid centre (km)
    [rotation_deg]                   ← optional rotation angle

Resistivity values are stored as ``ln(ρ)`` (LOGE).  LOG10 and LINEAR are also
supported on read, and converted to LOGE internally.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional, Tuple, Union

import numpy as np

from .base   import ModEmBase
from .config import ModEmConfig

PathLike = Union[str, Path]

__all__ = ["ModEmModel3D"]


# ---------------------------------------------------------------------------
# Low-level parser
# ---------------------------------------------------------------------------

def _parse_model3d(path: Path) -> dict:
    with path.open("r", errors="replace") as fh:
        lines = [ln.rstrip("\n") for ln in fh]

    N = len(lines)
    i = 0
    # skip blank / comment lines before the control line
    while i < N and (not lines[i].strip() or lines[i].strip().startswith("#")):
        i += 1

    ctrl  = lines[i].split()
    i    += 1
    nx    = int(ctrl[0])
    ny    = int(ctrl[1])
    nz    = int(ctrl[2])
    n_air = int(ctrl[3]) if len(ctrl) > 3 and ctrl[3].lstrip("-").isdigit() else 0
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

    x_widths = np.array(float_tokens[:nx],           dtype=float)
    y_widths = np.array(float_tokens[nx:nx + ny],     dtype=float)
    z_widths = np.array(float_tokens[nx + ny:target], dtype=float)

    # read rho grid: Nz * Ny rows of Nx values
    total = nx * ny * nz
    rho_flat: list[float] = []
    while i < N and len(rho_flat) < total:
        ln = lines[i].strip()
        i += 1
        if not ln or ln.startswith("#"):
            continue
        # stop at lines that look like "0 0 0" (centre coordinates) or "0" (rotation)
        parts = ln.split()
        if len(parts) <= 3 and all(_is_simple_numeric(p) for p in parts):
            # could be trailing centre/rotation — only stop if we already have
            # most of the data we need
            if len(rho_flat) >= total - nx:
                break
        for tok in parts:
            try:
                rho_flat.append(float(tok))
            except ValueError:
                pass

    rho_arr = np.array(rho_flat[:total], dtype=float).reshape(nz, ny, nx)

    if log_type == "LOG10":
        rho_loge = rho_arr * np.log(10.0)
    elif log_type == "LINEAR":
        with np.errstate(divide="ignore", invalid="ignore"):
            rho_loge = np.log(np.where(rho_arr > 0, rho_arr, np.nan))
    else:
        rho_loge = rho_arr   # already LOGE

    return {
        "nx": nx, "ny": ny, "nz": nz, "n_air": n_air,
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


# ---------------------------------------------------------------------------
# ModEmModel3D
# ---------------------------------------------------------------------------

class ModEmModel3D(ModEmBase):
    """ModEM 3-D model container.

    Attributes
    ----------
    x_widths : np.ndarray, shape (nx,)
        Cell widths in the X (easting) direction (m).
    y_widths : np.ndarray, shape (ny,)
        Cell widths in the Y (northing) direction (m).
    z_widths : np.ndarray, shape (nz,)
        Layer thicknesses (m), top to bottom.
    rho_loge : np.ndarray, shape (nz, ny, nx)
        Natural-log resistivity ln(Ω·m).
    n_air : int
        Number of air layers at the top of the grid (included in nz).
    log_type : str
        ``'LOGE'`` (default).
    """

    def __init__(self, config: Optional[ModEmConfig] = None, **kwargs):
        super().__init__(**kwargs)
        self.config:    ModEmConfig = config or ModEmConfig()
        self.x_widths:  np.ndarray  = np.array([])
        self.y_widths:  np.ndarray  = np.array([])
        self.z_widths:  np.ndarray  = np.array([])
        self.rho_loge:  np.ndarray  = np.empty((0, 0, 0))
        self.n_air:     int         = 0
        self.log_type:  str         = "LOGE"

    # ------------------------------------------------------------------
    # Derived
    # ------------------------------------------------------------------
    @property
    def nx(self) -> int:
        return len(self.x_widths)

    @property
    def ny(self) -> int:
        return len(self.y_widths)

    @property
    def nz(self) -> int:
        return len(self.z_widths)

    @property
    def x_nodes(self) -> np.ndarray:
        return np.concatenate([[0.0], np.cumsum(self.x_widths)])

    @property
    def y_nodes(self) -> np.ndarray:
        return np.concatenate([[0.0], np.cumsum(self.y_widths)])

    @property
    def z_nodes(self) -> np.ndarray:
        return np.concatenate([[0.0], np.cumsum(self.z_widths)])

    @property
    def rho_linear(self) -> np.ndarray:
        """Resistivity in Ω·m (linear)."""
        return np.exp(self.rho_loge)

    @property
    def shape(self) -> Tuple[int, int, int]:
        return (self.nz, self.ny, self.nx)

    # ------------------------------------------------------------------
    # Construction
    # ------------------------------------------------------------------
    @classmethod
    def halfspace(
        cls,
        data: "ModEmData",       # noqa: F821
        config: Optional[ModEmConfig] = None,
        **kwargs,
    ) -> "ModEmModel3D":
        """Build a uniform half-space starting model.

        Parameters
        ----------
        data : ModEmData
            Populated data object (``data.offsets`` and site coordinates required).
        config : ModEmConfig, optional

        Returns
        -------
        ModEmModel3D
        """
        cfg = config or ModEmConfig()
        obj = cls(config=cfg, **kwargs)

        # -- horizontal grids ------------------------------------------
        xs = np.sort(np.unique(data.x_coords)) if hasattr(data, "x_coords") else np.array([0.0])
        ys = np.sort(np.unique(data.y_coords)) if hasattr(data, "y_coords") else np.array([0.0])

        cell_h = cfg.cell_size_h
        n_pad  = cfg.n_padding_xy

        def _station_cells(coords: np.ndarray) -> list:
            if coords.size <= 1:
                return [cell_h]
            widths = []
            for i in range(len(coords) - 1):
                gap    = float(coords[i + 1] - coords[i])
                n_cell = max(1, round(gap / cell_h))
                widths.extend([gap / n_cell] * n_cell)
            return widths

        pad = [cell_h * float(2 ** (k + 1)) for k in range(n_pad)]

        x_station = _station_cells(xs)
        y_station = _station_cells(ys)

        x_w = np.array(list(reversed(pad)) + x_station + pad, dtype=float)
        y_w = np.array(list(reversed(pad)) + y_station + pad, dtype=float)

        # -- vertical grid ---------------------------------------------
        n_air    = cfg.n_airlayers
        n_active = cfg.nz
        air_h    = cfg.cell_size_v_top
        air_z    = [air_h] * n_air
        earth_z: list[float] = []
        thick = float(cfg.cell_size_v_top)
        for _ in range(n_active):
            earth_z.append(thick)
            thick *= cfg.depth_scale
        z_w = np.array(air_z + earth_z, dtype=float)

        nz_tot = len(z_w)
        ny_tot = len(y_w)
        nx_tot = len(x_w)

        rho_val  = float(np.log(cfg.initial_rho))
        rho_grid = np.full((nz_tot, ny_tot, nx_tot), rho_val, dtype=float)
        if n_air > 0:
            rho_grid[:n_air, :, :] = np.log(1e12)

        obj.x_widths = x_w
        obj.y_widths = y_w
        obj.z_widths = z_w
        obj.n_air    = n_air
        obj.rho_loge = rho_grid

        if obj.verbose:
            obj.logger.info(
                "ModEmModel3D.halfspace: %d×%d×%d grid, rho=%.1f Ω·m",
                nx_tot, ny_tot, nz_tot, cfg.initial_rho,
            )
        return obj

    # ------------------------------------------------------------------
    # I/O
    # ------------------------------------------------------------------
    @classmethod
    def read(cls, path: PathLike, **kwargs) -> "ModEmModel3D":
        """Parse an existing ModEM 3D model file (.ws)."""
        p = Path(path)
        if not p.exists():
            raise FileNotFoundError(f"ModEM 3D model file not found: {p}")
        d   = _parse_model3d(p)
        obj = cls(**kwargs)
        obj.x_widths = d["x_widths"]
        obj.y_widths = d["y_widths"]
        obj.z_widths = d["z_widths"]
        obj.rho_loge = d["rho_loge"]
        obj.n_air    = d["n_air"]
        obj.log_type = d["log_type"]
        if obj.verbose:
            obj.logger.info(
                "ModEmModel3D.read: %d×%d×%d from %s",
                obj.nx, obj.ny, obj.nz, p,
            )
        return obj

    def write(self, path: PathLike) -> Path:
        """Write model to *path* in ModEM WS (WinGLink) 3D format.

        Returns the resolved path.
        """
        p = Path(path)
        p.parent.mkdir(parents=True, exist_ok=True)

        def _floats(arr: np.ndarray, per_row: int = 12) -> str:
            rows = []
            for i in range(0, len(arr), per_row):
                rows.append("  " + "  ".join(f"{v:>12.4f}" for v in arr[i:i+per_row]))
            return "\n".join(rows) + "\n"

        lines: list[str] = [
            f"  {self.nx}  {self.ny}  {self.nz}  {self.n_air}  {self.log_type}\n",
        ]
        lines.append(_floats(self.x_widths))
        lines.append(_floats(self.y_widths))
        lines.append(_floats(self.z_widths))
        lines.append("\n")

        for iz in range(self.nz):
            for iy in range(self.ny):
                lines.append("  " + "  ".join(
                    f"{v:>12.5E}" for v in self.rho_loge[iz, iy, :]
                ) + "\n")

        with p.open("w") as fh:
            fh.writelines(lines)
        return p
