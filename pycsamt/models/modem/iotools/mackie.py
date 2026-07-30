# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Mackie-format model I/O — Python translation of the MATLAB ioAscii/ scripts.

Translates:
  read_mackie2d_model.m + readCond_2D.m   → :func:`read_mackie2d`
  write_mackie2d_model.m + writeCond_2D.m → :func:`write_mackie2d`
  read_mackie3d_model.m + readCond_3D.m   → :func:`read_mackie3d`
  write_mackie3d_model.m + writeCond_3D.m → :func:`write_mackie3d`

2D Mackie format
-----------------
::

    ny  nz  [LOGE]
    <ny y-widths (horizontal)>   -- 10 per line, %11.4G
    <nz z-widths (earth only)>   -- 10 per line, %11.4G
    1                            -- sentinel
    <nz rows of ny rho values>   -- 19 per line, %13.5E

3D Mackie format
-----------------
::

    nx  ny  nz  nzAir  VALUES  [LOGE]
    <nx x-widths>
    <ny y-widths>
    <nz z-widths>
    k [k2]               -- depth layer index, 1-based inclusive range
    <ny rows of nx rho>
    ...
    WINGLINK
    site
    1  1
    cx  cy  cz           -- origin in km
    rotation_deg

Notes
-----
* Widths are in **metres** (consistent with ``ModEmModel2D`` / ``ModEmModel3D``).
* Internally both model classes store ``rho_loge = ln(ρ)``; conversion from
  LINEAR and LOG10 file formats is performed automatically.
* ``readCond_2D`` adds 10 synthetic air layers above the earth grid; this
  behaviour is replicated in :func:`read_mackie2d`.
"""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path
from typing import Union

import numpy as np

PathLike = Union[str, Path]

__all__ = [
    "read_mackie2d",
    "write_mackie2d",
    "read_mackie3d",
    "write_mackie3d",
]

# -------------------------------------------------------------------------
# Constants
# -------------------------------------------------------------------------

_N_AIR_2D = 10  # air layers added by readCond_2D
_AIR_RHO = 1e12  # Ω·m for air cells
_AIR_LOGE = np.log(_AIR_RHO)  # ≈ 27.63
_AIR_THRESH = 20.0  # min loge value considered "air" for auto-detection
_VALUES_WRAP_WIDTH = 10  # widths per line (2D / 3D header widths)
_VALUES_WRAP_RHO_2D = 19  # rho values per line (2D body)


# -------------------------------------------------------------------------
# Internal helpers
# -------------------------------------------------------------------------


def _to_loge(arr: np.ndarray, log_type: str) -> np.ndarray:
    """Convert a raw rho array from *log_type* encoding to ln(ρ)."""
    if log_type == "LOG10":
        return arr * np.log(10.0)
    elif log_type == "LINEAR":
        safe = np.where(arr > 0, arr, np.finfo(float).tiny)
        return np.log(safe)
    else:  # LOGE
        return arr.copy()


def _from_loge(rho_loge: np.ndarray, log_type: str) -> np.ndarray:
    """Convert ln(ρ) to *log_type* encoding for writing."""
    if log_type == "LOG10":
        return rho_loge / np.log(10.0)
    elif log_type == "LINEAR":
        return np.exp(rho_loge)
    else:  # LOGE
        return rho_loge.copy()


def _detect_air_layers(rho_loge: np.ndarray) -> int:
    """Count consecutive top z-layers where every cell is clearly air."""
    n = 0
    for row in rho_loge:
        if float(np.min(row)) > _AIR_THRESH:
            n += 1
        else:
            break
    return n


def _write_widths(
    fh, widths: np.ndarray, per_line: int = _VALUES_WRAP_WIDTH
) -> None:
    """Write a sequence of widths, *per_line* per row, ``%11.4G`` format."""
    for j, v in enumerate(widths, 1):
        fh.write(f"{v:11.4G}")
        if j % per_line == 0 or j == len(widths):
            fh.write("\n")


# -------------------------------------------------------------------------
# 2D Mackie
# -------------------------------------------------------------------------


def _parse_mackie2d(path: Path) -> dict:
    """Low-level parser for the 2D Mackie file format."""
    with path.open("r", errors="replace") as fh:
        lines = [ln.rstrip("\n") for ln in fh]

    i = 0
    N = len(lines)
    while i < N and (not lines[i].strip() or lines[i].strip().startswith("#")):
        i += 1

    ctrl = lines[i].split()
    i += 1
    ny = int(ctrl[0])
    nz = int(ctrl[1])
    log_type = "LINEAR"
    for tok in ctrl:
        if tok.upper() in ("LOGE", "LOG10"):
            log_type = tok.upper()

    # flatten all remaining tokens (skip comment lines)
    all_tokens: list[str] = []
    for ln in lines[i:]:
        stripped = ln.strip()
        if stripped.startswith("#"):
            continue
        all_tokens.extend(stripped.split())

    t = 0
    y_widths = np.array(
        [float(all_tokens[t + k]) for k in range(ny)], dtype=float
    )
    t += ny
    z_earth = np.array(
        [float(all_tokens[t + k]) for k in range(nz)], dtype=float
    )
    t += nz
    t += 1  # skip sentinel integer

    rho_flat = np.array(
        [float(all_tokens[t + k]) for k in range(ny * nz)], dtype=float
    )
    # file layout: z outer, y inner → shape (nz, ny)
    rho_arr = rho_flat.reshape(nz, ny)
    rho_loge_earth = _to_loge(rho_arr, log_type)

    # add air layers — mirrors readCond_2D.m (nzAir = 10)
    n_air = _N_AIR_2D
    h0 = max(float(z_earth[0]), 10.0)
    dz_air_up = [h0 * (3.0**k) for k in range(n_air)]  # surface → up
    z_air = list(reversed(dz_air_up))  # top → surface
    z_full = np.array(z_air + list(z_earth), dtype=float)

    rho_air = np.full((n_air, ny), _AIR_LOGE, dtype=float)
    rho_loge = np.vstack([rho_air, rho_loge_earth])

    return {
        "x_widths": y_widths,  # horizontal dimension
        "z_widths": z_full,
        "rho_loge": rho_loge,
        "n_air": n_air,
        "log_type": log_type,
    }


def read_mackie2d(path: PathLike) -> ModEmModel2D:  # noqa: F821
    """Read a 2D Mackie resistivity model file.

    Equivalent to ``readCond_2D(fname)`` in MATLAB.  Ten synthetic air
    layers are prepended to the vertical grid (matching the MATLAB wrapper).

    Parameters
    ----------
    path : str or Path
        Path to the Mackie 2D model file (typically ``*.rho``).

    Returns
    -------
    ModEmModel2D
        Populated model object with ``n_air`` set to 10.
    """
    from ..model2d import ModEmModel2D

    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"Mackie 2D model file not found: {p}")

    d = _parse_mackie2d(p)
    obj = ModEmModel2D()
    obj.x_widths = d["x_widths"]
    obj.z_widths = d["z_widths"]
    obj.rho_loge = d["rho_loge"]
    obj.n_air = d["n_air"]  # dynamic attribute, not in base class
    obj.log_type = d["log_type"]  # type: ignore[attr-defined]
    return obj


def write_mackie2d(
    model: ModEmModel2D,  # noqa: F821
    path: PathLike,
    n_air: int | None = None,
    log_type: str = "LOGE",
) -> Path:
    """Write a ``ModEmModel2D`` in Mackie 2D format.

    Equivalent to ``writeCond_2D(fname, cond)`` in MATLAB.  Only the earth
    layers are written; air layers are stripped first.

    Parameters
    ----------
    model : ModEmModel2D
    path : str or Path
        Output file path.
    n_air : int, optional
        Number of air layers at the top of the model grid.  If *None*, uses
        ``model.n_air`` if present, otherwise auto-detects from rho values.
    log_type : {'LOGE', 'LINEAR', 'LOG10'}
        Encoding to use in the file.  Defaults to ``'LOGE'``.

    Returns
    -------
    Path
        Resolved output path.
    """
    if n_air is None:
        n_air = getattr(model, "n_air", None)
    if n_air is None:
        n_air = _detect_air_layers(model.rho_loge)

    rho_earth = model.rho_loge[n_air:]  # (nz_earth, nx)
    nz_earth, nx = rho_earth.shape
    z_earth = model.z_widths[n_air:]

    rho_out = _from_loge(rho_earth, log_type)

    p = Path(path)
    with p.open("w") as fh:
        # header: ny nz [LOGE]
        type_str = log_type if log_type != "LINEAR" else ""
        fh.write(f"{nx:3d} {nz_earth:3d} {type_str}\n")

        # y-widths (horizontal), 10 per line
        _write_widths(fh, model.x_widths)

        # z-widths (earth only), 10 per line
        _write_widths(fh, z_earth)

        # sentinel
        fh.write("1\n")

        # rho values: z outer, y inner; 19 per line with leading "  " per z-layer
        for iz in range(nz_earth):
            line = "  "
            for j, v in enumerate(rho_out[iz], 1):
                line += f"{v:13.5E}"
                if j % _VALUES_WRAP_RHO_2D == 0 or j == nx:
                    fh.write(line + "\n")
                    line = ""

    return p


# -------------------------------------------------------------------------
# 3D Mackie
# -------------------------------------------------------------------------


def _parse_mackie3d(path: Path) -> dict:
    """Low-level parser for the 3D Mackie file format."""
    with path.open("r", errors="replace") as fh:
        lines = [ln.rstrip("\n") for ln in fh]

    N = len(lines)
    i = 0
    while i < N and (not lines[i].strip() or lines[i].strip().startswith("#")):
        i += 1

    ctrl = lines[i].split()
    i += 1
    nx = int(ctrl[0])
    ny = int(ctrl[1])
    nz = int(ctrl[2])
    nz_air = int(ctrl[3]) if len(ctrl) > 3 else 0
    log_type = "LINEAR"
    for tok in ctrl:
        if tok.upper() in ("LOGE", "LOG10"):
            log_type = tok.upper()

    # collect width tokens: nx + ny + nz floats
    target = nx + ny + nz
    widths: list[float] = []
    while i < N and len(widths) < target:
        stripped = lines[i].strip()
        i += 1
        if not stripped or stripped.startswith("#"):
            continue
        for tok in stripped.split():
            try:
                widths.append(float(tok))
            except ValueError:
                pass
            if len(widths) == target:
                break

    x_widths = np.array(widths[:nx], dtype=float)
    y_widths = np.array(widths[nx : nx + ny], dtype=float)
    z_widths = np.array(widths[nx + ny :], dtype=float)

    # read blocks: each starts with "k1 [k2]" line then ny rows of nx values
    rho_arr = np.zeros((nz, ny, nx), dtype=float)
    assigned = [False] * nz

    while not all(assigned) and i < N:
        # skip blank lines
        while i < N and not lines[i].strip():
            i += 1
        if i >= N:
            break

        ln = lines[i].strip()
        # parse block index line: expect 1–2 integers
        parts = ln.split()
        ints: list[int] = []
        for p in parts:
            try:
                ints.append(int(p))
            except ValueError:
                break  # hit a non-integer → footer
        if not ints:
            break
        i += 1

        k1 = ints[0] - 1  # 0-indexed
        k2 = ints[1] - 1 if len(ints) >= 2 else k1

        # read ny rows of nx values; each row is on one line
        rows: list[list[float]] = []
        while len(rows) < ny and i < N:
            row_ln = lines[i].strip()
            i += 1
            if not row_ln:
                continue
            rows.append([float(v) for v in row_ln.split()])

        slice_arr = np.array(rows, dtype=float)  # (ny, nx)
        loge_slice = _to_loge(slice_arr, log_type)

        for k in range(k1, k2 + 1):
            if k < nz:
                rho_arr[k] = loge_slice
                assigned[k] = True

    # optional footer: WINGLINK … origin (km) … rotation
    origin = np.zeros(3, dtype=float)
    rotation = 0.0
    while i < N:
        ln = lines[i].strip()
        i += 1
        if not ln:
            continue
        parts = ln.split()
        if len(parts) == 3:
            try:
                vals = [float(p) for p in parts]
                origin = np.array(vals) * 1000.0  # km → m
                break
            except ValueError:
                pass

    while i < N:
        ln = lines[i].strip()
        i += 1
        if not ln:
            continue
        try:
            rotation = float(ln.split()[0])
            break
        except (ValueError, IndexError):
            pass

    return {
        "x_widths": x_widths,
        "y_widths": y_widths,
        "z_widths": z_widths,
        "rho_loge": rho_arr,
        "n_air": nz_air,
        "log_type": log_type,
        "origin": origin,
        "rotation": rotation,
    }


def read_mackie3d(path: PathLike) -> ModEmModel3D:  # noqa: F821
    """Read a 3D Mackie resistivity model file.

    Equivalent to ``readCond_3D(fname, format=1)`` in MATLAB.

    Parameters
    ----------
    path : str or Path
        Path to the Mackie 3D model file.

    Returns
    -------
    ModEmModel3D
        Populated model object.
    """
    from ..model3d import ModEmModel3D

    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"Mackie 3D model file not found: {p}")

    d = _parse_mackie3d(p)
    obj = ModEmModel3D()
    obj.x_widths = d["x_widths"]
    obj.y_widths = d["y_widths"]
    obj.z_widths = d["z_widths"]
    obj.rho_loge = d["rho_loge"]
    obj.n_air = d["n_air"]
    obj.log_type = d["log_type"]
    obj.origin = d["origin"]  # type: ignore[attr-defined]
    obj.rotation = d["rotation"]  # type: ignore[attr-defined]
    return obj


def write_mackie3d(
    model: ModEmModel3D,  # noqa: F821
    path: PathLike,
    log_type: str = "LOGE",
    origin: Sequence[float] = (0.0, 0.0, 0.0),
    rotation: float = 0.0,
) -> Path:
    """Write a ``ModEmModel3D`` in Mackie 3D format.

    Equivalent to ``writeCond_3D(fname, cond, format=1)`` in MATLAB.
    Each depth layer is written as a separate block with a single integer
    index on its own line.

    Parameters
    ----------
    model : ModEmModel3D
    path : str or Path
        Output file path.
    log_type : {'LOGE', 'LINEAR', 'LOG10'}
        Encoding to use.  Defaults to ``'LOGE'``.
    origin : sequence of 3 floats
        Grid centre ``(x, y, z)`` in **metres**.  Written as km.
    rotation : float
        Grid rotation angle in degrees.

    Returns
    -------
    Path
        Resolved output path.
    """
    nx, ny, nz = model.nx, model.ny, model.nz
    n_air = model.n_air

    rho_out = _from_loge(model.rho_loge, log_type)

    p = Path(path)
    with p.open("w") as fh:
        # header
        fh.write(f"{nx} {ny} {nz} {n_air} VALUES {log_type}\n")

        # widths: space-separated on one line each
        fh.write(" ".join(f"{v:G}" for v in model.x_widths) + "\n")
        fh.write(" ".join(f"{v:G}" for v in model.y_widths) + "\n")
        fh.write(" ".join(f"{v:G}" for v in model.z_widths) + "\n")

        # one block per depth layer
        for k in range(nz):
            fh.write(f"{k + 1}\n")
            for j in range(ny):
                row = rho_out[k, j, :]
                fh.write("".join(f"{v:13.5E}" for v in row) + "\n")

        # footer (WINGLINK compatibility)
        ox, oy, oz = [v / 1000.0 for v in origin]  # m → km
        fh.write("WINGLINK \n")
        fh.write("site \n")
        fh.write("1 1 \n")
        fh.write(f"{int(ox)} {int(oy)} {int(oz)}\n")
        fh.write(f"{int(rotation)}\n")

    return p
