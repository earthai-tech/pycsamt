# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Old-format impedance Z file I/O — Python translation of the MATLAB ioAscii/ scripts.

Translates (pre-2011 format, (c) Anna Kelbert 2008):
  readZ_3D_old.m  / writeZ_3D_old.m  → :func:`read_z3d_old`  / :func:`write_z3d_old`
  readZ_2D_old.m  / writeZ_2D_old.m  → :func:`read_z2d_old`  / :func:`write_z2d_old`
  writeZ_old2list_3D.m               → :func:`convert_z3d`
  writeZ_old2list_2D.m               → :func:`convert_z2d`

Also provides writers for the current ModEM list format (2011+):
  writeZ_3D.m (2011)                 → :func:`write_z3d_list`
  writeZ_2D.m (2011)                 → :func:`write_z2d_list`

Old 3D format
--------------
::

    Description: <info>
    Units: <units>
    Sign convention: <±1>

    <nTx>
    <T>  <nComp>  <nSites>
    <X1> ... <XN>
    <Y1> ... <YN>
    <Z1> ... <ZN>
    <comp1 comp1 comp2 comp2 ...>     -- nComp tokens (Re/Im pairs)
    <code>  <Re1> <Im1> <Re2> ...     -- data
              <err1> <err1> <err2> ... -- error
    ...

Old 2D format
--------------
Same structure but ``<T> <MODE> <nSites>`` (MODE = ``TE`` or ``TM``) and
2-column site locations (Y, Z only).

Current list format (3D)
-------------------------
::

    # <description>
    # Period(s) Code GG_Lat GG_Lon X(m) Y(m) Z(m) Component Real Imag Error
    > Full_Impedance (or Off_Diagonal_Impedance)
    > exp(-i\\omega t)
    > <units>
    > <orientation>
    > <origin_x> <origin_y>
    > <nTx> <nSites>
    <period> <code> <lat> <lon> <X> <Y> <Z> <comp> <Re> <Im> <err>
    ...
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Union

import numpy as np

PathLike = Union[str, Path]

__all__ = [
    "ZBlock",
    "ImpedanceFile",
    "read_z3d_old",
    "write_z3d_old",
    "read_z2d_old",
    "write_z2d_old",
    "write_z3d_list",
    "write_z2d_list",
    "convert_z3d",
    "convert_z2d",
]

# -------------------------------------------------------------------------
# Data structures
# -------------------------------------------------------------------------


@dataclass
class ZBlock:
    """One transmitter (period) block of impedance data.

    Attributes
    ----------
    period : float
        Period in seconds.
    site_names : list of str
        Site codes, length ``n_sites``.
    site_loc : ndarray, shape (n_sites, 3)
        Site locations in metres: ``[X, Y, Z]``.
    comp_names : list of str
        Complex component names, e.g. ``['ZXX', 'ZXY', 'ZYX', 'ZYY']``.
    Z : ndarray, shape (n_sites, n_comp), complex
        Impedance values.
    Zerr : ndarray, shape (n_sites, n_comp), float
        Impedance error estimates (standard deviations).
    mode : str
        ``'TE'``, ``'TM'``, or ``''`` (3-D).
    lat, lon : ndarray, shape (n_sites,) or None
        Geographic coordinates (optional).
    """

    period: float
    site_names: list[str]
    site_loc: np.ndarray  # (n_sites, 3) in m
    comp_names: list[str]
    Z: np.ndarray  # (n_sites, n_comp) complex
    Zerr: np.ndarray  # (n_sites, n_comp) float
    mode: str = ""
    lat: np.ndarray | None = None
    lon: np.ndarray | None = None

    @property
    def n_sites(self) -> int:
        return len(self.site_names)

    @property
    def n_comp(self) -> int:
        return len(self.comp_names)


@dataclass
class ImpedanceFile:
    """Container for an impedance data file (old or new format).

    Attributes
    ----------
    description : str
        Free-text description from the file header.
    units : str
        Impedance units, e.g. ``'[V/m]/[T]'`` or ``'Ohm'``.
    sign : int
        Time-variation sign convention: ``-1`` for ``exp(-iωt)``, ``+1`` otherwise.
    blocks : list of ZBlock
        One entry per period.
    origin : tuple of 3 floats
        Grid origin ``(x, y, z)`` in metres.
    orientation : float
        Grid orientation angle in degrees.
    """

    description: str = "ModEM impedance data"
    units: str = "[V/m]/[T]"
    sign: int = -1
    blocks: list[ZBlock] = field(default_factory=list)
    origin: tuple[float, float, float] = (0.0, 0.0, 0.0)
    orientation: float = 0.0

    @property
    def periods(self) -> np.ndarray:
        return np.array([b.period for b in self.blocks])

    @property
    def n_periods(self) -> int:
        return len(self.blocks)


# -------------------------------------------------------------------------
# Internal helpers
# -------------------------------------------------------------------------


def _read_lines(path: Path) -> list[str]:
    with path.open("r", errors="replace") as fh:
        return [ln.rstrip("\n") for ln in fh]


def _parse_header_3line(
    lines: list[str], i: int
) -> tuple[str, str, int, int]:
    """Parse the 3-line old-format header.  Returns (info, units, isign, next_i)."""
    info = (
        lines[i].split(":", 1)[1].strip()
        if ":" in lines[i]
        else lines[i].strip()
    )
    units = (
        lines[i + 1].split(":", 1)[1].strip()
        if ":" in lines[i + 1]
        else lines[i + 1].strip()
    )
    sign_str = (
        lines[i + 2].split(":", 1)[1].strip()
        if ":" in lines[i + 2]
        else lines[i + 2].strip()
    )
    try:
        isign = int(float(sign_str))
    except (ValueError, IndexError):
        isign = -1
    return info, units, isign, i + 3


def _skip_blank(lines: list[str], i: int) -> int:
    while i < len(lines) and not lines[i].strip():
        i += 1
    return i


def _read_n_floats(
    lines: list[str], i: int, n: int
) -> tuple[list[float], int]:
    """Collect exactly *n* float tokens from consecutive lines starting at *i*."""
    vals: list[float] = []
    while len(vals) < n and i < len(lines):
        for tok in lines[i].split():
            try:
                vals.append(float(tok))
            except ValueError:
                pass
        i += 1
    return vals, i


def _sign_str(sign: int) -> str:
    return r"exp(-i\omega t)" if sign == -1 else r"exp(+i\omega t)"


def _write_header(fh, imp: ImpedanceFile) -> None:
    fh.write(f"Description: {imp.description}\n")
    fh.write(f"Units: {imp.units}\n")
    fh.write(f"Sign convention: {imp.sign}\n\n")


# -------------------------------------------------------------------------
# Old 3D format
# -------------------------------------------------------------------------


def _parse_z3d_old_block(lines: list[str], i: int) -> tuple[ZBlock, int]:
    """Parse one transmitter block from the old 3D format.  Returns (block, next_i)."""
    i = _skip_blank(lines, i)
    parts = lines[i].split()
    i += 1
    period = float(parts[0])
    n_comp = int(parts[1])  # total real components = 2 × n_complex
    n_sites = int(parts[2])
    n_cplx = n_comp // 2

    # site locations: 3 rows of n_sites values (X, Y, Z each on its own line)
    xs, i = _read_n_floats(lines, i, n_sites)
    ys, i = _read_n_floats(lines, i, n_sites)
    zs, i = _read_n_floats(lines, i, n_sites)
    site_loc = np.column_stack([xs, ys, zs])  # (n_sites, 3)

    # component header: n_comp tokens; even-indexed give the component name
    # (old format repeats: ZXX ZXX ZXY ZXY ZYX ZYX ZYY ZYY)
    i = _skip_blank(lines, i)
    comp_tokens = lines[i].split()
    i += 1
    comp_names = [
        comp_tokens[k] for k in range(0, min(len(comp_tokens), n_comp), 2)
    ]

    # data: n_sites × 2 lines (data + error)
    site_names: list[str] = []
    Z = np.full((n_sites, n_cplx), np.nan, dtype=complex)
    Zerr = np.full((n_sites, n_cplx), np.nan, dtype=float)

    for k in range(n_sites):
        i = _skip_blank(lines, i)
        data_line = lines[i].split()
        i += 1
        site_names.append(data_line[0])
        d_vals = [float(v) for v in data_line[1:]]
        Z[k, :] = [
            d_vals[2 * c] + 1j * d_vals[2 * c + 1] for c in range(n_cplx)
        ]
        # error line (may start with whitespace-only padding)
        i = _skip_blank(lines, i)
        err_vals = [float(v) for v in lines[i].split()]
        i += 1
        Zerr[k, :] = [err_vals[2 * c] for c in range(n_cplx)]

    return ZBlock(
        period=period,
        site_names=site_names,
        site_loc=site_loc,
        comp_names=comp_names,
        Z=Z,
        Zerr=Zerr,
    ), i


def read_z3d_old(path: PathLike) -> ImpedanceFile:
    """Read an old-format (pre-2011) 3D impedance file.

    Equivalent to ``readZ_3D(cfile)`` / ``readZ_3D_old(cfile)`` in MATLAB
    (the 2008 version).

    Parameters
    ----------
    path : str or Path

    Returns
    -------
    ImpedanceFile
    """
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"Impedance file not found: {p}")

    lines = _read_lines(p)
    info, units, isign, i = _parse_header_3line(lines, 0)
    i = _skip_blank(lines, i)
    n_tx = int(lines[i].strip())
    i += 1

    blocks: list[ZBlock] = []
    for _ in range(n_tx):
        block, i = _parse_z3d_old_block(lines, i)
        blocks.append(block)

    # optional trailing origin
    origin = (0.0, 0.0, 0.0)
    while i < len(lines):
        ln = lines[i].strip()
        i += 1
        if ln.lower().startswith("origin"):
            try:
                vals = [float(v) for v in ln.split(":", 1)[1].split()]
                if len(vals) >= 3:
                    origin = (vals[0], vals[1], vals[2])
                elif len(vals) == 2:
                    origin = (vals[0], vals[1], 0.0)
            except (ValueError, IndexError):
                pass
            break

    return ImpedanceFile(
        description=info,
        units=units,
        sign=isign,
        blocks=blocks,
        origin=origin,
    )


def write_z3d_old(imp: ImpedanceFile, path: PathLike) -> Path:
    """Write an old-format (pre-2011) 3D impedance file.

    Equivalent to ``writeZ_3D(cfile, allData, info, units, isign)`` in MATLAB
    (the 2008 version).

    Parameters
    ----------
    imp : ImpedanceFile
    path : str or Path

    Returns
    -------
    Path
    """
    p = Path(path)
    with p.open("w") as fh:
        _write_header(fh, imp)
        fh.write(f"{imp.n_periods}\n")

        for blk in imp.blocks:
            n_comp_total = blk.n_comp * 2  # real components
            fh.write(
                f"{blk.period:12.6E} {n_comp_total:5d} {blk.n_sites:5d}\n"
            )
            # site locations: 3 rows (X, Y, Z)
            for col in range(3):
                row = " ".join(f"{v:12.3f}" for v in blk.site_loc[:, col])
                fh.write(row + "\n")
            # component header: each name repeated twice (Re Im pair)
            comp_hdr = "".join(
                f"{'':>10s} {c}{'':>10s} {c}" for c in blk.comp_names
            )
            fh.write(comp_hdr.rstrip() + "\n")
            # data and errors
            for k in range(blk.n_sites):
                # data line
                data_str = f"{blk.site_names[k]:>10s}"
                for c in range(blk.n_comp):
                    data_str += (
                        f"{blk.Z[k, c].real:15.6E} {blk.Z[k, c].imag:15.6E}"
                    )
                fh.write(data_str + "\n")
                # error line
                err_str = " " * 10
                for c in range(blk.n_comp):
                    err_str += (
                        f"{blk.Zerr[k, c]:15.6E} {blk.Zerr[k, c]:15.6E}"
                    )
                fh.write(err_str + "\n")

        ox, oy, oz = imp.origin
        if ox != 0.0 or oy != 0.0 or oz != 0.0:
            fh.write(f"\nOrigin: {ox:12.6f} {oy:12.6f} {oz:12.6f}\n")

    return p


# -------------------------------------------------------------------------
# Old 2D format
# -------------------------------------------------------------------------


def _parse_z2d_old_block(lines: list[str], i: int) -> tuple[ZBlock, int]:
    """Parse one transmitter block from the old 2D format."""
    i = _skip_blank(lines, i)
    parts = lines[i].split()
    i += 1
    period = float(parts[0])
    mode = parts[1].strip()  # 'TE' or 'TM'
    n_sites = int(parts[2])
    n_cplx = 1  # one complex component per TE or TM block

    # site locations: 2 rows of n_sites values (Y-offset, Z-depth)
    ys, i = _read_n_floats(lines, i, n_sites)
    zs, i = _read_n_floats(lines, i, n_sites)
    xs = [0.0] * n_sites
    site_loc = np.column_stack([xs, ys, zs])

    # component header (skip)
    i = _skip_blank(lines, i)
    i += 1

    site_names: list[str] = []
    Z = np.full((n_sites, n_cplx), np.nan, dtype=complex)
    Zerr = np.full((n_sites, n_cplx), np.nan, dtype=float)

    for k in range(n_sites):
        i = _skip_blank(lines, i)
        data_line = lines[i].split()
        i += 1
        site_names.append(data_line[0])
        Z[k, 0] = float(data_line[1]) + 1j * float(data_line[2])
        # error line
        i = _skip_blank(lines, i)
        err_vals = [float(v) for v in lines[i].split()]
        i += 1
        Zerr[k, 0] = err_vals[0] if err_vals else np.nan

    comp = ["ZTE"] if mode == "TE" else ["ZTM"]
    return ZBlock(
        period=period,
        site_names=site_names,
        site_loc=site_loc,
        comp_names=comp,
        Z=Z,
        Zerr=Zerr,
        mode=mode,
    ), i


def read_z2d_old(path: PathLike) -> ImpedanceFile:
    """Read an old-format (pre-2011) 2D impedance file.

    Equivalent to ``readZ_2D(cfile)`` (the 2008 version) in MATLAB.

    Parameters
    ----------
    path : str or Path

    Returns
    -------
    ImpedanceFile
    """
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"Impedance file not found: {p}")

    lines = _read_lines(p)
    info, units, isign, i = _parse_header_3line(lines, 0)
    i = _skip_blank(lines, i)
    n_tx = int(lines[i].strip())
    i += 1

    blocks: list[ZBlock] = []
    for _ in range(n_tx):
        block, i = _parse_z2d_old_block(lines, i)
        blocks.append(block)

    return ImpedanceFile(
        description=info, units=units, sign=isign, blocks=blocks
    )


def write_z2d_old(imp: ImpedanceFile, path: PathLike) -> Path:
    """Write an old-format (pre-2011) 2D impedance file.

    Equivalent to ``writeZ_2D(cfile, allData, ...)`` (the 2008 version).

    Parameters
    ----------
    imp : ImpedanceFile
    path : str or Path

    Returns
    -------
    Path
    """
    p = Path(path)
    with p.open("w") as fh:
        _write_header(fh, imp)
        fh.write(f"{imp.n_periods}\n")

        for blk in imp.blocks:
            mode = blk.mode or "TE"
            fh.write(f"{blk.period:12.6E} {mode:>5s} {blk.n_sites:5d}\n")
            # 2-column locations (Y, Z)
            fh.write(
                " ".join(
                    f"{blk.site_loc[k, 1]:12.3f}" for k in range(blk.n_sites)
                )
                + "\n"
            )
            fh.write(
                " ".join(
                    f"{blk.site_loc[k, 2]:12.3f}" for k in range(blk.n_sites)
                )
                + "\n"
            )
            # generic component header
            fh.write(f"{'':>12s} Re {'':>12s} Im\n")
            for k in range(blk.n_sites):
                data_str = f"{blk.site_names[k]:>10s}"
                data_str += (
                    f"{blk.Z[k, 0].real:15.6E} {blk.Z[k, 0].imag:15.6E}"
                )
                fh.write(data_str + "\n")
                err_str = " " * 10
                err_str += f"{blk.Zerr[k, 0]:15.6E} {blk.Zerr[k, 0]:15.6E}"
                fh.write(err_str + "\n")

    return p


# -------------------------------------------------------------------------
# Current list format writers
# -------------------------------------------------------------------------

_COMP_4 = ["ZXX", "ZXY", "ZYX", "ZYY"]
_COMP_2 = ["ZXY", "ZYX"]

_DATA_TYPE = {
    4: "Full_Impedance",
    2: "Off_Diagonal_Impedance",
}


def write_z3d_list(imp: ImpedanceFile, path: PathLike) -> Path:
    """Write impedance data in the current ModEM 3D list format.

    Equivalent to ``writeZ_3D(cfile, allData, ...)`` (the 2011 version).
    Each datum is one line:
    ``period code lat lon X Y Z comp Re Im err``

    Parameters
    ----------
    imp : ImpedanceFile
    path : str or Path

    Returns
    -------
    Path
    """
    p = Path(path)

    # Determine data type and component list from first block
    n_cplx = imp.blocks[0].n_comp if imp.blocks else 4
    data_type = _DATA_TYPE.get(n_cplx, "Full_Impedance")
    n_tx = imp.n_periods
    n_sites = max((b.n_sites for b in imp.blocks), default=0)
    ox, oy, _oz = imp.origin
    sign_s = _sign_str(imp.sign)

    with p.open("w") as fh:
        fh.write(f"# {imp.description}\n")
        fh.write(
            "# Period(s) Code GG_Lat GG_Lon X(m) Y(m) Z(m) Component Real Imag Error\n"
        )
        fh.write(f"> {data_type}\n")
        fh.write(f"> {sign_s}\n")
        fh.write(f"> {imp.units}\n")
        fh.write(f"> {imp.orientation:.2f}\n")
        fh.write(f"> {ox:.3f} {oy:.3f}\n")
        fh.write(f"> {n_tx} {n_sites}\n")

        for blk in imp.blocks:
            T = blk.period
            lat_arr = (
                blk.lat if blk.lat is not None else np.zeros(blk.n_sites)
            )
            lon_arr = (
                blk.lon if blk.lon is not None else np.zeros(blk.n_sites)
            )
            for k in range(blk.n_sites):
                code = blk.site_names[k]
                lat = float(lat_arr[k])
                lon = float(lon_arr[k])
                x, y, z = blk.site_loc[k]
                for c, cname in enumerate(blk.comp_names):
                    if np.isnan(blk.Z[k, c]):
                        continue
                    re = blk.Z[k, c].real
                    im = blk.Z[k, c].imag
                    err = float(blk.Zerr[k, c])
                    fh.write(
                        f"{T:12.6E} {code:s} {lat:8.3f} {lon:8.3f} "
                        f"{x:12.3f} {y:12.3f} {z:12.3f} "
                        f"{cname:s} {re:15.6E} {im:15.6E} {err:15.6E}\n"
                    )

    return p


def write_z2d_list(imp: ImpedanceFile, path: PathLike) -> Path:
    """Write impedance data in the current ModEM 2D list format.

    Equivalent to ``writeZ_2D(cfile, allData, ...)`` (the 2011 version).
    TE and TM blocks are written separately.

    Parameters
    ----------
    imp : ImpedanceFile
    path : str or Path

    Returns
    -------
    Path
    """
    p = Path(path)

    te_blocks = [b for b in imp.blocks if b.mode == "TE"]
    tm_blocks = [b for b in imp.blocks if b.mode == "TM"]
    ox, oy, _oz = imp.origin
    sign_s = _sign_str(imp.sign)

    def _write_mode_block(fh, blocks: list, data_type: str) -> None:
        n_tx = len(blocks)
        n_sites = max((b.n_sites for b in blocks), default=0)
        fh.write(f"# {imp.description}\n")
        fh.write(
            "# Period(s) Code GG_Lat GG_Lon X(m) Y(m) Z(m) Component Real Imag Error\n"
        )
        fh.write(f"> {data_type}\n")
        fh.write(f"> {sign_s}\n")
        fh.write(f"> {imp.units}\n")
        fh.write(f"> {imp.orientation:.2f}\n")
        fh.write(f"> {ox:.3f} {oy:.3f}\n")
        fh.write(f"> {n_tx} {n_sites}\n")
        for blk in blocks:
            T = blk.period
            lat_arr = (
                blk.lat if blk.lat is not None else np.zeros(blk.n_sites)
            )
            lon_arr = (
                blk.lon if blk.lon is not None else np.zeros(blk.n_sites)
            )
            for k in range(blk.n_sites):
                code = blk.site_names[k]
                lat = float(lat_arr[k])
                lon = float(lon_arr[k])
                _x, y, z = blk.site_loc[k]
                z_val = blk.Z[k, 0]
                if np.isnan(z_val):
                    continue
                err = float(blk.Zerr[k, 0])
                fh.write(
                    f"{T:12.6E} {code:s} {lat:8.3f} {lon:8.3f} "
                    f"{0.0:12.3f} {y:12.3f} {z:12.3f} "
                    f"ZXY {z_val.real:15.6E} {z_val.imag:15.6E} {err:15.6E}\n"
                )

    with p.open("w") as fh:
        if te_blocks:
            _write_mode_block(fh, te_blocks, "TE_Impedance")
        if tm_blocks:
            _write_mode_block(fh, tm_blocks, "TM_Impedance")

    return p


# -------------------------------------------------------------------------
# Conversion helpers
# -------------------------------------------------------------------------


def convert_z3d(old_path: PathLike, new_path: PathLike) -> Path:
    """Convert old 3D impedance format to current ModEM list format.

    Mirrors ``writeZ_old2list_3D.m``.

    Parameters
    ----------
    old_path : str or Path
        Input file in old (pre-2011) format.
    new_path : str or Path
        Output file in current ModEM list format.

    Returns
    -------
    Path
    """
    imp = read_z3d_old(old_path)
    return write_z3d_list(imp, new_path)


def convert_z2d(old_path: PathLike, new_path: PathLike) -> Path:
    """Convert old 2D impedance format to current ModEM list format.

    Mirrors ``writeZ_old2list_2D.m``.

    Parameters
    ----------
    old_path : str or Path
        Input file in old (pre-2011) format.
    new_path : str or Path
        Output file in current ModEM list format.

    Returns
    -------
    Path
    """
    imp = read_z2d_old(old_path)
    return write_z2d_list(imp, new_path)
