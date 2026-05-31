# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

"""
TEM data readers.

Currently implemented
---------------------
:func:`read_xyz`   Generic white-space / CSV file with columns
                   ``time  data  [error]``.

Planned (stubs raise ``NotImplementedError`` until implemented)
---------------------------------------------------------------
:func:`read_geosoft_dat`   Geosoft DAT format (Oasis Montaj)
:func:`read_amira`         AMIRA / EMIT .tem format
:func:`read_zonge`         Zonge GDP .avg / .tem format
:func:`read_walkttem`      WalkTEM / Aarhus Workbench format
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional, Sequence, Union

import numpy as np

from ._base import TEMSounding

__all__ = [
    "read_xyz",
    "read_geosoft_dat",
    "read_amira",
    "read_zonge",
    "read_walkttem",
]

Pathish = Union[str, Path]

# ---------------------------------------------------------------------------
# Generic XYZ / CSV reader  (fully implemented)
# ---------------------------------------------------------------------------

def read_xyz(
    path: Pathish,
    *,
    current: float,
    tx_area: Optional[float] = None,
    loop_side: Optional[float] = None,
    loop_radius: Optional[float] = None,
    data_type: str = "dBdt",
    rx_area: float = 1.0,
    rx_turns: int = 1,
    tx_turns: int = 1,
    offset: float = 0.0,
    station_name: str = "",
    x: float = 0.0,
    y: float = 0.0,
    elevation: float = 0.0,
    time_col: int = 0,
    data_col: int = 1,
    error_col: Optional[int] = None,
    skip_header: int = 0,
    delimiter: Optional[str] = None,
    time_unit: str = "s",
    data_unit: str = "SI",
) -> TEMSounding:
    r"""
    Read a TEM sounding from a plain-text XYZ / CSV file.

    The file may contain any number of columns; only ``time_col``
    and ``data_col`` (and optionally ``error_col``) are read.
    Comment lines starting with ``#`` or ``!`` are ignored.

    Parameters
    ----------
    path : str or Path
        Path to the text file.
    current : float
        Transmitter current in Amperes.
    tx_area : float, optional
        Transmitter loop area in m².  Supply this *or* ``loop_side``
        *or* ``loop_radius``.
    loop_side : float, optional
        Side of a square transmitter loop in metres.
    loop_radius : float, optional
        Radius of a circular transmitter loop in metres.
    data_type : str, optional
        Units / type of the ``data_col`` column.  One of
        ``"dBdt"``, ``"dHdt"``, ``"voltage"``,
        ``"normalized_voltage"``.  Default ``"dBdt"``.
    rx_area : float, optional
        Receiver coil area in m².  Default 1.0.
    rx_turns : int, optional
        Receiver coil turns.  Default 1.
    tx_turns : int, optional
        Transmitter turns.  Default 1.
    offset : float, optional
        Tx–Rx horizontal separation in m.  Default 0 (coincident).
    station_name : str, optional
        Name to store in the output :class:`TEMSounding`.
    x, y, elevation : float, optional
        Station coordinates in metres.
    time_col : int, optional
        Zero-based column index for time values.  Default 0.
    data_col : int, optional
        Zero-based column index for the decay data.  Default 1.
    error_col : int, optional
        Zero-based column index for uncertainty.  Default ``None``
        (no error column).
    skip_header : int, optional
        Number of header lines to skip before data rows.  Default 0.
    delimiter : str, optional
        Column separator.  ``None`` means any whitespace (default).
    time_unit : str, optional
        Unit of the time column: ``"s"`` (default), ``"ms"``, ``"us"``.
    data_unit : str, optional
        Unit of the data column: ``"SI"`` (no conversion), ``"nT/s"``
        (convert nT/s → T/s), ``"pT/s"``, ``"nV/Am2"``.

    Returns
    -------
    TEMSounding
        Populated :class:`~pycsamt.tdem.TEMSounding` instance.

    Examples
    --------
    A minimal two-column file (time in ms, dB/dt in nT/s):

    .. code-block:: text

        # TEM sounding, site A01
        0.010   1.24e+04
        0.018   8.91e+03
        0.031   5.23e+03

    >>> snd = read_xyz(
    ...     "site_A01.xyz",
    ...     current=8.0,
    ...     loop_side=100.0,
    ...     time_unit="ms",
    ...     data_unit="nT/s",
    ... )
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"TEM data file not found: {path}")

    rows: List[List[str]] = []
    with open(path, "r") as fh:
        for ln in fh:
            stripped = ln.strip()
            if not stripped or stripped[0] in ("#", "!"):
                continue
            rows.append(stripped.split(delimiter))

    if skip_header:
        rows = rows[skip_header:]

    if not rows:
        raise ValueError(f"No data rows found in {path}")

    def _col(rows, idx):
        try:
            return np.array([float(r[idx]) for r in rows])
        except (IndexError, ValueError) as exc:
            raise ValueError(
                f"Cannot read column {idx} from {path}: {exc}"
            ) from exc

    t = _col(rows, time_col)
    d = _col(rows, data_col)
    err = _col(rows, error_col) if error_col is not None else None

    # unit conversions
    _time_scale = {"s": 1.0, "ms": 1e-3, "us": 1e-6, "µs": 1e-6}
    if time_unit not in _time_scale:
        raise ValueError(f"time_unit must be one of {list(_time_scale)}, got '{time_unit}'")
    t = t * _time_scale[time_unit]

    _data_scale = {
        "SI": 1.0,
        "nT/s": 1e-9,
        "pT/s": 1e-12,
        "nV/Am2": 1e-9,
        "uV/Am2": 1e-6,
    }
    if data_unit not in _data_scale:
        raise ValueError(f"data_unit must be one of {list(_data_scale)}, got '{data_unit}'")
    scale = _data_scale[data_unit]
    d = d * scale
    if err is not None:
        err = err * scale

    return TEMSounding.from_arrays(
        t, d,
        current=current,
        tx_area=tx_area,
        loop_side=loop_side,
        loop_radius=loop_radius,
        data_type=data_type,
        tx_turns=tx_turns,
        rx_area=rx_area,
        rx_turns=rx_turns,
        offset=offset,
        station_name=station_name or path.stem,
        x=x,
        y=y,
        elevation=elevation,
        error=err,
    )


def read_xyz_multisite(
    path: Pathish,
    *,
    site_col: int = 0,
    time_col: int = 1,
    data_col: int = 2,
    error_col: Optional[int] = None,
    **kwargs,
) -> Dict[str, TEMSounding]:
    r"""
    Read multiple TEM soundings from a single XYZ file where the
    first column is a site identifier.

    Returns
    -------
    dict[str, TEMSounding]
        Mapping from site name to :class:`TEMSounding`.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(path)

    records: Dict[str, List] = {}
    with open(path, "r") as fh:
        for ln in fh:
            stripped = ln.strip()
            if not stripped or stripped[0] in ("#", "!"):
                continue
            parts = stripped.split()
            site = parts[site_col]
            if site not in records:
                records[site] = []
            records[site].append(parts)

    out = {}
    for site, rows in records.items():
        t = np.array([float(r[time_col]) for r in rows])
        d = np.array([float(r[data_col]) for r in rows])
        err = None
        if error_col is not None:
            err = np.array([float(r[error_col]) for r in rows])
        snd = TEMSounding.from_arrays(
            t, d,
            station_name=site,
            error=err,
            **kwargs,
        )
        out[site] = snd

    return out


# ---------------------------------------------------------------------------
# Geosoft DAT reader  (stub)
# ---------------------------------------------------------------------------

def read_geosoft_dat(
    path: Pathish,
    *,
    current: float,
    tx_area: Optional[float] = None,
    loop_side: Optional[float] = None,
    **kwargs,
) -> List[TEMSounding]:
    r"""
    Read TEM data from a Geosoft DAT file (Oasis Montaj format).

    .. note::
        **Not yet implemented.**  The Geosoft DAT format stores
        one line per station with columns for each time gate.
        Support will be added in a future release.

    Parameters
    ----------
    path : str or Path
        Path to the ``.dat`` file.
    current : float
        Transmitter current in Amperes.
    tx_area : float, optional
    loop_side : float, optional

    Raises
    ------
    NotImplementedError
        Always, until the reader is implemented.
    """
    raise NotImplementedError(
        "read_geosoft_dat is not yet implemented. "
        "Convert your Geosoft DAT to XYZ and use read_xyz() instead."
    )


# ---------------------------------------------------------------------------
# AMIRA / EMIT .tem reader  (stub)
# ---------------------------------------------------------------------------

def read_amira(
    path: Pathish,
    **kwargs,
) -> List[TEMSounding]:
    r"""
    Read TEM data from an AMIRA / EMIT ``.tem`` file.

    .. note::
        **Not yet implemented.**

    Raises
    ------
    NotImplementedError
        Always, until the reader is implemented.
    """
    raise NotImplementedError(
        "read_amira is not yet implemented. "
        "Export from EMIT to XYZ and use read_xyz() instead."
    )


# ---------------------------------------------------------------------------
# Zonge GDP reader  (stub)
# ---------------------------------------------------------------------------

def read_zonge(
    path: Pathish,
    **kwargs,
) -> List[TEMSounding]:
    r"""
    Read TEM data from a Zonge GDP ``.avg`` / ``.tem`` file.

    .. note::
        **Not yet implemented.**

    Raises
    ------
    NotImplementedError
        Always, until the reader is implemented.
    """
    raise NotImplementedError(
        "read_zonge is not yet implemented. "
        "Export from GDP to XYZ and use read_xyz() instead."
    )


# ---------------------------------------------------------------------------
# WalkTEM / Aarhus Workbench reader  (stub)
# ---------------------------------------------------------------------------

def read_walkttem(
    path: Pathish,
    **kwargs,
) -> List[TEMSounding]:
    r"""
    Read TEM data from a WalkTEM / Aarhus Workbench file.

    .. note::
        **Not yet implemented.**

    Raises
    ------
    NotImplementedError
        Always, until the reader is implemented.
    """
    raise NotImplementedError(
        "read_walkttem is not yet implemented."
    )
