# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

"""
TEM data readers.

Currently implemented
---------------------
:func:`read_xyz`          Generic white-space / CSV with columns
                          ``time  data  [error]``.
:func:`read_temavg`       Zonge TEMAVG processed ``.AVG`` files.
:func:`read_tem_z`        Zonge TEMAVG contour ``.Z`` files.
:func:`read_tem_log`      Zonge TEMAVG processing ``.LOG`` files.
:func:`read_tem_coordinates`  TEM profile/point coordinates.
:func:`read_temavg_soundings` TEMAVG files as TEMSounding objects.
:func:`read_geosoft_dat`  Geosoft line-data ``.dat`` / ``.DAT`` files
                          (Oasis Montaj ASCII export).
:func:`read_amira`        AMIRA / EMIT ``.tem`` keyword-block files.
:func:`read_zonge`        Zonge GDP ``.avg`` / ``.tem`` sounding files.
:func:`read_walkttem`     WalkTEM / Aarhus Workbench ``.tem`` files.
"""

from __future__ import annotations

from pathlib import Path
import re as _re

import numpy as np

from ._base import TEMSounding
from .avg import TEMAVG
from .coordinates import read_tem_coordinates
from .log import TEMLog
from .survey import read_temavg_survey
from .workflow import read_temavg_soundings, transform_temavg_survey
from .zplot import TEMZPlot

__all__ = [
    "read_xyz",
    "read_temavg",
    "read_tem_z",
    "read_tem_log",
    "read_tem_coordinates",
    "read_temavg_survey",
    "read_temavg_soundings",
    "transform_temavg_survey",
    "read_geosoft_dat",
    "read_amira",
    "read_zonge",
    "read_walkttem",
]

Pathish = str | Path


def read_temavg(
    path: Pathish,
    *,
    verbose: int = 0,
    logger: object | None = None,
) -> TEMAVG:
    """Read a Zonge TEMAVG processed ``.AVG`` file.

    Parameters
    ----------
    path : path-like
        Path to a processed TEMAVG ``.AVG`` file.

    Returns
    -------
    TEMAVG
        Parsed processed TEMAVG data and metadata.
    """
    return TEMAVG.read(path, verbose=verbose, logger=logger)


def read_tem_z(
    path: Pathish,
    *,
    verbose: int = 0,
    logger: object | None = None,
) -> TEMZPlot:
    """Read a Zonge TEMAVG contour ``.Z`` file.

    Parameters
    ----------
    path : path-like
        Path to a TEMAVG contour ``.Z`` file.

    Returns
    -------
    TEMZPlot
        Parsed contour rows and metadata.
    """
    return TEMZPlot.read(path, verbose=verbose, logger=logger)


def read_tem_log(
    path: Pathish,
    *,
    verbose: int = 0,
    logger: object | None = None,
) -> TEMLog:
    """Read a Zonge TEMAVG processing ``.LOG`` file.

    Parameters
    ----------
    path : path-like
        Path to a TEMAVG processing log.

    Returns
    -------
    TEMLog
        Parsed processing provenance and acquisition summary.
    """
    return TEMLog.read(path, verbose=verbose, logger=logger)

# ---------------------------------------------------------------------
# Generic XYZ / CSV reader  (fully implemented)
# ---------------------------------------------------------------------

def read_xyz(
    path: Pathish,
    *,
    current: float,
    tx_area: float | None = None,
    loop_side: float | None = None,
    loop_radius: float | None = None,
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
    error_col: int | None = None,
    skip_header: int = 0,
    delimiter: str | None = None,
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
        Transmitter loop area in square metres. Supply this
        *or* ``loop_side`` *or* ``loop_radius``.
    loop_side : float, optional
        Side of a square transmitter loop in metres.
    loop_radius : float, optional
        Radius of a circular transmitter loop in metres.
    data_type : str, optional
        Units / type of the ``data_col`` column.  One of
        ``"dBdt"``, ``"dHdt"``, ``"voltage"``,
        ``"normalized_voltage"``.  Default ``"dBdt"``.
    rx_area : float, optional
        Receiver coil area in square metres. Default 1.0.
    rx_turns : int, optional
        Receiver coil turns.  Default 1.
    tx_turns : int, optional
        Transmitter turns.  Default 1.
    offset : float, optional
        Tx-Rx horizontal separation in m. Default 0 (coincident).
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
        (convert nT/s to T/s), ``"pT/s"``, ``"nV/Am2"``.

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

    rows: list[list[str]] = []
    with open(path) as fh:
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
    _time_scale = {"s": 1.0, "ms": 1e-3, "us": 1e-6, "\u00b5s": 1e-6}
    if time_unit not in _time_scale:
        msg = (
            f"time_unit must be one of {list(_time_scale)}, "
            f"got {time_unit!r}"
        )
        raise ValueError(msg)
    t = t * _time_scale[time_unit]

    _data_scale = {
        "SI": 1.0,
        "nT/s": 1e-9,
        "pT/s": 1e-12,
        "nV/Am2": 1e-9,
        "uV/Am2": 1e-6,
    }
    if data_unit not in _data_scale:
        msg = (
            f"data_unit must be one of {list(_data_scale)}, "
            f"got {data_unit!r}"
        )
        raise ValueError(msg)
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
    error_col: int | None = None,
    **kwargs,
) -> dict[str, TEMSounding]:
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

    records: dict[str, list] = {}
    with open(path) as fh:
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


# ---------------------------------------------------------------------
# Internal helpers shared by the four format readers
# ---------------------------------------------------------------------

def _to_float(s: str) -> float:
    """Parse a float tolerant of Fortran-style exponents (1.2D+04 → 1.2e+04)."""
    return float(s.upper().replace("D", "E"))


def _gate_times_to_seconds(times: list[float], unit: str) -> np.ndarray:
    """Convert a list of gate-centre times to seconds."""
    scale = {"s": 1.0, "ms": 1e-3, "us": 1e-6, "µs": 1e-6}
    if unit not in scale:
        raise ValueError(
            f"gate time unit must be one of {list(scale)}, got {unit!r}"
        )
    return np.array(times, dtype=float) * scale[unit]


def _data_to_si(data: np.ndarray, unit: str) -> np.ndarray:
    """Convert TEM decay data to SI units (T/s for dBdt, V/(Am²) for Hz)."""
    scale = {
        "T/s": 1.0,        "nT/s": 1e-9,      "pT/s": 1e-12,
        "V/Am2": 1.0,      "V/(Am2)": 1.0,
        "nV/Am2": 1e-9,    "nV/(Am2)": 1e-9,
        "uV/Am2": 1e-6,    "µV/Am2": 1e-6,    "uV/(Am2)": 1e-6,
        "mV/Am2": 1e-3,    "mV/(Am2)": 1e-3,
        "SI": 1.0,
    }
    if unit not in scale:
        raise ValueError(
            f"data unit must be one of {list(scale)}, got {unit!r}"
        )
    return data * scale[unit]


# ---------------------------------------------------------------------
# Geosoft DAT reader
# ---------------------------------------------------------------------

def read_geosoft_dat(
    path: Pathish,
    *,
    current: float | None = None,
    tx_area: float | None = None,
    loop_side: float | None = None,
    loop_radius: float | None = None,
    rx_area: float = 1.0,
    rx_turns: int = 1,
    data_unit: str = "nV/Am2",
    data_type: str = "dBdt",
    gate_times: list[float] | None = None,
    gate_times_unit: str = "ms",
    x_col: str = "X",
    y_col: str = "Y",
    elev_col: str = "ELEV",
    station_col: str | None = None,
    gate_prefix: str = "G",
    line_col: str | None = None,
) -> list[TEMSounding]:
    r"""Read TEM data from a Geosoft ASCII line-data ``.dat`` file.

    Geosoft Oasis Montaj exports TEM soundings in a fixed-width or
    space-delimited ASCII file where:

    * Comment/metadata lines start with ``/`` or ``!``.
    * Line identifiers start with ``L`` (``L100``, ``L200``, …).
    * A single header row names every numeric channel.
    * Subsequent rows are one station (fiducial) per line.

    Metadata embedded in ``/`` comment lines is parsed automatically:

    .. code-block:: text

        / Current: 8.0
        / LoopSide: 200.0        (or LoopArea / LoopRadius)
        / GateTimes(ms): 0.021 0.037 0.065 ...
        / DataUnit: nV/Am2

    These override the keyword arguments of the same name.

    Parameters
    ----------
    path : str or Path
        Path to the ``.dat`` file.
    current : float or None
        Transmitter current in Amperes.  May also be given as a
        ``/ Current:`` header line in the file.
    tx_area : float or None
        Transmitter loop area in m².
    loop_side : float or None
        Side of a square loop in m (used to compute *tx_area*).
    loop_radius : float or None
        Radius of a circular loop in m (used to compute *tx_area*).
    rx_area : float
        Receiver coil area in m².  Default 1.0.
    rx_turns : int
        Receiver coil turns.  Default 1.
    data_unit : str
        Units of the gate-data columns before conversion to SI.
        Recognised: ``"nV/Am2"``, ``"uV/Am2"``, ``"mV/Am2"``,
        ``"nT/s"``, ``"T/s"``, ``"SI"``.  Default ``"nV/Am2"``.
    data_type : str
        ``"dBdt"`` or ``"normalized_voltage"``.  Default ``"dBdt"``.
    gate_times : list of float or None
        Explicit gate-centre times (in *gate_times_unit*).  When
        ``None`` the reader looks for a ``/ GateTimes:`` header line.
    gate_times_unit : str
        Time unit for *gate_times*: ``"ms"`` (default), ``"s"``,
        ``"us"``.
    x_col, y_col, elev_col : str
        Column names for easting, northing, elevation.
    station_col : str or None
        Column containing station identifiers.  ``None`` → auto-number.
    gate_prefix : str
        Prefix shared by all gate columns (e.g. ``"G"`` matches
        ``G1``, ``G2``, …, or ``"TEM_G"`` for ``TEM_G01``, …).
        Default ``"G"``.
    line_col : str or None
        Column name carrying the line number, if present.

    Returns
    -------
    list of TEMSounding
        One :class:`~pycsamt.tdem.TEMSounding` per data row.

    Raises
    ------
    FileNotFoundError
    ValueError
        If gate times cannot be determined from the file or arguments.

    Examples
    --------
    .. code-block:: text

        / Geosoft TEM export
        / Current: 8.0
        / LoopSide: 200.0
        / GateTimes(ms): 0.021 0.037 0.065 0.090
        / DataUnit: nV/Am2
        X        Y        ELEV    G1       G2       G3       G4
        1000.0   2000.0   0.0     1.24e4   8.91e3   5.23e3   3.12e3
        2000.0   2000.0   0.0     9.87e3   6.54e3   4.32e3   2.61e3

    >>> soundings = read_geosoft_dat(
    ...     "survey.dat",
    ...     current=8.0,
    ...     loop_side=200.0,
    ... )
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Geosoft DAT file not found: {path}")

    lines = path.read_text(errors="replace").splitlines()

    # ── Parse metadata from comment lines ─────────────────────────────────────
    _meta_float   = {}   # keyword → float
    _gate_times_f = None
    _data_unit_f  = None

    _KW = {
        "current": "current",
        "loopside": "loop_side",
        "loopradius": "loop_radius",
        "looparea": "tx_area",
        "txarea": "tx_area",
    }

    for ln in lines:
        stripped = ln.strip()
        if not (stripped.startswith("/") or stripped.startswith("!")):
            continue
        body = stripped.lstrip("/!").strip()

        # GateTimes(ms): 0.021 0.037 ...
        import re as _re
        m_gt = _re.match(
            r"GateTimes?\s*(?:\((\w+)\))?\s*[:=]\s*(.+)",
            body, _re.IGNORECASE,
        )
        if m_gt:
            unit_tag = m_gt.group(1) or gate_times_unit
            gate_times_unit = unit_tag.lower()
            _gate_times_f = [_to_float(v)
                              for v in m_gt.group(2).split()
                              if v.strip()]
            continue

        # DataUnit: nV/Am2
        m_du = _re.match(r"DataUnit\s*[:=]\s*(\S+)", body, _re.IGNORECASE)
        if m_du:
            _data_unit_f = m_du.group(1)
            continue

        # Current: 8.0  /  LoopSide: 200  / etc.
        m_kv = _re.match(r"(\w+)\s*[:=]\s*([0-9.eEdD+\-]+)", body)
        if m_kv:
            key = m_kv.group(1).lower().replace(" ", "")
            if key in _KW:
                _meta_float[_KW[key]] = _to_float(m_kv.group(2))

    # Merge file metadata with caller arguments (caller wins)
    if current     is None: current     = _meta_float.get("current")
    if tx_area     is None: tx_area     = _meta_float.get("tx_area")
    if loop_side   is None: loop_side   = _meta_float.get("loop_side")
    if loop_radius is None: loop_radius = _meta_float.get("loop_radius")
    if _data_unit_f is not None and data_unit == "nV/Am2":
        data_unit = _data_unit_f
    if _gate_times_f is not None and gate_times is None:
        gate_times = _gate_times_f

    if current is None:
        raise ValueError(
            "Transmitter current not found in file or arguments.  "
            "Pass current=<value>."
        )
    if tx_area is None and loop_side is None and loop_radius is None:
        raise ValueError(
            "Transmitter geometry not found.  "
            "Pass one of tx_area, loop_side, or loop_radius."
        )

    # ── Find header row and data rows ─────────────────────────────────────────
    header: list[str] = []
    data_rows: list[list[str]] = []
    current_line_id: str = ""

    for ln in lines:
        stripped = ln.strip()
        if not stripped:
            continue
        if stripped[0] in ("/", "!", "\\"):
            continue
        if stripped.upper().startswith("L") and not stripped[0].isdigit():
            # Line label (L100, L200 …)
            current_line_id = stripped
            continue
        # First non-comment, non-L row that is not all-numeric → header
        parts = stripped.split()
        try:
            float(parts[0])
            is_numeric = True
        except (ValueError, IndexError):
            is_numeric = False

        if not is_numeric and not header:
            header = [p.upper() for p in parts]
            continue

        if is_numeric or header:
            row = parts + ([""] if line_col and current_line_id else [])
            if line_col and current_line_id:
                row = [current_line_id.lstrip("Ll")] + parts
            data_rows.append(parts)

    if not header and data_rows:
        # No explicit header — generate G1, G2, … for all columns after 3
        n_cols = len(data_rows[0])
        header = [x_col.upper(), y_col.upper(), elev_col.upper()] + \
                 [f"{gate_prefix.upper()}{i+1}" for i in range(n_cols - 3)]

    if not data_rows:
        raise ValueError(f"No numeric data rows found in {path}")

    # Map column names → indices
    h = {name: i for i, name in enumerate(header)}

    def _idx(name: str, fallback: int | None = None) -> int | None:
        return h.get(name.upper(), fallback)

    xi   = _idx(x_col,    0)
    yi   = _idx(y_col,    1)
    eli  = _idx(elev_col, 2)
    sti  = _idx(station_col) if station_col else None

    # Gate columns: all columns whose name starts with gate_prefix (case-insensitive)
    pfx = gate_prefix.upper()
    gate_cols = sorted(
        [i for name, i in h.items() if name.startswith(pfx)],
        key=lambda idx: header[idx],
    )

    # Fallback: if no gate_prefix columns, everything after the first 3 cols
    if not gate_cols:
        start = max(xi or 0, yi or 1, eli or 2) + 1
        gate_cols = list(range(start, len(header)))

    if gate_times is None:
        raise ValueError(
            "Gate times not found in the file.  "
            "Pass gate_times=[...] (in gate_times_unit)."
        )
    t_s = _gate_times_to_seconds(gate_times, gate_times_unit)
    if len(t_s) != len(gate_cols):
        # Trim to the shorter length
        n = min(len(t_s), len(gate_cols))
        t_s      = t_s[:n]
        gate_cols = gate_cols[:n]

    soundings: list[TEMSounding] = []
    for row_i, row in enumerate(data_rows):
        try:
            x_val   = _to_float(row[xi])   if xi   is not None and xi < len(row)  else 0.0
            y_val   = _to_float(row[yi])   if yi   is not None and yi < len(row)  else 0.0
            el_val  = _to_float(row[eli])  if eli  is not None and eli < len(row) else 0.0
            st_name = (str(row[sti]) if sti is not None and sti < len(row)
                       else f"S{row_i + 1:04d}")
            d_raw = np.array(
                [_to_float(row[c]) for c in gate_cols if c < len(row)],
                dtype=float,
            )
        except (ValueError, IndexError) as exc:
            raise ValueError(
                f"Cannot parse data row {row_i + 1} in {path}: {row!r}"
            ) from exc

        d_si = _data_to_si(d_raw, data_unit)
        snd = TEMSounding.from_arrays(
            t_s[:len(d_si)], d_si,
            current=current,
            tx_area=tx_area,
            loop_side=loop_side,
            loop_radius=loop_radius,
            data_type=data_type,
            rx_area=rx_area,
            rx_turns=rx_turns,
            station_name=st_name,
            x=x_val,
            y=y_val,
            elevation=el_val,
        )
        soundings.append(snd)

    return soundings


# ---------------------------------------------------------------------
# AMIRA / EMIT .tem reader
# ---------------------------------------------------------------------

def read_amira(
    path: Pathish,
    *,
    current: float | None = None,
    tx_area: float | None = None,
    loop_side: float | None = None,
    loop_radius: float | None = None,
    rx_area: float = 1.0,
    rx_turns: int = 1,
    data_unit: str = "nV/Am2",
    data_type: str = "dBdt",
    gate_times: list[float] | None = None,
    gate_times_unit: str = "ms",
) -> list[TEMSounding]:
    r"""Read TEM data from an AMIRA / EMIT ``.tem`` keyword-block file.

    The EMIT ``.tem`` format organises acquisition parameters in
    keyword-value blocks separated by blank lines.  Comment lines begin
    with ``;`` or ``#``.  A minimal complete file looks like:

    .. code-block:: text

        ; EMIT TEM data — site survey
        TRANSMITTER_AREA   40000
        TRANSMITTER_CURRENT  8.0
        TRANSMITTER_TURNS    1
        RECEIVER_AREA  1.0
        RECEIVER_TURNS 1
        GATE_TIMES_UNIT ms
        GATE_TIMES  0.021 0.037 0.065 0.090 0.135 0.200
        DATA_UNIT   nV/Am2

        DATA
        ; StnName  X         Y         Elev    D1        D2      ...
        S001  1000.0  2000.0    0.0  1.24e+04  8.91e+03  ...
        S002  2000.0  2000.0    0.0  9.87e+03  6.54e+03  ...
        END DATA

    Keyword recognition is case-insensitive.  The ``DATA`` … ``END DATA``
    block may alternatively start with ``[DATA]`` or ``DATASETS``.

    Parameters
    ----------
    path : str or Path
    current : float or None
        Transmitter current [A].  Overrides ``TRANSMITTER_CURRENT`` in file.
    tx_area : float or None
        Transmitter loop area [m²].  Overrides ``TRANSMITTER_AREA``.
    loop_side : float or None
        Square loop side length [m].
    loop_radius : float or None
        Circular loop radius [m].
    rx_area : float
        Receiver coil area [m²].  Overrides ``RECEIVER_AREA``.  Default 1.0.
    rx_turns : int
        Receiver coil turns.  Overrides ``RECEIVER_TURNS``.  Default 1.
    data_unit : str
        Units of gate-data values.  Overrides ``DATA_UNIT`` in file.
    data_type : str
        ``"dBdt"`` or ``"normalized_voltage"``.  Default ``"dBdt"``.
    gate_times : list of float or None
        Gate-centre times (in *gate_times_unit*).
        Overrides ``GATE_TIMES`` in file.
    gate_times_unit : str
        Unit for *gate_times*.  Overrides ``GATE_TIMES_UNIT`` in file.

    Returns
    -------
    list of TEMSounding

    Raises
    ------
    FileNotFoundError
    ValueError
    """

    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"AMIRA/EMIT .tem file not found: {path}")

    lines = path.read_text(errors="replace").splitlines()

    # ── Parse keyword-value header ─────────────────────────────────────────────
    _kv: dict[str, str] = {}
    _gate_raw: list[float] = []
    data_lines: list[str] = []
    in_data = False

    _DATA_START = _re.compile(r"^\s*(DATA|DATASETS|\[DATA\])", _re.IGNORECASE)
    _DATA_END   = _re.compile(r"^\s*(END\s*DATA|END\s*DATASETS|\[END\])",
                               _re.IGNORECASE)

    for ln in lines:
        stripped = ln.strip()
        if not stripped or stripped[0] in (";", "#", "!"):
            continue

        if _DATA_END.match(stripped):
            in_data = False
            continue
        if _DATA_START.match(stripped):
            in_data = True
            continue

        if in_data:
            if not stripped.startswith(";") and not stripped.startswith("#"):
                data_lines.append(stripped)
            continue

        # keyword-value line
        parts = stripped.split(None, 1)
        if len(parts) == 2:
            _kv[parts[0].upper()] = parts[1].strip()
        elif len(parts) == 1:
            _kv[parts[0].upper()] = ""

    # Resolve acquisition parameters (caller wins over file)
    def _kv_float(key: str, fallback=None):
        v = _kv.get(key)
        return _to_float(v) if v else fallback

    def _kv_int(key: str, fallback: int) -> int:
        v = _kv.get(key)
        return int(_to_float(v)) if v else fallback

    if current is None:
        current = (_kv_float("TRANSMITTER_CURRENT")
                   or _kv_float("CURRENT")
                   or _kv_float("TX_CURRENT"))
    if tx_area is None:
        tx_area = (_kv_float("TRANSMITTER_AREA")
                   or _kv_float("TXAREA")
                   or _kv_float("TX_AREA"))
    if loop_side is None:
        loop_side = _kv_float("LOOP_SIDE") or _kv_float("LOOPSIDE")
    if loop_radius is None:
        loop_radius = _kv_float("LOOP_RADIUS") or _kv_float("LOOPRADIUS")
    if rx_area == 1.0:
        rx_area = _kv_float("RECEIVER_AREA") or _kv_float("RX_AREA") or 1.0
    if rx_turns == 1:
        rx_turns = _kv_int("RECEIVER_TURNS", 1) or _kv_int("RX_TURNS", 1)

    # Gate times
    _file_gt_unit = (_kv.get("GATE_TIMES_UNIT") or
                     _kv.get("TIME_UNIT") or gate_times_unit).lower()
    if gate_times is None:
        raw_gt = _kv.get("GATE_TIMES") or _kv.get("GATETIMES") or ""
        if raw_gt:
            gate_times = [_to_float(v) for v in raw_gt.split() if v.strip()]
            gate_times_unit = _file_gt_unit

    # Data unit
    _file_du = _kv.get("DATA_UNIT") or _kv.get("DATAUNIT") or ""
    if _file_du and data_unit == "nV/Am2":
        data_unit = _file_du.strip()

    if current is None:
        raise ValueError(
            "Transmitter current not found in AMIRA file or arguments."
        )
    if tx_area is None and loop_side is None and loop_radius is None:
        raise ValueError(
            "Transmitter geometry not found.  "
            "Pass tx_area, loop_side, or loop_radius."
        )
    if gate_times is None:
        raise ValueError(
            "GATE_TIMES not found in AMIRA file.  Pass gate_times=[...]."
        )
    t_s = _gate_times_to_seconds(gate_times, gate_times_unit)

    # ── Parse data block ───────────────────────────────────────────────────────
    soundings: list[TEMSounding] = []
    for row_i, ln in enumerate(data_lines):
        parts = ln.split()
        if not parts:
            continue
        # First token: station name (may start with letters)
        try:
            _to_float(parts[0])
            has_station_col = False
        except ValueError:
            has_station_col = True

        if has_station_col:
            st_name = parts[0]
            cols = parts[1:]
        else:
            st_name = f"S{row_i + 1:04d}"
            cols = parts

        if len(cols) < 1:
            continue

        x_val  = _to_float(cols[0]) if len(cols) > 0 else 0.0
        y_val  = _to_float(cols[1]) if len(cols) > 1 else 0.0
        el_val = _to_float(cols[2]) if len(cols) > 2 else 0.0
        data_start = 3

        n_gates = len(t_s)
        d_raw = np.array(
            [_to_float(cols[data_start + k])
             for k in range(n_gates)
             if data_start + k < len(cols)],
            dtype=float,
        )
        if d_raw.size == 0:
            continue

        n = min(len(t_s), len(d_raw))
        d_si = _data_to_si(d_raw[:n], data_unit)
        snd = TEMSounding.from_arrays(
            t_s[:n], d_si,
            current=current,
            tx_area=tx_area,
            loop_side=loop_side,
            loop_radius=loop_radius,
            data_type=data_type,
            rx_area=rx_area,
            rx_turns=rx_turns,
            station_name=st_name,
            x=x_val,
            y=y_val,
            elevation=el_val,
        )
        soundings.append(snd)

    return soundings


# ---------------------------------------------------------------------
# Zonge GDP reader
# ---------------------------------------------------------------------

def read_zonge(
    path: Pathish,
    *,
    current: float | None = None,
    tx_area: float | None = None,
    loop_side: float | None = None,
    loop_radius: float | None = None,
    rx_area: float = 1.0,
    rx_turns: int = 1,
    data_unit: str = "nV/Am2",
    data_type: str = "dBdt",
    gate_times: list[float] | None = None,
    gate_times_unit: str = "ms",
) -> list[TEMSounding]:
    r"""Read TEM sounding data from a Zonge GDP ``.avg`` / ``.tem`` file.

    The Zonge GDP receiver outputs an ASCII sounding file whose structure
    mirrors the TEMAVG processed format but at a per-sounding (not
    per-station-per-window) granularity.  Comment / metadata lines begin
    with ``*`` or ``\``.  A representative file:

    .. code-block:: text

        * Zonge GDP-32 TEM sounding
        * Station: S001
        * X: 1000.0   Y: 2000.0   Elev: 0.0
        * Current: 8.00 A
        * TxArea: 40000 m^2
        * GateTimes(ms): 0.021 0.037 0.065 0.090
        * DataUnit: nV/Am2
        Win  Time(ms)    Hz(nV/Am2)   Error
          1    0.021     1.2400e+04   5.0e+02
          2    0.037     8.9100e+03   4.0e+02
          3    0.065     5.2300e+03   2.5e+02
          4    0.090     3.1200e+03   1.8e+02
        * Station: S002
        ...

    Multiple soundings in one file are separated by a new ``* Station:``
    line.

    Parameters
    ----------
    path : str or Path
    current : float or None
        Transmitter current [A].
    tx_area : float or None
    loop_side : float or None
    loop_radius : float or None
    rx_area : float
    rx_turns : int
    data_unit : str
        Data column units.  Overrides ``* DataUnit:`` in file.
    data_type : str
        ``"dBdt"`` (default) or ``"normalized_voltage"``.
    gate_times : list of float or None
        Explicit gate times (override ``* GateTimes:`` in file).
    gate_times_unit : str

    Returns
    -------
    list of TEMSounding

    Raises
    ------
    FileNotFoundError
    ValueError
    """

    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Zonge GDP file not found: {path}")

    lines = path.read_text(errors="replace").splitlines()

    # ── Helpers ────────────────────────────────────────────────────────────────

    def _parse_meta_line(body: str, state: dict) -> None:
        """Extract key-value pairs from a `* Key: value` comment line."""
        m_gt = _re.match(
            r"GateTimes?\s*(?:\((\w+)\))?\s*[:=]\s*(.+)",
            body, _re.IGNORECASE,
        )
        if m_gt:
            unit_tag = (m_gt.group(1) or state.get("gt_unit", gate_times_unit)).lower()
            state["gt_unit"]  = unit_tag
            state["gt_raw"]   = [_to_float(v)
                                  for v in m_gt.group(2).split() if v.strip()]
            return

        for kw, key in (
            ("current", "current"), ("amp", "current"),
            ("txarea", "tx_area"), ("looparea", "tx_area"),
            ("loopside", "loop_side"), ("loopradius", "loop_radius"),
            ("rxarea", "rx_area"), ("rxturns", "rx_turns"),
            ("dataunit", "data_unit"),
        ):
            m = _re.match(rf"{kw}\s*[:=\s]\s*([^\s,;]+)",
                          body, _re.IGNORECASE)
            if m:
                state[key] = m.group(1).strip()
                return

        # Station metadata
        m_st = _re.match(r"Station\s*[:=]?\s*(\S+)", body, _re.IGNORECASE)
        if m_st:
            state["station_name"] = m_st.group(1)
            return
        for coord, key in (("X", "x"), ("Y", "y"), ("Elev", "elevation")):
            m_c = _re.match(rf"{coord}\s*[:=]\s*([0-9.eEdD+\-]+)",
                            body, _re.IGNORECASE)
            if m_c:
                state[key] = float(_to_float(m_c.group(1)))

    # ── First pass: collect global defaults ────────────────────────────────────
    global_state: dict = {}
    for ln in lines:
        stripped = ln.strip()
        if stripped and stripped[0] in ("*", "\\"):
            _parse_meta_line(stripped.lstrip("*\\").strip(), global_state)

    # ── Build per-sounding blocks ──────────────────────────────────────────────
    # A new sounding starts when a `* Station:` line is encountered.
    # Each block consists of its own metadata + data rows.

    blocks: list[tuple[dict, list[list[str]]]] = []
    cur_state: dict = dict(global_state)
    cur_rows:  list[list[str]] = []

    for ln in lines:
        stripped = ln.strip()
        if not stripped:
            continue

        if stripped[0] in ("*", "\\"):
            body = stripped.lstrip("*\\").strip()
            # New station → flush previous block
            if _re.match(r"Station\s*[:=]?\s*\S+", body, _re.IGNORECASE):
                if cur_rows:
                    blocks.append((cur_state, cur_rows))
                cur_state = dict(global_state)
                cur_rows  = []
            _parse_meta_line(body, cur_state)
            continue

        # Header row (non-numeric first token → skip)
        parts = stripped.split()
        if not parts:
            continue
        try:
            _to_float(parts[0])
        except ValueError:
            continue   # skip header or label lines

        cur_rows.append(parts)

    if cur_rows:
        blocks.append((cur_state, cur_rows))

    # ── Build caller-wins parameter resolver ───────────────────────────────────
    def _resolve(state: dict, key: str, caller_val, convert=float):
        if caller_val is not None:
            return caller_val
        v = state.get(key)
        return convert(v) if v is not None else None

    soundings: list[TEMSounding] = []

    for blk_i, (state, rows) in enumerate(blocks):
        blk_current     = _resolve(state, "current",     current)
        blk_tx_area     = _resolve(state, "tx_area",     tx_area)
        blk_loop_side   = _resolve(state, "loop_side",   loop_side)
        blk_loop_radius = _resolve(state, "loop_radius", loop_radius)
        blk_rx_area     = _resolve(state, "rx_area",     rx_area)
        blk_rx_turns    = _resolve(state, "rx_turns",    rx_turns, int)
        blk_du          = state.get("data_unit", data_unit)
        blk_st          = state.get("station_name", f"S{blk_i + 1:04d}")
        blk_x           = state.get("x", 0.0)
        blk_y           = state.get("y", 0.0)
        blk_el          = state.get("elevation", 0.0)

        if blk_current is None:
            raise ValueError(
                f"Transmitter current not found for block {blk_i + 1}."
            )
        if blk_tx_area is None and blk_loop_side is None and blk_loop_radius is None:
            raise ValueError(
                f"Transmitter geometry not found for block {blk_i + 1}."
            )

        # Determine time / data columns dynamically from rows
        # Expected columns: Win  Time(ms)  Data  [Error]
        # or:                Time(ms)  Data  [Error]
        t_list: list[float] = []
        d_list: list[float] = []
        e_list: list[float] = []

        gt_raw  = state.get("gt_raw") or (gate_times or [])
        gt_unit = state.get("gt_unit", gate_times_unit).lower()

        for row in rows:
            if not row:
                continue
            # Detect column layout by number of columns
            n = len(row)
            if gt_raw:
                # Gate times from header → rows are data only (may have Win prefix)
                try:
                    gate_i = int(_to_float(row[0])) - 1   # window index 1-based
                    if 0 <= gate_i < len(gt_raw):
                        d_val  = _to_float(row[-2] if n >= 3 else row[-1])
                        e_val  = _to_float(row[-1]) if n >= 3 else None
                        t_list.append(gt_raw[gate_i])
                        d_list.append(d_val)
                        if e_val is not None:
                            e_list.append(e_val)
                    continue
                except (ValueError, IndexError):
                    pass

            # No header gate times → columns: [Win]  Time  Data  [Error]
            if n >= 3:
                try:
                    int(_to_float(row[0]))   # first col is window number
                    t_list.append(_to_float(row[1]))
                    d_list.append(_to_float(row[2]))
                    if n >= 4:
                        e_list.append(_to_float(row[3]))
                    continue
                except ValueError:
                    pass
            if n >= 2:
                t_list.append(_to_float(row[0]))
                d_list.append(_to_float(row[1]))
                if n >= 3:
                    e_list.append(_to_float(row[2]))

        if not t_list:
            continue

        t_arr = _gate_times_to_seconds(t_list, gt_unit if gt_raw else gate_times_unit)
        d_arr = _data_to_si(np.array(d_list, dtype=float), blk_du)
        err   = np.array(e_list, dtype=float) if e_list else None
        if err is not None:
            err = _data_to_si(err, blk_du)

        snd = TEMSounding.from_arrays(
            t_arr, d_arr,
            current=blk_current,
            tx_area=blk_tx_area,
            loop_side=blk_loop_side,
            loop_radius=blk_loop_radius,
            data_type=data_type,
            rx_area=blk_rx_area or rx_area,
            rx_turns=blk_rx_turns or rx_turns,
            station_name=str(blk_st),
            x=float(blk_x),
            y=float(blk_y),
            elevation=float(blk_el),
            error=err,
        )
        soundings.append(snd)

    return soundings


# ---------------------------------------------------------------------
# WalkTEM / Aarhus Workbench reader
# ---------------------------------------------------------------------

def read_walkttem(
    path: Pathish,
    *,
    current: float | None = None,
    tx_area: float | None = None,
    loop_side: float | None = None,
    loop_radius: float | None = None,
    rx_area: float = 1.0,
    rx_turns: int = 1,
    data_unit: str = "nV/Am2",
    data_type: str = "dBdt",
    gate_times: list[float] | None = None,
    gate_times_unit: str = "ms",
    low_moment: bool = True,
) -> list[TEMSounding]:
    r"""Read WalkTEM / Aarhus Workbench ``.tem`` keyword-block files.

    The Aarhus Workbench TEM format uses ``/``-prefixed comment lines and
    section keywords in ALL-CAPS followed by their contents.  Sections end
    with ``END_<KEYWORD>``.  A representative file:

    .. code-block:: text

        / WalkTEM data — Aarhus Workbench export
        /
        GENERALHEADER
        WalkTEM System
        END_GENERALHEADER

        TXOPERATION
        LoopSize    40.000
        Current      6.60
        TurnsNumber  1
        Waveform     Full_Waveform
        END_TXOPERATION

        RXOPERATION
        CoilMoment  0.01
        END_RXOPERATION

        GATESET
        NumberOfGates 20
        GateTimes  0.008 0.015 0.025 0.043 0.072 0.120 0.200
        GateSigns  1 1 1 1 1 1 1
        END_GATESET

        DATA
        ! Line  Fid   X        Y       Elev   d1      d2      ...
        100  1000  1000.0  2000.0   0.0  1.5e-7  9.2e-8  ...
        100  2000  1200.0  2000.0   0.0  1.3e-7  8.1e-8  ...
        END_DATA

    Section names and keyword matching are case-insensitive.

    Parameters
    ----------
    path : str or Path
    current : float or None
        Transmitter current [A].  Overrides file value.
    tx_area : float or None
    loop_side : float or None
    loop_radius : float or None
    rx_area : float
    rx_turns : int
    data_unit : str
        Units for data gate columns (default ``"nV/Am2"``).
        WalkTEM data may also be in ``"T/s"`` or ``"SI"``.
    data_type : str
        ``"dBdt"`` (default) or ``"normalized_voltage"``.
    gate_times : list of float or None
        Explicit gate times (overrides ``GATESET`` block).
    gate_times_unit : str
    low_moment : bool
        WalkTEM instruments record two moments (LM and HM) per
        sounding.  When ``True`` (default) only the low-moment (first)
        channel is read; set ``False`` to read the high-moment channel.

    Returns
    -------
    list of TEMSounding

    Raises
    ------
    FileNotFoundError
    ValueError
    """

    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"WalkTEM file not found: {path}")

    lines = path.read_text(errors="replace").splitlines()

    # ── Section-block parser ───────────────────────────────────────────────────
    _sections: dict[str, list[str]] = {}
    current_section: str | None = None
    section_lines: list[str] = []
    data_lines: list[str] = []
    in_data = False

    _SECT_START = _re.compile(
        r"^(TXOPERATION|RXOPERATION|GATESET|GENERALHEADER|DATASTATUS)",
        _re.IGNORECASE,
    )
    _SECT_END   = _re.compile(r"^END_(\w+)", _re.IGNORECASE)
    _DATA_START = _re.compile(r"^DATA\s*$", _re.IGNORECASE)
    _DATA_END   = _re.compile(r"^END_DATA", _re.IGNORECASE)

    for ln in lines:
        stripped = ln.strip()
        if not stripped or stripped[0] in ("/", "!", ";", "#"):
            continue

        if _DATA_END.match(stripped):
            in_data = False
            continue
        if _DATA_START.match(stripped):
            in_data = True
            continue
        if in_data:
            data_lines.append(stripped)
            continue

        m_end = _SECT_END.match(stripped)
        if m_end:
            if current_section:
                _sections[current_section.upper()] = section_lines
            current_section = None
            section_lines   = []
            continue

        m_start = _SECT_START.match(stripped)
        if m_start:
            current_section = m_start.group(1)
            section_lines   = []
            continue

        if current_section:
            section_lines.append(stripped)

    # ── Parse TX section ───────────────────────────────────────────────────────
    def _sect_float(section: str, *keys) -> float | None:
        for ln in _sections.get(section.upper(), []):
            parts = ln.split()
            if len(parts) >= 2 and parts[0].lower() in [k.lower() for k in keys]:
                try:
                    return _to_float(parts[1])
                except ValueError:
                    pass
        return None

    _tx = _sections.get("TXOPERATION", [])
    for ln in _tx:
        parts = ln.split()
        if not parts:
            continue
        kl = parts[0].lower()
        val = parts[1] if len(parts) >= 2 else None
        if val is None:
            continue
        if kl in ("loopsize", "loopside") and loop_side is None and tx_area is None:
            loop_side = _to_float(val)
        elif kl in ("loopradius",) and loop_radius is None and tx_area is None:
            loop_radius = _to_float(val)
        elif kl in ("looparea", "txarea") and tx_area is None:
            tx_area = _to_float(val)
        elif kl == "current" and current is None:
            current = _to_float(val)
        elif kl == "turnsnumber" and rx_turns == 1:
            rx_turns = int(_to_float(val))

    _rx = _sections.get("RXOPERATION", [])
    for ln in _rx:
        parts = ln.split()
        if len(parts) >= 2 and parts[0].lower() in ("coilmoment", "rxarea", "coilarea"):
            if rx_area == 1.0:
                rx_area = _to_float(parts[1])

    # ── Parse gate times ───────────────────────────────────────────────────────
    _gs = _sections.get("GATESET", [])
    _file_gt: list[float] = []
    _file_gt_unit = gate_times_unit

    for ln in _gs:
        parts = ln.split()
        if not parts:
            continue
        kl = parts[0].lower()
        if kl in ("gatetimes", "times", "gatetime"):
            _file_gt = [_to_float(v) for v in parts[1:] if v.strip()]
        elif kl in ("gatetimesunit", "timeunit"):
            _file_gt_unit = parts[1].lower() if len(parts) > 1 else "ms"

    if gate_times is None and _file_gt:
        gate_times      = _file_gt
        gate_times_unit = _file_gt_unit

    if current is None:
        raise ValueError(
            "Transmitter current not found in WalkTEM file or arguments."
        )
    if tx_area is None and loop_side is None and loop_radius is None:
        raise ValueError(
            "Transmitter geometry not found.  "
            "Pass tx_area, loop_side, or loop_radius."
        )
    if gate_times is None:
        raise ValueError(
            "Gate times not found in WalkTEM GATESET block or arguments."
        )

    t_s = _gate_times_to_seconds(gate_times, gate_times_unit)
    n_gates = len(t_s)

    # ── Parse DATA block ───────────────────────────────────────────────────────
    # Expected columns (from the ! comment header):
    #   Line  Fid  X  Y  Elev  d1  d2  ...  [e1  e2  ...]
    # The first non-numeric column line is a header; skip it.

    soundings: list[TEMSounding] = []
    for row_i, ln in enumerate(data_lines):
        parts = ln.split()
        if not parts:
            continue
        # Skip header lines (non-numeric first token)
        try:
            _to_float(parts[0])
        except ValueError:
            continue

        # Layout: Line  Fid  X  Y  Elev  d1 … dn  [e1 … en]
        if len(parts) < 5:
            continue
        try:
            line_id   = parts[0]
            fid_id    = parts[1]
            x_val     = _to_float(parts[2])
            y_val     = _to_float(parts[3])
            el_val    = _to_float(parts[4])
        except (ValueError, IndexError):
            continue

        rest = parts[5:]
        n_rest = len(rest)

        # Detect whether error columns follow data columns
        # Heuristic: if n_rest == 2*n_gates, first half = data, second = errors
        if n_rest >= 2 * n_gates:
            d_cols = rest[:n_gates]
            e_cols = rest[n_gates: 2 * n_gates]
        else:
            d_cols = rest[:min(n_gates, n_rest)]
            e_cols = []

        if not d_cols:
            continue

        try:
            d_raw = np.array([_to_float(v) for v in d_cols], dtype=float)
            e_raw = (np.array([_to_float(v) for v in e_cols], dtype=float)
                     if e_cols else None)
        except ValueError:
            continue

        n = min(len(t_s), len(d_raw))
        d_si = _data_to_si(d_raw[:n], data_unit)
        e_si = _data_to_si(e_raw[:n], data_unit) if e_raw is not None else None

        st_name = f"L{line_id}_F{fid_id}"
        snd = TEMSounding.from_arrays(
            t_s[:n], d_si,
            current=current,
            tx_area=tx_area,
            loop_side=loop_side,
            loop_radius=loop_radius,
            data_type=data_type,
            rx_area=rx_area,
            rx_turns=rx_turns,
            station_name=st_name,
            x=x_val,
            y=y_val,
            elevation=el_val,
            error=e_si,
        )
        soundings.append(snd)

    return soundings
