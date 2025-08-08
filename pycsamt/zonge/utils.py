# -*- coding: utf-8 -*-
#       Author: LKouadio <etanoyau@gmail.com>
#       License: LGPL-3.0-or-later
"""pycsamt.zonge.utils
General‑purpose helpers for **Zonge** AVG / AMTAVG files and
accompanying *station* profiles.

The file now supports

* **Kind‑1** legacy whitespace tables.
* **Kind‑2** comma‑separated tables with leading metadata.

A tidy :class:`pandas.DataFrame` plus a *metadata* dict is returned
regardless of flavour.  Column names are normalised to a concise lower‑
case schema (``station, freq, emag, rho, phase, …``).
"""
from __future__ import annotations

from datetime import datetime 
from numbers import Integral
from pathlib import Path
import re
import io
from typing import (
    List, 
    Tuple, 
    Dict, 
    Optional, 
    Sequence,
    Any, 
    Iterable
)

import numpy as np
import pandas as pd

from ..gis.utils import ll_to_utm                # type: ignore
from ..exceptions import (
    AvgFileError, 
    AvgDataError, 
  )
from ..log.logger import get_logger


__all__ = [
    "load_avg", 
    "round_dipole_length", 
    "validate_stn_profile",
    "classify_avg_format", 
    "extract_core_columns",
    "number_stations", 
    "chunk_by_frequency"
]

logger = get_logger(__name__)


_RX_WS             = re.compile(r"\s+")
_RX_K2_HEADER      = re.compile(r"^Z\.mwgt\s*,", re.I)
_RX_K1_HEADER      = re.compile(r"^skp\s+Station", re.I)
_NUMERIC_REPLACE   = {"*": np.nan, "nan": np.nan, "NaN": np.nan,
                      "": np.nan}
_COL_MAP = {
    'Station': 'station', 'Stn': 'station', 'skp': 'station',
    'Freq': 'freq', 'Freq.': 'freq',
    'Comp': 'comp',
    'Tx.Amp': 'amps', 'Amps': 'amps',
    'E.mag': 'emag', 'Emag': 'emag',
    'E.phz': 'ephz', 'Ephz': 'ephz',
    'H.mag': 'hmag', 'B.mag': 'hmag',
    'Hphz':  'hphz', 'B.phz': 'hphz',
    'Resistivity': 'rho', 'ARes.mag': 'rho',
    'Phase': 'phase', 'Z.phz': 'phase',
}


_CANON2K2: dict[str, str] = {
    # survey logistics
    "station": "Station",
    "freq":    "Freq",
    "comp":    "Comp",
    # transmitter / field data
    "amps":  "Tx.Amp",
    "emag":  "E.mag",
    "ephz":  "E.phz",
    "hmag":  "B.mag",
    "hphz":  "B.phz",
    "zmag":  "Z.mag",
    "phase": "Z.phz",
    "rho":   "ARes.mag",
    # quality metrics
    "e.%err":   "E.%err",
    "e.perr":   "E.perr",
    "h.%err":   "B.%err",
    "h.perr":   "B.perr",
    "rho.%err": "ARes.%err",
    "phase.%err": "Z.perr",
    "z.%err":   "Z.%err",
    "z.perr":   "Z.perr",
    # ASTATIC weights (keep canonical labels)
    "z.mwgt": "Z.mwgt",
    "z.pwgt": "Z.pwgt",
    "e.wgt":  "E.wgt",
    "b.wgt":  "B.wgt",
}


def _to_float(val: str | float | int) -> float | np.floating:
    """Convert *val* to float while honouring project placeholders."""
    if isinstance(val, (float, int)):
        return float(val)
    txt = val.strip()
    if txt in _NUMERIC_REPLACE:
        return np.nan
    if txt.startswith('.'):
        txt = '0' + txt
    if txt.endswith('.'):
        txt = txt + '0'
    try:
        return float(txt)
    except ValueError:
        return np.nan

def classify_avg_format(lines: Sequence[str]) -> int:
    """Return **1** or **2** depending on AVG flavour.

    Parameters
    ----------
    lines : Sequence[str]
        Raw text lines read from the AVG file.

    Raises
    ------
    AvgFileError
        When the function cannot detect a valid header.
    """
    for ln in lines:
        if _RX_K2_HEADER.search(ln):
            logger.debug("Kind‑2 AVG detected")
            return 2
        if _RX_K1_HEADER.search(ln):
            logger.debug("Kind‑1 AVG detected")
            return 1
    raise AvgFileError("Unrecognised AVG header – cannot classify file")


def _parse_kind1(lines: Sequence[str]) -> pd.DataFrame:
    """Parse legacy fixed‑width (kind‑1) AVG table."""
    idx = next((i for i, ln in enumerate(lines)
                if _RX_K1_HEADER.search(ln)), None)
    if idx is None:
        raise AvgFileError("Header row not found in kind‑1 file")

    hdr_tokens = _RX_WS.sub(' ', lines[idx].strip()).split()
    data_rows: List[List[Any]] = []
    for ln in lines[idx + 1:]:
        if not ln.strip() or _RX_K1_HEADER.search(ln):
            break
        if ln.startswith(('\\', '$')):
            continue
        tokens = _RX_WS.sub(' ', ln.strip()).split()
        data_rows.append([_to_float(tk) if j >= 4 else tk
                          for j, tk in enumerate(tokens)])
    if not data_rows:
        raise AvgDataError("No data rows in kind‑1 file")
    df = pd.DataFrame(data_rows, columns=hdr_tokens)
    return _standardise_columns(df)


def _parse_kind2(
        lines: Sequence[str]) -> Tuple[pd.DataFrame, Dict[str, str]]:
    """Parse comma‑separated (kind‑2) AVG file.

    Returns a pair *(df, metadata)*.
    """
    meta: Dict[str, str] = {}
    block: List[str]      = []
    for ln in lines:
        if ln.startswith('$') and '=' in ln:
            key, val           = ln[1:].split('=', 1)
            meta[key.strip()]  = val.strip()
            continue
        if _RX_K2_HEADER.match(ln):
            block.append(ln)
            continue
        if block and ln.strip():
            block.append(ln)
        if block and not ln.strip():
            break
    if not block:
        raise AvgDataError("Data block missing in kind‑2 file")
    df = pd.read_csv(io.StringIO('\n'.join(block)), skipinitialspace=True)
    df = df.applymap(lambda v: _to_float(v) if isinstance(v, str) else v)
    return _standardise_columns(df), meta


def write_avg(
    core:  pd.DataFrame,
    extra: pd.DataFrame | None,
    meta:  dict[str, str] | None,
    path:  str | Path | None = None,     # ← path now optional
    *,
    stamp: bool = True,
    float_fmt: str = "%.6g",
    na_rep: str = "",
) -> Path:
    """
    Serialise *core* + *extra* to **kind-2 CSV**.

    Parameters
    ----------
    core, extra
        Frames returned by :func:`extract_core_columns`.
    meta
        Header ``$key = value`` pairs.
    path
        Destination.  When *None* a file called
        ``exported_kind2.avg`` is created in the **current working
        directory**.
    stamp
        Append a ``$Written`` UTC time-stamp.
    float_fmt, na_rep
        Forwarded to :pymeth:`pandas.DataFrame.to_csv`.

    Returns
    -------
    pathlib.Path
        Absolute path to the written file.
    """

    # 1) Resolve destination
    if path is None:
        path = Path.cwd() / "exported_kind2.avg"
    path = Path(path).expanduser().resolve()

    # 2) Build header block
    header: list[str] = []
    if meta:
        header.extend(f"${k} = {v}" for k, v in meta.items())
    if stamp:
        utc = datetime.utcnow().isoformat(timespec="seconds") + "Z"
        header.append(f"$Written = {utc}")
    header.append("")                       # blank line before table

    # 3) Assemble dataframe
    block = pd.concat([core, extra], axis=1) if extra is not None else core
    block = block.copy()

    # restore canonical → kind-2 column case
    rename_map = {c: _CANON2K2.get(c.lower(), c) for c in block.columns}
    block.rename(columns=rename_map, inplace=True)

    # 4) Dump to disk
    csv_txt = block.to_csv(index=False, float_format=float_fmt, na_rep=na_rep)
    path.write_text("\n".join(header) + csv_txt, encoding="utf8")

    logger.info("AVG written → %s", path)
    return path

def _standardise_columns(df: pd.DataFrame) -> pd.DataFrame:
    return df.rename(columns={c: _COL_MAP.get(c, c.lower()) for c in df.columns})


def load_avg(
    path: str | Path, *,
    ll_columns: Tuple[str, str] = ('latitude', 'longitude'),
    utm_zone:   Optional[int]   = None,
    inplace:    bool            = False
    ) -> Tuple[pd.DataFrame, Dict[str, str]]:
    """Read a Zonge AVG file.

    Parameters
    ----------
    path         : str | Path
        Filesystem location of the ``.avg`` file.
    ll_columns   : tuple(str, str), optional
        Column names holding lat / lon values if present.
    utm_zone     : int, optional
        Override automatic UTM zone detection.
    inplace      : bool, default *False*
        If *False* a copy is made before adding derived columns.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(path)
    lines = path.read_text(errors='replace').splitlines()
    kind  = classify_avg_format(lines)
    if kind == 1:
        df, meta = _parse_kind1(lines), {}
    else:
        df, meta = _parse_kind2(lines)
    if not inplace:
        df = df.copy()
    lat, lon = ll_columns
    if lat in df.columns and lon in df.columns:
        try:
            east, north    = ll_to_utm(df[lat].values, df[lon].values,
                                       zone=utm_zone)
            df['easting']  = east
            df['northing'] = north
        except Exception as exc:  # pragma: no cover
            logger.warning("Lat/Lon → UTM failed: %s", exc)
    return df, meta


# Station‑profile utilities
def validate_stn_profile(
    profile: Sequence[str],
    splitter: str | None = None
    ) -> Tuple[int, List[Tuple[str, int]]]:
    """Check *STN* file header and return *(score, index_list)*.

    Parameters
    ----------
    profile  : Sequence[str]
        Raw lines from the ``.stn`` file.
    splitter : str | None, default *None*
        Token delimiter.  Defaults to whitespace.
    """
    splitter        = splitter or ' '
    labels: set[str] = {
        'dot', 'station', 'sta',
        'e', 'east', 'easting',
        'n', 'north', 'northing',
        'h', 'elev', 'elevation',
        'lon', 'lat', 'utm_zone',
    }
    header_tokens   = [tk.strip().lower() for tk in
                       profile[0].split(splitter) if tk]
    matches         = [(tk, idx) for idx, tk in enumerate(header_tokens)
                       if tk in labels]
    return len(matches), matches


# Misc helpers (kept for backward compat.)

def round_dipole_length(length: float | int) -> float:
    """Round *length* to the nearest 5‑m increment."""
    length = float(length)
    mod    = length % 5
    if mod < 3:
        return length - mod
    if mod < 7:
        return length - mod + 5
    return length - mod + 10


# NEW v2 helpers
def extract_core_columns(
     df: pd.DataFrame,
     *,
    keep: Iterable[str] | None = None
    ) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Split a kind‑2 frame into **core** and **extra** columns.

    Parameters
    ----------
    df   : pandas.DataFrame
        Parsed DataFrame from :func:`load_avg`.
    keep : Iterable[str] | None, optional
        Canonical column names to retain.  Defaults to a
        built‑in minimal set.

    Returns
    -------
    core  : DataFrame
        Only the requested *keep* columns (plus *station* if absent).
    extra : DataFrame
        All remaining columns.
    """
    default_keep = {
        'station', 'freq', 'amps', 'emag', 'ephz',
        'hmag', 'hphz', 'rho', 'phase',
        'e.%err', 'e.perr', 'h.%err', 'h.perr',
        'rho.%err', 'phase.%err',
    }
    keep_set      = set(k.lower() for k in (keep or default_keep))
    cols_lower    = {c.lower(): c for c in df.columns}
    have_keep     = [cols_lower[c] for c in keep_set if c in cols_lower]
    if 'station' not in [c.lower() for c in have_keep]:
        have_keep.insert(0, cols_lower.get('station', df.columns[0]))
    core          = df[have_keep].copy()
    extra         = df.drop(columns=have_keep).copy()
    return core.reset_index(drop=True), extra.reset_index(drop=True)


def _block_to_dict(block: Sequence[str]) -> Dict[str, Any]:
    """Convert ``key=value`` lines to a dict (case‑insensitive keys)."""
    out: Dict[str, Any] = {}
    for ln in block:
        if '=' not in ln:
            continue
        k, v = ln.split('=', 1)
        out[k.strip().lower()] = v.strip()
    return out


def _dict_to_lines(data: Any) -> List[str]:
    """Serialise dict‑like *data* to ``key=value`` text lines."""
    if isinstance(data, str):
        import json
        data = json.loads(data)
    if not isinstance(data, dict):
        data = dict(data)
    return [f"{k}={v}\n" for k, v in data.items() if v is not None]


def number_stations(
        n_stations: int | Integral,
        n_freq:     int | Integral,
        *,
        prefix: str = "S"
    ) -> Tuple[List[str], List[str]]:
    """
    Generate station IDs and a frequency-expanded copy.

    Parameters
    ----------
    n_stations, n_freq
        Positive integers.  ``n_freq`` is the number of
        frequencies **per** station.
    prefix
        String prepended to every station label.

    Returns
    -------
    ids
        ``['S00', 'S01', …]`` up to ``n_stations – 1``.
    expanded
        Each ID repeated ``n_freq`` times (ordered like the
        original table: *all* rows for S00, then S01, …).
    """
    if n_stations < 1 or n_freq < 1:
        raise ValueError("n_stations and n_freq must be ≥ 1")
    ids       = [f"{prefix}{i:02d}" for i in range(int(n_stations))]
    expanded  = list(np.repeat(ids, int(n_freq)))
    return ids, expanded


def chunk_by_frequency(
    data: Sequence[Any] | np.ndarray,
    n_per_chunk: int | Integral,
    *,
    drop_remainder: bool = False
) -> List[np.ndarray]:
    """
    Split *data* into equally sized chunks (one per frequency).

    Parameters
    ----------
    data
        Any 1-D array-like object.  A copy is **not** made unless
        conversion is required.
    n_per_chunk
        Items per chunk (e.g. number of stations for that
        frequency).
    drop_remainder
        If *True*, discard a final partial chunk; otherwise it is
        returned as-is.

    Returns
    -------
    chunks
        List of ``numpy.ndarray`` slices.

    Examples
    --------
    >>> chunk_by_frequency([0, 1, 2, 3, 4], 2)
    [array([0, 1]), array([2, 3]), array([4])]
    """
    if n_per_chunk < 1:
        raise ValueError("n_per_chunk must be ≥ 1")
    arr   = np.asarray(data)
    total = arr.size
    idx   = np.arange(0, total, int(n_per_chunk))
    chunks: List[np.ndarray] = [
        arr[i : i + n_per_chunk] for i in idx
    ]
    if drop_remainder and chunks and chunks[-1].size < n_per_chunk:
        chunks.pop()
    return chunks


