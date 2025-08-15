# -*- coding: utf-8 -*-
#       Author: LKouadio <etanoyau@gmail.com>
#       License: LGPL-3.0 
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
import warnings
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
try:
    import xarray as xr  # type: ignore
except ImportError:  # pragma: no cover
    warnings.warn(
        "xarray is required for the package"
    )
from ..decorators import isdf 
from ..gis.utils import ll_to_utm # type: ignore
from ..exceptions import (
    AvgFileError, 
    AvgDataError, 
  )
from ..log.logger import get_logger
from ..utils.deps import ensure_pkg 

__all__ = [
    "load_avg", 
    "round_dipole_length", 
    "validate_stn_profile",
    "classify_avg_format", 
    "extract_core_columns",
    "number_stations", 
    "chunk_by_frequency", 
    "write_avg", 
    "to_xarray"
]

logger = get_logger(__name__)


_RX_WS             = re.compile(r"\s+")
# _RX_K2_HEADER      = re.compile(r"^Z\.mwgt\s*,", re.I)
# _RX_K1_HEADER      = re.compile(r"^skp\s+Station", re.I)

_RX_K2_HEADER = re.compile(r"^\s*Z\.mwgt\s*,", re.I)
_RX_K1_HEADER = re.compile(r"^\s*skp\s+Station", re.I)

_NUMERIC_REPLACE   = {"*": np.nan, "nan": np.nan, "NaN": np.nan,
                      "": np.nan}
_COMMENT_PREFIXES = ('\\', '/', '!', '"')

_COL_MAP = {
    # General & Legacy
    'Station': 'station', 
    'Stn': 'station',
    'skp': 'skp',
    'Freq': 'freq', 
    'Freq.': 'freq',
    'Comp': 'comp',
    'Amps': 'amps',
    'Emag': 'emag',
    'Ephz': 'ephz',
    'Hmag': 'hmag',
    'Hphz': 'hphz',
    'Resistivity': 'rho',
    'Phase': 'phase',
    'TMARES': 'rho_sc', 
    'SRES': 'rho_sc', 
    'TMARES/SRES': 'rho_sc',

    # Modern (CSAVGW) - overlaps are fine
    'Tx.Amp': 'amps',
    'E.mag': 'emag',
    'E.phz': 'ephz',
    'B.mag': 'hmag',  # Maps to hmag
    'B.phz': 'hphz',  # Maps to hphz
    'H.mag': 'hmag',
    'Z.mag': 'zmag',
    'Z.phz': 'phase',
    'ARes.mag': 'rho',
    'SRes': 'rho_sc',
    'E.%err': 'e.%err',
    'E.perr': 'e.perr',
    'B.%err': 'h.%err', # Maps to h.%err
    'B.perr': 'h.perr', # Maps to h.perr
    'Z.%err': 'z.%err',
    'Z.perr': 'z.perr',
    'ARes.%err': 'rho.%err',
    'E.wgt': 'e.wgt',
    'H.wgt': 'h.wgt',
    'Choer': 'coh',
    'Gdp.Blk': 'gdp_blk',
    'Gdp.Chn': 'gdp_chn',
    'Gdp.Time': 'gdp_time',
    '|Z|': 'zabs',
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
    """
    Return **1** or **2** depending on AVG flavour.

    This function uses a multi-pass approach for robustness. It
    first checks for explicit headers and then falls back to
    analyzing structural clues like keyword format and data
    delimiters.
    
    Parameters
    ----------
    lines : Sequence[str]
        Raw text lines read from the AVG file.

    Raises
    ------
    AvgFileError
        When the function cannot detect a valid header.
        
    """
    # First pass: Look for the most definitive header lines.
    # This is the fastest and most reliable method.
    for ln in lines:
        if _RX_K2_HEADER.search(ln):
            logger.debug(
                "Kind-2 AVG detected (found modern header)."
            )
            return 2
        if _RX_K1_HEADER.search(ln):
            logger.debug(
                "Kind-1 AVG detected (found legacy header)."
            )
            return 1

    # Second pass: If no header found, analyze file structure.
    has_commas = False
    has_dot_in_keyword = False

    for ln in lines:
        s = ln.strip()
        if not s:
            continue

        # Check for modern keyword style, e.g., '$Survey.Type='
        if s.startswith('$'):
            # This is a very strong indicator of a modern file.
            match = re.match(r"\$\s*(\w+\.\w+)\s*=", s)
            if match:
                has_dot_in_keyword = True
                break

        # Heuristic: check for comma-separated data.
        # Avoids matching headers with only one or two commas.
        if not s.startswith(('$', '\\', '/')) and s.count(',') > 2:
            has_commas = True

    if has_dot_in_keyword:
        logger.debug(
            "Kind-2 AVG detected (found dot-notation keywords)."
        )
        return 2

    if has_commas:
        logger.debug(
            "Kind-2 AVG detected (found comma-separated data)."
        )
        return 2

    # If no strong indicators are found, parsing cannot proceed.
    raise AvgFileError(
        "Unrecognised AVG header – cannot classify file"
    )

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
        data_rows.append([
            _to_float(tk) if j >= 4 else tk  # First 4 cols are non-numeric
            for j, tk in enumerate(tokens)
        ])
    if not data_rows:
        raise AvgDataError("No data rows in kind‑1 file")
    df = pd.DataFrame(data_rows, columns=hdr_tokens)
    return _standardise_columns(df)

def _is_comment(ln: str) -> bool:
    return bool(ln) and ln[0] in _COMMENT_PREFIXES


def _next_block(lines, i):
    """Yield (header_line_index, end_index_exclusive)."""
    n = len(lines)
    while i < n and not _RX_K2_HEADER.match(lines[i]):
        i += 1
    if i >= n:
        return None, n
    start = i
    i += 1
    while i < n:
        s = lines[i]
        if not s.strip():
            i += 1
            break
        # DO NOT break on comments — they'll be skipped later
        if _RX_K2_HEADER.match(s) or s.startswith('$'):
            break
        i += 1
    return (start, i)

def _parse_kind2(
    lines: Sequence[str]
) -> Tuple[pd.DataFrame, Dict[str, str]]:
    """
    Parse a modern CSAVGW (kind-2) AVG file that contains
    repeated CSV blocks, typically one per station/component.

    The function:
      * collects top-level $meta (once, before the first table),
      * collects per-block $meta (stamped onto the block),
      * reads each CSV block until a blank line, a new header,
        or the next $keyword line,
      * normalises column names to canonical lower-case,
      * stamps a 'station' column from $Rx.Stn (fallbacks to
        $Rx.GdpStn or $Stn.Beg),
      * returns a single tidy DataFrame plus a meta dict that
        includes a 'blocks' list holding each block's $meta.

    Parameters
    ----------
    lines : Sequence[str]
        Raw text lines read from the .avg file.

    Returns
    -------
    (df, meta) : (pandas.DataFrame, dict)
        df    : tidy concatenation of all blocks.
        meta  : top-level $meta and 'blocks' per-block $meta.

    Raises
    ------
    AvgDataError
        If no data blocks are found.
    """
    global_meta: Dict[str, str] = {}
    blocks_meta: List[Dict[str, str]] = []
    frames: List[pd.DataFrame] = []

    block_meta: Dict[str, str] = {}
    i, n = 0, len(lines)
    seen_table = False

    while i < n:
        ln = lines[i]

        # Skip any comment lines anywhere in the file.  CSAVGW
        # allows \, /, !, " to start comments.
        if _is_comment(ln):
            i += 1
            continue

        # Collect $key=value lines.  Before the first table, the
        # keys are also considered "global" survey/job config.
        if ln.startswith('$') and '=' in ln:
            key, val = ln[1:].split('=', 1)
            key = key.strip()
            val = val.strip()

            if not seen_table:
                global_meta[key] = val

            # Block-level $meta applies to the next table we hit.
            block_meta[key] = val
            i += 1
            continue

        # 3) If this line is the CSV header, parse *this* block
        if _RX_K2_HEADER.match(ln):
            # If we're here, try to locate the next table.  This can
            # be separated by blanks or interleaved comments/metadata.
            start, j = _next_block(lines, i)
            if start is None:
                break
    
            seen_table = True
    
            # Assemble header + rows for this block while skipping
            # inline comment lines that may appear among rows.
            table_txt = '\n'.join(
                [lines[start]] + [
                    s for s in lines[start + 1:j]
                    if not _is_comment(s)
                ]
            )
    
            # Parse CSV with forgiving whitespace.  Convert numeric
            # strings (including '*', '.5', '1.') with _to_float.
            dfb = pd.read_csv(
                io.StringIO(table_txt),
                skipinitialspace=True
            )
            dfb = dfb.applymap(
                lambda v: _to_float(v) if isinstance(v, str) else v
            )
    
            # Stamp station and a few helpful block-level fields as
            # columns.  Prefer client station number ($Rx.Stn).
            stn = (
                block_meta.get('Rx.Stn')
                or block_meta.get('Rx.GdpStn')
                or block_meta.get('Stn.Beg')
            )
            if stn is not None:
                try:
                    dfb['station'] = _to_float(stn)
                except Exception:
                    dfb['station'] = stn  # keep as text if odd
    
            # Component label and a couple of helpers can be handy for
            # QC.  They are optional and harmless if missing.
            if 'Rx.Cmp' in block_meta:
                dfb['comp'] = block_meta['Rx.Cmp']
            for k in ('Rx.Length', 'Rx.GdpStn'):
                if k in block_meta:
                    dfb[k.replace('.', '_').lower()] = block_meta[k]
    
            # Standardise to canonical lowercase names (e.g., ARes.mag
            # → 'rho', Z.phz → 'phase', etc.).
            dfb = _standardise_columns(dfb)
    
            # Keep this block and record its per-block metadata.
            frames.append(dfb)
            blocks_meta.append(dict(block_meta))
    
            # # Optional: keep "sticky" Rx.* meta for subsequent blocks that
            # # omit it; otherwise clear fully
            # sticky = ('Rx.Stn', 'Rx.GdpStn', 'Rx.Cmp', 'Rx.Length')
            # block_meta = {k: v for k, v in block_meta.items() if k in sticky}
        
            # Reset block meta and continue scanning from block end.
            block_meta.clear()
            i = j
            continue 
        
        # 4) Anything else (blank lines, stray text)
        # Non-meta, non-table line → just advance
        i += 1
    
    if not frames:
        raise AvgDataError("Data block(s) missing in kind-2 file")

    # Concatenate all blocks into a single tidy frame.
    df = pd.concat(frames, ignore_index=True)

    # Derive a convenient boolean selection flag from CSAVGW
    # weights (1 = keep, 0 = skip).  If weights are absent, the
    # column is simply not added.
    if 'z.mwgt' in df.columns or 'z.pwgt' in df.columns:
        mw = _get_weight_bool (df, 'z.mwgt')
        pw = _get_weight_bool (df, 'z.pwgt')
        # mw = df.get('z.mwgt', 1).fillna(1).astype(float) > 0
        # pw = df.get('z.pwgt', 1).fillna(1).astype(float) > 0
        df['use'] = mw & pw

    # Merge top-level meta with collected per-block meta.
    meta: Dict[str, Any] = {**global_meta, 'blocks': blocks_meta}
    return df, meta

def _get_weight_bool (df, comp ='z.mwgt'): 
    """ Get the bool weight for construction use. """
    val =df.get(comp, 1) 
    if isinstance (val, (float, int)): 
        return val  > 0
    return df.get(comp, 1).fillna(1).astype(float) > 0


def split_by_station(
    df: pd.DataFrame
) -> Dict[Any, pd.DataFrame]:
    """
    Split a tidy AVG DataFrame into per-station sub-frames.

    The splitter is robust to dtype quirks by forcing 'station'
    to numeric and normalising NumPy scalars into plain Python
    types for dict keys.

    Parameters
    ----------
    df : pandas.DataFrame
        Tidy table that includes a 'station' column.

    Returns
    -------
    dict[Any, pandas.DataFrame]
        Mapping of station id → sub-DataFrame (index reset).

    Raises
    ------
    AvgDataError
        If 'station' column is not present.
    """
    if 'station' not in df.columns:
        raise AvgDataError(
            "'station' column missing – cannot split"
        )

    # Coerce 'station' to numeric if needed to avoid object
    # mixes and to keep group keys consistent.
    if not np.issubdtype(df['station'].dtype, np.number):
        df = df.copy()
        df['station'] = pd.to_numeric(
            df['station'], errors='coerce'
        )

    out: Dict[Any, pd.DataFrame] = {}

    # Use dropna=False so NaN stations (if any) are still grouped
    # and visible to the caller.
    for stn, sub in df.groupby(
        'station', sort=True, dropna=False
    ):
        # Normalise potential NumPy scalar to a plain Python
        # number for stable dict keys and friendly equality.
        try:
            key = np.asarray(stn).item()
        except Exception:
            key = float(stn) if pd.notna(stn) else stn

        # If the key is a clean integral float (e.g., 25.0),
        # store it as an int for ergonomic lookups.
        if (
            isinstance(key, float)
            and pd.notna(key)
            and key.is_integer()
        ):
            key = int(key)

        out[key] = sub.reset_index(drop=True)

    return out


@ensure_pkg ('xarray', extra="xarray is required for to_xarray()")
@isdf
def to_xarray(
    df: pd.DataFrame,
    *,
    coords: Sequence[str] = ("station", "freq", "comp"),
    data_vars: Optional[Sequence[str]] = None,
    attrs: Optional[Dict[str, Any]] = None,
) -> "xr.Dataset":
    """
    Convert a tidy Zonge table to an :class:`xarray.Dataset`.

    The resulting dataset uses a multi-dimensional grid with
    coordinates (typically ``station × freq × comp``).  Ragged
    sampling across stations is handled implicitly by NaNs in
    the corresponding data variables.

    Parameters
    ----------
    df :
        Long / tidy :class:`pandas.DataFrame` as returned by
        :func:`load_avg` (kind-1 or kind-2).  Expected columns
        include at least a subset of ``station, freq, comp``
        and one or more numeric data columns such as
        ``emag, hmag, rho, phase, …``.
    coords :
        Names of the DataFrame columns to use as coordinates
        and dataset dimensions.  Columns not present in *df*
        are ignored.  The default (``station, freq, comp``)
        matches common CSAMT usage.
    data_vars :
        Names of columns to export as data variables.  When
        *None*, all numeric columns **except** those listed in
        *coords* are used.
    attrs :
        Mapping of global attributes to attach to the dataset.
        A typical value is the *meta* dict returned by
        :func:`load_avg`.  Keys like ``"blocks"`` (per-block
        stash) are ignored to keep attributes clean.

    Returns
    -------
    xr.Dataset
        Dataset with dimensions given by the intersection of
        *coords* and the columns present in *df*.  Coordinate
        ordering is preserved (``station`` → ``freq`` → ``comp``
        by default).

    Notes
    -----
    * If ``comp`` is missing, a default value of ``"ExHy"`` is
      injected so the *comp* dimension exists.
    * Duplicate rows with identical coordinates are averaged
      (numeric columns) to ensure a single value per cell.
    * Boolean columns are preserved as data variables.
    """
    # Work on a copy; we will normalise types and sort below.
    df = df.copy()

    # Ensure a 'comp' column exists so we always get a comp
    # dimension (kind-1 files often omit component labels).
    if "comp" not in df.columns:
        df["comp"] = "ExHy"

    # Determine which coord columns we actually have.
    idx_cols = [c for c in coords if c in df.columns]
    if not idx_cols:
        raise AvgDataError(
            "No coordinate columns found. Expected any of: "
            f"{coords!r}"
        )

    # Light type normalisation: make station/freq numeric when
    # possible; keep comp as string/categorical.
    if "station" in idx_cols:
        df["station"] = pd.to_numeric(
            df["station"], errors="ignore"
        )
    if "freq" in idx_cols:
        df["freq"] = pd.to_numeric(df["freq"], errors="ignore")

    # Provide a stable, interpretable order for 'comp'.  Keep a
    # canonical order first, then append any unexpected labels.
    if "comp" in idx_cols:
        canon = [
            "ExHy", "ExHx", "EyHx", "EyHy",
            "Zxx", "Zxy", "Zyx", "Zyy", "Zvec", "Zdet",
        ]
        present = (
            pd.Series(df["comp"].astype(str).unique())
            .tolist()
        )
        extras = [c for c in present if c not in canon]
        cats = canon + extras
        df["comp"] = pd.Categorical(
            df["comp"].astype(str),
            categories=cats,
            ordered=True,
        )

    # Decide which columns become data variables.  Default to
    # all numeric (including bool) excluding coordinate cols.
    if data_vars is None:
        num_like = df.select_dtypes(
            include=[np.number, "bool", "boolean"]
        ).columns.tolist()
        data_vars = [c for c in num_like if c not in idx_cols]

    if not data_vars:
        raise AvgDataError(
            "No data variables selected. Provide 'data_vars' "
            "or ensure df has numeric columns."
        )

    # Reduce duplicates: some files may contain repeated rows
    # for the same (station, freq, comp).  We average numeric
    # variables across duplicates to ensure uniqueness.
    dup_mask = df.duplicated(subset=idx_cols, keep=False)
    if bool(dup_mask.any()):
        logger.warning(
            "Duplicate coordinate rows found; averaging numeric "
            "columns over duplicates."
        )
        # Only aggregate what we plan to emit as variables.
        gb = df.groupby(idx_cols, sort=True, dropna=False)
        df_num = gb[data_vars].mean()
        df = df_num.reset_index()
    else:
        # Sort for predictable coordinate order.
        df = df.sort_values(idx_cols, kind="mergesort")

    # Build the Dataset.  MultiIndex → dense grid; raggedness
    # becomes NaNs where combinations are missing.
    ds = df.set_index(idx_cols)[data_vars].to_xarray()

    # Order dimensions as requested, dropping missing ones.
    dim_order = [d for d in coords if d in ds.dims]
    ds = ds.transpose(*dim_order)

    # Attach user attributes, filtering internal per-block stash.
    if attrs:
        clean = dict(attrs)
        clean.pop("blocks", None)
        ds.attrs.update(clean)

    return ds

def write_avg(
    core:  pd.DataFrame,
    extra: pd.DataFrame | None,
    meta:  dict[str, str] | None,
    path:  str | Path | None = None,
    *,
    stamp: bool = True,
    float_fmt: str = "%.6g",
    na_rep: str = "",
    header_spaces: bool = False,   # use $k=v by default
) -> Path:
    """
    Serialize to Zonge kind-2 (CSAVGW/ASTATIC) .avg.

    If the input frame has a 'station' column, data are written
    as multiple blocks (one CSV section per station) with
    $Rx.* lines preceding each block (recommended). Otherwise a
    single block is written (legacy/utility case).
    """
    # --- 0) destination 
    if path is None:
        path = Path.cwd() / "exported_kind2.avg"
    path = Path(path).expanduser().resolve()

    # --- 1) build header (global $meta) 
    # filter out non-$ keys like 'blocks'
    meta = dict(meta or {})
    meta.pop("blocks", None)

    eq = " = " if header_spaces else "="
    header_lines: list[str] = []
    for k, v in meta.items():
        header_lines.append(f"${k}{eq}{v}")
    if stamp:
        utc = datetime.utcnow().isoformat(timespec="seconds") + "Z"
        header_lines.append(f"$Written{eq}{utc}")

    # one blank before first table
    header_lines.append("")
    out_chunks: list[str] = ["\n".join(header_lines)]

    # --- 2) assemble data block(s) 
    block = pd.concat(
        [core, extra], axis=1) if extra is not None else core
    block = block.copy()

    # canonical → kind-2 casing (patch a couple of gaps)
    canon2k2 = dict(_CANON2K2)
    canon2k2.update({
        "h.wgt": "H.wgt",
        "rho_sc": "SRes",
        "coh": "Choer",
        "gdp_blk": "Gdp.Blk",
        "gdp_chn": "Gdp.Chn",
        "gdp_time": "Gdp.Time",
        "zabs": "|Z|",
    })
    rename_map = {
        c: canon2k2.get(c.lower(), c) for c in block.columns}
    block.rename(columns=rename_map, inplace=True)

    # Expected CSAVGW order; extras will be appended after these.
    ordered = [
        "Z.mwgt","Z.pwgt","Freq","Tx.Amp",
        "E.mag","E.phz","B.mag","B.phz",
        "Z.mag","Z.phz","ARes.mag","SRes",
        "E.wgt","H.wgt",
        "E.%err","E.perr","B.%err","B.perr",
        "Z.%err","Z.perr","ARes.%err",
        "Choer","Gdp.Blk","Gdp.Chn","Gdp.Time","|Z|",
    ]

    def _order_cols(df: pd.DataFrame) -> list[str]:
        present = [c for c in ordered if c in df.columns]
        extras  = [c for c in df.columns if c not in present]
        # Avoid writing a 'Station' column for kind-2 blocks
        extras  = [c for c in extras if c.lower() != "station"]
        return present + extras

    # Multi-station writer (group and stamp $Rx.*)
    if "station" in (c.lower() for c in block.columns):
        # Identify the actual column name (case may vary)
        stn_col = next(c for c in block.columns if c.lower() == "station")

        for stn, sub in block.groupby(stn_col, sort=True, dropna=False):
            # Defensive: skip NaN station group
            if pd.isna(stn):
                continue

            # Per-block $Rx.* (best-effort; pull from columns if present)
            rx_lines: list[str] = []
            rx_stn = int(stn) if float(stn).is_integer() else stn
            rx_lines.append(f"$Rx.Stn{eq}{rx_stn}")

            # Optional goodies
            if "rx_gdpstn" in block.columns:
                rx_lines.append(f"$Rx.GdpStn{eq}{block['rx_gdpstn'].iloc[0]}")
            if "rx_length" in block.columns:
                rx_lines.append(f"$Rx.Length{eq}{block['rx_length'].iloc[0]}")
            # If comp is uniform within the group, stamp it
            comp_col = next(
                (c for c in block.columns if c.lower()=="comp"), None)
            if comp_col is not None:
                vals = sub[comp_col].dropna().unique()
                if len(vals) == 1:
                    rx_lines.append(f"$Rx.Cmp{eq}{vals[0]}")

            out_chunks.append("\n".join(rx_lines))

            # Write the CSV table
            dfw = sub.drop(columns=[stn_col], errors="ignore")
            cols = _order_cols(dfw)
            csv_txt = dfw[cols].to_csv(
                index=False, float_format=float_fmt, na_rep=na_rep
            )
            out_chunks.append(csv_txt.rstrip("\n"))
            out_chunks.append("")  # blank line after the block
    else:
        # Single-block writer
        cols = _order_cols(block)
        csv_txt = block[cols].to_csv(
            index=False, float_format=float_fmt, na_rep=na_rep
        )
        out_chunks.append(csv_txt.rstrip("\n"))

    # --- 3) write to disk 
    path.write_text("\n".join(out_chunks), encoding="utf8")
    logger.info("AVG written → %s", path)
    return path


def _standardise_columns(df: pd.DataFrame) -> pd.DataFrame:
    return df.rename(
        columns={c: _COL_MAP.get(c, c.lower()) for c in df.columns})

def load_avg(
    path: str | Path,
    *,
    ll_columns: Tuple[str, str] = ('latitude', 'longitude'),
    utm_zone: Optional[int] = None,
    inplace: bool = False
) -> Tuple[pd.DataFrame, Dict[str, str]]:
    """
    Read a Zonge AVG file (legacy kind-1 or modern CSAVGW
    kind-2) and return a tidy DataFrame plus metadata.

    The function:
      * classifies the file by header,
      * parses kind-1 (fixed-width) via _parse_kind1,
      * parses kind-2 (CSV blocks) via _parse_kind2,
      * optionally adds UTM eastings/northings if lat/lon cols
        are present (column names configurable via ll_columns).

    Parameters
    ----------
    path : str | Path
        Filesystem path to the .avg file.
    ll_columns : (str, str), optional
        Column names containing latitude/longitude (deg).
    utm_zone : int, optional
        UTM zone override; if None, auto-detect.
    inplace : bool, default False
        If False, return a copy so callers can mutate safely.

    Returns
    -------
    (df, meta) : (pandas.DataFrame, dict)
        df    : tidy table with normalised column names.
        meta  : metadata dict.  For kind-2, includes a 'blocks'
                list with per-block metadata.

    Raises
    ------
    FileNotFoundError
        If the file path does not exist.
    AvgFileError / AvgDataError
        If classification or parsing fails.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(path)

    # Read text and classify by header patterns.
    lines = path.read_text(errors='replace').splitlines()
    kind = classify_avg_format(lines)

    # Dispatch to the appropriate parser.
    if kind == 1:
        df, meta = _parse_kind1(lines), {}
    else:
        df, meta = _parse_kind2(lines)

    # Copy unless caller explicitly wants "inplace".
    if not inplace:
        df = df.copy()

    # Optionally compute UTM if lat/lon columns are present.
    lat_col, lon_col = ll_columns
    if lat_col in df.columns and lon_col in df.columns:
        try:
            east, north = ll_to_utm(
                df[lat_col].values,
                df[lon_col].values,
                zone=utm_zone
            )
            df['easting'] = east
            df['northing'] = north
        except Exception as exc:  # pragma: no cover
            logger.warning(
                "Lat/Lon → UTM failed: %s",
                exc
            )

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

def _find_col(
    df: pd.DataFrame, candidates: Sequence[str]
) -> Optional[str]:
    """
    Return the first column name present in *df* among
    *candidates*.  Matching is case-insensitive and ignores
    whitespace.
    """
    low = {str(c).strip().lower(): c for c in df.columns}
    for want in candidates:
        key = str(want).strip().lower()
        if key in low:
            return low[key]
    return None


def _to_num(x: object) -> float:
    """
    Robust numeric coercion.  Strings like '', '*', 'NaN' become
    ``np.nan``.  Otherwise return a float when possible.
    """
    if x is None:
        return np.nan
    s = str(x).strip()
    if s in {"", "*", "nan", "NaN", "None", "null"}:
        return np.nan
    try:
        return float(s)
    except Exception:
        return np.nan


def _norm_comp(x: object) -> str:
    """
    Canonicalize component labels into the 2×2 slots used in the
    tensor layout.  Only a few common forms are normalized here;
    everything else is passed through as a string.

    Examples
    --------
    'exhy' → 'ExHy', 'EYHX' → 'EyHx', 'zxy' → 'Zxy'
    """
    if x is None:
        return "ExHy"
    s = str(x).strip()
    if not s:
        return "ExHy"

    s_up = s.upper()
    # zxx/zxy/zyx/zyy are accepted as-is (capitalized later)
    # classic galvanic pairs:
    if s_up in {"EXHY", "EX-HY", "E X H Y"}:
        return "ExHy"
    if s_up in {"EXHX", "EX-HX"}:
        return "ExHx"
    if s_up in {"EYHX", "EY-HX"}:
        return "EyHx"
    if s_up in {"EYHY", "EY-HY"}:
        return "EyHy"

    # impedance-like notations:
    if s_up in {"ZXX", "ZXY", "ZYX", "ZYY"}:
        return s_up.capitalize()

    # fallback: capitalize first, keep inner case
    return s[0:1].upper() + s[1:]
