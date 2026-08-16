#       Author: LKouadio <etanoyau@gmail.com>
#       License: LGPL-3.0

"""
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

import io
import re
import warnings  # noqa
from collections.abc import Iterable, Sequence
from datetime import datetime
from numbers import Integral
from pathlib import Path
from typing import (
    Any,
)

import numpy as np
import pandas as pd

try:
    import xarray as xr  # type: ignore
except ImportError:  # pragma: no cover
    pass
    # warnings.warn(
    #     "xarray is required for the package"
    # )
from ..compat.aliases import compat_alias
from ..decorators import isdf
from ..exceptions import (
    AvgDataError,
    AvgFileError,
    StationError,
)
from ..gis.utils import to_utm  # type: ignore
from ..log.logger import get_logger
from ..utils._dependency import import_optional_dependency
from .schema import (
    _CANON_TO_MODERN,
    _CANONICAL_MAP,
    _CSAVGW_ORDERED,
    _FLEXIBLE_LOOKUP,
    get_aliases,
)

__all__ = [
    "load_avg",
    "round_dipole_length",
    "validate_stn_profile",
    "classify_avg_format",
    "extract_core_columns",
    "number_stations",
    "chunk_by_frequency",
    "write_avg",
    "to_xarray",
    "to_numeric_if_possible",
]

logger = get_logger(__name__)


_RX_WS = re.compile(r"\s+")

_RX_K2_HEADER = re.compile(r"^\s*Z\.mwgt\s*,", re.I)
_RX_K1_HEADER = re.compile(r"^\s*skp\s+Station", re.I)

_NUMERIC_REPLACE = {"*": np.nan, "nan": np.nan, "NaN": np.nan, "": np.nan}
_COMMENT_PREFIXES = ("\\", "/", "!", '"')


def find_and_rename_column(
    df: pd.DataFrame, canonical_name: str
) -> pd.DataFrame:
    """
    Find a column by any of its aliases and rename it to the
    canonical name.
    """
    # Use the get_aliases function we already built
    aliases = get_aliases(canonical_name, kind="all")

    col_found = None
    for alias in aliases:
        if alias in df.columns:
            col_found = alias
            break

    if col_found and col_found != canonical_name:
        return df.rename(columns={col_found: canonical_name})

    return df


def _to_float(val: str | float | int) -> float | np.floating:
    """Convert *val* to float while honouring project placeholders."""
    if isinstance(val, (float, int)):
        return float(val)
    txt = val.strip()
    if txt in _NUMERIC_REPLACE:
        return np.nan
    if txt.startswith("."):
        txt = "0" + txt
    if txt.endswith("."):
        txt = txt + "0"
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
            logger.debug("Kind-2 AVG detected (found modern header).")
            return 2
        if _RX_K1_HEADER.search(ln):
            logger.debug("Kind-1 AVG detected (found legacy header).")
            return 1

    # Second pass: If no header found, analyze file structure.
    has_commas = False
    has_dot_in_keyword = False

    for ln in lines:
        s = ln.strip()
        if not s:
            continue

        # Check for modern keyword style, e.g., '$Survey.Type='
        if s.startswith("$"):
            # This is a very strong indicator of a modern file.
            match = re.match(r"\$\s*(\w+\.\w+)\s*=", s)
            if match:
                has_dot_in_keyword = True
                break

        # Heuristic: check for comma-separated data.
        # Avoids matching headers with only one or two commas.
        if not s.startswith(("$", "\\", "/")) and s.count(",") > 2:
            has_commas = True

    if has_dot_in_keyword:
        logger.debug("Kind-2 AVG detected (found dot-notation keywords).")
        return 2

    if has_commas:
        logger.debug("Kind-2 AVG detected (found comma-separated data).")
        return 2

    # If no strong indicators are found, parsing cannot proceed.
    raise AvgFileError("Unrecognised AVG header – cannot classify file")


def _parse_kind1(lines: Sequence[str]) -> pd.DataFrame:
    """Parse legacy fixed‑width (kind‑1) AVG table."""
    idx = next(
        (i for i, ln in enumerate(lines) if _RX_K1_HEADER.search(ln)), None
    )
    if idx is None:
        raise AvgFileError("Header row not found in kind‑1 file")

    hdr_tokens = _RX_WS.sub(" ", lines[idx].strip()).split()
    data_rows: list[list[Any]] = []
    for ln in lines[idx + 1 :]:
        if not ln.strip() or _RX_K1_HEADER.search(ln):
            break
        if ln.startswith(("\\", "$")):
            continue
        tokens = _RX_WS.sub(" ", ln.strip()).split()
        data_rows.append(
            [
                _to_float(tk) if j >= 4 else tk  # First 4 cols are non-numeric
                for j, tk in enumerate(tokens)
            ]
        )
    if not data_rows:
        raise AvgDataError("No data rows in kind‑1 file")

    # Some GDP/AMTAVG builds (e.g. AMTAVG 7.40) emit one or more extra
    # trailing data columns that are not listed in the "skp Station ..."
    # header row. When every row consistently carries the same surplus
    # of fields, extend the header with generic names instead of letting
    # pandas raise on the column-count mismatch.
    row_lengths = {len(row) for row in data_rows}
    if len(row_lengths) == 1:
        n_extra = row_lengths.pop() - len(hdr_tokens)
        if n_extra > 0:
            logger.debug(
                f"Kind-1 AVG data rows carry {n_extra} more field(s) than "
                "the header row; appending generic column name(s) for "
                "the surplus."
            )
            hdr_tokens = hdr_tokens + [
                f"extra_{i + 1}" for i in range(n_extra)
            ]

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
        if _RX_K2_HEADER.match(s) or s.startswith("$"):
            break
        i += 1
    return (start, i)


def _parse_kind2(lines: Sequence[str]) -> tuple[pd.DataFrame, dict[str, str]]:
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
    global_meta: dict[str, str] = {}
    blocks_meta: list[dict[str, str]] = []
    frames: list[pd.DataFrame] = []

    block_meta: dict[str, str] = {}
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
        if ln.startswith("$") and "=" in ln:
            key, val = ln[1:].split("=", 1)
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
            table_txt = "\n".join(
                [lines[start]]
                + [s for s in lines[start + 1 : j] if not _is_comment(s)]
            )

            # Parse CSV with forgiving whitespace.  Convert numeric
            # strings (including '*', '.5', '1.') with _to_float.
            dfb = pd.read_csv(io.StringIO(table_txt), skipinitialspace=True)
            # applymap was renamed to map in pandas ≥2.1
            _dfb_map = getattr(dfb, "map", None) or dfb.applymap
            dfb = _dfb_map(lambda v: _to_float(v) if isinstance(v, str) else v)

            # Stamp station and a few helpful block-level fields as
            # columns.  Prefer client station number ($Rx.Stn).
            stn = (
                block_meta.get("Rx.Stn")
                or block_meta.get("Rx.GdpStn")
                or block_meta.get("Stn.Beg")
            )
            if stn is not None:
                try:
                    dfb["station"] = _to_float(stn)
                except Exception:
                    dfb["station"] = stn  # keep as text if odd

            # Component label and a couple of helpers can be handy for
            # QC.  They are optional and harmless if missing.
            if "Rx.Cmp" in block_meta:
                dfb["comp"] = block_meta["Rx.Cmp"]
            for k in ("Rx.Length", "Rx.GdpStn"):
                if k in block_meta:
                    dfb[k.replace(".", "_").lower()] = block_meta[k]

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
    if "z_mwgt" in df.columns or "z_pwgt" in df.columns:
        mw = _get_weight_bool(df, "z_mwgt")
        pw = _get_weight_bool(df, "z_pwgt")
        # mw = df.get('z.mwgt', 1).fillna(1).astype(float) > 0
        # pw = df.get('z.pwgt', 1).fillna(1).astype(float) > 0
        df["use"] = mw & pw

    # Merge top-level meta with collected per-block meta.
    meta: dict[str, Any] = {**global_meta, "blocks": blocks_meta}
    return df, meta


def _get_weight_bool(df, comp="z.mwgt"):
    """Get the bool weight for construction use."""
    val = df.get(comp, 1)
    if isinstance(val, (float, int)):
        return val > 0
    return df.get(comp, 1).fillna(1).astype(float) > 0


def _standardise_columns(df: pd.DataFrame) -> pd.DataFrame:
    """
    Rename columns to a canonical schema using a two-pass strategy.
    """
    # Create a simple case-insensitive map for strict matching
    strict_lower_map = {k.lower(): v for k, v in _CANONICAL_MAP.items()}
    rename_dict = {}
    for col in df.columns:
        # Normalize the column name for flexible lookup
        norm_col = str(col).lower().strip()  # .replace(
        # '.', '').replace('_', '').replace('%', '')

        # 1. Try flexible lookup first for QC/weight columns
        if norm_col in _FLEXIBLE_LOOKUP:
            rename_dict[col] = _FLEXIBLE_LOOKUP[norm_col]
        # 2. Fallback to strict, case-insensitive lookup for all others
        elif str(col).lower() in strict_lower_map:
            rename_dict[col] = strict_lower_map[str(col).lower()]
        # 3. If no match, keep original name (or lowercase it)
        else:
            rename_dict[col] = str(col).lower()

    return df.rename(columns=rename_dict)


def to_numeric_if_possible(values):
    """Return numeric values when coercion succeeds, else original values.

    This avoids ``pd.to_numeric(errors="ignore")`` because Dask-backed
    pandas objects reject the ``"ignore"`` mode under newer CI stacks.
    """
    try:
        numeric = pd.to_numeric(values, errors="coerce")
    except Exception:
        return values

    try:
        original_na = pd.isna(values)
        converted_na = pd.isna(numeric)
        if bool((converted_na & ~original_na).any()):
            return values
    except Exception:
        return values
    return numeric


def split_by_station(df: pd.DataFrame) -> dict[Any, pd.DataFrame]:
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
    if "station" not in df.columns:
        raise AvgDataError("'station' column missing – cannot split")

    # Coerce 'station' to numeric if needed to avoid object
    # mixes and to keep group keys consistent.
    if not np.issubdtype(df["station"].dtype, np.number):
        df = df.copy()
        df["station"] = pd.to_numeric(df["station"], errors="coerce")

    out: dict[Any, pd.DataFrame] = {}

    # Use dropna=False so NaN stations (if any) are still grouped
    # and visible to the caller.
    for stn, sub in df.groupby("station", sort=True, dropna=False):
        # Normalise potential NumPy scalar to a plain Python
        # number for stable dict keys and friendly equality.
        try:
            key = np.asarray(stn).item()
        except Exception:
            key = float(stn) if pd.notna(stn) else stn

        # If the key is a clean integral float (e.g., 25.0),
        # store it as an int for ergonomic lookups.
        if isinstance(key, float) and pd.notna(key) and key.is_integer():
            key = int(key)

        out[key] = sub.reset_index(drop=True)

    return out


@isdf
def to_xarray(
    df: pd.DataFrame,
    *,
    coords: Sequence[str] = ("station", "freq", "comp"),
    data_vars: Sequence[str] | None = None,
    attrs: dict[str, Any] | None = None,
) -> xr.Dataset:
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
    import_optional_dependency(
        "xarray",
        extra="xarray is required for to_xarray()",
        errors="raise",
    )
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
            f"No coordinate columns found. Expected any of: {coords!r}"
        )

    # Light type normalisation: make station/freq numeric when
    # possible; keep comp as string/categorical.
    if "station" in idx_cols:
        df["station"] = to_numeric_if_possible(df["station"])
    if "freq" in idx_cols:
        df["freq"] = to_numeric_if_possible(df["freq"])

    # Provide a stable, interpretable order for 'comp'.  Keep a
    # canonical order first, then append any unexpected labels.
    if "comp" in idx_cols:
        canon = [
            "ExHy",
            "ExHx",
            "EyHx",
            "EyHy",
            "Zxx",
            "Zxy",
            "Zyx",
            "Zyy",
            "Zvec",
            "Zdet",
        ]
        present = pd.Series(df["comp"].astype(str).unique()).tolist()
        # Only categories actually present in the data: including
        # unused canonical labels here would make the later
        # set_index(...).to_xarray() densify a full cartesian grid
        # over all 10 canonical components, not just the ones the
        # file actually reports.
        cats = [c for c in canon if c in present] + [
            c for c in present if c not in canon
        ]
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
    core: pd.DataFrame,
    extra: pd.DataFrame | None,
    meta: dict[str, str] | None,
    path: str | Path | None = None,
    *,
    stamp: bool = True,
    float_fmt: str = "%.6g",
    na_rep: str = "*",
    header_spaces: bool = False,  # use $k=v by default
    banner_lines: Sequence[str] | None = None,
) -> Path:
    r"""Serialize a DataFrame to a Zonge kind-2 AVG file.

    This function serves as the core writer for creating modern,
    CSAVGW/ASTATIC-style ``.avg`` files. It takes a DataFrame
    with canonical column names, aggregates metadata, and formats
    the output into a structured text file with data blocks
    grouped by station.

    Parameters
    ----------
    core : pandas.DataFrame
        The main DataFrame containing the core measurement data.
        It is expected to have canonical column names (e.g.,
        'rho', 'phase', 'pc_emag').
    extra : pandas.DataFrame or None
        An optional DataFrame containing additional columns to be
        merged with the core data before writing.
    meta : dict, optional
        A dictionary of global metadata to be written as
        ``$keyword=value`` pairs in the file header.
    path : str or pathlib.Path, optional
        The output file path. If ``None``, a default filename
        like ``exported_kind2.avg`` is created in the current
        working directory.
    stamp : bool, default True
        If ``True``, a ``$Written=<timestamp>`` line is added to
        the header for provenance.
    float_fmt : str, default "%.6g"
        The format specifier for writing floating-point numbers
        in the data blocks.
    na_rep : str, default "*"
        The string representation for missing (NaN) values in the
        data blocks, conforming to the Zonge convention.
    header_spaces : bool, default False
        If ``True``, adds spaces around the equals sign in header
        keywords (e.g., ``$Key = Value``).
    banner_lines : sequence of str, optional
        A sequence of comment lines (starting with '\\') to be
        prepended to the file header, typically for hardware and
        processing information.

    Returns
    -------
    pathlib.Path
        The absolute path to the newly created ``.avg`` file.

    Notes
    -----
    This function is designed to produce clean, compliant, and
    human-readable kind-2 AVG files. Its key behaviors include:

    - **Column Renaming**: It uses the ``CANON_TO_MODERN_MAP`` from
      the :mod:`~.schema` module to automatically convert the
      internal canonical column names back to the standard modern
      format (e.g., 'rho' becomes 'ARes.mag').
    - **Block Grouping**: If a 'station' column is present, the
      function intelligently groups the data by station and writes
      a separate data block for each, preceded by its specific
      ``$Rx.*`` metadata.
    - **Smart Column Filtering**: It automatically detects and
      omits placeholder columns (like 'Choer', 'Gdp.Blk') if they
      contain no valid data, preventing empty columns from
      cluttering the output file. It also excludes internal helper
      columns from the final output.
    - **Aligned Formatting**: The function uses a custom CSV
      formatter to produce neatly aligned, fixed-width-like
      columns within the data blocks, matching the appearance of
      files generated by Zonge's proprietary software.

    See Also
    --------
    pycsamt.zonge.avg.BaseAVG.to_modern : The primary method that
        calls this function.
    pycsamt.zonge.utils.load_avg : The corresponding function for
        reading AVG files.
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
    header_lines: list[str] = list(banner_lines or [])
    if banner_lines:
        header_lines.append("")

    for k, v in meta.items():
        header_lines.append(f"${k}{eq}{v}")
    if stamp:
        utc = datetime.utcnow().isoformat(timespec="seconds") + "Z"
        header_lines.append(f"$Written{eq}{utc}")

    # one blank before first table
    header_lines.append("")
    out_chunks: list[str] = ["\n".join(header_lines)]

    # --- 2) assemble data block(s)
    block = pd.concat([core, extra], axis=1) if extra is not None else core
    block = block.copy()

    # Drop any "extra" columns that are completely empty
    extra_cols_to_check = ["coh", "gdp_blk", "gdp_chn", "gdp_time", "zabs"]
    for col in extra_cols_to_check:
        if col in block.columns and block[col].isnull().all():
            block = block.drop(columns=[col])

    # canonical → kind-2 casing (patch a couple of gaps)
    canon_to_modern = {
        k: v for k, v in _CANON_TO_MODERN.items() if k in block.columns
    }
    block.rename(columns=canon_to_modern, inplace=True)

    def _order_cols(df: pd.DataFrame) -> list[str]:
        # Expected CSAVGW order; extras will be appended after these.

        present = [c for c in _CSAVGW_ORDERED if c in df.columns]
        extras = [c for c in df.columns if c not in present]
        # Exclude all columns that are part of the block's metadata
        extras = [
            c
            for c in extras
            if c.lower() not in ("station", "comp", "rx_length", "rx_gdpstn")
        ]
        return present + extras

    # Identify the actual column name (case may vary)
    stn_col = next((c for c in block.columns if c.lower() == "station"), None)
    # Multi-station writer (group and stamp $Rx.*)
    # if "station" in (c.lower() for c in block.columns):
    #     # Identify the actual column name (case may vary)
    #     stn_col = next(c for c in block.columns if c.lower() == "station")
    if stn_col:
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
                (c for c in block.columns if c.lower() == "comp"), None
            )
            if comp_col is not None:
                vals = sub[comp_col].dropna().unique()
                if len(vals) == 1:
                    rx_lines.append(f"$Rx.Cmp{eq}{vals[0]}")

            out_chunks.append("\n".join(rx_lines))

            cols_to_write = _order_cols(sub)
            dfw = sub[cols_to_write]
            out_chunks.append(_format_aligned_csv(dfw, float_fmt, na_rep))

            out_chunks.append("")  # blank line after the block
    else:
        # Single-block writer
        cols = _order_cols(block)
        dfw = block[cols]
        out_chunks.append(_format_aligned_csv(dfw, float_fmt, na_rep))
    # --- 3) write to disk
    path.write_text("\n".join(out_chunks), encoding="utf8")
    logger.info("AVG written → %s", path)
    return path


def _format_aligned_csv(df: pd.DataFrame, float_fmt: str, na_rep: str) -> str:
    """Formats a DataFrame to a perfectly aligned CSV string."""
    # Convert all data to string representation first
    df_str = pd.DataFrame(index=df.index)
    for col in df.columns:
        if pd.api.types.is_float_dtype(df[col]):
            df_str[col] = df[col].apply(
                lambda x: na_rep if pd.isna(x) else float_fmt % x
            )
        else:
            df_str[col] = df[col].apply(
                lambda x: na_rep if pd.isna(x) else str(x)
            )

    # Calculate max width for each column
    widths = {
        col: max(df_str[col].str.len().max(), len(col))
        for col in df_str.columns
    }

    # Format header and rows
    header = ",".join(f"{col:<{widths[col]}}" for col in df.columns)
    rows = [header]
    for _, row in df_str.iterrows():
        rows.append(
            ",".join(f"{val:<{widths[col]}}" for col, val in row.items())
        )
    return "\n".join(rows)


def load_avg(
    path: str | Path,
    *,
    ll_columns: tuple[str, str] = ("latitude", "longitude"),
    utm_zone: int | None = None,
    inplace: bool = False,
) -> tuple[pd.DataFrame, dict[str, str]]:
    r"""Read a Zonge AVG file and return a tidy DataFrame and metadata.

    This function serves as the primary parser for both legacy
    (kind-1) and modern (kind-2) Zonge AVG files. It automatically
    detects the file format, parses the data accordingly, and
    standardizes all column names to a consistent, canonical
    schema.

    Parameters
    ----------
    path : str or pathlib.Path
        The filesystem path to the ``.avg`` file.
    ll_columns : tuple[str, str], default ('latitude', 'longitude')
        Column names in the source data that contain latitude and
        longitude values in decimal degrees. This is not a standard
        Zonge field but is supported for custom data formats.
    utm_zone : int, optional
        The UTM zone number to use for coordinate conversion. If
        ``None``, the zone is auto-detected from the longitude.
    inplace : bool, default False
        This parameter is deprecated and no longer has an effect,
        as the function always returns a new DataFrame.

    Returns
    -------
    df : pandas.DataFrame
        A tidy DataFrame where all column names have been
        standardized to the internal canonical schema (e.g.,
        'Resistivity' becomes 'rho').
    meta : dict
        A dictionary containing all header metadata from the file.
        For modern files, this includes a 'blocks' key with a
        list of per-station metadata blocks.

    Raises
    ------
    FileNotFoundError
        If the file specified by `path` does not exist.
    AvgFileError
        If the file format cannot be reliably classified as either
        legacy or modern.
    AvgDataError
        If parsing fails due to malformed or missing data within
        the file.

    Notes
    -----
    The standardization of column names is a key feature of this
    function. It ensures that all subsequent processing steps can
    rely on a consistent and predictable data structure, regardless
    of the input file's original format. This is achieved by using
    the ``_CANONICAL_MAP`` from the :mod:`~.schema` module.

    Examples
    --------
    >>> from pycsamt.zonge.utils import load_avg
    >>> # Load a modern AVG file
    >>> df, meta = load_avg("data/avg/K2.avg")
    >>> print(df.columns)
    Index(['z_mwgt', 'freq', ..., 'station', 'comp', 'use'], dtype='object')
    >>> print(meta["Survey.Type"])
    CSAMT

    See Also
    --------
    pycsamt.zonge.avg.AVG.from_file : The recommended high-level
        entry point for loading AVG data.
    pycsamt.zonge.utils.write_avg : The corresponding function for
        writing AVG files.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(path)

    # Read text and classify by header patterns.
    lines = path.read_text(errors="replace").splitlines()
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
            east, north, _ = to_utm(
                df[lat_col].values, df[lon_col].values, zone=utm_zone
            )
            df["easting"] = east
            df["northing"] = north
        except Exception as exc:  # pragma: no cover
            logger.warning("Lat/Lon → UTM failed: %s", exc)

    return df, meta


def read_stn(path: str | Path) -> pd.DataFrame:
    r"""
    Parse Zonge ``.stn`` files (legacy and extended forms).

    Supports:
    - legacy, space-delimited with four columns
      (``station easting northing elevation``)
    - CSV with quoted headers
    - extended CSV headers containing optional columns
      (e.g., heading, pitch, roll)
    - embedded first-data row at end of header line
      (e.g., ``...,roll -200,1472...``)

    Lines starting with ``!``, ``/``, ``\``, ``#``, or ``;``
    are treated as comments and skipped.

    Returns a raw DataFrame; column normalization is left to
    the caller (e.g., mapping ``dot/e/n/h`` to canonical
    names).
    """
    try:
        with open(Path(path), encoding="utf-8") as f:
            raw = f.read().splitlines()
    except Exception as exc:
        raise StationError(f"Cannot read STN file: {exc}") from exc

    # keep non-empty, non-comment lines; do not drop lines
    # beginning with quotes because quoted headers are valid
    def _is_comment(s: str) -> bool:
        s = s.lstrip()
        return (not s) or s[0] in {"!", "/", "\\", "#", ";"}

    lines: list[str] = [ln.strip() for ln in raw if not _is_comment(ln)]

    if not lines:
        raise StationError("Empty or comment-only STN file.")

    # detect header line: choose the last line that contains
    # any letter; numeric-only lines are data
    header_idx = -1
    for i, ln in enumerate(lines):
        if re.search(r"[A-Za-z]", ln):
            header_idx = i

    # split header lines that also contain data, e.g.:
    # "...,roll -200,1472..." → separate header/data
    if header_idx >= 0:
        m = re.match(
            r"^(?P<head>.*[A-Za-z].*?)\s+(?P<data>[-+0-9].*)$",
            lines[header_idx],
        )
        if m and ("," in m.group("data") or re.search(r"\s", m.group("data"))):
            lines[header_idx] = m.group("head")
            lines.insert(header_idx + 1, m.group("data"))

    # choose delimiter: comma preferred if present in header
    # or first data line; otherwise whitespace
    def _has_comma(s: str) -> bool:
        return "," in s

    delim: str
    if header_idx >= 0 and _has_comma(lines[header_idx]):
        delim = ","
    else:
        # look at first data line
        data_start = header_idx + 1 if header_idx >= 0 else 0
        if data_start < len(lines) and _has_comma(lines[data_start]):
            delim = ","
        else:
            delim = r"\s+"

    # build header
    if header_idx >= 0:
        raw_header = lines[header_idx]
        # strip stray HTML tags (e.g. a "<input type="hidden" />"
        # fragment accidentally pasted in from a browser-based export)
        # then common quotes (", “, ”) and extra whitespace
        cleaned = (
            re.sub(r"<[^>]*>", "", raw_header)
            .replace('"', "")
            .replace("“", "")
            .replace("”", "")
            .strip()
        )
        header = [h.strip() for h in re.split(delim, cleaned)]
        data_lines = lines[header_idx + 1 :]
    else:
        # no header: assume 4 legacy columns
        header = ["station", "easting", "northing", "elevation"]
        data_lines = lines

    if not data_lines:
        raise StationError("No data rows found in STN file.")

    # join data and parse with pandas; use python engine for
    # regex separators and ragged rows
    data_str = "\n".join(data_lines)
    try:
        df = pd.read_csv(
            io.StringIO(data_str),
            sep=delim,
            header=None,
            names=header,
            engine="python",
        )
    except Exception as exc:
        raise StationError(f"Failed to parse STN data: {exc}") from exc

    # drop fully-empty rows
    df = df.dropna(how="all")

    if df.empty:
        raise StationError("STN parse produced empty DataFrame.")

    return df


@compat_alias(
    "validate_stn_profile",
    since="2.0.0",
    remove_in="2.17.0",
    export=True,
    extra=("Use 'tensor2d'. Removal in v2.17.0."),
)
def detect_stn_header(
    profile: Sequence[str],
    splitter: str | None = None,
) -> tuple[int, list[tuple[str, int]]]:
    r"""
    Heuristically detect a ``.stn`` header and map token
    positions.

    Scans non-comment lines, tolerates CSV and whitespace,
    quoted labels (e.g., ``3x"dot3x"``), and embedded
    ``label=value`` segments. Returns a score (matched
    labels) and a list of ``(canonical_label, index)`` for
    the best header line found. If no header is detected,
    returns ``(0, [])``.

    Parameters
    ----------
    profile : sequence of str
        Raw lines of the STN file.
    splitter : str or None, default None
        Token delimiter. If ``None``, auto-detect per line
        (comma if present, else whitespace).

    Returns
    -------
    score : int
        Number of recognized header labels on the best line.
    matches : list of (str, int)
        Pairs of canonical label name and column index.
    """

    # comments and empties
    def _is_comment(s: str) -> bool:
        s = s.lstrip()
        return (not s) or s[0] in {"!", "/", "\\", "#", ";"}

    # alias map → canonical label
    aliases = {
        # station
        "dot": "station",
        "station": "station",
        "sta": "station",
        # easting
        "e": "easting",
        "east": "easting",
        "easting": "easting",
        # northing
        "n": "northing",
        "north": "northing",
        "northing": "northing",
        # elevation
        "h": "elevation",
        "elev": "elevation",
        "elevation": "elevation",
        # optional extras
        "line": "line",
        "lon": "lon",
        "lat": "lat",
        "utm_zone": "utm_zone",
        "heading": "heading",
        "pitch": "pitch",
        "roll": "roll",
    }

    def _split_line(ln: str) -> list[str]:
        # remove common quote chars
        s = ln.replace('"', "").replace("“", "").replace("”", "").strip()
        # strip label=value to label
        s = re.sub(r"\s*=\s*[^, \t]+", "", s)
        sep = "," if (splitter is None and "," in s) else (splitter or r"\s+")
        toks = re.split(sep, s)
        out: list[str] = []
        for t in toks:
            t = t.strip().lower()
            # keep letters and underscores only
            t = re.sub(r"[^a-z_]", "", t)
            if t:
                out.append(t)
        return out

    best_score = 0
    best_matches: list[tuple[str, int]] = []

    # evaluate each non-comment line; pick the best header
    for ln in profile:
        if _is_comment(ln):
            continue

        tokens = _split_line(ln)
        if not tokens:
            continue

        # header+data on same line → keep left of first digit
        m = re.search(r"[-+]?\d", ln)
        if m:
            # limit tokenization to the header part
            head = ln[: m.start()]
            tokens = _split_line(head) or tokens

        matches: list[tuple[str, int]] = []
        for idx, tk in enumerate(tokens):
            canon = aliases.get(tk)
            if canon is not None:
                matches.append((canon, idx))

        score = len(matches)
        if score > best_score or (score == best_score and score > 0):
            best_score = score
            best_matches = matches

    return best_score, best_matches


# Back-compat alias
validate_stn_profile = detect_stn_header


def round_dipole_length(length: float | int) -> float:
    """Round *length* to the nearest 5‑m increment."""
    length = float(length)
    mod = length % 5
    if mod < 3:
        return length - mod
    if mod < 7:
        return length - mod + 5
    return length - mod + 10


def extract_core_columns(
    df: pd.DataFrame, *, keep: Iterable[str] | None = None
) -> tuple[pd.DataFrame, pd.DataFrame]:
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
        "station",
        "freq",
        "amps",
        "emag",
        "ephz",
        "hmag",
        "hphz",
        "rho",
        "phase",
        "e.%err",
        "e.perr",
        "h.%err",
        "h.perr",
        "rho.%err",
        "phase.%err",
    }
    keep_set = set(k.lower() for k in (keep or default_keep))
    cols_lower = {c.lower(): c for c in df.columns}
    have_keep = [cols_lower[c] for c in keep_set if c in cols_lower]
    if "station" not in [c.lower() for c in have_keep]:
        have_keep.insert(0, cols_lower.get("station", df.columns[0]))
    core = df[have_keep].copy()
    extra = df.drop(columns=have_keep).copy()
    return core.reset_index(drop=True), extra.reset_index(drop=True)


def _block_to_dict(block: Sequence[str]) -> dict[str, Any]:
    """Convert ``key=value`` lines to a dict (case‑insensitive keys)."""
    out: dict[str, Any] = {}
    for ln in block:
        if "=" not in ln:
            continue
        k, v = ln.split("=", 1)
        out[k.strip().lower()] = v.strip()
    return out


def _dict_to_lines(data: Any) -> list[str]:
    """Serialise dict‑like *data* to ``key=value`` text lines."""
    if isinstance(data, str):
        import json

        data = json.loads(data)
    if not isinstance(data, dict):
        data = dict(data)
    return [f"{k}={v}\n" for k, v in data.items() if v is not None]


def number_stations(
    n_stations: int | Integral, n_freq: int | Integral, *, prefix: str = "S"
) -> tuple[list[str], list[str]]:
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
    ids = [f"{prefix}{i:02d}" for i in range(int(n_stations))]
    expanded = list(np.repeat(ids, int(n_freq)))
    return ids, expanded


def chunk_by_frequency(
    data: Sequence[Any] | np.ndarray,
    n_per_chunk: int | Integral,
    *,
    drop_remainder: bool = False,
) -> list[np.ndarray]:
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
    arr = np.asarray(data)
    total = arr.size
    idx = np.arange(0, total, int(n_per_chunk))
    chunks: list[np.ndarray] = [arr[i : i + n_per_chunk] for i in idx]
    if drop_remainder and chunks and chunks[-1].size < n_per_chunk:
        chunks.pop()
    return chunks


def _find_col(df: pd.DataFrame, candidates: Sequence[str]) -> str | None:
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


def _first_present(
    df: pd.DataFrame,
    candidates: Sequence[str],
) -> str | None:
    """
    Return the first column name found in *df*
    among *candidates*.

    Parameters
    ----------
    df
        Input tidy table.
    candidates
        Ordered list of legacy/modern aliases.

    Returns
    -------
    str or None
        The first matching column name, or *None* if absent.
    """

    cols_lc = {str(c).strip().lower(): c for c in df.columns}
    for a in candidates:
        a_lc = str(a).strip().lower()
        if a_lc in cols_lc:
            return cols_lc[a_lc]
    return None


def _to_numeric_percent(series: pd.Series) -> pd.Series:
    """
    Convert a percent-like column to float while tolerating blanks.

    Empty/asterisk tokens are mapped to NaN; strings like ``"5"`` or
    ``"5.0"`` become ``5.0``.
    """
    s = series.astype(str).str.strip().replace({"": np.nan, "*": np.nan})
    return pd.to_numeric(s, errors="coerce")


def _to_complex(x: Any) -> complex | float:
    """
    Robustly convert a value to a complex number.

    Handles strings, complex, float, and int types, returning
    np.nan for values that cannot be converted.
    """
    if x is None:
        return np.nan
    s = str(x).strip()
    if s in {"", "*", "nan", "NaN", "None", "null"}:
        return np.nan
    try:
        # The complex() constructor is robust and can handle
        # strings like '1+2j' as well as numeric types.
        return complex(s)
    except (ValueError, TypeError):
        return np.nan
