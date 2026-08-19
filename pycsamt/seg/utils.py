# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

"""
SEG-EDI helpers: minimal parsing and serialization utilities.

Utility helpers shared by *pycsamt.seg* sub-modules.
The public API is intentionally small – only helpers that are
useful **outside** the core package are exposed.

"""

from __future__ import annotations

import math
import re
import shutil
import sys
import warnings
from collections.abc import Iterable, Mapping, Sequence, Sized
from pathlib import Path
from typing import (
    TYPE_CHECKING,
    Any,
    TextIO,
    TypeVar,
    Union,
)

import numpy as np

from ..exceptions import FileHandlingError
from ..log.logger import get_logger

StrLike = Union[str, np.str_]
Container = Union[StrLike, Sequence[StrLike], np.ndarray]

# Generic TypeVar bound to *Edi* (forward-referenced to avoid
# import cycles when the module is imported very early).
if TYPE_CHECKING:
    from .edi import (
        Edi,  # only evaluated by Mypy/Pylance/pyright
    )

_Edi = TypeVar("_Edi", bound="Edi")


NUMBER_RE = re.compile(
    r"""
    ^[+\-]?(
        (\d+(\.\d*)?|\.\d+)([eEdD][+\-]?\d+)?   # float or sci
      | \d+                                     # int
    )$
    """,
    re.VERBOSE,
)

KVPAIR_RE = re.compile(
    r"""
    (?P<key>[A-Za-z_][A-Za-z0-9_]*)
    \s*=\s*
    (?P<val>
        "(?:[^"]*)"
      | '(?:[^']*)'
      | [^,\s]+
    )
    """,
    re.VERBOSE,
)

DMS_RE = re.compile(
    r"""
    ^\s*(?P<d>[-+]?\d+)\s*:\s*
        (?P<m>\d+)\s*:\s*
        (?P<s>\d+(?:\.\d*)?)\s*([NSEW])?\s*$
    """,
    re.VERBOSE | re.IGNORECASE,
)


logger = get_logger(__name__)


def _haversine(
    lat1: float,
    lon1: float,
    lat2: float,
    lon2: float,
    radius_km: float = 6_371.0,
) -> float:
    """Great-circle distance ( *km* )."""
    lat1, lon1, lat2, lon2 = map(math.radians, (lat1, lon1, lat2, lon2))
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    a = (
        math.sin(dlat / 2.0) ** 2.0
        + math.cos(lat1) * math.cos(lat2) * math.sin(dlon / 2.0) ** 2.0
    )
    return 2.0 * radius_km * math.asin(math.sqrt(a))


def _euclidean(
    lat1: float,
    lon1: float,
    lat2: float,
    lon2: float,
) -> float:
    """Rough planar distance in *deg* units."""
    return math.hypot(lat2 - lat1, lon2 - lon1)


_DIST_METRICS = {
    "haversine": _haversine,
    "cartesian": _euclidean,
    "euclidean": _euclidean,
}


def sort_edis_by_location(
    edi_objs: Iterable[_Edi],
    *,
    by: str = "index",
    method: str = "strict",
    metric: str = "cartesian",
    reference: tuple[float, float] | None = None,
    verbose: int = 0,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Order *EDI* sites according to geographic criteria.

    Parameters
    ----------
    edi_objs :
        Iterable of :class:`~pycsamt.seg.edi.Edi` objects.
    by :
        Sorting key:

        * ``"ll"`` or ``"lonlat"`` – ascending lon then lat;
        * ``"latlon"`` – ascending lat then lon;
        * ``"distance"`` – distance from *reference*;
        * ``"name"`` – basename of ``.edi`` file;
        * ``"dataid"`` – :attr:`~pycsamt.seg.heads.Head.dataid`;
        * ``"index"`` – first integer found in the filename
          (default, mimics legacy behaviour).
    method :
        ``"strict"`` (default) uses the *selected* key **only**.
        ``"naive"`` falls back to filename order when information
        is ambiguous (identical coordinates).
    metric :
        Distance metric when *by* is ``"distance"`` –
        one of ``"cartesian"`` or ``"haversine"``.
    reference :
        ``(lat, lon)`` of the origin for distance sorting.
        When *None* the first site becomes the origin.
    verbose :
        0 → silent, ≥1 prints a quick recap through
        :pyfunc:`~pycsamt.seg.stats.quick_edi_stats`.

    Returns
    -------
    edi_sorted, names :
        *edi_sorted* is a ``(N, )`` *numpy* array of *EDI* objects;
        *names* the corresponding filenames (without directory).

    Notes
    -----
    The function **never** mutates the input list.

    Examples
    --------
    >>> from pycsamt.seg.utils import sort_edis_by_location
    >>> edis, names = sort_edis_by_location(my_edi_list, by="lonlat")
    """
    edi_list = list(edi_objs)
    if not edi_list:
        raise FileHandlingError("no EDI objects provided")

    metric = str(metric).lower()
    if metric not in _DIST_METRICS:
        raise ValueError(f"unknown metric '{metric}'")

    by = str(by).lower()
    allowed = {
        "ll",
        "lonlat",
        "latlon",
        "distance",
        "name",
        "dataid",
        "index",
    }
    if by not in allowed:
        raise ValueError(f"'by' must be in {sorted(allowed)}, got '{by}'")

    ref_lat, ref_lon = (
        (
            reference
            if reference is not None
            else (edi_list[0].head.lat, edi_list[0].head.lon)
        )
        if by == "distance"
        else (None, None)
    )

    # pre-extract criteria
    recs: list[dict] = []
    _idx_re = re.compile(r"\d+")
    for edi in edi_list:
        file_name = Path(edi.path or "").name
        idx_match = _idx_re.search(file_name)
        recs.append(
            dict(
                obj=edi,
                name=file_name,
                dataid=edi.head.dataid or file_name,
                index=int(idx_match.group()) if idx_match else 0,
                lon=edi.head.lon,
                lat=edi.head.lat,
            )
        )

    # compute distance if requested
    if by == "distance":
        dist_fun = _DIST_METRICS[metric]
        for r in recs:
            r["distance"] = dist_fun(ref_lat, ref_lon, r["lat"], r["lon"])

    # define sorting key
    if by in {"ll", "lonlat"}:

        def key_fun(r):
            return (r["lon"], r["lat"])

    elif by == "latlon":

        def key_fun(r):
            return (r["lat"], r["lon"])

    else:

        def key_fun(r):
            return r[by]

    recs_sorted = sorted(recs, key=key_fun)

    if method != "strict" and by in {"ll", "lonlat", "latlon"}:
        # naive fallback: preserve original order for ties
        unique_keys = {key_fun(r) for r in recs_sorted}
        if len(unique_keys) < len(recs_sorted):
            recs_sorted = sorted(
                recs_sorted,
                key=lambda r: (key_fun(r), edi_list.index(r["obj"])),
            )

    edi_sorted = np.array([r["obj"] for r in recs_sorted], dtype=object)
    names = np.array([r["name"] for r in recs_sorted], dtype=str)

    if verbose:
        quick_edi_stats(
            total=len(edi_list),
            ok=len(edi_sorted),
            label="sort",
        )

    return edi_sorted, names


def _safe_number(s: str) -> int | float | str:
    """
    Try to parse a string as int/float (including E-format). If it
    fails, return the original string (stripped of quotes).
    """
    if s is None:
        return s
    s = s.strip()
    # strip quotes
    if (len(s) >= 2) and ((s[0] == s[-1]) and s[0] in ("'", '"')):
        s = s[1:-1]

    if NUMBER_RE.match(s):
        # Replace Fortran 'D' exponent with 'E'
        s_norm = s.replace("d", "e").replace("D", "E")
        try:
            f = float(s_norm)
            # Cast to int if integral
            return (
                int(f)
                if f.is_integer()
                and ("." not in s_norm and "e" not in s_norm.lower())
                else f
            )
        except ValueError:
            return s  # fall back to string
    return s


def _dms_to_deg(s: str) -> float | None:
    """
    Convert 'DD:MM:SS(.sss)' (optionally ending with N/S/E/W) to degrees.

    Returns ``None`` if the input does not look like DMS.
    """
    m = DMS_RE.match(s.strip())
    if not m:
        return None
    d = float(m.group("d"))
    mnt = float(m.group("m"))
    sec = float(m.group("s"))
    hemi = m.group(4)  # optional N/S/E/W
    val = abs(d) + mnt / 60.0 + sec / 3600.0
    if d < 0:
        val = -val
    if hemi:
        hemi = hemi.upper()
        if hemi in ("S", "W"):
            val = -abs(val)
        else:
            val = abs(val)
    return val


def dms_to_decimal_fallback(value: Any, hem: str | None = None) -> float:
    """Dependency-free DMS-string-to-decimal-degree conversion.

    Used only as a fallback when :func:`pycsamt.gis.utils.dms_to_decimal`
    cannot be imported (``seg/heads.py`` and ``seg/meas.py``), so it must
    not import anything beyond the standard library. Not a replacement
    for :func:`_dms_to_deg`, which serves a different lenient-parsing
    call site with its own ``None``-on-mismatch contract.
    """
    text = str(value).strip().upper()
    sign = -1.0 if text.startswith("-") else 1.0
    if text.endswith(("S", "W")):
        sign = -1.0
    text = text.rstrip("NSEW").lstrip("+-")
    parts = [part for part in re.split(r"[:\s]+", text) if part]
    if len(parts) == 1:
        return sign * float(parts[0])
    degree = float(parts[0])
    minute = float(parts[1]) if len(parts) > 1 else 0.0
    second = float(parts[2]) if len(parts) > 2 else 0.0
    return sign * (degree + minute / 60.0 + second / 3600.0)


def decimal_to_dms_fallback(value: Any, hem: str | None = None) -> str:
    """Dependency-free decimal-degree-to-DMS-string counterpart of
    :func:`dms_to_decimal_fallback`."""
    numeric = float(value)
    sign = "-" if numeric < 0 else ""
    numeric = abs(numeric)
    degree = int(numeric)
    rem = (numeric - degree) * 60.0
    minute = int(rem)
    second = (rem - minute) * 60.0
    return f"{sign}{degree:02d}:{minute:02d}:{second:08.5f}"


def parse_kv_pairs(s: str) -> dict[str, Any]:
    """
    Parse `KEY=VALUE` pairs from a single line.

    Handles quoted strings, numeric literals (including Fortran
    'D' exponents), and D:M:S lat/long. Returns a dictionary with
    best-effort typing.

    Examples
    --------
    >>> parse_kv_pairs("ID=10011.001 CHTYPE=HX X=0.0 Y=0.0 AZM=0")
    {'ID': 10011.001, 'CHTYPE': 'HX', 'X': 0.0, 'Y': 0.0, 'AZM': 0}
    """
    out: dict[str, Any] = {}
    for m in KVPAIR_RE.finditer(s):
        key = m.group("key")
        val_raw = m.group("val").strip()
        # DMS?
        as_dms = _dms_to_deg(val_raw)
        out[key] = _safe_number(val_raw) if as_dms is None else as_dms
    return out


# Public: measurement parser
def gather_measurement_key_value_with_str_parser(
    lines: Iterable[str],
) -> list[dict[str, Any]]:
    """
    Parse measurement definitions from lines.

    This recognizes blocks like::

        >HMEAS ID=... CHTYPE=HX X=... Y=... Z=... AZM=...
        >EMEAS ID=... CHTYPE=EX X=... Y=... X2=... Y2=...

    and returns a list of dicts with a mandatory ``"KIND"`` key
    (``"HMEAS"`` or ``"EMEAS"``) plus all parsed key/value pairs.

    Parameters
    ----------
    lines : iterable of str
        Raw lines from the EDI file.

    Returns
    -------
    list of dict
        Each dict contains ``"KIND"`` and the parsed keys.

    Notes
    -----
    - DMS coordinates (e.g. ``26:35:11.0``) are converted to
      decimal degrees.
    - Quoted strings are unquoted.
    - Numeric tokens (including Fortran ``D`` exponents) are
      converted to numbers.

    Examples
    --------
    >>> src = [
    ...     ">HMEAS ID=10011.001 CHTYPE=HX X=0.0 Y=0.0 AZM=0",
    ...     ">EMEAS ID=10014.001 CHTYPE=EX X=-10.0 Y=0.0 X2=10.0 Y2=0.0",
    ... ]
    >>> gather_measurement_key_value_with_str_parser(src)[0]["CHTYPE"]
    'HX'
    """
    out: list[dict[str, Any]] = []
    for raw in lines:
        s = raw.strip()
        if not s.startswith(">"):
            continue
        if s.upper().startswith(">HMEAS"):
            payload = s[6:].strip()  # after ">HMEAS"
            item = {"KIND": "HMEAS", **parse_kv_pairs(payload)}
            out.append(item)
        elif s.upper().startswith(">EMEAS"):
            payload = s[6:].strip()
            item = {"KIND": "EMEAS", **parse_kv_pairs(payload)}
            out.append(item)
    return out


# Formatting helpers for writer
def _format_block_numbers(
    arr: Sequence[float],
    *,
    per_line: int = 6,
    fmt: str = "{: .6E}",
    indent: int = 2,
) -> str:
    """
    Format a 1-D sequence of numbers into a block like the EDI
    examples (6 entries per line, scientific format).
    """
    buf: list[str] = []
    line: list[str] = []
    for i, v in enumerate(arr):
        line.append(fmt.format(float(v)))
        if (i + 1) % per_line == 0:
            buf.append(" " * indent + " ".join(line))
            line = []
    if line:
        buf.append(" " * indent + " ".join(line))
    return "\n".join(buf)


def quick_edi_stats(
    *,
    total: int,
    ok: int,
    label: str = "EDI",
    width: int | None = None,
) -> None:
    """One-line recap similar to legacy *show_stats*."""
    term = shutil.get_terminal_size((80, 20)).columns
    width = width or max(60, min(120, term))
    bar = "~" * width

    rate = (ok / total) * 100 if total else 0.0
    print(bar)
    print(
        f"{label:<15} read  : {ok:>6d}/{total:<6d}  —  success {rate:6.2f} %"
    )
    print(bar)


def _format_kv(k: str, v):
    kU = k.upper()
    # leave numerics bare; quote strings only if needed
    if isinstance(v, (int, float)):
        return f"{k}={v}"
    s = "" if v is None else str(v)
    # LAT/LONG/ELEV commonly numeric; if a string slipped in and
    # parses as a number, keep it bare.
    try:
        if kU in {"LAT", "LONG", "ELEV"}:
            fv = float(s.replace("D", "E").replace("d", "e"))
            return f"{k}={fv}"
    except Exception:
        pass
    # quote if whitespace or quotes present
    return f'{k}="{s.replace(chr(34), "")}"'


def _ensure_1d(a: Sequence[float] | np.ndarray) -> np.ndarray:
    v = np.asarray(a).ravel()
    return v.astype(float, copy=False)


def _quote(s: Any) -> str:
    if s is None:
        return '""'
    s = str(s)
    if any(ch.isspace() for ch in s) or any(ch in s for ch in ('"',)):
        return '"' + s.replace('"', "") + '"'
    return s


def minimum_parser_to_write_edi(obj: Mapping[str, Any]) -> str:
    """
    Write a minimal, valid SEG-EDI text from a structured mapping.

    This is deliberately small and conservative; it supports the
    headers and blocks present in your examples and is intended as
    the *last* step of the pipeline (the higher-level Edi class
    should prepare this structure).

    Expected structure (keys are optional unless noted):

    - ``head`` (mapping):
        Keys such as ``DATAID``, ``ACQBY``, ``ACQDATE``,
        ``FILEDATE``, ``PROSPECT``, ``LAT``, ``LONG``, ``ELEV``,
        ``STDVERS``, ``PROGVERS``, ``PROGDATE``, ``MAXSECT``,
        ``EMPTY``. Values are serialized with ``KEY=VALUE``
        (quoted if needed).
    - ``info`` (str or list[str]): free text under ``>INFO``.
    - ``definemeas`` (mapping):
        Keys like ``MAXCHAN``, ``MAXRUN``, ``MAXMEAS``, ``UNITS``,
        ``REFTYPE``, ``REFLAT``, ``REFLONG``, ``REFELEV``.
    - ``measurements`` (list[dict]): entries produced by
        :func:`gather_measurement_key_value_with_str_parser`.
        Each must include ``"KIND"``.
    - ``mtsect`` (mapping):
        Must include ``SECTID`` (str) and ``NFREQ`` (int).
        May include channel bindings (``HX``, ``HY``, ``HZ``,
        ``EX``, ``EY``, ``RX``, ``RY``).
    - ``freq`` (1-D array-like of float): frequency list.
    - ``zrot`` (1-D array-like of float): rotation angles.
    - ``impedance`` (mapping): optional numeric blocks with any of
        the standard keys (``ZXXR``, ``ZXXI``, ``ZXX.VAR``, ...).
    - ``resistivity`` (mapping): optional ``RHOXY``, ``RHOXY.ERR``,
        ``RHOYX``, ``RHOYX.ERR``, ``PHSXY``, ``PHSXY.ERR``,
        ``PHSYX``, ``PHSYX.ERR``.
    - ``coherence`` (list[tuple[str, Sequence[float]]]): optional
        pairs of a labelled preamble (e.g.
        ``'COH MEAS1=1 MEAS2=5 ROT=NORTH'``) and the numeric block.

    Parameters
    ----------
    obj : mapping
        Structured content as described above.

    Returns
    -------
    str
        Text of the EDI file.

    Notes
    -----
    - Numeric blocks are formatted at 6 values per line using
      scientific notation.
    - If counts are present, ``// N`` is emitted; otherwise they
      are inferred from the array length.
    - The function does **not** touch data ordering; it is up to
      the caller to ensure ``freq``, ``zrot`` and component
      arrays are aligned.

    Examples
    --------
    >>> edi_text = minimum_parser_to_write_edi(
    ...     {
    ...         "head": {
    ...             "DATAID": "E1_2",
    ...             "STDVERS": "SEG 1.0",
    ...             "EMPTY": 1e32,
    ...         },
    ...         "info": "Processed by pyCSAMT",
    ...         "definemeas": {"MAXCHAN": 16, "UNITS": "M"},
    ...         "measurements": [
    ...             {"KIND": "HMEAS", "ID": 1.001, "CHTYPE": "HX"}
    ...         ],
    ...         "mtsect": {"SECTID": "E1_2", "NFREQ": 2, "HX": 1.001},
    ...         "freq": [7e4, 5.88e4],
    ...         "zrot": [0.0, 0.0],
    ...     }
    ... )
    >>> edi_text.startswith(">HEAD")
    True
    """
    lines: list[str] = []

    # >HEAD
    head = dict(obj.get("head") or {})
    if head:
        lines.append(">HEAD")
        for k, v in head.items():
            lines.append("  " + _format_kv(k, v))

        lines.append("")

    # >INFO
    info = obj.get("info")
    if info is not None:
        lines.append(">INFO")
        if isinstance(info, str):
            # keep multiline as-is, indent each line
            for ln in info.splitlines():
                lines.append(f"  {ln.rstrip()}")
        else:
            for ln in info:
                lines.append(f"  {str(ln).rstrip()}")
        lines.append("")

    # >=DEFINEMEAS
    definemeas = dict(obj.get("definemeas") or {})
    meas = list(obj.get("measurements") or [])
    if definemeas or meas:
        lines.append(">=DEFINEMEAS")
        for k, v in definemeas.items():
            if isinstance(v, str):
                lines.append(f"  {k}={_quote(v)}")
            else:
                lines.append(f"  {k}={v}")
        if meas:
            lines.append(
                ">!***CHANNELS USING ORIGINAL SITE LAYOUT. FOR ROTATIONS SEE ZROT***!"
            )
            for m in meas:
                kind = str(m.get("KIND", "")).upper()
                kv = {k: v for k, v in m.items() if k != "KIND"}
                payload = " ".join(
                    f"{k}={_quote(v)}" if isinstance(v, str) else f"{k}={v}"
                    for k, v in kv.items()
                )
                lines.append(f">{kind}  {payload}")
        lines.append("")

    # >=MTSECT
    mtsect = dict(obj.get("mtsect") or {})
    if mtsect:
        lines.append(">=MTSECT")
        for k, v in mtsect.items():
            if isinstance(v, str):
                lines.append(f"  {k} =   {_quote(v)}")
            else:
                lines.append(f"  {k} =   {v}")
        lines.append("")

    # >FREQ
    freq = obj.get("freq")
    if freq is not None:
        fv = _ensure_1d(freq)
        lines.append(">!****FREQUENCIES****!")
        lines.append(f">FREQ  //{fv.size}")
        lines.append(_format_block_numbers(fv, indent=2))
        lines.append("")

    # >ZROT
    zrot = obj.get("zrot")
    if zrot is not None:
        zv = _ensure_1d(zrot)
        lines.append(">!****IMPEDANCE ROTATION ANGLES****!")
        lines.append(f">ZROT  //{zv.size}")
        lines.append(_format_block_numbers(zv, indent=2))
        lines.append("")

    # > IMPEDANCE blocks
    imp: Mapping[str, Any] = obj.get("impedance") or {}
    for key in (
        "ZXXR",
        "ZXXI",
        "ZXX.VAR",
        "ZXYR",
        "ZXYI",
        "ZXY.VAR",
        "ZYXR",
        "ZYXI",
        "ZYX.VAR",
        "ZYYR",
        "ZYYI",
        "ZYY.VAR",
    ):
        if key in imp and imp[key] is not None:
            arr = _ensure_1d(imp[key])
            lines.append(f">{key} ROT=ZROT  //{arr.size}")
            lines.append(_format_block_numbers(arr, indent=2))
            lines.append("")

    # > RESISTIVITY / PHASE blocks (optional)
    rho: Mapping[str, Any] = obj.get("resistivity") or {}
    for key in (
        "RHOROT",
        "RHOXY",
        "RHOXY.ERR",
        "RHOYX",
        "RHOYX.ERR",
        "PHSXY",
        "PHSXY.ERR",
        "PHSYX",
        "PHSYX.ERR",
    ):
        if key in rho and rho[key] is not None:
            arr = _ensure_1d(rho[key])
            # Only RHOROT is a pure angle list; others usually use ROT=RHOROT.
            if key == "RHOROT":
                lines.append(f">{key}   // {arr.size}")
            else:
                lines.append(f">{key}  ROT=RHOROT // {arr.size}")
            lines.append(_format_block_numbers(arr, indent=2))
            lines.append("")

    # > COH blocks (optional, generic label + data)
    for item in obj.get("coherence") or []:
        if not item:
            continue
        label, data = item
        arr = _ensure_1d(data)
        lines.append(f">{label} // {arr.size}")
        lines.append(_format_block_numbers(arr, indent=2))
        lines.append("")

    lines.append(">END")
    return "\n".join(lines)


def strip_item(
    item_to_clean: Container | None,
    item: str | None = None,
    multi_space: int = 12,
) -> Container | None:
    """
    Strip a token repeatedly from both ends of strings.

    This utility cleans leading/trailing *repetitions* of a token
    (default: whitespace) in a flexible way. It accepts a scalar
    string, a list of strings, or a NumPy array of strings and
    returns the same container type.

    If the input (or all items within) becomes empty after
    sanitization, ``None`` is returned and a warning is issued.
    This mirrors the legacy behavior where completely blank values
    are treated as missing.

    Parameters
    ----------
    item_to_clean : {str, list of str, np.ndarray of str}, optional
        The text or collection of texts to sanitize. If ``None``,
        returns ``None``.
    item : str, optional
        The token to strip from both ends. If ``None`` or a blank
        string, standard whitespace stripping (``str.strip()``) is
        used. For multi-character tokens (e.g., ``"//"``), the
        token is removed as repeated *whole* substrings rather
        than character sets.
    multi_space : int, default=12
        Maximum repetition count to consider at each end for a
        multi-character token (upper bound in the regex quantifier).
        Must be a positive integer.

    Returns
    -------
    {str, list of str, np.ndarray of str} or None
        Sanitized output preserving the input container type, or
        ``None`` if content is effectively empty after cleaning.

    Notes
    -----
    * For ``item`` that is ``None``/blank, this function applies
      ``str.strip()`` (whitespace). For a non-blank ``item``, it
      removes repeated occurrences of that *exact* token from both
      ends using a compiled regular expression anchored at the start
      and end of the string.
    * Returning ``None`` for fully empty results keeps backward
      compatibility with legacy usage that treats blank fields as
      missing.

    Examples
    --------
    >>> strip_item("    ss_data   ")
    'ss_data'
    >>> strip_item(["  a  ", "   b"], item=None)
    ['a', 'b']
    >>> arr = np.array(["////name////", "////x"], dtype="<U16")
    >>> strip_item(arr, item="//")
    array(['name', 'x'], dtype='<U16')
    >>> strip_item(["      ", "   "]) is None
    True
    """
    # Fast exits and validations
    if item_to_clean is None:
        return None

    if not isinstance(multi_space, int) or multi_space <= 0:
        raise TypeError(
            f"'multi_space' must be a positive integer, got {multi_space!r}"
        )

    # Normalize to a list of Python strings for processing; track original type
    input_is_array = isinstance(item_to_clean, np.ndarray)
    input_is_scalar = isinstance(item_to_clean, (str, np.str_))

    if input_is_scalar:
        data_list: list[str] = [str(item_to_clean)]
        orig_dtype = None
    elif input_is_array:
        # Ensure 1D iteration over elements (even for object dtype)
        arr = np.asarray(item_to_clean)
        orig_dtype = arr.dtype
        data_list = [str(x) if x is not None else "" for x in arr.tolist()]
    elif isinstance(item_to_clean, Sequence):
        data_list = [str(x) if x is not None else "" for x in item_to_clean]
        orig_dtype = None
    else:
        raise TypeError(
            "item_to_clean must be a string, a list/sequence of strings, "
            "or a NumPy array of strings."
        )

    # Build the stripper
    token = (item or "").strip("\n\r")
    cleaned: list[str] = []

    if token == "":
        # Whitespace strip
        for s in data_list:
            cleaned.append(s.strip())
    else:
        # Strip repeated *token* (not char-set) from both ends
        # ^(?:token){1,N} | (?:token){1,N}$
        pat = re.compile(
            rf"^(?:{re.escape(token)}){{1,{multi_space}}}|"
            rf"(?:{re.escape(token)}){{1,{multi_space}}}$"
        )
        for s in data_list:
            # Remove at both ends (pattern matches both anchors)
            ns = pat.sub("", s)
            # Also trim surrounding whitespace that may remain
            cleaned.append(ns.strip())

    # Decide if the result is effectively empty
    if all(c == "" for c in cleaned):
        warnings.warn(
            "No data found for sanitization; returning None.",
            RuntimeWarning,
            stacklevel=2,
        )
        return None

    # Reconstruct original container type
    if input_is_scalar:
        return cleaned[0]
    if input_is_array:
        # Preserve dtype if possible; fallback to unicode
        try:
            return np.array(cleaned, dtype=orig_dtype)
        except Exception:
            return np.array(cleaned)
    # Sequence -> list
    return cleaned


def _len(obj: Any) -> int:
    """
    Robust length helper.

    Returns ``len(obj)`` if *obj* is sized, otherwise
    treats *obj* as a scalar (length = 1).
    """
    return len(obj) if isinstance(obj, Sized) else 1


def show_edi_stats(
    collected: int | Iterable[Any],
    succeeded: int | Iterable[Any],
    *,
    failed: int | Iterable[Any] | None = None,
    elapsed: float | None = None,
    title: str = "EDI",
    width: int = 72,
    sep: str = "~",
    stream: TextIO | None = None,
) -> None:
    """
    Pretty print collection statistics.

    Parameters
    ----------
    collected
        Number of items *attempted* or an iterable of those
        items (length will be taken).
    succeeded
        Number of items *successfully* processed or an
        iterable (length will be taken).
    failed, optional
        Number of failed items or an iterable.  If *None*
        it is inferred: ``failed = collected - succeeded``.
    elapsed, optional
        Wall-clock time in *seconds*.
    title
        Short label for the object type (``"EDI"`` by
        default).
    width
        Total line width for the banner.
    sep
        Character used to draw the banner.
    stream, optional
        Output stream (defaults to :pydata:`sys.stdout`).
    """
    stream = stream or sys.stdout
    n_total = _len(collected)
    n_ok = _len(succeeded)
    n_fail = _len(failed) if failed is not None else n_total - n_ok
    rate = 0.0 if n_total == 0 else 100.0 * n_ok / n_total

    # Format banner
    bar = sep * width
    print(bar, file=stream)

    rows: list[str] = []
    rows.append(f"{'Collected':<14}: {n_total:>8,d}")
    rows.append(f"{'Succeeded':<14}: {n_ok:>8,d}")
    rows.append(f"{'Failed':<14}: {n_fail:>8,d}")
    rows.append(f"{'Success rate':<14}: {rate:>7.2f} %")
    if elapsed is not None:
        rows.append(f"{'Elapsed':<14}: {elapsed:>8.2f} s")

    # Centre title if it fits
    hdr = f"{title.upper()} STATISTICS"
    hdr_line = hdr.center(width)
    print(hdr_line, file=stream)

    for line in rows:
        print(line.center(width), file=stream)

    print(bar, file=stream)
