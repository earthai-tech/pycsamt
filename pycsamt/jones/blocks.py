# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, List, Optional, Sequence
from typing import Tuple
import math

import numpy as np

from .base import JComponentBase
from .config import ( 
    ENCODING_DEFAULT,
    DTYPE_SPECS,
    KIND_RHO_PHI, 
    KIND_COMPLEX_TF,
    RE_BLANK, 
    MISSING_FLOAT
)
# --- add at top with other imports ---
from .config import (
    RE_STATION,
    RE_DATATYPE_UNITS,
    RE_NPOINTS,
    RE_INFO,
    RE_COMMENT,
)
from .utils import parse_datatype_units, JParseError

from .utils import iter_lines
from .heads import Head

__all__ = [
    "RRow", "TFRow", "JBlock",
    "RBlock", "TFBlock", "JBlocks",
]

@dataclass
class RRow:
    period: float
    rho: float
    pha: float
    rhomax: float
    rhomin: float
    phamax: float
    phamin: float
    wrho: float
    wpha: float
    rej: bool = False

@dataclass
class TFRow:
    period: float
    real: float
    imag: float
    error: float
    weight: float
    rej: bool = False


class JBlock(JComponentBase):
    r"""
    Abstract base for J-format data blocks.

    A block ties a parsed header (:class:`Head`) to a homogeneous
    sequence of rows (:class:`RRow` or :class:`TFRow`). Subclasses
    implement the row parser, normalization, serialization, and
    simple QA summaries.

    Parameters
    ----------
    head : Head
        Parsed header triple (station, data-type, count).
    rows : sequence of rows
        Concrete row records for this block family.
    verbose : int, optional
        Verbosity for warnings during parsing and writing.

    Attributes
    ----------
    head : Head
        Block header with ``station``, ``kind``, ``comp``, ``n``.
    nrows : int
        Number of parsed data rows.
    station, kind, comp, units : str or None
        Convenience views onto ``head.dtype``.
    columns : tuple of str
        Column names for the structured numeric view.
    shape : (int, int)
        ``(nrows, ncols)`` derived from the current rows.

    Methods
    -------
    from_file(j_fn, *, verbose=0)
        Build a block from a file slice. The concrete subclass
        is chosen from the header's ``kind``.
    from_lines(lines, *, verbose=0)
        Parse from in-memory lines (header + body).
    read((head, body_lines))
        Populate the block from an existing :class:`Head` and
        the subsequent body lines. Marks the instance as read.
    write()
        Serialize header and all body rows into text.
    to_numpy()
        Structured array view with numeric normalization
        (period sign, missing sentinels, rejection flags).
    to_dataframe()
        Optional ``pandas`` representation (if available).
    qa_summary()
        Small dict of quick QA metrics for this block.

    Notes
    -----
    The base class is intentionally minimal to keep the parser
    light. Subclasses define the row regex, column names, and
    row-level normalization logic.

    Examples
    --------
    >>> # Pseudocode; subclass chosen from header
    >>> blk = JBlock.from_lines(['S01', 'RXY', '1', '...'])
    >>> blk.station, blk.nrows
    ('S01', 1)

    See Also
    --------
    RBlock, TFBlock : Concrete block families.
    Heads : Header-only convenience container.

    References
    ----------
    .. [1] A. G. Jones (1994). Magnetotelluric data file
       J-format, version 2.0.
    """



    def __init__(
        self,
        head: Optional[Head] = None,
        *,
        rows: Optional[Sequence[Any]] = None,
        verbose: int = 0,
    ) -> None:
        super().__init__(verbose=verbose)
        self.head: Head | None = head
        self.rows: List[Any] = list(rows or [])

    @classmethod
    def from_file(
        cls,
        j_fn: str | Path,
        *,
        verbose: int = 0,
    ) -> "JBlock":
        lines = list(iter_lines(j_fn, encoding=ENCODING_DEFAULT))
        return cls.from_lines(lines, verbose=verbose)

    @classmethod
    def from_lines(
        cls,
        lines: Sequence[str],
        *,
        verbose: int = 0,
    ) -> "JBlock":
        head, body = _extract_first_head_and_body(lines)
        blk = _block_factory(head, verbose=verbose)
        return blk.read((head, body))

    def read(
        self,
        data: Tuple[Head, Sequence[str]] | None = None,
    ) -> "JBlock":
        raise NotImplementedError("override in subclass")

    def write(self) -> List[str]:
        raise NotImplementedError("override in subclass")

    # ---- normalization & QA (common facade) ----
    def normalize(self) -> "JBlock":
        """
        Apply row-level normalization:

        - negative period => frequency (Hz) to period (s)
        - missing sentinel (-999) => NaN
        - set rejection flags per rules
        """
        return self  # implemented by subclasses

    def qa_summary(self) -> dict[str, Any]:
        """
        Quick QA metrics (counts and basic stats).
        """
        return {}  

    # ---- exports ----
    def to_numpy(self) -> dict[str, np.ndarray]:
        raise NotImplementedError("override in subclass")

    def to_dataframe(self):
        try:
            import pandas as pd  # type: ignore
        except Exception as exc:  # pragma: no cover
            raise ImportError(
                "pandas is required for to_dataframe()"
            ) from exc
        data = self.to_numpy()
        return pd.DataFrame(data)


    @property
    def nrows(self) -> int:
        return len(self.rows)

    @property
    def station(self) -> str | None:
        return ( 
            None if self.head is None 
            else self.head.station
            )
    
    @property
    def kind(self) -> str | None:
        return ( 
            None if self.head is None 
            else self.head.kind
            )
    
    @property
    def comp(self) -> str | None:
        return ( 
            None if self.head is None 
            else self.head.comp
        )
    
    @property
    def units(self) -> str | None:
        return ( 
            None if self.head is None 
            else self.head.units
        )
    
    @property
    def has_units(self) -> bool:
        return bool(self.units)
    
    @property
    def columns(self) -> tuple[str, ...]:
        # default; subclasses override if needed
        return tuple([])
    
    @property
    def shape(self) -> tuple[int, int]:
        return (self.nrows, len(self.columns))
    
    @property
    def periods(self):
        # cheap view; avoids pandas dep
        a = self.to_numpy()
        return a["period"]

    def __str__(self) -> str:
        st = self.head.station if self.head else None
        dt = getattr(self.head, "dtype", None)
        k = getattr(dt, "kind", None)
        c = getattr(dt, "comp", None)
        return (
            f"{self.__class__.__name__}(station={st!r}, "
            f"dtype={k}{c}, nrows={self.nrows})"
        )

class RBlock(JBlock):
    r"""
    Resistivity/phase (*R/S*) block implementation.

    Parses rows of the form ``period rho pha rhomax rhomin phamax
    phamin wrho wpha``. Applies J-format rules to normalize
    period sign and to mark rejected estimates.

    Parameters
    ----------
    head : Head
        Header whose kind is ``'R'`` or ``'S'``. The component
        encodes the tensor entry or average (e.g. ``RXY``, ``RDE``).
    rows : sequence of RRow, optional
        Prebuilt rows. Usually left to the parser.
    verbose : int, optional
        Verbosity for warnings during parsing.

    Attributes
    ----------
    columns : tuple of str
        ``('period', 'rho', 'pha', 'rhomax', 'rhomin',
        'phamax', 'phamin', 'wrho', 'wpha', 'rej')``.

    Notes
    -----
    Normalization rules:

    - If input period is negative, the value stores frequency in
      Hz and is converted to period (s).
    - ``rho < 0`` or ``wrho < 0`` marks the row as rejected.
    - ``-999.0`` sentinels are converted to ``NaN`` in numeric
      views.

    Examples
    --------
    >>> lines = [
    ...   'S01', 'RXY', '2',
    ...   '-1.0 100 45 110 90 50 40 1 1',
    ...   ' 2.0 -5  30  35 25 40 20 1 1',
    ... ]
    >>> blk = RBlock.from_lines(lines)
    >>> a = blk.to_numpy()
    >>> a['period'][0]  # 1/Hz -> s
    1.0
    >>> a['rej'][1]     # rho < 0
    True

    See Also
    --------
    TFBlock : Transfer-function block family.
    RRow : Single row record used by this block.

    References
    ----------
    .. [1] A. G. Jones (1994). Magnetotelluric data file
       J-format, version 2.0.
    """


    def read(
        self,
        data: Tuple[Head, Sequence[str]] | None = None,
    ) -> "RBlock":
        if data is None:
            raise ValueError("data is required")
        head, body = data
        spec = DTYPE_SPECS[head.dtype.kind]
        rows: List[RRow] = []
        n = head.n or 0
        i = 0
        for ln in body:
            if RE_BLANK.match(ln):
                continue
            m = spec.regex.match(ln)
            if not m:
                break
            d = m.groupdict()
            rows.append(
                RRow(
                    period=float(d["p"]),
                    rho=float(d["rho"]),
                    pha=float(d["pha"]),
                    rhomax=float(d["rhomax"]),
                    rhomin=float(d["rhomin"]),
                    phamax=float(d["phamax"]),
                    phamin=float(d["phamin"]),
                    wrho=float(d["wrho"]),
                    wpha=float(d["wpha"]),
                )
            )
            i += 1
            if n and i >= n:
                break
        self.head = head
        self.rows = rows
        self._mark_read(True)
        # apply default normalization
        return self.normalize()

    def normalize(self) -> "RBlock":
        for r in self.rows:
            # period: if negative => frequency (Hz)
            if r.period < 0:
                f = abs(r.period)
                r.period = _safe_div(1.0, f)
            # missing sentinels => NaN
            r.rho = _nanz(r.rho)
            r.pha = _nanz(r.pha)
            r.rhomax = _nanz(r.rhomax)
            r.rhomin = _nanz(r.rhomin)
            r.phamax = _nanz(r.phamax)
            r.phamin = _nanz(r.phamin)
            r.wrho = _nanz(r.wrho)
            r.wpha = _nanz(r.wpha)
            # rejection (rho<0 or wrho<0)
            rej = False
            if isinstance(r.rho, float) and not math.isnan(r.rho):
                rej = rej or (r.rho < 0.0)
            if isinstance(r.wrho, float) and not math.isnan(r.wrho):
                rej = rej or (r.wrho < 0.0)
            r.rej = bool(rej)
        return self

    def to_numpy(self) -> dict[str, np.ndarray]:
        n = self.nrows
        arr = {
            "period": np.empty(n, float),
            "rho": np.empty(n, float),
            "pha": np.empty(n, float),
            "rhomax": np.empty(n, float),
            "rhomin": np.empty(n, float),
            "phamax": np.empty(n, float),
            "phamin": np.empty(n, float),
            "wrho": np.empty(n, float),
            "wpha": np.empty(n, float),
            "rej": np.empty(n, bool),
        }
        for i, r in enumerate(self.rows):
            arr["period"][i] = r.period
            arr["rho"][i] = r.rho
            arr["pha"][i] = r.pha
            arr["rhomax"][i] = r.rhomax
            arr["rhomin"][i] = r.rhomin
            arr["phamax"][i] = r.phamax
            arr["phamin"][i] = r.phamin
            arr["wrho"][i] = r.wrho
            arr["wpha"][i] = r.wpha
            arr["rej"][i] = r.rej
        return arr

    def qa_summary(self) -> dict[str, Any]:
        a = self.to_numpy()
        n = self.nrows
        rej = int(a["rej"].sum())
        keep = n - rej
        finite = np.isfinite(a["rho"]).sum()
        pmin = np.nanmin(a["period"]) if n else np.nan
        pmax = np.nanmax(a["period"]) if n else np.nan
        return {
            "kind": "R/S",
            "nrows": n,
            "rejected": rej,
            "kept": keep,
            "finite_rho": int(finite),
            "period_min": float(pmin),
            "period_max": float(pmax),
        }

    def write(self) -> List[str]:
        if not self.head or self.head.n is None:
            raise ValueError("missing head")
        out: List[str] = []
        out.extend(self.head.write())
        for r in self.rows:
            out.append(
                (
                    "{: .6e} {: .6e} {: .6e} {: .6e} {: .6e} "
                    "{: .6e} {: .6e} {: .6e} {: .6e}"
                ).format(
                    r.period, r.rho, r.pha, r.rhomax, r.rhomin,
                    r.phamax, r.phamin, r.wrho, r.wpha
                ).strip()
            )
        return out

    @property
    def columns(self) -> tuple[str, ...]:
        return (
            "period", "rho", "pha", "rhomax", "rhomin",
            "phamax", "phamin", "wrho", "wpha", "rej",
        )


class TFBlock(JBlock):
    r"""
    Transfer-function (*Z/Q/C/T*) block implementation.

    Parses rows of the form ``period real imag error weight``.
    Applies J-format rules to normalize period sign and to mark
    rejected estimates.

    Parameters
    ----------
    head : Head
        Header whose kind is one of ``'Z'``, ``'Q'``, ``'C'``,
        or ``'T'``.
    rows : sequence of TFRow, optional
        Prebuilt rows. Usually left to the parser.
    verbose : int, optional
        Verbosity for warnings during parsing.

    Attributes
    ----------
    columns : tuple of str
        ``('period', 'real', 'imag', 'error', 'weight', 'rej')``.

    Notes
    -----
    Normalization rules:

    - If input period is negative, the value stores frequency in
      Hz and is converted to period (s).
    - ``weight < 0`` marks the row as rejected.
    - ``-999.0`` sentinels are converted to ``NaN`` in numeric
      views.

    Examples
    --------
    >>> lines = [
    ...   'S02', 'ZXY SI', '2',
    ...   '-1.0e+1 1.0 -2.0 0.1 1.0',
    ...   ' 5.0    -999 3.0 -999 -1.0',
    ... ]
    >>> blk = TFBlock.from_lines(lines)
    >>> a = blk.to_numpy()
    >>> round(a['period'][0], 3)
    0.1
    >>> a['rej'][1]     # weight < 0
    True

    See Also
    --------
    RBlock : Resistivity/phase block family.
    TFRow : Single row record used by this block.

    References
    ----------
    .. [1] A. G. Jones (1994). Magnetotelluric data file
       J-format, version 2.0.
    """
    def read(
        self,
        data: Tuple[Head, Sequence[str]] | None = None,
    ) -> "TFBlock":
        if data is None:
            raise ValueError("data is required")
        head, body = data
        spec = DTYPE_SPECS[head.dtype.kind]
        rows: List[TFRow] = []
        n = head.n or 0
        i = 0
        for ln in body:
            if RE_BLANK.match(ln):
                continue
            m = spec.regex.match(ln)
            if not m:
                break
            d = m.groupdict()
            rows.append(
                TFRow(
                    period=float(d["p"]),
                    real=float(d["real"]),
                    imag=float(d["imag"]),
                    error=float(d["error"]),
                    weight=float(d["weight"]),
                )
            )
            i += 1
            if n and i >= n:
                break
        self.head = head
        self.rows = rows
        self._mark_read(True)
        return self.normalize()

    def normalize(self) -> "TFBlock":
        for r in self.rows:
            if r.period < 0:
                f = abs(r.period)
                r.period = _safe_div(1.0, f)
            r.real = _nanz(r.real)
            r.imag = _nanz(r.imag)
            r.error = _nanz(r.error)
            r.weight = _nanz(r.weight)
            rej = False
            if isinstance(r.weight, float) and not math.isnan(r.weight):
                rej = rej or (r.weight < 0.0)
            r.rej = bool(rej)
        return self

    def to_numpy(self) -> dict[str, np.ndarray]:
        n = self.nrows
        arr = {
            "period": np.empty(n, float),
            "real": np.empty(n, float),
            "imag": np.empty(n, float),
            "error": np.empty(n, float),
            "weight": np.empty(n, float),
            "rej": np.empty(n, bool),
        }
        for i, r in enumerate(self.rows):
            arr["period"][i] = r.period
            arr["real"][i] = r.real
            arr["imag"][i] = r.imag
            arr["error"][i] = r.error
            arr["weight"][i] = r.weight
            arr["rej"][i] = r.rej
        return arr

    def qa_summary(self) -> dict[str, Any]:
        a = self.to_numpy()
        n = self.nrows
        rej = int(a["rej"].sum())
        keep = n - rej
        pmin = np.nanmin(a["period"]) if n else np.nan
        pmax = np.nanmax(a["period"]) if n else np.nan
        rms = np.sqrt(
            np.nanmean(np.square(a["error"]))
        ) if n else np.nan
        return {
            "kind": "TF",
            "nrows": n,
            "rejected": rej,
            "kept": keep,
            "period_min": float(pmin),
            "period_max": float(pmax),
            "error_rms": float(rms),
        }

    def write(self) -> List[str]:
        if not self.head or self.head.n is None:
            raise ValueError("missing head")
        out: List[str] = []
        out.extend(self.head.write())
        for r in self.rows:
            out.append(
                (
                    "{: .6e} {: .6e} {: .6e} {: .6e} {: .6e}"
                ).format(
                    r.period, r.real, r.imag, r.error, r.weight
                ).strip()
            )
        return out

    @property
    def columns(self) -> tuple[str, ...]:
        return (
            "period", "real", "imag", "error", "weight",
            "rej",
        )

class JBlocks(JComponentBase):
    r"""
    Container for a sequence of parsed J data blocks.

    Provides discovery from raw text, iteration, selection by
    family or component, serialization, and quick QA roll-ups
    across blocks.

    Parameters
    ----------
    blocks : sequence of JBlock, optional
        Prebuilt blocks. Usually created by the parser.
    verbose : int, optional
        Verbosity for warnings during parsing.

    Attributes
    ----------
    blocks : list of JBlock
        Parsed blocks in file order.
    n : int
        Number of blocks.
    stations : list of str
        All station identifiers seen across blocks.
    kinds : list of str
        The ``kind`` tokens across all blocks (e.g. ``'R'``, ``'Z'``).

    Methods
    -------
    from_file(j_fn, *, verbose=0)
        Parse all blocks found in a file.
    from_lines(lines, *, verbose=0)
        Parse blocks from an in-memory sequence of lines.
    read(lines)
        Replace current content with blocks parsed from lines.
    write()
        Serialize every block back to text.
    select(kind=None, comp=None)
        Return blocks filtered by family and/or component.
    period_range()
        Global period min/max across blocks (ignores ``NaN``).
    qa_summary()
        List of QA dicts, one per block.

    Notes
    -----
    Header discovery is tolerant to real-world quirks, including
    optional azimuth on the station line and the
    count-before-dtype pattern seen in some files.

    Examples
    --------
    >>> col = JBlocks.from_file('data/j/kb0-s001.txt')
    >>> col.n >= 1
    True
    >>> [b.station for b in col.select(kind='R')]
    ['KB0001', ...]

    See Also
    --------
    Head, Heads : Header parsers used during discovery.
    RBlock, TFBlock : Concrete block families contained here.

    References
    ----------
    .. [1] A. G. Jones (1994). Magnetotelluric data file
       J-format, version 2.0.
    """

    def __init__(
        self,
        blocks: Optional[Sequence[JBlock]] = None,
        *,
        verbose: int = 0,
    ) -> None:
        super().__init__(verbose=verbose)
        self.blocks: List[JBlock] = list(blocks or [])

    @property
    def n(self) -> int:
        return len(self.blocks)

    @classmethod
    def from_file(
        cls,
        j_fn: str | Path,
        *,
        verbose: int = 0,
    ) -> "JBlocks":
        lines = list(iter_lines(j_fn, encoding=ENCODING_DEFAULT))
        inst = cls.from_lines(lines, verbose=verbose)
        # inst is already marked by from_lines (see below)
        return inst

    @classmethod
    def from_lines(
        cls,
        lines: Sequence[str],
        *,
        verbose: int = 0,
    ) -> "JBlocks":
        blocks = _extract_all_blocks(lines, verbose=verbose)
        inst= cls(blocks=blocks, verbose=verbose)
        inst._mark_read(True)           # ← important
        return inst

    def read(
        self,
        lines: Sequence[str] | None = None,
    ) -> "JBlocks":
        if lines is None:
            raise ValueError("lines are required")
        blk = self.from_lines(lines, verbose=self.verbose)
        self.blocks = blk.blocks
        self._mark_read(True)
        return self

    def write(self) -> List[str]:
        out: List[str] = []
        for b in self.blocks:
            out.extend(b.write())
        return out

    # convenience exports
    def to_numpy(self) -> List[dict[str, np.ndarray]]:
        return [b.to_numpy() for b in self.blocks]

    def to_dataframe(self):
        try:
            import pandas as pd  # type: ignore
        except Exception as exc:  # pragma: no cover
            raise ImportError(
                "pandas is required for to_dataframe()"
            ) from exc
        dfs = [b.to_dataframe() for b in self.blocks]
        if not dfs:
            return pd.DataFrame()
        return pd.concat(dfs, ignore_index=True)

    def qa_summary(self) -> List[dict[str, Any]]:
        return [b.qa_summary() for b in self.blocks]

    @property
    def stations(self) -> list[str]:
        out: list[str] = []
        for b in self.blocks:
            if b.station is not None:
                out.append(b.station)
        return out
    
    @property
    def kinds(self) -> list[str]:
        out: list[str] = []
        for b in self.blocks:
            if b.kind is not None:
                out.append(b.kind)
        return out
    
    def select(
        self, *, kind: str | None = None, comp: str | None = None
    ) -> list["JBlock"]:
        out: list[JBlock] = []
        for b in self.blocks:
            if kind is not None and b.kind != kind:
                continue
            if comp is not None and b.comp != comp:
                continue
            out.append(b)
        return out
    
    def period_range(self) -> tuple[float, float] | None:
        lo, hi = None, None
        for b in self.blocks:
            a = b.to_numpy()
            if len(a["period"]) == 0:
                continue
            pmin = float(np.nanmin(a["period"]))
            pmax = float(np.nanmax(a["period"]))
            lo = pmin if lo is None else min(lo, pmin)
            hi = pmax if hi is None else max(hi, pmax)
        if lo is None or hi is None:
            return None
        return (lo, hi)
    

def _is_miss(x: float) -> bool:
    return x == MISSING_FLOAT


def _nanz(x: float) -> float:
    return np.nan if _is_miss(x) else x


def _safe_div(a: float, b: float) -> float:
    try:
        return a / b
    except Exception:
        return np.nan


def _index_of_subsequence(
    hay: Sequence[str], needle: Sequence[str]
) -> int:
    if not needle:
        return -1
    n = len(needle)
    for i in range(0, len(hay) - n + 1):
        ok = True
        for j in range(n):
            if hay[i + j].strip() != needle[j].strip():
                ok = False
                break
        if ok:
            return i
    return -1

def _block_factory(head: Head, *, verbose: int) -> "JBlock":
    k = head.dtype.kind
    if k in KIND_RHO_PHI:
        return RBlock(head=head, verbose=verbose)
    if k in KIND_COMPLEX_TF:
        return TFBlock(head=head, verbose=verbose)
    raise ValueError(f"unsupported kind: {k!r}")


def _extract_first_head_and_body(
    lines: Sequence[str],
) -> Tuple[Head, List[str]]:
    seq = list(lines)
    h = Head.from_lines(seq)
    trip = h.write()
    idx = _index_of_subsequence(seq, trip)
    if idx < 0:
        idx = 0
        while idx < len(seq) and trip[-1] != seq[idx].strip():
            idx += 1
    start = idx + len(trip)
    return h, seq[start:]



def _extract_all_blocks(
        lines: Sequence[str], *, verbose: int = 0
    ) -> List["JBlock"]:
    seq = list(lines)
    blocks: List[JBlock] = []

    while seq:
        # skip comments/info/blank
        i = 0
        while i < len(seq) and (RE_BLANK.match(seq[i]) or
                                seq[i].lstrip().startswith("#") or
                                seq[i].lstrip().startswith(">")):
            i += 1
        seq = seq[i:]
        if not seq:
            break

        # need a station
        if not RE_STATION.match(seq[0]):
            # drop one line and keep searching
            seq = seq[1:]
            continue

        # find dtype/count in either order
        j = 1
        while j < len(seq) and RE_BLANK.match(seq[j]):
            j += 1
        if j >= len(seq):
            break

        dtype_idx = count_idx = None
        try:
            parse_datatype_units(seq[j])       # station -> dtype -> count
            dtype_idx = j
            k = j + 1
            while k < len(seq) and RE_BLANK.match(seq[k]):
                k += 1
            if k < len(seq) and RE_NPOINTS.match(seq[k]):
                count_idx = k
        except JParseError:
            if RE_NPOINTS.match(seq[j]):       # station -> count -> dtype
                count_idx = j
                t = j + 1
                while t < len(seq) and not RE_STATION.match(seq[t]):
                    if not RE_BLANK.match(seq[t]):
                        try:
                            parse_datatype_units(seq[t])
                            dtype_idx = t
                            break
                        except JParseError:
                            pass
                    t += 1

        if dtype_idx is None or count_idx is None:
            # header incomplete, drop station and continue
            seq = seq[1:]
            continue

        # build head from raw header strings
        header = [seq[0], seq[dtype_idx], seq[count_idx]]
        head = Head(verbose=verbose).read(header)

        # body starts after count line; read up to n rows defensively
        body = seq[count_idx + 1:]
        blk = _block_factory(head, verbose=verbose).read((head, body))
        blocks.append(blk)

        # consume header (3) + declared rows (head.n) from seq
        rows = max(head.n or 0, 0)
        seq = body[rows:]

    return blocks


def _locate_header_indices(
    seq: Sequence[str], start: int
) -> tuple[int, int, int] | None:
    """
    Return (s_idx, d_idx, n_idx) for the next header triple starting
    at 'start', or None if none is found.

    Accepts both orders:
      station -> dtype -> count
      station -> count -> dtype (before next station/EOF)
    """
    n = len(seq)
    i = start

    # skip comments/info/blank first
    while i < n and (
        RE_COMMENT.match(seq[i]) or RE_INFO.match(seq[i]) or
        RE_BLANK.match(seq[i])
    ):
        i += 1

    # find station (must not be dtype/info/comment/blank)
    while i < n:
        if RE_STATION.match(seq[i]):
            if not RE_DATATYPE_UNITS.match(seq[i]):
                s = i
                break
        i += 1
    else:
        return None  # no station

    # first non-blank after station
    j = s + 1
    while j < n and RE_BLANK.match(seq[j]):
        j += 1
    if j >= n:
        return None

    # case A: dtype then count
    try:
        parse_datatype_units(seq[j])
        d = j
        k = d + 1
        while k < n and RE_BLANK.match(seq[k]):
            k += 1
        if k < n and RE_NPOINTS.match(seq[k]):
            return (s, d, k)
    except JParseError:
        pass

    # case B: count then dtype later (before next station/EOF)
    if RE_NPOINTS.match(seq[j]):
        k = j
        t = k + 1
        while t < n:
            if RE_STATION.match(seq[t]):
                break  # next header starts
            if RE_BLANK.match(seq[t]):
                t += 1
                continue
            try:
                parse_datatype_units(seq[t])
                d = t
                return (s, d, k)
            except JParseError:
                t += 1
                continue
        return None

    # unknown header after station; advance caller
    return None

TFRow.__doc__= r"""
One parsed transfer-function row of a J *Z/Q/C/T* block.

The row stores period (s), real/imag parts, a standard error,
a row weight, and a derived ``rej`` flag based on J-format
rules.

Parameters
----------
period : float
    Period in seconds. If input stored frequency (Hz) as a
    negative number, it is normalized to a positive period.
real, imag : float
    Real and imaginary parts of the transfer function.
error : float
    Standard error (format-specific). May be missing.
weight : float
    Row weight. Negative values mark the row as rejected.
rej : bool
    Convenience rejection flag derived from the rules above.

Notes
-----
The J specification allows several TF families (Z/Q/C/T). The
row schema stays the same; only the physical meaning differs.

Examples
--------
>>> tf = TFRow(
...     period=0.1, real=1.0, imag=-2.0,
...     error=0.1, weight=1.0, rej=False,
... )
>>> tf.period, tf.rej
(0.1, False)

See Also
--------
TFBlock : Sequence of :class:`TFRow` with I/O helpers.
RRow : Row model for resistivity/phase blocks.

References
----------
.. [1] A. G. Jones (1994). Magnetotelluric data file
   J-format, version 2.0.
"""

RRow.__doc__=r"""
One parsed resistivity/phase row of a J *R/S* block.

The row stores period (s), apparent resistivity (Ω·m), phase
(deg), their 1σ bounds, and per-column weights. A derived
boolean flag ``rej`` marks a row as rejected according to the
J-format rules.

Parameters
----------
period : float
    Period in seconds. If input stored frequency (Hz) as a
    negative number, it is normalized to a positive period.
rho : float
    Apparent resistivity. Values ``< 0`` mark the estimate as
    rejected.
pha : float
    Phase in degrees.
rhomax, rhomin : float
    Upper and lower 1σ bounds of resistivity. May be missing.
phamax, phamin : float
    Upper and lower 1σ bounds of phase. May be missing.
wrho, wpha : float
    Weights. A negative weight marks the estimate as rejected.
rej : bool
    Convenience rejection flag derived from the rules above.

Notes
-----
Missing numeric values are often written as ``-999.0`` in J
files. The higher-level APIs convert them to ``NaN`` for
numeric workflows.

Examples
--------
>>> row = RRow(
...     period=1.0, rho=100.0, pha=45.0,
...     rhomax=110.0, rhomin=90.0,
...     phamax=50.0, phamin=40.0,
...     wrho=1.0, wpha=1.0, rej=False,
... )
>>> row.rho, row.rej
(100.0, False)

See Also
--------
RBlock : Sequence of :class:`RRow` with I/O helpers.
TFRow : Row model for transfer-function blocks.

References
----------
.. [1] A. G. Jones (1994). Magnetotelluric data file
   J-format, version 2.0.
"""
