# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from ..log.logger import get_logger
from .base import JComponentBase
from .config import (
    DTYPE_SPECS,
    ENCODING_DEFAULT,
    KIND_COMPLEX_TF,
    KIND_RHO_PHI,
    MISSING_FLOAT,
    RE_BLANK,
    RE_COMMENT,
    RE_DATATYPE_UNITS,
    RE_INFO,
    RE_NPOINTS,
    RE_STATION,
)
from .heads import Head
from .utils import (
    JParseError,
    iter_lines,
    parse_datatype_units,
)

logger = get_logger(__name__)

__all__ = [
    "RRow",
    "TFRow",
    "JBlock",
    "RBlock",
    "TFBlock",
    "JBlocks",
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
    >>> blk = JBlock.from_lines(["S01", "RXY", "1", "..."])
    >>> blk.station, blk.nrows
    ('S01', 1)

    See Also
    --------
    RBlock, TFBlock : Concrete block families.
    Heads : Header-only convenience container.

    References
    ----------
    .. [JBlock-1] A. G. Jones (1994). Magnetotelluric data file
       J-format, version 2.0.
    """

    def __init__(
        self,
        head: Head | None = None,
        *,
        rows: Sequence[Any] | None = None,
        verbose: int = 0,
    ) -> None:
        super().__init__(verbose=verbose)
        self.head: Head | None = head
        self.rows: list[Any] = list(rows or [])

    @classmethod
    def from_file(
        cls,
        j_fn: str | Path,
        *,
        verbose: int = 0,
    ) -> JBlock:
        lines = list(iter_lines(j_fn, encoding=ENCODING_DEFAULT))
        return cls.from_lines(lines, verbose=verbose)

    @classmethod
    def from_lines(
        cls,
        lines: Sequence[str],
        *,
        verbose: int = 0,
    ) -> JBlock:
        head, body = _extract_first_head_and_body(lines)
        blk = _block_factory(head, verbose=verbose)
        return blk.read((head, body))

    def read(
        self,
        data: tuple[Head, Sequence[str]] | None = None,
    ) -> JBlock:
        raise NotImplementedError("override in subclass")

    def write(self) -> list[str]:
        raise NotImplementedError("override in subclass")

    # ---- normalization & QA (common facade) ----
    def normalize(self) -> JBlock:
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
            raise ImportError("pandas is required for to_dataframe()") from exc
        data = self.to_numpy()
        return pd.DataFrame(data)

    @property
    def nrows(self) -> int:
        return len(self.rows)

    @property
    def station(self) -> str | None:
        return None if self.head is None else self.head.station

    @property
    def kind(self) -> str | None:
        return None if self.head is None else self.head.kind

    @property
    def comp(self) -> str | None:
        return None if self.head is None else self.head.comp

    @property
    def units(self) -> str | None:
        return None if self.head is None else self.head.units

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
    ...     "S01",
    ...     "RXY",
    ...     "2",
    ...     "-1.0 100 45 110 90 50 40 1 1",
    ...     " 2.0 -5  30  35 25 40 20 1 1",
    ... ]
    >>> blk = RBlock.from_lines(lines)
    >>> a = blk.to_numpy()
    >>> a["period"][0]  # 1/Hz -> s
    1.0
    >>> a["rej"][1]  # rho < 0
    True

    See Also
    --------
    TFBlock : Transfer-function block family.
    RRow : Single row record used by this block.

    References
    ----------
    .. [RBlock-1] A. G. Jones (1994). Magnetotelluric data file
       J-format, version 2.0.
    """

    def read(
        self,
        data: tuple[Head, Sequence[str]] | None = None,
    ) -> RBlock:
        if data is None:
            raise ValueError("data is required")
        head, body = data
        spec = DTYPE_SPECS[head.dtype.kind]
        rows: list[RRow] = []
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

    def normalize(self) -> RBlock:
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

    def write(self) -> list[str]:
        if not self.head or self.head.n is None:
            raise ValueError("missing head")
        out: list[str] = []
        out.extend(self.head.write())
        for r in self.rows:
            out.append(
                (
                    f"{r.period: .6e} {r.rho: .6e} {r.pha: .6e} {r.rhomax: .6e} {r.rhomin: .6e} "
                    f"{r.phamax: .6e} {r.phamin: .6e} {r.wrho: .6e} {r.wpha: .6e}"
                ).strip()
            )
        return out

    @property
    def columns(self) -> tuple[str, ...]:
        return (
            "period",
            "rho",
            "pha",
            "rhomax",
            "rhomin",
            "phamax",
            "phamin",
            "wrho",
            "wpha",
            "rej",
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
    ...     "S02",
    ...     "ZXY SI",
    ...     "2",
    ...     "-1.0e+1 1.0 -2.0 0.1 1.0",
    ...     " 5.0    -999 3.0 -999 -1.0",
    ... ]
    >>> blk = TFBlock.from_lines(lines)
    >>> a = blk.to_numpy()
    >>> round(a["period"][0], 3)
    0.1
    >>> a["rej"][1]  # weight < 0
    True

    See Also
    --------
    RBlock : Resistivity/phase block family.
    TFRow : Single row record used by this block.

    References
    ----------
    .. [TFBlock-1] A. G. Jones (1994). Magnetotelluric data file
       J-format, version 2.0.
    """

    def read(
        self,
        data: tuple[Head, Sequence[str]] | None = None,
    ) -> TFBlock:
        if data is None:
            raise ValueError("data is required")
        head, body = data
        spec = DTYPE_SPECS[head.dtype.kind]
        rows: list[TFRow] = []
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

    def normalize(self) -> TFBlock:
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
        rms = np.sqrt(np.nanmean(np.square(a["error"]))) if n else np.nan
        return {
            "kind": "TF",
            "nrows": n,
            "rejected": rej,
            "kept": keep,
            "period_min": float(pmin),
            "period_max": float(pmax),
            "error_rms": float(rms),
        }

    def write(self) -> list[str]:
        if not self.head or self.head.n is None:
            raise ValueError("missing head")
        out: list[str] = []
        out.extend(self.head.write())
        for r in self.rows:
            out.append(
                (
                    f"{r.period: .6e} {r.real: .6e} {r.imag: .6e} {r.error: .6e} {r.weight: .6e}"
                ).strip()
            )
        return out

    @property
    def columns(self) -> tuple[str, ...]:
        return (
            "period",
            "real",
            "imag",
            "error",
            "weight",
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
    >>> col = JBlocks.from_file("data/j/kb0-s001.txt")
    >>> col.n >= 1
    True
    >>> [b.station for b in col.select(kind="R")]
    ['KB0001', ...]

    See Also
    --------
    Head, Heads : Header parsers used during discovery.
    RBlock, TFBlock : Concrete block families contained here.

    References
    ----------
    .. [JBlocks-1] A. G. Jones (1994). Magnetotelluric data file
       J-format, version 2.0.
    """

    def __init__(
        self,
        blocks: Sequence[JBlock] | None = None,
        *,
        verbose: int = 0,
    ) -> None:
        super().__init__(verbose=verbose)
        self.blocks: list[JBlock] = list(blocks or [])

    @property
    def n(self) -> int:
        return len(self.blocks)

    @classmethod
    def from_file(
        cls,
        j_fn: str | Path,
        *,
        verbose: int = 0,
    ) -> JBlocks:
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
    ) -> JBlocks:
        blocks = _extract_all_blocks(lines, verbose=verbose)
        inst = cls(blocks=blocks, verbose=verbose)
        inst._mark_read(True)  # ← important
        return inst

    def read(
        self,
        lines: Sequence[str] | None = None,
    ) -> JBlocks:
        if lines is None:
            raise ValueError("lines are required")
        blk = self.from_lines(lines, verbose=self.verbose)
        self.blocks = blk.blocks
        self._mark_read(True)
        return self

    def write(self) -> list[str]:
        out: list[str] = []
        for b in self.blocks:
            out.extend(b.write())
        return out

    # convenience exports
    def to_numpy(self) -> list[dict[str, np.ndarray]]:
        return [b.to_numpy() for b in self.blocks]

    def to_dataframe(self):
        try:
            import pandas as pd  # type: ignore
        except Exception as exc:  # pragma: no cover
            raise ImportError("pandas is required for to_dataframe()") from exc
        dfs = [b.to_dataframe() for b in self.blocks]
        if not dfs:
            return pd.DataFrame()
        return pd.concat(dfs, ignore_index=True)

    def qa_summary(self) -> list[dict[str, Any]]:
        return [b.qa_summary() for b in self.blocks]

    @property
    def stations(self) -> list[str]:
        out: list[str] = []
        for b in self.blocks:
            if b.station is not None:
                out.append(b.station)
        return out

    @property
    def station(self) -> str | None:
        """
        Return the single, unique station name for
        the block collection.
        Returns None if no blocks are present.
        """
        if not self.blocks:
            return None
        # All blocks should have the same
        # station name after correct parsing
        return self.blocks[0].station

    @property
    def kinds(self) -> list[str]:
        out: list[str] = []
        for b in self.blocks:
            if b.kind is not None:
                out.append(b.kind)
        return out

    def select(
        self, *, kind: str | None = None, comp: str | None = None
    ) -> list[JBlock]:
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


def _index_of_subsequence(hay: Sequence[str], needle: Sequence[str]) -> int:
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


def _block_factory(head: Head, *, verbose: int) -> JBlock:
    k = head.dtype.kind
    if k in KIND_RHO_PHI:
        return RBlock(head=head, verbose=verbose)
    if k in KIND_COMPLEX_TF:
        return TFBlock(head=head, verbose=verbose)
    raise ValueError(f"unsupported kind: {k!r}")


def _extract_first_head_and_body(
    lines: Sequence[str],
) -> tuple[Head, list[str]]:
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
) -> list[JBlock]:
    """
    Robustly parse all data blocks from a J-file for
    a single station.

    This parser finds the single station name for the
    file and then iteratively finds all subsequent data blocks,
    which are identified by a [data_type, count] pair.
    """
    seq = list(lines)
    blocks: list[JBlock] = []

    # 1. Find the first station line in the entire file. This will be the
    #    station for ALL subsequent blocks.
    station_line = None
    first_data_line_index = 0
    for i, line in enumerate(seq):
        # Skip over the initial comment and info blocks
        if (
            RE_COMMENT.match(line)
            or RE_INFO.match(line)
            or RE_BLANK.match(line)
        ):
            continue
        # The first non-header line must be the station
        if RE_STATION.match(line):
            station_line = line
            first_data_line_index = i + 1
            break

    if station_line is None:
        # If no station is found in the file, there are no blocks to parse.
        if verbose:
            logger.warning("No station line found. Cannot parse data blocks.")
        return []

    # 2. Now, loop through the rest of the
    # file from where we found the station.
    i = first_data_line_index
    while i < len(seq):
        # Find the next data type line, skipping any blank lines
        dtype_line = None
        dtype_idx = -1
        while i < len(seq):
            line = seq[i].strip()
            if RE_BLANK.match(line):
                i += 1
                continue
            try:
                parse_datatype_units(line)
                dtype_line = line
                dtype_idx = i
                break
            except JParseError:
                # This line is not a dtype, probably
                # leftover data or garbage. Skip it.
                i += 1

        if not dtype_line:
            break  # No more data type lines found, we are done.

        # 3. Find the next count line, skipping blank lines
        count_line = None
        count_idx = -1
        j = dtype_idx + 1
        while j < len(seq):
            line = seq[j].strip()
            if RE_BLANK.match(line):
                j += 1
                continue
            if RE_NPOINTS.match(line):
                count_line = line
                count_idx = j
                break
            # If we find something that is not
            # a count, the header is malformed.
            break

        if not count_line:
            # Could not find a count for the detected
            # data type, stop parsing.
            if verbose:
                logger.warning(
                    "Found data type '{dtype_line}'"
                    " but no following count. "
                    "Stopping parse."
                )
            break

        # 4. We have the full header. Create the block.
        header_lines = [station_line, dtype_line, count_line]
        head = Head(verbose=verbose).read(header_lines)

        body_start_index = count_idx + 1
        body_lines = seq[body_start_index:]

        # Create the block and let its `read` method consume the data rows
        blk = _block_factory(head, verbose=verbose).read((head, body_lines))
        blocks.append(blk)

        # 5. Advance the main loop index past the data we just consumed.
        rows_consumed = max(head.n or 0, 0)
        i = body_start_index + rows_consumed

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
        RE_COMMENT.match(seq[i])
        or RE_INFO.match(seq[i])
        or RE_BLANK.match(seq[i])
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


TFRow.__doc__ = r"""
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
.. [TFRow-1] A. G. Jones (1994). Magnetotelluric data file
   J-format, version 2.0.
"""

RRow.__doc__ = r"""
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
.. [RRow-1] A. G. Jones (1994). Magnetotelluric data file
   J-format, version 2.0.
"""
