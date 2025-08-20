# -*- coding: utf-8 -*-
#       Author: LKouadio <etanoyau@gmail.com>
#       License: LGPL-3.0-or-later
"""
pycsamt.zonge.base
==================

Foundational building blocks for the *Zonge* sub-package.

Goals
-----
- Provide a **single, flexible base** every component can
  inherit from (e.g., Station, Frequency, Resistivity…).
- Offer a robust **AVGFrame** holder for data + metadata.
- Define a clear **component contract**:
    • ``from_avg(...)``  → build from an AVG file/table
    • ``read(...)``       → parse & populate attributes
    • ``write()``         → reconstruct a textual AVG block
- Be agnostic to legacy/modern files.  Upstream loaders
  (e.g., ``utils.load_avg``) normalize columns; components
  consume the tidy frame.

Design notes
------------
- All public objects keep string widths short for readability.
- Columns are lower-case canonical names (e.g., ``'freq'``,
  ``'rho'``, ``'phase'``, ``'comp'``…).
- Components target *subsets* of the frame and may expose
  additional, derived attributes internally.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field, asdict
from datetime import datetime, timezone
from pathlib import Path
from typing import (
    Any,
    Dict,
    # Iterable,
    Mapping,
    MutableMapping,
    Optional,
    Sequence,
    Tuple,
    Union,
)

import numpy as np
import pandas as pd

from ..exceptions import AvgDataError 
from ..log.logger import get_logger
from ._transfer import LegacyAVGBase
from .utils import load_avg, classify_avg_format 

logger = get_logger(__name__)

__all__ = [
    "FieldAliases", "AvgRow", "AVGFrame" , "AVGComponentBase", 
    ]

@dataclass(slots=True)
class AVGFrame:
    """
    Tidy AVG table + side metadata.

    ``data``
        A *pandas* DataFrame with **canonical**, lower-case
        column names.  Examples: ``station``, ``freq``, ``comp``,
        ``emag``, ``hmag``, ``rho``, ``phase``…
    ``meta``
        Free-form metadata (e.g., ``$Survey.Type``, units…).
    ``source``
        Optional origin path (useful for provenance/logging).
    """

    data: pd.DataFrame
    meta: Dict[str, Any] = field(default_factory=dict)
    source: Optional[Path] = None

    @property
    def nrows(self) -> int:
        """Row count."""
        return int(len(self.data))

    @property
    def columns(self) -> Tuple[str, ...]:
        """Column labels (tuple for immutability)."""
        return tuple(map(str, self.data.columns))

    def copy(self) -> "AVGFrame":
        """Deep copy the frame (safe for mutation)."""
        return AVGFrame(
            data=self.data.copy(deep=True),
            meta=dict(self.meta),
            source=self.source,
        )

    def to_json(self, *, orient: str = "records",
                indent: int = 0) -> str:
        """Serialise *data* to JSON (metadata excluded)."""
        return self.data.to_json(orient=orient, indent=indent)

    def meta_as_json(self, *, indent: int = 0) -> str:
        """Serialise metadata to JSON."""
        return pd.io.json.dumps(self.meta, indent=indent)

    def asdict(self) -> Dict[str, Any]:
        """Plain dict for diagnostics / logging."""
        return {
            "data": pd.io.json.loads(self.to_json()),
            "meta": dict(self.meta),
            "source": str(self.source) if self.source else None,
        }

    def __str__(self) -> str:
        src = f", file='{self.source.name}'" if self.source else ""
        cols = ", ".join(self.columns[:6])
        tail = "…" if len(self.columns) > 6 else ""
        return (
            f"AVGFrame[{self.nrows}×{len(self.columns)}]"
            f"{src} cols=[{cols}{tail}]"
        )

    def __repr__(self) -> str:
        keys = list(self.meta)[:4]
        more = "…" if len(self.meta) > 4 else ""
        return (
            f"AVGFrame(nrows={self.nrows}, "
            f"ncols={len(self.columns)}, "
            f"meta_keys={keys}{more})"
        )

class FieldAliases:
    """
    Variants for frequently encountered labels.  Components may use
    these tuples to search/validate presence in foreign frames.

    Note: After parsing with ``load_avg``, your frame *should* be
    canonical already; these aliases mostly serve for validation or
    when consuming user-provided, partially-normalised tables.
    """

    # Spatial / survey
    station: Tuple[str, ...] = ("station", "stn")
    freq: Tuple[str, ...] = ("freq", "frequency")
    comp: Tuple[str, ...] = ("comp",)

    # Field magnitudes / phases
    amps: Tuple[str, ...] = ("amps", "tx.amp")
    emag: Tuple[str, ...] = ("emag", "e.mag")
    ephz: Tuple[str, ...] = ("ephz", "e.phz")
    hmag: Tuple[str, ...] = ("hmag", "b.mag")
    hphz: Tuple[str, ...] = ("hphz", "b.phz")

    # Impedance / rho / phase
    zmag: Tuple[str, ...] = ("zmag",)
    zabs: Tuple[str, ...] = ("zabs", "|z|")
    rho: Tuple[str, ...] = ("rho", "ares.mag")
    phase: Tuple[str, ...] = ("phase", "z.phz")

    # Quality / weights
    mwgt: Tuple[str, ...] = ("z.mwgt",)
    pwgt: Tuple[str, ...] = ("z.pwgt",)
    e_wgt: Tuple[str, ...] = ("e.wgt",)
    h_wgt: Tuple[str, ...] = ("h.wgt",)

    # Legacy skip flag
    skp: Tuple[str, ...] = ("skp",)

@dataclass(slots=True)
class AvgRow:
    """
    Minimal, format-agnostic representation of a single AVG row.
    Useful for tests, JSON export, or component round-trips.
    """

    station: Union[int, float]
    freq: float
    comp: str
    amps: Optional[float] = None
    emag: Optional[float] = None
    ephz: Optional[float] = None
    hmag: Optional[float] = None
    hphz: Optional[float] = None
    rho: Optional[float] = None
    phase: Optional[float] = None

    # Optional quality fields (kept simple)
    e_err: Optional[float] = None   # relative %
    e_perr: Optional[float] = None  # mrad
    h_err: Optional[float] = None
    h_perr: Optional[float] = None
    rho_err: Optional[float] = None
    z_perr: Optional[float] = None

    def __post_init__(self) -> None:
        self.comp = str(self.comp).strip() or "ExHy"

    def asdict(self) -> Dict[str, Any]:
        d = asdict(self)
        return d

    def __str__(self) -> str:
        r_fmt = (
            f"{self.rho:.1f}" if isinstance(self.rho, (int, float))
            and np.isfinite(self.rho) else "nan"
        )
        return (
            f"AvgRow(stn={self.station}, f={self.freq:g} Hz, "
            f"comp={self.comp}, rho={r_fmt})"
        )

    __repr__ = __str__


class AVGComponentBase(ABC):
    """
    Abstract base-class for a single AVG component.

    Sub-classes *must* implement:
      • ``read(df, meta)``  → mutate internal state from a tidy
        table and free-form metadata.
      • ``write()``         → return a sequence of **text** lines
        representing this component in an AVG-like block (e.g.,
        a header + CSV rows).  The base provides helpers.
      • ``from_avg(...)``   → convenience alternate constructor.

    Typical usage in a component:

    >>> class Station(AVGComponentBase):
    ...     required = {"station"}
    ...     provides = {"station"}
    ...     def read(self, df, meta):
    ...         self._frame = df[["station"]].dropna()
    ...         self._meta = dict(meta)
    ...     def write(self):
    ...         return self._write_csv_block(
    ...             cols=["station"], title="$Station Block"
    ...         )

    The base keeps payloads in a protected ``_frame`` (DataFrame)
    plus a free-form ``_meta`` dict.  A small JSON/str interface
    is provided for diagnostics.
    """

    # Sets listing which columns must be present / will be produced.
    required: set[str] = set()
    provides: set[str] = set()

    def __init__(
        self,
        data: Optional[pd.DataFrame] = None,
        meta: Optional[MutableMapping[str, Any]] = None,
        *,
        name: Optional[str] = None,
        verbose =0,  
    ) -> None:
        self._frame: pd.DataFrame = (
            data.copy(deep=True) if data is not None
            else pd.DataFrame()
        )
        self._meta: Dict[str, Any] = dict(meta or {})
        self._name: str = name or self.__class__.__name__
        self.verbose = verbose 
        self._logger = get_logger(self._name)

    @classmethod
    def from_avg(
        cls,
        avg: Union[
            str, Path, AVGFrame, pd.DataFrame,
            Tuple[pd.DataFrame, Mapping[str, Any]]
        ],
        *,
        meta: Optional[Mapping[str, Any]] = None,
        **kws
    ) -> "AVGComponentBase":
        """
        Build a component from a path / AVGFrame / dataframe.

        The method accepts:
          • path → uses ``utils.load_avg`` (modern or legacy)
          • ``AVGFrame`` → uses its data/meta directly
          • ``(df, meta)`` tuple
          • bare ``df`` + explicit ``meta`` kwarg
        """
        df: pd.DataFrame
        m: Mapping[str, Any]

        if isinstance(avg, (str, Path)):
            # When loading from a file, classify it first
            path = Path(avg)
            lines = path.read_text(errors="replace").splitlines()
            kind = classify_avg_format(lines)
            df, m = load_avg(path)

            # If it's legacy, transform it to modern structure
            if kind == 1:
                transformer = LegacyAVGBase()
                ds = transformer.from_dataframe(df, meta=m)
                df = ds.to_dataframe().reset_index()
                m = ds.attrs

            frame = AVGFrame(df, dict(m), path)

        elif isinstance(avg, AVGFrame):
            frame = avg
        elif isinstance(avg, tuple) and len(avg) == 2:
            df, m = avg
            frame = AVGFrame(df, dict(m))
        elif isinstance(avg, pd.DataFrame):
            frame = AVGFrame(avg, dict(meta or {}))
        else:
            raise TypeError(
                "from_avg expects Path|AVGFrame|DataFrame|"
                "(DataFrame, meta) tuple."
            )

        obj = cls()
        try:
            obj.read(frame.data, frame.meta)
        except TypeError:
            # Fallback for older components
            obj.read(frame.data)

        return obj
    
    @abstractmethod
    def read(
        self,
        source: pd.DataFrame,
        meta: Mapping[str, Any] | None = None,
    ) -> None:
        """
        Parse **source** and mutate component state.

        Implementations should:
          1) validate required columns using ``_require(...)``
          2) copy/select needed columns into ``self._frame``
          3) set/merge metadata into ``self._meta``
        """
        ...

    @abstractmethod
    def write(self) -> Sequence[str]:
        """
        Serialise the component to **text** lines.  Implementations
        may delegate to ``_write_csv_block`` for a consistent block
        format (header + CSV).  The base does not write to disk.
        """
        ...

    @property
    def frame(self) -> pd.DataFrame:
        """Read-only view of the component table."""
        return self._frame

    @property
    def meta(self) -> Dict[str, Any]:
        """Free-form metadata."""
        return self._meta

    @property
    def name(self) -> str:
        """Short human-readable component name."""
        return self._name

    @property
    def shape(self) -> Tuple[int, int]:
        """Table shape (rows, cols)."""
        return (int(len(self._frame)), int(self._frame.shape[1] or 0))

    def asdict(self, *, include_meta: bool = True) -> Dict[str, Any]:
        """Plain dict with data (+ meta optionally)."""
        d: Dict[str, Any] = {"data": self._frame.to_dict("list")}
        if include_meta:
            d["meta"] = dict(self._meta)
        return d

    def to_json(self, *, indent: int = 0) -> str:
        """JSON serialiser for diagnostics."""
        return pd.io.json.dumps(self.asdict(), indent=indent)

    def _require(self, *cols: str) -> None:
        """
        Ensure columns are present; raise a clear error otherwise.
        """
        missing = [c for c in cols if c not in self._frame.columns]
        if missing:
            raise AvgDataError(
                f"{self._name}: missing columns {missing}"
            )

    def _select(self, *cols: str) -> pd.DataFrame:
        """
        Return a *copy* with only selected columns (silently drops
        absent ones).  Useful when composing write-blocks.
        """
        keep = [c for c in cols if c in self._frame.columns]
        return self._frame.loc[:, keep].copy()


    def _write_csv_block(
        self,
        *,
        cols: Sequence[str],
        title: Optional[str] = None,
        float_fmt: str = "%.6g",
        na_rep: str = "",
        include_meta: bool = True,
        stamp: bool = True,
    ) -> Sequence[str]:
        """
        Build a *text* block with an optional title/meta preamble
        and a CSV table of the selected ``cols``.

        This does **not** write to disk; callers own the I/O.
        """
        lines: list[str] = []

        # Title / banner (optional but nice for humans)
        if title:
            lines.append(f"\\ {title}")

        # $key=value meta lines (optional)
        if include_meta and self._meta:
            for k, v in self._meta.items():
                # Avoid multi-line explosions in headers
                v_str = str(v).replace("\n", " ").strip()
                lines.append(f"${k}={v_str}")

        # UTC stamp for provenance if desired
        if stamp:
            ts = datetime.now(timezone.utc).strftime(
                "%Y-%m-%dT%H:%M:%SZ"
            )
            lines.append(f"$Written={ts}")

        # Guard against empty data
        if self._frame.empty:
            lines.append("")  # keep a blank for clean joins
            return lines

        # Emit CSV header + rows
        table = self._select(*cols)
        # Ensure stable column order as requested
        table = table.loc[:, list(cols)]
        csv = table.to_csv(
            index=False, float_format=float_fmt, na_rep=na_rep
        )
        # Separate header from CSV with a blank line for clarity
        lines.append("")
        lines.extend(csv.splitlines())
        return lines

    def __str__(self) -> str:
        r, c = self.shape
        cols = ", ".join(self._frame.columns[:6])
        tail = "…" if self._frame.shape[1] > 6 else ""
        return (
            f"{self._name}[{r}×{c}] "
            f"cols=[{cols}{tail}]"
        )

    __repr__ = __str__



from ._transfer import LegacyAVGBase # noqa 

__all__.extend(["LegacyAVGBase"])
