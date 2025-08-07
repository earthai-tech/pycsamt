# -*- coding: utf-8 -*-
#       Author: LKouadio <etanoyau@gmail.com>
#       License: LGPL‑3.0-or-later
"""pycsamt.zonge.base

A *single* abstract‑ish base container that all future AVG‑derived
objects will inherit from.  The class purposefully keeps the public API
minimal yet robust so that higher‑level helpers (tensors, station
objects, visualisers…) can rely on it without caring about the original
AVG flavour.

Light‑weight *data‑containers* for Zonge CSAMT/AMT AVG products.

Only the canonical, **format‑independent** fields are kept here.  The
module purposefully stays free of I/O logic – that lives in
:pyfile:`pycsamt.zonge.utils`.  These classes are designed so they can
be reused by higher‑level objects (e.g. *Station*, *Survey*, *Tensor*)
without dragging pandas around.

Notes
-----
* The class is *frozen* and uses *__slots__* – extending classes should
  declare their own extra slots.
* Input is expected to be the tidy frame returned by
  :func:`pycsamt.zonge.utils.load_avg`.
"""
from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field, asdict
from pathlib     import Path
from typing      import Dict, Any, Iterable, Sequence, Tuple, Optional
import json
from datetime import datetime 
from inspect import getmembers 

import numpy as np
import pandas as pd

from ..log.logger import get_logger 
logger = get_logger(__name__)


__all__ = [
    "FieldAliases", "AvgProperties", "BaseAVG" , "AVGComponentBase"]


@dataclass(slots=True, frozen=True)
class BaseAVG:  # noqa: D101 – documented above.
    data: pd.DataFrame
    meta: Dict[str, str] = field(default_factory=dict)
    source: Path | None   = None

    @property
    def nrows(self) -> int:  # noqa: D401 simple accessor
        return len(self.data)

    @property
    def columns(self) -> Sequence[str]:  # noqa: D401
        return tuple(self.data.columns)


    def __str__(self) -> str:  # noqa: D401 simple
        loc = f", file='{self.source.name}'" if self.source else ""
        return (f"AVGFrame[{self.nrows}×{len(self.columns)}]"
                f"{loc} cols={list(self.columns)}")

    def __repr__(self) -> str:  # noqa: D401 simple
        meta_keys = list(self.meta)[:3]
        extra     = "…" if len(self.meta) > 3 else ""
        return (f"BaseAVG(rows={self.nrows}, columns={len(self.columns)}, "
                f"meta_keys={meta_keys}{extra})")

    def to_json(self, *, orient: str = 'records', indent: int = 0) -> str:
        """Serialise the **data** frame to JSON (metadata excluded)."""
        return self.data.to_json(orient=orient, indent=indent)

    def meta_as_json(self, *, indent: int = 0) -> str:
        return json.dumps(self.meta, indent=indent)


    def __iter__(self) -> Iterable[dict]:  # iterate over rows as dicts
        for _, row in self.data.iterrows():
            yield row.to_dict()

    def asdict(self) -> Dict[str, Any]:
        """Whole object to a plain dict (including metadata)."""
        return {
            'data'  : json.loads(self.to_json()),
            'meta'  : self.meta,
            'source': str(self.source) if self.source else None,
        }

class FieldAliases:
    """Variants for frequently encountered column labels."""

    missing_values: tuple = (' ', '*', 'nan', 'NaN', None)

    # spatial
    longitude: Tuple[str, ...] = ('lon', 'longitude', 'LONG', 'LON')
    latitude:  Tuple[str, ...] = ('lat', 'latitude', 'LAT', 'LATITUDE')
    easting:   Tuple[str, ...] = ('e', 'east', 'easting', 'EASTING')
    northing:  Tuple[str, ...] = ('n', 'north', 'northing', 'NORTHING')
    elevation: Tuple[str, ...] = ('elev', 'elevation', 'ELEV', 'ELEVATION')

    # survey logistics
    station:   Tuple[str, ...] = ('sta', 'station', 'stn', 'dot', 'skp')
    freq:      Tuple[str, ...] = ('freq', 'frequency', 'Freq.')
    component: Tuple[str, ...] = ('comp', 'component')

    # field data
    amps:  Tuple[str, ...] = ('tx.amp', 'amps', 'Tx.Amp')
    emag:  Tuple[str, ...] = ('e.mag', 'emag', 'E.mag')
    ephz:  Tuple[str, ...] = ('e.phz', 'ephz', 'E.phz')
    hmag:  Tuple[str, ...] = ('h.mag', 'hmag', 'b.mag', 'B.mag')
    hphz:  Tuple[str, ...] = ('hphz', 'h.phz', 'b.phz', 'B.phz')
    rho:   Tuple[str, ...] = ('resistivity', 'rho', 'ares.mag')
    phase: Tuple[str, ...] = ('phase', 'z.phz')



@dataclass(slots=True)
class AvgProperties:
    """Minimal, *format‑agnostic* representation of an AVG row."""

    station: float | int
    freq:    float
    comp:    str
    amps:    float
    emag:    float
    ephz:    float
    hmag:    float
    hphz:    float
    rho:     float
    phase:   float

    # optional error bars 
    emag_err:   Optional[float] = None
    ephz_err:   Optional[float] = None
    hmag_err:   Optional[float] = None
    hphz_err:   Optional[float] = None
    rho_err:    Optional[float] = None
    phase_err:  Optional[float] = None

    # bookkeeping 
    source: Optional[Path] = field(default=None, repr=False)
    idx:    Optional[int]  = field(default=None, repr=False)


    def __post_init__(self):
        # ensure component code is always upper‑case for consistency
        self.comp = str(self.comp).strip().upper()

    # pretty string – concise
    def __str__(self) -> str:
        return (f"AVG(stn={self.station}, f={self.freq:g} Hz, "
                f"comp={self.comp}, rho={self.rho:.1f} Ω·m)")

    # repr – full but still readable within 62‑char width
    def __repr__(self) -> str:  # noqa: D401 simple.
        core = (
            f"station={self.station!r}, freq={self.freq!r}, "
            f"comp={self.comp!r}, emag={self.emag!r}, hmag={self.hmag!r}, "
            f"rho={self.rho!r}"
        )
        return f"AvgProperties({core})"

    def asdict(self, *, include_errors: bool = True) -> Dict[str, Any]:
        """Return a plain‑dict representation (for JSON, etc.)."""
        data = asdict(self)
        if not include_errors:
            for k in list(data.keys()):
                if k.endswith('_err'):
                    data.pop(k)
        # strip internal bookkeeping
        data.pop('source', None)
        data.pop('idx', None)
        return data

    def to_json(self, *, indent: int = 0) -> str:
        """Serialise to JSON (without path / row index)."""
        return json.dumps(self.asdict(), indent=indent)

    # allow dict‑like access – but *read‑only*
    def __getitem__(self, key: str) -> Any:  # type: ignore
        return getattr(self, key)

    def keys(self) -> Iterable[str]:  # noqa: D401 – iterator helper
        return (f.name for f in self.__dataclass_fields__.values())
    


class AVGComponentBase(ABC):
    """Abstract base for a single AVG data component."""

    __slots__ = ("_data", "_meta")

    def __init__(self,
                 data: np.ndarray | None = None,
                 meta: Dict[str, Any] | None = None) -> None:
        self._data: np.ndarray      = (
            np.asarray(data) if data is not None else np.empty(0)
        )
        self._meta: Dict[str, Any]  = meta or {}

    @abstractmethod
    def read(self, source: Sequence[str] | np.ndarray) -> None: ...

    @abstractmethod
    def write(self) -> Sequence[str]: ...

    @property
    def data(self) -> np.ndarray:                # read-only view
        return self._data

    @property
    def meta(self) -> Dict[str, Any]:
        return self._meta

    def asdict(self, *, include_meta: bool = True) -> Dict[str, Any]:
        d = {"data": self._data.tolist()}
        if include_meta:
            d["meta"] = self._meta
        return d

    # simple JSON serialiser
    def to_json(self, *, indent: int = 0) -> str:
        return json.dumps(self.asdict(), indent=indent)

    # printable representations
    def __str__(self) -> str:
        shape = "×".join(str(s) for s in self._data.shape)
        return f"{self.__class__.__name__}({shape})"

    __repr__ = __str__

    def to_frame(self) -> pd.DataFrame:
        """
        Return a *pandas* ``DataFrame`` representing the component.
        Sub-classes should override when a simple ``pd.Series`` cast
        is insufficient.
        """
        return pd.DataFrame({self.__class__.__name__.lower(): self._data})


# ------------------------------------------------------------------
# Canonical → kind-2 column label
_CANON2K2: dict[str, str] = {}

def _canon2k2(name: str) -> str:
    """
    Return the preferred *kind-2* label for a canonical ``name``.

    ``FieldAliases`` defines for each logical field a **tuple** of
    variants (first entry = Zonge’s original spelling).  Build a
    reverse lookup once at import time.
    """

    if not _CANON2K2:                                    # one-off build
        for attr, variants in getmembers(
                FieldAliases, lambda x: isinstance(x, tuple)):
            canon = attr.lower()
            _CANON2K2[canon] = variants[0]               # take the first variant
    return _CANON2K2.get(name.lower(), name)

def write_avg(
    core:     pd.DataFrame,
    extra:    pd.DataFrame | None,
    meta:     Dict[str, str] | None,
    path:     str | Path,
    *,
    stamp:    bool = True
    ) -> None:
    """
    Serialise *core* + *extra* columns back to an AVG text file
    using the modern **kind-2** comma-separated flavour.

    Parameters
    ----------
    core, extra
        DataFrames obtained from :func:`extract_core_columns`.
        *extra* can be *None* if you do not need ancillary fields.
    meta
        Header `key = value` dictionary.  Will be written at the
        top of the file (one per line, prefixed by ``$``).
    path
        Destination filename.
    stamp
        When *True* appends a ``$Written`` time-stamp.
    """
    path   = Path(path)
    header = []
    if meta:
        for k, v in meta.items():
            header.append(f"${k} = {v}")
    if stamp:
        ts = datetime.utcnow().isoformat(timespec="seconds") + "Z"
        header.append(f"$Written = {ts}")
    header.append("")                       # blank line before table

    block = pd.concat([core, extra], axis=1) if extra is not None else core
    block = block.copy()
    # restore original case for kind-2 standard columns
    # block.rename(columns={c: FieldAliases.get(c, c)  # type: ignore
    #                       for c in block.columns},
    #              inplace=True)
    
    # restore original *kind-2* labels
    block.rename(columns={col: _canon2k2(col) for col in block.columns},
                inplace=True)

    csv_txt = block.to_csv(index=False)
    path.write_text("\n".join(header) + csv_txt)
    logger.info("AVG written → %s", path)