# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

"""Base classes for SEG-EDI objects and components.

This module provides a small reusable base class that standardizes:
string representation, pretty printing of key/value lines in the
SEG-EDI style, and simple (de)serialization helpers.

Subclasses only need to:
- set ``_section`` (e.g., ``"HEAD"``, ``"INFO"``, ...)
- define typed attributes (via ``__annotations__`` or dataclass)
- optionally override ``validate`` or value formatting if needed.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from collections.abc import (
    Iterable,
    Iterator,
    Mapping,
    Sequence,
)
from dataclasses import fields as dc_fields
from dataclasses import is_dataclass
from pathlib import Path
from typing import (
    Any,
    TypeVar,
    Union,
)

import numpy as np

from ..log.logger import get_logger as _get_logger
from .validation import IsEdi

PathLike = Union[str, Path]


__all__ = [
    "Base",
    "BaseMixin",
    "SEGBase",
    "EDIComponentBase",
    "EdiFileBase",
]

_T = TypeVar("_T")


# Common, light mixin for facade
class Base:  # picked by seg.config (mixin discovery)
    INDENT: int = 2
    PER_LINE: int = 6
    FLOAT_FMT: str = "{: .6E}"

    def __init__(
        self,
        *args,
        verbose: int | bool = 0,
        logger=None,
        **kwargs,
    ) -> None:
        """
        Cooperative init: set verbosity and instance logger,
        then chain to next MRO parent.
        """
        try:
            super().__init__(*args, **kwargs)
        except Exception:
            # tolerate bases without cooperative __init__
            pass

        self.verbose: int = int(verbose)
        name = f"{self.__class__.__module__}.{self.__class__.__name__}"
        self._logger = logger if logger is not None else self._logger_factory(name)

    @staticmethod
    def _logger_factory(name: str):
        """Return a module/project logger (null-safe)."""
        if _get_logger is not None:
            return _get_logger(name)
        import logging as _logging

        lg = _logging.getLogger(name)
        if not lg.handlers:
            lg.addHandler(_logging.NullHandler())
        return lg


# Back-compat / alternate mixin names expected by config
class BaseMixin(Base):  # pragma: no cover
    pass


class SEGBase(Base):  # pragma: no cover
    pass


class EDIComponentBase(Base):
    _section: str | None = None
    _float_fmt: str = "{: .6E}"
    _indent: int = 2
    _show_in_str: Sequence[str] | None = None
    _repr_max_items: int = 6
    _repr_max_chars: int = 120

    def __init__(
        self,
        *args,
        verbose: int | bool = 0,
        logger=None,
        **kwargs,
    ) -> None:
        """
        Ensure all components get .verbose and ._logger.
        """
        super().__init__(
            *args,
            verbose=verbose,
            logger=logger,
            **kwargs,
        )

    # -------- repr/str ---------
    def __repr__(self) -> str:
        cls = self.__class__.__name__
        parts: list[str] = []
        chars = 0
        for i, (k, v) in enumerate(self._iter_kv(True)):
            if i >= self._repr_max_items:
                parts.append("…")
                break
            frag = f"{k}={self._repr_value(v)}"
            chars += len(frag)
            if chars > self._repr_max_chars:
                parts.append("…")
                break
            parts.append(frag)
        inside = ", ".join(parts)
        return f"{cls}({inside})"

    def __str__(self) -> str:
        if self._section:
            lines = [f">{self._section}"]
            lines.extend(self.to_lines(only=self._show_in_str))
            return "\n".join(lines)
        kv = " ".join(self.to_lines(only=self._show_in_str))
        return f"{self.__class__.__name__} {kv}"

    # public
    def to_dict(self, *, exclude_none: bool = True) -> dict[str, Any]:
        out: dict[str, Any] = {}
        for k, v in self._iter_kv(exclude_none):
            out[k] = v
        return out

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> EDIComponentBase:
        obj = cls()
        for k, v in (data or {}).items():
            try:
                setattr(obj, k, v)
            except Exception:  # pragma: no cover
                pass
        return obj

    def to_lines(
        self,
        *,
        only: Sequence[str] | None = None,
        indent: int | None = None,
    ) -> list[str]:
        pad = indent if indent is not None else self._indent
        kv = list(self._iter_kv(True))
        if only is not None:
            order = list(only)
            kv = [(k, v) for (k, v) in kv if k in order]
            kv.sort(key=lambda it: order.index(it[0]))
        return [self._format_kv(k, v, pad) for (k, v) in kv]

    def to_text(self) -> str:
        return str(self)

    def update(self, /, **kw: Any) -> EDIComponentBase:
        for k, v in kw.items():
            setattr(self, k, v)
        self.validate()
        return self

    def clone(self, /, **ov: Any) -> EDIComponentBase:
        new = self.__class__()
        for k, v in self._iter_kv(False):
            setattr(new, k, v)
        for k, v in ov.items():
            setattr(new, k, v)
        new.validate()
        return new

    # subclass hook
    def validate(self) -> None:
        return None

    #  internals
    @classmethod
    def _field_names(cls) -> list[str]:
        if is_dataclass(cls):
            return [f.name for f in dc_fields(cls)]
        ann = getattr(cls, "__annotations__", None) or {}
        if ann:
            return list(ann.keys())
        return []

    def _iter_kv(self, exclude_none: bool) -> Iterator[tuple[str, Any]]:
        seen = set()
        for name in self._field_names():
            if not hasattr(self, name):
                continue
            val = getattr(self, name)
            if exclude_none and val is None:
                continue
            seen.add(name)
            yield name, val
        for name, val in self.__dict__.items():
            if name.startswith("_") or name in seen:
                continue
            if exclude_none and val is None:
                continue
            yield name, val

    def _format_kv(self, key: str, value: Any, indent: int) -> str:
        return " " * indent + f"{key}={self._format_value(value)}"

    def _format_value(self, value: Any) -> str:
        if isinstance(value, str):
            v = value.strip()
            if (any(ch.isspace() for ch in v)) or ('"' in v):
                return '"' + v.replace('"', "") + '"'
            return v
        if isinstance(value, bool):
            return "1" if value else "0"
        if isinstance(value, float):
            if value != value:  # NaN
                return "NaN"
            return self._float_fmt.format(value)
        if isinstance(value, int):
            return str(value)
        if isinstance(value, (list, tuple)):
            return " ".join(self._format_value(v) for v in value)
        return str(value)

    @staticmethod
    def _repr_value(value: Any) -> str:
        if isinstance(value, str):
            v = value.replace('"', "")
            return f'"{v}"' if (len(v) > 20 or " " in v) else v
        if isinstance(value, float):
            return f"{value:.6g}"
        if isinstance(value, (list, tuple)):
            prev = ", ".join(
                (f"{x:.6g}" if isinstance(x, float) else repr(x)) for x in value[:4]
            )
            if len(value) > 4:
                prev += ", …"
            return f"[{prev}]"
        return repr(value)


class EdiFileBase(EDIComponentBase, ABC):
    path: Path | None = None
    strict_validate: bool = True
    block_size: int = 6
    number_fmt: str = " 15.6e"
    data_header_tpl: str = ">!****{title}****!\n"

    sections: dict[str, EDIComponentBase] = {}

    _log = Base._logger_factory(__name__)

    def __init__(
        self,
        *,
        path: PathLike | None = None,
        strict_validate: bool = True,
    ) -> None:
        self.path = Path(path) if path is not None else None
        self.strict_validate = bool(strict_validate)
        self.block_size = 6
        self.number_fmt = " 15.6e"
        self.data_header_tpl = ">!****{title}****!\n"
        self.sections = {}
        if self.path is not None and self.strict_validate:
            self._validate_path(self.path)

    # ---------- I/O ------------
    @classmethod
    def from_file(cls, file: PathLike, **kw: Any) -> EdiFileBase:
        inst = cls(path=Path(file), **kw)
        inst.read()
        return inst

    def write_file(
        self,
        file: PathLike | None = None,
        *,
        overwrite: bool = True,
    ) -> Path:
        target = Path(file) if file is not None else self._default_out()
        if target.exists() and not overwrite:
            raise FileExistsError(f"{target} exists and overwrite=False.")
        text = self.compose()
        if isinstance(text, (list, tuple)):
            text = "".join(text)
        target.parent.mkdir(parents=True, exist_ok=True)
        with target.open("w", encoding="utf-8", newline="\n") as f:
            f.write(text)
        self._log.info("Wrote EDI to %s", str(target))
        return target

    @abstractmethod
    def read(self) -> EdiFileBase:
        ...

    @abstractmethod
    def compose(self) -> str | list[str]:
        ...

    # ------- registry ----------
    def add_section(self, name: str, s: EDIComponentBase) -> None:
        key = self._normalize_section_name(name)
        self.sections[key] = s

    def get_section(self, name: str) -> EDIComponentBase | None:
        return self.sections.get(self._normalize_section_name(name))

    def remove_section(self, name: str) -> None:
        self.sections.pop(self._normalize_section_name(name), None)

    def iter_sections(
        self,
    ) -> Iterable[tuple[str, EDIComponentBase]]:
        return self.sections.items()

    # ------- proxies -----------
    @property
    def dataid(self) -> str | None:
        head = self.get_section("head")
        return getattr(head, "dataid", None) if head else None

    @property
    def lat(self) -> float | None:
        head = self.get_section("head")
        dm = self.get_section("definemeasurement")
        if head and getattr(head, "lat", None) is not None:
            return getattr(head, "lat", None)
        return getattr(dm, "reflat", None) if dm else None

    @property
    def lon(self) -> float | None:
        head = self.get_section("head")
        dm = self.get_section("definemeasurement")
        if head and getattr(head, "long", None) is not None:
            return getattr(head, "long", None)
        return getattr(dm, "reflong", None) if dm else None

    @property
    def elev(self) -> float | None:
        head = self.get_section("head")
        dm = self.get_section("definemeasurement")
        if head and getattr(head, "elev", None) is not None:
            return getattr(head, "elev", None)
        return getattr(dm, "refelev", None) if dm else None

    # ------- validation --------
    def _validate_path(self, file: Path) -> None:
        IsEdi._assert_edi(str(file), deep=self.strict_validate)

    # ------- defaults ----------
    def _default_out(self) -> Path:
        stem = self.dataid or "output"
        return Path(f"{stem}.edi")

    # ------- fmt helpers -------
    @staticmethod
    def normalize_section_title(title: str) -> str:
        t = (title or "").strip()
        return t.lstrip(">").upper()

    def format_block_header(self, title: str) -> str:
        return self.data_header_tpl.format(title=self.normalize_section_title(title))

    @staticmethod
    def format_section_head(head: str) -> str:
        h = head.strip()
        if not h.startswith(">"):
            h = ">" + h
        return h.upper() + ("\n" if not h.endswith("\n") else "")

    def format_kv(
        self,
        key: str,
        value: Any,
        *,
        quote: bool = False,
        width: int | None = None,
    ) -> str:
        k = key.strip().upper()
        if value is None:
            val = ""
        else:
            if quote and isinstance(value, str):
                val = f'"{value}"'
            else:
                val = str(value)
        if width is not None:
            val = f"{val:>{width}}"
        return f"  {k}={val}\n"

    def format_data_block(
        self,
        title: str,
        data: Iterable[float | int],
        *,
        per_line: int | None = None,
        header_comment: bool = True,
        count_comment: bool = True,
    ) -> list[str]:
        t = self.normalize_section_title(title)
        vals = list(data or [])
        n = len(vals)
        per = per_line or self.block_size
        out: list[str] = []
        if header_comment:
            out.append(self.format_block_header(t))
        if count_comment:
            out.append(self.format_section_head(f">{t}  //{n}"))
        else:
            out.append(self.format_section_head(f">{t}"))
        if n == 0:
            return out
        # lazy import to avoid cycles
        from .utils import _format_block_numbers as fmtnums

        out.append(fmtnums(vals, per_line=per, indent=2) + "\n")
        return out

    # ------- misc --------------
    @staticmethod
    def ensure_descending_frequency(
        freq: Iterable[float],
    ) -> tuple[list[float], list[int] | None]:
        f = list(freq or [])
        if not f:
            return f, None
        if len(f) < 2 or f[0] >= f[-1]:
            return f, None
        idx = list(range(len(f)))[::-1]
        return f[::-1], idx

    def to_dict(self) -> dict[str, Any]:
        return {
            "path": str(self.path) if self.path is not None else None,
            "strict_validate": self.strict_validate,
            "block_size": self.block_size,
            "number_fmt": self.number_fmt,
            "data_header_tpl": self.data_header_tpl,
            "sections": {
                k: (v.to_dict() if hasattr(v, "to_dict") else dict(v.__dict__))
                for k, v in self.sections.items()
            },
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> EdiFileBase:
        inst = cls(
            path=Path(data["path"]) if data.get("path") else None,
            strict_validate=bool(data.get("strict_validate", True)),
        )
        inst.block_size = int(data.get("block_size", inst.block_size))
        inst.number_fmt = str(data.get("number_fmt", inst.number_fmt))
        inst.data_header_tpl = str(data.get("data_header_tpl", inst.data_header_tpl))
        for k, v in (data.get("sections") or {}).items():
            comp = (
                EDIComponentBase.from_dict(v)
                if isinstance(v, Mapping)
                else EDIComponentBase()
            )
            inst.add_section(k, comp)
        return inst

    def summary(self) -> str:
        parts = [
            f"file={self.path.name if self.path else '<memory>'}",
            f"dataid={self.dataid or '-'}",
        ]
        loc = []
        if self.lat is not None:
            loc.append(f"lat={self.lat:.6f}")
        if self.lon is not None:
            loc.append(f"lon={self.lon:.6f}")
        if self.elev is not None:
            loc.append(f"elev={self.elev}")
        if loc:
            parts.append(" ".join(loc))
        parts.append(f"sections={list(self.sections.keys()) or '[]'}")
        return "Edi(" + ", ".join(parts) + ")"

    def __str__(self) -> str:  # pragma: no cover
        return self.summary()

    def __repr__(self) -> str:  # pragma: no cover
        cls = self.__class__.__name__
        return f"<{cls} path={self.path!r} " f"sections={list(self.sections.keys())!r}>"

    @staticmethod
    def _normalize_section_name(name: str) -> str:
        return (name or "").strip().lstrip(">=").lower()

    def compose_to_string(self) -> str:
        text = self.compose()
        return "".join(text) if isinstance(text, list) else str(text)


class SurveyBase:
    def __init__(self, *, verbose: int = 0) -> None:
        self.verbose = int(verbose)

    @staticmethod
    def _as_floats(v: Any) -> np.ndarray:
        try:
            a = np.asarray(v, float)
        except Exception:
            a = np.asarray([], float)
        return a

    @staticmethod
    def _rng(a: np.ndarray) -> tuple[float, float] | None:
        if a.size == 0 or not np.isfinite(a).any():
            return None
        return (float(np.nanmin(a)), float(np.nanmax(a)))

    @staticmethod
    def _fmt_rng(r: tuple[float, float] | None) -> str:
        if r is None:
            return "-"
        lo, hi = r
        return f"{lo:.3f}..{hi:.3f}"

    @staticmethod
    def _fmt_list(xs: Iterable[str], maxn: int = 6) -> str:
        vals = [str(x) for x in xs]
        if len(vals) <= maxn:
            return ", ".join(vals) if vals else "-"
        head = ", ".join(vals[: maxn - 1])
        return f"{head}, +{len(vals) - (maxn - 1)}"

    def _stations(self) -> list[str]:
        # try dedicated API, else rows, else empty
        if hasattr(self, "stations"):
            try:
                ss = self.stations
                if callable(ss):
                    return [str(s) for s in ss()]
                return [str(s) for s in ss]  # property
            except Exception:
                pass
        if hasattr(self, "as_table"):
            try:
                rows = self.as_table()
                return [str(r.get("station", "")) for r in rows]
            except Exception:
                pass
        return []

    def _lat(self) -> np.ndarray:
        if hasattr(self, "lat"):
            return self._as_floats(self.lat)
        return np.asarray([], float)

    def _lon(self) -> np.ndarray:
        if hasattr(self, "lon"):
            return self._as_floats(self.lon)
        return np.asarray([], float)

    def _elev(self) -> np.ndarray:
        if hasattr(self, "elev"):
            return self._as_floats(self.elev)
        if hasattr(self, "elevation"):
            return self._as_floats(self.elevation)
        return np.asarray([], float)

    def _distance(self) -> np.ndarray:
        if hasattr(self, "distance"):
            return self._as_floats(self.distance)
        return np.asarray([], float)

    def _azimuth(self) -> float | None:
        val = None
        if hasattr(self, "azimuth"):
            try:
                v = self.azimuth
                val = float(v)  # property or value
            except Exception:
                pass
        elif hasattr(self, "get_bearing"):
            try:
                v = self.get_bearing()
                val = None if v is None else float(v)
            except Exception:
                pass
        return val

    def _step(self) -> float | None:
        # prefer explicit method if it returns scalar
        if hasattr(self, "get_step"):
            try:
                v = self.get_step()
                a = self._as_floats(v)
                if a.ndim == 0:
                    return float(a)  # scalar
                if a.size > 1:
                    d = np.diff(a)
                    return float(np.median(d))
            except Exception:
                pass
        # fallback from distance
        d = self._distance()
        if d.size > 1:
            return float(np.median(np.diff(d)))
        return None

    def _n(self) -> int:
        # length priority: distance, lat, stations, len(self)
        d = self._distance().size
        if d:
            return int(d)
        la = self._lat().size
        if la:
            return int(la)
        st = len(self._stations())
        if st:
            return int(st)
        try:
            return int(len(self))  # type: ignore[arg-type]
        except Exception:
            return 0

    def summary_dict(self) -> dict[str, Any]:
        lat = self._lat()
        lon = self._lon()
        el = self._elev()
        di = self._distance()
        az = self._azimuth()
        st = self._stations()
        return {
            "n": self._n(),
            "stations": st,
            "lat_rng": self._rng(lat),
            "lon_rng": self._rng(lon),
            "elev_rng": self._rng(el),
            "dist_rng": self._rng(di),
            "azimuth": az,
            "step": self._step(),
        }

    def __repr__(self) -> str:  # pragma: no cover
        s = self.summary_dict()
        az = "-" if s["azimuth"] is None else f"{s['azimuth']:.1f}"
        st = self._fmt_list(s["stations"], maxn=4)
        dist = self._fmt_rng(s["dist_rng"])
        return (
            f"{self.__class__.__name__}("
            f"n={s['n']}, az={az}, dist={dist}, "
            f"stations=[{st}]"
            ")"
        )

    def __str__(self) -> str:  # pragma: no cover
        s = self.summary_dict()
        lines: list[str] = []
        lines.append(self.__class__.__name__)
        lines.append(f"  count   : {s['n']}")
        lines.append(f"  stations: {self._fmt_list(s['stations'])}")
        lines.append(f"  lat     : {self._fmt_rng(s['lat_rng'])}")
        lines.append(f"  lon     : {self._fmt_rng(s['lon_rng'])}")
        lines.append(f"  elev    : {self._fmt_rng(s['elev_rng'])}")
        lines.append(f"  dist    : {self._fmt_rng(s['dist_rng'])}")
        az = "-" if s["azimuth"] is None else f"{s['azimuth']:.2f}"
        st = "-" if s["step"] is None else f"{s['step']:.2f}"
        lines.append(f"  azimuth : {az} deg")
        lines.append(f"  step    : {st} m")
        return "\n".join(lines)

    def format_table(
        self,
        rows: list[dict[str, Any]] | None = None,
        *,
        cols: list[str] | None = None,
        max_rows: int | None = 24,
    ) -> str:
        if rows is None:
            if hasattr(self, "as_table"):
                try:
                    rows = self.as_table()
                except Exception:
                    rows = []
            else:
                rows = []
        if not rows:
            return "<empty>"
        if cols is None:
            cols = list(rows[0].keys())
        if max_rows is not None and len(rows) > max_rows:
            rows = rows[:max_rows]
        # widths
        w: dict[str, int] = {}
        for c in cols:
            w[c] = max(len(c), *(len(str(r.get(c, ""))) for r in rows))
        # header
        head = " | ".join(c.ljust(w[c]) for c in cols)
        sep = "-+-".join("-" * w[c] for c in cols)
        out = [head, sep]
        # body
        for r in rows:
            out.append(" | ".join(str(r.get(c, "")).ljust(w[c]) for c in cols))
        return "\n".join(out)
