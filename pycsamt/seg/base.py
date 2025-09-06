# -*- coding: utf-8 -*-
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
import datetime as _dt
import statistics as _stats
import textwrap as _tw

from pathlib import Path
from dataclasses import is_dataclass, fields as dc_fields
from typing import ( 
    Any, 
    Dict,
    Iterable, 
    List, 
    Optional, 
    Sequence, 
    Tuple, 
    Union, 
    Mapping, 
    Generic,
    TypeVar,
    Iterator,
    Callable 
)
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
        name = (
            f"{self.__class__.__module__}."
            f"{self.__class__.__name__}"
        )
        self._logger = (
            logger if logger is not None
            else self._logger_factory(name)
        )

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
    _section: Optional[str] = None
    _float_fmt: str = "{: .6E}"
    _indent: int = 2
    _show_in_str: Optional[Sequence[str]] = None
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
        parts: List[str] = []
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
    def to_dict(self, *, exclude_none: bool = True) -> Dict[str, Any]:
        out: Dict[str, Any] = {}
        for k, v in self._iter_kv(exclude_none):
            out[k] = v
        return out

    @classmethod
    def from_dict(
        cls, data: Mapping[str, Any]
    ) -> "EDIComponentBase":
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
        only: Optional[Sequence[str]] = None,
        indent: Optional[int] = None,
    ) -> List[str]:
        pad = indent if indent is not None else self._indent
        kv = list(self._iter_kv(True))
        if only is not None:
            order = list(only)
            kv = [(k, v) for (k, v) in kv if k in order]
            kv.sort(key=lambda it: order.index(it[0]))
        return [self._format_kv(k, v, pad) for (k, v) in kv]

    def to_text(self) -> str:
        return str(self)

    def update(self, /, **kw: Any) -> "EDIComponentBase":
        for k, v in kw.items():
            setattr(self, k, v)
        self.validate()
        return self

    def clone(self, /, **ov: Any) -> "EDIComponentBase":
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
    def _field_names(cls) -> List[str]:
        if is_dataclass(cls):
            return [f.name for f in dc_fields(cls)]
        ann = getattr(cls, "__annotations__", None) or {}
        if ann:
            return list(ann.keys())
        return []

    def _iter_kv(self, exclude_none: bool) -> Iterator[Tuple[str, Any]]:
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
                (f"{x:.6g}" if isinstance(x, float) else repr(x))
                for x in value[:4]
            )
            if len(value) > 4:
                prev += ", …"
            return f"[{prev}]"
        return repr(value)


class EdiFileBase(EDIComponentBase, ABC):
    path: Optional[Path] = None
    strict_validate: bool = True
    block_size: int = 6
    number_fmt: str = " 15.6e"
    data_header_tpl: str = ">!****{title}****!\n"

    sections: Dict[str, EDIComponentBase] = {}

    _log = Base._logger_factory(__name__)

    def __init__(
        self,
        *,
        path: Optional[PathLike] = None,
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
    def from_file(
        cls, file: PathLike, **kw: Any
    ) -> "EdiFileBase":
        inst = cls(path=Path(file), **kw)
        inst.read()
        return inst

    def write_file(
        self,
        file: Optional[PathLike] = None,
        *,
        overwrite: bool = True,
    ) -> Path:
        target = Path(file) if file is not None else self._default_out()
        if target.exists() and not overwrite:
            raise FileExistsError(
                f"{target} exists and overwrite=False."
            )
        text = self.compose()
        if isinstance(text, (list, tuple)):
            text = "".join(text)
        target.parent.mkdir(parents=True, exist_ok=True)
        with target.open("w", encoding="utf-8", newline="\n") as f:
            f.write(text)
        self._log.info("Wrote EDI to %s", str(target))
        return target

    @abstractmethod
    def read(self) -> "EdiFileBase":
        ...

    @abstractmethod
    def compose(self) -> Union[str, List[str]]:
        ...

    # ------- registry ----------
    def add_section(self, name: str, s: EDIComponentBase) -> None:
        key = self._normalize_section_name(name)
        self.sections[key] = s

    def get_section(self, name: str) -> Optional[EDIComponentBase]:
        return self.sections.get(self._normalize_section_name(name))

    def remove_section(self, name: str) -> None:
        self.sections.pop(self._normalize_section_name(name), None)

    def iter_sections(
        self,
    ) -> Iterable[Tuple[str, EDIComponentBase]]:
        return self.sections.items()

    # ------- proxies -----------
    @property
    def dataid(self) -> Optional[str]:
        head = self.get_section("head")
        return getattr(head, "dataid", None) if head else None

    @property
    def lat(self) -> Optional[float]:
        head = self.get_section("head")
        dm = self.get_section("definemeasurement")
        if head and getattr(head, "lat", None) is not None:
            return getattr(head, "lat", None)
        return getattr(dm, "reflat", None) if dm else None

    @property
    def lon(self) -> Optional[float]:
        head = self.get_section("head")
        dm = self.get_section("definemeasurement")
        if head and getattr(head, "long", None) is not None:
            return getattr(head, "long", None)
        return getattr(dm, "reflong", None) if dm else None

    @property
    def elev(self) -> Optional[float]:
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
        stem = (self.dataid or "output")
        return Path(f"{stem}.edi")

    # ------- fmt helpers -------
    @staticmethod
    def normalize_section_title(title: str) -> str:
        t = (title or "").strip()
        return t.lstrip(">").upper()

    def format_block_header(self, title: str) -> str:
        return self.data_header_tpl.format(
            title=self.normalize_section_title(title)
        )

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
        width: Optional[int] = None,
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
        data: Iterable[Union[float, int]],
        *,
        per_line: Optional[int] = None,
        header_comment: bool = True,
        count_comment: bool = True,
    ) -> List[str]:
        t = self.normalize_section_title(title)
        vals = list(data or [])
        n = len(vals)
        per = per_line or self.block_size
        out: List[str] = []
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

        out.append(
            fmtnums(vals, per_line=per, indent=2) + "\n"
        )
        return out

    # ------- misc --------------
    @staticmethod
    def ensure_descending_frequency(
        freq: Iterable[float],
    ) -> Tuple[List[float], Optional[List[int]]]:
        f = list(freq or [])
        if not f:
            return f, None
        if len(f) < 2 or f[0] >= f[-1]:
            return f, None
        idx = list(range(len(f)))[::-1]
        return f[::-1], idx

    def to_dict(self) -> Dict[str, Any]:
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
    def from_dict(cls, data: Mapping[str, Any]) -> "EdiFileBase":
        inst = cls(
            path=Path(data["path"]) if data.get("path") else None,
            strict_validate=bool(data.get("strict_validate", True)),
        )
        inst.block_size = int(data.get("block_size", inst.block_size))
        inst.number_fmt = str(data.get("number_fmt", inst.number_fmt))
        inst.data_header_tpl = str(
            data.get("data_header_tpl", inst.data_header_tpl)
        )
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
        return (
            f"<{cls} path={self.path!r} "
            f"sections={list(self.sections.keys())!r}>"
        )

    @staticmethod
    def _normalize_section_name(name: str) -> str:
        return (name or "").strip().lower()

    def compose_to_string(self) -> str:
        text = self.compose()
        return "".join(text) if isinstance(text, list) else str(text)



class CollectionBase(Generic[_T], Sequence[_T]):
    """
    Generic, immutable*¹* container.
    
    A thin, generic wrapper that groups together homogeneous
    Python objects and offers:

    * sequence-like behaviour (``len / iter / index``);
    * convenience helpers to add, filter, map or reduce items;
    * a compact, readable ``__repr__``/``__str__`` that never
      exceeds 62 characters per physical line;
    * a `.summary()` routine able to gather simple statistics
      (count / min / mean / max) for any numeric attribute or
      for values returned by a callable;
    * automatic timestamping → you always know when the
      collection was instantiated.

    All downstream domain-specific containers (e.g. an
    ``EdiCollection`` that handles many :class:`~pycsamt.seg.edi.Edi`
    objects) can build on that skeleton instead of re-implementing
    the plumbing each time.
    

    Parameters
    ----------
    items
        An iterable of objects to wrap.  If *None*, the
        collection starts empty.
    name
        Arbitrary label shown in the string representations.

    Notes
    -----
    *¹*  The underlying list is kept private.  Public mutators
    (`add`, `extend`, `clear`) always **replace** that internal
    list rather than mutating it in-place; this prevents hidden
    external references from becoming stale.
    """

    def __init__(
        self,
        items: Optional[Iterable[_T]] = None,
        *,
        name: Optional[str] = None,
    ) -> None:
        self._created: _dt.datetime = _dt.datetime.utcnow()
        self._name: str = name or self.__class__.__name__.lower()
        self._items: List[_T] = list(items) if items is not None else []

    #  basic sequence API #
    def __len__(self) -> int:  # noqa: D401  (pydocstyle: “missing period”)
        return len(self._items)

    def __iter__(self) -> Iterator[_T]:
        return iter(self._items)

    def __getitem__(self, key: int) -> _T:
        return self._items[key]

    def __contains__(self, item: object) -> bool:  # type: ignore[override]
        return item in self._items

    # pretty print 
    _WRAP = _tw.TextWrapper(
        width=62,
        subsequent_indent=" " * 4,
        break_on_hyphens=False,
        replace_whitespace=False,
    )

    def __repr__(self) -> str:  # noqa: D401
        cls = self.__class__.__name__
        return (
            f"<{cls} [{self._name!s}] "
            f"size={len(self)}, created={self._created.isoformat()}Z>"
        )

    def __str__(self) -> str:  # noqa: D401
        head = (
            f"{self._name} collection  •  {len(self)} item"
            f"{'' if len(self) == 1 else 's'}"
        )
        items = ", ".join(map(self._render_item, self._preview()))
        body = f"[{items}]" if items else "[empty]"
        return "\n".join(self._WRAP.wrap(head + "  " + body))


    # mutators – all return **self** to enable method-chaining
    def add(self, item: _T) -> "CollectionBase[_T]":
        """Append *item* and return *self*."""
        self._items = [*self._items, item]
        return self

    def extend(self, items: Iterable[_T]) -> "CollectionBase[_T]":
        """Append many *items* and return *self*."""
        self._items = [*self._items, *items]
        return self

    def clear(self) -> "CollectionBase[_T]":
        """Drop every stored object and return *self*."""
        self._items = []
        return self

    # higher order 
    def filter(
        self,
        predicate: Callable[[_T], bool],
        *,
        name: Optional[str] = None,
    ) -> "CollectionBase[_T]":
        """Return a **new** collection keeping objects that pass *predicate*."""
        return self.__class__(filter(predicate, self._items), name=name)

    def map(
        self,
        fn: Callable[[_T], Any],
        *,
        name: Optional[str] = None,
    ) -> "CollectionBase[Any]":
        """
        Return a **new** collection made of ``fn(item)`` for every
        element.  (Guarantees lazy *fn* evaluation through a
        generator comprehension.)
        """
        return CollectionBase(
            (fn(x) for x in self._items), name=name or self._name)

    # misc 
    def summary(
        self,
        attr: str | Callable[[_T], Any],
        *,
        numeric_only: bool = True,
    ) -> dict[str, Any]:
        """
        Compute simple stats on attribute *attr*.

        *attr* can be either:

        * the attribute **name** to pull from each object; or
        * a callable that returns the value to analyse.

        Returns a ``dict`` holding *count*, *min*, *mean*,
        *max* (keys absent if the data are not numeric).
        """
        if callable(attr):
            values = [attr(obj) for obj in self._items]
            name = getattr(attr, "__name__", "callable")
        else:
            values = [getattr(obj, attr, None) for obj in self._items]
            name = str(attr)

        out: dict[str, Any] = {"field": name, "count": len(values)}
        numeric_values = [v for v in values if isinstance(v, (int, float))]

        if numeric_values and (not numeric_only or len(numeric_values) == len(values)):
            out.update(
                min=min(numeric_values),
                mean=_stats.fmean(numeric_values),
                max=max(numeric_values),
            )

        return out

    #  internals #
    def _preview(self, n: int = 6) -> List[_T]:
        """Return *n* first items for rendering."""
        return self._items[: n - 1] + (["…"] if len(self) > n else [])

    @staticmethod
    def _render_item(obj: Any) -> str:
        """Render one object into a short single-line string."""
        if isinstance(obj, Path):
            return obj.name
        return str(obj)

