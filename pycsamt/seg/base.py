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
import logging
from pathlib import Path
from dataclasses import is_dataclass, fields as dc_fields
from dataclasses import dataclass, field  
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
try:
    from ..log.logger import get_logger as _get_logger
except Exception:  # pragma: no cover
    _get_logger = None

from .properties import IsEdi 

PathLike = Union[str, Path]


__all__ = [
    "CollectionBase", 
    "EDIComponentBase", 
    "EdiFileBase", 
  ]

_T = TypeVar("_T")


class EDIComponentBase:
    """Base class for SEG-EDI components.

    This mixin provides consistent ``__repr__`` and ``__str__``
    plus helpers to serialize component attributes into EDI-like
    ``KEY=VALUE`` lines.

    Subclasses may override:
    - ``_section``: name of the EDI block (e.g. ``"HEAD"``).
    - ``_float_fmt``: numeric format, default ``"%.6E"``.
    - ``_indent``: spaces for left indent when rendering lines.
    - ``_show_in_str``: ordered subset of keys to include in
      ``__str__``. If ``None``, all non-None fields are used.
    - ``validate``: to enforce invariants.

    Notes
    -----
    - Field order is derived from type hints (``__annotations__``)
      when available; otherwise it falls back to attribute order
      in ``__dict__`` at runtime.
    - ``None`` values are skipped in output lines.
    - Strings with spaces or quotes are double-quoted.
    - Floats are formatted with scientific notation by default.
    - ``__repr__`` is intentionally short (truncated) and geared
      for debugging, while ``__str__`` prints an EDI-style block.

    Examples
    --------
    >>> class Head(EDIComponentBase):
    ...     _section = "HEAD"
    ...     DATAID: str
    ...     STDVERS: str = "SEG 1.0"
    ...     EMPTY: float = 1.0e32
    ...
    >>> h = Head()
    >>> h.DATAID = "E1_2"
    >>> print(str(h)).startswith(">HEAD")
    True
    >>> "DATAID" in repr(h)
    True
    """

    # EDI rendering preferences (subclasses may override)
    _section: Optional[str] = None
    _float_fmt: str = "{:.6E}"
    _indent: int = 2
    _show_in_str: Optional[Sequence[str]] = None
    _repr_max_items: int = 6
    _repr_max_chars: int = 120

    # core protocol 
    def __repr__(self) -> str:
        cls = self.__class__.__name__
        items = self._iter_kv(exclude_none=True)
        # truncate for readability
        preview: List[str] = []
        total_chars = 0
        for i, (k, v) in enumerate(items):
            if i >= self._repr_max_items:
                preview.append("…")
                break
            frag = f"{k}={self._repr_value(v)}"
            total_chars += len(frag)
            if total_chars > self._repr_max_chars:
                preview.append("…")
                break
            preview.append(frag)
        inside = ", ".join(preview)
        return f"{cls}({inside})"

    def __str__(self) -> str:
        # If section is defined, render EDI block. Otherwise fallback
        # to a simple one-line KV view.
        if self._section:
            lines = [f">{self._section}"]
            for line in self.to_lines(only=self._show_in_str):
                lines.append(line)
            return "\n".join(lines)
        # fallback
        kv = " ".join(self.to_lines(only=self._show_in_str))
        return f"{self.__class__.__name__} {kv}"

    # public helpers 
    def to_dict(self, *, exclude_none: bool = True) -> Dict[str, Any]:
        """Return attributes as a dict.

        Parameters
        ----------
        exclude_none : bool, default=True
            If ``True``, drop keys whose values are ``None``.
        """
        out: Dict[str, Any] = {}
        for k, v in self._iter_kv(exclude_none=exclude_none):
            out[k] = v
        return out

    def to_lines(
        self,
        *,
        only: Optional[Sequence[str]] = None,
        indent: Optional[int] = None,
    ) -> List[str]:
        """Format attributes as EDI-like ``KEY=VALUE`` lines.

        Parameters
        ----------
        only : sequence of str, optional
            Subset of keys to include (order respected).
        indent : int, optional
            Number of leading spaces; defaults to ``self._indent``.
        """
        pad = indent if indent is not None else self._indent
        kv = list(self._iter_kv(exclude_none=True))
        if only is not None:
            order = list(only)
            # keep declared order; filter missing
            kv = [(k, v) for k, v in kv if k in order]
            # ensure explicit order
            kv.sort(key=lambda kv_: order.index(kv_[0]))
        return [self._format_kv(k, v, pad) for k, v in kv]

    def to_text(self) -> str:
        """Render the component as EDI block text.

        Returns
        -------
        str
            Block text. If ``_section`` is not set, a one-line
            ``KEY=VALUE`` sequence is returned.
        """
        return str(self)

    def update(self, /, **kwargs: Any) -> "EDIComponentBase":
        """Update attributes in place and revalidate.

        Returns
        -------
        self
        """
        for k, v in kwargs.items():
            setattr(self, k, v)
        self.validate()
        return self

    def clone(self, /, **overrides: Any) -> "EDIComponentBase":
        """Create a shallow copy with optional overrides."""
        new = self.__class__()
        for k, v in self._iter_kv(exclude_none=False):
            setattr(new, k, v)
        for k, v in overrides.items():
            setattr(new, k, v)
        new.validate()
        return new

    # subclass hooks 
    def validate(self) -> None:
        """Hook for subclasses to enforce invariants.

        The default implementation does nothing. Override to raise
        a ``ValueError`` (or a domain-specific error) when state
        is invalid.
        """
        return None

    # internals: formatting 
    @classmethod
    def _field_names(cls) -> List[str]:
        # Prefer dataclass field order if available
        if is_dataclass(cls):
            return [f.name for f in dc_fields(cls)]
        # fall back to annotations if present
        ann = getattr(cls, "__annotations__", None) or {}
        if ann:
            return list(ann.keys())
        # last resort: empty; subclasses without hints will rely
        # on instance __dict__ iteration order at runtime
        return []

    def _iter_kv(self, *, exclude_none: bool) -> Iterable[Tuple[str, Any]]:
        # Start from class-declared field order (dataclass or type hints)
        seen = set()
        for name in self._field_names():
            if not hasattr(self, name):
                continue
            val = getattr(self, name)
            if exclude_none and val is None:
                continue
            seen.add(name)
            yield name, val
        # Include any dynamically added attributes not declared in type hints
        for name, val in self.__dict__.items():
            if name.startswith("_") or name in seen:
                continue
            if exclude_none and val is None:
                continue
            yield name, val

    def _format_kv(self, key: str, value: Any, indent: int) -> str:
        return " " * indent + f"{key}={self._format_value(value)}"

    def _format_value(self, value: Any) -> str:
        # Strings: quote if spaces or quotes present.
        if isinstance(value, str):
            v = value.strip()
            if any(ch.isspace() for ch in v) or ('"' in v):
                return '"' + v.replace('"', "") + '"'
            return v
        # Booleans: render as 0/1 like many EDI writers do.
        if isinstance(value, bool):
            return "1" if value else "0"
        # Numbers: float -> scientific, int -> as-is
        if isinstance(value, float):
            if value != value:  # NaN
                return "NaN"
            return self._float_fmt.format(value)
        if isinstance(value, int):
            return str(value)
        # Sequences: join by space (subclasses may override)
        if isinstance(value, (list, tuple)):
            return " ".join(self._format_value(v) for v in value)
        # Fallback: default repr
        return str(value)

    @staticmethod
    def _repr_value(value: Any) -> str:
        if isinstance(value, str):
            v = value.replace('"', "")
            return f'"{v}"' if len(v) > 20 or " " in v else v
        if isinstance(value, float):
            return f"{value:.6g}"
        if isinstance(value, (list, tuple)):
            preview = ", ".join(
                (f"{x:.6g}" if isinstance(x, float) else repr(x)) for x in value[:4]
            )
            if len(value) > 4:
                preview += ", …"
            return f"[{preview}]"
        return repr(value)


def _logger_for(name: str) -> logging.Logger:
    if _get_logger is not None:
        return _get_logger(name)
    # Fallback to std logging (no-op if user doesn't configure it)
    lg = logging.getLogger(name)
    if not lg.handlers:
        lg.addHandler(logging.NullHandler())
    return lg


@dataclass
class EdiFileBase(EDIComponentBase, ABC):
    """
    A robust base for EDI file objects.

    Subclass this to create the concrete `Edi` class. It manages:
      • path & validation (via IsEdi)
      • section/component registry
      • pretty __repr__/__str__ and a short .summary()
      • read()/write() lifecycle hooks
      • composition helpers for blocks and fixed-width numeric data
      • metadata proxies (dataid, lat/lon/elev) if sections provide them
    """

    # --- File identity & config
    path: Optional[Path] = None
    strict_validate: bool = True
    block_size: int = 6                  # default numbers per line for data blocks
    number_fmt: str = " 15.6e"           # default float fmt for data blocks
    data_header_tpl: str = ">!****{title}****!\n"  # comment headers for blocks

    # --- Registry of sections (Head, Info, DefineMeasurement,
    # MTEMAP, Spectra, TimeSeries, Z, Tip…)
    sections: Dict[str, EDIComponentBase] = field(
        default_factory=dict, repr=False)

    # --- Internal logger
    _log: logging.Logger = field(
        default_factory=lambda: _logger_for(
            __name__), init=False, repr=False)


    # Construction & initialization
    def __post_init__(self) -> None:
        # Normalize path
        if isinstance(self.path, (str, Path)):
            self.path = Path(self.path) if self.path is not None else None

        # Optionally validate EDI file early
        if self.path is not None and self.strict_validate:
            self._validate_path(self.path)

    # Public API (to be used by concrete Edi subclass)
    @classmethod
    def from_file(cls, file: PathLike, **kwargs: Any) -> "EdiFileBase":
        """
        Build an instance from an on-disk EDI file and call .read().

        Subclasses should implement .read() to populate sections.
        """
        inst = cls(path=Path(file), **kwargs)
        inst.read()
        return inst

    def write_file(
            self, file: Optional[PathLike] = None, *, 
            overwrite: bool = True) -> Path:
        """
        Compose current content and write to disk.

        Subclasses must implement .compose() to return a string (or an
        iterable of lines). This base method handles path resolution,
        overwrite behavior, and encoding.
        """
        target = Path(file) if file is not None else self._default_output_path()
        if target.exists() and not overwrite:
            raise FileExistsError(f"{target} already exists and overwrite=False.")

        # Compose
        text = self.compose()
        if isinstance(text, (list, tuple)):
            text = "".join(text)

        # Ensure parent directory
        target.parent.mkdir(parents=True, exist_ok=True)
        with target.open("w", encoding="utf-8", newline="\n") as f:
            f.write(text)

        self._log.info("Wrote EDI to %s", str(target))
        return target


    @abstractmethod
    def read(self) -> "EdiFileBase":
        """
        Read from self.path (which should be set) and populate `sections`.
        Must return self.
        """
        ...

    @abstractmethod
    def compose(self) -> Union[str, List[str]]:
        """
        Generate the full EDI text (or list of lines) from the current `sections`.
        Must be implemented by the concrete Edi class.
        """
        ...

    # Section registry helpers
    def add_section(self, name: str, section: EDIComponentBase) -> None:
        """Register/replace a section object under a normalized key."""
        key = self._normalize_section_name(name)
        self.sections[key] = section

    def get_section(self, name: str) -> Optional[EDIComponentBase]:
        """Retrieve a section by name (normalized); returns None if missing."""
        return self.sections.get(self._normalize_section_name(name))

    def remove_section(self, name: str) -> None:
        """Remove a section by name; ignores missing keys."""
        self.sections.pop(self._normalize_section_name(name), None)

    def iter_sections(self) -> Iterable[Tuple[str, EDIComponentBase]]:
        """Iterate over (name, section) pairs in insertion order."""
        return self.sections.items()


    # Metadata proxies (if present)
    @property
    def dataid(self) -> Optional[str]:
        head = self.get_section("head")
        return getattr(head, "dataid", None) if head is not None else None

    @property
    def lat(self) -> Optional[float]:
        head = self.get_section("head")
        dm = self.get_section("definemeasurement")
        return (
            getattr(head, "lat", None)
            if head is not None and getattr(head, "lat", None) is not None
            else getattr(dm, "reflat", None) if dm is not None else None
        )

    @property
    def lon(self) -> Optional[float]:
        head = self.get_section("head")
        dm = self.get_section("definemeasurement")
        return (
            getattr(head, "long", None)
            if head is not None and getattr(head, "long", None) is not None
            else getattr(dm, "reflong", None) if dm is not None else None
        )

    @property
    def elev(self) -> Optional[float]:
        head = self.get_section("head")
        dm = self.get_section("definemeasurement")
        return (
            getattr(head, "elev", None)
            if head is not None and getattr(head, "elev", None) is not None
            else getattr(dm, "refelev", None) if dm is not None else None
        )

    # Validation
    def _validate_path(self, file: Path) -> None:
        """Use IsEdi to validate file existence and EDI-ness."""
        # Shallow validation (extension) if strict_validate is False
        IsEdi._assert_edi(str(file), deep=self.strict_validate)

    # Default output path (can be overridden)
    def _default_output_path(self) -> Path:
        """
        Resolve a reasonable default output file name.

        Uses dataid (if available) and today's date: <DATAID>_<YYYYMMDD>.edi
        Falls back to 'output_<timestamp>.edi'.
        """
        stem = (self.dataid or "output") + "_" + _dt.datetime.utcnow(
            ).strftime("%Y%m%d")
        return Path(f"{stem}.edi")

    # Formatting helpers
    @staticmethod
    def normalize_section_title(title: str) -> str:
        """Uppercase title, ensure no leading '>'."""
        t = (title or "").strip()
        return t.lstrip(">").upper()

    def format_block_header(self, title: str) -> str:
        """Return a 'comment' header line like '>!****TITLE****!\\n'."""
        return self.data_header_tpl.format(
            title=self.normalize_section_title(title))

    @staticmethod
    def format_section_head(head: str) -> str:
        """Return a section head marker (e.g., '>=MTSECT\\n' or '>INFO\\n')."""
        head = head.strip()
        if not head.startswith(">"):
            head = ">" + head
        return head.upper() + ("\n" if not head.endswith("\n") else "")

    def format_kv(
            self, key: str, value: Any, *, 
            quote: bool = False, 
            width: Optional[int] = None
            ) -> str:
        """
        Format '  KEY=VALUE\\n' with optional quoting and fixed width 
        right alignment.
        """
        key_up = key.strip().upper()
        if value is None:
            val = ""
        else:
            val = f'"{value}"' if quote and isinstance(value, str) else str(value)
        if width is not None:
            val = f"{val:>{width}}"
        return f"  {key_up}={val}\n"

    def format_data_block(
        self,
        title: str,
        data: Iterable[Union[float, int]],
        *,
        per_line: Optional[int] = None,
        num_fmt: Optional[str] = None,
        header_comment: bool = True,
        count_comment: bool = True,
    ) -> List[str]:
        """
        Format a numeric data block (e.g., FREQ, ZXXR, RHOXY) as a list of lines.

        Parameters
        ----------
        title : str
            The data block name (e.g., 'FREQ', 'ZXXR').
        data : Iterable[float|int]
            Numeric sequence to format.
        per_line : int
            Numbers per line (default to self.block_size).
        num_fmt : str
            Number format (default to self.number_fmt).
        header_comment : bool
            If True, include a '>!****TITLE****!' comment line.
        count_comment : bool
            If True, append '  //N' count on the first line after
            the title if EDI style needs it.

        Returns
        -------
        List[str] : lines ready to join/write.
        """
        title_norm = self.normalize_section_title(title)
        nums = list(data or [])
        n = len(nums)
        per_line = per_line or self.block_size
        num_fmt = num_fmt or self.number_fmt

        out: List[str] = []
        if header_comment:
            out.append(self.format_block_header(title_norm))

        # Title line (optional with count suffix)
        if count_comment:
            out.append(self.format_section_head(f">{title_norm}  //{n}"))
        else:
            out.append(self.format_section_head(f">{title_norm}"))

        # Body lines
        if n == 0:
            return out  # empty set represented by header only

        # Chunk and format numbers
        def _fmt(x: Union[float, int]) -> str:
            try:
                return f"{float(x):{num_fmt}}"
            except Exception:
                return f"{0.0:{num_fmt}}"

        for i in range(0, n, per_line):
            chunk = nums[i : i + per_line]
            out.append("".join(_fmt(v) for v in chunk) + "\n")

        return out


    # Frequency normalization helper
    @staticmethod
    def ensure_descending_frequency(
        freq: Iterable[float]
        ) -> Tuple[List[float], Optional[List[int]]]:
        """
        Return frequency in descending order. Also returns an index permutation
        if a reordering was needed; otherwise permutation is None.

        Useful so that Z/tipper arrays can be permuted consistently by caller.
        """
        f = list(freq or [])
        if not f:
            return f, None
        if len(f) < 2 or f[0] >= f[-1]:
            return f, None
        # ascending -> reverse
        idx = list(range(len(f)))[::-1]
        return f[::-1], idx

    # Serialization
    def to_dict(self) -> Dict[str, Any]:
        """Serialize this file object + sections."""
        d = {
            "path": str(self.path) if self.path is not None else None,
            "strict_validate": self.strict_validate,
            "block_size": self.block_size,
            "number_fmt": self.number_fmt,
            "data_header_tpl": self.data_header_tpl,
            "sections": {k: v.to_dict() if hasattr(
                v, "to_dict") else dict(
                    v.__dict__) for k, v in self.sections.items()},
        }
        return d

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> "EdiFileBase":
        """
        Rehydrate an instance from a dict produced by .to_dict().
        Sections are restored as generic EdiComponentBase instances unless
        the subclass overrides this method to dispatch to concrete types.
        """
        inst = cls(
            path=Path(data["path"]) if data.get("path") else None,
            strict_validate=bool(data.get("strict_validate", True)),
        )
        inst.block_size = int(data.get("block_size", inst.block_size))
        inst.number_fmt = str(data.get("number_fmt", inst.number_fmt))
        inst.data_header_tpl = str(data.get(
            "data_header_tpl", inst.data_header_tpl))

        sections = data.get("sections") or {}
        for k, v in sections.items():
            comp = EDIComponentBase.from_dict(v) if isinstance(
                v, Mapping) else EDIComponentBase()
            inst.add_section(k, comp)
        return inst


    # Representations
    def summary(self) -> str:
        """Short human-readable summary for __str__."""
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

    def __str__(self) -> str:  # pragma: no cover - simple formatting
        return self.summary()

    def __repr__(self) -> str:  # pragma: no cover - simple formatting
        cls = self.__class__.__name__
        return f"<{cls} path={self.path!r} sections={list(self.sections.keys())!r}>"


    # Utilities
    @staticmethod
    def _normalize_section_name(name: str) -> str:
        return (name or "").strip().lower()

    # Handy for tests or string output without touching disk
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

