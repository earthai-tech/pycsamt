# -*- coding: utf-8 -*-
"""
pycsamt.seg.collection_base
---------------------------

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
"""
from __future__ import annotations

import datetime as _dt
import statistics as _stats
import textwrap as _tw
from pathlib import Path
from typing import (
    Any,
    Callable,
    Generic,
    Iterable,
    Iterator,
    List,
    Optional,
    Sequence,
    TypeVar,
)

__all__ = ["CollectionBase"]

_T = TypeVar("_T")


class CollectionBase(Generic[_T], Sequence[_T]):
    """
    Generic, immutable*¹* container.

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
        return CollectionBase((fn(x) for x in self._items), name=name or self._name)

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

