"""General result container for public pyCSAMT APIs."""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

import pandas as pd

from ..property import PyCSAMTObject
from .frame import APIFrame, wrap_frame

__all__ = ["APIResult", "wrap_result"]


class APIResult(PyCSAMTObject):
    """Container for named result parts, including one or more tables."""

    __repr_fields__ = ("name", "kind", "keys")

    def __init__(
        self,
        /,
        **items: Any,
    ) -> None:
        self.name = str(items.pop("name", "result"))
        self.kind = items.pop("kind", None)
        self.meta = dict(items.pop("meta", {}) or {})
        self._items = dict(items)

    @property
    def keys(self) -> tuple[str, ...]:
        """Return stored item names."""
        return tuple(self._items.keys())

    def __getitem__(self, key: str) -> Any:
        return self._items[key]

    def __setitem__(self, key: str, value: Any) -> None:
        self._items[str(key)] = value

    def __contains__(self, key: object) -> bool:
        return key in self._items

    def __getattribute__(self, name: str) -> Any:
        if not name.startswith("_"):
            items = object.__getattribute__(
                self, "__dict__"
            ).get("_items")
            if isinstance(items, dict) and name in items:
                return items[name]
        return object.__getattribute__(self, name)

    def __iter__(self):
        return iter(self._items)

    def __len__(self) -> int:
        return len(self._items)

    def __dir__(self) -> list[str]:
        return sorted(set(super().__dir__()) | set(self._items.keys()))

    def items(self):
        """Return ``(name, value)`` pairs."""
        return self._items.items()

    def tables(self) -> dict[str, APIFrame]:
        """Return only table-like result parts."""
        items = object.__getattribute__(self, "_items")
        return {
            k: v for k, v in items.items()
            if isinstance(v, APIFrame)
        }

    def update_meta(self, /, **kwargs: Any) -> APIResult:
        """Update metadata in-place and return ``self``."""
        self.meta.update(kwargs)
        return self

    def _summary_text(self) -> str:
        """Return a compact result display."""
        lines = [f"APIResult: {self.name}"]
        if self.kind:
            lines.append(f"kind: {self.kind}")
        keys = tuple(object.__getattribute__(self, "_items").keys())
        lines.append(f"items: {', '.join(keys) if keys else '-'}")
        tables = APIResult.tables(self)
        if tables:
            shapes = ", ".join(
                f"{key}={value.shape[0]}x{value.shape[1]}"
                for key, value in tables.items()
            )
            lines.append(f"tables: {shapes}")
        return "\n".join(lines)

    def __str__(self) -> str:
        return self._summary_text()


def wrap_result(
    result: Mapping[str, Any],
    *,
    name: str = "result",
    kind: str | None = None,
    meta: Mapping[str, Any] | None = None,
    wrap_tables: bool = True,
) -> APIResult:
    """Wrap a mapping of result parts as an :class:`APIResult`."""
    items: dict[str, Any] = {}
    for key, value in result.items():
        if wrap_tables and isinstance(value, pd.DataFrame):
            items[str(key)] = wrap_frame(value, name=str(key), kind=kind)
        else:
            items[str(key)] = value
    return APIResult(name=name, kind=kind, meta=dict(meta or {}), **items)
