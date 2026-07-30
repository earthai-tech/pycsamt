"""Pandas-friendly dataframe result objects for public pyCSAMT APIs."""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from functools import wraps
from typing import Any, Callable

import numpy as np
import pandas as pd

from ..property import PyCSAMTObject

__all__ = [
    "APIFrame",
    "FrameProfile",
    "api_frame",
    "maybe_wrap_frame",
    "wrap_frame",
]


@dataclass(frozen=True)
class FrameProfile:
    """Small immutable profile describing a dataframe."""

    rows: int
    columns: int
    column_names: tuple[str, ...]
    numeric_columns: tuple[str, ...]
    missing_cells: int
    missing_fraction: float
    memory_bytes: int

    @classmethod
    def from_frame(cls, df: pd.DataFrame) -> FrameProfile:
        total = int(df.size)
        missing = int(df.isna().sum().sum()) if total else 0
        numeric = tuple(
            str(c) for c in df.select_dtypes(include="number").columns
        )
        try:
            memory = int(df.memory_usage(deep=True).sum())
        except Exception:
            memory = 0
        return cls(
            rows=int(df.shape[0]),
            columns=int(df.shape[1]),
            column_names=tuple(str(c) for c in df.columns),
            numeric_columns=numeric,
            missing_cells=missing,
            missing_fraction=(float(missing) / float(total)) if total else 0.0,
            memory_bytes=memory,
        )

    @property
    def shape(self) -> tuple[int, int]:
        """Return ``(rows, columns)``."""
        return self.rows, self.columns

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-friendly profile dictionary."""
        return {
            "rows": self.rows,
            "columns": self.columns,
            "column_names": list(self.column_names),
            "numeric_columns": list(self.numeric_columns),
            "missing_cells": self.missing_cells,
            "missing_fraction": self.missing_fraction,
            "memory_bytes": self.memory_bytes,
        }


class APIFrame(PyCSAMTObject):
    """A pyCSAMT dataframe view that keeps pandas behavior intact.

    ``APIFrame`` is a thin wrapper around :class:`pandas.DataFrame`.
    The dataframe remains available as :attr:`df`; pyCSAMT adds compact
    display, metadata, units, statistics, and convenient conversion helpers.
    """

    __repr_fields__ = ("name", "kind", "source", "shape")
    __repr_exclude__ = PyCSAMTObject.__repr_exclude__ | {"_df"}

    def __init__(
        self,
        data: Any = None,
        *,
        name: str | None = None,
        kind: str | None = None,
        source: Any = None,
        units: Mapping[str, str] | None = None,
        meta: Mapping[str, Any] | None = None,
        description: str | None = None,
        copy: bool = False,
        **frame_kwargs: Any,
    ) -> None:
        self._df = self._coerce_frame(data, copy=copy, **frame_kwargs)
        self.name = name or "dataframe"
        self.kind = kind
        self.source = source
        self.units = dict(units or {})
        self.meta = dict(meta or {})
        self.description = description

    @staticmethod
    def _coerce_frame(
        data: Any,
        *,
        copy: bool = False,
        **kwargs: Any,
    ) -> pd.DataFrame:
        if data is None:
            return pd.DataFrame(**kwargs)
        if isinstance(data, APIFrame):
            return data.df.copy() if copy else data.df
        if isinstance(data, pd.DataFrame):
            return data.copy() if copy else data
        return pd.DataFrame(data, **kwargs)

    @classmethod
    def from_records(
        cls,
        records: Iterable[Mapping[str, Any]],
        *,
        columns: Sequence[str] | None = None,
        **kwargs: Any,
    ) -> APIFrame:
        """Build an ``APIFrame`` from record dictionaries."""
        df = pd.DataFrame.from_records(records, columns=columns)
        return cls(df, **kwargs)

    @property
    def df(self) -> pd.DataFrame:
        """Return the underlying pandas dataframe."""
        return self._df

    @df.setter
    def df(self, value: Any) -> None:
        self._df = self._coerce_frame(value)

    @property
    def data(self) -> np.ndarray:
        """Return dataframe values as a NumPy array."""
        return self._df.to_numpy()

    @property
    def shape(self) -> tuple[int, int]:
        """Return dataframe shape."""
        return self._df.shape

    @property
    def columns(self) -> pd.Index:
        """Return dataframe columns."""
        return self._df.columns

    @property
    def stats(self) -> FrameProfile:
        """Return a compact dataframe profile."""
        return FrameProfile.from_frame(self._df)

    @property
    def schema(self) -> dict[str, str]:
        """Return column dtype names keyed by column name."""
        return {str(k): str(v) for k, v in self._df.dtypes.items()}

    def __len__(self) -> int:
        return len(self._df)

    def __iter__(self):
        return iter(self._df)

    def __contains__(self, key: object) -> bool:
        return key in self._df

    def __getitem__(self, key: Any) -> Any:
        return self._df.__getitem__(key)

    def __setitem__(self, key: Any, value: Any) -> None:
        if isinstance(key, tuple) and len(key) == 2:
            rows, column = key
            self._df.loc[rows, column] = value
            return
        if isinstance(key, str) and key in self._df.columns:
            self._df.loc[:, key] = value
            return
        self._df[key] = value

    def __array__(self, dtype: Any = None) -> np.ndarray:
        arr = self.data
        return arr.astype(dtype, copy=False) if dtype is not None else arr

    def __getattr__(self, name: str) -> Any:
        if name.startswith("_"):
            raise AttributeError(name)
        df = self.__dict__.get("_df")
        if df is not None and name in df.columns:
            return df[name]
        if df is not None and hasattr(df, name):
            return getattr(df, name)
        raise AttributeError(
            f"{self.__class__.__name__!s} object has no attribute {name!r}"
        )

    def __dir__(self) -> list[str]:
        base = set(super().__dir__())
        base.update(str(c) for c in self._df.columns)
        return sorted(base)

    def __repr__(self) -> str:
        return (
            f"APIFrame(name={self.name!r}, shape={self.shape}, "
            f"columns={list(map(str, self._df.columns[:5]))!r}"
            f"{', ...' if self._df.shape[1] > 5 else ''})"
        )

    def __str__(self) -> str:
        return self.summary()

    def summary(self, *, max_columns: int = 8) -> str:
        """Return a static display summary for printing."""
        profile = self.stats
        cols = list(profile.column_names[:max_columns])
        if len(profile.column_names) > max_columns:
            cols.append("...")
        lines = [f"APIFrame: {self.name}"]
        if self.kind:
            lines.append(f"kind: {self.kind}")
        lines.append(f"shape: {profile.rows} rows x {profile.columns} columns")
        lines.append(f"columns: {', '.join(cols) if cols else '-'}")
        lines.append(f"numeric: {len(profile.numeric_columns)} columns")
        lines.append(f"missing: {profile.missing_fraction:.1%}")
        if self.source is not None:
            lines.append(f"source: {self.source}")
        if self.description:
            lines.append(f"description: {self.description}")
        return "\n".join(lines)

    def profile(self) -> FrameProfile:
        """Return the same object as :attr:`stats`."""
        return self.stats

    def missing(self) -> pd.Series:
        """Return missing value counts by column."""
        return self._df.isna().sum()

    def numeric_stats(self, **kwargs: Any) -> pd.DataFrame:
        """Return pandas ``describe`` for numeric columns."""
        return self._df.describe(**kwargs)

    def to_pandas(self, *, copy: bool = False) -> pd.DataFrame:
        """Return the underlying dataframe, optionally copied."""
        return self._df.copy() if copy else self._df

    def to_numpy(self, *args: Any, **kwargs: Any) -> np.ndarray:
        """Return dataframe values as a NumPy array."""
        return self._df.to_numpy(*args, **kwargs)

    def to_dict(self, *args: Any, **kwargs: Any) -> dict:
        """Delegate to ``DataFrame.to_dict`` by default."""
        if not args and not kwargs:
            kwargs = {"orient": "list"}
        return self._df.to_dict(*args, **kwargs)

    def copy(self, *, deep: bool = True) -> APIFrame:
        """Return a copied view preserving metadata."""
        return self.with_df(self._df.copy(deep=deep))

    def with_df(self, df: Any, **overrides: Any) -> APIFrame:
        """Return a new ``APIFrame`` with another dataframe."""
        params = {
            "name": self.name,
            "kind": self.kind,
            "source": self.source,
            "units": self.units.copy(),
            "meta": self.meta.copy(),
            "description": self.description,
        }
        params.update(overrides)
        return APIFrame(df, **params)

    def update_meta(self, /, **kwargs: Any) -> APIFrame:
        """Update metadata in-place and return ``self``."""
        self.meta.update(kwargs)
        return self

    def set_units(self, /, **kwargs: str) -> APIFrame:
        """Update column units in-place and return ``self``."""
        self.units.update({str(k): str(v) for k, v in kwargs.items()})
        return self


def _default_wrap_frame(
    data: Any,
    *,
    name: str | None = None,
    kind: str | None = None,
    source: Any = None,
    units: Mapping[str, str] | None = None,
    meta: Mapping[str, Any] | None = None,
    description: str | None = None,
    copy: bool = False,
    **frame_kwargs: Any,
) -> APIFrame:
    """Wrap dataframe-like data as an :class:`APIFrame`."""
    if isinstance(data, APIFrame):
        if (
            any(
                v is not None
                for v in (name, kind, source, units, meta, description)
            )
            or copy
        ):
            return APIFrame(
                data.df,
                name=name or data.name,
                kind=kind if kind is not None else data.kind,
                source=source if source is not None else data.source,
                units=units or data.units,
                meta=meta or data.meta,
                description=(
                    description
                    if description is not None
                    else data.description
                ),
                copy=copy,
            )
        return data
    return APIFrame(
        data,
        name=name,
        kind=kind,
        source=source,
        units=units,
        meta=meta,
        description=description,
        copy=copy,
        **frame_kwargs,
    )


def wrap_frame(
    data: Any,
    *,
    name: str | None = None,
    kind: str | None = None,
    source: Any = None,
    units: Mapping[str, str] | None = None,
    meta: Mapping[str, Any] | None = None,
    description: str | None = None,
    copy: bool = False,
    **frame_kwargs: Any,
) -> Any:
    """Wrap dataframe-like data using the configured API view backend."""
    from .config import PYCSAMT_API_VIEW

    return PYCSAMT_API_VIEW.wrap_frame(
        data,
        name=name,
        kind=kind,
        source=source,
        units=units,
        meta=meta,
        description=description,
        copy=copy,
        **frame_kwargs,
    )


def maybe_wrap_frame(
    data: Any,
    *,
    api: bool | None = None,
    name: str | None = None,
    kind: str | None = None,
    source: Any = None,
    units: Mapping[str, str] | None = None,
    meta: Mapping[str, Any] | None = None,
    description: str | None = None,
    copy: bool = False,
    **frame_kwargs: Any,
) -> Any:
    """Conditionally wrap *data* as an :class:`APIFrame`.

    When *api* is ``None`` (the default), the global
    :data:`~pycsamt.api.view.config.PYCSAMT_API_VIEW` configuration decides.
    Pass ``api=True`` to force wrapping, or ``api=False`` to always return
    the raw dataframe regardless of the global setting.
    """
    if api is False:
        return data
    if api is True:
        # Explicit per-call override: always produce an APIFrame.
        return _default_wrap_frame(
            data,
            name=name,
            kind=kind,
            source=source,
            units=units,
            meta=meta,
            description=description,
            copy=copy,
            **frame_kwargs,
        )
    # api is None — defer to global config (may use custom wrapper or disable).
    from .config import PYCSAMT_API_VIEW

    if not PYCSAMT_API_VIEW.enabled():
        return data
    return wrap_frame(
        data,
        name=name,
        kind=kind,
        source=source,
        units=units,
        meta=meta,
        description=description,
        copy=copy,
        **frame_kwargs,
    )


def api_frame(
    _func: Callable[..., Any] | None = None,
    *,
    name: str | None = None,
    kind: str | None = None,
    source: Any = None,
    units: Mapping[str, str] | None = None,
    meta: Mapping[str, Any] | None = None,
    description: str | None = None,
    copy: bool = False,
) -> Callable[..., Any]:
    """Decorate a function so dataframe-like returns become ``APIFrame``."""

    def decorate(func: Callable[..., Any]) -> Callable[..., Any]:
        @wraps(func)
        def wrapper(*args: Any, **kwargs: Any) -> Any:
            out = func(*args, **kwargs)
            if isinstance(out, (APIFrame, pd.DataFrame)):
                return wrap_frame(
                    out,
                    name=name or getattr(func, "__name__", "dataframe"),
                    kind=kind,
                    source=source,
                    units=units,
                    meta=meta,
                    description=description,
                    copy=copy,
                )
            return out

        return wrapper

    if _func is not None:
        return decorate(_func)
    return decorate
