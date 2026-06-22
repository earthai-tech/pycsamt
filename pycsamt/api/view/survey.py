"""Public survey view over EDI collections."""

from __future__ import annotations

from typing import Any, Mapping, Sequence

import pandas as pd

from ..property import PyCSAMTObject
from .frame import APIFrame, wrap_frame

__all__ = ["APISurvey"]


class APISurvey(PyCSAMTObject):
    """Friendly public facade over an :class:`EDICollection`."""

    __repr_fields__ = ("name", "n_sites", "source")
    __repr_exclude__ = PyCSAMTObject.__repr_exclude__ | {
        "collection",
        "parser",
    }

    def __init__(
        self,
        collection: Any = None,
        *,
        name: str | None = None,
        source: Any = None,
        parser: Any = None,
        meta: Mapping[str, Any] | None = None,
    ) -> None:
        from ...seg.collection import EDICollection

        if collection is None:
            collection = EDICollection()
        if not isinstance(collection, EDICollection):
            collection = EDICollection(items=collection)
        self.collection = collection
        self.name = name or "survey"
        self.source = source
        self.parser = parser
        self.meta = dict(meta or {})

    @property
    def n_sites(self) -> int:
        """Return number of loaded sites."""
        return len(self.collection)

    @property
    def stations(self) -> list[str]:
        """Return station names."""
        return self.collection.stations()

    @property
    def paths(self) -> list[str]:
        """Return source paths for loaded stations."""
        return self.collection.paths

    @property
    def df(self) -> APIFrame:
        """Return the default summary table."""
        return self.summary()

    @property
    def data(self) -> Any:
        """Return the underlying EDI collection."""
        return self.collection

    def __len__(self) -> int:
        return len(self.collection)

    def __iter__(self):
        return iter(self.collection)

    def __getitem__(self, key: int | str):
        return self.collection[key]

    def __getattr__(self, name: str) -> Any:
        coll = self.__dict__.get("collection")
        if coll is not None and hasattr(coll, name):
            return getattr(coll, name)
        raise AttributeError(name)

    def __dir__(self) -> list[str]:
        coll = self.__dict__.get("collection")
        extra = set(dir(coll)) if coll is not None else set()
        return sorted(set(super().__dir__()) | extra)

    def summary(
        self,
        *,
        fields: Sequence[str] | None = None,
    ) -> APIFrame:
        """Return an editable APIFrame summary of loaded EDI files."""
        rows = self.collection.summary()
        df = pd.DataFrame(rows)
        if fields is not None:
            df = df.loc[:, list(fields)]
        return wrap_frame(
            df,
            name=f"{self.name}_summary",
            kind="edi.summary",
            source=self.source,
            meta={"n_sites": self.n_sites, **self.meta},
        )

    def errors(self) -> list[tuple[Any, BaseException]]:
        """Return parser errors from the last read operation."""
        if self.parser is None:
            return []
        errors = getattr(self.parser, "errors", None)
        return list(errors()) if callable(errors) else []

    def get_site(self, site: str, default: Any = None) -> Any:
        """Return one site by station, stem, filename, or path."""
        try:
            return self.collection._resolve(site)  # noqa: SLF001
        except Exception:
            return default

    def to_dataframe(self, *args: Any, api: bool = True, **kwargs: Any):
        """Return collection dataframe output, wrapped by default."""
        fn = getattr(self.collection, "to_dataframe", None)
        if not callable(fn):
            return self.summary()
        df = fn(*args, **kwargs)
        if not api:
            return df
        return wrap_frame(
            df,
            name=f"{self.name}_data",
            kind=kwargs.get("kind", None) or "edi.data",
            source=self.source,
            meta={"n_sites": self.n_sites, **self.meta},
        )

    def to_collection(self) -> Any:
        """Return the underlying collection."""
        return self.collection

    def update_meta(self, /, **kwargs: Any) -> "APISurvey":
        """Update metadata in-place and return ``self``."""
        self.meta.update(kwargs)
        return self

    def summary_text(self) -> str:
        """Return a compact survey display."""
        return self.__str__()

    def __str__(self) -> str:
        lines = [f"APISurvey: {self.name}"]
        lines.append(f"sites: {self.n_sites}")
        if self.stations:
            preview = ", ".join(self.stations[:8])
            if len(self.stations) > 8:
                preview += ", ..."
            lines.append(f"stations: {preview}")
        if self.source is not None:
            lines.append(f"source: {self.source}")
        errs = self.errors()
        if errs:
            lines.append(f"errors: {len(errs)}")
        return "\n".join(lines)
