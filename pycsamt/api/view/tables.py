"""Opt-in APIFrame table wrappers for common pyCSAMT summaries."""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any

from .frame import api_frame, wrap_frame

__all__ = [
    "geology_dataframe",
    "geology_table",
    "quality_dataframe",
    "quality_table",
    "sites_summary",
]


@api_frame(
    name="sites_summary",
    kind="edi.summary",
    description="Per-site EDI frequency, tipper, and coordinate summary.",
)
def sites_summary(
    sites: Any,
    *,
    fields: Sequence[str] = (
        "station",
        "n_freq",
        "has_tipper",
        "period_min",
        "period_max",
        "lat",
        "lon",
    ),
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> Any:
    """Return a public APIFrame wrapper around ``emtools`` site summary."""
    from ...emtools.inspect import (
        sites_summary as _sites_summary,
    )

    return _sites_summary(
        sites,
        fields=tuple(fields),
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )


@api_frame(
    name="quality_dataframe",
    kind="metadata.quality",
    description="Per-station data quality assessment table.",
)
def quality_dataframe(sites: Any) -> Any:
    """Return a public APIFrame wrapper around metadata quality rows."""
    from ...metadata.quality import (
        quality_dataframe as _quality_dataframe,
    )

    return _quality_dataframe(sites)


def quality_table(sites: Any) -> Any:
    """Alias for :func:`quality_dataframe`."""
    return quality_dataframe(sites)


def geology_dataframe(
    catalog: Any = None,
    *,
    name: str = "geology_catalog",
) -> Any:
    """Return a public APIFrame for a geology catalog.

    Parameters
    ----------
    catalog : object, optional
        Object exposing ``to_dataframe()``. If omitted, the package-level
        built-in geology catalog is used.
    name : str, default ``"geology_catalog"``
        Display name attached to the returned APIFrame.
    """
    if catalog is None:
        from ...metadata.geology import CATALOG

        catalog = CATALOG
    df = catalog.to_dataframe()
    return wrap_frame(
        df,
        name=name,
        kind="metadata.geology",
        source=getattr(catalog, "__class__", type(catalog)).__name__,
        description="Resistivity and lithology ranges by formation.",
    )


def geology_table(catalog: Any = None) -> Any:
    """Alias for :func:`geology_dataframe`."""
    return geology_dataframe(catalog)
