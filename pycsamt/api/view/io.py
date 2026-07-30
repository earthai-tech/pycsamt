"""Public readers for EDI files and site collections."""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path
from typing import Any

from .progress import iter_progress
from .survey import APISurvey

__all__ = ["read_edi", "read_edis", "read_sites"]


def read_edi(
    path: Any,
    *,
    verbose: int = 0,
    **kwargs: Any,
) -> Any:
    """Read one EDI file and return an :class:`EDIFile`."""
    from ...seg.edi import EDIFile

    return EDIFile(path, verbose=verbose, **kwargs)


def read_edis(
    sources: Any | Sequence[Any],
    *,
    recursive: bool = True,
    strict: bool = False,
    on_dup: str = "replace",
    progress: bool | str = "auto",
    leave: bool = False,
    verbose: int = 0,
) -> APISurvey:
    """Read many EDI files and return a public survey view."""
    from ...seg.cbase import CoreParser
    from ...seg.collection import EDICollection

    parser = CoreParser(
        recursive=recursive,
        strict=strict,
        on_dup=on_dup,
        verbose=verbose,
    )
    paths = list(
        parser._iter_edi_files(  # noqa: SLF001 - public facade over parser
            sources
            if isinstance(sources, list)
            else list(parser._iter_paths(sources))
        )
    )
    items: list[Any] = []
    parser.results.clear()
    parser._errors = getattr(parser, "_errors", [])  # noqa: SLF001
    for path in iter_progress(
        paths,
        enabled=progress,
        desc="Reading EDI",
        leave=leave,
        unit="file",
        total=len(paths),
    ):
        res = parser._read_one(Path(path))  # noqa: SLF001
        parser.results.append(res)
        if res.edi is not None:
            items.append(res.edi)

    by_station: dict[str, Any] = {}
    for ed in items:
        sid = getattr(ed, "station", None)
        if not sid:
            sid = parser._fast_station(ed.path) if ed.path else None
        sid = sid or (str(ed.path) if ed.path else "-")
        if sid in by_station and on_dup == "keep":
            continue
        by_station[sid] = ed

    collection = EDICollection(items=by_station.values(), verbose=verbose)
    survey = APISurvey(
        collection,
        name="edi_survey",
        source=sources,
        parser=parser,
        meta={
            "recursive": recursive,
            "strict": strict,
            "on_dup": on_dup,
        },
    )
    return survey


def read_sites(
    sources: Any | Sequence[Any],
    **kwargs: Any,
) -> APISurvey:
    """Alias for :func:`read_edis`."""
    return read_edis(sources, **kwargs)
