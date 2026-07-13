# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
# ruff: noqa: E501
"""Multi-line EDI ingestion for :mod:`pycsamt.map`.

A map-view survey is usually made of several profile lines, each
stored as its own EDI folder.  The renderers in :mod:`pycsamt.map`
already understand the per-station ``line`` attribute, but the
normalizer :func:`pycsamt.map._core.ensure_map_data` only consumes a
single source.  This module lifts the multi-line grouping logic out
of the GUI layer so scripts and notebooks can build one combined
:class:`~pycsamt.map.MapData` from many lines.

Accepted ``source`` forms
-------------------------
mapping
    ``{line_name: dir | iterable_of_edis}`` — explicit line layout.
sequence of directories
    ``[dir_a, dir_b, ...]`` — one line per directory (named after it).
sequence of EDI-like items / files
    Treated as a single ``"line"``.
single directory
    Grouped according to *detect* (see :func:`load_lines`).
single file
    A one-station line named after the file stem.
"""

from __future__ import annotations

import re
from collections.abc import Iterable, Mapping
from dataclasses import replace
from pathlib import Path
from typing import Any

from ._core import MapData, ensure_map_data

__all__ = ["load_lines", "resolve_line_groups"]

_EDI_GLOBS = ("*.edi", "*.EDI")
_ID_PREFIX = re.compile(r"^([A-Za-z]*\d+)")


def load_lines(
    source: Any,
    *,
    detect: str = "folder",
    recursive: bool = True,
    line_map: dict[str, Iterable[str]] | None = None,
    verbose: int = 0,
) -> MapData:
    """Build one combined :class:`MapData` from multiple lines.

    Parameters
    ----------
    source :
        Mapping, sequence, directory, or file (see module docstring).
        A :class:`MapData` is returned unchanged.
    detect : {"folder", "auto", "flat"}, default "folder"
        Line-grouping rule applied only when *source* is a single
        directory:

        ``folder``
            One line per immediate sub-folder holding EDI files
            (a flat folder becomes one line named after it).
        ``auto``
            Group by the station-ID prefix, e.g. ``22-001`` -> ``L22``.
        ``flat``
            Treat the whole directory as a single line.
    recursive :
        Recurse into sub-directories when discovering EDI files and
        when calling :func:`ensure_map_data`.
    line_map :
        Optional ``line -> station names`` mapping forwarded to
        :func:`ensure_map_data` for sources without embedded line
        metadata.
    verbose :
        Verbosity forwarded to :func:`ensure_map_data`.

    Returns
    -------
    MapData
        Stations from every line, re-indexed and tagged with their
        resolved line name; profiles are rebuilt per line.
    """
    if isinstance(source, MapData):
        return source

    groups = resolve_line_groups(
        source,
        detect=detect,
        recursive=recursive,
    )

    stations: list[Any] = []
    edis: list[Any] = []
    index = 0
    for line, line_source in groups.items():
        if not _has_items(line_source):
            continue
        line_data = ensure_map_data(
            line_source,
            recursive=recursive,
            line_map=line_map,
            verbose=verbose,
        )
        for station in line_data.stations:
            stations.append(replace(station, line=line, index=index))
            index += 1
        edis.extend(line_data.iter_edis())

    metadata = {
        "n_stations": len(stations),
        "n_lines": len({s.line for s in stations}),
        "lines": tuple(groups),
    }
    return MapData(
        sites=tuple(edis),
        stations=tuple(stations),
        profiles=(),
        metadata=metadata,
    )


def resolve_line_groups(
    source: Any,
    *,
    detect: str = "folder",
    recursive: bool = True,
) -> dict[str, Any]:
    """Return a ``{line_name: source}`` mapping for *source*.

    The values are passed straight to :func:`ensure_map_data`, so they
    may be directories, lists of file paths, or EDI-like iterables.
    This mirrors the ``groups`` store used by the GUI loader, making
    the two paths share one grouping contract.
    """
    if isinstance(source, Mapping):
        return {str(k): v for k, v in source.items()}

    if isinstance(source, (str, Path)):
        path = Path(source)
        if path.is_file():
            return {path.stem: [str(path)]}
        if not path.exists():
            raise FileNotFoundError(f"EDI source path not found: {path}")
        return _detect_dir_groups(
            path,
            detect=detect,
            recursive=recursive,
        )

    if isinstance(source, Iterable) and not isinstance(
        source,
        (str, bytes),
    ):
        items = list(source)
        if items and all(_is_dir(it) for it in items):
            return {Path(str(it)).name: str(it) for it in items}
        return {"line": items}

    raise TypeError(
        f"Unsupported EDI source for load_lines: {type(source)!r}."
    )


# ── directory grouping helpers ─────────────────────────


def _detect_dir_groups(
    path: Path,
    *,
    detect: str,
    recursive: bool,
) -> dict[str, Any]:
    mode = (detect or "folder").lower()
    if mode == "auto":
        return _group_by_id_prefix(path, recursive)
    if mode == "flat":
        return {path.name: _discover_edis(path, recursive)}
    if mode == "folder":
        return _group_by_subfolder(path, recursive)
    raise ValueError(
        f"Unknown detect mode {detect!r}. "
        "Expected 'folder', 'auto', or 'flat'."
    )


def _group_by_subfolder(
    path: Path,
    recursive: bool,
) -> dict[str, list[str]]:
    groups: dict[str, list[str]] = {}
    for edi in _discover_edis(path, recursive):
        parent = Path(edi).parent
        line = parent.name if parent != path else path.name
        groups.setdefault(line, []).append(edi)
    return groups or {path.name: []}


def _group_by_id_prefix(
    path: Path,
    recursive: bool,
) -> dict[str, list[str]]:
    groups: dict[str, list[str]] = {}
    for edi in _discover_edis(path, recursive):
        stem = Path(edi).stem
        match = _ID_PREFIX.match(stem)
        prefix = match.group(1) if match else stem[:4]
        if prefix.isdigit():
            prefix = f"L{prefix}"
        groups.setdefault(prefix, []).append(edi)
    return groups or {path.name: []}


def _discover_edis(path: Path, recursive: bool) -> list[str]:
    found: list[str] = []
    for pattern in _EDI_GLOBS:
        matches = path.rglob(pattern) if recursive else path.glob(pattern)
        found.extend(str(p) for p in matches)
    # de-duplicate (case-insensitive globs can overlap) and sort
    return sorted(dict.fromkeys(found))


def _is_dir(item: Any) -> bool:
    try:
        return Path(str(item)).is_dir()
    except (TypeError, ValueError):
        return False


def _has_items(line_source: Any) -> bool:
    if isinstance(line_source, (str, Path)):
        return True
    if isinstance(line_source, Iterable):
        return bool(list(line_source))
    return line_source is not None
