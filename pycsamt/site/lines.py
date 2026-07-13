# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
site/lines.py — survey-line detection and naming helpers.

Core functions used by both the web app and the desktop GUI so that
line-detection logic is not duplicated in the UI layer.
"""

from __future__ import annotations

import re
from collections.abc import Sequence
from pathlib import Path

__all__ = [
    "detect_lines_from_station_ids",
    "detect_lines_from_filenames",
    "apply_line_renames",
    "resolve_display_name",
    "pick_representative_stations",
]

_SEP = re.compile(r"[-_]")


def detect_lines_from_station_ids(
    station_ids: Sequence[str],
    *,
    min_stations_per_line: int = 1,
) -> dict[str, list[str]]:
    """
    Infer survey-line groupings from station ID strings.

    Strategy
    --------
    1. Split each ID on the first ``-`` or ``_`` separator.
    2. The leading token is the candidate line prefix.
    3. If the prefix is **purely numeric** → prepend ``'L'``
       (e.g. ``'18-001A'`` → line ``'L18'``).
    4. Non-numeric prefix used as-is, upper-cased
       (e.g. ``'WL07B'`` → line ``'WL07B'``).
    5. IDs with no separator → the whole ID becomes the key.
    6. If the resulting number of groups equals the number of stations
       (every ID ended up alone), the split is considered non-informative
       and all stations are returned as a single ``'auto'`` group.

    Parameters
    ----------
    station_ids : sequence of str
        Station identifiers, typically the ``ID`` column from
        ``station_records``.
    min_stations_per_line : int, default 1
        Groups smaller than this are merged into ``'unassigned'``.

    Returns
    -------
    dict[str, list[str]]
        Maps each inferred line name to the list of station IDs that
        belong to it, preserving insertion order.

    Examples
    --------
    >>> ids = ["18-001A", "18-003U", "22-345A", "22-007B"]
    >>> detect_lines_from_station_ids(ids)
    {'L18': ['18-001A', '18-003U'], 'L22': ['22-345A', '22-007B']}

    >>> ids = ["WL07A", "WL07B", "WL09A"]
    >>> detect_lines_from_station_ids(ids)
    {'WL07A': ['WL07A'], 'WL07B': ['WL07B'], 'WL09A': ['WL09A']}
    """
    groups: dict[str, list[str]] = {}
    for sid in station_ids:
        token = _SEP.split(str(sid), maxsplit=1)[0].strip()
        if not token:
            key = "unassigned"
        elif token.isdigit():
            key = f"L{token}"
        else:
            key = token.upper()
        groups.setdefault(key, []).append(str(sid))

    # Trivial split: every station in its own group → no useful line structure
    if len(groups) == len(station_ids) and len(station_ids) > 1:
        return {"auto": list(station_ids)}

    # Merge small groups
    if min_stations_per_line > 1:
        tiny = [
            k for k, v in groups.items() if len(v) < min_stations_per_line
        ]
        if tiny:
            unassigned = groups.setdefault("unassigned", [])
            for k in tiny:
                unassigned.extend(groups.pop(k))

    return groups


def detect_lines_from_filenames(
    file_paths: Sequence[str],
    *,
    use_parent: bool = False,
) -> dict[str, list[str]]:
    """
    Infer line groupings from EDI file paths.

    Parameters
    ----------
    file_paths : sequence of str
        Absolute or relative paths to EDI files.
    use_parent : bool, default False
        When True, use the immediate parent folder name as the line key
        (same as the *folder mode* in the web UI). When False, the file
        stem is parsed with the same prefix rule as
        :func:`detect_lines_from_station_ids`.

    Returns
    -------
    dict[str, list[str]]
        Maps each inferred line name to the file paths that belong to it.

    Examples
    --------
    >>> paths = ["/data/line18/S001.edi", "/data/line22/S002.edi"]
    >>> detect_lines_from_filenames(paths, use_parent=True)
    {'line18': ['/data/line18/S001.edi'], 'line22': ['/data/line22/S002.edi']}
    """
    groups: dict[str, list[str]] = {}
    for fp in file_paths:
        p = Path(fp)
        if use_parent:
            key = p.parent.name or "unassigned"
        else:
            stem = p.stem
            token = _SEP.split(stem, maxsplit=1)[0].strip()
            if not token:
                key = "unassigned"
            elif token.isdigit():
                key = f"L{token}"
            else:
                key = token.upper()
        groups.setdefault(key, []).append(str(fp))
    return groups


def apply_line_renames(
    station_records: list[dict],
    renames: dict[str, str],
) -> list[dict]:
    """
    Apply a ``{old_name: new_name}`` rename map to a list of station records.

    Parameters
    ----------
    station_records : list of dict
        Each dict must have a ``'Line'`` key.
    renames : dict[str, str]
        Mapping from current line name to desired display name.
        Entries where old == new are skipped.

    Returns
    -------
    list of dict
        A *new* list with the ``'Line'`` fields updated.  Records whose
        line is not in *renames* are returned unchanged.

    Examples
    --------
    >>> recs = [{"ID": "S1", "Line": "L18"}, {"ID": "S2", "Line": "L22"}]
    >>> apply_line_renames(recs, {"L18": "South Profile"})
    [{'ID': 'S1', 'Line': 'South Profile'}, {'ID': 'S2', 'Line': 'L22'}]
    """
    out = []
    for rec in station_records:
        new_rec = dict(rec)
        old = new_rec.get("Line", "")
        if old in renames and renames[old] != old:
            new_rec["Line"] = renames[old]
        out.append(new_rec)
    return out


def resolve_display_name(
    line_name: str,
    renames: dict[str, str] | None = None,
) -> str:
    """Return the display name for *line_name*, consulting *renames* first."""
    if renames and line_name in renames:
        return renames[line_name]
    return line_name


def pick_representative_stations(
    station_ids: Sequence[str],
    k: int,
) -> list[str]:
    """Return up to *k* evenly spaced station IDs from a sorted line.

    Used to build a compact ``profiles={line: [stations]}`` grouping for
    multi-profile ellipse-strip plots (see
    :func:`~pycsamt.emtools.tensor.plot_phase_tensor_strip_grid`) without
    drawing every single station in a long line.

    Parameters
    ----------
    station_ids : sequence of str
        Station IDs belonging to one survey line.
    k : int
        Maximum number of stations to keep.  The full (sorted) list is
        returned unchanged when it already has ``k`` or fewer entries.

    Examples
    --------
    >>> pick_representative_stations(
    ...     ["18-001A", "18-010U", "18-019U", "18-025A", "18-005U"], 3)
    ['18-001A', '18-010U', '18-025A']
    """
    names = sorted(station_ids)
    if len(names) <= k:
        return names
    step = (len(names) - 1) / (k - 1) if k > 1 else 0
    idx = sorted({round(i * step) for i in range(k)})
    return [names[i] for i in idx]
