# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Utility to find the most-recently modified MARE2DEM output file.

Port of ``m2d_getMostRecent.m``.

Useful when iterating over inversion runs to load the latest model
without hardcoding iteration numbers.
"""

from __future__ import annotations

from pathlib import Path

__all__ = ["get_most_recent"]


def get_most_recent(
    file_or_keyword: str | Path,
    pattern: str = "*.resistivity",
    *,
    search_dir: str | Path = ".",
) -> Path | None:
    """Return the most recently modified file matching *pattern*.

    Port of ``m2d_getMostRecent.m``.

    Parameters
    ----------
    file_or_keyword : str or path-like
        Either a literal file path, or one of the special keywords
        ``'lastiter'``, ``'last'``, or ``'newest'``.  When a keyword is
        given the function scans *search_dir* for files matching
        *pattern* and returns the most recently modified one.
    pattern : str, default "``*.resistivity``"
        Glob pattern used when a keyword is given.
    search_dir : path-like, default "."
        Directory to search when a keyword is given.

    Returns
    -------
    pathlib.Path or None
        Resolved file path, or ``None`` when the keyword is used but no
        matching files are found.

    Examples
    --------
    Load the latest inversion model:

    >>> from pycsamt.models.mare2dem.iotools.most_recent import get_most_recent
    >>> path = get_most_recent("newest", "*.resistivity", search_dir="./run")
    >>> path
    PosixPath('run/mare2dem.0020.resistivity')

    Use a literal path (pass-through):

    >>> get_most_recent("mare2dem.0010.resistivity")
    PosixPath('mare2dem.0010.resistivity')
    """
    keyword = str(file_or_keyword).strip().lower()
    if keyword not in ("lastiter", "last", "newest"):
        return Path(file_or_keyword)

    search = Path(search_dir)
    candidates = [
        p for p in search.glob(pattern) if not p.name.startswith(".")
    ]
    if not candidates:
        return None

    return max(candidates, key=lambda p: p.stat().st_mtime)
