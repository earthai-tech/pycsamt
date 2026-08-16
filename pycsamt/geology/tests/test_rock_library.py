# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for pycsamt.geology.rock_library — the built-in rock table itself."""

from __future__ import annotations

import pytest

from pycsamt.geology.rock_library import BUILTIN_ROCKS


def test_builtin_rocks_has_at_least_the_original_25_entries():
    assert len(BUILTIN_ROCKS) >= 25


def test_builtin_rocks_rows_are_well_formed():
    for row in BUILTIN_ROCKS:
        name, rho_min, rho_max, color, description, code, source = row
        assert isinstance(name, str) and name
        assert isinstance(rho_min, (int, float))
        assert isinstance(rho_max, (int, float))
        assert rho_min > 0
        assert rho_max > rho_min
        assert isinstance(color, str) and color.startswith("#")
        assert isinstance(description, str)
        assert isinstance(code, int) and code > 0
        assert isinstance(source, str) and source


def test_builtin_rocks_names_are_unique():
    names = [row[0] for row in BUILTIN_ROCKS]
    assert len(names) == len(set(names))


def test_builtin_rocks_codes_are_unique():
    codes = [row[5] for row in BUILTIN_ROCKS]
    assert len(codes) == len(set(codes))


def test_builtin_rocks_ranges_overlap_the_full_covered_span():
    # Design intent documented in the lithology guide: adjacent entries
    # should overlap or abut, so the covered span has no silent gaps a
    # user could fall into.
    by_min = sorted(BUILTIN_ROCKS, key=lambda row: row[1])
    running_max = by_min[0][2]
    for row in by_min[1:]:
        rho_min = row[1]
        assert rho_min <= running_max, (
            f"gap before {row[0]!r}: covered up to {running_max}, "
            f"next entry starts at {rho_min}"
        )
        running_max = max(running_max, row[2])


def test_rock_database_default_matches_builtin_rocks_length():
    from pycsamt.geology import RockDatabase

    db = RockDatabase.default()
    assert len(db) == len(BUILTIN_ROCKS)
