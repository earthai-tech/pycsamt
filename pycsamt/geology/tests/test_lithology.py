# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for pycsamt.geology.lithology — RockEntry, RockDatabase, Layer,
StratigraphicLog."""

from __future__ import annotations

import numpy as np
import pytest

from pycsamt.geology.lithology import (
    Layer,
    RockDatabase,
    RockEntry,
    StratigraphicLog,
)

# ─────────────────────────────────────────────────────────────────────────────
# RockEntry
# ─────────────────────────────────────────────────────────────────────────────


def test_rock_entry_rho_mid_is_geometric_mean():
    e = RockEntry(name="Sandstone", rho_min=50.0, rho_max=5000.0)
    assert e.rho_mid == pytest.approx(np.sqrt(50.0 * 5000.0))
    assert e.log_rho_mid == pytest.approx(np.log10(e.rho_mid))


def test_rock_entry_contains_is_inclusive_both_ends():
    e = RockEntry(name="Clay", rho_min=1.0, rho_max=20.0)
    assert e.contains(1.0) is True
    assert e.contains(20.0) is True
    assert e.contains(0.999) is False
    assert e.contains(20.001) is False


# ─────────────────────────────────────────────────────────────────────────────
# RockDatabase — construction
# ─────────────────────────────────────────────────────────────────────────────


def test_default_database_matches_builtin_rock_count():
    from pycsamt.geology.rock_library import BUILTIN_ROCKS

    db = RockDatabase.default()
    assert len(db) == len(BUILTIN_ROCKS)
    assert len(db.entries) == len(BUILTIN_ROCKS)
    assert len(db) >= 25  # never shrinks below the original built-in set


def test_entries_is_read_only_tuple_and_does_not_alias_internal_list():
    db = RockDatabase.default()
    entries = db.entries
    n = len(db)
    assert isinstance(entries, tuple)
    assert all(isinstance(e, RockEntry) for e in entries)
    # Mutating the returned tuple's backing list must not be possible /
    # must not affect the database.
    with pytest.raises(TypeError):
        entries[0] = entries[1]  # type: ignore[index]
    assert len(db) == n


def test_from_csv_round_trip(tmp_path):
    p = tmp_path / "rocks.csv"
    p.write_text(
        "name,rho_min,rho_max,color,description,code\n"
        "Peat,1,10,#000000,Organic soil,1\n"
        "Till,20,300,#888888,Glacial till,2\n"
    )
    db = RockDatabase.from_csv(p)
    assert len(db) == 2
    names = [e.name for e in db.entries]
    assert names == ["Peat", "Till"]
    assert db.entries[0].rho_min == 1.0
    assert db.entries[0].rho_max == 10.0
    assert db.entries[0].code == 1


def test_from_csv_defaults_when_optional_columns_missing(tmp_path):
    p = tmp_path / "rocks_min.csv"
    p.write_text("name,rho_min,rho_max\nPeat,1,10\n")
    db = RockDatabase.from_csv(p)
    assert db.entries[0].color == "#AAAAAA"
    assert db.entries[0].description == ""
    assert db.entries[0].code == 1  # 1-based fallback (i + 1)


# ─────────────────────────────────────────────────────────────────────────────
# RockDatabase — classification
# ─────────────────────────────────────────────────────────────────────────────


def test_classify_nearest_matches_default_database():
    db = RockDatabase.default()
    assert db.classify(250.0).name == "Granite (weathered)"
    assert db.classify(180.0).name == "Granite (weathered)"


def test_classify_overlap_returns_first_containing_entry_in_db_order():
    db = RockDatabase.default()
    result = db.classify(250.0, method="overlap")
    assert result.contains(250.0)
    # First entry (in insertion order) whose range brackets 250.0.
    containing = [e.name for e in db.entries if e.contains(250.0)]
    assert result.name == containing[0]


def test_classify_nearest_and_overlap_can_disagree():
    db = RockDatabase.default()
    nearest = db.classify(250.0, method="nearest")
    overlap = db.classify(250.0, method="overlap")
    assert nearest.name != overlap.name


def test_classify_nearest_extrapolates_beyond_database_coverage():
    # Below the lowest range (Sulfide ore body: 0.001-0.1 ohm.m) or above
    # the highest (Air / void: 1e6-1e12 ohm.m), "nearest" still returns an
    # entry, but that entry does not literally contain the query value.
    db = RockDatabase.default()
    below = db.classify(1e-6)
    above = db.classify(1e13)
    assert below.name == "Sulfide ore body"
    assert below.contains(1e-6) is False
    assert above.name == "Air / void"
    assert above.contains(1e13) is False


def test_classify_nan_or_nonpositive_returns_first_entry():
    db = RockDatabase.default()
    first = db.entries[0]
    assert db.classify(float("nan")).name == first.name
    assert db.classify(0.0).name == first.name
    assert db.classify(-5.0).name == first.name


def test_classify_column_uses_nearest_method_only():
    db = RockDatabase.default()
    rho_log10 = np.log10(np.array([250.0, 180.0]))
    result = db.classify_column(rho_log10)
    assert [e.name for e in result] == [
        db.classify(250.0, method="nearest").name,
        db.classify(180.0, method="nearest").name,
    ]


def test_database_repr():
    db = RockDatabase.default()
    assert repr(db) == f"RockDatabase({len(db)} entries)"


# ─────────────────────────────────────────────────────────────────────────────
# Layer
# ─────────────────────────────────────────────────────────────────────────────


def test_layer_thickness_and_rho_ohm_m():
    ly = Layer(top=10.0, bottom=25.0, rho_log10=2.0, lithology="Sandstone")
    assert ly.thickness == 15.0
    assert ly.rho_ohm_m == pytest.approx(100.0)


# ─────────────────────────────────────────────────────────────────────────────
# StratigraphicLog.from_column
# ─────────────────────────────────────────────────────────────────────────────


def _column():
    z = np.array([5.0, 15.0, 30.0, 55.0, 90.0])
    rho_ohm_m = np.array([420.0, 95.0, 42.0, 190.0, 1500.0])
    return z, np.log10(rho_ohm_m)


def test_from_column_default_database():
    z, rho_log10 = _column()
    log = StratigraphicLog.from_column("S02", 500.0, z, rho_log10)
    assert log.station_name == "S02"
    assert log.station_x == 500.0
    assert len(log.layers) >= 1
    assert log.layers[0].top == 0.0
    assert log.layers[-1].bottom == pytest.approx(z[-1] + (z[-1] - z[-2]) / 2)


def test_from_column_merges_adjacent_same_lithology_within_tolerance():
    z = np.array([10.0, 20.0, 30.0])
    # Same lithology bucket (both classify to the same entry) and close
    # in log-rho, so they must merge into a single layer.
    rho_log10 = np.array([1.7, 1.72, 1.71])
    log = StratigraphicLog.from_column(
        "S00", 0.0, z, rho_log10, merge_tolerance=0.2
    )
    assert len(log.layers) == 1


def test_from_column_splits_when_merge_tolerance_exceeded():
    z = np.array([10.0, 20.0])
    rho_log10 = np.array([1.0, 1.0])
    db = RockDatabase.default()
    entry = db.classify(10.0**1.0)
    # Force the same-lithology pair apart with an artificially tiny
    # tolerance; even identical values must not merge.
    log = StratigraphicLog.from_column(
        "S00", 0.0, z, rho_log10, db=db, merge_tolerance=-1.0
    )
    assert len(log.layers) == 2
    assert all(ly.lithology == entry.name for ly in log.layers)


def test_from_column_confidence_is_one_when_uniform():
    z = np.array([10.0, 20.0, 30.0])
    rho_log10 = np.array([1.7, 1.72, 1.71])
    log = StratigraphicLog.from_column("S00", 0.0, z, rho_log10)
    assert log.layers[0].confidence == pytest.approx(1.0)


def test_to_dict_roundtrip_shape():
    z, rho_log10 = _column()
    log = StratigraphicLog.from_column("S02", 500.0, z, rho_log10)
    d = log.to_dict()
    assert d["station_name"] == "S02"
    assert d["station_x"] == 500.0
    assert len(d["layers"]) == len(log.layers)
    assert set(d["layers"][0]) == {"top", "bottom", "rho_log10", "lithology"}


def test_to_dataframe_columns():
    pd = pytest.importorskip("pandas")
    z, rho_log10 = _column()
    log = StratigraphicLog.from_column("S02", 500.0, z, rho_log10)
    df = log.to_dataframe()
    assert isinstance(df, pd.DataFrame)
    assert list(df.columns) == [
        "station",
        "x",
        "top",
        "bottom",
        "thickness",
        "rho_log10",
        "rho_ohm_m",
        "lithology",
        "color",
        "confidence",
    ]
    assert len(df) == len(log.layers)


def test_log_repr():
    z, rho_log10 = _column()
    log = StratigraphicLog.from_column("S02", 500.0, z, rho_log10)
    assert repr(log).startswith("StratigraphicLog('S02', x=500.0 m,")
