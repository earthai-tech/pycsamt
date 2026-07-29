# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for DataController (Phase 2) — no Qt required.

Uses real EDI fixtures already checked into the repo (fast subsets of
``data/MT/kap03lmt_edis`` and ``data/AMT/WILLY_DATA``) rather than
inventing synthetic EDI content, per project convention.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

import pycsamt.site.base as site_base
from pycsamt.app.desktop.controllers.data_controller import DataController
from pycsamt.seg.edi import EDIFile
from pycsamt.site.base import Sites

REPO_ROOT = Path(__file__).resolve().parents[3]
KAP03_DIR = REPO_ROOT / "data" / "MT" / "kap03lmt_edis"
WILLY_DIR = REPO_ROOT / "data" / "AMT" / "WILLY_DATA" / "L18PLT"


@pytest.fixture
def kap03_paths() -> list[Path]:
    """Three real KAP TIPPER EDIs — no LAT/LONG in their headers."""
    paths = sorted(KAP03_DIR.glob("*.edi"))[:3]
    assert len(paths) == 3
    return paths


@pytest.fixture
def willy_paths() -> list[Path]:
    """Three real WILLY AMT EDIs — carry real LAT/LONG/ELEV."""
    paths = sorted(WILLY_DIR.glob("*.edi"))[:3]
    assert len(paths) == 3
    return paths


# ── Construction ─────────────────────────────────────────────────────────


def test_creates_with_no_callback():
    ctrl = DataController()
    assert ctrl.is_loaded() is False
    assert ctrl.sites is None
    assert ctrl.dataframe.empty
    assert list(ctrl.dataframe.columns) == DataController.STATION_COLUMNS


def test_station_columns_constant():
    assert DataController.STATION_COLUMNS == [
        "ID",
        "Line",
        "Latitude",
        "Longitude",
        "Elevation",
        "N_freq",
        "Tipper",
    ]


def test_creates_with_progress_callback_stored():
    calls = []
    cb = calls.append
    ctrl = DataController(progress_callback=cb)
    assert ctrl._progress_cb is cb


# ── load() ───────────────────────────────────────────────────────────────


def test_load_returns_sites_instance(kap03_paths):
    ctrl = DataController()
    result = ctrl.load(kap03_paths)
    assert isinstance(result, Sites)
    assert ctrl.sites is result


def test_load_sets_is_loaded_true(kap03_paths):
    ctrl = DataController()
    assert ctrl.is_loaded() is False
    ctrl.load(kap03_paths)
    assert ctrl.is_loaded() is True


def test_load_accepts_str_paths(kap03_paths):
    ctrl = DataController()
    ctrl.load([str(p) for p in kap03_paths])
    assert len(ctrl.dataframe) == 3


def test_load_builds_dataframe_with_expected_rows(kap03_paths):
    ctrl = DataController()
    ctrl.load(kap03_paths)
    df = ctrl.dataframe
    assert list(df.columns) == DataController.STATION_COLUMNS
    assert len(df) == 3
    assert set(df["ID"]) == {p.stem for p in kap03_paths}
    # kap03lmt EDIs are all TIPPER surveys.
    assert df["Tipper"].all()
    assert (df["N_freq"] > 0).all()


def test_load_empty_paths_list_does_not_raise():
    ctrl = DataController()
    sites = ctrl.load([])
    assert isinstance(sites, Sites)
    assert ctrl.dataframe.empty


def test_load_path_to_line_mapping(kap03_paths):
    ctrl = DataController()
    path_to_line = {
        str(kap03_paths[0]): "L1",
        str(kap03_paths[1]): "L2",
    }
    ctrl.load(kap03_paths, path_to_line=path_to_line)
    df = ctrl.dataframe.set_index("ID")
    assert df.loc[kap03_paths[0].stem, "Line"] == "L1"
    assert df.loc[kap03_paths[1].stem, "Line"] == "L2"
    # Unmapped station falls back to the em-dash placeholder.
    assert df.loc[kap03_paths[2].stem, "Line"] == "—"


def test_load_progress_callback_reaches_100(kap03_paths):
    calls = []
    ctrl = DataController(progress_callback=calls.append)
    ctrl.load(kap03_paths)
    # 3 files: 30%, 60%, 90% during the loop, then 100% at the end.
    assert calls == [30, 60, 90, 100]


def test_load_without_callback_does_not_raise(kap03_paths):
    ctrl = DataController(progress_callback=None)
    ctrl.load(kap03_paths)  # must not raise
    assert ctrl.is_loaded()


def test_load_real_latlon_populated(willy_paths):
    ctrl = DataController()
    ctrl.load(willy_paths)
    df = ctrl.dataframe
    assert df["Latitude"].notna().all()
    assert df["Longitude"].notna().all()
    assert df["Elevation"].notna().all()


# ── _build_dataframe fallback branches ──────────────────────────────────


def test_build_dataframe_sites_none_returns_empty_frame():
    ctrl = DataController()
    ctrl._sites = None
    df = ctrl._build_dataframe()
    assert df.empty
    assert list(df.columns) == DataController.STATION_COLUMNS


def test_build_dataframe_converts_via_to_sites_when_no_as_list(kap03_paths):
    # A plain list of EDIFile objects has no .as_list(), so
    # _build_dataframe() must route it through site_base.to_sites().
    edis = [EDIFile(p) for p in kap03_paths[:2]]
    ctrl = DataController()
    ctrl._sites = edis
    assert not hasattr(ctrl._sites, "as_list")
    df = ctrl._build_dataframe()
    assert isinstance(ctrl._sites, Sites)  # converted in place
    assert len(df) == 2


def test_build_dataframe_to_sites_failure_returns_empty_frame(monkeypatch, kap03_paths):
    def boom(*_a, **_k):
        raise RuntimeError("cannot coerce")

    monkeypatch.setattr(site_base, "to_sites", boom)

    ctrl = DataController()
    ctrl._sites = object()  # no as_list -> forces the to_sites() path
    df = ctrl._build_dataframe()
    assert df.empty
    assert list(df.columns) == DataController.STATION_COLUMNS


def test_build_dataframe_summary_fallback_hits_location_attribute_bug(
    monkeypatch, willy_paths
):
    """
    Document a real bug (not fixed here, per task constraints):

    When ``Sites.get()``/``Site.summary()`` fails, the fallback in
    ``_build_dataframe`` reads ``edi.get_section("head").Location``
    and accesses ``loc.lat`` / ``loc.lon`` / ``loc.elev``. But
    ``pycsamt.loc.Location`` exposes ``latitude`` / ``longitude`` /
    ``elevation`` — there is no ``.lat``/``.lon``/``.elev`` attribute.
    So whenever the primary summary() path fails but the EDI header
    is readable, this "fallback" itself raises AttributeError instead
    of degrading gracefully, and the exception is not caught anywhere
    up the call chain (unlike the inner get_section()/Location lookup,
    which *is* guarded).
    """
    monkeypatch.setattr(Sites, "get", lambda self, name: None)

    ctrl = DataController()
    with pytest.raises(AttributeError, match="lat"):
        ctrl.load(willy_paths[:1])


def test_build_dataframe_double_fallback_yields_nan_row(monkeypatch, willy_paths):
    """
    When *both* site.summary() and the Location lookup fail (header
    section present but its cached ``.Location`` missing), the code
    degrades to a NaN-filled row instead of crashing, because that
    inner failure *is* guarded by its own try/except.
    """
    edi = EDIFile(willy_paths[0])
    head = edi.get_section("head")
    del head.__dict__["Location"]

    sites = Sites([edi])
    monkeypatch.setattr(Sites, "get", lambda self, name: None)

    ctrl = DataController()
    ctrl._sites = sites
    df = ctrl._build_dataframe()

    assert len(df) == 1
    row = df.iloc[0]
    assert row["ID"] == edi.station
    assert pd.isna(row["Latitude"])
    assert pd.isna(row["Longitude"])
    assert pd.isna(row["Elevation"])
    assert row["N_freq"] == edi.n_freq
    assert bool(row["Tipper"]) == bool(edi.has_tipper)


# ── dataframe / sites properties ────────────────────────────────────────


def test_dataframe_is_empty_before_load():
    ctrl = DataController()
    assert ctrl.dataframe.empty


def test_dataframe_returns_independent_copy(kap03_paths):
    ctrl = DataController()
    ctrl.load(kap03_paths)
    df1 = ctrl.dataframe
    df1.loc[0, "ID"] = "MUTATED"
    df2 = ctrl.dataframe
    assert df2.loc[0, "ID"] != "MUTATED"


def test_sites_property_before_and_after_load(kap03_paths):
    ctrl = DataController()
    assert ctrl.sites is None
    ctrl.load(kap03_paths)
    assert ctrl.sites is not None
    assert isinstance(ctrl.sites, Sites)


# ── filter_by_ids ────────────────────────────────────────────────────────


def test_filter_by_ids_before_load_returns_empty():
    ctrl = DataController()
    result = ctrl.filter_by_ids(["kap103"])
    assert result.empty
    assert list(result.columns) == DataController.STATION_COLUMNS


def test_filter_by_ids_returns_matching_rows(kap03_paths):
    ctrl = DataController()
    ctrl.load(kap03_paths)
    wanted = [kap03_paths[0].stem, kap03_paths[2].stem]
    result = ctrl.filter_by_ids(wanted)
    assert set(result["ID"]) == set(wanted)
    assert len(result) == 2


def test_filter_by_ids_no_match_returns_empty(kap03_paths):
    ctrl = DataController()
    ctrl.load(kap03_paths)
    result = ctrl.filter_by_ids(["does-not-exist"])
    assert result.empty


def test_filter_by_ids_resets_index(kap03_paths):
    ctrl = DataController()
    ctrl.load(kap03_paths)
    result = ctrl.filter_by_ids([kap03_paths[2].stem])
    assert list(result.index) == [0]


# ── station_coords ───────────────────────────────────────────────────────


def test_station_coords_columns(willy_paths):
    ctrl = DataController()
    ctrl.load(willy_paths)
    coords = ctrl.station_coords()
    assert list(coords.columns) == ["ID", "Latitude", "Longitude"]


def test_station_coords_drops_nan_rows(kap03_paths):
    # kap03lmt EDIs carry no LAT/LONG -> summary() yields NaN lat/lon,
    # so station_coords() must drop every row.
    ctrl = DataController()
    ctrl.load(kap03_paths)
    assert ctrl.dataframe["Latitude"].isna().all()
    coords = ctrl.station_coords()
    assert coords.empty


def test_station_coords_keeps_rows_with_real_latlon(willy_paths):
    ctrl = DataController()
    ctrl.load(willy_paths)
    coords = ctrl.station_coords()
    assert len(coords) == 3
    assert set(coords["ID"]) == {p.stem for p in willy_paths}


def test_station_coords_mixed_nan_and_valid_rows(kap03_paths, willy_paths):
    ctrl = DataController()
    ctrl.load([kap03_paths[0], willy_paths[0]])
    coords = ctrl.station_coords()
    # Only the WILLY station (has real lat/lon) should survive.
    assert list(coords["ID"]) == [willy_paths[0].stem]


# ── is_loaded ─────────────────────────────────────────────────────────────


def test_is_loaded_false_initially():
    assert DataController().is_loaded() is False


def test_is_loaded_true_after_load(kap03_paths):
    ctrl = DataController()
    ctrl.load(kap03_paths)
    assert ctrl.is_loaded() is True
