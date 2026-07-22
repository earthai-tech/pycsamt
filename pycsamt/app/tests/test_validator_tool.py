# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Tests for pycsamt.app.desktop.tools.validator_tool (EDIValidatorDialog).

Real data
---------
data/AMT/WILLY_DATA/L18PLT/ — WILLY AMT EDIs (first 5 files used for speed).

Strategy
--------
* ``_check_site`` is exercised directly with small duck-typed mock objects
  so every branch (pass / warn / fail, coordinate fallbacks, station-name
  fallbacks, Z-tensor shape edge cases) is deterministic and fast.
* It is also exercised against real unwrapped EDI objects for an
  integration sanity check.
* ``_swap_excluded`` is tested directly with a handful of index pairs.
* The full ``EDIValidatorDialog`` is driven end to end (construction,
  ``_run_checks``/``_fill_table``, item editing, selection, rename,
  delete, exclude, reorder, context menu, apply, export) against a small
  real-data survey, with modal dialogs neutralised via monkeypatching.
"""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("PySide6")

from PySide6.QtCore import QItemSelectionModel, QPoint
from PySide6.QtWidgets import QFileDialog, QInputDialog, QMenu, QMessageBox

from pycsamt.app.desktop.tools import validator_tool as _validator_tool_mod
from pycsamt.app.desktop.tools.validator_tool import (
    _COL_STATION,
    _COLS,
    _EXCLUDED,
    _GREEN,
    _RED,
    _YELLOW,
    _check_site,
    _icon,
    _swap_excluded,
    EDIValidatorDialog,
)

# ── Paths ─────────────────────────────────────────────────────────────────────

_ROOT = Path(__file__).parents[3]  # pycsamt/
_WILLY_DIR = _ROOT / "data" / "AMT" / "WILLY_DATA" / "L18PLT"
_HAS_WILLY = _WILLY_DIR.exists() and any(_WILLY_DIR.glob("*.edi"))


@pytest.fixture(scope="session")
def willy_edi_paths():
    if not _HAS_WILLY:
        pytest.skip("WILLY L18PLT data not available")
    return [str(p) for p in sorted(_WILLY_DIR.glob("*.edi"))[:5]]


@pytest.fixture(scope="session")
def willy_sites(willy_edi_paths):
    """5-station real Sites collection, loaded once per session."""
    pytest.importorskip("pycsamt.emtools")
    from pycsamt.emtools import ensure_sites

    return ensure_sites(willy_edi_paths)


@pytest.fixture()
def willy_items(willy_edi_paths):
    """Fresh list of *unwrapped* raw EDI objects (not shared across tests)."""
    from pycsamt.emtools import ensure_sites
    from pycsamt.emtools._core import _unwrap

    sites = ensure_sites(willy_edi_paths)
    return [_unwrap(s) for s in sites]


# ── Mock duck-typed site objects for _check_site ───────────────────────────────


class _Z:
    def __init__(self, z=None, z_err=None, freq=None):
        self.z = z
        self.z_err = z_err
        self.freq = freq


class _Site:
    """Minimal duck-typed stand-in for an EDI/Site object."""

    def __init__(self, **kw):
        for k, v in kw.items():
            setattr(self, k, v)


def _full_z(n=5):
    z = np.ones((n, 2, 2), dtype=complex)
    z_err = np.ones((n, 2, 2)) * 0.1
    freq = np.linspace(1.0, 100.0, n)
    return z, z_err, freq


# ── _icon ────────────────────────────────────────────────────────────────────


class TestIcon:
    def test_existing_icon_returns_nonnull(self):
        icon = _icon("recompute")
        assert not icon.isNull()

    def test_missing_icon_returns_empty(self):
        icon = _icon("definitely_not_a_real_icon_xyz")
        assert icon.isNull()


# ── _swap_excluded ───────────────────────────────────────────────────────────


class TestSwapExcluded:
    def test_neither_excluded_unchanged(self):
        assert _swap_excluded(set(), 0, 1) == set()

    def test_both_excluded_unchanged(self):
        assert _swap_excluded({0, 1}, 0, 1) == {0, 1}

    def test_only_a_excluded_moves_to_b(self):
        assert _swap_excluded({2}, 2, 3) == {3}

    def test_only_b_excluded_moves_to_a(self):
        assert _swap_excluded({3}, 2, 3) == {2}

    def test_unrelated_indices_untouched(self):
        assert _swap_excluded({0, 5}, 2, 3) == {0, 5}

    def test_does_not_mutate_input(self):
        original = {2}
        _swap_excluded(original, 2, 3)
        assert original == {2}


# ── _check_site — mock-driven branch coverage ──────────────────────────────────


class TestCheckSiteMocks:
    def test_full_pass(self):
        z, z_err, freq = _full_z()
        ed = _Site(station="ST1", lat=1.0, lon=2.0, Z=_Z(z, z_err, freq))
        res = _check_site(ed)
        assert res["Station"] == "ST1"
        assert res["Lat"] == "1.0000"
        assert res["Lon"] == "2.0000"
        assert res["Z ok"] == "YES"
        assert res["Errors"] == "YES"
        assert res["N freq"] == "5"
        assert res["NaN rows"] == "0"
        assert res["Status"] == "PASS"

    def test_missing_coords_warns(self):
        z, z_err, freq = _full_z()
        ed = _Site(station="ST2", Z=_Z(z, z_err, freq))
        res = _check_site(ed)
        assert res["Lat"] == "MISSING"
        assert res["Lon"] == "MISSING"
        assert res["Status"] == "WARN"

    def test_coords_via_header_fallback(self):
        class Header:
            lat = 10.0
            lon = 20.0

        z, z_err, freq = _full_z()
        ed = _Site(station="ST3", Header=Header(), Z=_Z(z, z_err, freq))
        res = _check_site(ed)
        assert res["Lat"] == "10.0000"
        assert res["Lon"] == "20.0000"
        assert res["Status"] == "PASS"

    def test_station_name_via_id_fallback(self):
        ed = _Site(id="ID123")
        res = _check_site(ed)
        assert res["Station"] == "ID123"

    def test_station_name_via_site_dict_fallback(self):
        ed = _Site(Site={"Name": "SiteXYZ"}, id="fallback-id")
        res = _check_site(ed)
        assert res["Station"] == "SiteXYZ"

    def test_station_name_unknown_defaults_to_qmark(self):
        ed = _Site()
        res = _check_site(ed)
        assert res["Station"] == "?"

    def test_no_z_object_fails(self):
        ed = _Site(station="ST4", lat=0.0, lon=0.0)
        res = _check_site(ed)
        assert res["Z ok"] == "NO"
        assert res["Status"] == "FAIL"

    def test_z_zero_size_fails(self):
        ed = _Site(
            station="ST5",
            lat=0.0,
            lon=0.0,
            Z=_Z(z=np.array([]), z_err=None, freq=None),
        )
        res = _check_site(ed)
        assert res["Z ok"] == "NO"
        assert res["Status"] == "FAIL"

    def test_z_object_via_edi_attribute_fallback(self):
        z, z_err, freq = _full_z()
        inner = _Site(edi=_Site(Z=_Z(z, z_err, freq)))
        # 'Z' not present directly on ed, only on ed.edi
        res = _check_site(inner)
        assert res["Z ok"] == "YES"

    def test_z_non_3d_is_partial(self):
        z = np.array([1 + 1j, 2 + 2j, 3 + 3j])  # ndim == 1
        ed = _Site(
            station="ST6",
            lat=0.0,
            lon=0.0,
            Z=_Z(z=z, z_err=np.ones(3), freq=np.array([1.0, 2.0, 3.0])),
        )
        res = _check_site(ed)
        assert res["Z ok"] == "PARTIAL"
        assert res["Status"] == "WARN"
        # NaN-rows check only runs for ndim == 3 arrays
        assert res["NaN rows"] == "—"

    def test_z_partial_xy_only(self):
        n = 5
        z = np.zeros((n, 2, 2), dtype=complex)
        z[:, 0, 1] = 1 + 1j  # xy finite
        z[:, 1, 0] = np.nan  # yx all NaN
        ed = _Site(
            station="ST7",
            lat=0.0,
            lon=0.0,
            Z=_Z(z=z, z_err=np.ones((n, 2, 2)), freq=np.linspace(1, 10, n)),
        )
        res = _check_site(ed)
        assert res["Z ok"] == "PARTIAL"
        assert res["Status"] == "WARN"
        assert res["NaN rows"] == "0"

    def test_nan_rows_counted(self):
        n = 5
        z = np.ones((n, 2, 2), dtype=complex)
        z[2, :, :] = np.nan  # entire row fully NaN
        ed = _Site(
            station="ST8",
            lat=0.0,
            lon=0.0,
            Z=_Z(z=z, z_err=np.ones((n, 2, 2)), freq=np.linspace(1, 10, n)),
        )
        res = _check_site(ed)
        assert res["Z ok"] == "YES"  # other rows keep xy/yx finite
        assert res["NaN rows"] == "1"
        assert res["Status"] == "WARN"

    def test_no_error_estimates_warns(self):
        n = 5
        z = np.ones((n, 2, 2), dtype=complex)
        ed = _Site(
            station="ST9",
            lat=0.0,
            lon=0.0,
            Z=_Z(z=z, z_err=None, freq=np.linspace(1, 10, n)),
        )
        res = _check_site(ed)
        assert res["Errors"] == "NO"
        assert res["Status"] == "WARN"

    def test_all_nan_ndarray_errors_flagged_no(self):
        n = 5
        z = np.ones((n, 2, 2), dtype=complex)
        z_err = np.full((n, 2, 2), np.nan)
        ed = _Site(
            station="ST10",
            lat=0.0,
            lon=0.0,
            Z=_Z(z=z, z_err=z_err, freq=np.linspace(1, 10, n)),
        )
        res = _check_site(ed)
        assert res["Errors"] == "NO"

    def test_all_nan_list_errors_not_detected(self):
        """
        Observation: the "all-NaN errors" check only triggers for objects
        exposing a ``.size`` attribute (ndarrays). A plain Python list of
        NaNs bypasses ``hasattr(z_err, "size")`` and is reported as "YES"
        even though every value is NaN. Locking in current behaviour here;
        see final report for discussion (not fixed, per instructions).
        """
        n = 5
        z = np.ones((n, 2, 2), dtype=complex)
        z_err = [float("nan")] * n
        ed = _Site(
            station="ST11",
            lat=0.0,
            lon=0.0,
            Z=_Z(z=z, z_err=z_err, freq=np.linspace(1, 10, n)),
        )
        res = _check_site(ed)
        assert res["Errors"] == "YES"

    def test_few_frequencies_warns(self):
        n = 5
        z = np.ones((n, 2, 2), dtype=complex)
        ed = _Site(
            station="ST12",
            lat=0.0,
            lon=0.0,
            Z=_Z(z=z, z_err=np.ones((n, 2, 2)), freq=np.array([1.0, 2.0])),
        )
        res = _check_site(ed)
        assert res["N freq"] == "2"
        assert res["Status"] == "WARN"

    # ── Status-escalates-from-PASS branches ─────────────────────────────────
    #
    # NOTE: every mock test above (and most below) uses lat=0.0, lon=0.0.
    # Because of the falsy-0.0 coordinate bug documented in
    # test_zero_frequencies_fails_when_status_pass, that already downgrades
    # Status to WARN before the Z/errors/frequency checks run, so the
    # `and result["Status"] == "PASS"` escalation guards inside those
    # checks are never actually exercised by the tests above. These use a
    # non-zero coordinate so Status is still "PASS" entering each check.

    def test_partial_z_warns_from_pass(self):
        n = 5
        z = np.zeros((n, 2, 2), dtype=complex)
        z[:, 0, 1] = 1 + 1j  # xy finite
        z[:, 1, 0] = np.nan  # yx all NaN -> partial
        ed = _Site(
            station="ST17",
            lat=1.0,
            lon=2.0,
            Z=_Z(z=z, z_err=np.ones((n, 2, 2)), freq=np.linspace(1, 10, n)),
        )
        res = _check_site(ed)
        assert res["Z ok"] == "PARTIAL"
        assert res["Status"] == "WARN"

    def test_nan_rows_warns_from_pass(self):
        n = 5
        z = np.ones((n, 2, 2), dtype=complex)
        z[2, :, :] = np.nan  # one entire row NaN, others fully finite (Z ok=YES)
        ed = _Site(
            station="ST18",
            lat=1.0,
            lon=2.0,
            Z=_Z(z=z, z_err=np.ones((n, 2, 2)), freq=np.linspace(1, 10, n)),
        )
        res = _check_site(ed)
        assert res["Z ok"] == "YES"
        assert res["NaN rows"] == "1"
        assert res["Status"] == "WARN"

    def test_no_error_estimates_warns_from_pass(self):
        n = 5
        z = np.ones((n, 2, 2), dtype=complex)
        ed = _Site(
            station="ST19",
            lat=1.0,
            lon=2.0,
            Z=_Z(z=z, z_err=None, freq=np.linspace(1, 10, n)),
        )
        res = _check_site(ed)
        assert res["Errors"] == "NO"
        assert res["Status"] == "WARN"

    def test_all_nonpositive_freqs_warns_from_pass(self):
        """All frequencies filtered out (fa.size == 0) still hits the
        len(fa) < 3 escalation, skipping the T-min/T-max computation."""
        n = 5
        z = np.ones((n, 2, 2), dtype=complex)
        ed = _Site(
            station="ST20",
            lat=1.0,
            lon=2.0,
            Z=_Z(
                z=z,
                z_err=np.ones((n, 2, 2)),
                freq=np.array([0.0, -1.0, -2.0]),
            ),
        )
        res = _check_site(ed)
        assert res["N freq"] == "0"
        assert res["T min (s)"] == "—"
        assert res["T max (s)"] == "—"
        assert res["Status"] == "WARN"

    def test_negative_and_zero_freqs_filtered(self):
        n = 5
        z = np.ones((n, 2, 2), dtype=complex)
        ed = _Site(
            station="ST13",
            lat=0.0,
            lon=0.0,
            Z=_Z(
                z=z,
                z_err=np.ones((n, 2, 2)),
                freq=np.array([0.0, -1.0, 5.0, 10.0]),
            ),
        )
        res = _check_site(ed)
        assert res["N freq"] == "2"  # only 5.0 and 10.0 survive fa > 0

    def test_zero_frequencies_fails_when_status_pass(self):
        # NOTE: lat/lon of exactly 0.0 is a separate real bug in
        # ``_check_site``: it reads coordinates via
        # ``getattr(ed, "lat", None) or getattr(ed, "latitude", None)``,
        # and 0.0 is falsy in Python, so a legitimate equatorial/prime-
        # meridian station is misread as missing and Status is prematurely
        # downgraded to WARN before the frequency check ever runs. Using a
        # non-zero coordinate here so this test actually exercises its
        # intended branch (Status starting at PASS); not fixed in
        # production, per instructions.
        n = 5
        z = np.ones((n, 2, 2), dtype=complex)
        ed = _Site(
            station="ST14",
            lat=10.0,
            lon=20.0,
            Z=_Z(z=z, z_err=np.ones((n, 2, 2)), freq=None),
        )
        res = _check_site(ed)
        assert res["N freq"] == "0"
        assert res["Status"] == "FAIL"

    def test_zero_frequencies_does_not_escalate_existing_warn(self):
        """
        Observation: Status is only escalated to FAIL on the zero-frequency
        branch when it is still "PASS" at that point (``if status == PASS``
        guard). If an earlier, milder check (e.g. missing coordinates) has
        already downgraded Status to "WARN", a station with *zero usable
        frequencies* — arguably at least as severe as "no Z" — is never
        escalated past "WARN". This looks like a real severity-ordering
        bug: the final Status can understate the worst problem found,
        depending on which check runs first. Locking in current (probably
        unintended) behaviour; not fixed, per instructions.
        """
        n = 5
        z = np.ones((n, 2, 2), dtype=complex)
        ed = _Site(
            # lat/lon intentionally omitted -> Status becomes WARN first
            station="ST15",
            Z=_Z(z=z, z_err=np.ones((n, 2, 2)), freq=None),
        )
        res = _check_site(ed)
        assert res["Lat"] == "MISSING"
        assert res["N freq"] == "0"
        assert res["Status"] == "WARN"  # not escalated to FAIL

    def test_t_min_max_computed_from_periods(self):
        n = 3
        z = np.ones((n, 2, 2), dtype=complex)
        ed = _Site(
            station="ST16",
            lat=0.0,
            lon=0.0,
            Z=_Z(
                z=z,
                z_err=np.ones((n, 2, 2)),
                freq=np.array([1.0, 10.0, 100.0]),
            ),
        )
        res = _check_site(ed)
        assert res["T min (s)"] == "0.01"
        assert res["T max (s)"] == "1"

    def test_check_site_uncaught_exception_propagates(self):
        """
        A property raising a non-AttributeError inside the Z accessor is
        NOT swallowed by getattr(..., default) and propagates out of
        _check_site itself (caught one level up, inside the dialog's
        _run_checks loop — see TestDialogRunChecks).
        """

        class RaisingZ:
            @property
            def z(self):
                raise RuntimeError("boom-z-access")

            z_err = None
            freq = None

        ed = _Site(station="BAD1", Z=RaisingZ())
        with pytest.raises(RuntimeError, match="boom-z-access"):
            _check_site(ed)


# ── _check_site — real EDI data ─────────────────────────────────────────────────


class TestCheckSiteRealData:
    def test_real_station_reports_expected_fields(self, willy_items):
        res = _check_site(willy_items[0])
        assert res["Station"] == "18-001A"
        # This dataset has no coordinates attached to the raw EDI objects.
        assert res["Lat"] == "MISSING"
        assert res["Lon"] == "MISSING"
        assert res["Z ok"] == "YES"
        assert res["Errors"] == "YES"
        assert int(res["N freq"]) > 0
        assert res["Status"] == "WARN"  # due to missing coordinates only

    def test_all_real_stations_checkable(self, willy_items):
        for ed in willy_items:
            res = _check_site(ed)
            assert res["Status"] in ("PASS", "WARN", "FAIL")
            assert res["Station"]


# ── Dialog helpers ───────────────────────────────────────────────────────────


def _select_rows(table, rows):
    table.clearSelection()
    sm = table.selectionModel()
    for r in rows:
        sm.select(
            table.model().index(r, 0),
            QItemSelectionModel.SelectionFlag.Select
            | QItemSelectionModel.SelectionFlag.Rows,
        )


# ── Dialog construction ──────────────────────────────────────────────────────


class TestDialogConstruction:
    def test_construct_with_real_sites(self, qapp, willy_sites):
        dlg = EDIValidatorDialog(willy_sites)
        assert len(dlg._items) == 5
        assert len(dlg._results) == 5
        assert dlg._table.rowCount() == 5
        assert dlg._apply_btn.isEnabled()
        assert "5 stations" in dlg._summary_lbl.text()

    def test_construct_with_none_sites(self, qapp):
        dlg = EDIValidatorDialog(None)
        assert dlg._summary_lbl.text() == "No survey loaded."
        assert dlg._items == []
        assert not dlg._apply_btn.isEnabled()

    def test_construct_with_non_iterable_sites(self, qapp):
        dlg = EDIValidatorDialog(object())
        assert dlg._items == []
        assert dlg._table.rowCount() == 0
        assert not dlg._apply_btn.isEnabled()

    def test_construct_with_plain_list_of_mocks(self, qapp):
        z, z_err, freq = _full_z()
        sites = [
            _Site(station="A", lat=1.0, lon=1.0, Z=_Z(z, z_err, freq)),
            _Site(station="B", lat=2.0, lon=2.0, Z=_Z(z, z_err, freq)),
        ]
        dlg = EDIValidatorDialog(sites)
        assert len(dlg._items) == 2
        assert dlg._results[0]["Status"] == "PASS"

    def test_per_item_exception_caught_and_marked_fail(self, qapp):
        class RaisingZ:
            @property
            def z(self):
                raise RuntimeError("boom-integration")

            z_err = None
            freq = None

        bad = _Site(station="BAD1", Z=RaisingZ())
        dlg = EDIValidatorDialog([bad])
        assert dlg._results[0]["Status"] == "FAIL"
        assert "boom-integration" in dlg._results[0]["NaN rows"]

    def test_open_recompute_requested_signal_exists(self, qapp, willy_sites):
        dlg = EDIValidatorDialog(willy_sites)
        received = []
        dlg.open_recompute_requested.connect(lambda: received.append(True))
        dlg.open_recompute_requested.emit()
        assert received == [True]

    def test_modified_sites_initially_none(self, qapp, willy_sites):
        dlg = EDIValidatorDialog(willy_sites)
        assert dlg.modified_sites is None


def _boom(*_a, **_k):
    raise RuntimeError("boom")


class TestRunChecksIterFallback:
    """_run_checks falls back to plain list()/[obj] when the
    ``pycsamt.emtools._core._iter_items`` helper itself raises."""

    def test_iter_items_raises_falls_back_to_list(self, qapp, monkeypatch):
        import pycsamt.emtools._core as core_mod

        monkeypatch.setattr(core_mod, "_iter_items", _boom)
        z, z_err, freq = _full_z()
        sites = [_Site(station="A", lat=1.0, lon=1.0, Z=_Z(z, z_err, freq))]
        dlg = EDIValidatorDialog(sites)
        assert len(dlg._items) == 1
        assert dlg._results[0]["Station"] == "A"

    def test_iter_items_and_list_fallback_both_raise_shows_message(
        self, qapp, monkeypatch
    ):
        import pycsamt.emtools._core as core_mod

        monkeypatch.setattr(core_mod, "_iter_items", _boom)

        class _BadIter:
            def __iter__(self):
                raise RuntimeError("also boom")

        dlg = EDIValidatorDialog(_BadIter())
        assert dlg._summary_lbl.text() == "Cannot iterate over loaded data."
        assert dlg._items == []

    def test_unwrap_raises_for_pycsamt_item_keeps_wrapped_item(
        self, qapp, monkeypatch, willy_sites
    ):
        """The per-item `_unwrap` call is only attempted when the item's
        type module contains "pycsamt"; force it to raise and confirm the
        loop degrades to using the still-wrapped item instead of aborting."""
        import pycsamt.emtools._core as core_mod

        monkeypatch.setattr(core_mod, "_unwrap", _boom)
        dlg = EDIValidatorDialog(willy_sites)
        assert len(dlg._items) == 5
        assert len(dlg._results) == 5


# ── _run_checks / _fill_table coloring ──────────────────────────────────────


class TestFillTableColors:
    def test_warn_rows_are_yellow(self, qapp, willy_sites):
        dlg = EDIValidatorDialog(willy_sites)
        # All 5 WILLY stations are WARN (missing coordinates only).
        for r in range(dlg._table.rowCount()):
            item = dlg._table.item(r, _COL_STATION)
            assert item.background().color().name() == _YELLOW.name()

    def test_excluded_rows_are_gray_and_struck_out(self, qapp, willy_sites):
        dlg = EDIValidatorDialog(willy_sites)
        dlg._excluded = {0}
        dlg._fill_table()
        item = dlg._table.item(0, _COL_STATION)
        assert item.background().color().name() == _EXCLUDED.name()
        assert item.font().strikeOut()

    def test_pass_row_is_green(self, qapp):
        z, z_err, freq = _full_z()
        ed = _Site(station="A", lat=1.0, lon=1.0, Z=_Z(z, z_err, freq))
        dlg = EDIValidatorDialog([ed])
        item = dlg._table.item(0, _COL_STATION)
        assert item.background().color().name() == _GREEN.name()

    def test_fail_row_is_red(self, qapp):
        ed = _Site(station="A")  # no Z at all -> FAIL
        dlg = EDIValidatorDialog([ed])
        item = dlg._table.item(0, _COL_STATION)
        assert item.background().color().name() == _RED.name()

    def test_only_station_column_editable(self, qapp, willy_sites):
        dlg = EDIValidatorDialog(willy_sites)
        from PySide6.QtCore import Qt as _Qt

        station_item = dlg._table.item(0, _COL_STATION)
        other_item = dlg._table.item(0, 1)
        assert bool(station_item.flags() & _Qt.ItemFlag.ItemIsEditable)
        assert not bool(other_item.flags() & _Qt.ItemFlag.ItemIsEditable)

    def test_summary_counts(self, qapp, willy_sites):
        dlg = EDIValidatorDialog(willy_sites)
        assert "5 stations" in dlg._summary_lbl.text()
        assert "5 warn" in dlg._summary_lbl.text()


# ── _on_item_changed ─────────────────────────────────────────────────────────


class TestOnItemChanged:
    def test_editing_station_cell_updates_results(self, qapp, willy_sites):
        dlg = EDIValidatorDialog(willy_sites)
        item = dlg._table.item(0, _COL_STATION)
        item.setText("  Renamed-Station  ")
        # itemChanged is connected -> _on_item_changed should have fired
        assert dlg._results[0]["Station"] == "Renamed-Station"

    def test_editing_non_station_cell_ignored(self, qapp, willy_sites):
        dlg = EDIValidatorDialog(willy_sites)
        before = dict(dlg._results[0])
        item = dlg._table.item(0, 1)  # Lat column, non-editable but callable directly
        dlg._on_item_changed(item)
        assert dlg._results[0] == before

    def test_populating_guard_suppresses_updates(self, qapp, willy_sites):
        dlg = EDIValidatorDialog(willy_sites)
        dlg._populating = True
        item = dlg._table.item(0, _COL_STATION)
        before = dlg._results[0]["Station"]
        item.setText("Should-Not-Apply")
        assert dlg._results[0]["Station"] == before
        dlg._populating = False

    def test_out_of_range_row_is_noop(self, qapp, willy_sites):
        dlg = EDIValidatorDialog(willy_sites)
        item = dlg._table.item(0, _COL_STATION)
        # Simulate a stale item pointing past the current results list.
        orig_row = item.row
        item.row = lambda: 999
        try:
            dlg._on_item_changed(item)  # must not raise
        finally:
            item.row = orig_row


# ── _selected_rows ───────────────────────────────────────────────────────────


class TestSelectedRows:
    def test_no_selection(self, qapp, willy_sites):
        dlg = EDIValidatorDialog(willy_sites)
        assert dlg._selected_rows() == []

    def test_single_selection(self, qapp, willy_sites):
        dlg = EDIValidatorDialog(willy_sites)
        _select_rows(dlg._table, [2])
        assert dlg._selected_rows() == [2]

    def test_multiple_selection_sorted(self, qapp, willy_sites):
        dlg = EDIValidatorDialog(willy_sites)
        _select_rows(dlg._table, [3, 1])
        assert dlg._selected_rows() == [1, 3]


# ── _on_rename ────────────────────────────────────────────────────────────────


class TestOnRename:
    def test_no_selection_does_nothing(self, qapp, willy_sites, monkeypatch):
        dlg = EDIValidatorDialog(willy_sites)
        calls = []
        monkeypatch.setattr(
            QInputDialog,
            "getText",
            staticmethod(lambda *a, **k: calls.append(1) or ("X", True)),
        )
        dlg._on_rename()
        assert calls == []

    def test_valid_rename_applies(self, qapp, willy_sites, monkeypatch):
        dlg = EDIValidatorDialog(willy_sites)
        _select_rows(dlg._table, [0])
        monkeypatch.setattr(
            QInputDialog,
            "getText",
            staticmethod(lambda *a, **k: ("NewName", True)),
        )
        dlg._on_rename()
        assert dlg._results[0]["Station"] == "NewName"

    def test_cancelled_rename_ignored(self, qapp, willy_sites, monkeypatch):
        dlg = EDIValidatorDialog(willy_sites)
        _select_rows(dlg._table, [0])
        before = dlg._results[0]["Station"]
        monkeypatch.setattr(
            QInputDialog,
            "getText",
            staticmethod(lambda *a, **k: ("NewName", False)),
        )
        dlg._on_rename()
        assert dlg._results[0]["Station"] == before

    def test_empty_rename_ignored(self, qapp, willy_sites, monkeypatch):
        dlg = EDIValidatorDialog(willy_sites)
        _select_rows(dlg._table, [0])
        before = dlg._results[0]["Station"]
        monkeypatch.setattr(
            QInputDialog,
            "getText",
            staticmethod(lambda *a, **k: ("   ", True)),
        )
        dlg._on_rename()
        assert dlg._results[0]["Station"] == before

    def test_unchanged_name_ignored(self, qapp, willy_sites, monkeypatch):
        dlg = EDIValidatorDialog(willy_sites)
        _select_rows(dlg._table, [0])
        old_name = dlg._results[0]["Station"]
        monkeypatch.setattr(
            QInputDialog,
            "getText",
            staticmethod(lambda *a, **k: (old_name, True)),
        )
        dlg._on_rename()
        assert dlg._results[0]["Station"] == old_name


# ── _on_delete ────────────────────────────────────────────────────────────────


class TestOnDelete:
    def test_no_selection_does_nothing(self, qapp, willy_sites):
        dlg = EDIValidatorDialog(willy_sites)
        dlg._on_delete()
        assert len(dlg._items) == 5

    def test_confirmed_single_delete(self, qapp, willy_sites):
        # no_modal_dialogs autouse fixture makes QMessageBox.question -> Yes
        dlg = EDIValidatorDialog(willy_sites)
        _select_rows(dlg._table, [1])
        dlg._on_delete()
        assert len(dlg._items) == 4
        assert len(dlg._results) == 4
        assert dlg._table.rowCount() == 4

    def test_confirmed_multi_delete(self, qapp, willy_sites):
        dlg = EDIValidatorDialog(willy_sites)
        _select_rows(dlg._table, [0, 2, 4])
        dlg._on_delete()
        assert len(dlg._items) == 2

    def test_cancelled_delete_keeps_all(self, qapp, willy_sites, monkeypatch):
        dlg = EDIValidatorDialog(willy_sites)
        _select_rows(dlg._table, [0])
        monkeypatch.setattr(
            QMessageBox,
            "question",
            staticmethod(
                lambda *a, **k: QMessageBox.StandardButton.Cancel
            ),
        )
        dlg._on_delete()
        assert len(dlg._items) == 5

    def test_delete_message_truncates_long_selection(
        self, qapp, willy_sites, monkeypatch
    ):
        dlg = EDIValidatorDialog(willy_sites)
        _select_rows(dlg._table, [0, 1, 2, 3, 4])
        captured = {}

        def _fake_question(*a, **k):
            captured["text"] = a[2] if len(a) > 2 else k.get("text", "")
            return QMessageBox.StandardButton.Yes

        monkeypatch.setattr(
            QMessageBox, "question", staticmethod(_fake_question)
        )
        dlg._on_delete()
        assert "more" in captured["text"]
        assert len(dlg._items) == 0

    def test_excluded_indices_remapped_after_delete(self, qapp, willy_sites):
        dlg = EDIValidatorDialog(willy_sites)
        dlg._excluded = {1, 3}  # station indices 1 and 3 excluded
        _select_rows(dlg._table, [0])  # delete row 0
        dlg._on_delete()
        # after removing row 0, old row1 -> new row0, old row3 -> new row2
        assert dlg._excluded == {0, 2}

    def test_deleted_index_dropped_from_excluded(self, qapp, willy_sites):
        dlg = EDIValidatorDialog(willy_sites)
        dlg._excluded = {1}
        _select_rows(dlg._table, [1])
        dlg._on_delete()
        assert dlg._excluded == set()


# ── _on_toggle_exclude ────────────────────────────────────────────────────────


class TestOnToggleExclude:
    def test_no_selection_does_nothing(self, qapp, willy_sites):
        dlg = EDIValidatorDialog(willy_sites)
        dlg._on_toggle_exclude()
        assert dlg._excluded == set()

    def test_toggle_excludes_then_includes(self, qapp, willy_sites):
        dlg = EDIValidatorDialog(willy_sites)
        _select_rows(dlg._table, [1])
        dlg._on_toggle_exclude()
        assert dlg._excluded == {1}
        _select_rows(dlg._table, [1])
        dlg._on_toggle_exclude()
        assert dlg._excluded == set()

    def test_mixed_selection_excludes_all(self, qapp, willy_sites):
        dlg = EDIValidatorDialog(willy_sites)
        dlg._excluded = {0}
        _select_rows(dlg._table, [0, 1])
        dlg._on_toggle_exclude()
        assert dlg._excluded == {0, 1}


# ── _on_move_up / _on_move_down ──────────────────────────────────────────────


class TestOnMove:
    def test_move_up_first_row_is_noop(self, qapp, willy_sites):
        dlg = EDIValidatorDialog(willy_sites)
        names_before = [r["Station"] for r in dlg._results]
        _select_rows(dlg._table, [0])
        dlg._on_move_up()
        assert [r["Station"] for r in dlg._results] == names_before

    def test_move_down_last_row_is_noop(self, qapp, willy_sites):
        dlg = EDIValidatorDialog(willy_sites)
        n = len(dlg._results)
        names_before = [r["Station"] for r in dlg._results]
        _select_rows(dlg._table, [n - 1])
        dlg._on_move_down()
        assert [r["Station"] for r in dlg._results] == names_before

    def test_move_up_swaps_with_previous(self, qapp, willy_sites):
        dlg = EDIValidatorDialog(willy_sites)
        names_before = [r["Station"] for r in dlg._results]
        _select_rows(dlg._table, [1])
        dlg._on_move_up()
        names_after = [r["Station"] for r in dlg._results]
        assert names_after[0] == names_before[1]
        assert names_after[1] == names_before[0]
        assert names_after[2:] == names_before[2:]

    def test_move_down_swaps_with_next(self, qapp, willy_sites):
        dlg = EDIValidatorDialog(willy_sites)
        names_before = [r["Station"] for r in dlg._results]
        _select_rows(dlg._table, [1])
        dlg._on_move_down()
        names_after = [r["Station"] for r in dlg._results]
        assert names_after[1] == names_before[2]
        assert names_after[2] == names_before[1]

    def test_move_up_carries_excluded_flag(self, qapp, willy_sites):
        dlg = EDIValidatorDialog(willy_sites)
        dlg._excluded = {2}
        _select_rows(dlg._table, [2])
        dlg._on_move_up()
        assert dlg._excluded == {1}

    def test_move_down_carries_excluded_flag(self, qapp, willy_sites):
        dlg = EDIValidatorDialog(willy_sites)
        dlg._excluded = {2}
        _select_rows(dlg._table, [2])
        dlg._on_move_down()
        assert dlg._excluded == {3}

    def test_move_up_multi_row_selection(self, qapp, willy_sites):
        dlg = EDIValidatorDialog(willy_sites)
        names_before = [r["Station"] for r in dlg._results]
        _select_rows(dlg._table, [1, 2])
        dlg._on_move_up()
        names_after = [r["Station"] for r in dlg._results]
        # both rows 1 and 2 shift up by one; row 0 pushed down
        assert names_after[0] == names_before[1]
        assert names_after[1] == names_before[2]
        assert names_after[2] == names_before[0]

    def test_move_down_multi_row_selection(self, qapp, willy_sites):
        dlg = EDIValidatorDialog(willy_sites)
        n = len(dlg._results)
        names_before = [r["Station"] for r in dlg._results]
        # rows [n-2, n-1] would include the last row, which is a no-op
        # per _on_move_down's own `rows[-1] >= n - 1` guard (see
        # test_move_down_last_row_is_noop) -- select one row short of
        # that so the cascade actually runs.
        _select_rows(dlg._table, [n - 3, n - 2])
        dlg._on_move_down()
        names_after = [r["Station"] for r in dlg._results]
        # both selected rows shift down by one; the row below them is
        # pushed up (mirrors test_move_up_multi_row_selection)
        assert names_after[-1] == names_before[-2]
        assert names_after[-2] == names_before[-3]


# ── _show_context_menu ────────────────────────────────────────────────────────


class TestShowContextMenu:
    def test_no_selection_menu_not_shown(self, qapp, willy_sites, monkeypatch):
        dlg = EDIValidatorDialog(willy_sites)
        calls = []

        class _FakeMenu(QMenu):
            def exec(self, *a, **k):
                calls.append(1)
                return None

        monkeypatch.setattr(_validator_tool_mod, "QMenu", _FakeMenu)
        dlg._show_context_menu(QPoint(0, 0))
        assert calls == []

    def test_selection_shows_menu_with_exclude_label(
        self, qapp, willy_sites, monkeypatch
    ):
        # NOTE: QMenu.exec is a compiled/Shiboken-bound method — reassigning
        # it via monkeypatch.setattr(QMenu, "exec", ...) is silently ignored
        # (the real modal exec() still runs, which hangs forever offscreen
        # since nothing ever dismisses it). Subclassing and swapping the
        # module-level QMenu name is the only override that actually takes.
        dlg = EDIValidatorDialog(willy_sites)
        _select_rows(dlg._table, [0])
        captured = {}

        class _FakeMenu(QMenu):
            def exec(self, *a, **k):
                captured["actions"] = [act.text() for act in self.actions()]
                return None

        monkeypatch.setattr(_validator_tool_mod, "QMenu", _FakeMenu)
        dlg._show_context_menu(QPoint(1, 1))
        assert any("Exclude" in t for t in captured["actions"])
        assert any("Rename" in t for t in captured["actions"])
        assert any("Delete" in t for t in captured["actions"])
        assert any("Move Up" in t for t in captured["actions"])
        assert any("Move Down" in t for t in captured["actions"])

    def test_all_excluded_selection_shows_include_label(
        self, qapp, willy_sites, monkeypatch
    ):
        dlg = EDIValidatorDialog(willy_sites)
        dlg._excluded = {0}
        _select_rows(dlg._table, [0])
        captured = {}

        class _FakeMenu(QMenu):
            def exec(self, *a, **k):
                captured["actions"] = [act.text() for act in self.actions()]
                return None

        monkeypatch.setattr(_validator_tool_mod, "QMenu", _FakeMenu)
        dlg._show_context_menu(QPoint(1, 1))
        assert any("Include" in t for t in captured["actions"])


# ── _on_apply ─────────────────────────────────────────────────────────────────


class TestOnApply:
    def test_apply_without_exclusions(self, qapp, willy_sites):
        from PySide6.QtWidgets import QDialog

        dlg = EDIValidatorDialog(willy_sites)
        dlg._on_apply()
        assert dlg.modified_sites is not None
        assert len(dlg.modified_sites) == 5
        assert not dlg._apply_btn.isEnabled()
        assert dlg.result() == QDialog.DialogCode.Accepted

    def test_apply_with_exclusions_filters_list(self, qapp, willy_sites):
        dlg = EDIValidatorDialog(willy_sites)
        _select_rows(dlg._table, [0, 1])
        dlg._on_toggle_exclude()
        dlg._on_apply()
        assert len(dlg.modified_sites) == 3
        assert "excluded" in dlg._summary_lbl.text()

    def test_apply_pushes_rename_to_underlying_object(
        self, qapp, willy_sites
    ):
        dlg = EDIValidatorDialog(willy_sites)
        dlg._results[0]["Station"] = "RenamedByApply"
        dlg._on_apply()
        assert getattr(dlg._items[0], "station", None) == "RenamedByApply"

    def test_apply_with_empty_name_does_not_touch_object(self, qapp):
        ed = _Site(station="Orig")
        dlg = EDIValidatorDialog([ed])
        dlg._results[0]["Station"] = ""
        dlg._on_apply()
        assert ed.station == "Orig"

    def test_apply_setattr_failure_on_one_attr_is_swallowed(self, qapp):
        """A read-only ``station`` property raises on assignment; the
        per-attr try/except in _on_apply must swallow it and move on
        instead of aborting the whole Apply."""

        class _ReadOnlyStation:
            @property
            def station(self):
                return "OLD"

        ed = _ReadOnlyStation()
        dlg = EDIValidatorDialog([ed])
        dlg._results[0]["Station"] = "NewName"
        dlg._on_apply()  # must not raise
        assert ed.station == "OLD"  # setattr silently failed
        assert dlg.modified_sites == [ed]


# ── _on_export ────────────────────────────────────────────────────────────────


class TestOnExport:
    def test_no_results_noop(self, qapp, monkeypatch):
        dlg = EDIValidatorDialog(None)
        calls = []
        monkeypatch.setattr(
            QFileDialog,
            "getSaveFileName",
            staticmethod(lambda *a, **k: calls.append(1) or ("", "")),
        )
        dlg._on_export()
        assert calls == []

    def test_cancel_save_dialog_no_file(self, qapp, willy_sites, monkeypatch, tmp_path):
        dlg = EDIValidatorDialog(willy_sites)
        monkeypatch.setattr(
            QFileDialog,
            "getSaveFileName",
            staticmethod(lambda *a, **k: ("", "CSV (*.csv)")),
        )
        dlg._on_export()
        assert not any(tmp_path.iterdir())

    def test_successful_export_writes_csv(
        self, qapp, willy_sites, monkeypatch, tmp_path
    ):
        out = tmp_path / "report.csv"
        monkeypatch.setattr(
            QFileDialog,
            "getSaveFileName",
            staticmethod(lambda *a, **k: (str(out), "CSV (*.csv)")),
        )
        dlg = EDIValidatorDialog(willy_sites)
        dlg._on_export()
        assert out.exists()
        with open(out, newline="", encoding="utf-8") as fh:
            rows = list(csv.DictReader(fh))
        assert len(rows) == 5
        assert set(_COLS + ["Excluded"]).issubset(rows[0].keys())
        assert all(r["Excluded"] == "no" for r in rows)
        assert "saved" in dlg._summary_lbl.text()

    def test_export_marks_excluded_rows(
        self, qapp, willy_sites, monkeypatch, tmp_path
    ):
        out = tmp_path / "report2.csv"
        monkeypatch.setattr(
            QFileDialog,
            "getSaveFileName",
            staticmethod(lambda *a, **k: (str(out), "CSV (*.csv)")),
        )
        dlg = EDIValidatorDialog(willy_sites)
        dlg._excluded = {0}
        dlg._on_export()
        with open(out, newline="", encoding="utf-8") as fh:
            rows = list(csv.DictReader(fh))
        assert rows[0]["Excluded"] == "yes"
        assert rows[1]["Excluded"] == "no"

    def test_export_error_reports_message(
        self, qapp, willy_sites, monkeypatch, tmp_path
    ):
        # Point the export at a directory path -> open() raises IsADirectoryError/
        # PermissionError, exercised through the except branch.
        target_dir = tmp_path / "a_directory"
        target_dir.mkdir()
        monkeypatch.setattr(
            QFileDialog,
            "getSaveFileName",
            staticmethod(lambda *a, **k: (str(target_dir), "CSV (*.csv)")),
        )
        dlg = EDIValidatorDialog(willy_sites)
        dlg._on_export()
        assert "Export error" in dlg._summary_lbl.text()
