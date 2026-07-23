# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for SettingsController — no Qt required.

``SettingsController`` is a thin, Qt-free wrapper around several
process-wide ``PYCSAMT_*`` singletons:

* :data:`pycsamt.api.control.PYCSAMT_CONTROL`             (view controls)
* :data:`pycsamt.api.station.PYCSAMT_STATION_RENDERING`   (station markers)
* :data:`pycsamt.api.section.PYCSAMT_SECTION`             (section/axis style)
* :data:`pycsamt.topo.PYCSAMT_TOPO`                       (topography)
* :data:`pycsamt.api.style.PYCSAMT_STYLE`                 (display, reset only)
* :data:`pycsamt.api.interp.PYCSAMT_INTERP`               (interpretation, reset only)

All of these are module-level singletons that persist for the whole
pytest session, so every test that mutates them must leave them exactly
as it found them.  We rely on ``SettingsController.snapshot()`` /
``.restore()`` themselves as the isolation mechanism (autouse fixture
below) — they cover every field any ``apply_*`` method can write, so
round-tripping through them is sufficient to undo any test's mutations.
``PYCSAMT_STYLE`` / ``PYCSAMT_INTERP`` are only ever touched via
``.reset()`` (idempotent — always returns package defaults), so they
need no separate snapshot.
"""

from __future__ import annotations

import json
import sys

import pytest

from pycsamt.app.desktop.controllers.settings_controller import (
    SettingsController,
)

# ── Isolation ────────────────────────────────────────────────────────────


@pytest.fixture(autouse=True)
def _isolate_singletons():
    """Snapshot every PYCSAMT_* field this controller can write, restore after."""
    ctrl = SettingsController()
    snap = ctrl.snapshot()
    yield
    ctrl.restore(snap)


@pytest.fixture()
def ctrl():
    return SettingsController()


def _break_import(monkeypatch, module_name):
    """Force ``import module_name`` (or ``from module_name import x``) to raise."""
    monkeypatch.setitem(sys.modules, module_name, None)


# ── snapshot() ───────────────────────────────────────────────────────────


def test_snapshot_has_expected_top_level_keys(ctrl):
    snap = ctrl.snapshot()
    assert set(snap.keys()) == {
        "view_controls",
        "station",
        "section",
        "topography",
    }


def test_snapshot_reflects_current_singleton_state(ctrl):
    ctrl.apply_view_controls(rho_view="linear", phase_unit="radian")
    ctrl.apply_station(side="bottom", max_labels=5)
    ctrl.apply_section(y_direction="up", station_side="bottom")
    ctrl.apply_topography(enabled=True, exaggeration=3.5)

    snap = ctrl.snapshot()
    assert snap["view_controls"]["rho_view"] == "linear"
    assert snap["view_controls"]["phase_unit"] == "radian"
    assert snap["station"]["side"] == "bottom"
    assert snap["station"]["max_labels"] == 5
    assert snap["section"]["y_direction"] == "up"
    assert snap["section"]["station_side"] == "bottom"
    assert snap["topography"]["enabled"] is True
    assert snap["topography"]["exaggeration"] == 3.5


def test_snapshot_phase_range_is_a_list(ctrl):
    # snapshot() explicitly does list(C.phase.range) for JSON-friendliness.
    snap = ctrl.snapshot()
    assert isinstance(snap["view_controls"]["phase_range"], list)


def test_snapshot_swallows_import_error_per_section(monkeypatch, ctrl):
    _break_import(monkeypatch, "pycsamt.api.control")
    snap = ctrl.snapshot()
    assert "view_controls" not in snap
    # other sections still populated
    assert "station" in snap
    assert "section" in snap
    assert "topography" in snap


def test_snapshot_swallows_import_error_station(monkeypatch, ctrl):
    _break_import(monkeypatch, "pycsamt.api.station")
    snap = ctrl.snapshot()
    assert "station" not in snap
    assert "view_controls" in snap
    assert "section" in snap
    assert "topography" in snap


def test_snapshot_swallows_import_error_section(monkeypatch, ctrl):
    _break_import(monkeypatch, "pycsamt.api.section")
    snap = ctrl.snapshot()
    assert "section" not in snap
    assert "view_controls" in snap
    assert "station" in snap
    assert "topography" in snap


def test_snapshot_swallows_import_error_topography(monkeypatch, ctrl):
    _break_import(monkeypatch, "pycsamt.topo")
    snap = ctrl.snapshot()
    assert "topography" not in snap
    assert "view_controls" in snap
    assert "station" in snap
    assert "section" in snap


# ── restore() ────────────────────────────────────────────────────────────


def test_restore_round_trip(ctrl):
    baseline = ctrl.snapshot()

    ctrl.apply_view_controls(rho_view="linear")
    ctrl.apply_station(side="bottom")
    ctrl.apply_section(y_direction="up")
    ctrl.apply_topography(enabled=True)

    ctrl.restore(baseline)
    assert ctrl.snapshot() == baseline


def test_restore_empty_dict_is_noop(ctrl):
    baseline = ctrl.snapshot()
    ctrl.restore({})
    assert ctrl.snapshot() == baseline


def test_restore_partial_dict_only_touches_that_section(ctrl):
    baseline = ctrl.snapshot()
    ctrl.apply_view_controls(rho_view="linear")
    ctrl.apply_topography(enabled=True)

    # restore only topography; view_controls edit must survive
    ctrl.restore({"topography": baseline["topography"]})

    snap = ctrl.snapshot()
    assert snap["view_controls"]["rho_view"] == "linear"
    assert snap["topography"]["enabled"] == baseline["topography"]["enabled"]


# ── apply_view_controls() ───────────────────────────────────────────────


def test_apply_view_controls_all_fields(ctrl):
    from pycsamt.api.control import PYCSAMT_CONTROL as C

    ctrl.apply_view_controls(
        rho_view="linear",
        phase_range=[-90.0, 90.0],
        phase_unit="radian",
        phase_wrap=False,
        x_view="frequency",
    )
    assert C.rho.view == "linear"
    assert C.phase.range == (-90.0, 90.0)
    assert isinstance(C.phase.range, tuple)
    assert C.phase.unit == "radian"
    assert C.phase.wrap is False
    assert C.x.view == "frequency"


def test_apply_view_controls_partial_kwargs_leaves_rest_untouched(ctrl):
    from pycsamt.api.control import PYCSAMT_CONTROL as C

    before_x_view = C.x.view
    ctrl.apply_view_controls(rho_view="linear")
    assert C.rho.view == "linear"
    assert C.x.view == before_x_view


def test_apply_view_controls_phase_wrap_bool_coercion(ctrl):
    from pycsamt.api.control import PYCSAMT_CONTROL as C

    ctrl.apply_view_controls(phase_wrap=0)
    assert C.phase.wrap is False
    ctrl.apply_view_controls(phase_wrap=1)
    assert C.phase.wrap is True


def test_apply_view_controls_import_error_is_swallowed(monkeypatch, ctrl):
    from pycsamt.api.control import PYCSAMT_CONTROL as C

    before = C.rho.view
    _break_import(monkeypatch, "pycsamt.api.control")
    ctrl.apply_view_controls(rho_view="linear")  # must not raise
    assert C.rho.view == before


def test_apply_view_controls_bad_value_aborts_remaining_kwargs(ctrl):
    """A ValueError on one field silently aborts the *rest* of the call.

    apply_view_controls wraps its whole body in one try/except, so if an
    earlier kwarg's coercion raises, later kwargs in the same call are
    never applied (documented here, not fixed — see final report).
    """
    from pycsamt.api.control import PYCSAMT_CONTROL as C

    before_x_view = C.x.view
    # phase_range is processed before x_view; tuple() of a non-iterable
    # int raises TypeError, which the broad except swallows.
    ctrl.apply_view_controls(phase_range=123, x_view="frequency")
    assert C.x.view == before_x_view  # never reached


# ── apply_station() ──────────────────────────────────────────────────────


def test_apply_station_all_fields(ctrl):
    from pycsamt.api.station import PYCSAMT_STATION_RENDERING as SR

    ctrl.apply_station(
        side="bottom",
        show_markers=False,
        marker_symbol="o",
        marker_size=50.0,
        marker_offset=2.0,
        max_labels=8,
    )
    ps = SR.pseudosection
    assert ps.side == "bottom"
    assert ps.show_markers is False
    assert ps.marker.marker == "o"
    assert ps.marker.size == 50.0
    assert ps.marker.offset == 2.0
    assert ps.max_labels == 8


def test_apply_station_partial_kwargs_leaves_rest_untouched(ctrl):
    from pycsamt.api.station import PYCSAMT_STATION_RENDERING as SR

    before_symbol = SR.pseudosection.marker.marker
    ctrl.apply_station(side="bottom")
    assert SR.pseudosection.side == "bottom"
    assert SR.pseudosection.marker.marker == before_symbol


def test_apply_station_numeric_string_coercion(ctrl):
    from pycsamt.api.station import PYCSAMT_STATION_RENDERING as SR

    ctrl.apply_station(marker_size="40.5", marker_offset="1.5", max_labels="9")
    ps = SR.pseudosection
    assert ps.marker.size == 40.5
    assert ps.marker.offset == 1.5
    assert ps.max_labels == 9


def test_apply_station_import_error_is_swallowed(monkeypatch, ctrl):
    from pycsamt.api.station import PYCSAMT_STATION_RENDERING as SR

    before = SR.pseudosection.side
    _break_import(monkeypatch, "pycsamt.api.station")
    ctrl.apply_station(side="bottom")  # must not raise
    assert SR.pseudosection.side == before


def test_apply_station_bad_value_aborts_remaining_kwargs(ctrl):
    from pycsamt.api.station import PYCSAMT_STATION_RENDERING as SR

    before_max_labels = SR.pseudosection.max_labels
    ctrl.apply_station(marker_size="not-a-number", max_labels=99)
    # marker_size coercion raised first -> max_labels never applied
    assert SR.pseudosection.max_labels == before_max_labels


# ── apply_section() ──────────────────────────────────────────────────────


_APPLIED_PRESETS = ("pseudosection", "dashboard", "compact", "publication", "dynamic")


def test_apply_section_writes_all_applied_presets(ctrl):
    from pycsamt.api.section import PYCSAMT_SECTION as SEC

    ctrl.apply_section(y_direction="up", station_side="bottom")
    for name in _APPLIED_PRESETS:
        ax = SEC.style_for(name).axis
        assert ax.y_direction == "up"
        assert ax.station_side == "bottom"


def test_apply_section_does_not_touch_inversion_preset(ctrl):
    """apply_section only writes the 'pseudosection-family' presets per its
    docstring; 'inversion' is deliberately excluded (not in the hard-coded
    tuple), even though it is a valid PYCSAMT_SECTION preset elsewhere."""
    from pycsamt.api.section import PYCSAMT_SECTION as SEC

    before = SEC.inversion.axis.y_direction
    ctrl.apply_section(y_direction="up")
    assert SEC.inversion.axis.y_direction == before


def test_apply_section_partial_kwargs(ctrl):
    from pycsamt.api.section import PYCSAMT_SECTION as SEC

    before_side = SEC.pseudosection.axis.station_side
    ctrl.apply_section(y_direction="up")
    assert SEC.pseudosection.axis.y_direction == "up"
    assert SEC.pseudosection.axis.station_side == before_side


def test_apply_section_station_side_only(ctrl):
    # Exercise the other half of the "y_direction"/"station_side" branch
    # pair inside the per-preset loop.
    from pycsamt.api.section import PYCSAMT_SECTION as SEC

    before_y = SEC.pseudosection.axis.y_direction
    ctrl.apply_section(station_side="bottom")
    assert SEC.pseudosection.axis.station_side == "bottom"
    assert SEC.pseudosection.axis.y_direction == before_y


def test_apply_section_import_error_is_swallowed(monkeypatch, ctrl):
    from pycsamt.api.section import PYCSAMT_SECTION as SEC

    before = SEC.pseudosection.axis.y_direction
    _break_import(monkeypatch, "pycsamt.api.section")
    ctrl.apply_section(y_direction="up")  # must not raise
    assert SEC.pseudosection.axis.y_direction == before


def test_apply_section_isolates_per_preset_failures(monkeypatch, ctrl):
    """One failing preset (style_for raising) must not block the others."""
    from pycsamt.api.section import PYCSAMT_SECTION as SEC

    orig_style_for = SEC.style_for

    def flaky_style_for(preset):
        if preset == "compact":
            raise RuntimeError("boom")
        return orig_style_for(preset)

    monkeypatch.setattr(SEC, "style_for", flaky_style_for)
    before_compact = SEC.compact.axis.y_direction

    ctrl.apply_section(y_direction="up")

    assert SEC.pseudosection.axis.y_direction == "up"
    assert SEC.dashboard.axis.y_direction == "up"
    assert SEC.publication.axis.y_direction == "up"
    assert SEC.dynamic.axis.y_direction == "up"
    assert SEC.compact.axis.y_direction == before_compact  # untouched


# ── apply_topography() ───────────────────────────────────────────────────


def test_apply_topography_all_fields(ctrl):
    from pycsamt.topo import PYCSAMT_TOPO as T

    ctrl.apply_topography(enabled=True, exaggeration=2.5, marker_pad=0.05)
    assert T.enabled is True
    assert T.exaggeration == 2.5
    assert T.marker_pad_fraction == 0.05


def test_apply_topography_partial_kwargs(ctrl):
    from pycsamt.topo import PYCSAMT_TOPO as T

    before_exag = T.exaggeration
    ctrl.apply_topography(enabled=True)
    assert T.enabled is True
    assert T.exaggeration == before_exag


def test_apply_topography_bool_coercion(ctrl):
    from pycsamt.topo import PYCSAMT_TOPO as T

    ctrl.apply_topography(enabled=1)
    assert T.enabled is True
    ctrl.apply_topography(enabled=0)
    assert T.enabled is False


def test_apply_topography_import_error_is_swallowed(monkeypatch, ctrl):
    from pycsamt.topo import PYCSAMT_TOPO as T

    before = T.enabled
    _break_import(monkeypatch, "pycsamt.topo")
    ctrl.apply_topography(enabled=not before)  # must not raise
    assert T.enabled == before


def test_apply_topography_bad_value_aborts_remaining_kwargs(ctrl):
    from pycsamt.topo import PYCSAMT_TOPO as T

    before_pad = T.marker_pad_fraction
    ctrl.apply_topography(exaggeration="not-a-number", marker_pad=0.9)
    assert T.marker_pad_fraction == before_pad  # never reached


# ── reset_tab() ──────────────────────────────────────────────────────────


def test_reset_tab_view_controls(ctrl):
    from pycsamt.api.control import PYCSAMT_CONTROL as C

    ctrl.apply_view_controls(rho_view="linear", phase_unit="radian")
    ctrl.reset_tab("view_controls")
    assert C.rho.view == "log10"
    assert C.phase.unit == "degree"
    assert C.phase.range == (-180.0, 180.0)
    assert C.phase.wrap is True
    assert C.x.view == "log10_period"


def test_reset_tab_pseudosections(ctrl):
    from pycsamt.api.section import PYCSAMT_SECTION as SEC
    from pycsamt.api.station import PYCSAMT_STATION_RENDERING as SR

    ctrl.apply_station(side="bottom", marker_size=99.0)
    ctrl.apply_section(y_direction="up", station_side="bottom")
    ctrl.reset_tab("pseudosections")

    assert SR.pseudosection.side == "top"
    assert SR.pseudosection.marker.size == 34.0
    assert SEC.pseudosection.axis.y_direction == "down"
    assert SEC.pseudosection.axis.station_side == "top"


def test_reset_tab_topography(ctrl):
    from pycsamt.topo import PYCSAMT_TOPO as T

    ctrl.apply_topography(enabled=True, exaggeration=5.0)
    ctrl.reset_tab("topography")
    assert T.enabled is False
    assert T.exaggeration == 1.0
    assert T.marker_pad_fraction == 0.015


def test_reset_tab_display(ctrl):
    # PYCSAMT_STYLE has no apply_* counterpart in this controller; dirty it
    # directly to confirm reset_tab("display") really reaches
    # PYCSAMT_STYLE.reset().
    from pycsamt.api.style import PYCSAMT_STYLE as S

    S.rose.bar_style = "flat"
    ctrl.reset_tab("display")
    assert S.rose.bar_style == "gradient"


def test_reset_tab_interpretation(ctrl):
    # Same rationale as above, for PYCSAMT_INTERP.
    from pycsamt.api.interp import PYCSAMT_INTERP as I

    I.default.section.wt_color = "chartreuse"
    ctrl.reset_tab("interpretation")
    assert I.default.section.wt_color == "deepskyblue"


def test_reset_tab_unknown_tab_is_noop(ctrl):
    from pycsamt.api.control import PYCSAMT_CONTROL as C

    ctrl.apply_view_controls(rho_view="linear")
    ctrl.reset_tab("not_a_real_tab")  # must not raise, must not reset anything
    assert C.rho.view == "linear"


def test_reset_tab_import_error_is_swallowed(monkeypatch, ctrl):
    from pycsamt.api.control import PYCSAMT_CONTROL as C

    ctrl.apply_view_controls(rho_view="linear")
    _break_import(monkeypatch, "pycsamt.api.control")
    ctrl.reset_tab("view_controls")  # must not raise
    assert C.rho.view == "linear"  # reset never happened


def test_reset_tab_pseudosections_import_error_is_swallowed(monkeypatch, ctrl):
    from pycsamt.api.station import PYCSAMT_STATION_RENDERING as SR

    ctrl.apply_station(side="bottom")
    _break_import(monkeypatch, "pycsamt.api.section")
    ctrl.reset_tab("pseudosections")  # must not raise
    # section import broke before either .reset() call executed
    assert SR.pseudosection.side == "bottom"


# ── reset_all() ──────────────────────────────────────────────────────────


def test_reset_all_reverts_every_singleton(ctrl):
    from pycsamt.api.control import PYCSAMT_CONTROL as C
    from pycsamt.api.section import PYCSAMT_SECTION as SEC
    from pycsamt.api.station import PYCSAMT_STATION_RENDERING as SR
    from pycsamt.topo import PYCSAMT_TOPO as T

    ctrl.apply_view_controls(rho_view="linear")
    ctrl.apply_station(side="bottom")
    ctrl.apply_section(y_direction="up")
    ctrl.apply_topography(enabled=True)

    ctrl.reset_all()

    assert C.rho.view == "log10"
    assert SR.pseudosection.side == "top"
    assert SEC.pseudosection.axis.y_direction == "down"
    assert T.enabled is False


# ── save() / load() ──────────────────────────────────────────────────────


def test_save_writes_json_matching_snapshot(monkeypatch, ctrl, tmp_path):
    target = tmp_path / "nested" / "settings.json"
    monkeypatch.setattr(SettingsController, "SETTINGS_PATH", target)

    ctrl.apply_view_controls(rho_view="linear")
    ctrl.save()

    assert target.exists()
    on_disk = json.loads(target.read_text())
    assert on_disk == ctrl.snapshot()


def test_save_with_explicit_path(ctrl, tmp_path):
    target = tmp_path / "profile.json"
    ctrl.apply_topography(enabled=True)
    ctrl.save(target)

    assert target.exists()
    on_disk = json.loads(target.read_text())
    assert on_disk["topography"]["enabled"] is True


def test_save_accepts_string_path(ctrl, tmp_path):
    target = tmp_path / "profile_str.json"
    ctrl.save(str(target))
    assert target.exists()


def test_load_round_trip(ctrl, tmp_path):
    target = tmp_path / "profile.json"
    ctrl.apply_view_controls(rho_view="linear")
    ctrl.apply_station(side="bottom")
    ctrl.save(target)
    saved_snap = ctrl.snapshot()

    # dirty the singletons differently, then load back
    ctrl.apply_view_controls(rho_view="rho")
    ctrl.apply_station(side="top")

    ok = ctrl.load(target)
    assert ok is True
    assert ctrl.snapshot() == saved_snap


def test_load_missing_file_returns_false(ctrl, tmp_path):
    missing = tmp_path / "does_not_exist.json"
    assert ctrl.load(missing) is False


def test_load_corrupt_json_returns_false(ctrl, tmp_path):
    bad = tmp_path / "corrupt.json"
    bad.write_text("{not valid json::")
    before = ctrl.snapshot()

    assert ctrl.load(bad) is False
    assert ctrl.snapshot() == before  # nothing mutated


def test_load_default_path(monkeypatch, ctrl, tmp_path):
    target = tmp_path / "default_settings.json"
    monkeypatch.setattr(SettingsController, "SETTINGS_PATH", target)

    ctrl.apply_topography(enabled=True, exaggeration=4.0)
    ctrl.save()

    ctrl.apply_topography(enabled=False, exaggeration=1.0)
    assert ctrl.load() is True

    from pycsamt.topo import PYCSAMT_TOPO as T

    assert T.enabled is True
    assert T.exaggeration == 4.0


def test_load_missing_default_path_returns_false(monkeypatch, ctrl, tmp_path):
    target = tmp_path / "never_saved.json"
    monkeypatch.setattr(SettingsController, "SETTINGS_PATH", target)
    assert ctrl.load() is False
