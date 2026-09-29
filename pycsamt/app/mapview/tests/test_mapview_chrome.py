# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for pycsamt.app.mapview.callbacks.chrome — seed/theme/help/dock."""

from __future__ import annotations

import numpy as np
import pytest

pytest.importorskip("dash", reason="dash required")
pytest.importorskip("dash_bootstrap_components", reason="dbc required")

from pycsamt.app.mapview.callbacks import chrome as chrome_mod
from pycsamt.map._core import MapData, StationRecord
from pycsamt.map.view import MapView


def _capture():
    captured: dict = {}

    class _App:
        def callback(self, *a, **k):
            def deco(fn):
                captured[fn.__name__] = fn
                return fn

            return deco

        def clientside_callback(self, *a, **k):
            captured.setdefault("_clientside_calls", 0)
            captured["_clientside_calls"] += 1

    chrome_mod.register_chrome(_App())
    return captured


def _map_data() -> MapData:
    class _Z:
        def __init__(self):
            self.freq = [10.0, 1.0]
            self.resistivity = np.ones((2, 2, 2)) * 100.0
            self.phase = np.ones((2, 2, 2)) * 45.0

    class _Edi:
        def __init__(self, station):
            self.station = station
            self.Z = _Z()

    return MapData(
        sites=[_Edi("S00"), _Edi("S01")],
        stations=(
            StationRecord("S00", 1.0, 2.0, 10.0, "L1", 0),
            StationRecord("S01", 1.1, 2.1, 20.0, "L1", 1),
        ),
        metadata={},
    )


class TestRegisterChrome:
    def test_register_chrome_is_callable(self):
        from pycsamt.app.mapview.callbacks.chrome import register_chrome

        assert callable(register_chrome)

    def test_expected_outputs_wired(self):
        from pycsamt.app.mapview._ids import IDs
        from pycsamt.app.mapview.app import create_app

        app = create_app()
        cb_outputs = str(app.callback_map)
        assert IDs.STORE_DATA in cb_outputs
        assert IDs.STORE_THEME in cb_outputs
        assert IDs.MODAL_HELP in cb_outputs

    def test_resize_snippet_dispatches_resize_event(self):
        from pycsamt.app.mapview.callbacks.chrome import _RESIZE

        assert "dispatchEvent" in _RESIZE
        assert "resize" in _RESIZE

    def test_registers_expected_callbacks_and_three_clientside(self):
        captured = _capture()
        assert {"adopt_seed", "toggle_theme", "toggle_help"} <= set(captured)
        assert captured["_clientside_calls"] == 3


# ---------------------------------------------------------------------------
# _register_seed / adopt_seed
# ---------------------------------------------------------------------------


def test_adopt_seed_no_seed_returns_no_update(monkeypatch):
    captured = _capture()
    monkeypatch.setattr(chrome_mod, "take_seed", lambda: None)
    store, badge, cls = captured["adopt_seed"]("sess-1")
    assert store is chrome_mod.no_update
    assert badge is chrome_mod.no_update
    assert cls is chrome_mod.no_update


def test_adopt_seed_no_session_id_returns_no_update_even_with_seed(monkeypatch):
    captured = _capture()
    view = MapView(_map_data())
    monkeypatch.setattr(chrome_mod, "take_seed", lambda: view)
    store, badge, cls = captured["adopt_seed"](None)
    assert store is chrome_mod.no_update
    assert badge is chrome_mod.no_update
    assert cls is chrome_mod.no_update


def test_adopt_seed_with_seed_and_session_writes_store_and_badge(monkeypatch):
    captured = _capture()
    view = MapView(_map_data())
    monkeypatch.setattr(chrome_mod, "take_seed", lambda: view)
    seen: dict = {}
    monkeypatch.setattr(
        chrome_mod, "set_view", lambda sid, v: seen.update(session_id=sid, view=v)
    )
    store, badge, cls = captured["adopt_seed"]("sess-1")
    assert seen == {"session_id": "sess-1", "view": view}
    assert store["n_stations"] == 2
    assert "2 stations" in badge
    assert "1 line(s)" in badge
    assert cls == "mv-data-badge visible"


# ---------------------------------------------------------------------------
# _register_theme / toggle_theme
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "current, new_theme, icon",
    [
        ("light", "dark", "bi bi-sun"),
        ("dark", "light", "bi bi-moon-stars"),
        (None, "dark", "bi bi-sun"),
        ("", "dark", "bi bi-sun"),
    ],
)
def test_toggle_theme_flips_and_sets_icon(current, new_theme, icon):
    captured = _capture()
    theme, cls = captured["toggle_theme"](1, current)
    assert theme == new_theme
    assert cls == icon


# ---------------------------------------------------------------------------
# _register_help / toggle_help
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("is_open", [True, False])
def test_toggle_help_flips_state(is_open):
    captured = _capture()
    assert captured["toggle_help"](1, None, is_open) is (not is_open)


def test_toggle_help_close_button_also_flips_state():
    captured = _capture()
    assert captured["toggle_help"](None, 1, False) is True
    assert captured["toggle_help"](None, 1, True) is False
