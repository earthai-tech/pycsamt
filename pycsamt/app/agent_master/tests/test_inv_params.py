# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for pycsamt.app.agent_master.callbacks.inv_params."""

from __future__ import annotations

import pytest

pytest.importorskip("dash", reason="dash required")
pytest.importorskip("dash_bootstrap_components", reason="dbc required")


class TestIsPinnOrHybrid:
    def test_pinn_keywords_match(self):
        from pycsamt.app.agent_master.callbacks.inv_params import (
            _is_pinn_or_hybrid,
        )

        assert _is_pinn_or_hybrid("run a PINN inversion")
        assert _is_pinn_or_hybrid("physics-informed neural network")
        assert _is_pinn_or_hybrid("no training data required")

    def test_hybrid_keywords_match(self):
        from pycsamt.app.agent_master.callbacks.inv_params import (
            _is_pinn_or_hybrid,
        )

        assert _is_pinn_or_hybrid("use the hybrid mode")
        assert _is_pinn_or_hybrid("two-stage warm start")

    def test_unrelated_text_does_not_match(self):
        from pycsamt.app.agent_master.callbacks.inv_params import (
            _is_pinn_or_hybrid,
        )

        assert not _is_pinn_or_hybrid("run qc on the data")

    def test_case_insensitive(self):
        from pycsamt.app.agent_master.callbacks.inv_params import (
            _is_pinn_or_hybrid,
        )

        assert _is_pinn_or_hybrid("PHYSICS INFORMED model")


class TestRegisterInvParams:
    def test_register_is_callable(self):
        from pycsamt.app.agent_master.callbacks.inv_params import (
            register_inv_params,
        )

        assert callable(register_inv_params)

    def test_expected_outputs_wired(self):
        from pycsamt.app.agent_master._ids import IDs
        from pycsamt.app.agent_master.app import create_app

        app = create_app()
        cb_outputs = str(app.callback_map)
        assert IDs.STORE_INV_CONFIG in cb_outputs
        assert IDs.INV_PANEL_2D in cb_outputs
        assert IDs.INV_PANEL_HYBRID in cb_outputs


def _unwrap(entry):
    fn = entry["callback"]
    return getattr(fn, "__wrapped__", fn)


def _find(agent_app, input_id, output_hint):
    matches = [
        k
        for k, entry in agent_app.callback_map.items()
        if entry["inputs"]
        and entry["inputs"][0]["id"] == input_id
        and output_hint in k
    ]
    assert len(matches) == 1, (input_id, output_hint, matches)
    return _unwrap(agent_app.callback_map[matches[0]])


class TestToggleDimPanels:
    def _fn(self, agent_app):
        from pycsamt.app.agent_master._ids import IDs

        return _find(agent_app, IDs.INV_DIM, IDs.INV_PANEL_2D)

    def test_dim_1_hides_both_panels(self, agent_app):
        fn = self._fn(agent_app)
        d2, d3 = fn(1)
        assert d2 == {"display": "none"}
        assert d3 == {"display": "none"}

    def test_dim_2_shows_2d_only(self, agent_app):
        fn = self._fn(agent_app)
        d2, d3 = fn(2)
        assert d2 == {"display": "block"}
        assert d3 == {"display": "none"}

    def test_dim_3_shows_both(self, agent_app):
        fn = self._fn(agent_app)
        d2, d3 = fn(3)
        assert d2 == {"display": "block"}
        assert d3 == {"display": "block"}


class TestToggleHybridPanel:
    def _fn(self, agent_app):
        from pycsamt.app.agent_master._ids import IDs

        return _find(agent_app, IDs.INV_MODE, IDs.INV_PANEL_HYBRID)

    def test_hybrid_mode_shows_panel(self, agent_app):
        fn = self._fn(agent_app)
        assert fn("hybrid") == {"display": "block"}

    def test_other_mode_hides_panel(self, agent_app):
        fn = self._fn(agent_app)
        assert fn("pinn") == {"display": "none"}


class TestConfirmInvParams:
    def _fn(self, agent_app):
        from pycsamt.app.agent_master._ids import IDs

        return _find(agent_app, IDs.BTN_INV_CONFIRM, IDs.STORE_INV_CONFIG)

    def test_prevent_update_without_click(self, agent_app):
        from dash.exceptions import PreventUpdate

        fn = self._fn(agent_app)
        with pytest.raises(PreventUpdate):
            fn(None, "pinn", 1, "mt1d", 10, 2000, 500, 0.01, 0.01, 0.005, 0.005, 5000, "")

    def test_validation_errors_returned(self, agent_app):
        from dash import no_update

        fn = self._fn(agent_app)
        cfg, msg, is_open = fn(
            1, "pinn", 1, "mt1d", 1, 0, 0, 0, 0.01, 0.005, 0.005, 5000, ""
        )
        assert cfg is no_update
        assert is_open is no_update
        text = str(msg)
        assert "n_layers must be >= 2" in text
        assert "depth_max must be positive" in text
        assert "epochs must be >= 1" in text
        assert "Learning rate must be > 0" in text

    def test_valid_params_saved(self, agent_app):
        fn = self._fn(agent_app)
        cfg, msg, is_open = fn(
            1, "hybrid", 3, "mt3d", 12, 3000.0, 800, 0.02,
            0.02, 0.01, 0.01, 6000.0, "ckpt.pt",
        )
        assert cfg == {
            "mode": "hybrid",
            "dim": 3,
            "solver": "mt3d",
            "n_layers": 12,
            "depth_max": 3000.0,
            "epochs": 800,
            "lr": 0.02,
            "smoothness_weight": 0.02,
            "lateral_weight": 0.01,
            "graph_weight": 0.01,
            "radius": 6000.0,
            "checkpoint": "ckpt.pt",
        }
        assert is_open is False
        assert "Parameters saved" in str(msg)

    def test_valid_params_apply_defaults(self, agent_app):
        fn = self._fn(agent_app)
        cfg, _msg, _is_open = fn(
            1, None, None, None, 5, 100.0, 10, 0.001,
            None, None, None, None, None,
        )
        assert cfg["mode"] == "pinn"
        assert cfg["dim"] == 1
        assert cfg["solver"] == "mt1d"
        assert cfg["smoothness_weight"] == 0.01
        assert cfg["lateral_weight"] == 0.005
        assert cfg["graph_weight"] == 0.005
        assert cfg["radius"] == 5000.0
        assert cfg["checkpoint"] == ""


class TestLoadInvForm:
    def _fn(self, agent_app):
        from pycsamt.app.agent_master._ids import IDs

        return _find(agent_app, IDs.CANVAS_INV_PARAMS, IDs.INV_MODE)

    def test_prevent_update_when_closed(self, agent_app):
        from dash.exceptions import PreventUpdate

        fn = self._fn(agent_app)
        with pytest.raises(PreventUpdate):
            fn(False, {"mode": "pinn"})

    def test_prevent_update_without_config(self, agent_app):
        from dash.exceptions import PreventUpdate

        fn = self._fn(agent_app)
        with pytest.raises(PreventUpdate):
            fn(True, None)

    def test_loads_config_values(self, agent_app):
        fn = self._fn(agent_app)
        cfg = {
            "mode": "hybrid",
            "dim": 2,
            "solver": "mt2d",
            "n_layers": 8,
            "depth_max": 1500.0,
            "epochs": 300,
            "lr": 0.005,
            "smoothness_weight": 0.02,
            "lateral_weight": 0.01,
            "graph_weight": 0.01,
            "radius": 4000.0,
            "checkpoint": "c.pt",
        }
        result = fn(True, cfg)
        assert result == (
            "hybrid", 2, "mt2d", 8, 1500.0, 300, 0.005,
            0.02, 0.01, 0.01, 4000.0, "c.pt",
        )

    def test_loads_defaults_for_missing_keys(self, agent_app):
        fn = self._fn(agent_app)
        result = fn(True, {"unrelated": True})
        assert result == (
            "pinn", 1, "mt1d", 10, 2000.0, 500, 0.01,
            0.01, 0.005, 0.005, 5000.0, "",
        )
