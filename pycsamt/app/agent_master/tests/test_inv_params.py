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
