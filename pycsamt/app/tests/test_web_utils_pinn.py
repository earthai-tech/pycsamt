# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Unit tests for pycsamt.app.web.utils_pinn (PINN/Hybrid web bridge)."""

from __future__ import annotations

import base64
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

import numpy as np
import pandas as pd
import pytest


# ── decode_npz_checkpoint ────────────────────────────────────────────────────


def test_decode_npz_checkpoint_writes_and_returns_path():
    from pycsamt.app.web.utils_pinn import decode_npz_checkpoint

    raw = b"not a real npz, just some bytes \x00\x01\x02"
    content_b64 = "data:application/octet-stream;base64," + base64.b64encode(
        raw
    ).decode("ascii")

    path = decode_npz_checkpoint(content_b64)
    try:
        assert path.endswith(".npz")
        with open(path, "rb") as fh:
            assert fh.read() == raw
    finally:
        import os

        os.remove(path)


# ── session_to_obs_1d / session_to_obs_2d ───────────────────────────────────


def test_session_to_obs_1d_raises_when_no_session_data():
    from pycsamt.app.web.utils_pinn import session_to_obs_1d

    with patch("pycsamt.app.web.cache.cache_get", return_value=None):
        with pytest.raises(ValueError, match="No data loaded"):
            session_to_obs_1d("sess-1")


def test_session_to_obs_1d_filters_by_stations():
    from pycsamt.app.web.utils_pinn import session_to_obs_1d

    obs_all = [
        SimpleNamespace(name="S1"),
        SimpleNamespace(name="S2"),
        SimpleNamespace(name="S3"),
    ]
    with patch("pycsamt.app.web.cache.cache_get", return_value=object()):
        with patch(
            "pycsamt.ai.inversion._sites_bridge.sites_to_obs_1d",
            return_value=obs_all,
        ) as mock_bridge:
            obs = session_to_obs_1d(
                "sess-1", stations=["S2"], comp="yx"
            )

    assert [o.name for o in obs] == ["S2"]
    assert mock_bridge.call_args.kwargs["comp"] == "yx"


def test_session_to_obs_1d_raises_when_filter_matches_none():
    from pycsamt.app.web.utils_pinn import session_to_obs_1d

    obs_all = [SimpleNamespace(name="S1")]
    with patch("pycsamt.app.web.cache.cache_get", return_value=object()):
        with patch(
            "pycsamt.ai.inversion._sites_bridge.sites_to_obs_1d",
            return_value=obs_all,
        ):
            with pytest.raises(ValueError, match="None of the selected"):
                session_to_obs_1d("sess-1", stations=["ZZZ"])


def test_session_to_obs_2d_raises_when_no_session_data():
    from pycsamt.app.web.utils_pinn import session_to_obs_2d

    with patch("pycsamt.app.web.cache.cache_get", return_value=None):
        with pytest.raises(ValueError, match="No data loaded"):
            session_to_obs_2d("sess-1")


def test_session_to_obs_2d_filters_by_stations():
    from pycsamt.app.web.utils_pinn import session_to_obs_2d

    obs_all = [SimpleNamespace(name="A"), SimpleNamespace(name="B")]
    with patch("pycsamt.app.web.cache.cache_get", return_value=object()):
        with patch(
            "pycsamt.ai.inversion._sites_bridge.sites_to_obs_2d",
            return_value=obs_all,
        ) as mock_bridge:
            obs = session_to_obs_2d(
                "sess-1", stations=["B"], comp_te="xy", comp_tm="yx"
            )

    assert [o.name for o in obs] == ["B"]
    assert mock_bridge.call_args.kwargs["comp_te"] == "xy"
    assert mock_bridge.call_args.kwargs["comp_tm"] == "yx"


# ── build_pinn_inv ───────────────────────────────────────────────────────────


def test_build_pinn_inv_1d_passes_expected_kwargs():
    from pycsamt.app.web.utils_pinn import build_pinn_inv

    obs = [SimpleNamespace(name="S1")]
    sentinel = MagicMock(name="PINNInverter1D-instance")
    with patch(
        "pycsamt.ai.inversion.PINNInverter1D", return_value=sentinel
    ) as ctor:
        result = build_pinn_inv(obs, "1d", solver="mt1d", comp="yx")

    assert result is sentinel
    ctor.assert_called_once()
    args, kwargs = ctor.call_args
    assert args[0] is obs
    assert kwargs["solver"] == "mt1d"
    assert kwargs["comp"] == "yx"


def test_build_pinn_inv_2d_passes_expected_kwargs():
    from pycsamt.app.web.utils_pinn import build_pinn_inv

    obs = [SimpleNamespace(name="S1"), SimpleNamespace(name="S2")]
    sentinel = MagicMock(name="PINNInverter2D-instance")
    with patch(
        "pycsamt.ai.inversion.PINNInverter2D", return_value=sentinel
    ) as ctor:
        result = build_pinn_inv(obs, "2d", mode="tm", epochs=50)

    assert result is sentinel
    kwargs = ctor.call_args.kwargs
    assert kwargs["mode"] == "tm"
    assert kwargs["epochs"] == 50


def test_build_pinn_inv_3d_builds_station_coords():
    from pycsamt.app.web.utils_pinn import build_pinn_inv

    obs = [SimpleNamespace(name=f"S{i}") for i in range(4)]
    sentinel = MagicMock(name="PINNInverter3D-instance")
    with patch(
        "pycsamt.ai.inversion.PINNInverter3D", return_value=sentinel
    ) as ctor:
        result = build_pinn_inv(obs, "3d", station_spacing=250.0)

    assert result is sentinel
    coords = ctor.call_args.kwargs["station_coords"]
    assert coords.shape == (4, 2)
    np.testing.assert_allclose(coords[:, 0], [0.0, 250.0, 500.0, 750.0])
    np.testing.assert_allclose(coords[:, 1], [0.0, 0.0, 0.0, 0.0])


def test_build_pinn_inv_device_auto_becomes_none():
    from pycsamt.app.web.utils_pinn import build_pinn_inv

    obs = [SimpleNamespace(name="S1")]
    with patch(
        "pycsamt.ai.inversion.PINNInverter1D", return_value=MagicMock()
    ) as ctor:
        build_pinn_inv(obs, "1d", device="auto")

    assert ctor.call_args.kwargs["device"] is None


# ── build_hybrid_inv ─────────────────────────────────────────────────────────


def test_build_hybrid_inv_1d_passes_checkpoint_and_kwargs():
    from pycsamt.app.web.utils_pinn import build_hybrid_inv

    obs = [SimpleNamespace(name="S1")]
    sentinel = MagicMock(name="HybridInverter1D-instance")
    with patch(
        "pycsamt.ai.inversion.HybridInverter1D", return_value=sentinel
    ) as ctor:
        result = build_hybrid_inv(
            obs, "1d", "/tmp/ckpt.npz", max_iter=99
        )

    assert result is sentinel
    args, kwargs = ctor.call_args
    assert args[0] is obs
    assert args[1] == "/tmp/ckpt.npz"
    assert kwargs["max_iter"] == 99


def test_build_hybrid_inv_2d_passes_expected_kwargs():
    from pycsamt.app.web.utils_pinn import build_hybrid_inv

    obs = [SimpleNamespace(name="S1")]
    checkpoint = MagicMock(name="fitted-inverter")
    sentinel = MagicMock(name="HybridInverter2D-instance")
    with patch(
        "pycsamt.ai.inversion.HybridInverter2D", return_value=sentinel
    ) as ctor:
        result = build_hybrid_inv(obs, "2d", checkpoint, n_layers=7)

    assert result is sentinel
    args, kwargs = ctor.call_args
    assert args[1] is checkpoint
    assert kwargs["n_layers"] == 7
    # epochs is aliased from max_iter for the 2-D/3-D constructors
    assert kwargs["epochs"] == 150


def test_build_hybrid_inv_device_auto_becomes_none():
    from pycsamt.app.web.utils_pinn import build_hybrid_inv

    obs = [SimpleNamespace(name="S1")]
    with patch(
        "pycsamt.ai.inversion.HybridInverter1D", return_value=MagicMock()
    ) as ctor:
        build_hybrid_inv(obs, "1d", "ckpt", device="auto")

    assert ctor.call_args.kwargs["device"] is None


def test_build_hybrid_inv_3d_builds_station_coords():
    from pycsamt.app.web.utils_pinn import build_hybrid_inv

    obs = [SimpleNamespace(name=f"S{i}") for i in range(3)]
    sentinel = MagicMock(name="HybridInverter3D-instance")
    with patch(
        "pycsamt.ai.inversion.HybridInverter3D", return_value=sentinel
    ) as ctor:
        result = build_hybrid_inv(
            obs, "3d", "ckpt", station_spacing=100.0
        )

    assert result is sentinel
    coords = ctor.call_args.kwargs["station_coords"]
    assert coords.shape == (3, 2)
    np.testing.assert_allclose(coords[:, 0], [0.0, 100.0, 200.0])


# ── _extract_loss_history ────────────────────────────────────────────────────


def test_extract_loss_history_prefers_convergence_curve():
    from pycsamt.app.web.utils_pinn import _extract_loss_history

    inv = SimpleNamespace(
        convergence_curve=lambda: pd.DataFrame({"loss": [3.0, 2.0, 1.0]}),
    )
    assert _extract_loss_history(inv) == [3.0, 2.0, 1.0]


def test_extract_loss_history_falls_back_to_loss_curves():
    from pycsamt.app.web.utils_pinn import _extract_loss_history

    df = pd.DataFrame(
        {
            "epoch": [1, 1, 2, 2],
            "loss": [4.0, 2.0, 3.0, 1.0],
        }
    )
    inv = SimpleNamespace(loss_curves=lambda: df)
    history = _extract_loss_history(inv)
    assert history == [3.0, 2.0]


def test_extract_loss_history_falls_back_to_history_attr():
    from pycsamt.app.web.utils_pinn import _extract_loss_history

    inv = SimpleNamespace(_history=[5.0, 4.0, 3.0])
    assert _extract_loss_history(inv) == [5.0, 4.0, 3.0]


def test_extract_loss_history_returns_empty_when_nothing_available():
    from pycsamt.app.web.utils_pinn import _extract_loss_history

    assert _extract_loss_history(SimpleNamespace()) == []


# ── plot_pinn_convergence ────────────────────────────────────────────────────


def test_plot_pinn_convergence_with_history_returns_data_uri():
    from pycsamt.app.web.utils_pinn import plot_pinn_convergence

    inv = SimpleNamespace(_history=[3.0, 2.0, 1.0])
    src = plot_pinn_convergence(inv, dark=True, label="Hybrid")

    assert src.startswith("data:image/png;base64,")
    base64.b64decode(src.split(",", 1)[1])  # must decode cleanly


def test_plot_pinn_convergence_without_history_returns_placeholder():
    from pycsamt.app.web.utils_pinn import plot_pinn_convergence

    src = plot_pinn_convergence(SimpleNamespace(), dark=False)

    assert src.startswith("data:image/png;base64,")


# ── plot_pinn_data_fit ───────────────────────────────────────────────────────


def _residuals_df():
    return pd.DataFrame(
        {
            "station": ["S1", "S1", "S2"],
            "freq": [10.0, 1.0, 10.0],
            "rho_obs": [100.0, 120.0, 90.0],
            "rho_pred": [98.0, 118.0, 91.0],
            "phase_obs": [45.0, 40.0, 44.0],
            "phase_pred": [46.0, 39.0, 43.0],
        }
    )


def test_plot_pinn_data_fit_with_df_returns_data_uri():
    from pycsamt.app.web.utils_pinn import plot_pinn_data_fit

    src = plot_pinn_data_fit(
        None, "S1", dark=True, df=_residuals_df()
    )

    assert src.startswith("data:image/png;base64,")
    base64.b64decode(src.split(",", 1)[1])


def test_plot_pinn_data_fit_empty_station_returns_placeholder():
    from pycsamt.app.web.utils_pinn import plot_pinn_data_fit

    src = plot_pinn_data_fit(
        None, "NOPE", dark=False, df=_residuals_df()
    )

    assert src.startswith("data:image/png;base64,")


def test_plot_pinn_data_fit_residuals_error_returns_placeholder():
    from pycsamt.app.web.utils_pinn import plot_pinn_data_fit

    def _boom():
        raise RuntimeError("no residuals")

    inv = SimpleNamespace(residuals=_boom)
    src = plot_pinn_data_fit(inv, "S1", dark=True)

    assert src.startswith("data:image/png;base64,")


# ── plot_pinn_section ────────────────────────────────────────────────────────


def _model_1d(n_layers=3):
    return SimpleNamespace(
        resistivity=np.linspace(50.0, 500.0, n_layers),
        thickness=np.full(n_layers, 100.0),
    )


def test_plot_pinn_section_1d():
    from pycsamt.app.web.utils_pinn import plot_pinn_section

    inv = SimpleNamespace(
        stations=["S1", "S2"],
        predict=lambda: [_model_1d(), _model_1d()],
    )
    src = plot_pinn_section(inv, "1d", dark=True)

    assert src.startswith("data:image/png;base64,")
    base64.b64decode(src.split(",", 1)[1])


def test_plot_pinn_section_2d():
    from pycsamt.app.web.utils_pinn import plot_pinn_section

    mat = np.random.default_rng(0).uniform(1.0, 3.0, size=(5, 4))
    inv = SimpleNamespace(
        stations=["S1", "S2", "S3", "S4"],
        resistivity_section=lambda as_log10=True: mat,
    )
    src = plot_pinn_section(inv, "2d", dark=False)

    assert src.startswith("data:image/png;base64,")


def test_plot_pinn_section_3d():
    from pycsamt.app.web.utils_pinn import plot_pinn_section

    mat = np.random.default_rng(0).uniform(1.0, 3.0, size=(5, 3))
    inv = SimpleNamespace(
        stations=["S1", "S2", "S3"],
        resistivity_volume=lambda as_log10=True: mat,
    )
    src = plot_pinn_section(inv, "3d", dark=True)

    assert src.startswith("data:image/png;base64,")


# ── pinn_stats_div ───────────────────────────────────────────────────────────


def test_pinn_stats_div_basic_rows():
    from pycsamt.app.web.utils_pinn import pinn_stats_div

    inv = SimpleNamespace(
        n_sites=12, n_layers=8, depth_max=2000.0, epochs=300, lr=0.01
    )
    div = pinn_stats_div(inv, "1d", elapsed_s=12.3, label="PINN")

    table = div.children
    rows = table.children
    # 7 base rows, no "Final loss" row since there's no loss history
    assert len(rows) == 7
    keys = [row.children[0].children for row in rows]
    assert keys == [
        "Track",
        "Stations",
        "Layers",
        "Depth max (m)",
        "Epochs",
        "Learning rate",
        "Elapsed (s)",
    ]


def test_pinn_stats_div_includes_final_loss_when_history_present():
    from pycsamt.app.web.utils_pinn import pinn_stats_div

    inv = SimpleNamespace(
        n_sites=5,
        n_layers=4,
        depth_max=1000.0,
        max_iter=150,
        lr=0.005,
        _history=[2.0, 1.5, 1.0],
    )
    div = pinn_stats_div(inv, "2d", label="Hybrid")

    rows = div.children.children
    assert len(rows) == 8
    assert rows[-1].children[0].children == "Final loss"
    assert rows[-1].children[1].children == f"{1.0:.4e}"
