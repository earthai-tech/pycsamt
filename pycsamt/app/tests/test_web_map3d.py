"""Regression tests for the web 3D map figure builders."""

from __future__ import annotations

import numpy as np


def test_map3d_scene_uses_light_background_when_requested():
    from pycsamt.app.web.callbacks.map3d import _dark_scene

    layout = _dark_scene(dark=False)

    assert layout["scene"]["bgcolor"] == "#ffffff"
    assert layout["paper_bgcolor"] == "#eff1f5"
    assert layout["scene"]["xaxis"]["tickfont"]["color"] == "#6c6f85"


def test_empty_map3d_figure_uses_transparent_canvas():
    from pycsamt.app.web.pages.map3d import _empty_3d_fig

    fig = _empty_3d_fig(dark=False)

    assert fig.layout.paper_bgcolor == "rgba(0,0,0,0)"
    assert fig.layout.plot_bgcolor == "rgba(0,0,0,0)"
    assert fig.layout.scene.bgcolor == "rgba(0,0,0,0)"
    assert fig.layout.scene.xaxis.visible is False


def test_survey_line_profiles_require_named_lines():
    from pycsamt.app.web.callbacks.map3d import (
        _has_survey_line_profiles,
    )

    assert not _has_survey_line_profiles(None)
    assert not _has_survey_line_profiles({"line_counts": {"L1": 4}})
    assert _has_survey_line_profiles({"line_counts": {"L1": 4, "L2": 5}})


def test_inversion_result_converts_to_linear_map3d_profile():
    from pycsamt.app.web.callbacks.map3d import (
        _profiles_from_inversion_result,
    )
    from pycsamt.inversion.results import InversionResult

    result = InversionResult(
        method="mt",
        dimension="2d",
        backend="builtin",
        model={
            "rho_2d": np.array([[2.0, 2.2], [2.4, 2.6]]),
            "x_centers": np.array([0.0, 500.0]),
            "z_centers": np.array([50.0, 150.0]),
        },
    )

    profiles = _profiles_from_inversion_result(result)
    profile = next(iter(profiles.values()))

    assert profile["x"].tolist() == [0.0, 500.0]
    assert profile["z"].tolist() == [50.0, 150.0]
    np.testing.assert_allclose(profile["rho"], 10.0 ** result.model["rho_2d"])


def test_block_fig_drops_non_finite_cells_and_clamps_iso_band():
    from pycsamt.app.web.callbacks.map3d import (
        _build_block_fig,
    )

    x = np.linspace(-100.0, 300.0, 5)
    y = np.array([0.0, 1000.0])
    z = np.linspace(0.0, 500.0, 4)

    base = np.linspace(10.0, 10_000.0, x.size * y.size * z.size)
    rho = base.reshape(y.size, x.size, z.size)
    rho[0, 0, 0] = np.nan
    rho[1, 3, 2] = np.inf

    fig = _build_block_fig(
        x_arr=x,
        y_arr=y,
        z_arr=z,
        rho_3d=rho,
        cmap="Viridis",
        vmin=1,
        vmax=100_000,
        opacity=0.7,
        isovalue_lo=0.0,
        isovalue_hi=5.0,
        n_surfaces=8,
        clip_x=1.0,
        clip_y=1.0,
        clip_z=1.0,
        contours=True,
        title="",
        rho_lo=1,
        rho_hi=100_000,
        depth_lo=0,
        depth_hi=500,
    )

    assert fig.data
    trace = fig.data[0]
    values = np.asarray(trace.value, dtype=float)

    assert trace.type == "volume"
    assert values.size > 0
    assert np.isfinite(values).all()
    assert trace.isomin >= values.min()
    assert trace.isomax <= values.max()
    assert trace.isomax > trace.isomin


def test_block_fig_clips_axes_from_axis_bounds_not_axis_maximum():
    from pycsamt.app.web.callbacks.map3d import (
        _build_block_fig,
    )

    x = np.linspace(-500.0, 500.0, 6)
    y = np.array([-750.0, 750.0])
    z = np.array([0.0, 100.0, 250.0])
    rho = np.full((y.size, x.size, z.size), 100.0)
    rho[:, x >= 0, :] = 1000.0

    fig = _build_block_fig(
        x_arr=x,
        y_arr=y,
        z_arr=z,
        rho_3d=rho,
        cmap="Viridis",
        vmin=1,
        vmax=10_000,
        opacity=0.8,
        isovalue_lo=1.0,
        isovalue_hi=3.0,
        n_surfaces=4,
        clip_x=0.5,
        clip_y=0.5,
        clip_z=1.0,
        contours=False,
        title="",
        rho_lo=1,
        rho_hi=10_000,
        depth_lo=0,
        depth_hi=250,
    )

    trace = fig.data[0]
    assert max(trace.x) <= 0.0
    assert max(trace.y) <= 0.0
    assert len(trace.value) > 0


def test_block_fig_uses_finite_profile_neighbour_when_other_line_is_missing():
    from pycsamt.app.web.callbacks.map3d import (
        _build_block_fig,
    )

    x = np.array([0.0, 100.0])
    y = np.array([0.0, 1000.0])
    z = np.array([0.0, 250.0])
    rho = np.array(
        [
            [[100.0, 120.0], [200.0, 240.0]],
            [[np.nan, np.nan], [300.0, 360.0]],
        ]
    )

    fig = _build_block_fig(
        x_arr=x,
        y_arr=y,
        z_arr=z,
        rho_3d=rho,
        cmap="Viridis",
        vmin=1,
        vmax=10_000,
        opacity=0.8,
        isovalue_lo=1.0,
        isovalue_hi=3.0,
        n_surfaces=4,
        clip_x=1.0,
        clip_y=1.0,
        clip_z=1.0,
        contours=False,
        title="",
        rho_lo=1,
        rho_hi=10_000,
        depth_lo=0,
        depth_hi=250,
    )

    values = np.asarray(fig.data[0].value, dtype=float)
    assert values.size > rho.size
    assert np.isfinite(values).all()
