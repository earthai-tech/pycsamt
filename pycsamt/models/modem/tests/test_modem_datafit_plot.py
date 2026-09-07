"""Tests for PlotDataFit and PlotMisfitMap (models/modem/plot.py).

These use the bundled compact Willy 27-frequency ModEM sample when
present, which carries an observed data file plus numbered
``Modular_NLCG_NNN.dat`` response files — enough for a real
observed-vs-model comparison. All tests skip when the sample is absent.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

matplotlib = pytest.importorskip("matplotlib")
matplotlib.use("Agg")
import matplotlib.figure  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

_SAMPLE = (
    Path(__file__).parents[4]
    / "data"
    / "modem"
    / "willy_27freq_watex_line02_sample"
)

pytestmark = pytest.mark.skipif(
    not (_SAMPLE.exists() and any(_SAMPLE.glob("*.dat"))),
    reason="compact Willy ModEM sample not available",
)


@pytest.fixture(scope="module")
def result():
    from pycsamt.models.modem.results import InversionResult

    return InversionResult(
        _SAMPLE, load_models=False, load_covariance=False
    )


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------


def test_present_components_subset_of_canonical(result):
    from pycsamt.models.modem.plot import _RESP_COMPS, _present_components

    comps = _present_components(result.data_obs)
    assert comps
    assert set(comps).issubset(set(_RESP_COMPS))


def test_station_misfit_matches_pooled_residual(result):
    from pycsamt.models.modem.plot import (
        _present_components,
        _rms,
        _station_misfit,
        _z_residuals,
        _collect_z_rows,
    )

    comps = _present_components(result.data_obs)
    name = result.data_obs.site_names[0]
    rms, per_comp, n = _station_misfit(
        result.data_obs, result.data_pred, name, comps
    )
    # recompute independently
    pooled = []
    for c in comps:
        o = _collect_z_rows(result.data_obs, name, c)
        p = _collect_z_rows(
            result.data_pred, name, c, filter_masked=False
        )
        r = _z_residuals(o, p)
        if r.size:
            pooled.append(r)
    ref = _rms(np.concatenate(pooled)) if pooled else float("nan")
    assert np.isclose(rms, ref, equal_nan=True)
    assert n >= 0
    assert set(per_comp) == set(comps)


# ---------------------------------------------------------------------------
# PlotDataFit
# ---------------------------------------------------------------------------


def test_datafit_returns_figure(result):
    from pycsamt.models.modem.plot import PlotDataFit

    fig = PlotDataFit(result=result, max_stations=2).plot()
    assert isinstance(fig, matplotlib.figure.Figure)
    plt.close(fig)


def test_datafit_drops_absent_components(result):
    from pycsamt.models.modem.plot import PlotDataFit, _present_components

    n_present = len(_present_components(result.data_obs))
    fig = PlotDataFit(result=result, stations=[result.data_obs.site_names[0]]).plot()
    # 2 rows (rho, phase) per present component, single station
    assert len(fig.axes) == 2 * n_present
    plt.close(fig)


def test_datafit_component_override(result):
    from pycsamt.models.modem.plot import PlotDataFit

    fig = PlotDataFit(
        result=result,
        stations=[result.data_obs.site_names[0]],
        components=["ZXY"],
    ).plot()
    assert len(fig.axes) == 2
    plt.close(fig)


def test_datafit_accepts_data_objects_without_result(result):
    from pycsamt.models.modem.plot import PlotDataFit

    fig = PlotDataFit(
        data_obs=result.data_obs,
        data_pred=result.data_pred,
        max_stations=2,
    ).plot()
    assert isinstance(fig, matplotlib.figure.Figure)
    plt.close(fig)


def test_datafit_stacked_results_have_band_labels(result):
    from pycsamt.models.modem.plot import PlotDataFit

    fig = PlotDataFit(
        results=[result, result],
        stations=[result.data_obs.site_names[0]],
    ).plot()
    texts = {t.get_text() for t in fig.texts} | {
        a.get_text()
        for ax in fig.axes
        for a in ax.texts
    }
    assert "a)" in texts and "b)" in texts
    plt.close(fig)


def test_datafit_no_data_raises(tmp_path):
    from pycsamt.models.modem.plot import PlotDataFit
    from pycsamt.models.modem.results import InversionResult

    r = InversionResult(tmp_path)
    with pytest.raises(ValueError, match="data_obs"):
        PlotDataFit(result=r).plot()


# ---------------------------------------------------------------------------
# PlotMisfitMap
# ---------------------------------------------------------------------------


def test_misfit_map_returns_figure(result):
    from pycsamt.models.modem.plot import PlotMisfitMap

    fig = PlotMisfitMap(result=result).plot()
    assert isinstance(fig, matplotlib.figure.Figure)
    plt.close(fig)


def test_misfit_map_overall_rms_in_title_matches_log(result):
    from pycsamt.models.modem.plot import PlotMisfitMap

    fig = PlotMisfitMap(result=result).plot()
    blob = " ".join(
        [t.get_text() for t in fig.texts]
        + [a.get_text() for ax in fig.axes for a in ax.texts]
        + [ax.get_title() for ax in fig.axes]
    )
    assert "overall RMS" in blob
    val = float(blob.split("overall RMS")[1].split()[0])
    if np.isfinite(result.final_rms):
        assert abs(val - result.final_rms) < 0.15
    plt.close(fig)


def test_misfit_map_by_component(result):
    from pycsamt.models.modem.plot import PlotMisfitMap, _present_components

    n = len(_present_components(result.data_obs))
    fig = PlotMisfitMap(result=result, by_component=True).plot()
    # one visible axes per component (+ its colorbar axes)
    visible = [ax for ax in fig.axes if ax.get_visible() and ax.has_data()]
    assert len(visible) >= n
    plt.close(fig)


def test_misfit_map_line_labels_present(result):
    from pycsamt.models.modem.plot import PlotMisfitMap, _survey_line_of

    lines = {_survey_line_of(n) for n in result.data_obs.site_names}
    fig = PlotMisfitMap(result=result, show_line_labels=True).plot()
    drawn = {
        a.get_text().lstrip("L")
        for ax in fig.axes
        for a in ax.texts
    }
    if len(lines) >= 2:
        assert lines & {d.lstrip("L") for d in drawn}
    plt.close(fig)


def test_misfit_map_rotation_switches_to_local_km(result):
    from pycsamt.models.modem.plot import PlotMisfitMap

    fig = PlotMisfitMap(result=result, rotate_deg=-90, aspect="auto").plot()
    ax = fig.axes[0]
    assert "km" in ax.get_xlabel()
    plt.close(fig)


def test_misfit_map_requires_predicted(tmp_path):
    from pycsamt.models.modem.plot import PlotMisfitMap
    from pycsamt.models.modem.results import InversionResult

    r = InversionResult(tmp_path)
    with pytest.raises(ValueError):
        PlotMisfitMap(result=r).plot()
