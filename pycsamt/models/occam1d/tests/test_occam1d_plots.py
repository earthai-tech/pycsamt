"""Visual smoke tests for the Occam1D model-image products."""

import matplotlib
import matplotlib.pyplot as plt

from pycsamt.models.occam1d import Occam1DResult
from pycsamt.models.occam1d.tests.test_occam1d_results import _make_run

matplotlib.use("Agg")


def test_all_main_plots_return_figures(tmp_path):
    run = tmp_path / "run"
    _make_run(run)
    result = Occam1DResult(run)
    figures = [
        result.plot_model(),
        result.plot_response(),
        result.plot_convergence(),
        result.plot_summary(),
    ]
    assert all(hasattr(figure, "savefig") for figure in figures)
    for figure in figures:
        plt.close(figure)


def test_save_main_images_uses_stable_names(tmp_path):
    run = tmp_path / "run"
    _make_run(run)
    paths = Occam1DResult(run).save_main_images(tmp_path / "model-image")
    assert set(paths) == {"model", "response", "convergence", "summary"}
    assert all(path.is_file() and path.stat().st_size > 0 for path in paths.values())
    assert paths["summary"].name == "S01_occam1d_summary.png"


def test_plot_style_customizes_artists_and_hides_model_legend(tmp_path):
    run = tmp_path / "run"
    _make_run(run)
    result = Occam1DResult(run)

    response = result.plot_response(
        style_overrides={
            "observed__marker": "x",
            "predicted__color": "purple",
            "response_legend": False,
        }
    )
    model = result.plot_model(
        style_overrides={"model__color": "green"}
    )

    assert response.axes[0].lines[0].get_marker() == "x"
    assert any(line.get_color() == "purple" for line in response.axes[0].lines)
    assert response.axes[0].get_legend() is None
    assert model.axes[0].lines[0].get_color() == "green"
    assert model.axes[0].get_legend() is None
    plt.close(response)
    plt.close(model)
