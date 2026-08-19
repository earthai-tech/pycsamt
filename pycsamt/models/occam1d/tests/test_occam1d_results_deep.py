"""Deep contracts for Occam1D run aggregation and visualization."""

import json
import shutil

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pytest

from pycsamt.models.occam1d import (
    Occam1DResult,
    Occam1DResultFiles,
    Occam1DResultState,
)
from pycsamt.models.occam1d.plot import PlotModel
from pycsamt.models.occam1d.tests.test_occam1d_results import _make_run

matplotlib.use("Agg")


def test_result_inventory_uses_signatures_not_required_names(tmp_path):
    run = tmp_path / "run"
    _make_run(run)
    (run / "Occam1DData").rename(run / "sounding.native")
    (run / "Occam1DModel").rename(run / "layers.native")
    result = Occam1DResult(run)
    assert result.files.data.name == "sounding.native"
    assert result.files.model.name == "layers.native"


def test_result_exposes_immutable_provenance_and_state(tmp_path):
    run = tmp_path / "run"
    _make_run(run)
    result = Occam1DResult(run)
    assert result.result_state is Occam1DResultState.READY
    assert result.is_ready
    assert isinstance(result.files, Occam1DResultFiles)
    assert result.available_iterations == (1, 2)
    with pytest.raises(Exception):
        result.files.data = run / "other"


def test_result_rejects_ambiguous_required_file_in_strict_mode(tmp_path):
    run = tmp_path / "run"
    _make_run(run)
    shutil.copy2(run / "Occam1DData", run / "data.backup")
    with pytest.raises(ValueError, match="Multiple Occam1D data"):
        Occam1DResult(run)


def test_non_strict_result_records_ambiguous_file_warning(tmp_path):
    run = tmp_path / "run"
    _make_run(run)
    shutil.copy2(run / "Occam1DData", run / "data.backup")
    result = Occam1DResult(run, strict=False)
    assert result.is_ready
    assert any("Multiple Occam1D data" in item for item in result.warnings)


def test_result_rejects_unavailable_iteration_with_inventory(tmp_path):
    run = tmp_path / "run"
    _make_run(run)
    with pytest.raises(ValueError, match="available iterations: 1, 2"):
        Occam1DResult(run, iteration=9)


def test_result_can_load_startup_only_build(tmp_path):
    run = tmp_path / "run"
    _make_run(run)
    for path in run.glob("ITER*"):
        path.unlink()
    (run / "RESP_02.resp").unlink()
    result = Occam1DResult(run)
    assert result.iteration == 0
    assert result.iteration_file is None
    assert result.response is None
    assert np.isnan(result.final_rms)


def test_metadata_export_is_json_safe_and_has_file_provenance(tmp_path):
    run = tmp_path / "run"
    _make_run(run)
    result = Occam1DResult(run)
    target = result.export_metadata(tmp_path / "metadata.json")
    values = json.loads(target.read_text(encoding="utf8"))
    assert values["result_state"] == "ready"
    assert values["files"]["iteration"].endswith("ITER_02")
    assert values["selected_iteration"] == 2


def test_export_prefix_is_sanitized(tmp_path):
    run = tmp_path / "run"
    _make_run(run)
    outputs = Occam1DResult(run).export_text(
        tmp_path / "text",
        prefix="site 01/unsafe",
    )
    assert all(path.name.startswith("site_01_unsafe") for path in outputs.values())


def test_plot_model_uses_inversion_station_marker(tmp_path):
    run = tmp_path / "run"
    _make_run(run)
    figure = Occam1DResult(run).plot_model()
    try:
        collections = figure.axes[0].collections
        assert collections
        face = collections[-1].get_facecolor()[0, :3]
        np.testing.assert_allclose(face, 0.0)
    finally:
        plt.close(figure)


def test_plot_validation_rejects_invalid_configuration(tmp_path):
    run = tmp_path / "run"
    _make_run(run)
    result = Occam1DResult(run)
    with pytest.raises(ValueError, match="dpi"):
        PlotModel(result, dpi=0)
    with pytest.raises(ValueError, match="figsize"):
        PlotModel(result, figsize=(4,))
    with pytest.raises(ValueError, match="depth_max"):
        result.plot_model(depth_max=-1)


def test_response_plot_requires_matching_response(tmp_path):
    run = tmp_path / "run"
    _make_run(run)
    result = Occam1DResult(run, iteration=1)
    assert result.response is None
    with pytest.raises(RuntimeError, match="No response"):
        result.plot_response()
