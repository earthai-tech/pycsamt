"""Tests for OccamModel.read()."""

import math
from pathlib import Path

import pytest

DATA_DIR = Path(__file__).parents[4] / "data" / "occam2D"
MODEL_FILE = DATA_DIR / "Occam2DModel"

pytestmark = pytest.mark.skipif(
    not MODEL_FILE.exists(), reason="sample occam2D data not available"
)


@pytest.fixture(scope="module")
def model():
    from pycsamt.models.occam2d.model import OccamModel

    return OccamModel.read(MODEL_FILE)


def test_model_format(model):
    assert model.format_str == "OCCAM2MTMOD_1.0"


def test_model_name(model):
    assert "MTPY" in model.name.upper() or "MODEL" in model.name.upper()


def test_model_mesh_file(model):
    assert model.mesh_file == "Occam2DMesh"


def test_model_mesh_type(model):
    assert model.mesh_type == "PW2D"


def test_model_binding_offset(model):
    assert math.isclose(model.binding_offset, 0.0, abs_tol=1e-9)


def test_model_n_layers(model):
    assert model.n_layers == 26


def test_model_layers_length(model):
    assert len(model.layers) == 26


def test_model_layer0_n_merge(model):
    assert model.layers[0]["n_merge"] == 2


def test_model_layer0_n_cols(model):
    assert model.layers[0]["n_cols"] == 283


def test_model_layer0_params_length(model):
    assert len(model.layers[0]["params"]) == 283


def test_model_layer0_boundary_edges(model):
    # First and last column of every layer should be a boundary cell (code 7)
    params = model.layers[0]["params"]
    assert params[0] == 7
    assert params[-1] == 7


def test_model_last_layer_n_merge(model):
    assert model.layers[-1]["n_merge"] == 5


def test_model_last_layer_n_cols(model):
    assert model.layers[-1]["n_cols"] == 3


def test_model_n_params(model):
    # Sum of n_cols across all layers = startup Param Count
    assert model.n_params == 3752


def test_model_n_free_params(model):
    # Non-boundary cells
    assert model.n_free_params == model.n_params - 52  # 52 boundary cells


def test_model_n_exceptions(model):
    assert model.n_exceptions == 0


def test_model_statics_file(model):
    assert "none" in model.statics_file.lower()


def test_model_prejudice_file(model):
    assert "none" in model.prejudice_file.lower()


# -----------------------------------------------------------------------
# Defaults
# -----------------------------------------------------------------------


def test_model_defaults():
    from pycsamt.models.occam2d.model import OccamModel

    m = OccamModel()
    assert m.n_layers == 0
    assert m.n_params == 0
    assert m.layers == []


# -----------------------------------------------------------------------
# Error handling
# -----------------------------------------------------------------------


def test_missing_file_raises():
    from pycsamt.models.occam2d.model import OccamModel

    with pytest.raises(FileNotFoundError):
        OccamModel.read("/nonexistent/Occam2DModel")


def test_wrong_format_raises(tmp_path):
    from pycsamt.models.occam2d.model import OccamModel

    bad = tmp_path / "BadModel"
    bad.write_text("Format: WRONG_FORMAT\nModel Name: test\n")
    with pytest.raises(ValueError, match="OCCAM2MTMOD"):
        OccamModel.read(bad)


# -----------------------------------------------------------------------
# Roughness-penalty exceptions
# -----------------------------------------------------------------------


def _small_model():
    from pycsamt.models.occam2d.model import OccamModel

    m = OccamModel()
    m.n_layers = 2
    m.layers = [
        {"n_merge": 1, "n_cols": 2, "params": [2, 2]},
        {"n_merge": 1, "n_cols": 2, "params": [2, 2]},
    ]
    return m


def test_model_exceptions_default_empty():
    from pycsamt.models.occam2d.model import OccamModel

    assert OccamModel().exceptions == []


def test_model_write_no_exceptions_writes_zero_count(tmp_path):
    model = _small_model()
    path = model.write(tmp_path / "Occam2DModel")
    text = path.read_text()
    assert "NO. EXCEPTIONS:   0" in text


def test_model_exceptions_roundtrip(tmp_path):
    from pycsamt.models.occam2d.model import OccamModel

    model = _small_model()
    model.exceptions = [(1, 2, 0.3), (2, 4, 1.0)]
    path = model.write(tmp_path / "Occam2DModel")

    restored = OccamModel.read(path)
    assert restored.n_exceptions == -2
    assert len(restored.exceptions) == 2
    for (bi, bj, expen), (rbi, rbj, rexpen) in zip(
        model.exceptions, restored.exceptions
    ):
        assert (bi, bj) == (rbi, rbj)
        assert math.isclose(expen, rexpen, rel_tol=1e-6)


def test_model_exceptions_written_after_layer_blocks(tmp_path):
    model = _small_model()
    model.exceptions = [(1, 3, 0.0)]
    path = model.write(tmp_path / "Occam2DModel")
    text = path.read_text()
    header_index = text.index("NO. EXCEPTIONS:")
    row_index = text.index("1 3 0")
    assert row_index > header_index


def test_model_exceptions_rejects_self_link(tmp_path):
    model = _small_model()
    model.exceptions = [(1, 1, 0.5)]
    with pytest.raises(ValueError, match="itself"):
        model.write(tmp_path / "Occam2DModel")


def test_model_exceptions_rejects_negative_expen(tmp_path):
    model = _small_model()
    model.exceptions = [(1, 2, -0.5)]
    with pytest.raises(ValueError, match="non-negative"):
        model.write(tmp_path / "Occam2DModel")
