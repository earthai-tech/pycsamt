"""Phase 4 tests — InversionResult."""

from pathlib import Path

import numpy as np
import pytest

_EX2D = (
    Path(__file__).parents[4] / "ModEMv626" / "ModEM" / "examples" / "2D_MT" / "BLOCK2"
)
_EX3D = (
    Path(__file__).parents[4] / "ModEMv626" / "ModEM" / "examples" / "3D_MT" / "BLOCK2"
)

pytestmark = pytest.mark.skipif(
    not _EX2D.exists(), reason="ModEMv626 example data not available"
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def result_2d():
    from pycsamt.models.modem.results import InversionResult

    return InversionResult(_EX2D)


@pytest.fixture(scope="module")
def result_3d():
    from pycsamt.models.modem.results import InversionResult

    return InversionResult(_EX3D)


# ---------------------------------------------------------------------------
# 2D result
# ---------------------------------------------------------------------------


def test_result_2d_mode(result_2d):
    assert result_2d.mode == "2d"


def test_result_2d_log_loaded(result_2d):
    assert result_2d.log is not None
    assert result_2d.log.n_iter > 0


def test_result_2d_n_iter(result_2d):
    assert result_2d.n_iter > 0


def test_result_2d_final_rms(result_2d):
    assert np.isfinite(result_2d.final_rms)
    assert result_2d.final_rms > 0


def test_result_2d_rms_history(result_2d):
    h = result_2d.rms_history
    assert len(h) == result_2d.n_iter
    assert np.all(np.isfinite(h))


def test_result_2d_best_rms(result_2d):
    assert result_2d.best_rms <= result_2d.final_rms + 1e-9


def test_result_2d_model_initial_loaded(result_2d):
    assert result_2d.model_initial is not None
    assert result_2d.model_initial.nx > 0


def test_result_2d_model_final_loaded(result_2d):
    assert result_2d.model_final is not None
    assert result_2d.model_final.nx > 0


def test_result_2d_model_initial_type(result_2d):
    from pycsamt.models.modem.model2d import ModEmModel2D

    assert isinstance(result_2d.model_initial, ModEmModel2D)


def test_result_2d_models_dict(result_2d):
    assert "m0" in result_2d.models


def test_result_2d_data_obs_loaded(result_2d):
    assert result_2d.data_obs is not None
    assert result_2d.data_obs.n_sites > 0


def test_result_2d_repr(result_2d):
    r = repr(result_2d)
    assert "2d" in r
    assert "final_rms" in r


# ---------------------------------------------------------------------------
# 3D result
# ---------------------------------------------------------------------------


def test_result_3d_mode(result_3d):
    assert result_3d.mode == "3d"


def test_result_3d_log_loaded(result_3d):
    assert result_3d.log is not None


def test_result_3d_model_final_type(result_3d):
    from pycsamt.models.modem.model3d import ModEmModel3D

    assert isinstance(result_3d.model_final, ModEmModel3D)


def test_result_3d_model_shape(result_3d):
    m = result_3d.model_final
    assert m.rho_loge.ndim == 3


def test_result_3d_control_loaded(result_3d):
    assert result_3d.control is not None
    assert result_3d.control.max_iterations > 0


def test_result_3d_covariance_loaded(result_3d):
    assert result_3d.covariance is not None
    assert result_3d.covariance.nz_earth > 0


# ---------------------------------------------------------------------------
# Error handling
# ---------------------------------------------------------------------------


def test_result_missing_workdir_raises():
    from pycsamt.models.modem.results import InversionResult

    with pytest.raises(FileNotFoundError):
        InversionResult("/no/such/dir")


def test_result_empty_workdir(tmp_path):
    from pycsamt.models.modem.results import InversionResult

    r = InversionResult(tmp_path)
    assert r.mode == "unknown"
    assert r.log is None
    assert r.model_initial is None
    assert r.n_iter == 0
    assert not np.isfinite(r.final_rms)
