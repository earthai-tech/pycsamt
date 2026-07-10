"""Phase 1 tests — ModEmConfig, validation predicates, and _source layout."""

from pathlib import Path

import pytest

_SOURCE = Path(__file__).parents[1] / "_source"
_EXAMPLES_2D = Path(__file__).parents[4] / "ModEMv626" / "ModEM" / "examples" / "2D_MT" / "BLOCK2"
_EXAMPLES_3D = Path(__file__).parents[4] / "ModEMv626" / "ModEM" / "examples" / "3D_MT" / "BLOCK2"

pytestmark_examples = pytest.mark.skipif(
    not _EXAMPLES_2D.exists(),
    reason="ModEMv626 example data not available",
)


# ---------------------------------------------------------------------------
# ModEmConfig
# ---------------------------------------------------------------------------

class TestModEmConfig:
    def test_default_mode(self):
        from pycsamt.models.modem.config import ModEmConfig
        cfg = ModEmConfig()
        assert cfg.mode == "3d"

    def test_is_3d_true(self):
        from pycsamt.models.modem.config import ModEmConfig
        assert ModEmConfig(mode="3d").is_3d is True

    def test_is_3d_false(self):
        from pycsamt.models.modem.config import ModEmConfig
        assert ModEmConfig(mode="2d").is_3d is False

    def test_binary_name_3d(self):
        from pycsamt.models.modem.config import ModEmConfig
        assert ModEmConfig(mode="3d").binary_name == "Mod3DMT"

    def test_binary_name_2d(self):
        from pycsamt.models.modem.config import ModEmConfig
        assert ModEmConfig(mode="2d").binary_name == "Mod2DMT"

    def test_custom_error_floor(self):
        from pycsamt.models.modem.config import ModEmConfig
        cfg = ModEmConfig(error_floor_z=0.1)
        assert cfg.error_floor_z == pytest.approx(0.1)

    def test_custom_initial_rho(self):
        from pycsamt.models.modem.config import ModEmConfig
        cfg = ModEmConfig(initial_rho=300.0)
        assert cfg.initial_rho == pytest.approx(300.0)

    def test_mpi_defaults(self):
        from pycsamt.models.modem.config import ModEmConfig
        cfg = ModEmConfig()
        assert cfg.use_mpi is False
        assert cfg.n_procs == 4
        assert cfg.mpi_command == "mpirun"

    def test_file_name_defaults(self):
        from pycsamt.models.modem.config import ModEmConfig
        cfg = ModEmConfig()
        assert cfg.data_file.endswith(".dat")
        assert cfg.model_file.endswith(".rho")
        assert cfg.covariance_file.endswith(".cov")
        assert cfg.control_file.endswith(".inv")

    def test_write_and_read_python_template(self, tmp_path):
        from pycsamt.models.modem.config import ModEmConfig
        path = ModEmConfig.write_template(tmp_path / "modem_config.py")
        text = path.read_text()
        assert "CONFIG = {" in text
        assert "Workflow dimensionality" in text
        cfg = ModEmConfig.from_file(path)
        assert cfg.mode == "3d"

    def test_write_and_read_json_template(self, tmp_path):
        from pycsamt.models.modem.config import ModEmConfig
        path = ModEmConfig.write_template(tmp_path / "modem_config.json")
        text = path.read_text()
        assert '"_schema"' in text
        assert '"config"' in text
        cfg = ModEmConfig.from_file(path)
        assert cfg.binary_name == "Mod3DMT"

    def test_write_and_read_yaml_template(self, tmp_path):
        from pycsamt.models.modem.config import ModEmConfig
        path = ModEmConfig.write_template(tmp_path / "modem_config.yml")
        text = path.read_text()
        assert "# ---- Dimensionality ----" in text
        cfg = ModEmConfig.from_file(path)
        assert cfg.n_procs == 4

    def test_to_template_uses_instance_values(self, tmp_path):
        from pycsamt.models.modem.config import ModEmConfig
        cfg = ModEmConfig(mode="2d", n_procs=12)
        path = cfg.to_template(tmp_path / "custom.py")
        loaded = ModEmConfig.read(path)
        assert loaded.mode == "2d"
        assert loaded.n_procs == 12

    def test_from_file_rejects_unknown_keys(self, tmp_path):
        from pycsamt.models.modem.config import ModEmConfig
        path = tmp_path / "bad.py"
        path.write_text("CONFIG = {'mode': '3d', 'bad_key': 1}\n")
        with pytest.raises(ValueError, match="bad_key"):
            ModEmConfig.from_file(path)


# ---------------------------------------------------------------------------
# Validation — against real example files
# ---------------------------------------------------------------------------

@pytestmark_examples
class TestValidationWithExamples:
    def test_detect_2d_data_file(self):
        from pycsamt.models.modem.validation import (
            is_data_file,
        )
        assert is_data_file(_EXAMPLES_2D / "d0_JT.dat") is True

    def test_detect_3d_data_file(self):
        from pycsamt.models.modem.validation import (
            is_data_file,
        )
        assert is_data_file(_EXAMPLES_3D / "d0.dat") is True

    def test_detect_2d_model_file(self):
        from pycsamt.models.modem.validation import (
            is_model_2d_file,
            is_model_3d_file,
        )
        assert is_model_2d_file(_EXAMPLES_2D / "m0.rho") is True
        assert is_model_3d_file(_EXAMPLES_2D / "m0.rho") is False

    def test_detect_3d_model_file(self):
        from pycsamt.models.modem.validation import (
            is_model_2d_file,
            is_model_3d_file,
        )
        assert is_model_3d_file(_EXAMPLES_3D / "m0.ws") is True
        assert is_model_2d_file(_EXAMPLES_3D / "m0.ws") is False

    def test_detect_covariance_file(self):
        from pycsamt.models.modem.validation import (
            is_covariance_file,
        )
        assert is_covariance_file(_EXAMPLES_3D / "example.cov") is True

    def test_detect_control_file(self):
        from pycsamt.models.modem.validation import (
            is_control_file,
        )
        assert is_control_file(_EXAMPLES_3D / "block2.inv") is True

    def test_detect_log_file_2d(self):
        from pycsamt.models.modem.validation import (
            is_log_file,
        )
        assert is_log_file(_EXAMPLES_2D / "Modular_NLCG.log") is True

    def test_data_is_not_model(self):
        from pycsamt.models.modem.validation import (
            is_model_file,
        )
        assert is_model_file(_EXAMPLES_2D / "d0_JT.dat") is False

    def test_model_is_not_data(self):
        from pycsamt.models.modem.validation import (
            is_data_file,
        )
        assert is_data_file(_EXAMPLES_2D / "m0.rho") is False

    def test_detect_file_type_data(self):
        from pycsamt.models.modem.validation import (
            ModEmFileType,
            detect_file_type,
        )
        assert detect_file_type(_EXAMPLES_3D / "d0.dat") == ModEmFileType.DATA

    def test_detect_file_type_model_3d(self):
        from pycsamt.models.modem.validation import (
            ModEmFileType,
            detect_file_type,
        )
        assert detect_file_type(_EXAMPLES_3D / "m0.ws") == ModEmFileType.MODEL_3D

    def test_detect_file_type_model_2d(self):
        from pycsamt.models.modem.validation import (
            ModEmFileType,
            detect_file_type,
        )
        assert detect_file_type(_EXAMPLES_2D / "m0.rho") == ModEmFileType.MODEL_2D

    def test_detect_file_type_covariance(self):
        from pycsamt.models.modem.validation import (
            ModEmFileType,
            detect_file_type,
        )
        assert detect_file_type(_EXAMPLES_3D / "example.cov") == ModEmFileType.COVARIANCE

    def test_detect_file_type_control(self):
        from pycsamt.models.modem.validation import (
            ModEmFileType,
            detect_file_type,
        )
        assert detect_file_type(_EXAMPLES_3D / "block2.inv") == ModEmFileType.CONTROL

    def test_detect_file_type_log(self):
        from pycsamt.models.modem.validation import (
            ModEmFileType,
            detect_file_type,
        )
        assert detect_file_type(_EXAMPLES_2D / "Modular_NLCG.log") == ModEmFileType.LOG


class TestValidationMissingFile:
    def test_missing_file_data(self):
        from pycsamt.models.modem.validation import (
            is_data_file,
        )
        assert is_data_file("/no/such/file.dat") is False

    def test_missing_file_model(self):
        from pycsamt.models.modem.validation import (
            is_model_file,
        )
        assert is_model_file("/no/such/file.rho") is False

    def test_missing_file_covariance(self):
        from pycsamt.models.modem.validation import (
            is_covariance_file,
        )
        assert is_covariance_file("/no/such/file.cov") is False

    def test_unknown_type(self):
        import os
        import tempfile

        from pycsamt.models.modem.validation import (
            ModEmFileType,
            detect_file_type,
        )
        with tempfile.NamedTemporaryFile(suffix=".txt", delete=False) as f:
            f.write(b"hello world\n")
            name = f.name
        try:
            assert detect_file_type(name) == ModEmFileType.UNKNOWN
        finally:
            os.unlink(name)


# ---------------------------------------------------------------------------
# _source layout
# ---------------------------------------------------------------------------

class TestSourceLayout:
    def test_source_2d_makefile_exists(self):
        assert (_SOURCE / "2D" / "Makefile").exists()

    def test_source_3d_makefile_exists(self):
        assert (_SOURCE / "3D" / "Makefile").exists()

    def test_install_script_exists(self):
        assert (_SOURCE / "install_sources.sh").exists()

    def test_readme_exists(self):
        assert (_SOURCE / "README.txt").exists()

    def test_makefile_2d_has_gfortran(self):
        content = (_SOURCE / "2D" / "Makefile").read_text()
        assert "gfortran" in content

    def test_makefile_3d_has_mpi_target(self):
        content = (_SOURCE / "3D" / "Makefile").read_text()
        assert "mpi" in content and "mpif90" in content
