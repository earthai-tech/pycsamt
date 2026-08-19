"""Behavioral contract tests for the Occam1D object foundation."""

import io
from unittest.mock import Mock

import pytest

from pycsamt.models.occam1d import Occam1DBase, Occam1DObjectState


def _logger():
    logger = Mock()
    for name in ("debug", "info", "warning", "error", "exception"):
        setattr(logger, name, Mock())
    return logger


def test_new_object_has_explicit_lifecycle_and_diagnostics():
    base = Occam1DBase(logger=_logger(), metadata={"station": "S01"})
    assert base.state is Occam1DObjectState.NEW
    assert base.is_valid and not base.is_bound
    assert base.diagnostics()["metadata"] == {"station": "S01"}


def test_path_binding_and_clearing_updates_lifecycle(tmp_path):
    base = Occam1DBase(logger=_logger())
    bound = base._bind_path(tmp_path / "run" / "Startup")
    assert bound.is_absolute()
    assert base.path == bound
    assert base.state is Occam1DObjectState.BOUND
    base.path = None
    assert base.path is None
    assert base.state is Occam1DObjectState.NEW


def test_input_and_output_path_boundaries(tmp_path):
    base = Occam1DBase(logger=_logger())
    source = tmp_path / "Data"
    source.write_text("Format: EMData_1.1\n", encoding="utf8")
    assert base._require_input_file(source) == source.resolve()
    output = base._prepare_output_file(tmp_path / "nested" / "Model")
    assert output.parent.is_dir()
    with pytest.raises(IsADirectoryError):
        base._prepare_output_file(tmp_path)
    with pytest.raises(FileNotFoundError):
        base._require_input_file(tmp_path / "missing")


def test_verbosity_is_independent_from_logging():
    logger = _logger()
    stream = io.StringIO()
    quiet = Occam1DBase(verbose=0, logger=logger, stream=stream)
    assert quiet._emit("hidden") is False
    assert stream.getvalue() == ""
    logger.info.assert_not_called()

    visible = Occam1DBase(verbose=2, logger=logger, stream=stream)
    assert visible._emit("stage", level=1) is True
    assert visible._emit("detail", level=2) is True
    assert stream.getvalue().splitlines() == ["stage", "detail"]
    logger.info.assert_not_called()


def test_warning_provenance_is_deduplicated_and_logged():
    logger = _logger()
    base = Occam1DBase(logger=logger)
    base.add_warning("phase uncertainty was floored")
    base.add_warning("phase uncertainty was floored")
    assert base.warnings == ("phase uncertainty was floored",)
    logger.warning.assert_called_once()


def test_mark_invalid_preserves_reason_and_diagnostics():
    logger = _logger()
    base = Occam1DBase(logger=logger)
    base.mark_invalid("parameter count mismatch")
    assert base.state is Occam1DObjectState.INVALID
    assert not base.is_valid
    assert base.diagnostics()["state_reason"] == "parameter count mismatch"
    logger.error.assert_called_once()


def test_invalid_state_requires_explicit_recovery(tmp_path):
    base = Occam1DBase(logger=_logger())
    base.mark_invalid("bad format")
    base._bind_path(tmp_path / "Data")
    assert base.state is Occam1DObjectState.INVALID
    assert base.state_reason == "bad format"
    base.mark_valid()
    assert base.state is Occam1DObjectState.BOUND
    assert base.state_reason is None


@pytest.mark.parametrize("verbose", [-1, 1.5, "1", None])
def test_invalid_verbosity_is_rejected(verbose):
    error = ValueError if verbose == -1 else TypeError
    with pytest.raises(error):
        Occam1DBase(verbose=verbose, logger=_logger())


def test_invalid_logger_metadata_stream_and_unknown_kwargs_are_rejected():
    with pytest.raises(TypeError, match="logger"):
        Occam1DBase(logger=object())
    with pytest.raises(TypeError, match="metadata"):
        Occam1DBase(logger=_logger(), metadata=[])
    with pytest.raises(TypeError, match="stream"):
        Occam1DBase(logger=_logger(), stream=object())
    with pytest.raises(TypeError, match="Unexpected"):
        Occam1DBase(logger=_logger(), mystery=True)


def test_metadata_view_is_read_only_and_source_mapping_is_copied():
    source = {"survey": "A"}
    base = Occam1DBase(logger=_logger(), metadata=source)
    source["survey"] = "B"
    assert base.metadata_view["survey"] == "A"
    with pytest.raises(TypeError):
        base.metadata_view["survey"] = "C"
