from __future__ import annotations

import pytest

from pycsamt.api.pipe.config import (
    PYCSAMT_PIPE,
    PipelineAPIConfig,
    configure_pipe,
    reset_pipe,
)


@pytest.fixture(autouse=True)
def _restore_global_pipe_config():
    """PYCSAMT_PIPE is a process-wide singleton; reset it afterwards so
    configure_pipe calls never leak into other tests."""
    yield
    PYCSAMT_PIPE.reset()


def test_configure_sets_multiple_attributes():
    cfg = PipelineAPIConfig()
    cfg.configure(plot_dpi=300, plot_fmt="pdf")
    assert cfg.plot_dpi == 300
    assert cfg.plot_fmt == "pdf"


def test_configure_rejects_unknown_key():
    cfg = PipelineAPIConfig()
    with pytest.raises(AttributeError, match="Unknown pipe config key"):
        cfg.configure(bogus_key=1)


def test_context_reverts_after_block():
    cfg = PipelineAPIConfig()
    original_dpi = cfg.plot_dpi
    with cfg.context(plot_dpi=600, plot_fmt="svg") as ctx:
        assert ctx is cfg
        assert cfg.plot_dpi == 600
        assert cfg.plot_fmt == "svg"
    assert cfg.plot_dpi == original_dpi
    assert cfg.plot_fmt == "png"


def test_context_reverts_on_exception():
    cfg = PipelineAPIConfig()
    original = cfg.output_root
    with pytest.raises(ValueError):
        with cfg.context(output_root="tmp"):
            raise ValueError("boom")
    assert cfg.output_root == original


def test_reset_restores_defaults():
    cfg = PipelineAPIConfig()
    cfg.configure(plot_dpi=999, show_progress=False)
    cfg.reset()
    assert cfg.plot_dpi == 150
    assert cfg.show_progress is True


def test_clone_returns_independent_deep_copy():
    cfg = PipelineAPIConfig()
    clone = cfg.clone()
    clone.plot_dpi = 999
    assert cfg.plot_dpi == 150
    assert clone is not cfg


def test_repr_lists_every_field():
    cfg = PipelineAPIConfig()
    text = repr(cfg)
    assert text.startswith("PipelineAPIConfig")
    assert "output_root" in text
    assert "plot_dpi" in text


def test_module_level_configure_and_reset_pipe():
    configure_pipe(plot_dpi=222)
    assert PYCSAMT_PIPE.plot_dpi == 222
    reset_pipe()
    assert PYCSAMT_PIPE.plot_dpi == 150
