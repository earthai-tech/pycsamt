from __future__ import annotations

import pytest

from pycsamt.api.view.progress import (
    PYCSAMT_PROGRESS,
    ProgressAPIConfig,
    ProgressConfig,
    configure_progress,
    get_progress_bar,
    iter_progress,
    progress_enabled,
    reset_progress,
)


@pytest.fixture(autouse=True)
def _reset():
    reset_progress()
    yield
    reset_progress()


@pytest.fixture
def cfg() -> ProgressAPIConfig:
    return ProgressAPIConfig()


def test_progress_config_dataclass_defaults():
    pc = ProgressConfig()
    assert pc.enabled == "auto"
    assert pc.unit == "it"
    assert pc.leave is False


def test_defaults(cfg):
    assert cfg.force_verbose is None
    assert cfg.style == "ascii"
    assert cfg.min_interval == 0.1
    assert cfg.log_every is None


def test_configure_rejects_bad_style(cfg):
    with pytest.raises(ValueError):
        cfg.configure(style="fancy")


def test_configure_rejects_unknown_key(cfg):
    with pytest.raises(AttributeError):
        cfg.configure(not_a_key=1)


def test_configure_sets_valid_keys(cfg):
    cfg.configure(style="unicode", min_interval=0.5, log_every=10)
    assert cfg.style == "unicode"
    assert cfg.min_interval == 0.5
    assert cfg.log_every == 10


def test_context_restores_after_block(cfg):
    cfg.configure(style="unicode")
    with cfg.context(style="ascii") as ctx:
        assert ctx is cfg
        assert cfg.style == "ascii"
    assert cfg.style == "unicode"


def test_context_restores_on_exception(cfg):
    cfg.configure(style="unicode")
    with pytest.raises(RuntimeError):
        with cfg.context(style="ascii"):
            raise RuntimeError("boom")
    assert cfg.style == "unicode"


def test_reset_restores_defaults(cfg):
    cfg.configure(style="unicode", min_interval=0.9)
    cfg.reset()
    assert cfg.style == "ascii"
    assert cfg.min_interval == 0.1


def test_resolve_force_verbose_overrides_call_site(cfg):
    cfg.configure(force_verbose=0)
    assert cfg.resolve(True) == 0
    assert cfg.resolve(2) == 0


def test_resolve_uses_default_verbose_when_none(cfg):
    cfg.configure(default_verbose=2)
    assert cfg.resolve(None) == 2


def test_summary_and_repr_list_all_keys(cfg):
    text = cfg.summary()
    assert "ProgressAPIConfig" in text
    assert "style" in text
    assert "min_interval" in text
    assert repr(cfg) == text


def test_module_level_configure_and_reset():
    configure_progress(style="unicode")
    assert PYCSAMT_PROGRESS.style == "unicode"
    reset_progress()
    assert PYCSAMT_PROGRESS.style == "ascii"


# ─────────────────────────────────────────────────────────────────────────
# get_progress_bar
# ─────────────────────────────────────────────────────────────────────────


def test_get_progress_bar_respects_forced_silence(capsys):
    configure_progress(force_verbose=0)
    with get_progress_bar(total=2, desc="x", verbose=True) as bar:
        for _ in range(2):
            bar.update(1)
    assert capsys.readouterr().out == ""


def test_get_progress_bar_uses_log_every_override(capsys):
    configure_progress(force_verbose=2)
    with get_progress_bar(total=1, desc="x", verbose=True, log_every=1) as bar:
        bar.update(1)


# ─────────────────────────────────────────────────────────────────────────
# progress_enabled
# ─────────────────────────────────────────────────────────────────────────


def test_progress_enabled_bool_passthrough():
    assert progress_enabled(True) is True
    assert progress_enabled(False) is False


@pytest.mark.parametrize("value", ["1", "true", "yes", "on", "TRUE"])
def test_progress_enabled_truthy_strings(value):
    assert progress_enabled(value) is True


@pytest.mark.parametrize("value", ["0", "false", "no", "off", "none", "FALSE"])
def test_progress_enabled_falsy_strings(value):
    assert progress_enabled(value) is False


def test_progress_enabled_auto_checks_stderr_isatty(monkeypatch):
    import sys

    monkeypatch.setattr(sys.stderr, "isatty", lambda: True, raising=False)
    assert progress_enabled("auto") is True
    monkeypatch.setattr(sys.stderr, "isatty", lambda: False, raising=False)
    assert progress_enabled("auto") is False


def test_progress_enabled_rejects_unknown_string():
    with pytest.raises(ValueError):
        progress_enabled("bogus")


# ─────────────────────────────────────────────────────────────────────────
# iter_progress
# ─────────────────────────────────────────────────────────────────────────


def test_iter_progress_disabled_yields_plain_items():
    items = list(iter_progress([1, 2, 3], enabled=False))
    assert items == [1, 2, 3]


def test_iter_progress_enabled_wraps_with_progress_bar(capsys):
    items = list(
        iter_progress([1, 2, 3], enabled=True, desc="test", total=3)
    )
    assert items == [1, 2, 3]
