from __future__ import annotations

import pytest

from pycsamt.api.view.config import (
    PYCSAMT_API_VIEW,
    APIViewConfig,
    configure_api_view,
    reset_api_view,
)


@pytest.fixture
def cfg() -> APIViewConfig:
    return APIViewConfig()


def test_default_backend_is_pycsamt(cfg):
    assert cfg.backend == "pycsamt"
    assert cfg.enabled() is True


def test_configure_backend_true_false_and_none(cfg):
    cfg.configure(backend=True)
    assert cfg.backend == "pycsamt"
    cfg.configure(backend=False)
    assert cfg.backend == "pandas"
    assert cfg.enabled() is False
    cfg.configure(backend=None)
    assert cfg.backend == "pandas"


def test_configure_backend_string_variants(cfg):
    for name in ("api", "apiframe", "pycsamt", "view", "PYCSAMT"):
        cfg.configure(backend=name)
        assert cfg.backend == "pycsamt"
    for name in ("false", "none", "off", "pandas", "raw", "RAW"):
        cfg.configure(backend=name)
        assert cfg.backend == "pandas"


def test_configure_backend_rejects_unknown_string(cfg):
    with pytest.raises(ValueError):
        cfg.configure(backend="bogus")


def test_configure_via_enabled_alias(cfg):
    cfg.configure(enabled=False)
    assert cfg.backend == "pandas"


def test_configure_via_wrapper_alias(cfg):
    def custom(data, **meta):
        return data

    cfg.configure(wrapper=custom)
    assert cfg.backend is custom


def test_configure_rejects_unknown_keys(cfg):
    with pytest.raises(ValueError):
        cfg.configure(bogus_key=1)


def test_context_restores_backend_after_block(cfg):
    cfg.configure(backend=False)
    with cfg.context(backend=True) as ctx:
        assert ctx is cfg
        assert cfg.backend == "pycsamt"
    assert cfg.backend == "pandas"


def test_context_restores_backend_on_exception(cfg):
    cfg.configure(backend=False)
    with pytest.raises(RuntimeError):
        with cfg.context(backend=True):
            raise RuntimeError("boom")
    assert cfg.backend == "pandas"


def test_context_with_no_overrides(cfg):
    with cfg.context():
        assert cfg.backend == "pycsamt"


def test_reset_restores_default(cfg, monkeypatch):
    monkeypatch.delenv("PYCSAMT_API_VIEW", raising=False)
    cfg.configure(backend=False)
    cfg.reset()
    assert cfg.backend == "pycsamt"


def test_wrap_frame_pandas_backend_returns_raw(cfg):
    import pandas as pd

    cfg.configure(backend=False)
    df = pd.DataFrame({"a": [1]})
    assert cfg.wrap_frame(df) is df


def test_wrap_frame_callable_backend_calls_wrapper(cfg):
    calls = []

    def custom(data, **meta):
        calls.append(meta)
        return "wrapped"

    cfg.configure(wrapper=custom)
    result = cfg.wrap_frame({"a": 1}, name="x")
    assert result == "wrapped"
    assert calls[0]["name"] == "x"


def test_wrap_frame_pycsamt_backend_uses_default_wrap(cfg):
    import pandas as pd

    from pycsamt.api.view.frame import APIFrame

    df = pd.DataFrame({"a": [1]})
    result = cfg.wrap_frame(df, name="demo")
    assert isinstance(result, APIFrame)
    assert result.name == "demo"


def test_summary_and_repr_string_backend(cfg):
    text = cfg.summary()
    assert text == "APIViewConfig(backend='pycsamt')"
    assert repr(cfg) == text


def test_summary_and_repr_callable_backend_named_function(cfg):
    def my_wrapper(data, **meta):
        return data

    cfg.configure(wrapper=my_wrapper)
    text = cfg.summary()
    assert "my_wrapper" in text
    assert repr(cfg) == text


def test_summary_callable_backend_without_name_uses_class_name(cfg):
    class CallableWrapper:
        def __call__(self, data, **meta):
            return data

    instance = CallableWrapper()
    cfg.configure(wrapper=instance)
    text = cfg.summary()
    assert "CallableWrapper" in text


def test_module_level_configure_and_reset_singleton():
    try:
        configure_api_view(backend=False)
        assert PYCSAMT_API_VIEW.backend == "pandas"
        reset_api_view()
        assert PYCSAMT_API_VIEW.backend == "pycsamt"
    finally:
        reset_api_view()
