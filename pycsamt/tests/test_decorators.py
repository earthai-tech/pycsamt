from __future__ import annotations

import subprocess

import numpy as np
import pandas as pd
import pytest

from pycsamt.decorators import (
    Deprecated,
    GdalDataCheck,
    ReplaceBy,
    check_empty,
    ensure_fit,
    has_fit,
    isdf,
    noop,
)

# --------------------------------------------------------------------------
# has_fit
# --------------------------------------------------------------------------


def test_has_fit_noop_when_fit_already_defined():
    @has_fit()
    class WithFit:
        def fit(self, x):
            return x

    assert WithFit().fit(3) == 3


def test_has_fit_idempotent_double_decoration():
    @has_fit()
    class Loader:
        def read(self, src):
            return f"read:{src}"

    # Applying again should be a no-op due to the _has_fit_alias guard.
    Loader2 = has_fit()(Loader)
    assert Loader2 is Loader
    assert Loader().fit("a") == "read:a"


def test_has_fit_wraps_read_as_fit():
    @has_fit()
    class Loader:
        def read(self, src, **kw):
            """Read docstring."""
            return f"reading {src!r}"

    loader = Loader()
    assert loader.fit("dataset.avg") == "reading 'dataset.avg'"
    assert loader.fit.__name__ == "read"


def test_has_fit_raise_when_neither_fit_nor_read():
    with pytest.raises(AttributeError, match="neither 'fit' nor 'read'"):

        @has_fit("raise")
        class Empty:
            pass


def test_has_fit_warn_installs_noop_and_warns():
    with pytest.warns(RuntimeWarning, match="neither 'fit' nor 'read'"):

        @has_fit("warn")
        class Empty:
            pass

    assert Empty().fit(1, x=2) is None


def test_has_fit_ignore_installs_silent_noop():
    @has_fit("ignore")
    class Empty:
        pass

    assert Empty().fit() is None


# --------------------------------------------------------------------------
# noop
# --------------------------------------------------------------------------


def test_noop_returns_object_unchanged():
    def original(x):
        return x + 1

    decorated = noop()(original)
    assert decorated is original

    decorated2 = noop("not implemented yet")(original)
    assert decorated2 is original


# --------------------------------------------------------------------------
# Deprecated
# --------------------------------------------------------------------------


def test_deprecated_requires_nonempty_reason():
    with pytest.raises(ValueError):
        Deprecated("")
    with pytest.raises(ValueError):
        Deprecated("   ")
    with pytest.raises(ValueError):
        Deprecated(None)  # type: ignore[arg-type]


def test_deprecated_wraps_and_warns():
    @Deprecated("use new_func instead")
    def old_func(a, b=2):
        """Original doc."""
        return a + b

    assert "[DEPRECATED] use new_func instead" in old_func.__doc__
    assert "Original doc." in old_func.__doc__

    with pytest.warns(DeprecationWarning, match="old_func is deprecated"):
        result = old_func(1, b=3)
    assert result == 4


def test_deprecated_handles_missing_original_doc():
    @Deprecated("gone")
    def bare(a):
        return a

    assert bare.__doc__.startswith("[DEPRECATED] gone")
    with pytest.warns(DeprecationWarning):
        assert bare(5) == 5


# --------------------------------------------------------------------------
# GdalDataCheck
# --------------------------------------------------------------------------


@pytest.fixture(autouse=True)
def _reset_gdal_check_state():
    GdalDataCheck._checked = False
    GdalDataCheck._available = False
    yield
    GdalDataCheck._checked = False
    GdalDataCheck._available = False


def test_gdal_data_check_uses_existing_env_var(monkeypatch, tmp_path):
    monkeypatch.setenv("GDAL_DATA", str(tmp_path))

    @GdalDataCheck()
    def fn():
        return "ok"

    assert fn() == "ok"
    assert GdalDataCheck._available is True


def test_gdal_data_check_locates_via_gdal_config(monkeypatch, tmp_path):
    monkeypatch.delenv("GDAL_DATA", raising=False)

    def fake_run(*a, **k):
        return subprocess.CompletedProcess(args=a, returncode=0, stdout=str(tmp_path))

    monkeypatch.setattr(subprocess, "run", fake_run)

    @GdalDataCheck()
    def fn():
        return "ok"

    assert fn() == "ok"
    assert GdalDataCheck._available is True
    import os

    assert os.environ.get("GDAL_DATA") == str(tmp_path)


def test_gdal_data_check_missing_raises_when_configured(monkeypatch, tmp_path):
    monkeypatch.delenv("GDAL_DATA", raising=False)
    monkeypatch.delenv("PYCSAMT_DOCS_BUILD", raising=False)

    def fake_run(*a, **k):
        raise FileNotFoundError("no gdal-config")

    monkeypatch.setattr(subprocess, "run", fake_run)

    @GdalDataCheck(raise_on_missing=True)
    def fn():
        return "ok"

    with pytest.raises(ImportError):
        fn()


def test_gdal_data_check_second_instantiation_skips_relookup(monkeypatch, tmp_path):
    monkeypatch.setenv("GDAL_DATA", str(tmp_path))
    GdalDataCheck()
    assert GdalDataCheck._checked is True

    # A second instantiation must not re-run the lookup (branch where
    # `_checked` is already True short-circuits `_locate_gdal_data`).
    def boom():
        raise AssertionError("should not be called again")

    monkeypatch.setattr(GdalDataCheck, "_locate_gdal_data", boom)
    GdalDataCheck()  # no raise


def test_gdal_data_check_config_returns_invalid_dir(monkeypatch):
    monkeypatch.delenv("GDAL_DATA", raising=False)

    def fake_run(*a, **k):
        return subprocess.CompletedProcess(
            args=a, returncode=0, stdout="/not/a/real/dir"
        )

    monkeypatch.setattr(subprocess, "run", fake_run)

    @GdalDataCheck(raise_on_missing=False)
    def fn():
        return "ok"

    assert fn() == "ok"
    assert GdalDataCheck._available is False


def test_gdal_data_check_missing_but_not_raising(monkeypatch):
    monkeypatch.delenv("GDAL_DATA", raising=False)
    monkeypatch.setenv("PYCSAMT_DOCS_BUILD", "1")

    def fake_run(*a, **k):
        raise FileNotFoundError("no gdal-config")

    monkeypatch.setattr(subprocess, "run", fake_run)

    @GdalDataCheck(raise_on_missing=False)
    def fn():
        return "still runs"

    assert fn() == "still runs"
    assert GdalDataCheck._available is False


# --------------------------------------------------------------------------
# ReplaceBy
# --------------------------------------------------------------------------


def test_replace_by_default_reason_and_delegation():
    def new_func(a, b):
        return a * b

    @ReplaceBy(new_func)
    def old_func(a, b):
        raise AssertionError("should not be called")

    with pytest.warns(DeprecationWarning, match="Use new_func instead"):
        assert old_func(3, 4) == 12


def test_replace_by_custom_reason():
    def new_func(a):
        return a + 1

    @ReplaceBy(new_func, reason="deprecated in v2")
    def old_func(a):
        raise AssertionError("should not be called")

    with pytest.warns(DeprecationWarning, match="deprecated in v2"):
        assert old_func(1) == 2


# --------------------------------------------------------------------------
# isdf
# --------------------------------------------------------------------------


def test_isdf_no_parameters_passthrough():
    @isdf
    def fn():
        return "no-params"

    assert fn() == "no-params"


def test_isdf_converts_data_param_to_dataframe():
    @isdf
    def process(data, columns=None):
        return data

    out = process(np.array([[1, 2], [3, 4]]), columns=["a", "b"])
    assert isinstance(out, pd.DataFrame)
    assert list(out.columns) == ["a", "b"]


def test_isdf_columns_as_string_becomes_list():
    @isdf
    def process(data, columns=None):
        return data

    out = process(np.array([[1], [2]]), columns="only")
    assert isinstance(out, pd.DataFrame)
    assert list(out.columns) == ["only"]


def test_isdf_columns_mismatch_ignores_columns():
    @isdf
    def process(data, columns=None):
        return data

    # Two columns of data but only one column name -> mismatch branch.
    out = process(np.array([[1, 2], [3, 4]]), columns=["only"])
    assert isinstance(out, pd.DataFrame)
    assert list(out.columns) == ["0", "1"]


def test_isdf_leaves_existing_dataframe_untouched():
    @isdf
    def process(data):
        return data

    df = pd.DataFrame({"x": [1, 2]})
    out = process(df)
    assert out is df


def test_isdf_none_data_is_untouched():
    @isdf
    def process(data=None):
        return data

    assert process() is None


def test_isdf_first_positional_param_used_when_no_data_name():
    @isdf
    def process(values, columns=None):
        return values

    out = process(np.array([[1, 2]]), columns=["a", "b"])
    assert isinstance(out, pd.DataFrame)


def test_isdf_works_on_bound_method():
    class Proc:
        @isdf
        def process(self, data, columns=None):
            return data

    out = Proc().process(np.array([[1, 2]]), columns=["a", "b"])
    assert isinstance(out, pd.DataFrame)
    assert list(out.columns) == ["a", "b"]


def test_isdf_conversion_failure_raises_value_error():
    @isdf
    def process(data, columns=None):
        return data

    class Unconvertible:
        shape = None

        def __iter__(self):
            raise TypeError("cannot iterate")

    with pytest.raises(ValueError, match="Unable to convert"):
        process(Unconvertible())


# --------------------------------------------------------------------------
# check_empty
# --------------------------------------------------------------------------


def test_check_empty_no_parens_raises_on_empty_positional():
    @check_empty
    def process_data(data, *args, **kwargs):
        return data

    with pytest.raises(ValueError, match="considered empty"):
        process_data([])

    assert process_data([1, 2]) == [1, 2]


def test_check_empty_no_parens_checks_first_kwarg():
    @check_empty
    def process_data(data=None):
        return data

    with pytest.raises(ValueError):
        process_data(data=[])

    assert process_data(data=[1]) == [1]


def test_check_empty_no_parens_no_args_no_kwargs():
    @check_empty
    def process_data():
        return "ok"

    assert process_data() == "ok"


def test_check_empty_no_parens_none_allowed_by_default():
    @check_empty
    def process_data(data=None):
        return data

    assert process_data(None) is None


def test_check_empty_no_parens_non_len_value_is_skipped():
    @check_empty
    def process_data(data):
        return data

    assert process_data(5) == 5


def test_check_empty_with_parens_checks_only_named_params():
    @check_empty(params=["x"])
    def load_data(x, y=None):
        return x, y

    with pytest.raises(ValueError):
        load_data([], y=[])

    assert load_data([1], y=[]) == ([1], [])


def test_check_empty_with_parens_warn_mode():
    @check_empty(params=["x"], error="warn")
    def load_data(x):
        return x

    with pytest.warns(UserWarning, match="considered empty"):
        out = load_data([])
    assert out == []


def test_check_empty_with_parens_ignore_mode():
    @check_empty(params=["x"], error="ignore")
    def load_data(x):
        return x

    assert load_data([]) == []


def test_check_empty_allow_none_false_raises_on_none():
    @check_empty(params=["x"], allow_none=False)
    def load_data(x):
        return x

    with pytest.raises(ValueError):
        load_data(None)


def test_check_empty_none_as_empty_true_raises_on_none():
    @check_empty(params=["x"], none_as_empty=True)
    def load_data(x):
        return x

    with pytest.raises(ValueError):
        load_data(None)


def test_check_empty_with_parens_skips_absent_param_name():
    @check_empty(params=["missing"])
    def load_data(x):
        return x

    # "missing" is not a real parameter of load_data, so the check is
    # skipped entirely and the call proceeds normally.
    assert load_data([]) == []


def test_check_empty_with_parens_no_params_list():
    @check_empty()
    def load_data(x):
        return x

    # params is None -> nothing is checked, even an empty value passes.
    assert load_data([]) == []


# --------------------------------------------------------------------------
# ensure_fit
# --------------------------------------------------------------------------


def test_ensure_fit_noop_when_fit_defined():
    @ensure_fit()
    class WithFit:
        def fit(self, x):
            return x * 2

    assert WithFit().fit(3) == 6


def test_ensure_fit_aliases_read_as_fit():
    @ensure_fit()
    class Loader:
        def read(self, src):
            return f"read:{src}"

    assert Loader().fit("a") == "read:a"


def test_ensure_fit_raise_when_neither():
    with pytest.raises(AttributeError, match="neither 'fit' nor 'read'"):

        @ensure_fit("raise")
        class Empty:
            pass


def test_ensure_fit_warn_installs_noop():
    with pytest.warns(RuntimeWarning, match="neither 'fit' nor 'read'"):

        @ensure_fit("warn")
        class Empty:
            pass

    assert Empty().fit() is None


def test_ensure_fit_ignore_installs_silent_noop():
    @ensure_fit("ignore")
    class Empty:
        pass

    assert Empty().fit(1, 2, x=3) is None
