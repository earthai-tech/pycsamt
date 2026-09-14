from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from pycsamt.api.view.config import reset_api_view
from pycsamt.api.view.frame import (
    APIFrame,
    FrameProfile,
    api_frame,
    maybe_wrap_frame,
    wrap_frame,
)


@pytest.fixture(autouse=True)
def _reset_config():
    reset_api_view()
    yield
    reset_api_view()


# ─────────────────────────────────────────────────────────────────────────
# FrameProfile
# ─────────────────────────────────────────────────────────────────────────


def test_frame_profile_from_frame_basic():
    df = pd.DataFrame({"a": [1.0, None], "b": ["x", "y"]})
    profile = FrameProfile.from_frame(df)
    assert profile.rows == 2
    assert profile.columns == 2
    assert profile.numeric_columns == ("a",)
    assert profile.missing_cells == 1
    assert profile.missing_fraction == 0.25
    assert profile.shape == (2, 2)


def test_frame_profile_empty_dataframe_has_zero_missing_fraction():
    df = pd.DataFrame()
    profile = FrameProfile.from_frame(df)
    assert profile.missing_fraction == 0.0


def test_frame_profile_to_dict():
    df = pd.DataFrame({"a": [1, 2]})
    d = FrameProfile.from_frame(df).to_dict()
    assert d["rows"] == 2
    assert d["column_names"] == ["a"]
    assert d["numeric_columns"] == ["a"]


def test_frame_profile_memory_usage_failure_falls_back_to_zero(monkeypatch):
    df = pd.DataFrame({"a": [1, 2]})

    def _raise(*a, **k):
        raise RuntimeError("boom")

    monkeypatch.setattr(pd.DataFrame, "memory_usage", _raise)
    profile = FrameProfile.from_frame(df)
    assert profile.memory_bytes == 0


# ─────────────────────────────────────────────────────────────────────────
# APIFrame construction / coercion
# ─────────────────────────────────────────────────────────────────────────


def test_api_frame_none_data_builds_empty_dataframe():
    obj = APIFrame(None)
    assert obj.df.empty


def test_api_frame_from_another_api_frame_shares_by_default():
    base = APIFrame({"a": [1]})
    wrapped = APIFrame(base)
    assert wrapped.df is base.df


def test_api_frame_from_another_api_frame_copy_true():
    base = APIFrame({"a": [1]})
    wrapped = APIFrame(base, copy=True)
    assert wrapped.df is not base.df
    assert wrapped.df.equals(base.df)


def test_api_frame_from_dataframe_shares_by_default():
    df = pd.DataFrame({"a": [1]})
    obj = APIFrame(df)
    assert obj.df is df


def test_api_frame_from_dataframe_copy_true():
    df = pd.DataFrame({"a": [1]})
    obj = APIFrame(df, copy=True)
    assert obj.df is not df


def test_api_frame_from_dict_builds_dataframe():
    obj = APIFrame({"a": [1, 2]})
    assert list(obj.df.columns) == ["a"]


def test_api_frame_from_records():
    obj = APIFrame.from_records(
        [{"a": 1, "b": 2}, {"a": 3, "b": 4}], name="records",
    )
    assert obj.name == "records"
    assert obj.df.shape == (2, 2)


def test_api_frame_default_name():
    assert APIFrame({"a": [1]}).name == "dataframe"


# ─────────────────────────────────────────────────────────────────────────
# Properties / dunder protocol
# ─────────────────────────────────────────────────────────────────────────


def test_api_frame_df_setter_coerces():
    obj = APIFrame({"a": [1]})
    obj.df = {"b": [2]}
    assert list(obj.df.columns) == ["b"]


def test_api_frame_data_and_shape_and_columns():
    obj = APIFrame({"a": [1, 2], "b": [3, 4]})
    assert obj.data.shape == (2, 2)
    assert obj.shape == (2, 2)
    assert list(obj.columns) == ["a", "b"]


def test_api_frame_schema():
    obj = APIFrame({"a": [1], "b": ["x"]})
    schema = obj.schema
    assert schema["a"] == "int64"
    assert schema["b"] == "object"


def test_api_frame_len_iter_contains():
    obj = APIFrame({"a": [1, 2, 3]})
    assert len(obj) == 3
    assert list(iter(obj)) == ["a"]
    assert "a" in obj
    assert "z" not in obj


def test_api_frame_getitem():
    obj = APIFrame({"a": [1, 2]})
    assert obj["a"].tolist() == [1, 2]


def test_api_frame_setitem_tuple_key():
    obj = APIFrame({"a": [1, 2]})
    obj[0, "a"] = 99
    assert obj.df.loc[0, "a"] == 99


def test_api_frame_setitem_existing_column_key():
    obj = APIFrame({"a": [1, 2]})
    obj["a"] = [9, 9]
    assert obj["a"].tolist() == [9, 9]


def test_api_frame_setitem_new_column_key():
    obj = APIFrame({"a": [1, 2]})
    obj["b"] = [5, 6]
    assert obj["b"].tolist() == [5, 6]


def test_api_frame_array_protocol():
    obj = APIFrame({"a": [1, 2]})
    arr = np.asarray(obj)
    assert arr.tolist() == [[1], [2]]
    arr2 = np.array(obj, dtype=float)
    assert arr2.dtype == np.float64


def test_api_frame_getattr_column_access():
    obj = APIFrame({"a": [1, 2]})
    assert obj.a.tolist() == [1, 2]


def test_api_frame_getattr_delegates_to_dataframe_method():
    obj = APIFrame({"a": [1, 2]})
    assert obj.sum()["a"] == 3


def test_api_frame_getattr_raises_for_unknown():
    obj = APIFrame({"a": [1]})
    with pytest.raises(AttributeError):
        obj.totally_unknown_attr
    with pytest.raises(AttributeError):
        obj._totally_private


def test_api_frame_dir_includes_columns():
    obj = APIFrame({"my_col": [1]})
    assert "my_col" in dir(obj)


def test_api_frame_repr_truncates_many_columns():
    data = {f"c{i}": [1] for i in range(8)}
    obj = APIFrame(data, name="wide")
    text = repr(obj)
    assert "..." in text
    assert "wide" in text


def test_api_frame_repr_no_truncation_for_few_columns():
    obj = APIFrame({"a": [1]}, name="narrow")
    text = repr(obj)
    assert "..." not in text


# ─────────────────────────────────────────────────────────────────────────
# summary / profile / missing / numeric_stats / conversions
# ─────────────────────────────────────────────────────────────────────────


def test_summary_includes_kind_source_and_description():
    obj = APIFrame(
        {"a": [1, None], "b": [2, 3]},
        name="s",
        kind="k",
        source="src",
        description="desc",
    )
    text = obj.summary(max_columns=1)
    assert "APIFrame: s" in text
    assert "kind: k" in text
    assert "..." in text  # truncated column list
    assert "source: src" in text
    assert "description: desc" in text


def test_summary_no_columns_shows_dash():
    obj = APIFrame(None)
    text = obj.summary()
    assert "columns: -" in text


def test_str_uses_summary():
    obj = APIFrame({"a": [1]}, name="x")
    assert str(obj) == obj.summary()


def test_profile_matches_stats():
    obj = APIFrame({"a": [1]})
    assert obj.profile() == obj.stats


def test_missing_counts():
    obj = APIFrame({"a": [1, None, 3]})
    assert obj.missing()["a"] == 1


def test_numeric_stats_delegates_to_describe():
    obj = APIFrame({"a": [1, 2, 3]})
    desc = obj.numeric_stats()
    assert "mean" in desc.index


def test_to_pandas_copy_and_share():
    obj = APIFrame({"a": [1]})
    assert obj.to_pandas() is obj.df
    copied = obj.to_pandas(copy=True)
    assert copied is not obj.df


def test_to_numpy():
    obj = APIFrame({"a": [1, 2]})
    assert obj.to_numpy().tolist() == [[1], [2]]


def test_to_dict_default_orient_list():
    obj = APIFrame({"a": [1, 2]})
    assert obj.to_dict() == {"a": [1, 2]}


def test_to_dict_explicit_kwargs():
    obj = APIFrame({"a": [1]})
    assert obj.to_dict("records") == [{"a": 1}]


def test_copy_preserves_metadata():
    obj = APIFrame(
        {"a": [1]}, name="n", kind="k", source="s",
        units={"a": "m"}, meta={"x": 1}, description="d",
    )
    clone = obj.copy()
    assert clone is not obj
    assert clone.df is not obj.df
    assert clone.name == "n"
    assert clone.units == {"a": "m"}
    assert clone.meta == {"x": 1}
    assert clone.description == "d"


def test_with_df_overrides():
    obj = APIFrame({"a": [1]}, name="n")
    other = obj.with_df({"b": [2]}, name="override")
    assert other.name == "override"
    assert list(other.df.columns) == ["b"]


def test_update_meta_and_set_units_return_self():
    obj = APIFrame({"a": [1]})
    result = obj.update_meta(x=1).set_units(a="ohm.m")
    assert result is obj
    assert obj.meta == {"x": 1}
    assert obj.units == {"a": "ohm.m"}


# ─────────────────────────────────────────────────────────────────────────
# wrap_frame / maybe_wrap_frame / api_frame decorator
# ─────────────────────────────────────────────────────────────────────────


def test_wrap_frame_uses_global_config():
    result = wrap_frame({"a": [1]}, name="global")
    assert isinstance(result, APIFrame)
    assert result.name == "global"


def test_maybe_wrap_frame_api_false_returns_raw():
    df = pd.DataFrame({"a": [1]})
    assert maybe_wrap_frame(df, api=False) is df


def test_maybe_wrap_frame_api_true_forces_wrap_even_when_disabled():
    from pycsamt.api.view.config import PYCSAMT_API_VIEW

    PYCSAMT_API_VIEW.configure(backend=False)
    try:
        df = pd.DataFrame({"a": [1]})
        result = maybe_wrap_frame(df, api=True, name="forced")
        assert isinstance(result, APIFrame)
        assert result.name == "forced"
    finally:
        reset_api_view()


def test_maybe_wrap_frame_api_true_with_custom_wrapper_uses_wrap_frame():
    from pycsamt.api.view.config import PYCSAMT_API_VIEW

    calls = []

    def custom(data, **meta):
        calls.append(meta)
        return data

    PYCSAMT_API_VIEW.configure(wrapper=custom)
    try:
        df = pd.DataFrame({"a": [1]})
        maybe_wrap_frame(df, api=True, name="via-wrapper")
        assert calls and calls[0]["name"] == "via-wrapper"
    finally:
        reset_api_view()


def test_maybe_wrap_frame_api_none_respects_disabled_global():
    from pycsamt.api.view.config import PYCSAMT_API_VIEW

    PYCSAMT_API_VIEW.configure(backend=False)
    try:
        df = pd.DataFrame({"a": [1]})
        assert maybe_wrap_frame(df) is df
    finally:
        reset_api_view()


def test_maybe_wrap_frame_api_none_wraps_when_enabled():
    result = maybe_wrap_frame(pd.DataFrame({"a": [1]}), name="enabled")
    assert isinstance(result, APIFrame)


def test_default_wrap_frame_returns_same_instance_without_overrides():
    base = APIFrame({"a": [1]})
    result = maybe_wrap_frame(base, api=True)
    assert result is base


def test_default_wrap_frame_rewraps_when_overrides_given():
    base = APIFrame({"a": [1]}, name="base", kind="k1")
    result = maybe_wrap_frame(base, api=True, name="override")
    assert result is not base
    assert result.name == "override"
    assert result.kind == "k1"


def test_default_wrap_frame_rewraps_when_copy_true():
    base = APIFrame({"a": [1]}, name="base")
    result = maybe_wrap_frame(base, api=True, copy=True)
    assert result is not base
    assert result.df is not base.df


def test_api_frame_decorator_no_parens():
    @api_frame
    def build():
        return pd.DataFrame({"a": [1]})

    result = build()
    assert isinstance(result, APIFrame)
    assert result.name == "build"


def test_api_frame_decorator_with_args():
    @api_frame(name="custom", kind="k")
    def build():
        return pd.DataFrame({"a": [1]})

    result = build()
    assert result.name == "custom"
    assert result.kind == "k"


def test_api_frame_decorator_passes_through_non_dataframe():
    @api_frame
    def build():
        return "not a frame"

    assert build() == "not a frame"
