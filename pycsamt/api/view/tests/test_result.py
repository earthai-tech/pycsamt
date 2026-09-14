from __future__ import annotations

import pandas as pd
import pytest

from pycsamt.api.view.config import reset_api_view
from pycsamt.api.view.frame import APIFrame
from pycsamt.api.view.result import APIResult, wrap_result


@pytest.fixture(autouse=True)
def _reset_config():
    reset_api_view()
    yield
    reset_api_view()


def test_default_name_and_kind():
    res = APIResult()
    assert res.name == "result"
    assert res.kind is None
    assert res.meta == {}
    assert res.keys == ()


def test_items_stored_and_accessible_as_attributes():
    res = APIResult(name="r", kind="k", a=1, b="two")
    assert res.name == "r"
    assert res.kind == "k"
    assert res.keys == ("a", "b")
    assert res.a == 1
    assert res.b == "two"


def test_getitem_setitem_contains():
    res = APIResult(a=1)
    assert res["a"] == 1
    res["b"] = 2
    assert res["b"] == 2
    assert "b" in res
    assert "z" not in res


def test_iter_and_len():
    res = APIResult(a=1, b=2)
    assert set(iter(res)) == {"a", "b"}
    assert len(res) == 2


def test_dir_includes_items():
    res = APIResult(a=1)
    assert "a" in dir(res)


def test_items_method():
    res = APIResult(a=1, b=2)
    assert dict(res.items()) == {"a": 1, "b": 2}


def test_tables_filters_api_frame_values():
    frame = APIFrame({"x": [1]})
    res = APIResult(a=frame, b="not a table")
    tables = res.tables()
    assert tables == {"a": frame}


def test_update_meta_returns_self():
    res = APIResult()
    result = res.update_meta(x=1)
    assert result is res
    assert res.meta == {"x": 1}


def test_str_without_kind_or_tables():
    res = APIResult(name="plain", a=1)
    text = str(res)
    assert "APIResult: plain" in text
    assert "kind" not in text
    assert "items: a" in text
    assert "tables" not in text


def test_str_with_kind_and_no_items():
    res = APIResult(name="n", kind="k")
    text = str(res)
    assert "kind: k" in text
    assert "items: -" in text


def test_str_with_tables_shows_shapes():
    frame = APIFrame({"x": [1, 2], "y": [3, 4]})
    res = APIResult(name="n", summary=frame)
    text = str(res)
    assert "tables: summary=2x2" in text


def test_wrap_result_wraps_dataframes_by_default():
    res = wrap_result(
        {"summary": pd.DataFrame({"a": [1]}), "note": "ok"},
        name="run",
        kind="qc",
    )
    assert isinstance(res, APIResult)
    assert isinstance(res.summary, APIFrame)
    assert res.summary.kind == "qc"
    assert res.note == "ok"


def test_wrap_result_wrap_tables_false_keeps_raw_dataframe():
    res = wrap_result(
        {"summary": pd.DataFrame({"a": [1]})}, wrap_tables=False,
    )
    assert isinstance(res.summary, pd.DataFrame)
    assert not isinstance(res.summary, APIFrame)


def test_wrap_result_forwards_meta():
    res = wrap_result({"a": 1}, meta={"owner": "me"})
    assert res.meta == {"owner": "me"}
