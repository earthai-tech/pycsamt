# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
from __future__ import annotations

import re

import pytest

from pycsamt.utils.text import (
    fmt_text,
    listing_items_format,
    smart_format,
    str2columns,
)

# ------------------------------- smart_format ------------------------------


def test_smart_format_non_iterable_quoted_and_unquoted():
    assert smart_format(42) == "42"
    assert smart_format("x") == "'x'"
    assert smart_format("x", quote=False) == "x"


def test_smart_format_empty_list():
    assert smart_format([]) == ""


def test_smart_format_single_item():
    assert smart_format(["a"]) == "'a'"


def test_smart_format_two_items_default_and_custom_choice():
    assert smart_format(["a", "b"]) == "'a' and 'b'"
    assert smart_format(["a", "b"], choice="or") == "'a' or 'b'"


def test_smart_format_more_than_two_oxford_comma_true():
    assert smart_format(["a", "b", "c"]) == "'a', 'b', and 'c'"


def test_smart_format_more_than_two_oxford_comma_false():
    assert smart_format(["a", "b", "c"], oxford_comma=False) == (
        "'a', 'b' and 'c'"
    )


def test_smart_format_quote_false():
    assert smart_format(["a", "b", "c"], quote=False) == "a, b, and c"


# -------------------------------- fmt_text ----------------------------------


def test_fmt_text_dict_scalar_values():
    out = fmt_text({"a": 1, "b": 2})
    assert "a" in out and "1" in out
    lines = out.splitlines()
    assert lines[0] == lines[2] == lines[-1]


def test_fmt_text_dict_iterable_values():
    out = fmt_text({"a": [1, 2], "b": [3, 4]})
    assert "a" in out
    assert "3" in out


def test_fmt_text_dict_mixed_scalar_and_iterable_values():
    out = fmt_text({"a": 1, "b": [2, 3, 4]})
    lines = out.splitlines()
    assert lines[0] == lines[2] == lines[-1]
    assert "b" in out and "4" in out


def test_fmt_text_list_of_rows():
    out = fmt_text([[1, 2], [3, 4]])
    assert "1" in out and "4" in out


def test_fmt_text_ragged_list_of_rows():
    out = fmt_text([[1, 2], [3, 4, 5]])
    lines = out.splitlines()
    assert lines[0] == lines[2] == lines[-1]
    assert "5" in out


def test_fmt_text_single_flat_row():
    out = fmt_text([1, 2, 3])
    lines = out.splitlines()
    header = lines[1]
    body = lines[3]
    assert "1" in body and "2" in body and "3" in body
    assert "0" in header and "1" in header and "2" in header


def test_fmt_text_headers_shorter_than_ncols_padded():
    out = fmt_text([[1, 2, 3]], headers=["only"])
    header_line = out.splitlines()[1]
    assert "only" in header_line


def test_fmt_text_scalar_non_iterable_input():
    out = fmt_text(42)
    assert "42" in out


def test_fmt_text_custom_col_sep_inline_width():
    out = fmt_text([[1, 2]], headers=["a", "b"], col_sep=" ~ ", inline="=")
    lines = out.splitlines()
    assert lines[0][0] == "="
    assert "~" in lines[1]


# ------------------------------- str2columns --------------------------------


def test_str2columns_default_pattern():
    assert str2columns("hello, world! foo-bar") == [
        "hello",
        "world",
        "foo",
        "bar",
    ]


def test_str2columns_custom_regex_pattern_object():
    splitter = re.compile(r"\s+")
    assert str2columns("a b  c", regex=splitter) == ["a", "b", "c"]


def test_str2columns_custom_pattern_string():
    assert str2columns("a-b_c", pattern=r"[-_]") == ["a", "b", "c"]


def test_str2columns_lower_true():
    assert str2columns("Hello_World", lower=True) == ["hello", "world"]


def test_str2columns_unique_true_preserves_order():
    assert str2columns("a,a,b,b,c", unique=True) == ["a", "b", "c"]


def test_str2columns_strip_chars():
    assert str2columns("  .a.  .b.  ", strip_chars=".") == ["a", "b"]


# --------------------------- listing_items_format ---------------------------


def test_listing_items_format_verbose_true_prints_and_returns_none(capsys):
    result = listing_items_format(
        ["a", "b"], begin_text="Items:", end_text="Done", verbose=True
    )
    assert result is None
    captured = capsys.readouterr()
    assert "Items:" in captured.out
    assert "a" in captured.out
    assert "Done" in captured.out


def test_listing_items_format_verbose_false_returns_string():
    result = listing_items_format(["a", "b"], verbose=False)
    assert isinstance(result, str)
    assert "a" in result and "b" in result


def test_listing_items_format_begin_end_text():
    result = listing_items_format(
        ["a"], begin_text="Start", end_text="End", verbose=False
    )
    assert result.startswith("Start")
    assert result.endswith("End")


def test_listing_items_format_enumerate_items_true_vs_false():
    enumerated = listing_items_format(["a", "b"], verbose=False)
    assert "1." in enumerated and "2." in enumerated

    bulleted = listing_items_format(
        ["a", "b"], verbose=False, enumerate_items=False, bullet="*"
    )
    assert "*" in bulleted
    assert "1." not in bulleted


def test_listing_items_format_non_iterable_single_item_wrapped():
    result = listing_items_format(5, verbose=False)
    assert "5" in result
    assert result.count(".") == 1


def test_listing_items_format_broken_iterable_falls_back_to_wrap():
    class Broken:
        def __iter__(self):
            raise RuntimeError("boom")

        def __str__(self):
            return "BROKEN"

    result = listing_items_format(Broken(), verbose=False)
    assert "BROKEN" in result
