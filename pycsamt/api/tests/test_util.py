from __future__ import annotations

import sys
import types
import warnings

import numpy as np
import pandas as pd
import pytest

from pycsamt.api import util as u


# ─────────────────────────────────────────────────────────────────────────
# TerminalSize
# ─────────────────────────────────────────────────────────────────────────


def test_get_terminal_size_uses_shutil_success():
    cols, rows = u.TerminalSize.get_terminal_size()
    assert isinstance(cols, int) and isinstance(rows, int)


def test_get_terminal_size_dispatches_to_windows_helper(monkeypatch):
    def _raise(*a, **k):
        raise OSError("no tty")

    monkeypatch.setattr(u.shutil, "get_terminal_size", _raise)
    monkeypatch.setattr(u.os, "name", "nt")
    monkeypatch.setattr(
        u.TerminalSize, "_get_terminal_size_windows", staticmethod(lambda: (1, 2))
    )
    assert u.TerminalSize.get_terminal_size() == (1, 2)


def test_get_terminal_size_dispatches_to_unix_helper(monkeypatch):
    def _raise(*a, **k):
        raise OSError("no tty")

    monkeypatch.setattr(u.shutil, "get_terminal_size", _raise)
    monkeypatch.setattr(u.os, "name", "posix")
    monkeypatch.setattr(
        u.TerminalSize, "_get_terminal_size_unix", staticmethod(lambda: (3, 4))
    )
    assert u.TerminalSize.get_terminal_size() == (3, 4)


def test_get_terminal_size_falls_back_to_default_for_unknown_os(monkeypatch):
    def _raise(*a, **k):
        raise OSError("no tty")

    monkeypatch.setattr(u.shutil, "get_terminal_size", _raise)
    monkeypatch.setattr(u.os, "name", "java")
    assert u.TerminalSize.get_terminal_size() == u.TerminalSize.DEFAULT_SIZE


@pytest.mark.skipif(
    sys.platform != "win32", reason="ctypes.windll only exists on Windows"
)
def test_get_terminal_size_windows_falls_back_without_real_console():
    # Under pytest there is normally no attached console, so this
    # exercises the "res falsy / exception" fallback path for real.
    size = u.TerminalSize._get_terminal_size_windows()
    assert size == u.TerminalSize.DEFAULT_SIZE


@pytest.mark.skipif(
    sys.platform != "win32", reason="ctypes.windll only exists on Windows"
)
def test_get_terminal_size_windows_success(monkeypatch):
    import struct
    from ctypes import windll

    packed = struct.pack("hhhhHhhhhhh", 0, 0, 0, 0, 0, 0, 0, 79, 24, 0, 0)

    def fake_get_std_handle(_which):
        return 1

    def fake_get_console_screen_buffer_info(_handle, csbi):
        csbi.raw = packed
        return 1

    monkeypatch.setattr(
        windll.kernel32, "GetStdHandle", fake_get_std_handle,
    )
    monkeypatch.setattr(
        windll.kernel32,
        "GetConsoleScreenBufferInfo",
        fake_get_console_screen_buffer_info,
    )
    size = u.TerminalSize._get_terminal_size_windows()
    assert size == (80, 25)


def _install_fake_posix_modules(monkeypatch, *, ioctl_result):
    fcntl_mod = types.ModuleType("fcntl")
    termios_mod = types.ModuleType("termios")
    termios_mod.TIOCGWINSZ = 0x5413

    def fake_ioctl(fd, request, buf):
        if ioctl_result is None:
            raise OSError("no tty")
        return ioctl_result

    fcntl_mod.ioctl = fake_ioctl
    monkeypatch.setitem(sys.modules, "fcntl", fcntl_mod)
    monkeypatch.setitem(sys.modules, "termios", termios_mod)


def test_get_terminal_size_unix_success_via_fd(monkeypatch):
    import struct

    packed = struct.pack("hh", 40, 120)  # (rows, cols)
    _install_fake_posix_modules(monkeypatch, ioctl_result=packed)
    size = u.TerminalSize._get_terminal_size_unix()
    assert size == (120, 40)


def test_get_terminal_size_unix_success_via_ctermid(monkeypatch, tmp_path):
    import struct

    packed = struct.pack("hh", 40, 120)  # (rows, cols)
    _install_fake_posix_modules(monkeypatch, ioctl_result=None)
    ctermid_file = tmp_path / "fake_tty"
    ctermid_file.write_bytes(b"")
    monkeypatch.setattr(
        u.os, "ctermid", lambda: str(ctermid_file), raising=False,
    )

    def fake_ioctl(fd, request, buf):
        if fd in (0, 1, 2):
            raise OSError("no tty")
        return packed

    sys.modules["fcntl"].ioctl = fake_ioctl
    size = u.TerminalSize._get_terminal_size_unix()
    assert size == (120, 40)


def test_get_terminal_size_unix_falls_back_to_default(monkeypatch):
    _install_fake_posix_modules(monkeypatch, ioctl_result=None)
    monkeypatch.setattr(
        u.os,
        "ctermid",
        lambda: (_ for _ in ()).throw(AttributeError("no ctermid")),
        raising=False,
    )
    size = u.TerminalSize._get_terminal_size_unix()
    assert size == u.TerminalSize.DEFAULT_SIZE


# ─────────────────────────────────────────────────────────────────────────
# format_value / apply_precision / validate_precision
# ─────────────────────────────────────────────────────────────────────────


def test_format_value_integer_and_float():
    assert u.format_value(123) == "123"
    assert u.format_value(123.456789) == "123.4568"


def test_format_value_non_numeric_passthrough():
    assert u.format_value("hello") == "hello"


def test_apply_precision_integer_returned_as_int():
    assert u.apply_precision(np.int32(456), 2) == 456


def test_apply_precision_float_rounded_only_when_needed():
    assert u.apply_precision(123.4, 2) == 123.4
    assert u.apply_precision(123.456789, 2) == 123.46


def test_apply_precision_non_numeric_passthrough():
    assert u.apply_precision("abc") == "abc"


def test_validate_precision_accepts_valid_values():
    assert u.validate_precision(3) == 3
    assert u.validate_precision(3.0) == 3
    assert u.validate_precision(None) == 4  # falsy -> default 4


def test_validate_precision_rejects_negative():
    with pytest.raises(ValueError):
        u.validate_precision(-1)


def test_validate_precision_rejects_non_numeric():
    with pytest.raises(ValueError):
        u.validate_precision("three")


# ─────────────────────────────────────────────────────────────────────────
# parse_component_kind
# ─────────────────────────────────────────────────────────────────────────


_PC_LIST = [
    ("pc1", ["f1", "f2"], [0.8, 0.5]),
    ("pc2", ["f1", "f2"], [0.6, 0.4]),
]


def test_parse_component_kind_extracts_valid_component():
    names, values = u.parse_component_kind(_PC_LIST, "pc1")
    assert names == ["f1", "f2"]
    assert values == [0.8, 0.5]


def test_parse_component_kind_out_of_range_raises():
    with pytest.raises(ValueError):
        u.parse_component_kind(_PC_LIST, "pc9")


def test_parse_component_kind_no_digit_raises():
    with pytest.raises(ValueError):
        u.parse_component_kind(_PC_LIST, "pc")


# ─────────────────────────────────────────────────────────────────────────
# find_maximum_table_width
# ─────────────────────────────────────────────────────────────────────────


def test_find_maximum_table_width_finds_longest_header():
    summary = "Title\n====\nrow\n======\n"
    assert u.find_maximum_table_width(summary) == 6


def test_find_maximum_table_width_no_headers_returns_zero():
    assert u.find_maximum_table_width("no headers here") == 0


# ─────────────────────────────────────────────────────────────────────────
# format_text
# ─────────────────────────────────────────────────────────────────────────


def test_format_text_wraps_with_key():
    text = "a " * 40
    out = u.format_text(text, key="Note", key_length=10, max_char_text=30)
    lines = out.split("\n")
    assert lines[0].startswith("Note")
    assert len(lines) > 1


def test_format_text_no_key_but_key_length_given():
    out = u.format_text("short text", key=None, key_length=5)
    assert out.startswith("     ")


def test_format_text_no_key_no_key_length():
    out = u.format_text("short text", key=None, key_length=None)
    assert out == "short text"


def test_format_text_key_length_none_uses_key_length():
    out = u.format_text("value", key="K", key_length=None)
    assert out.startswith("K : ")


def test_format_text_force_break_when_no_space_found():
    text = "x" * 100
    out = u.format_text(text, key="K", key_length=3, max_char_text=20)
    assert "\n" in out


def test_format_text_add_frame_lines():
    out = u.format_text(
        "hello", key=None, key_length=None,
        add_frame_lines=True, max_char_text=10,
    )
    lines = out.split("\n")
    assert lines[0] == "=" * 10
    assert lines[-1] == "=" * 10
    assert "hello" in out


def test_format_text_narrow_width_with_default_key_length_terminates():
    # Regression test: max_char_text smaller than the default key_length
    # padding used to make the internal wrap width negative, which left
    # `text` unchanged after slicing and spun forever. It must now always
    # terminate, and every character of the input must still appear.
    out = u.format_text("hello", add_frame_lines=True, max_char_text=10)
    assert "".join(ch for ch in out if ch.isalpha()) == "hello"


# ─────────────────────────────────────────────────────────────────────────
# get_frame_chars
# ─────────────────────────────────────────────────────────────────────────


@pytest.mark.parametrize(
    "char,expected",
    [
        ("[", ("]", "[", "]")),
        ("{", ("}", "{", "}")),
        ("(", (")", "(", ")")),
        ("<", (">", "<", ">")),
    ],
)
def test_get_frame_chars_known(char, expected):
    assert u.get_frame_chars(char) == expected


def test_get_frame_chars_unknown_falls_back_to_dots():
    assert u.get_frame_chars("?") == (".", ".", ".")


# ─────────────────────────────────────────────────────────────────────────
# format_cell
# ─────────────────────────────────────────────────────────────────────────


def test_format_cell_truncates_long_text():
    out = u.format_cell("a very long string value", 10)
    assert out.endswith("...")
    assert len(out) == 10


def test_format_cell_right_aligns_with_max_width():
    out = u.format_cell("ab", 10, max_width=6)
    assert out == "    ab"


def test_format_cell_no_max_width_returns_as_is():
    assert u.format_cell("ab", 10) == "ab"


# ─────────────────────────────────────────────────────────────────────────
# get_table_width_from
# ─────────────────────────────────────────────────────────────────────────


def test_get_table_width_from_no_border_warns():
    with pytest.warns(UserWarning):
        result = u.get_table_width_from("no borders here", error="warn")
    assert result is None


def test_get_table_width_from_no_border_raises():
    with pytest.raises(ValueError):
        u.get_table_width_from("no borders here", error="raise")


def test_get_table_width_from_no_border_ignore():
    assert u.get_table_width_from("no borders here", error="ignore") is None


def test_get_table_width_from_shallow_check_first():
    text = "=====\nrow\n===\n"
    assert (
        u.get_table_width_from(text, deep_check=False, check_first=True) == 5
    )


def test_get_table_width_from_shallow_check_last():
    text = "=====\nrow\n===\n"
    assert (
        u.get_table_width_from(text, deep_check=False, check_first=False)
        == 3
    )


def test_get_table_width_from_deep_check_strategies():
    text = "=====\nrow\n===\n"
    assert u.get_table_width_from(text, width_strategy="max") == 5
    assert u.get_table_width_from(text, width_strategy="min") == 3
    assert u.get_table_width_from(text, width_strategy="average") == 4


# ─────────────────────────────────────────────────────────────────────────
# generate_legend
# ─────────────────────────────────────────────────────────────────────────


def test_generate_legend_rejects_non_dict_custom_markers():
    with pytest.raises(TypeError):
        u.generate_legend(custom_markers=["not", "a", "dict"])


def test_generate_legend_default_hides_diagonal():
    text = u.generate_legend()
    assert "Diagonal" not in text


def test_generate_legend_shows_diagonal_when_requested():
    text = u.generate_legend(hide_diag=False)
    assert "Diagonal" in text


def test_generate_legend_custom_markers_merge_with_defaults():
    text = u.generate_legend(custom_markers={"++": "Custom positive"})
    normalized = " ".join(text.split())
    assert "Custom positive" in normalized


def test_generate_legend_removes_placeholder_when_falsy():
    text = u.generate_legend(no_corr_placeholder="")
    assert "Non-correlated" not in text


# ─────────────────────────────────────────────────────────────────────────
# to_snake_case
# ─────────────────────────────────────────────────────────────────────────


def test_to_snake_case_standard_camel():
    assert u.to_snake_case("CamelCaseName") == "camel_case_name"


def test_to_snake_case_standard_with_symbols():
    assert u.to_snake_case("Hello-World!!") == "hello_world"


def test_to_snake_case_soft_mode():
    assert u.to_snake_case("  Hello   World  ", mode="soft") == "hello_world"


# ─────────────────────────────────────────────────────────────────────────
# generate_column_name_mapping
# ─────────────────────────────────────────────────────────────────────────


def test_generate_column_name_mapping():
    mapping = u.generate_column_name_mapping(["ColumnOne", "ColumnTwo"])
    assert mapping == {
        "column_one": "ColumnOne",
        "column_two": "ColumnTwo",
    }


# ─────────────────────────────────────────────────────────────────────────
# series_to_dataframe
# ─────────────────────────────────────────────────────────────────────────


def test_series_to_dataframe_string_index():
    series = pd.Series(data=[1, 2, 3], index=["a", "b", "c"])
    df = u.series_to_dataframe(series)
    assert list(df.columns) == ["a", "b", "c"]
    assert df.iloc[0].tolist() == [1, 2, 3]


def test_series_to_dataframe_numeric_index_converted_to_str():
    series = pd.Series(data=[4, 5, 6], index=[10, 20, 30])
    df = u.series_to_dataframe(series)
    assert list(df.columns) == ["10", "20", "30"]


def test_series_to_dataframe_rejects_non_series():
    with pytest.raises(TypeError):
        u.series_to_dataframe([1, 2, 3])


# ─────────────────────────────────────────────────────────────────────────
# get_table_size / get_terminal_size (module function)
# ─────────────────────────────────────────────────────────────────────────


def test_get_table_size_auto_uses_terminal_width(monkeypatch):
    monkeypatch.setattr(u, "get_terminal_size", lambda: (100, 30))
    assert u.get_table_size(width="auto") == 100


def test_get_table_size_explicit_within_terminal_width(monkeypatch):
    monkeypatch.setattr(u, "get_terminal_size", lambda: (100, 30))
    assert u.get_table_size(width=50) == 50


def test_get_table_size_explicit_exceeds_terminal_width_warns(monkeypatch):
    monkeypatch.setattr(u, "get_terminal_size", lambda: (100, 30))
    with pytest.warns(UserWarning):
        result = u.get_table_size(width=200, error="warn")
    assert result == 200


def test_get_table_size_explicit_exceeds_width_no_warning_when_ignored(
    monkeypatch,
):
    monkeypatch.setattr(u, "get_terminal_size", lambda: (100, 30))
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        result = u.get_table_size(width=200, error="ignore")
    assert result == 200


def test_get_table_size_invalid_width_raises(monkeypatch):
    monkeypatch.setattr(u, "get_terminal_size", lambda: (100, 30))
    with pytest.raises(ValueError):
        u.get_table_size(width="not-a-number")


def test_get_table_size_return_height(monkeypatch):
    monkeypatch.setattr(u, "get_terminal_size", lambda: (100, 30))
    assert u.get_table_size(width="auto", return_height=True) == (100, 30)


def test_module_get_terminal_size_returns_two_ints():
    cols, rows = u.get_terminal_size()
    assert isinstance(cols, int) and isinstance(rows, int)


# ─────────────────────────────────────────────────────────────────────────
# to_camel_case
# ─────────────────────────────────────────────────────────────────────────


def test_to_camel_case_with_explicit_delimiter():
    assert u.to_camel_case("outlier_results", "_") == "OutlierResults"


def test_to_camel_case_auto_detects_space_delimiter():
    assert u.to_camel_case("outlier results") == "OutlierResults"


def test_to_camel_case_auto_detects_underscore_delimiter():
    assert u.to_camel_case("data_science_rocks") == "DataScienceRocks"


def test_to_camel_case_both_space_and_underscore_present():
    assert u.to_camel_case("data_science rocks") == "DataScienceRocks"


def test_to_camel_case_no_delimiter_single_word():
    assert u.to_camel_case("outlierresults") == "Outlierresults"


def test_to_camel_case_use_regex():
    assert u.to_camel_case("multi@var_analysis", use_regex=True) == (
        "MultiVarAnalysis"
    )


def test_to_camel_case_already_camelcase_returned_as_is():
    assert u.to_camel_case("OutlierResults") == "OutlierResults"
    assert u.to_camel_case("BoxFormatter") == "BoxFormatter"


# ─────────────────────────────────────────────────────────────────────────
# beautify_dict
# ─────────────────────────────────────────────────────────────────────────


def test_beautify_dict_rejects_non_dict():
    with pytest.raises(TypeError):
        u.beautify_dict(["not", "a", "dict"])


def test_beautify_dict_empty_dict():
    out = u.beautify_dict({})
    assert out.startswith("{")
    assert out.endswith("}")


def test_beautify_dict_basic_formatting():
    d = {"b": "second", "a": "first"}
    out = u.beautify_dict(d, space=2, max_char=100)
    assert "a" in out and "b" in out
    assert out.index("'a'".strip("'")) < out.index("second")


def test_beautify_dict_truncates_long_values():
    d = {"a": "x" * 100}
    out = u.beautify_dict(d, max_char=10)
    assert "..." in out


def test_beautify_dict_with_key_indents_and_truncates_lines():
    d = {"a": "value_one", "bb": "value_two_longer"}
    out = u.beautify_dict(d, key="MyDict", max_char=15)
    assert out.startswith("MyDict : {")


def test_beautify_dict_auto_max_char(monkeypatch):
    monkeypatch.setattr(u, "get_terminal_size", lambda: (100, 30))
    out = u.beautify_dict({"a": "short"})
    assert "short" in out


# ─────────────────────────────────────────────────────────────────────────
# remove_extra_spaces
# ─────────────────────────────────────────────────────────────────────────


def test_remove_extra_spaces():
    text = "this is      text that    have   extra          space."
    assert u.remove_extra_spaces(text) == (
        "this is text that have extra space."
    )


# ─────────────────────────────────────────────────────────────────────────
# format_iterable
# ─────────────────────────────────────────────────────────────────────────


def test_format_iterable_numeric_list():
    out = u.format_iterable([1, 2, 3])
    assert out.startswith("list (min=1")
    assert "len=3" in out


def test_format_iterable_numeric_ndarray():
    out = u.format_iterable(np.array([1.0, 2.0, 3.0]))
    assert out.startswith("ndarray (")
    assert "shape=(3,)" in out


def test_format_iterable_non_numeric_ndarray():
    out = u.format_iterable(np.array(["a", "b"]))
    assert out.startswith("ndarray (")
    assert "min" not in out


def test_format_iterable_numeric_series():
    out = u.format_iterable(pd.Series([1.0, 2.0, 3.0]))
    assert out.startswith("Series (")
    assert "len=3" in out


def test_format_iterable_object_series():
    out = u.format_iterable(pd.Series(["a", "b"]))
    assert out.startswith("Series (len=2")


def test_format_iterable_numeric_dataframe():
    df = pd.DataFrame({"a": [1.0, 2.0], "b": [3.0, 4.0]})
    out = u.format_iterable(df)
    assert out.startswith("DataFrame (")
    assert "n_rows=2" in out


def test_format_iterable_dataframe_without_numeric_columns():
    df = pd.DataFrame({"a": ["x", "y"]})
    out = u.format_iterable(df)
    assert out.startswith("DataFrame (n_rows=2")


def test_format_iterable_fallback_to_str():
    assert u.format_iterable({"a": 1}) == str({"a": 1})


# ─────────────────────────────────────────────────────────────────────────
# format_dict_result
# ─────────────────────────────────────────────────────────────────────────


def test_format_dict_result_basic_and_message():
    d = {"key1": "short value", "key2": 42}
    out = u.format_dict_result(d, dict_name="D", include_message=True)
    assert out.startswith("D({")
    assert "Use <D.key>" in out


def test_format_dict_result_truncates_long_values():
    d = {"key1": "x" * 100}
    out = u.format_dict_result(d, max_char=10)
    assert "..." in out


def test_format_dict_result_no_message_by_default():
    out = u.format_dict_result({"a": 1})
    assert "Use <" not in out


def test_format_dict_result_metaclass_name_access_failure_falls_back():
    class Meta(type):
        @property
        def __name__(cls):  # noqa: N807
            raise RuntimeError("no class name")

    class Weird(metaclass=Meta):
        pass

    instance = Weird()
    instance.__name__ = "fallback_name"
    out = u.format_dict_result({"weird": instance})
    assert "fallback_name" in out


# ─────────────────────────────────────────────────────────────────────────
# count_functions
# ─────────────────────────────────────────────────────────────────────────


@pytest.fixture
def fixture_module(tmp_path, monkeypatch):
    src = '''
def public_func():
    def _nested_local():
        pass
    return _nested_local


def _private_func():
    pass


class PublicClass:
    pass


class _PrivateClass:
    pass
'''
    mod_path = tmp_path / "pycsamt_util_fixture_mod.py"
    mod_path.write_text(src, encoding="utf-8")
    monkeypatch.syspath_prepend(str(tmp_path))
    import sys as _sys

    _sys.modules.pop("pycsamt_util_fixture_mod", None)
    yield "pycsamt_util_fixture_mod"
    _sys.modules.pop("pycsamt_util_fixture_mod", None)


def test_count_functions_default_public_top_level_only(fixture_module):
    result = u.count_functions(fixture_module, return_counts=False)
    assert result == ["public_func"]


def test_count_functions_include_private(fixture_module):
    result = u.count_functions(
        fixture_module, return_counts=False, include_private=True,
    )
    assert result == ["_private_func", "public_func"]


def test_count_functions_include_local(fixture_module):
    result = u.count_functions(
        fixture_module,
        return_counts=False,
        include_private=True,
        include_local=True,
    )
    assert "_nested_local" in result


def test_count_functions_include_class_counts(fixture_module):
    count = u.count_functions(
        fixture_module, include_class=True, return_counts=True,
    )
    assert count == 2  # public_func + PublicClass


def test_count_functions_include_class_and_private_listing(fixture_module):
    result = u.count_functions(
        fixture_module,
        include_class=True,
        return_counts=False,
        include_private=True,
    )
    assert set(result) >= {
        "PublicClass", "_PrivateClass", "public_func", "_private_func",
    }


# ─────────────────────────────────────────────────────────────────────────
# round_numeric_values
# ─────────────────────────────────────────────────────────────────────────


def test_round_numeric_values_rounds_floats_only():
    df = pd.DataFrame(
        {
            "A": [1.12345, 2.6789, 3],
            "B": [4, 5.98765, "text"],
        }
    )
    out = u.round_numeric_values(df, precision=2)
    assert out["A"][0] == 1.12
    assert out["B"][1] == 5.99
    assert out["B"][2] == "text"
    assert out["A"][2] == 3  # int untouched
