# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for pycsamt.app.converter.settings."""

from __future__ import annotations

import pytest

pytest.importorskip("PySide6", reason="PySide6 required")

from pycsamt.app.converter.settings import (  # noqa: E402
    ConverterSettings,
    _coerce,
    load_settings,
    save_settings,
)


def test_converter_settings_defaults():
    s = ConverterSettings()
    assert s.default_epsg is None
    assert s.default_utm_zone == ""
    assert s.station_z_convention == "auto"
    assert s.air_threshold_ohm_m == 1e8
    assert s.on_loss == "warn"
    assert s.prefer_spectra is True
    assert s.xml_strict is True
    assert s.default_pcsf_format == "pcsf"
    assert s.log10_view is False
    assert s.open_output_folder_after_convert is False
    assert s.overwrite_without_asking is False
    assert s.theme == "system"
    assert s.last_dirs == {}


def test_converter_settings_last_dirs_not_shared_between_instances():
    # __post_init__ must give each instance its own dict, not a shared
    # mutable default.
    a = ConverterSettings()
    b = ConverterSettings()
    a.last_dirs["inversion"] = "/some/path"
    assert b.last_dirs == {}


def test_load_settings_returns_defaults_when_unset(qapp, isolated_qsettings):
    s = load_settings()
    assert s == ConverterSettings()


def test_save_and_load_settings_roundtrip(qapp, isolated_qsettings):
    original = ConverterSettings(
        default_epsg=32650,
        default_utm_zone="48N",
        station_z_convention="elevation",
        air_threshold_ohm_m=5e7,
        on_loss="raise",
        prefer_spectra=False,
        xml_strict=False,
        default_pcsf_format="pcsm",
        log10_view=True,
        open_output_folder_after_convert=True,
        overwrite_without_asking=True,
        theme="dark",
        last_dirs={"inversion": "/a/b", "batch": "/c/d"},
    )
    save_settings(original)
    loaded = load_settings()
    assert loaded == original


def test_save_settings_persists_across_qsettings_instances(qapp, isolated_qsettings):
    save_settings(ConverterSettings(default_epsg=4326))
    # A second, independent load_settings() call must see the same file.
    loaded_again = load_settings()
    assert loaded_again.default_epsg == 4326


def test_load_settings_partial_write_falls_back_to_defaults(qapp, isolated_qsettings):
    from PySide6.QtCore import QSettings

    from pycsamt.app.converter import settings as settings_mod

    qs = settings_mod._qsettings()
    qs.setValue("default_epsg", 111)
    qs.sync()
    loaded = load_settings()
    assert loaded.default_epsg == 111
    # Everything else must fall back to the dataclass default.
    assert loaded.station_z_convention == "auto"
    assert loaded.theme == "system"


# ── _coerce ──────────────────────────────────────────────────────────────


def test_coerce_none_returns_default():
    assert _coerce(None, 42) == 42
    assert _coerce(None, "x") == "x"


@pytest.mark.parametrize(
    "raw,expected",
    [
        (True, True),
        (False, False),
        ("true", True),
        ("True", True),
        ("1", True),
        ("yes", True),
        ("on", True),
        ("false", False),
        ("0", False),
        ("no", False),
        ("garbage", False),
    ],
)
def test_coerce_bool(raw, expected):
    assert _coerce(raw, False) is expected


def test_coerce_int_valid_and_invalid():
    assert _coerce("42", 0) == 42
    assert _coerce("not-an-int", 7) == 7
    assert _coerce(None, 7) == 7


def test_coerce_float_valid_and_invalid():
    assert _coerce("1.5", 0.0) == 1.5
    assert _coerce("nope", 2.5) == 2.5


def test_coerce_passthrough_for_other_types():
    assert _coerce("48N", "") == "48N"
    assert _coerce({"a": 1}, {}) == {"a": 1}
