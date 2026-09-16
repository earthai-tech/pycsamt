# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for pycsamt.app.converter.pages.settings_page.SettingsPage."""

from __future__ import annotations

import pytest

pytest.importorskip("PySide6", reason="PySide6 required")

from pycsamt.app.converter.pages.settings_page import SettingsPage  # noqa: E402
from pycsamt.app.converter.settings import ConverterSettings, load_settings  # noqa: E402


@pytest.fixture
def page(qapp, isolated_qsettings):
    p = SettingsPage()
    yield p
    p.close()


def test_page_loads_defaults_on_construction(page):
    assert page._epsg.value() == 0
    assert page._utm_zone.text() == ""
    assert page._station_z.currentText() == "auto"
    assert page._on_loss.currentText() == "warn"
    assert page._prefer_spectra.isChecked() is True
    assert page._xml_strict.isChecked() is True
    assert page._log10_view.isChecked() is False
    assert page._open_after.isChecked() is False
    assert page._overwrite.isChecked() is False


def test_page_loads_persisted_values(qapp, isolated_qsettings):
    from pycsamt.app.converter.settings import save_settings

    save_settings(
        ConverterSettings(
            default_epsg=32650,
            default_utm_zone="48N",
            station_z_convention="elevation",
            air_threshold_ohm_m=5e7,
            on_loss="raise",
            prefer_spectra=False,
            xml_strict=False,
            log10_view=True,
            open_output_folder_after_convert=True,
            overwrite_without_asking=True,
        )
    )
    p = SettingsPage()
    assert p._epsg.value() == 32650
    assert p._utm_zone.text() == "48N"
    assert p._station_z.currentText() == "elevation"
    assert p._air_threshold.value() == 5e7
    assert p._on_loss.currentText() == "raise"
    assert p._prefer_spectra.isChecked() is False
    assert p._xml_strict.isChecked() is False
    assert p._log10_view.isChecked() is True
    assert p._open_after.isChecked() is True
    assert p._overwrite.isChecked() is True
    p.close()


def test_save_writes_values_and_updates_status(page):
    page._epsg.setValue(4326)
    page._utm_zone.setText("31N")
    page._station_z.setCurrentText("depth_down")
    page._air_threshold.setValue(1e6)
    page._on_loss.setCurrentText("ignore")
    page._prefer_spectra.setChecked(False)
    page._xml_strict.setChecked(False)
    page._log10_view.setChecked(True)
    page._open_after.setChecked(True)
    page._overwrite.setChecked(True)

    page._on_save()

    assert "Settings saved." in page._status.text()

    reloaded = load_settings()
    assert reloaded.default_epsg == 4326
    assert reloaded.default_utm_zone == "31N"
    assert reloaded.station_z_convention == "depth_down"
    assert reloaded.air_threshold_ohm_m == 1e6
    assert reloaded.on_loss == "ignore"
    assert reloaded.prefer_spectra is False
    assert reloaded.xml_strict is False
    assert reloaded.log10_view is True
    assert reloaded.open_output_folder_after_convert is True
    assert reloaded.overwrite_without_asking is True


def test_save_with_epsg_zero_persists_as_none(page):
    page._epsg.setValue(0)
    page._on_save()
    reloaded = load_settings()
    assert reloaded.default_epsg is None
