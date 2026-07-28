# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for AboutDialog, _HeroBanner, and _LinkButton."""

from __future__ import annotations

import sys
from unittest import mock

import pytest

pytest.importorskip("PySide6", reason="PySide6 required")

from PySide6.QtWidgets import QDialog, QPushButton

from pycsamt.app.desktop.dialogs import about_dialog as about_mod
from pycsamt.app.desktop.dialogs.about_dialog import (
    AboutDialog,
    _HeroBanner,
    _LinkButton,
)

# ── _HeroBanner ───────────────────────────────────────────────────────────


def test_hero_banner_creates_with_real_logo(qapp):
    banner = _HeroBanner()
    assert banner._renderer is not None
    banner.close()


def test_hero_banner_paints_without_error(qapp):
    banner = _HeroBanner()
    banner.resize(400, 175)
    pix = banner.grab()  # forces a real paintEvent pass
    assert not pix.isNull()
    banner.close()


def test_hero_banner_missing_logo_file(qapp, monkeypatch, tmp_path):
    monkeypatch.setattr(about_mod, "_LOGO_SVG", tmp_path / "nope.svg")
    banner = _HeroBanner()
    assert banner._renderer is None
    banner.resize(400, 175)
    pix = banner.grab()
    assert not pix.isNull()
    banner.close()


def test_hero_banner_invalid_svg_content(qapp, monkeypatch, tmp_path):
    bad_svg = tmp_path / "bad.svg"
    bad_svg.write_text("not an svg at all")
    monkeypatch.setattr(about_mod, "_LOGO_SVG", bad_svg)
    banner = _HeroBanner()
    assert banner._renderer is None  # r.isValid() is False -> stays None
    banner.close()


def test_hero_banner_renderer_construction_raises(qapp, monkeypatch):
    monkeypatch.setattr(
        "PySide6.QtSvg.QSvgRenderer",
        mock.Mock(side_effect=RuntimeError("boom")),
    )
    banner = _HeroBanner()  # must not raise
    assert banner._renderer is None
    banner.close()


def test_hero_banner_paint_version_import_failure(qapp, monkeypatch):
    banner = _HeroBanner()
    banner.resize(400, 175)
    with mock.patch.dict(sys.modules, {"pycsamt": None}):
        pix = banner.grab()  # exercises the `except Exception: ver = "2.0.0"` path
        assert not pix.isNull()
    banner.close()


def test_hero_banner_zero_width_no_division_error(qapp):
    banner = _HeroBanner()
    banner.resize(0, 175)
    pix = banner.grab()
    assert pix is not None
    banner.close()


# ── _LinkButton ───────────────────────────────────────────────────────────


def test_link_button_creates(qapp):
    btn = _LinkButton("X", "Label", "https://example.com/")
    assert "Label" in btn.text()
    btn.close()


def test_link_button_click_opens_url(qapp, monkeypatch):
    from PySide6.QtGui import QDesktopServices

    calls = []
    monkeypatch.setattr(
        QDesktopServices,
        "openUrl",
        staticmethod(lambda url: calls.append(url.toString())),
    )
    btn = _LinkButton("X", "Docs", "https://pycsamt.org/")
    btn.click()
    assert calls == ["https://pycsamt.org/"]
    btn.close()


# ── AboutDialog ───────────────────────────────────────────────────────────


@pytest.fixture
def dlg(qapp):
    d = AboutDialog()
    yield d
    d.close()


def test_creates(dlg):
    assert dlg.windowTitle() == "About pycsamt"


def test_fixed_width(dlg):
    assert dlg.maximumWidth() == 530


def test_has_close_button_that_accepts(dlg):
    buttons = dlg.findChildren(QPushButton)
    close_btn = next(b for b in buttons if b.text() == "Close")
    close_btn.click()
    assert dlg.result() == QDialog.DialogCode.Accepted


def test_has_two_link_buttons(dlg):
    links = dlg.findChildren(_LinkButton)
    assert len(links) == 2


def test_has_hero_banner(dlg):
    banners = dlg.findChildren(_HeroBanner)
    assert len(banners) == 1


def test_metadata_version_import_failure(qapp):
    with mock.patch.dict(sys.modules, {"pycsamt": None}):
        d = AboutDialog()  # exercises the `except Exception: ver = "2.0.0"` path
        assert d is not None
        d.close()


def test_parent_wiring(qapp):
    from PySide6.QtWidgets import QWidget

    parent = QWidget()
    d = AboutDialog(parent=parent)
    assert d.parent() is parent
    d.close()
    parent.close()
