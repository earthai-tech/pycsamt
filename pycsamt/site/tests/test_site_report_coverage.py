"""Behavior and fallback coverage for site reports."""

from types import SimpleNamespace

import numpy as np
import pytest

import pycsamt.site.report as report_module
from pycsamt.site.report import SiteReport, SitesReport


class FakeSite:
    def __init__(self, name="S01", populated=True):
        self.name = name
        self.coords = (5.25, -3.75, 120.0)
        self.freq = np.array([0.001, 1.0, 2000.0]) if populated else None
        self.z = (
            np.array(
                [
                    [[1 + 1j, 2 + 2j], [3 + 3j, 4 + 4j]],
                    [[np.nan + 1j, 2 + 2j], [3 + 3j, 4 + 4j]],
                    [[1 + 1j, 2 + 2j], [3 + 3j, 4 + 4j]],
                ]
            )
            if populated
            else None
        )
        self.rho = np.arange(1.0, 13.0).reshape(3, 2, 2) if populated else None
        self.phase = np.arange(12.0).reshape(3, 2, 2) if populated else None
        self.tipper = np.ones((3, 2)) if populated else None

    def has_component(self, name):
        return self.z is not None and name != "zyy"

    def to_dataframe(self, kind, api=None):
        return (kind, api)


def test_format_and_array_helpers_cover_boundaries():
    assert len(report_module._bar(-1, 4)) == 4
    assert len(report_module._bar(2, 4)) == 4
    assert report_module._fmt_rho(None, None)
    assert "10" in report_module._fmt_rho(10, 2)
    assert report_module._fmt_phi(None, None)
    assert "10.0" in report_module._fmt_phi(10, 2)
    for value in (None, 2000, 2, 0.2, 0.0002):
        assert report_module._fmt_freq(value)
    assert report_module._check(True) != report_module._check(False)
    assert report_module._safe_arr(None) is None
    assert report_module._safe_arr([]) is None
    assert report_module._safe_arr([1 + 2j]).tolist() == [1.0]


@pytest.mark.parametrize("helper", [report_module._rho_stats, report_module._phase_stats])
def test_stats_helpers_accept_layouts_and_bad_inputs(helper):
    cube = np.arange(1.0, 13.0).reshape(3, 2, 2)
    matrix = cube.reshape(3, 4)
    assert helper(cube, 1)[0] is not None
    assert helper(matrix, 1)[0] is not None
    assert helper(np.array([1.0]), 0) == (None, None)
    assert helper(np.array([["bad"]]), 0) == (None, None)


def test_quality_and_component_fallbacks():
    z = np.ones((2, 4), dtype=complex)
    z[0, 1] = complex(np.nan, 1)
    assert report_module._quality_pct(z, 1) == 0.5
    assert report_module._quality_pct(np.ones((2, 2, 2)), 1) == 1.0
    assert report_module._quality_pct(None, 0) is None
    assert report_module._quality_pct(np.ones(2), 0) is None
    assert report_module._quality_pct(np.array([["bad"]]), 0) is None

    fallback = SimpleNamespace(z=np.ones((1, 2, 2)))
    assert report_module._has_component(fallback, "zxy") is True
    assert report_module._has_component(SimpleNamespace(z=None), "zxy") is False


def test_site_report_populated_and_plain_output(monkeypatch, capsys):
    site = FakeSite()
    report = SiteReport(site)
    stats = report.to_dict()
    assert stats["nfreq"] == 3
    assert stats["has_tipper"] is True
    assert stats["quality"]["Zxx"] == pytest.approx(2 / 3)
    assert "S01" in report.summary()
    assert repr(report) == report.summary()
    assert report.to_dataframe("z", api=False) == ("z", False)
    assert any("Coordinates" in line for line in report._plain_lines())

    monkeypatch.setattr(report_module, "_RICH", False)
    report.report(detail=True)
    assert "Site: S01" in capsys.readouterr().out
    assert report_module._console() is None


def test_site_report_handles_missing_and_broken_site():
    site = FakeSite(populated=False)
    site.coords = property(lambda self: (_ for _ in ()).throw(ValueError()))
    report = SiteReport(site)
    assert report.to_dict()["nfreq"] == 0
    assert any("Coordinates" in line for line in report._plain_lines())


def test_sites_report_statistics_exports_and_plain_output(monkeypatch, capsys):
    sites = [FakeSite("A"), FakeSite("B", populated=False)]
    report = SitesReport(sites)
    assert report.to_dict()[0]["name"] == "A"
    assert len(report.to_dataframe(api=False)) == 2
    assert "2 stations" in report.summary()
    assert repr(report) == report.summary()
    assert any("Survey: 2" in line for line in report._plain_lines(report._records))

    monkeypatch.setattr(report_module, "_RICH", False)
    report.report(top=1, detail=True)
    assert "Survey: 2" in capsys.readouterr().out


@pytest.mark.skipif(not report_module._RICH, reason="rich is optional")
def test_rich_renderers_accept_populated_and_empty_reports():
    site_report = SiteReport(FakeSite())
    site_report._rich_report(detail=True)
    sites_report = SitesReport([FakeSite("A"), FakeSite("B", False)])
    sites_report._rich_report(sites_report._records[:1], detail=True)
    SitesReport([])._rich_report([], detail=True)


def test_empty_sites_report_has_safe_summary():
    report = SitesReport([])
    assert report.to_dict() == []
    assert "0 stations" in report.summary()
