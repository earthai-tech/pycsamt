# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Focused coverage tests for pycsamt.models.mare2dem.edi.

Exercises `_latlon` header/attribute fallbacks, the `stations_from_edi`
skip paths (missing Z, missing coords, mismatched frequency tables,
missing error blocks), the confidence-weighting helpers in isolation,
and `make_mt_data_from_edi`'s pass-through to the ZMM writer.
"""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from pycsamt.models.mare2dem import edi as edi_mod


# ==========================================================================
# _latlon
# ==========================================================================


class TestLatLon:
    def test_header_section_wins(self):
        head = SimpleNamespace(lat=10.0, long=20.0, elev=300.0)
        ed = SimpleNamespace(get_section=lambda *_: head)
        assert edi_mod._latlon(ed) == (10.0, 20.0)

    def test_falls_back_to_plain_attrs_when_no_header(self):
        # No get_section -> get_coords returns NaN, so the function must
        # fall back to probing lat/lon attributes directly (real bug fix:
        # get_coords never returns None, it returns a NaN-filled _Coord,
        # so the old "if c is not None" guard never reached this branch).
        ed = SimpleNamespace(lat=1.5, lon=2.5)
        assert edi_mod._latlon(ed) == (1.5, 2.5)

    def test_falls_back_to_latitude_longitude_attrs(self):
        ed = SimpleNamespace(latitude=3.5, longitude=4.5)
        assert edi_mod._latlon(ed) == (3.5, 4.5)

    def test_falls_back_to_edi_wrapper_attrs(self):
        ed = SimpleNamespace(edi=SimpleNamespace(lat=5.0, lon=6.0))
        assert edi_mod._latlon(ed) == (5.0, 6.0)

    def test_non_numeric_attrs_are_skipped(self):
        ed = SimpleNamespace(lat="north", lon=6.0, edi=SimpleNamespace(lat=7.0, lon=8.0))
        assert edi_mod._latlon(ed) == (7.0, 8.0)

    def test_nothing_available_returns_nan(self):
        ed = SimpleNamespace()
        lat, lon = edi_mod._latlon(ed)
        assert np.isnan(lat) and np.isnan(lon)

    def test_get_coords_exception_falls_back_to_attrs(self):
        def _boom(*_a, **_kw):
            raise RuntimeError("no head section")

        ed = SimpleNamespace(get_section=_boom, lat=9.0, lon=10.0)
        assert edi_mod._latlon(ed) == (9.0, 10.0)

    def test_get_coords_itself_raising_falls_back_to_attrs(self, monkeypatch):
        # get_coords normally swallows its own internal errors (returns a
        # NaN-filled _Coord); force the imported name itself to raise so the
        # outer except in _latlon is exercised.
        from pycsamt.site import utils as site_utils

        def _raise(*_a, **_kw):
            raise RuntimeError("boom")

        monkeypatch.setattr(site_utils, "get_coords", _raise)
        ed = SimpleNamespace(lat=11.0, lon=12.0)
        assert edi_mod._latlon(ed) == (11.0, 12.0)


# ==========================================================================
# stations_from_edi
# ==========================================================================


def _edi_like(
    station,
    lat,
    lon,
    freqs,
    zxy,
    zyx,
    z_err=None,
):
    n = len(freqs)
    z = np.zeros((n, 2, 2), dtype=np.complex128)
    z[:, 0, 1] = zxy
    z[:, 1, 0] = zyx
    kw = dict(z=z, freq=np.asarray(freqs, dtype=float))
    if z_err is not None:
        kw["z_err"] = z_err
    return SimpleNamespace(
        station=station,
        lat=lat,
        lon=lon,
        Z=SimpleNamespace(**kw),
    )


class TestStationsFromEdi:
    def _patch_sites(self, monkeypatch, stations):
        from pycsamt.emtools import _core

        monkeypatch.setattr(_core, "ensure_sites", lambda source, **_: source)
        return stations

    def test_basic_two_stations(self, monkeypatch):
        freqs = [10.0, 1.0]
        s1 = _edi_like("S1", 1.0, 2.0, freqs, [1 + 1j, 2 + 2j], [1 - 1j, 2 - 2j])
        s2 = _edi_like("S2", 3.0, 4.0, freqs, [1 + 1j, 2 + 2j], [1 - 1j, 2 - 2j])
        self._patch_sites(monkeypatch, [s1, s2])

        out = edi_mod.stations_from_edi([s1, s2])
        assert [st.name for st in out] == ["S1", "S2"]
        assert out[0].latitude == 1.0 and out[0].longitude == 2.0
        # periods sorted ascending
        assert np.all(np.diff(out[0].periods) > 0)
        # default rel error applied (no z_err on fixture)
        np.testing.assert_allclose(out[0].phase_te_se, np.degrees(0.05))

    def test_missing_z_block_is_skipped(self, monkeypatch):
        freqs = [10.0, 1.0]
        good = _edi_like("GOOD", 1.0, 2.0, freqs, [1 + 1j, 2 + 2j], [1 - 1j, 2 - 2j])
        bad = SimpleNamespace(station="BAD", lat=5.0, lon=6.0)  # no Z at all
        self._patch_sites(monkeypatch, [good, bad])

        out = edi_mod.stations_from_edi([good, bad])
        assert [st.name for st in out] == ["GOOD"]

    def test_missing_coords_is_skipped(self, monkeypatch):
        freqs = [10.0, 1.0]
        good = _edi_like("GOOD", 1.0, 2.0, freqs, [1 + 1j, 2 + 2j], [1 - 1j, 2 - 2j])
        nocoord = _edi_like("NOCOORD", float("nan"), float("nan"), freqs, [1j], [1j])
        self._patch_sites(monkeypatch, [good, nocoord])

        out = edi_mod.stations_from_edi([good, nocoord])
        assert [st.name for st in out] == ["GOOD"]

    def test_mismatched_frequency_table_is_skipped(self, monkeypatch):
        ref = _edi_like(
            "REF", 1.0, 2.0, [10.0, 1.0], [1 + 1j, 2 + 2j], [1 - 1j, 2 - 2j]
        )
        odd = _edi_like(
            "ODD",
            3.0,
            4.0,
            [10.0, 5.0, 1.0],
            [1 + 1j, 1 + 1j, 2 + 2j],
            [1 - 1j, 1 - 1j, 2 - 2j],
        )
        self._patch_sites(monkeypatch, [ref, odd])

        out = edi_mod.stations_from_edi([ref, odd])
        assert [st.name for st in out] == ["REF"]

    def test_uses_stored_z_err_when_present(self, monkeypatch):
        freqs = [10.0, 1.0]
        z_err = np.zeros((2, 2, 2))
        z_err[:, 0, 1] = 0.5  # 50% relative error on Zxy
        z_err[:, 1, 0] = 0.2
        s = _edi_like(
            "S1", 1.0, 2.0, freqs, [1 + 0j, 2 + 0j], [1 + 0j, 2 + 0j], z_err=z_err
        )
        self._patch_sites(monkeypatch, [s])

        out = edi_mod.stations_from_edi([s])
        # |zxy| = [1, 2] (already ascending-period order), z_err Zxy = 0.5
        rel_te = np.array([0.5 / 1.0, 0.5 / 2.0])
        np.testing.assert_allclose(out[0].phase_te_se, np.degrees(rel_te))

    def test_degenerate_stored_error_falls_back_to_default(self, monkeypatch):
        freqs = [10.0, 1.0]
        z_err = np.zeros((2, 2, 2))  # all zero -> degenerate, falls back
        s = _edi_like(
            "S1", 1.0, 2.0, freqs, [1 + 0j, 2 + 0j], [1 + 0j, 2 + 0j], z_err=z_err
        )
        self._patch_sites(monkeypatch, [s])

        out = edi_mod.stations_from_edi([s], default_rel_error=0.1)
        np.testing.assert_allclose(out[0].phase_te_se, np.degrees(0.1))

    def test_unwrap_failure_falls_back_to_raw_ed(self, monkeypatch):
        from pycsamt.emtools import _core

        freqs = [10.0, 1.0]
        s = _edi_like("S1", 1.0, 2.0, freqs, [1 + 1j, 2 + 2j], [1 - 1j, 2 - 2j])
        self._patch_sites(monkeypatch, [s])

        def _raise(_ed):
            raise RuntimeError("cannot unwrap")

        monkeypatch.setattr(_core, "_unwrap", _raise)
        out = edi_mod.stations_from_edi([s])
        assert [st.name for st in out] == ["S1"]

    def test_all_skipped_raises_value_error(self, monkeypatch):
        bad = SimpleNamespace(station="BAD")
        self._patch_sites(monkeypatch, [bad])
        with pytest.raises(ValueError, match="No EDI station"):
            edi_mod.stations_from_edi([bad])

    def test_confidence_weighting_pipeline_end_to_end(self, monkeypatch):
        freqs = [10.0, 1.0]
        s = _edi_like("S1", 1.0, 2.0, freqs, [1 + 1j], [1 - 1j])
        s.Z.z = np.zeros((2, 2, 2), dtype=np.complex128)
        s.Z.z[:, 0, 1] = [1 + 1j, 2 + 2j]
        s.Z.z[:, 1, 0] = [1 - 1j, 2 - 2j]
        self._patch_sites(monkeypatch, [s])
        monkeypatch.setattr(
            edi_mod,
            "_frequency_confidence_lookup",
            lambda *_a, **_kw: {"S1": (np.array([10.0, 1.0]), np.array([0.5, 1.0]))},
        )
        out = edi_mod.stations_from_edi(
            [s], confidence_weighting=True, confidence_min=0.05, confidence_power=1.0
        )
        # freq=1.0 (index 1 after sort by ascending period == descending freq)
        # has CR=1.0 -> factor 1 -> unaffected; freq=10 has CR=0.5 -> factor 2
        assert out[0].phase_te_se[0] > out[0].phase_te_se[1]


# ==========================================================================
# _frequency_confidence_lookup
# ==========================================================================


class TestFrequencyConfidenceLookup:
    def test_none_table_returns_empty(self, monkeypatch):
        from pycsamt.emtools import qc

        monkeypatch.setattr(qc, "frequency_confidence_table", lambda *a, **kw: None)
        out = edi_mod._frequency_confidence_lookup(
            [], method="composite", weights=None
        )
        assert out == {}

    def test_empty_table_returns_empty(self, monkeypatch):
        from pycsamt.emtools import qc

        monkeypatch.setattr(
            qc, "frequency_confidence_table", lambda *a, **kw: pd.DataFrame()
        )
        out = edi_mod._frequency_confidence_lookup(
            [], method="composite", weights=None
        )
        assert out == {}

    def test_groups_and_filters_by_station(self, monkeypatch):
        from pycsamt.emtools import qc

        df = pd.DataFrame(
            {
                "station": ["S1", "S1", "S1", "S2", "S3"],
                "frequency_hz": [10.0, 1.0, -1.0, 5.0, -3.0],
                "confidence": [0.5, 1.5, 0.9, -0.2, 0.7],
            }
        )
        monkeypatch.setattr(qc, "frequency_confidence_table", lambda *a, **kw: df)
        out = edi_mod._frequency_confidence_lookup(
            [], method="composite", weights=None
        )
        # S3's only row has a negative frequency -> entirely filtered out,
        # so it must not appear in the lookup at all.
        assert set(out) == {"S1", "S2"}
        freq_s1, conf_s1 = out["S1"]
        # the negative-frequency row must be excluded
        assert -1.0 not in freq_s1
        assert len(freq_s1) == 2
        # confidence clipped to [0, 1]
        assert conf_s1.max() <= 1.0
        freq_s2, conf_s2 = out["S2"]
        assert conf_s2[0] == 0.0  # -0.2 clipped to 0


# ==========================================================================
# _station_confidence_values
# ==========================================================================


class TestStationConfidenceValues:
    def test_missing_station_returns_ones(self):
        out = edi_mod._station_confidence_values({}, "S1", np.array([1.0, 2.0]))
        np.testing.assert_allclose(out, [1.0, 1.0])

    def test_exact_match_used(self):
        lookup = {"S1": (np.array([1.0, 10.0]), np.array([0.3, 0.9]))}
        out = edi_mod._station_confidence_values(
            lookup, "S1", np.array([10.0, 1.0])
        )
        np.testing.assert_allclose(out, [0.9, 0.3])

    def test_non_matching_frequency_stays_one(self):
        lookup = {"S1": (np.array([1.0]), np.array([0.3]))}
        out = edi_mod._station_confidence_values(lookup, "S1", np.array([999.0]))
        np.testing.assert_allclose(out, [1.0])

    def test_invalid_frequency_skipped(self):
        lookup = {"S1": (np.array([1.0]), np.array([0.3]))}
        out = edi_mod._station_confidence_values(
            lookup, "S1", np.array([0.0, float("nan")])
        )
        np.testing.assert_allclose(out, [1.0, 1.0])


# ==========================================================================
# _confidence_error_factor
# ==========================================================================


class TestConfidenceErrorFactor:
    def test_basic_scaling(self):
        out = edi_mod._confidence_error_factor(
            np.array([1.0, 0.5, 0.25]), confidence_min=0.1, confidence_power=1.0
        )
        np.testing.assert_allclose(out, [1.0, 2.0, 4.0])

    def test_floor_applied(self):
        out = edi_mod._confidence_error_factor(
            np.array([0.001]), confidence_min=0.1, confidence_power=1.0
        )
        np.testing.assert_allclose(out, [10.0])

    def test_power_zero_is_neutral(self):
        out = edi_mod._confidence_error_factor(
            np.array([0.2]), confidence_min=0.05, confidence_power=0.0
        )
        np.testing.assert_allclose(out, [1.0])

    def test_nonfinite_replaced_with_one(self):
        out = edi_mod._confidence_error_factor(
            np.array([float("nan")]), confidence_min=0.1, confidence_power=1.0
        )
        np.testing.assert_allclose(out, [1.0])


# ==========================================================================
# make_mt_data_from_edi
# ==========================================================================


class TestMakeMtDataFromEdi:
    def test_passes_through_to_writer(self, monkeypatch, tmp_path):
        captured = {}

        def _fake_stations_from_edi(source, **kw):
            captured["stations_kw"] = kw
            return ["STATION_SENTINEL"]

        def _fake_writer(stations, out_file, **kw):
            captured["stations"] = stations
            captured["out_file"] = out_file
            captured["writer_kw"] = kw
            return "EMDATA_SENTINEL"

        monkeypatch.setattr(edi_mod, "stations_from_edi", _fake_stations_from_edi)
        monkeypatch.setattr(
            edi_mod, "make_mt_data_from_stations", _fake_writer
        )

        out = edi_mod.make_mt_data_from_edi(
            "some/source",
            tmp_path / "out.emdata",
            error_floor_te=0.1,
            error_floor_tm=0.2,
            confidence_weighting=True,
            confidence_min=0.2,
            utm_zone="19N",
        )

        assert out == "EMDATA_SENTINEL"
        assert captured["stations"] == ["STATION_SENTINEL"]
        assert captured["out_file"] == tmp_path / "out.emdata"
        assert captured["writer_kw"]["error_floor_te"] == 0.1
        assert captured["writer_kw"]["error_floor_tm"] == 0.2
        assert captured["writer_kw"]["utm_zone"] == "19N"
        assert captured["stations_kw"]["confidence_weighting"] is True
        assert captured["stations_kw"]["confidence_min"] == 0.2
