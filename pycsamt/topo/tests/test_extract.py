# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for pycsamt.topo.extract using real WILLY AMT EDI data.

The two WILLY profiles committed to the repository are exercised.
Tests verify:
  * elevation extraction from EDIFile list and EDICollection
  * monotonic along-profile chainage
  * station naming
  * fallback behaviour for zero / missing elevation
"""

from __future__ import annotations

import glob
import os

import numpy as np
import pytest

from pycsamt.seg.collection import EDICollection
from pycsamt.seg.edi import EDIFile
from pycsamt.topo.extract import (
    extract_chainage,
    extract_elevation,
    extract_station_names,
    has_elevation,
)

# ---------------------------------------------------------------------------
# Paths to real data
# ---------------------------------------------------------------------------

_DATA_ROOT = os.path.join(
    os.path.dirname(__file__),
    "..",
    "..",
    "..",
    "data",
    "AMT",
    "WILLY_DATA",
)
_CSAMT_ROOT = os.path.join(
    os.path.dirname(__file__), "..", "..", "..", "data", "CSAMT"
)


def _load_profile(profile: str) -> list[EDIFile]:
    edi_dir = os.path.join(_DATA_ROOT, profile)
    paths = sorted(glob.glob(os.path.join(edi_dir, "*.edi")))
    if not paths:
        pytest.skip(f"WILLY {profile} EDI files not found at {edi_dir}")
    return [EDIFile(p) for p in paths]


def _load_csamt() -> list[EDIFile]:
    paths = sorted(glob.glob(os.path.join(_CSAMT_ROOT, "*.edi")))
    if not paths:
        pytest.skip(f"CSAMT EDI files not found at {_CSAMT_ROOT}")
    return [EDIFile(p) for p in paths]


# Shared profiles
L18 = None
COLL18 = None


def _get_l18():
    global L18, COLL18
    if L18 is None:
        L18 = _load_profile("L18PLT")
        COLL18 = EDICollection(L18)
    return L18, COLL18


# ---------------------------------------------------------------------------
# has_elevation
# ---------------------------------------------------------------------------


class TestHasElevation:
    def test_true_for_edi_file_list(self):
        edis, _ = _get_l18()
        assert has_elevation(edis) is True

    def test_true_for_edi_collection(self):
        _, coll = _get_l18()
        assert has_elevation(coll) is True

    def test_true_for_single_edi(self):
        edis, _ = _get_l18()
        assert has_elevation(edis[0]) is True

    def test_false_for_empty_list(self):
        assert has_elevation([]) is False

    def test_false_warns_for_zero_elevation(self):
        """A list of EDIFile objects with no lat/lon/elev should warn."""

        class _FakeEDI:
            """Minimal stand-in with zero elevation and no head."""

            station = "FAKE001"

        # extract_elevation emits the warning; has_elevation suppresses it
        with pytest.warns(UserWarning, match="elevation"):
            extract_elevation([_FakeEDI()])


# ---------------------------------------------------------------------------
# extract_elevation
# ---------------------------------------------------------------------------


class TestExtractElevation:
    def test_length_matches_station_count(self):
        edis, _ = _get_l18()
        assert len(extract_elevation(edis)) == len(edis)

    def test_first_station_value(self):
        """18-001A is documented at 99 m a.s.l."""
        edis, _ = _get_l18()
        assert extract_elevation(edis)[0] == pytest.approx(99.0)

    def test_all_positive(self):
        edis, _ = _get_l18()
        assert np.all(extract_elevation(edis) > 0)

    def test_range_l18(self):
        """L18PLT elevation range is 37–144 m."""
        edis, _ = _get_l18()
        elev = extract_elevation(edis)
        assert elev.min() >= 30.0
        assert elev.max() <= 200.0

    def test_collection_equals_list(self):
        edis, coll = _get_l18()
        np.testing.assert_array_equal(
            extract_elevation(edis),
            extract_elevation(coll),
        )

    def test_dtype_float(self):
        edis, _ = _get_l18()
        assert extract_elevation(edis).dtype == np.float64

    def test_no_nan_values(self):
        edis, _ = _get_l18()
        assert np.all(np.isfinite(extract_elevation(edis)))

    def test_empty_list_returns_empty(self):
        arr = extract_elevation([])
        assert len(arr) == 0

    def test_all_committed_profiles(self):
        for prof in ("L18PLT", "L22PLT"):
            edis = _load_profile(prof)
            elev = extract_elevation(edis)
            assert len(elev) == len(edis)
            assert np.all(elev > 0), f"{prof}: found zero/negative elevations"

    def test_zero_elev_emits_warning(self):
        """Synthetic EDI with elev=0 should trigger UserWarning."""

        class _ZeroElev:
            class _Head:
                elev = 0.0
                lat = 0.0
                lon = 0.0

            Head = _Head()

        with pytest.warns(UserWarning, match="elevation"):
            extract_elevation([_ZeroElev()])


# ---------------------------------------------------------------------------
# >INFO-section fallback (real CSAMT data: LATITUDE=/LONGITUDE=/ELEVATION=
# recorded under >INFO instead of the standard >HEAD keys)
# ---------------------------------------------------------------------------


class TestInfoSectionFallback:
    """Real Tongkeng CSAMT delivery: no LAT=/LONG=/ELEV= under >HEAD at
    all, but real values under >INFO, already parsed onto
    ``pycsamt.seg.heads.Info``'s own ``.latitude``/``.longitude``/
    ``.elevation`` attributes -- confirmed directly against the raw EDI
    text (``grep ELEVATION data/CSAMT/csa000.edi`` shows ``ELEVATION=573.4``).
    A standards-compliant fallback, not a delivery-specific hack: the
    SEG-EDI spec permits location metadata in ``>INFO``.
    """

    def test_has_elevation_true_via_info_section(self):
        edis = _load_csamt()
        assert has_elevation(edis) is True

    def test_extract_elevation_matches_raw_info_block(self):
        """csa000/csa450 documented directly in the raw EDI text."""
        edis = _load_csamt()
        elev = extract_elevation(edis)
        assert elev[0] == pytest.approx(573.4)
        assert elev[-1] == pytest.approx(510.3)

    def test_extract_elevation_no_zero_shadowing(self):
        """A >HEAD Location defaulting to elevation=0.0 (because the file
        never set it) must not shadow a real, non-zero >INFO value.
        """
        edis = _load_csamt()
        elev = extract_elevation(edis)
        assert np.all(elev > 0.0)

    def test_extract_chainage_via_info_latlon(self):
        """extract_chainage also needs >INFO's latitude/longitude to work
        for this delivery -- station spacing here is a real, known 50 m.
        """
        edis = _load_csamt()
        chain_m = extract_chainage(edis) * 1000.0
        spacings = np.diff(chain_m)
        assert np.allclose(spacings, 50.0, atol=5.0)

    def test_info_section_does_not_override_a_real_head_value(self):
        """A HEAD-provided non-zero elevation must win outright -- >INFO
        is a fallback for missing/zero HEAD data, not an override.
        """

        class _Info:
            elevation = 999.0
            latitude = 99.0
            longitude = 99.0

        class _Location:
            elevation = 250.0

        class _Head:
            Location = _Location()

        class _EDI:
            Head = _Head()

            def get_section(self, name):
                return _Info() if name == "info" else None

        elev = extract_elevation([_EDI()])
        assert elev[0] == pytest.approx(250.0)

    def test_info_section_used_when_head_elevation_is_zero_default(self):
        """A HEAD Location left at its 0.0 default (never set by the
        file) must fall through to a real >INFO value.
        """

        class _Info:
            elevation = 573.4

        class _Location:
            elevation = 0.0

        class _Head:
            Location = _Location()

        class _EDI:
            Head = _Head()

            def get_section(self, name):
                return _Info() if name == "info" else None

        elev = extract_elevation([_EDI()])
        assert elev[0] == pytest.approx(573.4)

    def test_falls_back_to_zero_when_info_also_missing(self):
        """No HEAD, no INFO: still degrades to the documented all-zero /
        warning behaviour rather than raising.
        """

        class _EDI:
            pass

        with pytest.warns(UserWarning, match="elevation"):
            elev = extract_elevation([_EDI()])
        assert elev[0] == 0.0


# ---------------------------------------------------------------------------
# extract_chainage
# ---------------------------------------------------------------------------


class TestExtractChainage:
    def test_length_matches_station_count(self):
        edis, _ = _get_l18()
        chain = extract_chainage(edis)
        assert len(chain) == len(edis)

    def test_starts_at_zero(self):
        edis, _ = _get_l18()
        assert extract_chainage(edis)[0] == pytest.approx(0.0)

    def test_monotonically_non_decreasing(self):
        edis, _ = _get_l18()
        chain = extract_chainage(edis)
        assert np.all(np.diff(chain) >= 0.0)

    def test_strictly_increasing_for_distinct_stations(self):
        """All L18PLT stations are at distinct positions."""
        edis, _ = _get_l18()
        chain = extract_chainage(edis)
        assert np.all(np.diff(chain) > 0.0)

    def test_total_profile_length_reasonable(self):
        """L18PLT spans ~2.4 km along profile."""
        edis, _ = _get_l18()
        chain = extract_chainage(edis)
        assert 1.5 < chain[-1] < 5.0, f"unexpected total: {chain[-1]:.2f} km"

    def test_dtype_float(self):
        edis, _ = _get_l18()
        assert extract_chainage(edis).dtype == np.float64

    def test_collection_equals_list(self):
        edis, coll = _get_l18()
        np.testing.assert_allclose(
            extract_chainage(edis),
            extract_chainage(coll),
            rtol=1e-9,
        )

    def test_empty_list_returns_empty(self):
        assert len(extract_chainage([])) == 0

    def test_single_station_returns_zero(self):
        edis, _ = _get_l18()
        chain = extract_chainage([edis[0]])
        assert chain.shape == (1,)
        assert chain[0] == pytest.approx(0.0)

    def test_all_committed_profiles_have_reasonable_lengths(self):
        """Every WILLY profile should span at least 1 km and at most 15 km.

        L22PLT has one large gap (~2.1 km) due to a station with a large
        lon offset, giving a total chainage of ~10.5 km.
        """
        for prof in ("L18PLT", "L22PLT"):
            edis = _load_profile(prof)
            chain = extract_chainage(edis)
            total_km = chain[-1]
            assert (
                1.0 < total_km < 15.0
            ), f"{prof}: total chainage {total_km:.2f} km out of expected range"

    def test_inter_station_spacing_realistic(self):
        """WILLY typical (median) inter-station spacing is ~80–200 m.

        A few near-duplicate stations in L18PLT produce sub-10 m minimum
        spacings — median is used to characterise the typical pattern.
        """
        edis, _ = _get_l18()
        chain = extract_chainage(edis)
        spacings_m = np.diff(chain) * 1000.0
        # The majority of spacings are in the 50–200 m range
        assert np.median(spacings_m) > 50.0
        assert np.median(spacings_m) < 200.0
        assert spacings_m.max() < 500.0


# ---------------------------------------------------------------------------
# extract_station_names
# ---------------------------------------------------------------------------


class TestExtractStationNames:
    def test_length_matches_station_count(self):
        edis, _ = _get_l18()
        assert len(extract_station_names(edis)) == len(edis)

    def test_first_station_name(self):
        edis, _ = _get_l18()
        names = extract_station_names(edis)
        assert names[0] == "23-18-001A"

    def test_all_strings(self):
        edis, _ = _get_l18()
        names = extract_station_names(edis)
        assert all(isinstance(n, str) for n in names)

    def test_no_empty_strings(self):
        edis, _ = _get_l18()
        names = extract_station_names(edis)
        assert all(len(n) > 0 for n in names)

    def test_unique_within_profile(self):
        edis, _ = _get_l18()
        names = extract_station_names(edis)
        assert len(names) == len(set(names)), "Duplicate station names in L18PLT"

    def test_collection_equals_list(self):
        edis, coll = _get_l18()
        assert extract_station_names(edis) == extract_station_names(coll)

    def test_fallback_for_missing_name(self):
        """If dataid is missing, should fall back to S000-style names."""

        class _NoName:
            class Head:
                dataid = None
                station = None
                lat = 1.0
                lon = 1.0
                elev = 100.0

        names = extract_station_names([_NoName()])
        assert names == ["S000"]

    def test_all_committed_profiles(self):
        for prof in ("L18PLT", "L22PLT"):
            edis = _load_profile(prof)
            names = extract_station_names(edis)
            assert len(names) == len(edis)
            assert all(isinstance(n, str) and len(n) > 0 for n in names)


# ---------------------------------------------------------------------------
# Edge cases / robustness
# ---------------------------------------------------------------------------


class TestEdgeCases:
    def test_extract_elevation_single_edi(self):
        edis, _ = _get_l18()
        elev = extract_elevation(edis[0])
        assert elev.shape == (1,)
        assert elev[0] == pytest.approx(99.0)

    def test_extract_chainage_two_stations(self):
        edis, _ = _get_l18()
        chain = extract_chainage(edis[:2])
        assert chain[0] == pytest.approx(0.0)
        assert chain[1] > 0.0

    def test_old_head_attribute_style(self):
        """Verify compatibility with old MTpy-style .Head attribute."""

        class _OldHead:
            elev = 250.0
            lat = 10.0
            lon = 20.0
            dataid = "OLD001"

        class _OldEDI:
            Head = _OldHead()

        assert has_elevation([_OldEDI()])
        elev = extract_elevation([_OldEDI()])
        assert elev[0] == pytest.approx(250.0)
        names = extract_station_names([_OldEDI()])
        assert names[0] == "OLD001"

    def test_edi_collection_all_profiles(self):
        """Load each profile into EDICollection and extract successfully."""
        for prof in ("L18PLT", "L22PLT"):
            edis = _load_profile(prof)
            coll = EDICollection(edis)
            elev = extract_elevation(coll)
            chain = extract_chainage(coll)
            assert len(elev) == len(edis)
            assert len(chain) == len(edis)
            assert has_elevation(coll)
