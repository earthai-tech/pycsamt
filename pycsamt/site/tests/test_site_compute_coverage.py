# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Branch-level coverage for ``pycsamt.site.compute``.

Complements ``test_site_compute.py`` (which exercises the public API
against a real simulated EDI) with targeted tests against small
duck-typed stubs for the defensive branches and private helpers that
a well-formed EDI never triggers: missing/invalid ``Z``/tipper
sections, empty-frequency arrays, NaN-poisoned bands, and the
``Sites``-as-input path through ``_as_sites_iter``.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from pycsamt.seg.edi import EDIFile
from pycsamt.site import compute as cmp
from pycsamt.site import edit as ed
from pycsamt.site.base import Sites


class _EdiNoZ:
    """Duck-typed EDI-like object exposing ``Z=None``."""

    Z = None

    def get_section(self, name=None):
        return None


class _ZHolderRaisingThenValid:
    """``Z`` object whose ``z`` attribute raises; ``impedance`` is valid."""

    def __init__(self, impedance):
        self.impedance = impedance

    @property
    def z(self):
        raise RuntimeError("boom")


class _ZHolderWrongShapeThenValid:
    """``z`` has the wrong shape; ``impedance`` is a valid (n, 2, 2)."""

    def __init__(self, bad, good):
        self.z = bad
        self.impedance = good


class _ZHolderAllNone:
    z = None
    impedance = None
    _z = None


class _EdiWithZ:
    def __init__(self, zobj, freq=None):
        self.Z = zobj
        if freq is not None:
            self.Z.freq = freq

    def get_section(self, name=None):
        return None


def _load_edi(p):
    return EDIFile(p)


# ---------------------------------------------------------------------------
# strike_estimate — missing-Z / method branches
# ---------------------------------------------------------------------------


def test_strike_estimate_returns_nan_when_z_missing():
    ang = cmp.strike_estimate(_EdiNoZ(), method="swift")
    assert isinstance(ang, float)
    assert math.isnan(ang)


def test_strike_estimate_phase_diff_method(tmp_path, simulated_edi):
    edi = _load_edi(simulated_edi)
    edi = ed.fill_missing(edi, how="zero", components=("Z",), inplace=False)
    ang = cmp.strike_estimate(edi, method="phase_diff")
    assert ang in (0.0, 90.0)


def test_strike_estimate_unknown_method_falls_back_to_swift(simulated_edi):
    edi = _load_edi(simulated_edi)
    edi = ed.fill_missing(edi, how="zero", components=("Z",), inplace=False)
    ang_default = cmp.strike_estimate(edi, method="swift")
    ang_unknown = cmp.strike_estimate(edi, method="bogus")
    assert ang_default == pytest.approx(ang_unknown)


def test_phase_diff_theta_nan_when_all_values_nan():
    Z = np.full((3, 2, 2), np.nan, dtype=complex)
    assert math.isnan(cmp._phase_diff_theta(Z))


def test_swift_cost_and_theta_with_empty_frequency_axis():
    Z = np.empty((0, 2, 2), dtype=complex)
    assert cmp._swift_cost(Z, 10.0) == float("inf")
    assert math.isnan(cmp._swift_theta(Z))


# ---------------------------------------------------------------------------
# res_at_freq — missing data and non-positive frequency
# ---------------------------------------------------------------------------


def test_res_at_freq_returns_nan_row_when_z_missing():
    out = cmp.res_at_freq(_EdiNoZ(), 100.0, how="nearest")
    assert math.isnan(out["res_xy"])
    assert math.isnan(out["res_yx"])
    assert math.isnan(out["f_used"])


def test_rho_from_z_non_positive_or_nonfinite_frequency():
    assert math.isnan(cmp._rho_from_z(1 + 1j, 0.0))
    assert math.isnan(cmp._rho_from_z(1 + 1j, -5.0))
    assert math.isnan(cmp._rho_from_z(1 + 1j, float("nan")))
    assert cmp._rho_from_z(1 + 1j, 10.0) > 0.0


def test_interp_returns_nan_on_shape_error():
    out = cmp._interp(np.array([1.0, 2.0]), np.array([[1.0, 2.0], [3.0, 4.0]]), 1.5)
    assert math.isnan(out)


# ---------------------------------------------------------------------------
# phase_slope — missing data, empty band, and NaN-poisoned fits
# ---------------------------------------------------------------------------


def test_phase_slope_returns_nan_row_when_z_missing():
    out = cmp.phase_slope(_EdiNoZ(), band=(1.0, 100.0))
    assert math.isnan(out["slope_xy"])
    assert math.isnan(out["slope_yx"])


def test_phase_slope_returns_nan_when_band_has_no_overlap(simulated_edi):
    edi = _load_edi(simulated_edi)
    edi = ed.fill_missing(edi, how="zero", components=("Z",), inplace=False)
    out = cmp.phase_slope(edi, band=(1e6, 2e6))
    assert math.isnan(out["slope_xy"])
    assert math.isnan(out["slope_yx"])


def test_phase_slope_nan_when_fit_raises(monkeypatch):
    # An infinite frequency inside the band feeds an infinite value into
    # np.polyfit's design matrix, which raises LinAlgError internally --
    # unlike an all-NaN *phase* series, which polyfit tolerates and just
    # returns NaN coefficients without raising (a NaN frequency instead
    # of infinite would simply fail the band mask and be dropped). A
    # real ``Z`` section refuses to store a non-finite frequency, so a
    # lightweight stub is used to force the array through untouched.
    Z = np.array(
        [[[1 + 1j, 1 + 1j], [1 + 1j, 1 + 1j]]] * 2,
        dtype=complex,
    )

    class _E:
        Z = object()  # only needs to look edi-like for is_edi_file()

        def get_section(self, name=None):
            return None

    stub = _E()
    monkeypatch.setattr(cmp, "get_freq", lambda _e: np.array([1.0, np.inf]))
    monkeypatch.setattr(cmp, "_get_z", lambda _e: Z)

    out = cmp.phase_slope(stub, band=(1.0, float("inf")))
    assert math.isnan(out["slope_xy"])
    assert math.isnan(out["slope_yx"])


# ---------------------------------------------------------------------------
# tipper_magnitude — per_freq branches and the empty-input dead ends
# ---------------------------------------------------------------------------


class _NeitherEdiNorIterable:
    """Not edi-like and not iterable -> _as_sites_iter yields nothing."""


def test_tipper_magnitude_no_rows_summary_returns_nan_dict():
    out = cmp.tipper_magnitude(_NeitherEdiNorIterable(), per_freq=False)
    assert out == {"mean": np.nan, "median": np.nan, "max": np.nan}


def test_tipper_magnitude_no_rows_per_freq_returns_empty_arrays():
    out = cmp.tipper_magnitude(_NeitherEdiNorIterable(), per_freq=True)
    assert out["freq"].size == 0
    assert out["mag"].size == 0


def test_tipper_magnitude_per_freq_single_site(tmp_path, simulated_edi):
    edi = _load_edi(simulated_edi)
    edi = ed.fill_missing(edi, how="zero", components=("Tip",), inplace=False)
    out = cmp.tipper_magnitude(edi, per_freq=True)
    assert set(out.keys()) == {"freq", "mag"}
    assert len(out["freq"]) == len(out["mag"])
    assert len(out["freq"]) > 0


def test_tipper_magnitude_per_freq_multi_site_skips_missing(
    tmp_path, simulated_edi
):
    p1 = tmp_path / "TA.edi"
    p1.write_text(simulated_edi.read_text(encoding="utf-8"), encoding="utf-8")
    e1 = _load_edi(p1)
    e1 = ed.fill_missing(e1, how="zero", components=("Tip",), inplace=False)

    e2 = _EdiNoZ()  # no tipper at all -> "continue" branch under per_freq=True

    df = cmp.tipper_magnitude([e1, e2], per_freq=True, api=False)
    assert set(df["station"].unique()) == {cmp.station_name(e1)}


# ---------------------------------------------------------------------------
# _get_z — attribute lookup, exceptions, and shape mismatches
# ---------------------------------------------------------------------------


def test_get_z_returns_none_when_no_z_attribute():
    assert cmp._get_z(object()) is None


def test_get_z_returns_none_when_all_candidates_are_none():
    class _E:
        Z = _ZHolderAllNone()

    assert cmp._get_z(_E()) is None


def test_get_z_recovers_after_attribute_access_raises():
    good = np.zeros((2, 2, 2), complex)

    class _E:
        Z = _ZHolderRaisingThenValid(good)

    arr = cmp._get_z(_E())
    assert arr is not None
    assert arr.shape == (2, 2, 2)


def test_get_z_skips_wrong_shape_and_uses_next_candidate():
    bad = np.zeros((2, 3), float)
    good = np.zeros((2, 2, 2), complex)

    class _E:
        Z = _ZHolderWrongShapeThenValid(bad, good)

    arr = cmp._get_z(_E())
    assert arr is not None
    assert arr.shape == (2, 2, 2)


# ---------------------------------------------------------------------------
# _as_sites_iter — Sites container (``as_list``) path
# ---------------------------------------------------------------------------


def test_as_sites_iter_uses_as_list_for_sites_container(tmp_path, simulated_edi):
    p1 = tmp_path / "N01.edi"
    p2 = tmp_path / "N02.edi"
    p1.write_text(simulated_edi.read_text(encoding="utf-8"), encoding="utf-8")
    p2.write_text(simulated_edi.read_text(encoding="utf-8"), encoding="utf-8")
    sites = Sites([_load_edi(p1), _load_edi(p2)])

    pairs = list(cmp._as_sites_iter(sites))
    assert len(pairs) == 2
    assert all(isinstance(name, str) for name, _ in pairs)


def test_strike_estimate_accepts_sites_container(tmp_path, simulated_edi):
    p1 = tmp_path / "M01.edi"
    p2 = tmp_path / "M02.edi"
    p1.write_text(simulated_edi.read_text(encoding="utf-8"), encoding="utf-8")
    p2.write_text(simulated_edi.read_text(encoding="utf-8"), encoding="utf-8")
    e1 = ed.fill_missing(
        _load_edi(p1), how="zero", components=("Z",), inplace=False
    )
    e2 = ed.fill_missing(
        _load_edi(p2), how="zero", components=("Z",), inplace=False
    )
    sites = Sites([e1, e2])
    df = cmp.strike_estimate(sites, method="swift", api=False)
    assert len(df) == 2


# ---------------------------------------------------------------------------
# _tip_arr — tipper-holder discovery, exceptions, and Z-fallback shapes
# ---------------------------------------------------------------------------


def test_tip_arr_none_when_nothing_present():
    assert cmp._tip_arr(_EdiNoZ()) is None


def test_tip_arr_recovers_after_attribute_access_raises_then_uses_tipper():
    class _TipObj:
        @property
        def tipper(self):
            raise RuntimeError("boom")

        _tipper = np.zeros((3, 2), complex)

    class _E:
        Tip = _TipObj()

    arr = cmp._tip_arr(_E())
    assert arr is not None
    assert arr.shape == (3, 2)


def test_tip_arr_widens_1d_tipper():
    class _TipObj:
        tipper = np.zeros(4, complex)

    class _E:
        T = _TipObj()

    arr = cmp._tip_arr(_E())
    assert arr.shape == (4, 1)


def test_tip_arr_widens_3d_tipper_with_singleton_axis():
    class _TipObj:
        tipper = np.zeros((5, 1, 2), complex)

    class _E:
        TIP = _TipObj()

    arr = cmp._tip_arr(_E())
    assert arr.shape == (5, 2)


@pytest.mark.parametrize("attr_name", ["tipper", "tip", "_tipper"])
def test_tip_arr_falls_back_to_z_tipper_2d(attr_name):
    class _Z:
        pass

    z = _Z()
    setattr(z, attr_name, np.zeros((3, 2), complex))

    class _E:
        Z = z

    arr = cmp._tip_arr(_E())
    assert arr.shape == (3, 2)


def test_tip_arr_falls_back_to_z_tipper_3d_singleton():
    class _Z:
        tipper = np.zeros((2, 1, 2), complex)

    class _E:
        Z = _Z()

    arr = cmp._tip_arr(_E())
    assert arr.shape == (2, 2)


def test_tip_arr_falls_back_to_z_tipper_1d():
    class _Z:
        tipper = np.zeros(6, complex)

    class _E:
        Z = _Z()

    arr = cmp._tip_arr(_E())
    assert arr.shape == (6, 1)


def test_tip_arr_none_when_z_present_but_no_tipper_like_attrs():
    class _Z:
        pass

    class _E:
        Z = _Z()

    assert cmp._tip_arr(_E()) is None
