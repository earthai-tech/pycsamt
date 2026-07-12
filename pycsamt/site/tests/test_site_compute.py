# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from pycsamt.api import APIFrame, reset_api_view
from pycsamt.seg.edi import EDIFile
from pycsamt.site import compute as cmp
from pycsamt.site import edit as ed


def _load_edi(p: Path) -> EDIFile:
    return EDIFile(p)  # type: ignore


def _dup_edi(tmp_path: Path, src: Path, stem: str) -> Path:
    dst = tmp_path / f"{stem}.edi"
    dst.write_text(src.read_text(encoding="utf-8"),
                   encoding="utf-8")
    return dst


def _mk_two_edifiles(
    tmp_path: Path,
    simulated_edi: Path,
    n1: str = "C01",
    n2: str = "C02",
) -> tuple[EDIFile, EDIFile]:
    p1 = _dup_edi(tmp_path, simulated_edi, n1)
    p2 = _dup_edi(tmp_path, simulated_edi, n2)
    return _load_edi(p1), _load_edi(p2)


def test_strike_estimate_single_and_multi(
    tmp_path: Path, simulated_edi: Path
) -> None:
    ed1 = _load_edi(simulated_edi)
    # Ensure Z has the right shape for deterministic result
    ed1 = ed.fill_missing(ed1, how="zero",
                          components=("Z",), inplace=False)

    ang = cmp.strike_estimate(ed1, method="swift")
    assert isinstance(ang, float)
    assert 0.0 <= ang < 180.0
    # or: assert np.isfinite(ang)

    # For a zero tensor, first minimum is at 0 deg
    # Force Z to exact zeros so the minimum is at 0 deg
    Z = getattr(ed1, "Z", None)
    if Z is not None:
        f = cmp.get_freq(ed1)  # or use the shared get_freq
        n = len(f) if f is not None else 0
        if n > 0:
            Z.z = np.zeros((n, 2, 2), complex)
    ang = cmp.strike_estimate(ed1, method="swift")
    assert ang == pytest.approx(0.0, abs=1e-6)

    e1, e2 = _mk_two_edifiles(tmp_path, simulated_edi, "S01",
                              "S02")
    e1 = ed.fill_missing(e1, how="zero",
                         components=("Z",), inplace=False)
    e2 = ed.fill_missing(e2, how="zero",
                         components=("Z",), inplace=False)

    df = cmp.strike_estimate([e1, e2], method="swift", api=False)
    assert isinstance(df, pd.DataFrame)
    assert list(df.columns) == ["station", "method",
                                "theta_deg"]
    assert len(df) == 2
    assert np.all(np.isfinite(df["theta_deg"].values))


def test_compute_dataframe_results_can_return_api_frame(
    tmp_path: Path, simulated_edi: Path
) -> None:
    reset_api_view()
    e1, e2 = _mk_two_edifiles(tmp_path, simulated_edi, "A11",
                              "A12")
    e1 = ed.fill_missing(e1, how="zero",
                         components=("Z",), inplace=False)
    e2 = ed.fill_missing(e2, how="zero",
                         components=("Z",), inplace=False)
    sites = [e1, e2]

    strike = cmp.strike_estimate(sites, api=True)
    rho = cmp.res_at_freq(sites, 150.0, api=True)
    slope = cmp.phase_slope(sites, (90.0, 210.0), api=True)
    tipper = cmp.tipper_magnitude(sites, api=True)

    assert isinstance(strike, APIFrame)
    assert strike.kind == "site.compute.strike"
    assert isinstance(rho, APIFrame)
    assert rho.kind == "site.compute.resistivity"
    assert isinstance(slope, APIFrame)
    assert slope.kind == "site.compute.phase_slope"
    assert isinstance(tipper, APIFrame)
    assert tipper.kind == "site.compute.tipper"


def test_res_at_freq_nearest_and_interp(
    tmp_path: Path, simulated_edi: Path
) -> None:
    ed1 = _load_edi(simulated_edi)
    # Fill Z so Zxy, Zyx exist and are zeros
    ed1 = ed.fill_missing(ed1, how="zero",
                          components=("Z",), inplace=False)

    # Nearest mode: f_used must be one from the header
    out_near = cmp.res_at_freq(ed1, 150.0, how="nearest")
    assert set(out_near.keys()) == {"res_xy", "res_yx",
                                    "f_used"}
    assert out_near["f_used"] in {100.0, 200.0}
    # Off-diagonals are 0 -> apparent resistivities are 0
    assert out_near["res_xy"] == pytest.approx(0.0)
    assert out_near["res_yx"] == pytest.approx(0.0)

    # Interp mode: f_used equals query freq
    out_int = cmp.res_at_freq(ed1, 150.0, how="interp")
    assert out_int["f_used"] == pytest.approx(150.0)
    assert out_int["res_xy"] == pytest.approx(0.0)
    assert out_int["res_yx"] == pytest.approx(0.0)

    # Multi-site path returns a DataFrame
    e1, e2 = _mk_two_edifiles(tmp_path, simulated_edi, "R01",
                              "R02")
    e1 = ed.fill_missing(e1, how="zero",
                         components=("Z",), inplace=False)
    e2 = ed.fill_missing(e2, how="zero",
                         components=("Z",), inplace=False)

    df = cmp.res_at_freq([e1, e2], 150.0, how="interp", api=False)
    assert isinstance(df, pd.DataFrame)
    assert list(df.columns) == ["station", "res_xy", "res_yx",
                                "f_used"]
    assert len(df) == 2
    assert np.allclose(df["f_used"].values,
                       np.array([150.0, 150.0]))


def test_phase_slope_band(tmp_path: Path,
                          simulated_edi: Path) -> None:
    ed1 = _load_edi(simulated_edi)
    # Zero Z => phase is 0 deg over the band -> slope 0
    ed1 = ed.fill_missing(ed1, how="zero",
                          components=("Z",), inplace=False)

    out = cmp.phase_slope(ed1, band=(90.0, 210.0))
    assert set(out.keys()) == {"slope_xy", "slope_yx"}
    assert out["slope_xy"] == pytest.approx(0.0, abs=1e-6)
    assert out["slope_yx"] == pytest.approx(0.0, abs=1e-6)

    # Multi-site variant
    e1, e2 = _mk_two_edifiles(tmp_path, simulated_edi, "P01",
                              "P02")
    e1 = ed.fill_missing(e1, how="zero",
                         components=("Z",), inplace=False)
    e2 = ed.fill_missing(e2, how="zero",
                         components=("Z",), inplace=False)
    df = cmp.phase_slope([e1, e2], band=(90.0, 210.0), api=False)
    assert isinstance(df, pd.DataFrame)
    assert list(df.columns) == ["station", "slope_xy",
                                "slope_yx"]
    assert len(df) == 2
    assert np.allclose(df["slope_xy"].values,
                       np.zeros(2), atol=1e-6)
    assert np.allclose(df["slope_yx"].values,
                       np.zeros(2), atol=1e-6)


def test_tipper_magnitude_presence_and_stats(
    tmp_path: Path, simulated_edi: Path
) -> None:
    # No tipper present -> NaNs for summary
    ed1 = _load_edi(simulated_edi)
    out = cmp.tipper_magnitude(ed1, per_freq=False)
    assert set(out.keys()) == {"mean", "median", "max"}
    assert math.isnan(out["mean"])
    assert math.isnan(out["median"])
    assert math.isnan(out["max"])

    # Allocate an all-zero tipper -> stats become zeros
    e1, e2 = _mk_two_edifiles(tmp_path, simulated_edi, "T01",
                              "T02")
    e1 = ed.fill_missing(e1, how="zero",
                         components=("Tip",), inplace=False)
    e2 = ed.fill_missing(e2, how="zero",
                         components=("Tip",), inplace=False)

    # (n_freq, 1, 2) zero-filled tippers are accepted since the
    # tipper-shape widening in site.compute._tip_arr
    s = cmp.tipper_magnitude(e1, per_freq=False)
    assert s["mean"] == 0.0
    assert s["median"] == 0.0
    assert s["max"] == 0.0

    # Create a tipper holder and attach zero arrays
    f = cmp.get_freq(e1)
    n = len(f) if f is not None else 0
    class _Tip: pass
    tp = _Tip()
    tp._tipper = np.zeros((n, 2), complex)
    e1.Tip = tp

    # Now fill_missing can operate on the existing tipper
    e1 = ed.fill_missing(e1, how="zero", components=("Tip",),
                         inplace=True)

    s = cmp.tipper_magnitude(e1, per_freq=False)
    assert s["mean"] == pytest.approx(0.0)
    assert s["median"] == pytest.approx(0.0)
    assert s["max"] == pytest.approx(0.0)
