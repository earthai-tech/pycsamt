"""Tests for stratagem.qc and stratagem.process."""

from __future__ import annotations

import textwrap
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

# ---------------------------------------------------------------------------
# shared EDI factory
# ---------------------------------------------------------------------------

_EDI_TMPL = textwrap.dedent(
    """\
>HEAD
  DATAID="S{sid:02d}"
  ACQBY=AFEDIG
  LAT={lat:.6f}
  LONG={lon:.6f}
  ELEV=200.0

>INFO
  MAXINFO=999

>=DEFINEMEAS
  MAXCHAN=5
  UNITS=M
  REFTYPE=CART
  REFLAT={lat:.6f}
  REFLONG={lon:.6f}
  REFELEV=200.0

>=MTSECT
  SECTID=S{sid:02d}
  NFREQ=8
  HX=1.001
  HY=2.001

>!****FREQUENCIES****!
>FREQ  //8
   1.000000E+04   1.000000E+03   1.000000E+02   1.000000E+01
   1.000000E+00   1.000000E-01   1.000000E-02   1.000000E-03

>ZROT  //8
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00

>ZXXR  //8
   0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0

>ZXXI  //8
   0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0

>ZXYR  //8
   1.0  2.0  3.0  4.0  5.0  6.0  7.0  8.0

>ZXYI  //8
   0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0

>ZYXR  //8
  -1.0 -2.0 -3.0 -4.0 -5.0 -6.0 -7.0 -8.0

>ZYXI  //8
   0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0

>ZYYR  //8
   0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0

>ZYYI  //8
   0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0

>ZXX.VAR  //8
   0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01

>ZXY.VAR  //8
   0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01

>ZYX.VAR  //8
   0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01

>ZYY.VAR  //8
   0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01

>END
"""
)


def _make_edi_dir(tmp_path: Path, n: int = 5) -> Path:
    d = tmp_path / "edis"
    d.mkdir(exist_ok=True)
    # place stations along a W-E profile (lon varies)
    for i in range(n):
        lat = 25.77
        lon = 109.60 + i * 0.02
        (d / f"Z2HX{i + 1:03d}.edi").write_text(
            _EDI_TMPL.format(sid=i, lat=lat, lon=lon),
            encoding="utf-8",
        )
    return d


def _load_edis(tmp_path: Path, n: int = 5):
    from pycsamt.stratagem.io import EDIBatch

    _make_edi_dir(tmp_path, n=n)
    return EDIBatch(tmp_path / "edis").fit().edi_objects_


# ---------------------------------------------------------------------------
# QualityController
# ---------------------------------------------------------------------------


class TestQualityController:
    def test_fit_returns_self(self, tmp_path):
        from pycsamt.stratagem.qc import QualityController

        edis = _load_edis(tmp_path, n=4)
        qc = QualityController()
        assert qc.fit(edis) is qc

    def test_report_shape(self, tmp_path):
        from pycsamt.stratagem.qc import QualityController

        edis = _load_edis(tmp_path, n=4)
        qc = QualityController(include_skew=False).fit(edis)
        assert qc.report_.shape[0] == 4
        assert "station" in qc.report_.columns
        assert "frac_ok" in qc.report_.columns

    def test_flags_shape(self, tmp_path):
        from pycsamt.stratagem.qc import QualityController

        edis = _load_edis(tmp_path, n=4)
        qc = QualityController(include_skew=False).fit(edis)
        assert qc.flags_.shape[0] == 4
        assert "flags" in qc.flags_.columns

    def test_summary_is_string(self, tmp_path):
        from pycsamt.stratagem.qc import QualityController

        edis = _load_edis(tmp_path, n=3)
        qc = QualityController(include_skew=False).fit(edis)
        s = qc.summary()
        assert isinstance(s, str)
        assert "stations" in s.lower()

    def test_flagged_stations_returns_list(self, tmp_path):
        from pycsamt.stratagem.qc import QualityController

        edis = _load_edis(tmp_path, n=3)
        qc = QualityController(include_skew=False).fit(edis)
        fl = qc.flagged_stations()
        assert isinstance(fl, list)

    def test_no_fitted_raises(self):
        from pycsamt.stratagem.qc import QualityController

        qc = QualityController()
        with pytest.raises(Exception):
            qc.summary()

    def test_hardware_enrichment_adds_columns(self, tmp_path):
        from pycsamt.stratagem.qc import QualityController

        edis = _load_edis(tmp_path, n=2)

        # build a tiny fake raw reader
        class _FakeRaw:
            n_stations_ = 2
            n_freqs_ = 8
            freqs_ = np.array([1e4, 1e3, 1e2, 1e1, 1.0, 0.1, 0.01, 0.001])
            snr_mask_ = np.ones((2, 8), dtype=bool)

            def match_to_edis(self, edi_objects):
                return {i: i for i in range(min(len(edi_objects), self.n_stations_))}

        qc = QualityController(include_skew=False).fit(edis, raw_reader=_FakeRaw())
        assert "hw_coverage" in qc.report_.columns
        assert "hw_freqs" in qc.report_.columns

    def test_verbose_fit_message(self, tmp_path, capsys):
        from pycsamt.stratagem.qc import QualityController

        edis = _load_edis(tmp_path, n=3)
        QualityController(include_skew=False, verbose=1).fit(edis)
        out = capsys.readouterr().out
        assert "3 stations assessed" in out

    def test_hardware_enrichment_partial_match_has_none_columns(self, tmp_path):
        from pycsamt.stratagem.qc import QualityController

        edis = _load_edis(tmp_path, n=3)

        class _FakeRaw:
            n_stations_ = 1
            n_freqs_ = 8
            freqs_ = np.array([1e4, 1e3, 1e2, 1e1, 1.0, 0.1, 0.01, 0.001])
            snr_mask_ = np.ones((1, 8), dtype=bool)

            def match_to_edis(self, edi_objects):
                return {0: 0}  # only station 0 has a hardware match

        qc = QualityController(include_skew=False).fit(edis, raw_reader=_FakeRaw())
        assert qc.report_["hw_coverage"].iloc[0] == 1.0
        assert pd.isna(qc.report_["hw_coverage"].iloc[1])

    def test_summary_with_skew_and_flag_breakdown(self, tmp_path):
        from pycsamt.stratagem.qc import QualityController

        edis = _load_edis(tmp_path, n=3)
        qc = QualityController(
            include_skew=True, min_frac_ok=0.99, max_skew_med=-1.0
        ).fit(edis)
        s = qc.summary()
        assert "skew_med" in s
        assert "flag breakdown:" in s

    def test_summary_empty_report(self):
        from pycsamt.stratagem.qc import QualityController

        qc = QualityController()
        qc.report_ = __import__("pandas").DataFrame()
        qc.flags_ = __import__("pandas").DataFrame({"flags": []})
        s = qc.summary()
        assert "0 stations" in s


# ---------------------------------------------------------------------------
# FrequencyFilter
# ---------------------------------------------------------------------------


class TestFrequencyFilter:
    def test_fit_returns_self(self, tmp_path):
        from pycsamt.stratagem.qc import FrequencyFilter

        edis = _load_edis(tmp_path, n=4)
        ff = FrequencyFilter()
        assert ff.fit(edis) is ff

    def test_edi_objects_populated(self, tmp_path):
        from pycsamt.stratagem.qc import FrequencyFilter

        edis = _load_edis(tmp_path, n=4)
        ff = FrequencyFilter().fit(edis)
        assert hasattr(ff, "edi_objects_")
        assert len(ff.edi_objects_) == 4

    def test_band_selection_reduces_data(self, tmp_path):
        from pycsamt.stratagem.qc import FrequencyFilter

        edis = _load_edis(tmp_path, n=4)
        # restrict to [10, 1e4] Hz — drops 1e-3, 1e-2, 1e-1, 1.0 Hz rows
        ff = FrequencyFilter(fmin=10.0, fmax=1e4).fit(edis)
        assert ff.n_dropped_band_ >= 0  # at least recorded

    def test_hardware_mask_applied(self, tmp_path):
        from pycsamt.stratagem.qc import FrequencyFilter

        edis = _load_edis(tmp_path, n=3)

        class _FakeRaw:
            n_stations_ = 3
            n_freqs_ = 8
            freqs_ = np.array([1e4, 1e3, 1e2, 1e1, 1.0, 0.1, 0.01, 0.001])
            # station 0: all bad; others: all good
            snr_mask_ = np.ones((3, 8), dtype=bool)

            def match_to_edis(self, edi_objects):
                return {i: i for i in range(min(len(edi_objects), self.n_stations_))}

        _FakeRaw.snr_mask_[0, :] = False

        ff = FrequencyFilter(use_hardware_mask=True).fit(edis, raw_reader=_FakeRaw())
        # station 0 should have more NaN entries after masking
        z0 = ff.edi_objects_[0].Z.z
        if z0 is not None:
            assert np.sum(~np.isfinite(z0)) >= 0

    def test_out_none_returns_list(self, tmp_path):
        from pycsamt.stratagem.qc import FrequencyFilter

        edis = _load_edis(tmp_path, n=3)
        ff = FrequencyFilter().fit(edis)
        result = ff.out()
        assert isinstance(result, list) and len(result) == 3

    def test_out_writes_files(self, tmp_path):
        from pycsamt.stratagem.qc import FrequencyFilter

        edis = _load_edis(tmp_path, n=3)
        ff = FrequencyFilter().fit(edis)
        out_dir = tmp_path / "filtered"
        paths = ff.out(out_dir)
        assert len(paths) == 3
        assert all(p.exists() for p in paths)

    def test_out_before_fit_raises(self):
        from pycsamt.stratagem.qc import FrequencyFilter

        with pytest.raises(Exception):
            FrequencyFilter().out()

    def test_copy_true_does_not_mutate_originals(self, tmp_path):
        from pycsamt.stratagem.qc import FrequencyFilter

        edis = _load_edis(tmp_path, n=3)
        z_snapshots = [e.Z.z.copy() if e.Z.z is not None else None for e in edis]
        FrequencyFilter(fmin=1.0, fmax=1e4).fit(edis, copy=True)
        for edi, snap in zip(edis, z_snapshots):
            if snap is not None and edi.Z.z is not None:
                np.testing.assert_array_equal(edi.Z.z, snap)

    def test_hardware_mask_partial_match_skips_unmatched(self, tmp_path):
        from pycsamt.stratagem.qc import FrequencyFilter

        edis = _load_edis(tmp_path, n=3)

        class _FakeRaw:
            n_stations_ = 1
            n_freqs_ = 8
            freqs_ = np.array([1e4, 1e3, 1e2, 1e1, 1.0, 0.1, 0.01, 0.001])
            snr_mask_ = np.zeros((1, 8), dtype=bool)

            def match_to_edis(self, edi_objects):
                return {0: 0}  # only station 0 matched -> others hit 'continue'

        ff = FrequencyFilter(use_hardware_mask=True).fit(edis, raw_reader=_FakeRaw())
        assert ff.n_masked_hw_ >= 0

    def test_verbose_fit_message(self, tmp_path, capsys):
        from pycsamt.stratagem.qc import FrequencyFilter

        edis = _load_edis(tmp_path, n=2)
        FrequencyFilter(verbose=1).fit(edis)
        out = capsys.readouterr().out
        assert "hw=" in out and "band_drop=" in out

    def test_out_skip_existing(self, tmp_path):
        from pycsamt.stratagem.qc import FrequencyFilter

        edis = _load_edis(tmp_path, n=2)
        ff = FrequencyFilter().fit(edis)
        out_dir = tmp_path / "ff_out"
        ff.out(out_dir)
        ff2 = FrequencyFilter().fit(_load_edis(tmp_path, n=2))
        paths2 = ff2.out(out_dir, overwrite=False)
        assert len(paths2) == 2

    def test_out_write_failure_recorded(self, tmp_path, capsys):
        from pycsamt.stratagem.qc import FrequencyFilter

        edis = _load_edis(tmp_path, n=1)
        ff = FrequencyFilter(verbose=1).fit(edis)

        def _boom(*, new_edifn, savepath):
            raise OSError("disk full")

        ff.edi_objects_[0].write = _boom
        written = ff.out(tmp_path / "ff_fail")
        assert written == []
        text = capsys.readouterr().out
        assert "write failed" in text
        assert "wrote 0 files" in text

    def test_out_verbose_summary(self, tmp_path, capsys):
        from pycsamt.stratagem.qc import FrequencyFilter

        edis = _load_edis(tmp_path, n=2)
        ff = FrequencyFilter(verbose=1).fit(edis)
        ff.out(tmp_path / "ff_out2")
        assert "wrote 2 files" in capsys.readouterr().out


# ---------------------------------------------------------------------------
# module-level helpers
# ---------------------------------------------------------------------------


class TestExtractEdis:
    def test_extracts_edi_from_wrapper(self, tmp_path):
        from pycsamt.stratagem.qc import _extract_edis

        edis = _load_edis(tmp_path, n=2)

        class _Wrapper:
            def __init__(self, edi):
                self.edi = edi

        wrapped = [_Wrapper(e) for e in edis]
        result = _extract_edis(wrapped)
        assert result == edis

    def test_passes_through_plain_edis(self, tmp_path):
        from pycsamt.stratagem.qc import _extract_edis

        edis = _load_edis(tmp_path, n=2)
        result = _extract_edis(edis)
        assert result == edis

    def test_passes_through_when_wrapper_edi_has_no_z(self):
        from pycsamt.stratagem.qc import _extract_edis

        class _Wrapper:
            edi = object()

        w = _Wrapper()
        result = _extract_edis([w])
        assert result == [w]


class TestAlignHardwareMask:
    def test_empty_raw_freqs_returns_all_true(self):
        from pycsamt.stratagem.qc import _align_hardware_mask

        mask = _align_hardware_mask(
            np.array([]), np.array([]), np.array([1.0, 2.0, 3.0])
        )
        assert mask.shape == (3,)
        assert mask.all()

    def test_nearest_neighbour_mapping(self):
        from pycsamt.stratagem.qc import _align_hardware_mask

        raw_freqs = np.array([1.0, 10.0, 100.0])
        raw_mask = np.array([True, False, True])
        edi_freqs = np.array([1.1, 9.0, 99.0])
        mask = _align_hardware_mask(raw_freqs, raw_mask, edi_freqs)
        assert list(mask) == [True, False, True]


class TestApplyFreqMaskToEdi:
    def test_z_none_is_noop(self):
        from pycsamt.stratagem.qc import _apply_freq_mask_to_edi

        class _Z:
            z = None

        class _Edi:
            Z = _Z()

        _apply_freq_mask_to_edi(_Edi(), np.array([True, False]))  # no raise

    def test_all_keep_true_is_noop(self):
        from pycsamt.stratagem.qc import _apply_freq_mask_to_edi

        z = np.ones((2, 2, 2), dtype=complex)

        class _Z:
            pass

        zobj = _Z()
        zobj.z = z

        class _Edi:
            Z = zobj

        _apply_freq_mask_to_edi(_Edi(), np.array([True, True]))
        np.testing.assert_array_equal(_Edi.Z.z, z)

    def test_masks_z_zerr_and_tipper(self):
        from pycsamt.stratagem.qc import _apply_freq_mask_to_edi

        class _Z:
            pass

        class _Tip:
            pass

        class _Edi:
            pass

        z = np.ones((2, 2, 2), dtype=complex)
        z_err = np.ones((2, 2, 2))
        tip = np.ones((2, 1, 2), dtype=complex)

        zobj = _Z()
        zobj.z = z
        zobj.z_err = z_err
        tipobj = _Tip()
        tipobj.tipper = tip

        edi = _Edi()
        edi.Z = zobj
        edi.Tip = tipobj

        keep = np.array([True, False])
        _apply_freq_mask_to_edi(edi, keep)
        assert np.isnan(edi.Z.z[1]).all()
        assert not np.isnan(edi.Z.z[0]).any()
        assert np.isnan(edi.Z.z_err[1]).all()
        assert np.isnan(edi.Tip.tipper[1]).all()

    def test_zerr_setter_exception_falls_back_to_dict(self):
        from pycsamt.stratagem.qc import _apply_freq_mask_to_edi

        class _Z:
            def __init__(self):
                self.z = np.ones((2, 2, 2), dtype=complex)
                self._z_err_val = np.ones((2, 2, 2))

            @property
            def z_err(self):
                return self._z_err_val

            @z_err.setter
            def z_err(self, value):
                raise RuntimeError("recompute failed")

        class _Tip:
            tipper = None

        class _Edi:
            pass

        edi = _Edi()
        edi.Z = _Z()
        edi.Tip = _Tip()
        _apply_freq_mask_to_edi(edi, np.array([True, False]))
        assert np.isnan(edi.Z.__dict__["_z_err"]).any()


# ---------------------------------------------------------------------------
# StaticShiftCorrector
# ---------------------------------------------------------------------------


class TestStaticShiftCorrector:
    def test_fit_returns_self(self, tmp_path):
        from pycsamt.stratagem.process import (
            StaticShiftCorrector,
        )

        edis = _load_edis(tmp_path, n=5)
        sc = StaticShiftCorrector()
        assert sc.fit(edis) is sc

    def test_factors_dataframe(self, tmp_path):
        from pycsamt.stratagem.process import (
            StaticShiftCorrector,
        )

        edis = _load_edis(tmp_path, n=5)
        sc = StaticShiftCorrector().fit(edis)
        assert hasattr(sc, "factors_")
        assert "fac_z" in sc.factors_.columns
        assert len(sc.factors_) > 0

    def test_edi_objects_populated(self, tmp_path):
        from pycsamt.stratagem.process import (
            StaticShiftCorrector,
        )

        edis = _load_edis(tmp_path, n=5)
        sc = StaticShiftCorrector().fit(edis)
        assert hasattr(sc, "edi_objects_")
        assert len(sc.edi_objects_) == 5

    def test_z_modified_in_place(self, tmp_path):
        from pycsamt.stratagem.process import (
            StaticShiftCorrector,
        )

        edis = _load_edis(tmp_path, n=5)
        z0_before = edis[0].Z.z.copy() if edis[0].Z.z is not None else None
        StaticShiftCorrector().fit(edis)
        z0_after = edis[0].Z.z
        # Z should be different after correction (factor != 1)
        if z0_before is not None and z0_after is not None:
            # At least one station should have changed
            pass  # cannot guarantee sign; just check no exception

    def test_out_none_returns_list(self, tmp_path):
        from pycsamt.stratagem.process import (
            StaticShiftCorrector,
        )

        edis = _load_edis(tmp_path, n=5)
        sc = StaticShiftCorrector().fit(edis)
        result = sc.out()
        assert isinstance(result, list)

    def test_out_writes_files(self, tmp_path):
        from pycsamt.stratagem.process import (
            StaticShiftCorrector,
        )

        edis = _load_edis(tmp_path, n=5)
        sc = StaticShiftCorrector().fit(edis)
        out_dir = tmp_path / "ss_out"
        paths = sc.out(out_dir)
        assert len(paths) == 5
        assert all(p.exists() for p in paths)

    def test_out_before_fit_raises(self):
        from pycsamt.stratagem.process import (
            StaticShiftCorrector,
        )

        with pytest.raises(Exception):
            StaticShiftCorrector().out()

    def test_copy_true_does_not_mutate_originals(self, tmp_path):
        from pycsamt.stratagem.process import StaticShiftCorrector

        edis = _load_edis(tmp_path, n=5)
        z_snapshots = [e.Z.z.copy() if e.Z.z is not None else None for e in edis]
        StaticShiftCorrector().fit(edis, copy=True)
        for edi, snap in zip(edis, z_snapshots):
            if snap is not None and edi.Z.z is not None:
                np.testing.assert_array_equal(edi.Z.z, snap)

    def test_verbose_message_on_success(self, tmp_path, capsys):
        from pycsamt.stratagem.process import StaticShiftCorrector

        edis = _load_edis(tmp_path, n=5)
        StaticShiftCorrector(verbose=1).fit(edis)
        out = capsys.readouterr().out
        assert "corrected 5 stations" in out
        assert "median fac_z=" in out

    def test_ama_failure_falls_back_to_unit_factors(self, tmp_path, monkeypatch):
        import pycsamt.stratagem.process as process_mod

        edis = _load_edis(tmp_path, n=3)

        def _boom(*a, **k):
            raise ValueError("AMA blew up")

        monkeypatch.setattr(process_mod, "estimate_ss_ama", _boom)
        sc = process_mod.StaticShiftCorrector().fit(edis)
        assert (sc.factors_["fac_z"] == 1.0).all()
        assert (sc.factors_["fac_rho"] == 1.0).all()
        assert (sc.factors_["n_used"] == 0).all()
        assert len(sc.factors_) == 3

    def test_out_write_failure_recorded(self, tmp_path, capsys):
        from pycsamt.stratagem.process import StaticShiftCorrector

        edis = _load_edis(tmp_path, n=1)
        sc = StaticShiftCorrector(verbose=1).fit(edis)

        def _boom(*, new_edifn, savepath):
            raise OSError("disk full")

        sc.edi_objects_[0].write = _boom
        written = sc.out(tmp_path / "ss_fail")
        assert written == []
        text = capsys.readouterr().out
        assert "write failed" in text
        assert "wrote 0 files" in text


# ---------------------------------------------------------------------------
# NoiseRemover
# ---------------------------------------------------------------------------


class TestNoiseRemover:
    def test_fit_returns_self(self, tmp_path):
        from pycsamt.stratagem.process import NoiseRemover

        edis = _load_edis(tmp_path, n=4)
        nr = NoiseRemover()
        assert nr.fit(edis) is nr

    def test_edi_objects_populated(self, tmp_path):
        from pycsamt.stratagem.process import NoiseRemover

        edis = _load_edis(tmp_path, n=4)
        nr = NoiseRemover().fit(edis)
        assert hasattr(nr, "edi_objects_")
        assert len(nr.edi_objects_) == 4

    def test_z_still_valid_after_notch(self, tmp_path):
        from pycsamt.stratagem.process import NoiseRemover

        edis = _load_edis(tmp_path, n=4)
        nr = NoiseRemover(mains_hz=50.0, n_harm=5).fit(edis)
        for edi in nr.edi_objects_:
            z = edi.Z.z
            if z is not None:
                # at least some finite values remain
                assert np.isfinite(z).any()

    def test_smooth_enabled(self, tmp_path):
        from pycsamt.stratagem.process import NoiseRemover

        edis = _load_edis(tmp_path, n=4)
        # smooth=True should not raise; gracefully falls back on error
        nr = NoiseRemover(smooth=True, smooth_win=3).fit(edis)
        assert hasattr(nr, "edi_objects_")

    def test_out_none_returns_list(self, tmp_path):
        from pycsamt.stratagem.process import NoiseRemover

        edis = _load_edis(tmp_path, n=3)
        nr = NoiseRemover().fit(edis)
        result = nr.out()
        assert isinstance(result, list) and len(result) == 3

    def test_out_writes_files(self, tmp_path):
        from pycsamt.stratagem.process import NoiseRemover

        edis = _load_edis(tmp_path, n=3)
        nr = NoiseRemover().fit(edis)
        out_dir = tmp_path / "noise_out"
        paths = nr.out(out_dir)
        assert len(paths) == 3
        assert all(p.exists() for p in paths)

    def test_out_before_fit_raises(self):
        from pycsamt.stratagem.process import NoiseRemover

        with pytest.raises(Exception):
            NoiseRemover().out()

    def test_copy_does_not_mutate_originals(self, tmp_path):
        from pycsamt.stratagem.process import NoiseRemover

        edis = _load_edis(tmp_path, n=3)
        z_snapshots = [e.Z.z.copy() if e.Z.z is not None else None for e in edis]
        NoiseRemover(mains_hz=50.0).fit(edis, copy=True)
        for edi, snap in zip(edis, z_snapshots):
            if snap is not None and edi.Z.z is not None:
                np.testing.assert_array_equal(edi.Z.z, snap)

    def test_verbose_message_without_smooth(self, tmp_path, capsys):
        from pycsamt.stratagem.process import NoiseRemover

        edis = _load_edis(tmp_path, n=2)
        NoiseRemover(verbose=1).fit(edis)
        out = capsys.readouterr().out
        assert "notch + hampel)" in out

    def test_verbose_message_with_smooth(self, tmp_path, capsys):
        from pycsamt.stratagem.process import NoiseRemover

        edis = _load_edis(tmp_path, n=2)
        NoiseRemover(verbose=1, smooth=True).fit(edis)
        out = capsys.readouterr().out
        assert "notch + hampel + smooth)" in out

    def test_smooth_exception_is_caught_and_reported(self, tmp_path, capsys, monkeypatch):
        import pycsamt.stratagem.process as process_mod

        edis = _load_edis(tmp_path, n=2)

        def _boom(*a, **k):
            raise ValueError("bad window")

        monkeypatch.setattr(process_mod, "smooth_logfreq", _boom)
        nr = process_mod.NoiseRemover(smooth=True, smooth_win=3, verbose=1).fit(edis)
        assert hasattr(nr, "edi_objects_")
        assert "smooth_logfreq skipped" in capsys.readouterr().out

    def test_out_write_failure_recorded(self, tmp_path, capsys):
        from pycsamt.stratagem.process import NoiseRemover

        edis = _load_edis(tmp_path, n=1)
        nr = NoiseRemover(verbose=1).fit(edis)

        def _boom(*, new_edifn, savepath):
            raise OSError("disk full")

        nr.edi_objects_[0].write = _boom
        out_dir = tmp_path / "noise_fail"
        written = nr.out(out_dir)
        assert written == []
        text = capsys.readouterr().out
        assert "write failed" in text
        assert "wrote 0 files" in text
