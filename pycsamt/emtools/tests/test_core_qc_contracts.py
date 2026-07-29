from __future__ import annotations

from types import SimpleNamespace

import matplotlib
import numpy as np
import pandas as pd
import pytest

matplotlib.use("Agg")

from pycsamt.emtools import _core, qc, ss


class _FakeZ:
    def __init__(self, z, freq, z_err=None):
        self.z = np.asarray(z, dtype=np.complex128)
        self.freq = np.asarray(freq, dtype=float)
        if z_err is not None:
            self.z_err = np.asarray(z_err, dtype=float)


class _FakeTipper:
    def __init__(self, tipper, freq):
        self.tipper = np.asarray(tipper, dtype=np.complex128)
        self.freq = np.asarray(freq, dtype=float)


class _FakeSite:
    def __init__(self, station, z, freq, *, z_err=None, tipper=None):
        self.station = station
        self.Z = _FakeZ(z, freq, z_err=z_err)
        self.freq = np.asarray(freq, dtype=float)
        if tipper is not None:
            self.Tipper = _FakeTipper(tipper, freq)


def _z_tensor(freq, *, rho=100.0, z_err_scale=0.05):
    z_abs = np.sqrt(5.0 * np.asarray(freq, dtype=float) * rho)
    z = np.zeros((len(freq), 2, 2), dtype=np.complex128)
    z[:, 0, 1] = z_abs * (1.0 + 1.0j) / np.sqrt(2.0)
    z[:, 1, 0] = -z_abs * (1.0 + 1.0j) / np.sqrt(2.0)
    z_err = np.abs(z) * z_err_scale
    return z, z_err


def _site(station="S00", *, n_freq=4, rho=100.0, z_err_scale=0.05):
    freq = np.logspace(0.0, 3.0, n_freq)
    z, z_err = _z_tensor(freq, rho=rho, z_err_scale=z_err_scale)
    tipper = np.zeros((n_freq, 2), dtype=np.complex128)
    return _FakeSite(station, z, freq, z_err=z_err, tipper=tipper)


def test_ensure_sites_rejects_none_before_loader():
    with pytest.raises(ValueError, match="got None"):
        _core.ensure_sites(None)


def test_ensure_sites_forwards_default_global_order(monkeypatch):
    class _OrderedSites:
        def ordered(self, by):
            self.by = by
            return self

        def __len__(self):
            return 1

    fake = _OrderedSites()
    monkeypatch.setattr("pycsamt.site.base.to_sites", lambda *args, **kwargs: fake)
    monkeypatch.setattr("pycsamt.site.base.Sites", _OrderedSites)

    out = _core.ensure_sites("ignored")

    assert out is fake
    assert out.by is None


def test_qc_table_has_stable_columns_and_skips_invalid_z(monkeypatch):
    valid = _site("S01")
    invalid = SimpleNamespace(station="BAD", Z=SimpleNamespace(z=np.ones((2, 3))))
    skew = pd.DataFrame(
        {
            "station": ["S01", "S01"],
            "beta": [2.0, -4.0],
        }
    )

    monkeypatch.setattr(qc, "ensure_sites", lambda sites, **_: sites)
    monkeypatch.setattr(qc, "build_phase_tensor_table", lambda *_, **__: skew)

    table = qc.build_qc_table([valid, invalid], include_skew=True)

    assert list(table.columns) == [
        "station",
        "n_freq",
        "n_ok",
        "frac_ok",
        "n_tip",
        "n_tip_ok",
        "snr_med",
        "pmin",
        "pmax",
        "skew_med",
        "skew_iqr",
    ]
    assert table["station"].tolist() == ["S01"]
    assert table.loc[0, "n_freq"] == 4
    assert table.loc[0, "n_ok"] == 4
    assert table.loc[0, "n_tip_ok"] == 4
    assert table.loc[0, "frac_ok"] == pytest.approx(1.0)
    assert table.loc[0, "skew_med"] == pytest.approx(3.0)


def test_qc_flags_report_low_coverage_and_low_snr(monkeypatch):
    site = _site("S02", z_err_scale=10.0)
    site.Z.z[1, 0, 1] = np.nan
    site.Z.z[1, 1, 0] = np.nan

    monkeypatch.setattr(qc, "ensure_sites", lambda sites, **_: sites)
    monkeypatch.setattr(
        qc,
        "build_phase_tensor_table",
        lambda *_, **__: pd.DataFrame({"station": ["S02"], "beta": [1.0]}),
    )

    flags = qc.qc_flags([site], min_frac_ok=0.9, min_snr_med=2.0)

    assert flags.loc[0, "station"] == "S02"
    assert set(flags.loc[0, "flags"].split(",")) == {
        "low_coverage",
        "low_snr",
    }


def test_plot_ss_summary_uses_synthetic_arrays():
    freqs = np.array([1.0, 10.0, 100.0])
    before = np.array([[2.0, 2.1, 2.2], [2.3, 2.4, 2.5]])
    after = before - 0.2

    fig = ss.plot_ss_summary(
        before,
        after,
        freqs=freqs,
        station_labels=["S01", "S02"],
    )

    assert len(fig.axes) >= 4
    assert fig.axes[0].get_title().lower().startswith("(a)")
    assert fig.axes[2].get_title().lower().startswith("(c)")


def test_plot_ss_1d_curves_falls_back_when_station_selection_is_unknown():
    freqs = np.array([1.0, 10.0, 100.0])
    before = np.array([[2.0, 2.1, 2.2], [2.3, 2.4, 2.5]])
    after = before - 0.2

    fig = ss.plot_ss_1d_curves(
        before,
        after,
        freqs=freqs,
        stations=["missing"],
        station_labels=["S01", "S02"],
        n_cols=1,
    )

    assert fig.axes[0].get_title() == "S01"
    assert {line.get_label() for line in fig.axes[0].lines} >= {
        "before",
        "after",
    }
