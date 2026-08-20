"""Extra coverage tests for pycsamt.emtools.tf (transfer-function / tipper).

This module targets the large swaths of ``tf.py`` left uncovered by
``test_emtools_tf.py`` -- primarily the map-view / section / polar /
rose / spectra-wrapper / multi-period functions that consume real
tipper data.

Fixture design note (historical)
---------------------------------
Two upstream bugs -- now fixed -- used to make it hard to reach the
real tipper/2-D-map code paths in this module through a normal
``Sites``-wrapped collection:

- ``pycsamt.emtools._core._get_t_block`` used to probe lowercase
  ``"tipper"`` on ``ed`` directly, which matched ``Site.tipper``
  (already a plain ``(n, 2)`` ndarray) before falling through to
  ``ed.edi.Tipper`` (the real wrapper). The fake ``Tipper`` fixture
  below used to store a **pre-transposed** ``(2, n)`` array to exploit
  an accidental ``ndarray.T`` access and cancel the effect back out.
  ``_get_t_block`` now probes only exact-cased wrapper names
  (``"Tipper"``/``"Tip"``) on ``ed``, mirroring ``_get_z_block``'s
  ``"Z"``-only first pass, so the fixture stores tipper data in its
  natural ``(n, 2)`` orientation again.
- ``_station_xy`` (defined in ``tf.py`` itself) used to check
  ``getattr(ed, "east", None)`` etc. directly on ``ed`` with no
  ``ed.edi`` fallback, so genuine ``Site``-wrapped collections (which
  only expose geography via ``.coords`` -> ``(lat, lon, elev)``)
  always fell back to profile/index mode. ``_station_xy`` now also
  falls back to ``ed.coords`` (as x=lon, y=lat) when the flat
  east/north/x/y/lon/lat attributes aren't present. A couple of tests
  below still monkeypatch ``tf._station_xy`` directly where an exact,
  hand-picked 2-D layout is needed regardless of ``.coords``.
"""

from __future__ import annotations

import matplotlib
import numpy as np
import pytest

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from pycsamt.emtools import tf as tfm
from pycsamt.emtools.tf import (
    plot_induction_arrows,
    plot_induction_convention,
    plot_induction_map,
    plot_induction_map_from_spectra,
    plot_induction_multiperiod_map,
    plot_induction_rose,
    plot_induction_rose_from_spectra,
    plot_induction_section,
    plot_tipper_hodograms,
    plot_tipper_polar,
    plot_tipper_polar_from_spectra,
)

# ─────────────────────────────────────────────────────────────────────────────
# Shared fixtures
# ─────────────────────────────────────────────────────────────────────────────


class _FakeZ:
    def __init__(self, z, freq):
        self.z = np.asarray(z, dtype=complex)
        self.freq = np.asarray(freq, dtype=float)


class _FakeTipper:
    """Tipper sub-object, exposed via ``ed.edi.Tipper`` once Site-wrapped."""

    def __init__(self, tipper, freq):
        self.tipper = np.asarray(tipper, dtype=complex)  # (n, 2)
        self.freq = np.asarray(freq, dtype=float)


class _FakeSite:
    def __init__(self, station, z, freq, *, tipper=None, east=None, north=None):
        self.station = station
        self.Z = _FakeZ(z, freq)
        self.freq = np.asarray(freq, dtype=float)
        if tipper is not None:
            self.Tipper = _FakeTipper(tipper, freq)
        if east is not None:
            self.east = float(east)
            self.north = float(north)

    def get_section(self, *_, **__):
        return None


def _freqs(n: int = 12, f_lo: float = 0.1, f_hi: float = 1e4) -> np.ndarray:
    return np.logspace(np.log10(f_lo), np.log10(f_hi), n)


def _iso_z(freqs: np.ndarray, rho: float = 100.0) -> np.ndarray:
    amp = np.sqrt(5.0 * freqs * rho)
    z = np.zeros((freqs.size, 2, 2), dtype=complex)
    z[:, 0, 1] = amp * (1 + 1j) / np.sqrt(2)
    z[:, 1, 0] = -amp * (1 + 1j) / np.sqrt(2)
    return z


def _tipper_arr(
    freqs: np.ndarray, amp: float = 0.15, phase_span: float = np.pi / 4
) -> np.ndarray:
    """Synthetic, correctly-oriented (n, 2) complex tipper array."""
    t = np.zeros((freqs.size, 2), dtype=complex)
    t[:, 0] = amp * np.exp(1j * np.linspace(0, phase_span, freqs.size))
    t[:, 1] = amp * 0.5 * np.exp(1j * np.linspace(0, phase_span * 0.6, freqs.size))
    return t


def _site(
    name: str,
    n: int = 12,
    *,
    with_tipper: bool = True,
    east=None,
    north=None,
    amp: float = 0.15,
) -> _FakeSite:
    fr = _freqs(n)
    tip = _tipper_arr(fr, amp=amp) if with_tipper else None
    return _FakeSite(name, _iso_z(fr), fr, tipper=tip, east=east, north=north)


def _profile(n_sites: int = 4, *, n: int = 12, with_tipper: bool = True):
    return [
        _site(
            f"S{i:02d}",
            n=n,
            east=i * 200.0,
            north=0.0,
            with_tipper=with_tipper,
            amp=0.1 + 0.02 * i,
        )
        for i in range(n_sites)
    ]


# ── Sanity check that the workaround actually restores real tipper data ──


def test_fixture_reaches_real_tipper_path():
    """Guard test: confirms plot_tipper_hodograms sees real tipper data."""
    sites = [_site("S00", n=8)]
    fig = plot_tipper_hodograms(sites)
    plt.close("all")
    # With real tipper data the figure has two axes (Tx, Ty panels)
    # rather than the single "no tipper" placeholder axis.
    assert len(fig.axes) == 2


# ── Fake Spectra duck types for the *_from_spectra wrappers ──


class _FakeSpectraTip:
    def __init__(self, tipper):
        # shape (nf, 1, 2) mirroring the real Tipper container so that
        # ``tip.tipper[:, 0, :]`` yields (nf, 2).
        self.tipper = np.asarray(tipper, dtype=complex)[:, np.newaxis, :]
        self.freq = None  # freq is read off the returned tip object


class _FakeSpectraTipWithFreq(_FakeSpectraTip):
    def __init__(self, tipper, freq):
        super().__init__(tipper)
        self.freq = np.asarray(freq, dtype=float)


class _FakeSpectraObj:
    """Duck-typed ``Spectra`` stand-in exposing only ``to_Z``."""

    def __init__(self, tipper=None, freq=None, name="sp0"):
        self.name = name
        if tipper is None:
            self._tip = None
        else:
            self._tip = _FakeSpectraTipWithFreq(tipper, freq)

    def to_Z(self, *, estimate_error=False):
        return None, self._tip


def _spectra_dict(n_sites=3, n=10, with_tipper=True):
    fr = _freqs(n)
    out = {}
    for i in range(n_sites):
        tip = _tipper_arr(fr, amp=0.1 + 0.03 * i) if with_tipper else None
        out[f"SP{i:02d}"] = _FakeSpectraObj(tip, fr, name=f"SP{i:02d}")
    return out, fr


# ─────────────────────────────────────────────────────────────────────────────
# plot_induction_map
# ─────────────────────────────────────────────────────────────────────────────


class TestPlotInductionMap:
    def test_basic_profile(self):
        sites = _profile(4)
        ax = plot_induction_map(sites, period=1.0)
        plt.close("all")
        assert ax is not None

    def test_no_tipper_data(self):
        sites = _profile(3, with_tipper=False)
        ax = plot_induction_map(sites, period=1.0)
        plt.close("all")
        assert "no tipper" in ax.texts[0].get_text()

    def test_show_real_only(self):
        sites = _profile(3)
        ax = plot_induction_map(sites, period=1.0, show_imag=False)
        plt.close("all")
        assert ax is not None

    def test_show_imag_only(self):
        sites = _profile(3)
        ax = plot_induction_map(sites, period=1.0, show_real=False)
        plt.close("all")
        assert ax is not None

    def test_no_arrows_no_legend(self):
        sites = _profile(3)
        ax = plot_induction_map(sites, period=1.0, show_real=False, show_imag=False)
        plt.close("all")
        assert ax is not None

    def test_clim_and_no_colorbar(self):
        sites = _profile(3)
        ax = plot_induction_map(sites, period=1.0, clim=(0.0, 1.0), show_colorbar=False)
        plt.close("all")
        assert ax is not None

    def test_single_site_scale_fallback(self):
        sites = [_site("S00")]
        ax = plot_induction_map(sites, period=1.0)
        plt.close("all")
        assert ax is not None

    def test_station_labels_off(self):
        sites = _profile(6)
        ax = plot_induction_map(sites, period=1.0, station_labels=False)
        plt.close("all")
        assert ax is not None

    def test_external_ax_and_title(self):
        sites = _profile(3)
        fig, ax_in = plt.subplots()
        ax = plot_induction_map(sites, period=1.0, ax=ax_in, title="Custom title")
        plt.close("all")
        assert ax.get_title() == "Custom title"

    def test_explicit_scale_and_reference_arrow(self):
        sites = _profile(3)
        ax = plot_induction_map(sites, period=1.0, scale=5.0, reference_arrow=0.3)
        plt.close("all")
        assert ax is not None

    def test_true_2d_map_branch(self, monkeypatch):
        """Force _station_xy to report real 2-D coordinates.

        Exercises the "map" branch (vs. profile) of the axis-label
        logic and the equal-aspect branch of ``_set_map_aspect``,
        which are unreachable via the public API given the
        ``_station_xy`` limitation documented at the top of this file.
        """
        sites = _profile(4)
        xy_map = {f"S{i:02d}": (float(i) * 100.0, float(i) * 137.0) for i in range(4)}

        def _fake_xy(ed, i):
            from pycsamt.emtools._core import _name

            return xy_map[_name(ed, i)]

        monkeypatch.setattr(tfm, "_station_xy", _fake_xy)
        ax = plot_induction_map(sites, period=1.0)
        plt.close("all")
        assert ax.get_ylabel().startswith("y")


# ─────────────────────────────────────────────────────────────────────────────
# plot_induction_section
# ─────────────────────────────────────────────────────────────────────────────


class TestPlotInductionSection:
    def test_basic_abs(self):
        sites = _profile(4, n=16)
        ax = plot_induction_section(sites)
        plt.close("all")
        assert ax is not None

    def test_component_real_and_imag(self):
        sites = _profile(3, n=16)
        ax_r = plot_induction_section(sites, component="real")
        plt.close("all")
        ax_i = plot_induction_section(sites, component="imag")
        plt.close("all")
        assert ax_r is not None
        assert ax_i is not None

    def test_no_tipper(self):
        sites = _profile(3, with_tipper=False)
        ax = plot_induction_section(sites)
        plt.close("all")
        assert "no tipper" in ax.texts[0].get_text()

    def test_single_period(self):
        sites = _profile(3, n=16)
        ax = plot_induction_section(sites, n_periods=1)
        plt.close("all")
        assert ax is not None

    def test_custom_clim_and_cmap(self):
        sites = _profile(3, n=16)
        ax = plot_induction_section(sites, clim=(0.0, 1.0), cmap="viridis")
        plt.close("all")
        assert ax is not None

    def test_external_ax(self):
        sites = _profile(3, n=16)
        fig, ax_in = plt.subplots()
        ax = plot_induction_section(sites, ax=ax_in)
        plt.close("all")
        assert ax is ax_in

    def test_accepts_airborne_sites(self):
        from pycsamt.airborne import NavigationTrack
        from pycsamt.airborne.site import AirborneSites
        from pycsamt.airborne.ztem import build_ztem_line

        freqs = np.array([30.0, 90.0, 270.0])
        n = 5
        nav = NavigationTrack(
            sample_ids=tuple(f"S{i}" for i in range(n)),
            easting=np.arange(n, dtype=float) * 100.0,
            northing=np.zeros(n),
        )
        tipper = np.zeros((n, freqs.size, 2), dtype=complex)
        tipper[..., 0] = 0.1 + 0.05j
        tipper[..., 1] = 0.05 - 0.02j
        line = build_ztem_line("L001", nav, tipper, frequency=freqs)
        asites = AirborneSites.from_line(line, technology="ztem")

        ax = plot_induction_section(asites)
        plt.close("all")
        assert ax is not None
        assert "no tipper" not in (
            ax.texts[0].get_text() if ax.texts else ""
        )

    def test_custom_title(self):
        sites = _profile(3, n=16)
        ax = plot_induction_section(sites, title="My section")
        plt.close("all")
        assert ax.get_title() == "My section"


# ─────────────────────────────────────────────────────────────────────────────
# plot_induction_convention (extra branches beyond test_emtools_tf.py)
# ─────────────────────────────────────────────────────────────────────────────


class TestPlotInductionConventionExtra:
    def test_no_tipper(self):
        sites = _profile(3, with_tipper=False)
        axs = plot_induction_convention(sites)
        plt.close("all")
        assert axs.shape == (2, 2)

    def test_single_site(self):
        sites = [_site("S00")]
        axs = plot_induction_convention(sites)
        plt.close("all")
        assert axs.shape == (2, 2)

    def test_station_labels_off(self):
        sites = _profile(6)
        axs = plot_induction_convention(sites, station_labels=False)
        plt.close("all")
        assert axs.shape == (2, 2)

    def test_explicit_scale(self):
        sites = _profile(3)
        axs = plot_induction_convention(sites, scale=2.0)
        plt.close("all")
        assert axs.shape == (2, 2)

    def test_true_2d_map_branch(self, monkeypatch):
        sites = _profile(4)
        xy_map = {f"S{i:02d}": (float(i) * 50.0, float(i) * 90.0) for i in range(4)}

        def _fake_xy(ed, i):
            from pycsamt.emtools._core import _name

            return xy_map[_name(ed, i)]

        monkeypatch.setattr(tfm, "_station_xy", _fake_xy)
        axs = plot_induction_convention(sites)
        plt.close("all")
        assert axs.shape == (2, 2)


# ─────────────────────────────────────────────────────────────────────────────
# plot_tipper_polar (extra branches)
# ─────────────────────────────────────────────────────────────────────────────


class TestPlotTipperPolarExtra:
    def test_real_component_data(self):
        sites = [_site("S00", n=16)]
        ax = plot_tipper_polar(sites, component="real")
        plt.close("all")
        assert ax is not None

    def test_imag_component(self):
        sites = [_site("S00", n=16)]
        ax = plot_tipper_polar(sites, component="imag")
        plt.close("all")
        assert ax is not None

    def test_abs_component(self):
        sites = [_site("S00", n=16)]
        ax = plot_tipper_polar(sites, component="abs")
        plt.close("all")
        assert ax is not None

    def test_custom_style_kwargs(self):
        sites = [_site("S00", n=16)]
        ax = plot_tipper_polar(sites, cmap="plasma", lw=2.0, alpha=0.5, title="polar")
        plt.close("all")
        assert ax.get_title() == "polar"

    def test_external_ax(self):
        sites = [_site("S00", n=16)]
        fig = plt.figure()
        ax_in = fig.add_subplot(111, polar=True)
        ax = plot_tipper_polar(sites, ax=ax_in)
        plt.close("all")
        assert ax is ax_in

    def test_single_frequency(self):
        sites = [_site("S00", n=1)]
        ax = plot_tipper_polar(sites)
        plt.close("all")
        assert ax is not None


# ─────────────────────────────────────────────────────────────────────────────
# plot_induction_rose (extra branches)
# ─────────────────────────────────────────────────────────────────────────────


class TestPlotInductionRoseExtra:
    def test_component_variants(self):
        sites = _profile(3, n=16)
        for comp in ("real", "imag", "abs"):
            ax = plot_induction_rose(sites, component=comp)
            plt.close("all")
            assert ax is not None

    def test_pband_filter(self):
        sites = _profile(3, n=16)
        ax = plot_induction_rose(sites, pband=(1.0, 100.0))
        plt.close("all")
        assert ax is not None

    def test_pband_excludes_everything(self):
        sites = _profile(3, n=16)
        ax = plot_induction_rose(sites, pband=(1e9, 2e9))
        plt.close("all")
        assert "no arrows" in ax.get_title()

    def test_custom_nbins(self):
        sites = _profile(3, n=16)
        ax = plot_induction_rose(sites, nbins=12)
        plt.close("all")
        assert ax is not None

    def test_minimal_style_solid_no_mean(self):
        sites = _profile(3, n=16)
        ax = plot_induction_rose(sites, style="minimal")
        plt.close("all")
        assert ax is not None

    def test_external_ax(self):
        sites = _profile(3, n=16)
        fig = plt.figure()
        ax_in = fig.add_subplot(111, projection="polar")
        ax = plot_induction_rose(sites, ax=ax_in)
        plt.close("all")
        assert ax is ax_in

    def test_custom_title(self):
        sites = _profile(3, n=16)
        ax = plot_induction_rose(sites, title="Rose title")
        plt.close("all")
        assert ax.get_title() == "Rose title"


# ─────────────────────────────────────────────────────────────────────────────
# plot_induction_arrows (extra branches beyond test_emtools_tf.py)
# ─────────────────────────────────────────────────────────────────────────────


class TestPlotInductionArrowsExtra:
    def test_real_convention(self):
        sites = _profile(3)
        ax = plot_induction_arrows(sites, periods=[1.0], convention="real")
        plt.close("all")
        assert ax is not None

    def test_imag_convention(self):
        sites = _profile(3)
        ax = plot_induction_arrows(sites, periods=[1.0], convention="imag")
        plt.close("all")
        assert ax is not None

    def test_no_normalize(self):
        sites = _profile(3)
        ax = plot_induction_arrows(sites, periods=[1.0], normalize=False)
        plt.close("all")
        assert ax is not None

    def test_no_strike_ticks(self):
        sites = _profile(3)
        ax = plot_induction_arrows(sites, periods=[1.0], strike_ticks=False)
        plt.close("all")
        assert ax is not None

    def test_true_2d_map_branch(self, monkeypatch):
        sites = _profile(4)
        xy_map = {f"S{i:02d}": (float(i) * 60.0, float(i) * 45.0) for i in range(4)}

        def _fake_xy(ed, i):
            from pycsamt.emtools._core import _name

            return xy_map[_name(ed, i)]

        monkeypatch.setattr(tfm, "_station_xy", _fake_xy)
        ax = plot_induction_arrows(sites, periods=[1.0])
        plt.close("all")
        assert ax.get_xlabel().startswith("East")


# ─────────────────────────────────────────────────────────────────────────────
# Spectra-direct wrappers
# ─────────────────────────────────────────────────────────────────────────────


class TestPlotInductionMapFromSpectra:
    def test_basic_dict_with_2d_coords(self):
        sp, fr = _spectra_dict(4)
        coords = {
            name: (float(i) * 100.0, float(i) * 33.0) for i, name in enumerate(sp)
        }
        ax = plot_induction_map_from_spectra(
            sp, period=1.0 / fr[len(fr) // 2], coords=coords
        )
        plt.close("all")
        assert ax is not None

    def test_list_input_auto_coords(self):
        sp_dict, fr = _spectra_dict(3)
        sp_list = list(sp_dict.values())
        ax = plot_induction_map_from_spectra(sp_list, period=1.0)
        plt.close("all")
        assert ax is not None

    def test_no_tipper(self):
        sp, _ = _spectra_dict(3, with_tipper=False)
        ax = plot_induction_map_from_spectra(sp, period=1.0)
        plt.close("all")
        assert "no tipper" in ax.texts[0].get_text()

    def test_show_real_only_and_show_imag_only(self):
        sp, _ = _spectra_dict(3)
        ax1 = plot_induction_map_from_spectra(sp, period=1.0, show_imag=False)
        plt.close("all")
        ax2 = plot_induction_map_from_spectra(sp, period=1.0, show_real=False)
        plt.close("all")
        assert ax1 is not None
        assert ax2 is not None

    def test_station_labels_off_and_external_ax(self):
        sp, _ = _spectra_dict(6)
        fig, ax_in = plt.subplots()
        ax = plot_induction_map_from_spectra(
            sp, period=1.0, ax=ax_in, station_labels=False
        )
        plt.close("all")
        assert ax is ax_in

    def test_single_site_scale_fallback(self):
        sp, _ = _spectra_dict(1)
        ax = plot_induction_map_from_spectra(sp, period=1.0)
        plt.close("all")
        assert ax is not None

    def test_custom_title(self):
        sp, _ = _spectra_dict(3)
        ax = plot_induction_map_from_spectra(sp, period=1.0, title="Spectra map")
        plt.close("all")
        assert ax.get_title() == "Spectra map"


class TestPlotTipperPolarFromSpectra:
    def test_basic(self):
        sp, _ = _spectra_dict(1)
        ax = plot_tipper_polar_from_spectra(sp)
        plt.close("all")
        assert ax is not None

    def test_component_variants(self):
        sp, _ = _spectra_dict(1)
        for comp in ("real", "imag", "abs"):
            ax = plot_tipper_polar_from_spectra(sp, component=comp)
            plt.close("all")
            assert ax is not None

    def test_no_tipper(self):
        sp, _ = _spectra_dict(1, with_tipper=False)
        ax = plot_tipper_polar_from_spectra(sp)
        plt.close("all")
        assert ax.get_title() == "no tipper"

    def test_no_tipper_external_ax(self):
        sp, _ = _spectra_dict(1, with_tipper=False)
        fig = plt.figure()
        ax_in = fig.add_subplot(111, polar=True)
        ax = plot_tipper_polar_from_spectra(sp, ax=ax_in)
        plt.close("all")
        assert ax is ax_in

    def test_external_ax_with_tipper(self):
        sp, _ = _spectra_dict(1)
        fig = plt.figure()
        ax_in = fig.add_subplot(111, polar=True)
        ax = plot_tipper_polar_from_spectra(sp, ax=ax_in)
        plt.close("all")
        assert ax is ax_in

    def test_single_spectra_object_dict(self):
        sp, _ = _spectra_dict(1)
        ax = plot_tipper_polar_from_spectra(sp, title="One site")
        plt.close("all")
        assert ax.get_title() == "One site"


class TestPlotInductionRoseFromSpectra:
    def test_basic(self):
        sp, _ = _spectra_dict(4)
        ax = plot_induction_rose_from_spectra(sp)
        plt.close("all")
        assert ax is not None

    def test_component_variants(self):
        sp, _ = _spectra_dict(3)
        for comp in ("real", "imag", "abs"):
            ax = plot_induction_rose_from_spectra(sp, component=comp)
            plt.close("all")
            assert ax is not None

    def test_pband_filter(self):
        sp, fr = _spectra_dict(3)
        per = 1.0 / fr
        ax = plot_induction_rose_from_spectra(
            sp, pband=(float(per.min()), float(per.max()))
        )
        plt.close("all")
        assert ax is not None

    def test_no_arrows(self):
        sp, _ = _spectra_dict(3, with_tipper=False)
        ax = plot_induction_rose_from_spectra(sp)
        plt.close("all")
        assert "no arrows" in ax.get_title()

    def test_custom_style_and_nbins(self):
        sp, _ = _spectra_dict(3)
        ax = plot_induction_rose_from_spectra(sp, style="minimal", nbins=8)
        plt.close("all")
        assert ax is not None

    def test_external_ax(self):
        sp, _ = _spectra_dict(3)
        fig = plt.figure()
        ax_in = fig.add_subplot(111, projection="polar")
        ax = plot_induction_rose_from_spectra(sp, ax=ax_in)
        plt.close("all")
        assert ax is ax_in


# ─────────────────────────────────────────────────────────────────────────────
# plot_induction_multiperiod_map
# ─────────────────────────────────────────────────────────────────────────────


class TestPlotTipperHodogramsExtra:
    def test_explicit_bands(self):
        sites = [_site("S00", n=16)]
        fig = plot_tipper_hodograms(sites, bands=[(0.01, 1.0), (1.0, 100.0)])
        plt.close("all")
        assert len(fig.axes) == 2

    def test_normalize(self):
        sites = [_site("S00", n=16)]
        fig = plot_tipper_hodograms(sites, normalize=True)
        plt.close("all")
        assert len(fig.axes) == 2


class TestThinLabelIndicesViaManyStations:
    def test_many_stations_thin_labels(self):
        sites = _profile(40)
        ax = plot_induction_map(sites, period=1.0, station_labels=True)
        plt.close("all")
        assert ax is not None

    def test_many_stations_convention_thin_labels(self):
        sites = _profile(40)
        axs = plot_induction_convention(sites, station_labels=True)
        plt.close("all")
        assert axs.shape == (2, 2)


class TestPlotInductionMultiperiodMap:
    def test_basic_default_background(self):
        sites = _profile(4)
        fig, axs = plot_induction_multiperiod_map(sites, periods=[1.0, 10.0])
        plt.close("all")
        assert len(axs) == 2

    def test_no_sites(self):
        fig, axs = plot_induction_multiperiod_map([], periods=[1.0])
        plt.close("all")
        assert axs.size == 1
        assert "no sites" in axs[0].texts[0].get_text()

    def test_single_site(self):
        sites = [_site("S00")]
        fig, axs = plot_induction_multiperiod_map(sites, periods=[1.0])
        plt.close("all")
        assert len(axs) == 1

    def test_explicit_tipper_data_matching_shape(self):
        sites = _profile(3)
        tipper_data = {
            1.0: np.array(
                [[0.1 + 0.05j, 0.02 + 0.01j] for _ in range(3)],
                dtype=complex,
            )
        }
        fig, axs = plot_induction_multiperiod_map(
            sites, periods=[1.0], tipper_data=tipper_data
        )
        plt.close("all")
        assert len(axs) == 1

    def test_explicit_tipper_data_mismatched_shape_falls_back(self):
        sites = _profile(3)
        tipper_data = {1.0: np.array([[0.1 + 0.05j, 0.02 + 0.01j]], dtype=complex)}
        fig, axs = plot_induction_multiperiod_map(
            sites, periods=[1.0], tipper_data=tipper_data
        )
        plt.close("all")
        assert len(axs) == 1

    def test_custom_background_and_extent(self):
        sites = _profile(3)
        bg = np.random.default_rng(0).normal(size=(20, 25))
        ext = (0.0, 600.0, -50.0, 50.0)
        fig, axs = plot_induction_multiperiod_map(
            sites,
            periods=[1.0],
            background=bg,
            background_extent=ext,
            background_clim=(-1.0, 1.0),
        )
        plt.close("all")
        assert len(axs) == 1

    def test_wiese_convention(self):
        sites = _profile(3)
        fig, axs = plot_induction_multiperiod_map(
            sites, periods=[1.0], convention="wiese"
        )
        plt.close("all")
        assert len(axs) == 1

    def test_bg_colorbar_sides(self):
        sites = _profile(3)
        for side in ("right", "left", "top", "bottom"):
            fig, axs = plot_induction_multiperiod_map(
                sites, periods=[1.0], bg_colorbar_side=side
            )
            plt.close("all")
            assert len(axs) == 1

    def test_bg_colorbar_side_invalid_raises(self):
        sites = _profile(3)
        with pytest.raises(ValueError):
            plot_induction_multiperiod_map(
                sites, periods=[1.0], bg_colorbar_side="diagonal"
            )
        plt.close("all")

    def test_show_background_cbar_false(self):
        sites = _profile(3)
        fig, axs = plot_induction_multiperiod_map(
            sites, periods=[1.0], show_background_cbar=False
        )
        plt.close("all")
        assert len(axs) == 1

    def test_station_labels_and_show_stations_false(self):
        sites = _profile(5)
        fig, axs = plot_induction_multiperiod_map(
            sites, periods=[1.0], station_labels=True
        )
        plt.close("all")
        fig2, axs2 = plot_induction_multiperiod_map(
            sites, periods=[1.0], show_stations=False
        )
        plt.close("all")
        assert len(axs) == 1
        assert len(axs2) == 1

    def test_annotations(self):
        sites = _profile(3)
        annotations = {
            "Fault A": (100.0, 0.0),
            "Fault B": (300.0, 0.0, {"color": "red", "fontsize": 10}),
        }
        fig, axs = plot_induction_multiperiod_map(
            sites, periods=[1.0], annotations=annotations
        )
        plt.close("all")
        assert len(axs) == 1

    def test_reference_panel_and_label(self):
        sites = _profile(3)
        fig, axs = plot_induction_multiperiod_map(
            sites,
            periods=[1.0, 10.0],
            reference_panel=1,
            reference_label="Ref",
        )
        plt.close("all")
        assert len(axs) == 2

    def test_external_axes(self):
        sites = _profile(3)
        fig_in, axs_in = plt.subplots(2, 1)
        fig, axs = plot_induction_multiperiod_map(
            sites, periods=[1.0, 10.0], axes=axs_in
        )
        plt.close("all")
        assert len(axs) == 2

    def test_title_and_panel_labels(self):
        sites = _profile(3)
        fig, axs = plot_induction_multiperiod_map(
            sites,
            periods=[1.0, 10.0],
            title="Multi-period map",
            panel_labels=["(A)", "(B)"],
        )
        plt.close("all")
        assert fig._suptitle.get_text() == "Multi-period map"

    def test_no_tipper_zero_arrows(self):
        sites = _profile(3, with_tipper=False)
        fig, axs = plot_induction_multiperiod_map(sites, periods=[1.0])
        plt.close("all")
        assert len(axs) == 1

    def test_many_periods(self):
        sites = _profile(3)
        fig, axs = plot_induction_multiperiod_map(
            sites, periods=[0.1, 1.0, 10.0, 100.0, 1000.0]
        )
        plt.close("all")
        assert len(axs) == 5
