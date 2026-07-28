"""
Additional coverage tests for pycsamt.emtools.ss.

Complements test_emtools_ss_ns.py by exercising the static-shift
estimators (AMA/LOESS/bilateral/reference-median), the apply/correct
helpers, the QC plotting functions (delta pseudosection, station
curves, delta profile, publication-quality comparison/summary plots,
one-shot QC wrappers, radar), and additional detect_near_surface /
plot_ns_detection branches (ordering helpers, skew/period masking,
wrapped Site objects).
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from pycsamt.api import APIFrame, reset_api_view
from pycsamt.emtools.ss import (
    _scale_site_Z,
    apply_ss_factors,
    correct_ss_ama,
    detect_near_surface,
    estimate_ss_ama,
    estimate_ss_bilateral,
    estimate_ss_loess,
    estimate_ss_refmedian,
    plot_ns_detection,
    plot_ss_1d_curves,
    plot_ss_comparison_psection,
    plot_ss_delta_profile,
    plot_ss_delta_psection,
    plot_ss_radar,
    plot_ss_station_curves,
    plot_ss_summary,
    ss_comparison_psection,
    ss_qc_profile,
    ss_qc_psection,
    ss_qc_station_curves,
)

# ----------------------------- fixtures ----------------------------------- #


class _FakeZ:
    def __init__(self, z, freq):
        self.z = np.asarray(z, dtype=complex)
        self.freq = np.asarray(freq, dtype=float)


class _FakeSite:
    """Minimal EDI-like object (duck-type compatible with is_edi_file)."""

    def __init__(self, station, z, freq, lon=None, lat=None, coords=None):
        self.station = station
        self.Z = _FakeZ(z, freq)
        self.freq = np.asarray(freq, dtype=float)
        if lon is not None:
            self.lon = float(lon)
        if lat is not None:
            self.lat = float(lat)
        if coords is not None:
            self.coords = coords

    def get_section(self, *_, **__):
        return None


class _WrappedSite:
    """Sites-level ``Site`` wrapper: real EDI lives at ``.edi``."""

    def __init__(self, edi):
        self.edi = edi
        self.station = "wrapped-should-not-be-used"

    def get_section(self, *_, **__):
        return None


class _NoZSite:
    """Passes ``is_edi_file`` duck-typing (has ``Z`` + ``get_section``) but
    carries no usable impedance data -- exercises the ``Z is None`` guard
    branches sprinkled through ss.py's per-site loops."""

    def __init__(self, station):
        self.station = station
        self.Z = None

    def get_section(self, *_, **__):
        return None


# pycsamt unit convention: ρ_a = 0.2 * |Z|² / f  =>  |Z| = sqrt(5 f ρ)
def _make_z(freqs: np.ndarray, rho) -> np.ndarray:
    freqs = np.asarray(freqs, dtype=float)
    rho_arr = np.broadcast_to(np.asarray(rho, dtype=float), freqs.shape)
    z_abs = np.sqrt(5.0 * freqs * rho_arr)
    z = np.zeros((freqs.size, 2, 2), dtype=complex)
    z[:, 0, 1] = z_abs * (1 + 1j) / np.sqrt(2)
    z[:, 1, 0] = -z_abs * (1 + 1j) / np.sqrt(2)
    return z


def _site(
    station: str,
    rho: float = 100.0,
    n: int = 10,
    f_lo: float = -2,
    f_hi: float = 3,
    **kw,
) -> _FakeSite:
    fr = np.logspace(f_lo, f_hi, n)
    return _FakeSite(station, _make_z(fr, rho), fr, **kw)


def _profile(n: int = 8, rho: float = 100.0, n_freq: int = 10, **kw):
    return [_site(f"S{i:02d}", rho=rho, n=n_freq, **kw) for i in range(n)]


def _shifted_profile(n: int = 8, base_rho: float = 100.0, shift: float = 0.3):
    """Profile where every other station is uniformly shifted."""
    out = []
    for i in range(n):
        rho = base_rho * (10.0**shift) if i % 3 == 0 else base_rho
        out.append(_site(f"S{i:02d}", rho=rho, n=10))
    return out


# ═══════════════════════════════════════════════════════════════════════════
# _site_coords / _site_order_key (via sort_by="lon"/"lat")
# ═══════════════════════════════════════════════════════════════════════════


def test_estimate_ss_ama_sort_by_lon_with_mixed_coord_sources():
    sites = [
        _site("A", lon=1.0, lat=1.0),
        _site("B", coords=(2.0, 2.0)),  # (lat, lon) fallback
        _site("C"),  # neither -> (1, inf, name) fallback
        _site("D", coords=(None, None)),  # exception branch -> None
        _site("E", lon=0.5, lat=0.5),
    ]
    tbl = estimate_ss_ama(sites, sort_by="lon", half_window=2)
    assert set(tbl["station"]).issubset({"A", "B", "C", "D", "E"})


def test_estimate_ss_ama_sort_by_lat_with_mixed_coord_sources():
    sites = [
        _site("A", lon=1.0, lat=3.0),
        _site("B", coords=(2.0, 2.0)),
        _site("C"),
        _site("D", coords=(None, None)),
        _site("E", lon=0.5, lat=0.1),
    ]
    tbl = estimate_ss_ama(sites, sort_by="lat", half_window=2)
    assert set(tbl["station"]).issubset({"A", "B", "C", "D", "E"})


def test_estimate_ss_ama_sort_by_name():
    sites = _profile(5)
    tbl = estimate_ss_ama(sites, sort_by="name", half_window=2)
    assert list(tbl["station"]) == sorted(tbl["station"])


# ═══════════════════════════════════════════════════════════════════════════
# estimate_ss_ama — main computational path
# ═══════════════════════════════════════════════════════════════════════════


@pytest.mark.parametrize("weights", ["tri", "gauss", "uniform"])
def test_estimate_ss_ama_weight_schemes(weights):
    sites = _shifted_profile(9)
    tbl = estimate_ss_ama(sites, half_window=2, weights=weights)
    assert len(tbl) > 0
    assert {"delta_log10_rho", "fac_rho", "fac_z", "n_used"}.issubset(tbl.columns)


@pytest.mark.parametrize("robust_freq", ["median", "mean"])
@pytest.mark.parametrize("robust_overall", ["median", "mean"])
def test_estimate_ss_ama_robust_combinations(robust_freq, robust_overall):
    sites = _shifted_profile(7)
    tbl = estimate_ss_ama(
        sites,
        half_window=2,
        robust_freq=robust_freq,
        robust_overall=robust_overall,
    )
    assert len(tbl) > 0
    assert np.isfinite(tbl["delta_log10_rho"]).any()


def test_estimate_ss_ama_pband_filters_frequencies():
    sites = _profile(6, n_freq=12)
    tbl = estimate_ss_ama(sites, half_window=2, pband=(1e-2, 1.0))
    assert len(tbl) > 0


def test_estimate_ss_ama_max_skew_none_skips_masking():
    sites = _shifted_profile(6)
    tbl = estimate_ss_ama(sites, half_window=2, max_skew=None)
    assert len(tbl) > 0


def test_estimate_ss_ama_max_skew_excludes_all_frequencies():
    """An impossible max_skew (<0) masks every frequency for every site.

    Exercises the ``fr1.size == 0: continue`` branch while still routing
    through the skew-masking code path (pt is non-empty, sdf is non-empty).
    """
    sites = _profile(5)
    tbl = estimate_ss_ama(sites, half_window=2, max_skew=-1.0)
    assert tbl.empty or len(tbl) >= 0  # degrades gracefully, no crash


def test_estimate_ss_ama_api_wrap_nonempty():
    reset_api_view()
    sites = _shifted_profile(6)
    viewed = estimate_ss_ama(sites, half_window=2, api=True)
    assert isinstance(viewed, APIFrame)
    assert viewed.kind == "emtools.ss.ama"
    assert not viewed.df.empty


# ═══════════════════════════════════════════════════════════════════════════
# apply_ss_factors / correct_ss_ama
# ═══════════════════════════════════════════════════════════════════════════


def test_apply_ss_factors_dataframe_form():
    sites = [_site("S1"), _site("S2")]
    factors = pd.DataFrame({"station": ["S1", "S2"], "fac_z": [1.5, 0.8]})
    corr = apply_ss_factors(sites, factors, key="fac_z", inplace=False)
    items = list(corr) if not hasattr(corr, "items") else corr
    assert items is not None


def test_apply_ss_factors_bad_dataframe_raises():
    sites = [_site("S1")]
    bad = pd.DataFrame({"foo": [1], "bar": [2]})
    with pytest.raises(ValueError):
        apply_ss_factors(sites, bad, key="fac_z", inplace=False)


def test_apply_ss_factors_unknown_station_defaults_to_one():
    sites = [_site("S1"), _site("S2")]
    corr = apply_ss_factors(sites, {"NOT_S1": 5.0}, key="fac_z", inplace=False)
    assert corr is not None


def test_correct_ss_ama_end_to_end():
    sites = _shifted_profile(8)
    corrected = correct_ss_ama(sites, half_window=2, inplace=False)
    assert corrected is not None
    # corrected sites should still carry finite impedance
    from pycsamt.emtools._core import _get_z_block, _iter_items

    n_finite = 0
    for ed in _iter_items(corrected):
        _, z, _ = _get_z_block(ed)
        if z is not None:
            n_finite += int(np.isfinite(z).sum())
    assert n_finite > 0


def test_correct_ss_ama_inplace_true():
    sites = _shifted_profile(6)
    corrected = correct_ss_ama(sites, half_window=2, inplace=True)
    assert corrected is not None


# ═══════════════════════════════════════════════════════════════════════════
# estimate_ss_loess
# ═══════════════════════════════════════════════════════════════════════════


@pytest.mark.parametrize("poly", [0, 1])
@pytest.mark.parametrize("summary", ["median", "mean"])
def test_estimate_ss_loess_main_path(poly, summary):
    sites = _shifted_profile(9)
    tbl = estimate_ss_loess(sites, half_window=3, poly=poly, it=2, summary=summary)
    assert len(tbl) > 0
    assert {"station", "delta_log10_rho", "fac_rho", "fac_z", "n_used"} == (
        set(tbl.columns)
    )


def test_estimate_ss_loess_it_one_iteration():
    sites = _shifted_profile(7)
    tbl = estimate_ss_loess(sites, half_window=2, it=1)
    assert len(tbl) > 0


def test_estimate_ss_loess_pband_and_skew():
    sites = _profile(6, n_freq=12)
    tbl = estimate_ss_loess(sites, half_window=2, pband=(1e-2, 10.0), max_skew=None)
    assert len(tbl) > 0


def test_estimate_ss_loess_single_station():
    tbl = estimate_ss_loess([_site("Solo")])
    assert len(tbl) == 1


def test_estimate_ss_loess_api_wrap():
    reset_api_view()
    sites = _shifted_profile(6)
    viewed = estimate_ss_loess(sites, api=True)
    assert isinstance(viewed, APIFrame)
    assert viewed.kind == "emtools.ss.loess"


# ═══════════════════════════════════════════════════════════════════════════
# estimate_ss_bilateral
# ═══════════════════════════════════════════════════════════════════════════


@pytest.mark.parametrize("summary", ["median", "mean"])
def test_estimate_ss_bilateral_main_path(summary):
    sites = _shifted_profile(9)
    tbl = estimate_ss_bilateral(sites, half_window=3, summary=summary)
    assert len(tbl) > 0


def test_estimate_ss_bilateral_explicit_sigmas():
    sites = _shifted_profile(7)
    tbl = estimate_ss_bilateral(sites, half_window=2, sig_dist=1.5, sig_val=0.2)
    assert len(tbl) > 0


def test_estimate_ss_bilateral_single_station():
    tbl = estimate_ss_bilateral([_site("Solo")])
    assert len(tbl) == 1


def test_estimate_ss_bilateral_api_wrap():
    reset_api_view()
    sites = _shifted_profile(6)
    viewed = estimate_ss_bilateral(sites, api=True)
    assert isinstance(viewed, APIFrame)
    assert viewed.kind == "emtools.ss.bilateral"


# ═══════════════════════════════════════════════════════════════════════════
# estimate_ss_refmedian
# ═══════════════════════════════════════════════════════════════════════════


@pytest.mark.parametrize("summary", ["median", "mean"])
def test_estimate_ss_refmedian_main_path(summary):
    sites = _shifted_profile(8)
    tbl = estimate_ss_refmedian(sites, summary=summary)
    assert len(tbl) > 0


def test_estimate_ss_refmedian_smooth_sites_noop_branch():
    sites = _shifted_profile(6)
    tbl = estimate_ss_refmedian(sites, smooth_sites=3)
    assert len(tbl) > 0


def test_estimate_ss_refmedian_single_station():
    tbl = estimate_ss_refmedian([_site("Solo")])
    assert len(tbl) == 1


def test_estimate_ss_refmedian_api_wrap():
    reset_api_view()
    sites = _shifted_profile(6)
    viewed = estimate_ss_refmedian(sites, api=True)
    assert isinstance(viewed, APIFrame)
    assert viewed.kind == "emtools.ss.refmedian"


# ═══════════════════════════════════════════════════════════════════════════
# plot_ss_delta_psection / plot_ss_station_curves / plot_ss_delta_profile
# ═══════════════════════════════════════════════════════════════════════════


@pytest.fixture(autouse=True, scope="module")
def _agg_backend():
    import matplotlib

    matplotlib.use("Agg")
    yield


def _before_after(n=6, shift=0.25):
    before = _profile(n, rho=100.0, n_freq=10)
    after = []
    for i, s in enumerate(before):
        rho = 100.0 * (10.0**shift if i == 0 else 1.0)
        after.append(_site(s.station, rho=rho, n=10))
    return before, after


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_plot_ss_delta_psection_main_path():
    import matplotlib.pyplot as plt

    before, after = _before_after(6)
    ax = plot_ss_delta_psection(before, after)
    assert ax is not None
    plt.close("all")


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_plot_ss_delta_psection_period_axis_and_vlim():
    import matplotlib.pyplot as plt

    before, after = _before_after(5)
    ax = plot_ss_delta_psection(
        before, after, axis_y="period", vlim=0.4, pband=(1e-2, 1e2)
    )
    assert ax is not None
    plt.close("all")


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_plot_ss_delta_psection_no_overlap():
    import matplotlib.pyplot as plt

    before = [_site("A")]
    after = [_site("Z")]
    ax = plot_ss_delta_psection(before, after)
    assert ax is not None
    plt.close("all")


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_plot_ss_station_curves_default_station():
    import matplotlib.pyplot as plt

    before, after = _before_after(4)
    ax = plot_ss_station_curves(before, after)
    assert ax is not None
    plt.close("all")


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_plot_ss_station_curves_named_station_and_pband():
    import matplotlib.pyplot as plt

    before, after = _before_after(4)
    ax = plot_ss_station_curves(before, after, station="S02", pband=(1e-2, 1e2))
    assert ax is not None
    plt.close("all")


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_plot_ss_station_curves_no_common_stations():
    import matplotlib.pyplot as plt

    ax = plot_ss_station_curves([_site("A")], [_site("Z")])
    assert ax is not None
    plt.close("all")


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_plot_ss_station_curves_station_not_found():
    import matplotlib.pyplot as plt

    before, after = _before_after(3)
    ax = plot_ss_station_curves(before, after, station="does-not-exist")
    assert ax is not None
    plt.close("all")


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_plot_ss_delta_profile_main_path():
    import matplotlib.pyplot as plt

    before, after = _before_after(6)
    ax = plot_ss_delta_profile(before, after, robust="mean")
    assert ax is not None
    plt.close("all")


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_plot_ss_delta_profile_no_overlap():
    import matplotlib.pyplot as plt

    ax = plot_ss_delta_profile([_site("A")], [_site("Z")])
    assert ax is not None
    plt.close("all")


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_plot_ss_delta_profile_existing_axes():
    import matplotlib.pyplot as plt

    _, ax_in = plt.subplots()
    before, after = _before_after(4)
    ax_out = plot_ss_delta_profile(before, after, ax=ax_in)
    assert ax_out is ax_in
    plt.close("all")


# ═══════════════════════════════════════════════════════════════════════════
# Publication-quality comparison / 1-D curves / summary (array-based)
# ═══════════════════════════════════════════════════════════════════════════


def _synthetic_arrays(n_st=6, n_f=10, seed=0):
    rng = np.random.default_rng(seed)
    freqs = np.logspace(-2, 3, n_f)
    logRho_before = rng.normal(2.0, 0.25, size=(n_st, n_f))
    shift = rng.normal(0.0, 0.2, size=(n_st, 1))
    logRho_after = logRho_before + shift
    labels = [f"S{i:02d}" for i in range(n_st)]
    return logRho_before, logRho_after, freqs, labels


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_plot_ss_comparison_psection_with_delta():
    import matplotlib.pyplot as plt

    b, a, freqs, labels = _synthetic_arrays()
    fig = plot_ss_comparison_psection(
        b, a, freqs=freqs, station_labels=labels, show_delta=True
    )
    assert fig is not None
    assert len(fig.axes) >= 3
    plt.close("all")


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_plot_ss_comparison_psection_no_delta():
    import matplotlib.pyplot as plt

    b, a, freqs, labels = _synthetic_arrays()
    fig = plot_ss_comparison_psection(
        b, a, freqs=freqs, station_labels=labels, show_delta=False
    )
    assert fig is not None
    plt.close("all")


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_plot_ss_comparison_psection_explicit_clim_and_axes():
    import matplotlib.pyplot as plt

    b, a, freqs, labels = _synthetic_arrays()
    fig, axes = plt.subplots(3, 1, figsize=(8, 10))
    out = plot_ss_comparison_psection(
        b,
        a,
        freqs=freqs,
        station_labels=labels,
        clim=(1.0, 3.0),
        delta_vlim=0.5,
        period_up=False,
        suptitle="test",
        axes=axes,
    )
    assert out is fig
    plt.close("all")


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_plot_ss_1d_curves_default_selection():
    import matplotlib.pyplot as plt

    b, a, freqs, labels = _synthetic_arrays(n_st=5, n_f=8)
    fig = plot_ss_1d_curves(b, a, freqs=freqs, station_labels=labels)
    assert fig is not None
    plt.close("all")


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_plot_ss_1d_curves_int_station_selection():
    import matplotlib.pyplot as plt

    b, a, freqs, labels = _synthetic_arrays(n_st=6, n_f=8)
    fig = plot_ss_1d_curves(
        b, a, freqs=freqs, station_labels=labels, stations=[0, 2, 5]
    )
    assert fig is not None
    plt.close("all")


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_plot_ss_1d_curves_str_station_selection_with_fallback():
    import matplotlib.pyplot as plt

    b, a, freqs, labels = _synthetic_arrays(n_st=4, n_f=8)
    # "does-not-exist" fails lookup -> triggers `idx` empty -> fallback range
    fig = plot_ss_1d_curves(
        b, a, freqs=freqs, station_labels=labels, stations=["does-not-exist"]
    )
    assert fig is not None
    plt.close("all")


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_plot_ss_1d_curves_max_stations_cap():
    import matplotlib.pyplot as plt

    b, a, freqs, labels = _synthetic_arrays(n_st=20, n_f=6)
    fig = plot_ss_1d_curves(
        b,
        a,
        freqs=freqs,
        station_labels=labels,
        max_stations=5,
        n_cols=3,
        log_period=False,
        show_shift_annotation=False,
        show_grid=False,
        title="cap test",
    )
    assert fig is not None
    plt.close("all")


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_plot_ss_summary_default():
    import matplotlib.pyplot as plt

    b, a, freqs, labels = _synthetic_arrays(n_st=5, n_f=8)
    fig = plot_ss_summary(b, a, freqs=freqs, station_labels=labels)
    assert fig is not None
    assert len(fig.axes) >= 4
    plt.close("all")


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_plot_ss_summary_mean_shift_and_explicit_axes():
    import matplotlib.pyplot as plt

    b, a, freqs, labels = _synthetic_arrays(n_st=5, n_f=8)
    fig = plt.figure(figsize=(10, 12))
    ax_before = fig.add_subplot(3, 2, 1)
    ax_after = fig.add_subplot(3, 2, 2)
    ax_delta = fig.add_subplot(3, 1, 2)
    ax_bar = fig.add_subplot(3, 1, 3)
    out = plot_ss_summary(
        b,
        a,
        freqs=freqs,
        station_labels=labels,
        shift_robust="mean",
        axes=(ax_before, ax_after, ax_delta, ax_bar),
    )
    assert out is fig
    plt.close("all")


# ═══════════════════════════════════════════════════════════════════════════
# One-shot QC wrappers (sites in)
# ═══════════════════════════════════════════════════════════════════════════


@pytest.mark.parametrize("method", ["ama", "loess", "bilateral", "refmedian"])
@pytest.mark.filterwarnings("ignore::UserWarning")
def test_ss_qc_psection_all_methods(method):
    import matplotlib.pyplot as plt

    sites = _shifted_profile(7)
    ax = ss_qc_psection(sites, method=method, half_window=2)
    assert ax is not None
    plt.close("all")


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_ss_qc_psection_return_sites():
    import matplotlib.pyplot as plt

    sites = _shifted_profile(6)
    ax, corrected = ss_qc_psection(
        sites, method="ama", half_window=2, return_sites=True
    )
    assert ax is not None
    assert corrected is not None
    plt.close("all")


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_ss_qc_station_curves_default():
    import matplotlib.pyplot as plt

    sites = _shifted_profile(6)
    ax, corrected = ss_qc_station_curves(
        sites, method="loess", half_window=2, return_sites=True
    )
    assert ax is not None
    assert corrected is not None
    plt.close("all")


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_ss_qc_profile_default():
    import matplotlib.pyplot as plt

    sites = _shifted_profile(6)
    ax = ss_qc_profile(sites, method="bilateral", half_window=2)
    assert ax is not None
    plt.close("all")


def test_ss_qc_wrappers_unknown_method_raises():
    sites = _shifted_profile(4)
    with pytest.raises(ValueError):
        ss_qc_psection(sites, method="not-a-method")


@pytest.mark.parametrize("method", ["ama", "refmedian"])
@pytest.mark.filterwarnings("ignore::UserWarning")
def test_ss_comparison_psection_methods(method):
    import matplotlib.pyplot as plt

    sites = _shifted_profile(6)
    fig = ss_comparison_psection(sites, method=method, half_window=2)
    assert fig is not None
    plt.close("all")


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_ss_comparison_psection_return_sites_and_no_delta():
    import matplotlib.pyplot as plt

    sites = _shifted_profile(6)
    fig, corrected = ss_comparison_psection(
        sites,
        method="ama",
        half_window=2,
        show_delta=False,
        return_sites=True,
    )
    assert fig is not None
    assert corrected is not None
    plt.close("all")


# ═══════════════════════════════════════════════════════════════════════════
# plot_ss_radar
# ═══════════════════════════════════════════════════════════════════════════


@pytest.mark.parametrize("rotate", ["pt", "deg", "none"])
@pytest.mark.filterwarnings("ignore::UserWarning")
def test_plot_ss_radar_rotation_modes(rotate):
    import matplotlib.pyplot as plt

    sites = _profile(4)
    ax = plot_ss_radar(sites, rotate=rotate, rotate_deg=20.0)
    assert ax is not None
    plt.close("all")


@pytest.mark.parametrize("theta_axis", ["logperiod", "period", "freq"])
@pytest.mark.filterwarnings("ignore::UserWarning")
def test_plot_ss_radar_theta_axes(theta_axis):
    import matplotlib.pyplot as plt

    sites = _profile(3)
    ax = plot_ss_radar(sites, theta_axis=theta_axis)
    assert ax is not None
    plt.close("all")


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_plot_ss_radar_radial_linear_and_fill():
    import matplotlib.pyplot as plt

    sites = _profile(3)
    ax = plot_ss_radar(sites, radial="rho", fill_between=True)
    assert ax is not None
    plt.close("all")


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_plot_ss_radar_named_station_and_pband():
    import matplotlib.pyplot as plt

    sites = _profile(4)
    ax = plot_ss_radar(sites, station="S02", pband=(1e-2, 1e2))
    assert ax is not None
    plt.close("all")


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_plot_ss_radar_empty_sites():
    import matplotlib.pyplot as plt

    ax = plot_ss_radar([])
    assert ax is not None
    plt.close("all")


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_plot_ss_radar_station_not_found():
    import matplotlib.pyplot as plt

    sites = _profile(2)
    ax = plot_ss_radar(sites, station="nonexistent")
    assert ax is not None
    plt.close("all")


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_plot_ss_radar_reuses_polar_axes():
    import matplotlib.pyplot as plt

    _, polar_ax = plt.subplots(subplot_kw={"polar": True})
    sites = _profile(2)
    ax_out = plot_ss_radar(sites, ax=polar_ax)
    assert ax_out is polar_ax
    plt.close("all")


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_plot_ss_radar_converts_non_polar_axes():
    import matplotlib.pyplot as plt

    _, cart_ax = plt.subplots()
    sites = _profile(2)
    ax_out = plot_ss_radar(sites, ax=cart_ax)
    assert getattr(ax_out, "name", "") == "polar"
    plt.close("all")


# ═══════════════════════════════════════════════════════════════════════════
# detect_near_surface — extra branches
# ═══════════════════════════════════════════════════════════════════════════


def test_detect_near_surface_pband_filter():
    sites = _profile(5, n_freq=14)
    df = detect_near_surface(sites, pband=(1e-2, 10.0))
    assert len(df) == 5


def test_detect_near_surface_max_skew_none():
    sites = _profile(5)
    df = detect_near_surface(sites, max_skew=None)
    assert len(df) == 5


def test_unwrap_ns_returns_inner_edi_when_present():
    """``_unwrap_ns`` should follow ``.edi`` when it exposes ``Z``."""
    from pycsamt.emtools.ss import _unwrap_ns

    inner = _site("W00")
    wrapped = _WrappedSite(inner)
    assert _unwrap_ns(wrapped) is inner


def test_unwrap_ns_passthrough_without_edi():
    from pycsamt.emtools.ss import _unwrap_ns

    plain = _site("W01")
    assert _unwrap_ns(plain) is plain


def test_detect_near_surface_single_frequency_per_band():
    """Single point on each side of f_split -> slopes undefined (NaN)."""
    fr = np.array([0.5, 2.0])
    site = _FakeSite("S1", _make_z(fr, 100.0), fr)
    df = detect_near_surface([site], f_split=1.0)
    assert len(df) == 1
    row = df.iloc[0]
    assert np.isnan(row["slope_hf"]) or np.isnan(row["slope_lf"])


def test_detect_near_surface_api_wrap():
    reset_api_view()
    sites = _profile(5)
    viewed = detect_near_surface(sites, api=True)
    assert isinstance(viewed, APIFrame)
    assert viewed.kind == "emtools.ss.near_surface"
    assert not viewed.df.empty


# ═══════════════════════════════════════════════════════════════════════════
# plot_ns_detection — extra branches
# ═══════════════════════════════════════════════════════════════════════════


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_plot_ns_detection_mixed_profile_all_types():
    import matplotlib.pyplot as plt

    rng = np.random.default_rng(7)
    sites = []
    # clean
    sites += [_site(f"C{i}") for i in range(3)]
    # static shift
    sites += [_site(f"T{i}", rho=100.0 * 10**0.4) for i in range(2)]
    # near-surface (HF scatter)
    for i in range(2):
        fr = np.logspace(-2, 3, 12)
        log_rho = np.full(12, 2.0)
        hf = fr >= 1.0
        log_rho[hf] += rng.normal(0, 1.3, size=hf.sum())
        sites.append(_FakeSite(f"N{i}", _make_z(fr, 10.0**log_rho), fr))
    ax = plot_ns_detection(sites, ns_threshold=1.5, ss_threshold=0.2)
    assert ax is not None
    plt.close("all")


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_plot_ns_detection_existing_axes_reused():
    import matplotlib.pyplot as plt

    _, ax_in = plt.subplots()
    sites = _profile(4)
    ax_out = plot_ns_detection(sites, ax=ax_in)
    assert ax_out is ax_in
    plt.close("all")


# ═══════════════════════════════════════════════════════════════════════════
# _scale_site_Z — direct unit tests (guard / exception branches)
# ═══════════════════════════════════════════════════════════════════════════


def test_scale_site_z_no_z_block_is_noop():
    site = _NoZSite("S1")
    _scale_site_Z(site, 2.0)  # must not raise
    assert site.Z is None


def test_scale_site_z_non_numeric_factor_is_noop():
    site = _site("S1")
    before = site.Z.z.copy()
    _scale_site_Z(site, "not-a-number")
    np.testing.assert_array_equal(site.Z.z, before)


def test_scale_site_z_scales_z_err_when_present():
    site = _site("S1")
    z = site.Z.z
    site.Z.z_err = np.abs(z) * 0.05  # matching-shape error array
    before_err = site.Z.z_err.copy()
    _scale_site_Z(site, 2.0)
    np.testing.assert_allclose(site.Z.z_err, before_err * 2.0)


# ═══════════════════════════════════════════════════════════════════════════
# "no usable Z" branches across per-site loops (_prep_lr_curves,
# detect_near_surface, plot_ss_station_curves, plot_ss_radar,
# ss_comparison_psection)
# ═══════════════════════════════════════════════════════════════════════════


def test_estimate_ss_loess_skips_site_with_no_z():
    sites = [_NoZSite("BAD")] + _shifted_profile(6)
    tbl = estimate_ss_loess(sites, half_window=2)
    assert "BAD" not in set(tbl["station"])
    assert len(tbl) > 0


def test_prep_lr_curves_pband_excludes_one_station_only():
    """One station's whole band falls outside pband -> its `continue`
    branch fires while the rest of the profile still yields rows."""
    sites = _profile(4, n_freq=10, f_lo=-2, f_hi=3)
    outlier = _site("OUT", n=6, f_lo=4, f_hi=5)  # far outside pband below
    tbl = estimate_ss_bilateral(sites + [outlier], half_window=2, pband=(1e-2, 1.0))
    assert "OUT" not in set(tbl["station"])


def test_detect_near_surface_skips_site_with_no_z():
    sites = [_NoZSite("BAD")] + _profile(4)
    df = detect_near_surface(sites)
    assert "BAD" not in set(df["station"])
    assert len(df) == 4


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_plot_ss_station_curves_no_z_blocks():
    import matplotlib.pyplot as plt

    before = [_NoZSite("S1")]
    after = [_NoZSite("S1")]
    ax = plot_ss_station_curves(before, after)
    assert ax is not None
    plt.close("all")


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_plot_ss_radar_no_z():
    import matplotlib.pyplot as plt

    ax = plot_ss_radar([_NoZSite("S1")])
    assert ax is not None
    plt.close("all")


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_ss_comparison_psection_no_z_data():
    import matplotlib.pyplot as plt

    sites = [_NoZSite("S1"), _NoZSite("S2")]
    fig = ss_comparison_psection(sites, method="ama", half_window=1)
    assert fig is not None
    plt.close("all")
