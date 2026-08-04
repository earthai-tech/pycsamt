"""Tests for pycsamt.emtools.tensor"""

from __future__ import annotations

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from pycsamt.emtools.tensor import (
    antisymmetrize,
    balance_offdiag,
    build_phase_tensor_table,
    invert,
    rotate,
    sigma_clip_z,
)

# ─────────────────────────────────────────────────────────────────────────────
# Shared helpers
# ─────────────────────────────────────────────────────────────────────────────


class _FakeZ:
    def __init__(self, z, freq):
        self.z = np.asarray(z, dtype=complex)
        self.freq = np.asarray(freq, dtype=float)


class _FakeSite:
    def __init__(self, station, z, freq):
        self.station = station
        self.Z = _FakeZ(z, freq)
        self.freq = np.asarray(freq, dtype=float)

    def get_section(self, *_, **__):
        return None


def _freqs(n: int = 10, f_lo: float = 1.0, f_hi: float = 1e4) -> np.ndarray:
    return np.logspace(np.log10(f_lo), np.log10(f_hi), n)


def _iso_z(freqs: np.ndarray, rho: float = 100.0) -> np.ndarray:
    amp = np.sqrt(5.0 * freqs * rho)
    z = np.zeros((freqs.size, 2, 2), dtype=complex)
    z[:, 0, 1] = amp * (1 + 1j) / np.sqrt(2)
    z[:, 1, 0] = -amp * (1 + 1j) / np.sqrt(2)
    return z


def _3d_z(freqs: np.ndarray, skew_frac: float = 0.4) -> np.ndarray:
    z = _iso_z(freqs)
    amp = np.abs(z[:, 0, 1])
    z[:, 0, 0] = skew_frac * amp * (0.6 + 0.8j)
    z[:, 1, 0] = -skew_frac * amp * (0.5 + 0.7j)
    return z


def _site(name: str, z=None, n: int = 10) -> _FakeSite:
    fr = _freqs(n)
    if z is None:
        z = _iso_z(fr)
    return _FakeSite(name, z, fr)


# ─────────────────────────────────────────────────────────────────────────────
# build_phase_tensor_table
# ─────────────────────────────────────────────────────────────────────────────


class TestBuildPhaseTensorTable:
    def test_returns_dataframe(self):
        import pandas as pd

        sites = [_site("S00")]
        df = build_phase_tensor_table(sites)
        assert isinstance(df, pd.DataFrame)

    def test_expected_columns(self):
        sites = [_site("S00")]
        df = build_phase_tensor_table(sites)
        for col in ("station", "freq", "period", "beta", "ellipt", "theta"):
            assert col in df.columns

    def test_row_count(self):
        n = 8
        fr = _freqs(n)
        sites = [_FakeSite("S00", _iso_z(fr), fr)]
        df = build_phase_tensor_table(sites)
        assert len(df) == n

    def test_multiple_sites(self):
        sites = [_site(f"S{i:02d}") for i in range(4)]
        df = build_phase_tensor_table(sites)
        assert df["station"].nunique() == 4

    def test_empty_input_empty_df(self):
        import pandas as pd

        df = build_phase_tensor_table([])
        assert isinstance(df, pd.DataFrame)

    def test_period_positive(self):
        sites = [_site("S00")]
        df = build_phase_tensor_table(sites)
        assert (df["period"].dropna() > 0).all()

    def test_skips_non_finite_impedance_rows(self):
        fr = _freqs(5)
        z = _iso_z(fr)
        z[2, 0, 1] = np.nan + 0j
        df = build_phase_tensor_table([_FakeSite("S00", z, fr)])
        assert len(df) == 4
        assert np.isfinite(df["freq"].to_numpy(float)).all()

    def test_ellipt_between_0_1(self):
        sites = [_site("S00")]
        df = build_phase_tensor_table(sites)
        elp = df["ellipt"].dropna()
        assert (elp >= 0.0).all()
        assert (elp <= 1.0).all()


# ─────────────────────────────────────────────────────────────────────────────
# rotate
# ─────────────────────────────────────────────────────────────────────────────


class TestRotate:
    def test_returns_sites(self):
        from pycsamt.site.base import Sites

        sites = [_site("S00")]
        result = rotate(sites, 30.0)
        assert isinstance(result, Sites)

    def test_site_count_preserved(self):
        sites = [_site(f"S{i:02d}") for i in range(3)]
        result = rotate(sites, 45.0)
        assert sum(1 for _ in result) == 3

    def test_zero_rotation_unchanged(self):
        """Rotating by 0° must not change |Z_xy|."""
        fr = _freqs(6)
        z0 = _iso_z(fr)
        sites = [_FakeSite("S00", z0.copy(), fr)]
        result = rotate(sites, 0.0)
        from pycsamt.emtools._core import (
            _get_z_block,
            _iter_items,
        )

        ed = next(_iter_items(result))
        _, z1, _ = _get_z_block(ed)
        if z1 is not None:
            np.testing.assert_allclose(
                np.abs(z1[:, 0, 1]), np.abs(z0[:, 0, 1]), rtol=1e-6
            )

    def test_rotation_360_returns_to_start(self):
        fr = _freqs(6)
        z0 = _iso_z(fr)
        sites = [_FakeSite("S00", z0.copy(), fr)]
        result = rotate(sites, 360.0)
        from pycsamt.emtools._core import (
            _get_z_block,
            _iter_items,
        )

        ed = next(_iter_items(result))
        _, z1, _ = _get_z_block(ed)
        if z1 is not None:
            np.testing.assert_allclose(
                np.abs(z1[:, 0, 1]), np.abs(z0[:, 0, 1]), rtol=1e-5
            )

    def test_empty_sites(self):
        from pycsamt.site.base import Sites

        result = rotate([], 30.0)
        assert isinstance(result, Sites)


# ─────────────────────────────────────────────────────────────────────────────
# antisymmetrize
# ─────────────────────────────────────────────────────────────────────────────


class TestAntisymmetrize:
    def test_returns_sites(self):
        from pycsamt.site.base import Sites

        sites = [_site("S00")]
        result = antisymmetrize(sites)
        assert isinstance(result, Sites)

    def test_offdiag_antisymmetry_enforced(self):
        """After antisymmetrize, |Zxy| == |Zyx| (average mode)."""
        fr = _freqs(6)
        z0 = _3d_z(fr)  # asymmetric off-diagonal
        sites = [_FakeSite("S00", z0.copy(), fr)]
        result = antisymmetrize(sites, how="rms")
        from pycsamt.emtools._core import (
            _get_z_block,
            _iter_items,
        )

        ed = next(_iter_items(result))
        _, z1, _ = _get_z_block(ed)
        if z1 is not None:
            np.testing.assert_allclose(
                np.abs(z1[:, 0, 1]), np.abs(z1[:, 1, 0]), rtol=1e-6
            )

    def test_empty_sites(self):
        from pycsamt.site.base import Sites

        result = antisymmetrize([])
        assert isinstance(result, Sites)


# ─────────────────────────────────────────────────────────────────────────────
# invert
# ─────────────────────────────────────────────────────────────────────────────


class TestInvert:
    def test_returns_sites(self):
        from pycsamt.site.base import Sites

        sites = [_site("S00")]
        result = invert(sites)
        assert isinstance(result, Sites)

    def test_double_invert_recovers_original(self):
        """Z^{-1-1} ≈ Z (up to numerical tolerance)."""
        fr = _freqs(6)
        z0 = _iso_z(fr)
        sites = [_FakeSite("S00", z0.copy(), fr)]
        twice = invert(invert(sites))
        from pycsamt.emtools._core import (
            _get_z_block,
            _iter_items,
        )

        ed = next(_iter_items(twice))
        _, z2, _ = _get_z_block(ed)
        if z2 is not None:
            np.testing.assert_allclose(
                np.abs(z2[:, 0, 1]), np.abs(z0[:, 0, 1]), rtol=1e-5
            )

    def test_empty_sites(self):
        from pycsamt.site.base import Sites

        result = invert([])
        assert isinstance(result, Sites)


# ─────────────────────────────────────────────────────────────────────────────
# sigma_clip_z
# ─────────────────────────────────────────────────────────────────────────────


class TestSigmaClipZ:
    def test_returns_sites(self):
        from pycsamt.site.base import Sites

        sites = [_site("S00")]
        result = sigma_clip_z(sites, sigma=3.0)
        assert isinstance(result, Sites)

    def test_preserves_clean_data(self):
        """For a smooth synthetic Z, sigma clipping should not nan everything."""
        fr = _freqs(10)
        z0 = _iso_z(fr)
        sites = [_FakeSite("S00", z0.copy(), fr)]
        result = sigma_clip_z(sites, sigma=3.0)
        from pycsamt.emtools._core import (
            _get_z_block,
            _iter_items,
        )

        ed = next(_iter_items(result))
        _, z1, _ = _get_z_block(ed)
        if z1 is not None:
            # at least half the rows should remain finite
            finite_frac = np.isfinite(z1[:, 0, 1]).mean()
            assert finite_frac > 0.5


# ─────────────────────────────────────────────────────────────────────────────
# balance_offdiag
# ─────────────────────────────────────────────────────────────────────────────


class TestBalanceOffdiag:
    def test_returns_sites(self):
        from pycsamt.site.base import Sites

        sites = [_site("S00")]
        result = balance_offdiag(sites)
        assert isinstance(result, Sites)

    def test_offdiag_equal_magnitude_after_balance(self):
        fr = _freqs(6)
        z0 = _3d_z(fr)  # |Zxy| != |Zyx|
        sites = [_FakeSite("S00", z0.copy(), fr)]
        result = balance_offdiag(sites, mode="avgabs")
        from pycsamt.emtools._core import (
            _get_z_block,
            _iter_items,
        )

        ed = next(_iter_items(result))
        _, z1, _ = _get_z_block(ed)
        if z1 is not None:
            np.testing.assert_allclose(
                np.abs(z1[:, 0, 1]), np.abs(z1[:, 1, 0]), rtol=1e-6
            )

    def test_empty_sites(self):
        from pycsamt.site.base import Sites

        result = balance_offdiag([])
        assert isinstance(result, Sites)


# ─────────────────────────────────────────────────────────────────────────────
# plot_strike_director_field
# ─────────────────────────────────────────────────────────────────────────────


class TestPlotStrikeDirectorField:
    def _sites(self, n_sta: int = 5, n_freq: int = 12):
        return [_site(f"S{i:02d}", _3d_z(_freqs(n_freq)), n_freq) for i in range(n_sta)]

    def test_returns_axes(self):
        from pycsamt.emtools.tensor import (
            plot_strike_director_field,
        )

        ax = plot_strike_director_field(self._sites())
        assert isinstance(ax, plt.Axes)
        plt.close("all")

    def test_exported_from_package(self):
        from pycsamt.emtools import plot_strike_director_field

        assert callable(plot_strike_director_field)

    def test_empty_input_no_crash(self):
        from pycsamt.emtools.tensor import (
            plot_strike_director_field,
        )

        ax = plot_strike_director_field([])
        assert isinstance(ax, plt.Axes)
        plt.close("all")

    def test_draws_one_director_per_cell(self):
        from pycsamt.emtools.tensor import (
            plot_strike_director_field,
        )

        n_sta, n_freq = 4, 10
        ax = plot_strike_director_field(self._sites(n_sta, n_freq), streamlines=False)
        # exactly one quiver collection holding n_sta * n_freq directors
        quivers = [c for c in ax.collections if c.__class__.__name__ == "Quiver"]
        assert len(quivers) == 1
        assert quivers[0].N == n_sta * n_freq
        plt.close("all")

    def test_streamlines_toggle(self):
        from pycsamt.emtools.tensor import (
            plot_strike_director_field,
        )

        ax = plot_strike_director_field(self._sites(), streamlines=True)
        assert isinstance(ax, plt.Axes)
        plt.close("all")

    def test_color_and_length_options(self):
        from pycsamt.emtools.tensor import (
            plot_strike_director_field,
        )

        ax = plot_strike_director_field(
            self._sites(),
            color_by="ellipt",
            length_by=None,
            streamlines=False,
            show_legend=False,
        )
        assert isinstance(ax, plt.Axes)
        plt.close("all")

    def test_accepts_external_axes(self):
        from pycsamt.emtools.tensor import (
            plot_strike_director_field,
        )

        fig, ax0 = plt.subplots()
        ax = plot_strike_director_field(self._sites(), ax=ax0, streamlines=False)
        assert ax is ax0
        plt.close("all")


def test_phase_tensor_psection_uses_dynamic_frame_and_skew_limits():
    """Legends must not create empty period bands or fix beta to +/-3."""
    from pycsamt.emtools.tensor import plot_phase_tensor_psection

    fr = _freqs(12, f_lo=1e-4, f_hi=1e-1)
    sites = [
        _FakeSite("S00", _3d_z(fr, skew_frac=0.7), fr),
        _FakeSite("S01", _3d_z(fr, skew_frac=1.1), fr),
    ]
    ax = plot_phase_tensor_psection(
        sites, c_by="beta", period_up=False, recursive=False
    )
    data_y = np.log10(1.0 / fr)
    ylim = ax.get_ylim()
    assert min(ylim) >= data_y.min() - 0.5
    assert max(ylim) <= data_y.max() + 0.5

    cbar = ax.figure.axes[-1]
    mappables = [item for item in cbar.collections if item.get_clim() != (None, None)]
    assert mappables
    vmin, vmax = mappables[0].get_clim()
    assert not np.isclose(vmax, 3.0)
    assert not np.isclose(vmin, -3.0)
    assert np.isclose(abs(vmin), abs(vmax))
    plt.close("all")


def test_phase_tensor_psection_shape_mode_and_explicit_skew_clip():
    """Shape mode keeps cells visible while allowing a publication beta scale."""
    from matplotlib.patches import Ellipse
    from pycsamt.emtools.tensor import plot_phase_tensor_psection

    fr = _freqs(8, f_lo=1e-3, f_hi=1.0)
    sites = [_FakeSite("S00", _3d_z(fr, skew_frac=1.0), fr)]
    ax = plot_phase_tensor_psection(
        sites,
        c_by="beta",
        normalise_by="shape",
        min_aspect=0.12,
        clim=(-3.0, 3.0),
        recursive=False,
    )
    ellipses = [patch for patch in ax.patches if isinstance(patch, Ellipse)]
    assert len(ellipses) >= len(fr)
    renderer = ax.figure.canvas.get_renderer()
    for ellipse in ellipses[: len(fr)]:
        bounds = ellipse.get_window_extent(renderer).bounds
        width, height = bounds[2], bounds[3]
        assert min(width, height) / max(width, height) >= 0.11

    cbar = ax.figure.axes[-1]
    mappables = [item for item in cbar.collections if item.get_clim() != (None, None)]
    assert mappables
    assert np.allclose(mappables[0].get_clim(), (-3.0, 3.0))
    plt.close("all")


def test_phase_tensor_psection_segmented_colors_and_artist_kwargs():
    """Discrete skew classes and ellipse/colorbar kwargs are user-controlled."""
    from matplotlib.colors import BoundaryNorm, to_rgba
    from matplotlib.patches import Ellipse
    from pycsamt.emtools.tensor import plot_phase_tensor_psection

    fr = _freqs(6, f_lo=1e-3, f_hi=1.0)
    ax = plot_phase_tensor_psection(
        [_FakeSite("S00", _3d_z(fr, skew_frac=1.0), fr)],
        c_by="beta",
        normalise_by="shape",
        color_mode="segmented",
        segment_bounds=(-3.0, 3.0),
        segment_colors=("navy", "white", "firebrick"),
        ellipse_kws={"edgecolor": "lime", "linewidth": 0.7},
        cb_kws={"size": "4%", "pad": 0.1},
        recursive=False,
    )
    ellipses = [patch for patch in ax.patches if isinstance(patch, Ellipse)]
    assert ellipses
    assert np.allclose(ellipses[0].get_edgecolor()[:3], to_rgba("lime")[:3])
    assert any(
        isinstance(collection.norm, BoundaryNorm)
        for collection in ax.figure.axes[-1].collections
    )
    assert [tick.get_text() for tick in ax.figure.axes[-1].get_yticklabels()] == [
        "< -3", "-3 to 3", "> 3"
    ]
    plt.close("all")
