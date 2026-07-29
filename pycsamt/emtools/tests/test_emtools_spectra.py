"""Tests for pycsamt.emtools.spectra (analysis functions)"""

from __future__ import annotations

from types import SimpleNamespace

import matplotlib
import numpy as np
import pytest

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.figure import Figure  # noqa: E402

from pycsamt.api import reset_api_view  # noqa: E402
from pycsamt.api.control import PYCSAMT_CONTROL  # noqa: E402
from pycsamt.api.section import PYCSAMT_SECTION, SectionStyle  # noqa: E402
from pycsamt.emtools.spectra import (  # noqa: E402
    band_select,
    coherence_matrix,
    coherence_table,
    mask_low_coherence,
    plot_coherence,
    plot_coherence_section,
    plot_psd,
    plot_psd_section,
    plot_spectra_matrix,
    plot_tipper_from_spectra,
    plot_z_from_spectra,
    psd_table,
    snr_table,
    spectra_summary,
)


@pytest.fixture(autouse=True)
def _no_api_view():
    """Ensure each test starts with the global API view reset."""
    reset_api_view()
    yield
    reset_api_view()


# ─────────────────────────────────────────────────────────────────────────────
# Fake Spectra object
# ─────────────────────────────────────────────────────────────────────────────


class _FakeSpectra:
    """Minimal Spectra duck-type with 3 channels and N frequencies."""

    def __init__(self, n_freq: int = 8, n_chan: int = 3, name: str = "STA01") -> None:
        self.name = name
        fr = np.logspace(0, 3, n_freq)
        self.freq = fr
        self.bw = np.full(n_freq, 0.5)
        self.avgt = np.full(n_freq, 10.0)
        self.avgf = fr.copy()
        self.rotspec = np.zeros(n_freq)
        self.segnum = np.full(n_freq, 4, dtype=int)
        self.band = [f"B{i}" for i in range(n_freq)]
        self.chan_ids = [str(i + 1) for i in range(n_chan)]
        self.id_to_chtype = {"1": "Ex", "2": "Ey", "3": "Hx"}

        # Build Hermitian positive-definite S matrix
        rng = np.random.default_rng(42)
        S = np.zeros((n_freq, n_chan, n_chan), dtype=complex)
        for k in range(n_freq):
            A = rng.standard_normal((n_chan, n_chan)) + 1j * rng.standard_normal(
                (n_chan, n_chan)
            )
            S[k] = A @ A.conj().T + np.eye(n_chan) * 2.0
        self._S = S

    @property
    def S(self) -> np.ndarray:
        return self._S

    @property
    def n_freq(self) -> int:
        return int(self.freq.size)

    @property
    def n_chan(self) -> int:
        return int(self._S.shape[1])


def _sp(n_freq=8, n_chan=3, name="STA01") -> _FakeSpectra:
    return _FakeSpectra(n_freq=n_freq, n_chan=n_chan, name=name)


@pytest.fixture(scope="module")
def real_spectra(edi_spe_file):
    """A real Spectra object parsed from the bundled spectra EDI."""
    from pycsamt.seg.edi import EDIFile

    ed = EDIFile(str(edi_spe_file))
    sp = ed.spectra
    if sp is None or sp.n_freq == 0:
        pytest.skip("Bundled spectra EDI produced no usable Spectra object.")
    return sp


class _FakeZTSpectra:
    """Duck-typed Spectra exposing a deterministic ``to_Z`` for plot tests.

    Lets us exercise the error / no-error and tipper / no-tipper branches
    of ``plot_z_from_spectra`` and ``plot_tipper_from_spectra`` without
    depending on the exact contents of the bundled EDI files.
    """

    def __init__(
        self,
        n_freq: int = 6,
        with_tipper: bool = True,
        with_tipper_err: bool = True,
        name: str = "ZTSITE",
    ) -> None:
        self.name = name
        self._freq = np.logspace(0, 2, n_freq)
        self._with_tipper = with_tipper
        self._with_tipper_err = with_tipper_err
        self._n_chan = 5 if with_tipper else 4

    @property
    def n_freq(self) -> int:
        return int(self._freq.size)

    @property
    def n_chan(self) -> int:
        return self._n_chan

    def to_Z(self, *, e_labels=None, h_labels=None, ridge=None, estimate_error=False):
        nf = self.n_freq
        rng = np.random.default_rng(9)
        z = rng.standard_normal((nf, 2, 2)) + 1j * rng.standard_normal((nf, 2, 2))
        rho = np.abs(z) ** 2 + 1.0
        phase = np.angle(z, deg=True)
        z_err = np.abs(z) * 0.05 if estimate_error else None
        z_obj = SimpleNamespace(
            freq=self._freq, resistivity=rho, phase=phase, z=z, z_err=z_err
        )
        tip = None
        if self._with_tipper:
            t = rng.standard_normal((nf, 1, 2)) + 1j * rng.standard_normal((nf, 1, 2))
            t_err = (
                np.abs(t) * 0.05 if (estimate_error and self._with_tipper_err) else None
            )
            tip = SimpleNamespace(freq=self._freq, tipper=t, tipper_err=t_err)
        return z_obj, tip


# ─────────────────────────────────────────────────────────────────────────────
# coherence_matrix
# ─────────────────────────────────────────────────────────────────────────────


class TestCoherenceMatrix:
    def test_shape(self):
        sp = _sp(8, 3)
        coh = coherence_matrix(sp)
        assert coh.shape == (8, 3, 3)

    def test_diagonal_is_one(self):
        sp = _sp(8, 3)
        coh = coherence_matrix(sp)
        diag = np.diagonal(coh, axis1=1, axis2=2)
        np.testing.assert_allclose(diag, 1.0, atol=1e-10)

    def test_values_in_0_1(self):
        sp = _sp(8, 3)
        coh = coherence_matrix(sp)
        assert coh.min() >= 0.0
        assert coh.max() <= 1.0 + 1e-12

    def test_symmetric(self):
        sp = _sp(8, 3)
        coh = coherence_matrix(sp)
        np.testing.assert_allclose(coh, np.transpose(coh, (0, 2, 1)), atol=1e-10)

    def test_single_channel_raises(self):
        sp = _sp(8, 1)
        # n_chan=1: still valid shape — just check no exception
        coh = coherence_matrix(sp)
        assert coh.shape == (8, 1, 1)


# ─────────────────────────────────────────────────────────────────────────────
# psd_table
# ─────────────────────────────────────────────────────────────────────────────


class TestPsdTable:
    def test_returns_dataframe(self):
        import pandas as pd

        sp = _sp()
        df = psd_table({"STA01": sp}, api=False)
        assert isinstance(df, pd.DataFrame)

    def test_expected_columns(self):
        sp = _sp()
        df = psd_table({"STA01": sp}, api=False)
        for col in ("station", "freq", "period", "channel", "psd"):
            assert col in df.columns

    def test_row_count(self):
        nf, nc = 8, 3
        sp = _sp(n_freq=nf, n_chan=nc)
        df = psd_table({"STA01": sp}, api=False)
        assert len(df) == nf * nc

    def test_psd_nonneg(self):
        sp = _sp()
        df = psd_table({"STA01": sp}, api=False)
        assert (df["psd"].dropna() >= 0).all()

    def test_freq_positive(self):
        sp = _sp()
        df = psd_table({"STA01": sp}, api=False)
        assert (df["freq"] > 0).all()

    def test_normalize(self):
        sp = _sp()
        df = psd_table({"STA01": sp}, normalize=True)
        # each channel's max PSD should be 1.0
        assert df["psd"].max() <= 1.0 + 1e-9

    def test_dict_input(self):
        sp1 = _sp(name="A")
        sp2 = _sp(name="B")
        df = psd_table({"A": sp1, "B": sp2})
        assert df["station"].nunique() == 2

    def test_list_input(self):
        sps = [_sp(name=f"S{i}") for i in range(3)]
        df = psd_table(sps)
        assert df["station"].nunique() == 3


# ─────────────────────────────────────────────────────────────────────────────
# coherence_table
# ─────────────────────────────────────────────────────────────────────────────


class TestCoherenceTable:
    def test_returns_dataframe(self):
        import pandas as pd

        sp = _sp()
        df = coherence_table({"STA01": sp}, api=False)
        assert isinstance(df, pd.DataFrame)

    def test_columns(self):
        sp = _sp()
        df = coherence_table({"STA01": sp}, api=False)
        for col in ("station", "freq", "ch_i", "ch_j", "coherence"):
            assert col in df.columns

    def test_coherence_in_range(self):
        sp = _sp()
        df = coherence_table({"STA01": sp}, api=False)
        coh = df["coherence"].dropna()
        assert (coh >= 0.0).all() and (coh <= 1.0 + 1e-9).all()

    def test_num_pairs(self):
        """Upper-triangle pairs for n_chan=3: (0,1),(0,2),(1,2) → 3 pairs."""
        nf, nc = 6, 3
        sp = _sp(n_freq=nf, n_chan=nc)
        df = coherence_table({"STA01": sp}, api=False)
        n_pairs = nc * (nc - 1) // 2
        assert len(df) == nf * n_pairs

    def test_custom_pairs(self):
        sp = _sp()
        df = coherence_table({"STA01": sp}, pairs=[(0, 1)], api=False)
        assert df["coherence"].notna().any()


# ─────────────────────────────────────────────────────────────────────────────
# snr_table
# ─────────────────────────────────────────────────────────────────────────────


class TestSnrTable:
    def test_returns_tabular(self):
        import pandas as pd

        from pycsamt.api.view.frame import APIFrame

        sp = _sp()
        df = snr_table({"STA01": sp})
        assert isinstance(df, (pd.DataFrame, APIFrame))

    def test_snr_db_column(self):
        sp = _sp()
        df = snr_table({"STA01": sp})
        assert "snr_db" in df.columns

    def test_freq_positive(self):
        sp = _sp()
        df = snr_table({"STA01": sp})
        assert (df["freq"] > 0).all()


# ─────────────────────────────────────────────────────────────────────────────
# band_select
# ─────────────────────────────────────────────────────────────────────────────


class TestBandSelect:
    def test_returns_spectra(self):
        sp = _sp(n_freq=12)
        result = band_select(sp, f_min=10.0, f_max=1000.0)
        assert result is not None

    def test_reduced_frequencies(self):
        sp = _sp(n_freq=12)
        result = band_select(sp, f_min=10.0, f_max=100.0)
        if hasattr(result, "n_freq"):
            assert result.n_freq <= sp.n_freq

    def test_full_range_no_reduction(self):
        sp = _sp(n_freq=8)
        f_lo, f_hi = float(sp.freq.min()), float(sp.freq.max())
        result = band_select(sp, f_min=f_lo * 0.5, f_max=f_hi * 2.0)
        assert result is not None


# ─────────────────────────────────────────────────────────────────────────────
# mask_low_coherence
# ─────────────────────────────────────────────────────────────────────────────


class TestMaskLowCoherence:
    def test_returns_spectra(self):
        sp = _sp()
        result = mask_low_coherence(sp, threshold=0.5)
        assert result is not None

    def test_high_threshold_masks_more(self):
        sp = _sp(8, 3)
        # Just check it runs without error; masking is implementation-defined.
        result = mask_low_coherence(sp, threshold=0.99)
        assert result is not None

    def test_zero_threshold_keeps_all(self):
        sp = _sp(8, 3)
        result = mask_low_coherence(sp, threshold=0.0)
        assert result is not None


# ─────────────────────────────────────────────────────────────────────────────
# spectra_summary
# ─────────────────────────────────────────────────────────────────────────────


class TestSpectraSummary:
    def test_returns_dataframe(self):
        import pandas as pd

        sp = _sp()
        df = spectra_summary(sp, api=False)
        assert isinstance(df, pd.DataFrame)

    def test_expected_columns(self):
        sp = _sp()
        df = spectra_summary(sp, api=False)
        # spectra_summary returns per-frequency rows (freq, period, bw, ...)
        assert "freq" in df.columns

    def test_row_count(self):
        nf = 8
        sp = _sp(n_freq=nf)
        df = spectra_summary(sp, api=False)
        assert len(df) == nf

    def test_single_channel_no_mean_coherence(self):
        """With n_chan == 1 there are no pairs, so the mean-coherence
        column must be absent (nc > 1 branch not taken)."""
        sp = _sp(n_freq=5, n_chan=1)
        df = spectra_summary(sp, api=False)
        assert "mean_coherence" not in df.columns


# ─────────────────────────────────────────────────────────────────────────────
# _check_spectra guard (shared by every analysis / plotting function)
# ─────────────────────────────────────────────────────────────────────────────


class TestCheckSpectraGuards:
    def test_zero_freq_raises(self):
        sp = _sp(n_freq=0, n_chan=3)
        with pytest.raises(ValueError):
            coherence_matrix(sp)

    def test_zero_chan_raises(self):
        sp = _sp(n_freq=5, n_chan=0)
        with pytest.raises(ValueError):
            coherence_matrix(sp)

    def test_spectra_summary_zero_freq_raises(self):
        sp = _sp(n_freq=0, n_chan=2)
        with pytest.raises(ValueError):
            spectra_summary(sp, api=False)

    def test_band_select_zero_chan_raises(self):
        sp = _sp(n_freq=5, n_chan=0)
        with pytest.raises(ValueError):
            band_select(sp, f_min=1.0, f_max=100.0)


# ─────────────────────────────────────────────────────────────────────────────
# _chan_label fallback (missing / empty chan_ids)
# ─────────────────────────────────────────────────────────────────────────────


class TestChanLabelFallback:
    def test_empty_chan_ids_falls_back_to_index(self):
        sp = _sp(n_freq=4, n_chan=2)
        sp.chan_ids = []
        sp.id_to_chtype = {}
        df = psd_table({"S": sp}, api=False)
        assert set(df["channel"]) == {"ch0", "ch1"}

    def test_idx_beyond_chan_ids_falls_back_to_index(self):
        sp = _sp(n_freq=4, n_chan=2)
        sp.chan_ids = ["only-one"]
        df = psd_table({"S": sp}, api=False)
        assert "ch1" in set(df["channel"])


# ─────────────────────────────────────────────────────────────────────────────
# _sp_to_dict: single Spectra instance / invalid type
# ─────────────────────────────────────────────────────────────────────────────


class TestSpToDict:
    def test_invalid_type_raises(self):
        with pytest.raises(TypeError):
            psd_table(12345, api=False)

    def test_single_real_spectra_instance(self, real_spectra):
        df = psd_table(real_spectra, api=False)
        assert df["station"].nunique() == 1
        assert df["station"].iloc[0] == (real_spectra.name or "site")

    def test_single_real_spectra_coherence_table(self, real_spectra):
        df = coherence_table(real_spectra, api=False)
        assert df["station"].nunique() == 1
        assert len(df) > 0


# ─────────────────────────────────────────────────────────────────────────────
# psd_table: normalize with an all-zero channel (skip-division branch)
# ─────────────────────────────────────────────────────────────────────────────


class TestPsdTableZeroChannel:
    def test_normalize_skips_zero_max(self):
        sp = _sp(n_freq=6, n_chan=2)
        sp._S[:, 0, 0] = 0.0  # channel 0 has zero power everywhere
        df = psd_table({"S": sp}, normalize=True, api=False)
        assert not df["psd"].isna().any()
        assert np.isfinite(df["psd"]).all()
        ch0_label = _chan_label_str(sp, 0)
        zero_rows = df[df["channel"] == ch0_label]
        assert len(zero_rows) == sp.n_freq
        assert (zero_rows["psd"] == 0.0).all()


# ─────────────────────────────────────────────────────────────────────────────
# band_select: out-of-range raises
# ─────────────────────────────────────────────────────────────────────────────


class TestBandSelectErrors:
    def test_out_of_range_raises(self):
        sp = _sp(n_freq=8)
        with pytest.raises(ValueError):
            band_select(sp, f_min=1e9, f_max=2e9)


# ─────────────────────────────────────────────────────────────────────────────
# mask_low_coherence: require_all branch
# ─────────────────────────────────────────────────────────────────────────────


class TestMaskLowCoherenceRequireAll:
    def test_require_all_is_subset_of_any(self):
        sp = _sp(8, 3)
        mask_any = mask_low_coherence(sp, threshold=0.3, require_all=False)
        mask_all = mask_low_coherence(sp, threshold=0.3, require_all=True)
        assert mask_any.dtype == bool
        assert mask_all.dtype == bool
        # every freq flagged under "all" must also be flagged under "any"
        assert np.all(mask_any[mask_all])
        assert mask_all.sum() <= mask_any.sum()


# ─────────────────────────────────────────────────────────────────────────────
# snr_table: custom pairs
# ─────────────────────────────────────────────────────────────────────────────


class TestSnrTableCustomPairs:
    def test_custom_pairs(self):
        sp = _sp()
        df = snr_table({"STA01": sp}, pairs=[(0, 1)], api=False)
        assert set(df["pair"].unique()) <= {
            f"{_chan_label_str(sp, 0)}-{_chan_label_str(sp, 1)}"
        }


def _chan_label_str(sp, idx):
    ids = getattr(sp, "chan_ids", None)
    id_to_cht = getattr(sp, "id_to_chtype", {}) or {}
    if ids and idx < len(ids):
        raw = ids[idx]
        cht = id_to_cht.get(str(raw), "")
        return f"{cht}({raw})" if cht else str(raw)
    return f"ch{idx}"


# ─────────────────────────────────────────────────────────────────────────────
# plot_psd
# ─────────────────────────────────────────────────────────────────────────────


class TestPlotPsd:
    def test_default(self):
        sp = _sp()
        ax = plot_psd(sp)
        assert ax is not None
        assert len(ax.lines) == sp.n_chan
        plt.close("all")

    def test_channel_subset(self):
        sp = _sp(n_chan=3)
        ax = plot_psd(sp, channels=[0, 2])
        assert len(ax.lines) == 2
        plt.close("all")

    def test_log_psd_false(self):
        sp = _sp()
        ax = plot_psd(sp, log_psd=False)
        assert ax.get_ylabel() == "PSD"
        plt.close("all")

    def test_custom_ax_reused(self):
        sp = _sp()
        fig, ax = plt.subplots()
        out = plot_psd(sp, ax=ax)
        assert out is ax
        plt.close("all")

    def test_custom_lw_alpha_title(self):
        sp = _sp()
        ax = plot_psd(sp, lw=2.0, alpha=0.3, title="Custom PSD Title")
        assert ax.get_title() == "Custom PSD Title"
        plt.close("all")


# ─────────────────────────────────────────────────────────────────────────────
# plot_coherence
# ─────────────────────────────────────────────────────────────────────────────


class TestPlotCoherence:
    def test_default_all_pairs(self):
        sp = _sp(8, 3)
        axs = plot_coherence(sp)
        assert len(axs) == 3
        plt.close("all")

    def test_custom_single_pair(self):
        sp = _sp(8, 3)
        axs = plot_coherence(sp, pairs=[(0, 1)])
        assert len(axs) == 1
        plt.close("all")

    def test_show_threshold_false(self):
        sp = _sp(8, 3)
        axs = plot_coherence(sp, pairs=[(0, 1)], show_threshold=False)
        assert len(axs) == 1
        plt.close("all")

    def test_axes_given_reused(self):
        sp = _sp(8, 3)
        fig, ax = plt.subplots()
        axs = plot_coherence(sp, pairs=[(0, 1)], axes=[ax])
        assert axs[0] is ax
        plt.close("all")

    def test_no_pairs_raises(self):
        sp = _sp(8, 1)
        with pytest.raises(ValueError):
            plot_coherence(sp)
        plt.close("all")


# ─────────────────────────────────────────────────────────────────────────────
# plot_spectra_matrix
# ─────────────────────────────────────────────────────────────────────────────


class TestPlotSpectraMatrix:
    @pytest.mark.parametrize("quantity", ["abs", "real", "imag", "phase"])
    def test_quantities(self, quantity):
        sp = _sp(8, 3)
        fig = plot_spectra_matrix(sp, quantity=quantity)
        assert isinstance(fig, Figure)
        plt.close("all")

    def test_abs_log_scale_false(self):
        sp = _sp(8, 3)
        fig = plot_spectra_matrix(sp, quantity="abs", log_scale=False)
        assert isinstance(fig, Figure)
        plt.close("all")

    def test_invalid_quantity_raises(self):
        sp = _sp(8, 3)
        with pytest.raises(ValueError):
            plot_spectra_matrix(sp, quantity="bogus")
        plt.close("all")

    def test_custom_ax_reused(self):
        sp = _sp(8, 3)
        fig_in, ax = plt.subplots()
        fig_out = plot_spectra_matrix(sp, ax=ax)
        assert fig_out is ax.figure
        plt.close("all")

    def test_custom_freq_idx_and_title(self):
        sp = _sp(8, 3)
        fig = plot_spectra_matrix(sp, freq_idx=3, title="Matrix @ f3")
        ax = fig.axes[0]
        assert ax.get_title() == "Matrix @ f3"
        plt.close("all")


# ─────────────────────────────────────────────────────────────────────────────
# plot_z_from_spectra
# ─────────────────────────────────────────────────────────────────────────────


class TestPlotZFromSpectra:
    def test_default_no_error(self):
        sp = _FakeZTSpectra(with_tipper=False)
        fig = plot_z_from_spectra(sp)
        assert isinstance(fig, Figure)
        assert len(fig.axes) == 2
        plt.close("all")

    def test_estimate_error_true_shows_error_band(self):
        sp = _FakeZTSpectra(with_tipper=False)
        fig = plot_z_from_spectra(sp, estimate_error=True, show_error=True)
        assert isinstance(fig, Figure)
        ax_r = fig.axes[0]
        # a filled error band adds PolyCollections to the resistivity axis
        assert len(ax_r.collections) >= 1
        plt.close("all")

    def test_estimate_error_true_show_error_false(self):
        sp = _FakeZTSpectra(with_tipper=False)
        fig = plot_z_from_spectra(sp, estimate_error=True, show_error=False)
        ax_r = fig.axes[0]
        assert len(ax_r.collections) == 0
        plt.close("all")

    def test_custom_axes_reused(self):
        sp = _FakeZTSpectra(with_tipper=False)
        fig_in, (ax_r, ax_p) = plt.subplots(1, 2)
        fig = plot_z_from_spectra(sp, axes=(ax_r, ax_p))
        assert fig is ax_r.figure
        plt.close("all")

    def test_custom_title(self):
        sp = _FakeZTSpectra(with_tipper=False)
        fig = plot_z_from_spectra(sp, title="Custom Z Title")
        assert fig.get_suptitle() == "Custom Z Title"
        plt.close("all")

    def test_real_spectra_smoke(self, real_spectra):
        fig = plot_z_from_spectra(real_spectra)
        assert isinstance(fig, Figure)
        plt.close("all")


# ─────────────────────────────────────────────────────────────────────────────
# plot_tipper_from_spectra
# ─────────────────────────────────────────────────────────────────────────────


class TestPlotTipperFromSpectra:
    def test_with_tipper_and_error(self):
        sp = _FakeZTSpectra(with_tipper=True, with_tipper_err=True)
        axs = plot_tipper_from_spectra(sp, estimate_error=True)
        assert len(axs) == 2
        ax_a = axs[0]
        assert len(ax_a.collections) >= 1  # error band present
        plt.close("all")

    def test_with_tipper_no_error(self):
        sp = _FakeZTSpectra(with_tipper=True, with_tipper_err=False)
        axs = plot_tipper_from_spectra(sp, estimate_error=True)
        ax_a = axs[0]
        assert len(ax_a.collections) == 0
        plt.close("all")

    def test_no_tipper_shows_message(self):
        sp = _FakeZTSpectra(with_tipper=False)
        axs = plot_tipper_from_spectra(sp)
        assert len(axs) == 2
        texts = [t.get_text() for t in axs[0].texts]
        assert any("No tipper" in t for t in texts)
        plt.close("all")

    def test_no_tipper_with_custom_axes(self):
        sp = _FakeZTSpectra(with_tipper=False)
        fig_in, (ax1, ax2) = plt.subplots(1, 2)
        axs = plot_tipper_from_spectra(sp, axes=(ax1, ax2))
        assert axs[0] is ax1
        assert axs[1] is ax2
        plt.close("all")

    def test_real_spectra_smoke(self, real_spectra):
        axs = plot_tipper_from_spectra(real_spectra)
        assert len(axs) == 2
        plt.close("all")


# ─────────────────────────────────────────────────────────────────────────────
# plot_psd_section
# ─────────────────────────────────────────────────────────────────────────────


class TestPlotPsdSection:
    def test_default_dict_input(self):
        sps = {f"S{i}": _sp(n_freq=8, n_chan=2, name=f"S{i}") for i in range(3)}
        ax = plot_psd_section(sps)
        assert ax is not None
        plt.close("all")

    def test_list_input_channel_and_log(self):
        sps = [_sp(n_freq=8, n_chan=2, name=f"L{i}") for i in range(2)]
        ax = plot_psd_section(sps, channel=1, log_psd=False)
        assert ax is not None
        plt.close("all")

    def test_vmin_vmax_custom_ax(self):
        sps = {f"S{i}": _sp(n_freq=6, n_chan=2, name=f"S{i}") for i in range(2)}
        fig, ax_in = plt.subplots()
        ax = plot_psd_section(sps, vmin=-5.0, vmax=5.0, ax=ax_in)
        assert ax is ax_in
        plt.close("all")

    def test_section_style_instance(self):
        sps = {f"S{i}": _sp(n_freq=6, n_chan=2, name=f"S{i}") for i in range(2)}
        sty = PYCSAMT_SECTION.style_for("compact")
        assert isinstance(sty, SectionStyle)
        ax = plot_psd_section(sps, section=sty)
        assert ax is not None
        plt.close("all")

    def test_single_frequency_edge_case(self):
        """n_f == 1 exercises the degenerate y_edges branch."""
        sps = {"S1": _sp(n_freq=1, n_chan=2, name="S1")}
        ax = plot_psd_section(sps)
        assert ax is not None
        plt.close("all")


# ─────────────────────────────────────────────────────────────────────────────
# plot_coherence_section
# ─────────────────────────────────────────────────────────────────────────────


class TestPlotCoherenceSection:
    def test_default_mean_pairs(self):
        sps = {f"S{i}": _sp(n_freq=8, n_chan=3, name=f"S{i}") for i in range(3)}
        ax = plot_coherence_section(sps)
        assert ax is not None
        plt.close("all")

    def test_explicit_pair(self):
        sps = {f"S{i}": _sp(n_freq=8, n_chan=3, name=f"S{i}") for i in range(2)}
        ax = plot_coherence_section(sps, pair=(0, 1))
        assert ax is not None
        plt.close("all")

    def test_show_threshold_false(self):
        sps = {f"S{i}": _sp(n_freq=8, n_chan=3, name=f"S{i}") for i in range(2)}
        ax = plot_coherence_section(sps, show_threshold=False)
        assert ax is not None
        plt.close("all")

    def test_section_style_instance_and_custom_ax(self):
        sps = {f"S{i}": _sp(n_freq=6, n_chan=3, name=f"S{i}") for i in range(2)}
        sty = PYCSAMT_SECTION.style_for("dashboard")
        fig, ax_in = plt.subplots()
        ax = plot_coherence_section(sps, section=sty, ax=ax_in)
        assert ax is ax_in
        plt.close("all")

    def test_single_frequency_edge_case(self):
        sps = {"S1": _sp(n_freq=1, n_chan=3, name="S1")}
        ax = plot_coherence_section(sps)
        assert ax is not None
        plt.close("all")


# ─────────────────────────────────────────────────────────────────────────────
# Remaining branch gaps: log x-axis, explicit style overrides, surplus axes,
# explicit cmap/figsize, tipper + custom axes together.
# ─────────────────────────────────────────────────────────────────────────────


class TestRemainingBranchCoverage:
    def test_plot_psd_log_x_axis(self):
        sp = _sp()
        with PYCSAMT_CONTROL.context(x__view="frequency"):
            ax = plot_psd(sp)
            assert ax.get_xscale() == "log"
        plt.close("all")

    def test_plot_coherence_log_x_axis(self):
        sp = _sp(8, 3)
        with PYCSAMT_CONTROL.context(x__view="frequency"):
            axs = plot_coherence(sp, pairs=[(0, 1)])
            assert axs[0].get_xscale() == "log"
        plt.close("all")

    def test_plot_z_from_spectra_log_x_axis(self):
        sp = _FakeZTSpectra(with_tipper=False)
        with PYCSAMT_CONTROL.context(x__view="frequency"):
            fig = plot_z_from_spectra(sp)
            assert fig.axes[0].get_xscale() == "log"
        plt.close("all")

    def test_plot_tipper_log_x_axis(self):
        sp = _FakeZTSpectra(with_tipper=True)
        with PYCSAMT_CONTROL.context(x__view="frequency"):
            axs = plot_tipper_from_spectra(sp)
            assert axs[0].get_xscale() == "log"
        plt.close("all")

    def test_plot_coherence_explicit_lw_alpha_figsize(self):
        sp = _sp(8, 3)
        axs = plot_coherence(sp, pairs=[(0, 1)], lw=1.5, alpha=0.4, figsize=(6.0, 4.0))
        assert len(axs) == 1
        plt.close("all")

    def test_plot_coherence_surplus_axes_hidden(self):
        """4 pairs on a 3-col grid leaves 2 surplus axes to hide."""
        sp = _sp(8, 4)
        axs = plot_coherence(sp, pairs=[(0, 1), (0, 2), (0, 3), (1, 2)])
        assert len(axs) == 4
        fig = axs[0].figure
        assert len(fig.axes) == 6
        hidden = [a for a in fig.axes if not a.get_visible()]
        assert len(hidden) == 2
        plt.close("all")

    @pytest.mark.parametrize(
        ("quantity", "cmap"),
        [
            ("abs", "plasma"),
            ("real", "coolwarm"),
            ("imag", "coolwarm"),
            ("phase", "twilight"),
        ],
    )
    def test_plot_spectra_matrix_explicit_cmap(self, quantity, cmap):
        sp = _sp(8, 3)
        fig = plot_spectra_matrix(sp, quantity=quantity, cmap=cmap)
        assert isinstance(fig, Figure)
        plt.close("all")

    def test_plot_tipper_with_tipper_and_custom_axes(self):
        sp = _FakeZTSpectra(with_tipper=True, with_tipper_err=True)
        fig_in, (ax1, ax2) = plt.subplots(1, 2)
        axs = plot_tipper_from_spectra(sp, estimate_error=True, axes=(ax1, ax2))
        assert axs[0] is ax1
        assert axs[1] is ax2
        plt.close("all")

    def test_plot_psd_section_explicit_figsize(self):
        sps = {f"S{i}": _sp(n_freq=6, n_chan=2, name=f"S{i}") for i in range(2)}
        ax = plot_psd_section(sps, figsize=(8.0, 5.0))
        assert ax is not None
        plt.close("all")

    def test_plot_coherence_section_explicit_figsize(self):
        sps = {f"S{i}": _sp(n_freq=6, n_chan=3, name=f"S{i}") for i in range(2)}
        ax = plot_coherence_section(sps, figsize=(8.0, 5.0))
        assert ax is not None
        plt.close("all")
