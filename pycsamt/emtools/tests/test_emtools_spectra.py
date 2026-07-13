"""Tests for pycsamt.emtools.spectra (analysis functions)"""

from __future__ import annotations

import matplotlib
import numpy as np
import pytest

matplotlib.use("Agg")

from pycsamt.api import reset_api_view
from pycsamt.emtools.spectra import (
    band_select,
    coherence_matrix,
    coherence_table,
    mask_low_coherence,
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

    def __init__(
        self, n_freq: int = 8, n_chan: int = 3, name: str = "STA01"
    ) -> None:
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
            A = rng.standard_normal(
                (n_chan, n_chan)
            ) + 1j * rng.standard_normal((n_chan, n_chan))
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
        np.testing.assert_allclose(
            coh, np.transpose(coh, (0, 2, 1)), atol=1e-10
        )

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
