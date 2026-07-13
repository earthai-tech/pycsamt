from __future__ import annotations

import numpy as np
import pandas as pd

from pycsamt.api import APIFrame, reset_api_view
from pycsamt.emtools import dimensionality as dim
from pycsamt.emtools import spectra


class _FakeSpectra:
    def __init__(self) -> None:
        self.name = "STA01"
        self.freq = np.asarray([10.0, 1.0], dtype=float)
        self.bw = np.asarray([0.1, 0.2], dtype=float)
        self.avgt = np.asarray([8.0, 9.0], dtype=float)
        self.avgf = np.asarray([10.0, 1.0], dtype=float)
        self.rotspec = np.asarray([0.0, 0.0], dtype=float)
        self.segnum = np.asarray([4, 5], dtype=int)
        self.band = ["A", "B"]
        self.chan_ids = ["1", "2", "3"]
        self.id_to_chtype = {"1": "Ex", "2": "Ey", "3": "Hx"}
        self._S = np.asarray(
            [
                [
                    [4.0, 1.0 + 0.0j, 0.5 + 0.0j],
                    [1.0 + 0.0j, 9.0, 1.5 + 0.0j],
                    [0.5 + 0.0j, 1.5 + 0.0j, 16.0],
                ],
                [
                    [5.0, 1.2 + 0.0j, 0.4 + 0.0j],
                    [1.2 + 0.0j, 8.0, 1.0 + 0.0j],
                    [0.4 + 0.0j, 1.0 + 0.0j, 14.0],
                ],
            ],
            dtype=complex,
        )

    @property
    def S(self) -> np.ndarray:
        return self._S

    @property
    def n_freq(self) -> int:
        return int(self.freq.size)

    @property
    def n_chan(self) -> int:
        return int(self.S.shape[1])


def _feature_table() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "station": ["S1", "S1"],
            "freq": [10.0, 1.0],
            "period": [0.1, 1.0],
            "beta_abs": [1.0, 5.0],
            "ellipt_abs": [0.1, 0.3],
            "logrho_det": [2.0, 2.2],
            "phi_det": [45.0, 43.0],
            "tip_amp": [0.0, 0.1],
        }
    )


def test_spectra_tables_api_flag_overrides_global():
    reset_api_view()
    sp = _FakeSpectra()

    plain = spectra.psd_table({"STA01": sp}, api=False)
    viewed = spectra.psd_table({"STA01": sp}, api=True)

    assert isinstance(plain, pd.DataFrame)
    assert isinstance(viewed, APIFrame)
    assert viewed.kind == "emtools.spectra.psd"
    assert viewed.df.equals(plain)


def test_spectra_related_tables_return_api_frame():
    reset_api_view()
    sp = _FakeSpectra()

    coh = spectra.coherence_table({"STA01": sp}, api=True)
    snr = spectra.snr_table({"STA01": sp}, api=True)
    summary = spectra.spectra_summary(sp, api=True)

    assert isinstance(coh, APIFrame)
    assert coh.kind == "emtools.spectra.coherence"
    assert isinstance(snr, APIFrame)
    assert snr.kind == "emtools.spectra.snr"
    assert "snr_db" in snr.columns
    assert isinstance(summary, APIFrame)
    assert summary.kind == "emtools.spectra.summary"


def test_classify_dimensionality_api_flag_overrides_global(monkeypatch):
    reset_api_view()
    monkeypatch.setattr(
        dim, "phase_features_table", lambda *a, **k: _feature_table()
    )

    plain = dim.classify_dimensionality(["S1"], api=False)
    viewed = dim.classify_dimensionality(["S1"], api=True)

    assert isinstance(plain, pd.DataFrame)
    assert isinstance(viewed, APIFrame)
    assert viewed.kind == "emtools.dimensionality.classification"
    assert viewed.df.equals(plain)


def test_encode_dimensionality_missing_model_can_return_api_frame():
    reset_api_view()

    viewed = dim.encode_dimensionality(["S1"], {}, api=True)

    assert isinstance(viewed, APIFrame)
    assert viewed.kind == "emtools.dimensionality.encoding"
    assert viewed.df.empty
