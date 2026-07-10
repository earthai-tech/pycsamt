from __future__ import annotations

import numpy as np

from pycsamt.emtools._core import _get_z_block
from pycsamt.emtools.gb import (
    GroomBaileyResult,
    apply_groom_bailey,
    groom_bailey_decomposition,
    groom_bailey_table,
)


class _FakeZ:
    def __init__(self, z, freq):
        self.z = np.asarray(z, dtype=np.complex128)
        self.freq = np.asarray(freq, dtype=float)

    def compute_resistivity_phase(self):
        return None


class _FakeSite:
    def __init__(self, station, z, freq):
        self.station = station
        self.Z = _FakeZ(z, freq)
        self.freq = np.asarray(freq, dtype=float)

    def get_section(self, *_, **__):
        return None


def _freqs(n: int = 12) -> np.ndarray:
    return np.logspace(0.0, 3.0, n)


def _regional_2d(freq: np.ndarray) -> np.ndarray:
    z = np.zeros((freq.size, 2, 2), dtype=np.complex128)
    amp = np.sqrt(5.0 * freq * 100.0)
    z[:, 0, 1] = amp * (1.0 + 1.0j) / np.sqrt(2.0)
    z[:, 1, 0] = -1.4 * amp * (0.8 + 1.2j) / np.sqrt(2.08)
    return z


def _distorted_site(D: np.ndarray) -> _FakeSite:
    fr = _freqs()
    z = D[None, :, :] @ _regional_2d(fr)
    return _FakeSite("S0", z, fr)


def _diag_ratio(z: np.ndarray) -> float:
    diag = np.sqrt(np.abs(z[:, 0, 0]) ** 2 + np.abs(z[:, 1, 1]) ** 2)
    off = np.sqrt(np.abs(z[:, 0, 1]) ** 2 + np.abs(z[:, 1, 0]) ** 2)
    return float(np.nanmedian(diag / np.maximum(off, 1e-24)))


def test_groom_bailey_table_fits_distorted_2d_tensor():
    D = np.array([[1.0, 0.25], [-0.15, 1.1]], dtype=float)
    site = _distorted_site(D)

    table = groom_bailey_table([site], min_freq=4, robust=False)

    assert list(table["station"]) == ["S0"]
    assert table.loc[0, "status"] == "ok"
    assert table.loc[0, "n_freq"] == 12
    assert table.loc[0, "rms_fit"] < 1e-8
    assert table.loc[0, "diagonal_ratio_after"] < 1e-8
    assert table.loc[0, "diagonal_ratio_before"] > 0.05
    assert np.isfinite(table.loc[0, "twist_deg"])
    assert np.isfinite(table.loc[0, "shear"])
    assert np.isfinite(table.loc[0, "anisotropy"])


def test_apply_groom_bailey_removes_diagonal_leakage():
    D = np.array([[1.0, 0.25], [-0.15, 1.1]], dtype=float)
    site = _distorted_site(D)
    before = _diag_ratio(site.Z.z)

    corrected = apply_groom_bailey([site], robust=False)
    item = next(iter(corrected))
    _, z_after, _ = _get_z_block(item)
    after = _diag_ratio(z_after)

    assert before > 0.05
    assert after < 1e-8


def test_groom_bailey_decomposition_result_container():
    D = np.array([[1.0, 0.2], [-0.1, 1.0]], dtype=float)
    result = groom_bailey_decomposition([_distorted_site(D)], apply=True)

    assert isinstance(result, GroomBaileyResult)
    assert result.applied is True
    assert result.n_station == 1
    assert "GroomBaileyResult" in result.summary()


def test_groom_bailey_table_reports_insufficient_frequencies():
    D = np.eye(2)
    site = _distorted_site(D)
    site.Z.z = site.Z.z[:2]
    site.Z.freq = site.Z.freq[:2]

    table = groom_bailey_table([site], min_freq=4)

    assert table.loc[0, "status"] == "insufficient_frequencies"
    assert table.loc[0, "n_freq"] == 2
