import numpy as np
import pytest

from pycsamt.seg.edi import EDIFile
from pycsamt.seg.ops import synthesize_spectra_from_z
from pycsamt.seg.spectra import Spectra, spectra_from_Z


def _rand_cx(shape, scale=1.0, rng=None):
    rng = np.random.default_rng() if rng is None else rng
    a = rng.standard_normal(shape) + 1j * rng.standard_normal(
        shape
    )
    return a * scale


def _make_spd_s_hh(nf, rng=None):
    rng = np.random.default_rng() if rng is None else rng
    A = _rand_cx((nf, 2, 2), rng=rng)
    # S = A A^H + eps*I, per frequency
    S = np.einsum("fij,fkj->fik", A, np.conjugate(A))
    eps = 1e-6
    S[:, 0, 0] += eps
    S[:, 1, 1] += eps
    return S


def test_synthesize_spectra_from_z_synth_no_hz():
    nf = 8
    rng = np.random.default_rng(123)
    Z = _rand_cx((nf, 2, 2), rng=rng)
    S_HH = _make_spd_s_hh(nf, rng=rng)

    Sfull, order = synthesize_spectra_from_z(
        Z,
        S_HH=S_HH,
        include_hz=False,
        chan_order=("HX", "HY", "EX", "EY"),
    )
    assert Sfull.shape == (nf, 4, 4)
    assert list(order) == ["HX", "HY", "EX", "EY"]

    # Hermitian check
    Sswap = np.swapaxes(np.conjugate(Sfull), -1, -2)
    assert np.allclose(Sfull, Sswap, rtol=1e-10, atol=1e-10)

    idx = {k: i for i, k in enumerate(order)}
    # References
    S_EH_ref = np.einsum("fij,fjk->fik", Z, S_HH)
    S_EE_ref = np.einsum(
        "fij,fjk,flk->fil", Z, S_HH, np.conjugate(Z)
    )

    # E/H blocks
    assert np.allclose(
        Sfull[:, idx["EX"], idx["HX"]], S_EH_ref[:, 0, 0]
    )
    assert np.allclose(
        Sfull[:, idx["EX"], idx["HY"]], S_EH_ref[:, 0, 1]
    )
    assert np.allclose(
        Sfull[:, idx["EY"], idx["HX"]], S_EH_ref[:, 1, 0]
    )
    assert np.allclose(
        Sfull[:, idx["EY"], idx["HY"]], S_EH_ref[:, 1, 1]
    )

    # E/E blocks
    assert np.allclose(
        Sfull[:, idx["EX"], idx["EX"]], S_EE_ref[:, 0, 0]
    )
    assert np.allclose(
        Sfull[:, idx["EX"], idx["EY"]], S_EE_ref[:, 0, 1]
    )
    assert np.allclose(
        Sfull[:, idx["EY"], idx["EY"]], S_EE_ref[:, 1, 1]
    )
    # symmetry for E/E off-diagonal
    assert np.allclose(
        Sfull[:, idx["EY"], idx["EX"]],
        np.conjugate(S_EE_ref[:, 0, 1]),
    )


def test_synthesize_spectra_from_z_synth_with_hz():
    nf = 6
    rng = np.random.default_rng(321)
    Z = _rand_cx((nf, 2, 2), rng=rng)
    S_HH = _make_spd_s_hh(nf, rng=rng)
    T = _rand_cx((nf, 1, 2), rng=rng)

    Sfull, order = synthesize_spectra_from_z(
        Z,
        S_HH=S_HH,
        tipper=T,
        include_hz=True,
        chan_order=("HX", "HY", "HZ", "EX", "EY"),
    )
    assert Sfull.shape == (nf, 5, 5)
    assert list(order) == ["HX", "HY", "HZ", "EX", "EY"]

    idx = {k: i for i, k in enumerate(order)}
    S_ZH_ref = np.einsum("fik,fkj->fij", T, S_HH)
    # HZ/HX
    assert np.allclose(
        Sfull[:, idx["HZ"], idx["HX"]], S_ZH_ref[:, 0, 0]
    )
    # HX/HZ conjugate
    assert np.allclose(
        Sfull[:, idx["HX"], idx["HZ"]],
        np.conjugate(S_ZH_ref[:, 0, 0]),
    )


def test_spectra_from_Z_function_and_class_synthetic():
    nf = 7
    rng = np.random.default_rng(111)
    Z = _rand_cx((nf, 2, 2), rng=rng)
    S_HH = _make_spd_s_hh(nf, rng=rng)
    T = _rand_cx((nf, 1, 2), rng=rng)

    class DummyZ:
        pass

    z_obj = DummyZ()
    z_obj.z = Z
    z_obj.freq = np.geomspace(0.1, 100.0, nf)
    z_obj.name = "siteX"
    z_obj.verbose = 0

    sp_fun = spectra_from_Z(
        z_obj,
        S_HH=S_HH,
        tipper=T,
        include_hz=True,
        chan_order=("HX", "HY", "HZ", "EX", "EY"),
    )
    assert sp_fun.n_freq == nf
    assert sp_fun.S.shape == (nf, 5, 5)
    assert sp_fun.chan_ids[:3] == ["HX", "HY", "HZ"]

    sp_cls = Spectra.from_Z(
        z_obj,
        S_HH=S_HH,
        tipper=T,
        include_hz=True,
        chan_order=("HX", "HY", "HZ", "EX", "EY"),
    )
    assert np.allclose(sp_fun.freq, sp_cls.freq)
    assert np.allclose(sp_fun.S, sp_cls.S)


@pytest.mark.usefixtures("edi_spe_file")
def test_realfile_from_spectra_roundtrip(edi_spe_file):
    ed = EDIFile(edi_spe_file)
    sp = ed.spectra
    assert sp is not None
    assert sp.n_freq >= 1

    # Recover Z from spectra (mapping auto via sect)
    Zhat, That = sp.to_Z()
    assert Zhat.n_freq == sp.n_freq
    # Build spectra back from Z with unit power
    sp2 = Spectra.from_Z(
        Zhat,
        include_hz=That is not None,
        tipper=That if That is not None else None,
        chan_order=("HX", "HY", "EX", "EY"),
    )
    assert sp2.n_freq == sp.n_freq
    assert sp2.S.shape[1] in (4, 5)

    # Hermitian sanity on synthesized spectra
    Sswap = np.swapaxes(np.conjugate(sp2.S), -1, -2)
    assert np.allclose(sp2.S, Sswap, rtol=1e-10, atol=1e-10)


def test_synthesize_spectra_from_z_errors():
    nf = 4
    Z = np.zeros((nf, 2, 2), complex)

    with pytest.raises(ValueError):
        synthesize_spectra_from_z(
            Z,
            include_hz=True,
            chan_order=("HX", "HY", "HZ", "EX", "EY"),
        )

    with pytest.raises(ValueError):
        bad_S_HH = np.zeros((nf, 3, 3), complex)
        synthesize_spectra_from_z(Z, S_HH=bad_S_HH)
