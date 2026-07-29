"""
Tests for pycsamt.emtools.legacy — v2 compatibility.

All functions are imported via the emtools namespace and exercised
with v2 :class:`~pycsamt.z.z.Z` objects (the primary data carrier
after the Edi→Sites refactor).  Plotting tests use the Agg backend.
"""

from __future__ import annotations

import warnings

import numpy as np
import pytest

# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────


def _freqs(n: int = 6, f_lo: float = 1.0, f_hi: float = 1e4) -> np.ndarray:
    return np.logspace(np.log10(f_lo), np.log10(f_hi), n)


def _make_z(freqs: np.ndarray, rho: float = 100.0):
    """Minimal v2 Z object with rho_a=rho, phase=45°."""
    from pycsamt.z.z import Z

    z_abs = np.sqrt(5.0 * freqs * rho)
    za = np.zeros((len(freqs), 2, 2), dtype=complex)
    za[:, 0, 1] = z_abs * (1 + 1j) / np.sqrt(2)
    za[:, 1, 0] = -z_abs * (1 + 1j) / np.sqrt(2)
    return Z(z_array=za, freq=freqs)


def _z_list(n: int = 4, rho: float = 100.0) -> list:
    fr = _freqs()
    return [_make_z(fr, rho * (1 + 0.1 * i)) for i in range(n)]


# ─────────────────────────────────────────────────────────────────────────────
# Import tests — all legacy symbols must be importable
# ─────────────────────────────────────────────────────────────────────────────


def test_import_legacy_module():
    from pycsamt.emtools import legacy  # noqa: F401


def test_all_symbols_importable():
    from pycsamt.emtools.legacy import (
        align_tensor,
        check_em_kind,
        compute_qc,
        export_edis,
        extract_z_list,
        full_freq,
        parse_tensor,
        plot_confidence,
        plot_lcurve,
        plot_station_tensors,
        plot_strike,
        plot_tensors,
        tensor2d,
        wrap_phase,
    )

    for fn in (
        check_em_kind,
        extract_z_list,
        parse_tensor,
        compute_qc,
        full_freq,
        tensor2d,
        align_tensor,
        export_edis,
        plot_confidence,
        plot_strike,
        plot_tensors,
        plot_station_tensors,
        wrap_phase,
        plot_lcurve,
    ):
        assert callable(fn)


def test_deprecated_aliases_importable():
    """Aliases installed by install_compat_aliases must exist."""
    import pycsamt.emtools.legacy as leg

    for name in (
        "get_full_frequency",
        "qc",
        "plot_tensors2",
        "plot_l_curve",
        "get2dtensor",
    ):
        assert hasattr(leg, name), f"missing alias: {name!r}"


def test_deprecated_aliases_emit_future_warning():
    import pycsamt.emtools.legacy as leg

    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        # calling a deprecated alias (get_full_frequency = full_freq)
        zs = _z_list(2)
        leg.get_full_frequency(zs)
    categories = [str(x.category) for x in w]
    assert any(
        "FutureWarning" in c for c in categories
    ), "No FutureWarning raised for deprecated alias"


# ─────────────────────────────────────────────────────────────────────────────
# check_em_kind
# ─────────────────────────────────────────────────────────────────────────────


def test_check_em_kind_z_list():
    from pycsamt.emtools.legacy import check_em_kind

    assert check_em_kind(_z_list()) == "Z"


def test_check_em_kind_single_z():
    from pycsamt.emtools.legacy import check_em_kind

    z = _make_z(_freqs())
    assert check_em_kind([z]) == "Z"


def test_check_em_kind_empty_raises():
    from pycsamt.emtools.legacy import check_em_kind

    with pytest.raises(TypeError):
        check_em_kind([])


def test_check_em_kind_string_raises():
    from pycsamt.emtools.legacy import check_em_kind

    with pytest.raises(TypeError):
        check_em_kind("notalist")


def test_check_em_kind_mixed_raises():
    """A mix of Z and a plain object should raise."""
    from pycsamt.emtools.legacy import check_em_kind
    from pycsamt.exceptions import EMError

    z = _make_z(_freqs())
    with pytest.raises((EMError, Exception)):
        check_em_kind([z, object()])


# ─────────────────────────────────────────────────────────────────────────────
# extract_z_list
# ─────────────────────────────────────────────────────────────────────────────


def test_extract_z_list_returns_list():
    from pycsamt.emtools.legacy import extract_z_list

    zs = _z_list(3)
    result = extract_z_list(zs)
    assert isinstance(result, list)
    assert len(result) == 3


def test_extract_z_list_identity_for_z():
    from pycsamt.emtools.legacy import extract_z_list
    from pycsamt.z.z import Z

    zs = _z_list(2)
    result = extract_z_list(zs)
    assert all(isinstance(r, Z) for r in result)


# ─────────────────────────────────────────────────────────────────────────────
# parse_tensor
# ─────────────────────────────────────────────────────────────────────────────


@pytest.mark.parametrize(
    "token,expected_name,expected_comp",
    [
        ("zxy", "z", "xy"),
        ("zyx", "z", "yx"),
        ("resxy", "resistivity", "xy"),
        ("phaseyx", "phase", "yx"),
        ("freq", "_freq", None),
        ("rhoaxy", "resistivity", "xy"),
    ],
)
def test_parse_tensor_valid_tokens(token, expected_name, expected_comp):
    from pycsamt.emtools.legacy import parse_tensor

    name, comp = parse_tensor(token)
    assert name == expected_name
    assert comp == expected_comp


def test_parse_tensor_explicit_tensor_comp():
    from pycsamt.emtools.legacy import parse_tensor

    name, comp = parse_tensor(tensor="res", component="yx")
    assert name == "resistivity"
    assert comp == "yx"


def test_parse_tensor_invalid_raises():
    from pycsamt.emtools.legacy import parse_tensor

    with pytest.raises(ValueError):
        parse_tensor("unknownfoo")


def test_parse_tensor_one_of_pair_raises():
    from pycsamt.emtools.legacy import parse_tensor
    from pycsamt.exceptions import EMError

    with pytest.raises((EMError, ValueError)):
        parse_tensor(tensor="res")


# ─────────────────────────────────────────────────────────────────────────────
# full_freq
# ─────────────────────────────────────────────────────────────────────────────


def test_full_freq_returns_array():
    from pycsamt.emtools.legacy import full_freq

    zs = _z_list(4)
    f = full_freq(zs)
    assert isinstance(f, np.ndarray)
    assert f.ndim == 1
    assert f.size > 0


def test_full_freq_shape_matches_longest_site():
    from pycsamt.emtools.legacy import full_freq

    f_long = _freqs(8)
    f_short = _freqs(4)
    zs = [_make_z(f_long), _make_z(f_short)]
    f = full_freq(zs)
    assert len(f) == 8


def test_full_freq_to_log10():
    from pycsamt.emtools.legacy import full_freq

    zs = _z_list(3)
    flog = full_freq(zs, to_log10=True)
    assert np.all(np.isfinite(flog))


# ─────────────────────────────────────────────────────────────────────────────
# tensor2d
# ─────────────────────────────────────────────────────────────────────────────


def test_tensor2d_shape():
    from pycsamt.emtools.legacy import tensor2d

    n_sta = 5
    n_freq = 6
    fr = _freqs(n_freq)
    zs = [_make_z(fr) for _ in range(n_sta)]
    mat = tensor2d(zs, tensor="res", component="xy")
    assert mat.shape == (n_freq, n_sta)


def test_tensor2d_return_freqs():
    from pycsamt.emtools.legacy import tensor2d

    zs = _z_list(3)
    mat, fr = tensor2d(zs, tensor="res", component="xy", return_freqs=True)
    assert isinstance(fr, np.ndarray)
    assert mat.shape[0] == len(fr)


def test_tensor2d_positive_resistivity():
    from pycsamt.emtools.legacy import tensor2d

    zs = _z_list(4, rho=100.0)
    mat = tensor2d(zs, tensor="res", component="xy")
    finite = mat[np.isfinite(mat)]
    assert (finite > 0).all()


def test_tensor2d_z_component():
    from pycsamt.emtools.legacy import tensor2d

    zs = _z_list(3)
    mat = tensor2d(zs, tensor="z", component="xy", kind="modulus")
    assert mat.shape[1] == 3
    assert (mat[np.isfinite(mat)] >= 0).all()


def test_tensor2d_mixed_freq_fills_nan():
    """Site with a strict subset of frequencies → NaN in the gaps."""
    from pycsamt.emtools.legacy import tensor2d

    # f6 is the reference; f3 takes every other point (a clean subset)
    f6 = _freqs(6)
    f3 = f6[::2]  # indices 0, 2, 4 — guaranteed subset of f6
    zs = [_make_z(f6), _make_z(f3)]
    mat = tensor2d(zs, tensor="res", component="xy")
    assert mat.shape[1] == 2
    assert np.isnan(mat[:, 1]).any()  # gaps at indices 1, 3, 5


# ─────────────────────────────────────────────────────────────────────────────
# align_tensor
# ─────────────────────────────────────────────────────────────────────────────


def test_align_tensor_no_gaps():
    from pycsamt.emtools.legacy import align_tensor

    fr = _freqs(5)
    z = np.ones(5, dtype=complex)
    result = align_tensor(fr, fr, z)
    assert np.allclose(result, z)


def test_align_tensor_with_gaps():
    from pycsamt.emtools.legacy import align_tensor

    fr = _freqs(5)
    site_fr = fr[[0, 2, 4]]  # skip indices 1 and 3
    z = np.array([1 + 0j, 2 + 0j, 3 + 0j])
    result = align_tensor(fr, site_fr, z)
    assert len(result) == 5
    assert np.isnan(result[1].real)
    assert np.isnan(result[3].real)
    assert result[0] == pytest.approx(1.0)
    assert result[2] == pytest.approx(2.0)
    assert result[4] == pytest.approx(3.0)


def test_align_tensor_out_of_grid_raises():
    """Frequencies not in ref_freq violate the contract → EMError."""
    from pycsamt.emtools.legacy import align_tensor
    from pycsamt.exceptions import EMError

    ref_fr = _freqs(4)
    site_fr = np.array([99999.0])  # not in ref
    z = np.array([1.0 + 0j])
    with pytest.raises(EMError):
        align_tensor(ref_fr, site_fr, z)


# ─────────────────────────────────────────────────────────────────────────────
# compute_qc
# ─────────────────────────────────────────────────────────────────────────────


def test_compute_qc_returns_rate():
    from pycsamt.emtools.legacy import compute_qc

    zs = _z_list(4)
    (rate,) = compute_qc(zs)
    assert 0.0 <= rate <= 1.0


def test_compute_qc_full_data_rate_one():
    from pycsamt.emtools.legacy import compute_qc

    zs = _z_list(5, rho=100.0)
    (rate,) = compute_qc(zs)
    assert rate == pytest.approx(1.0)


def test_compute_qc_return_freq():
    from pycsamt.emtools.legacy import compute_qc

    zs = _z_list(4)
    rate, fr = compute_qc(zs, return_freq=True)
    assert isinstance(fr, np.ndarray)
    assert fr.size > 0


def test_compute_qc_return_qco():
    from pycsamt.api.bunch import Bunch
    from pycsamt.emtools.legacy import compute_qc

    zs = _z_list(3)
    rep = compute_qc(zs, return_qco=True)
    assert isinstance(rep, Bunch)
    assert hasattr(rep, "rate_")
    assert hasattr(rep, "freqs_")
    assert hasattr(rep, "data_")


def test_compute_qc_tolerance_filter():
    """Strict tolerance drops more frequencies than loose tolerance."""
    from pycsamt.emtools.legacy import compute_qc

    # Build a 2-site collection where site B has only a subset of freqs
    f_long = _freqs(8)
    f_sub = f_long[::2]  # 4-element strict subset → 4 NaN rows in site B
    zs = [_make_z(f_long), _make_z(f_sub)]
    # tol=0.01 is very strict (drop any row with >1 % missing)
    rate_strict, fr_strict = compute_qc(zs, tol=0.01, return_freq=True)
    # tol=0.99 keeps almost everything
    rate_loose, fr_loose = compute_qc(zs, tol=0.99, return_freq=True)
    assert len(fr_strict) <= len(fr_loose)


# ─────────────────────────────────────────────────────────────────────────────
# wrap_phase  (pure math — no EDI dependency)
# ─────────────────────────────────────────────────────────────────────────────


def test_wrap_phase_default_0_360():
    from pycsamt.emtools.legacy import wrap_phase

    x = np.array([-360.0, 0.0, 180.0, 360.0, 540.0])
    out = wrap_phase(x)
    assert ((out >= 0) & (out < 360)).all()


def test_wrap_phase_symmetric_range():
    from pycsamt.emtools.legacy import wrap_phase

    x = np.array([-270.0, -90.0, 0.0, 90.0, 270.0])
    out = wrap_phase(x, value_range=(-180.0, 180.0))
    assert ((out >= -180) & (out < 180)).all()


def test_wrap_phase_scalar_range():
    from pycsamt.emtools.legacy import wrap_phase

    x = np.array([0.0, 90.0, 180.0, 270.0])
    out = wrap_phase(x, value_range=180.0)
    assert ((out >= 0) & (out < 180)).all()


def test_wrap_phase_invalid_mod_base():
    from pycsamt.emtools.legacy import wrap_phase

    with pytest.raises(ValueError):
        wrap_phase(np.array([1.0]), mod_base=45)


def test_wrap_phase_180_base():
    from pycsamt.emtools.legacy import wrap_phase

    x = np.array([0.0, 90.0, 180.0, 270.0])
    out = wrap_phase(x, mod_base=180)
    assert ((out >= 0) & (out < 180)).all()


# ─────────────────────────────────────────────────────────────────────────────
# plot_lcurve  (pure visualisation, no EDI dependency)
# ─────────────────────────────────────────────────────────────────────────────


@pytest.mark.filterwarnings("ignore")
def test_plot_lcurve_returns_axes():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from pycsamt.emtools.legacy import plot_lcurve

    n = 10
    rms = np.linspace(1.0, 0.1, n)
    rough = np.linspace(0.1, 10.0, n)
    lam = np.logspace(2, -1, n)
    ax = plot_lcurve(rms, rough, lam)
    assert ax is not None
    plt.close("all")


@pytest.mark.filterwarnings("ignore")
def test_plot_lcurve_no_lambda():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from pycsamt.emtools.legacy import plot_lcurve

    rms = np.linspace(1.0, 0.1, 8)
    rough = np.linspace(0.1, 8.0, 8)
    ax = plot_lcurve(rms, rough)
    assert ax is not None
    plt.close("all")


# ─────────────────────────────────────────────────────────────────────────────
# plot_confidence (smoke test)
# ─────────────────────────────────────────────────────────────────────────────


@pytest.mark.filterwarnings("ignore")
def test_plot_confidence_1d_smoke():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from pycsamt.emtools.legacy import plot_confidence
    from pycsamt.z.z import Z

    fr = _freqs(6)
    err = np.zeros((6, 2, 2))
    err[:, 0, 1] = 0.05  # small constant error
    zs = []
    for _ in range(5):
        z_abs = np.sqrt(5.0 * fr * 100.0)
        za = np.zeros((6, 2, 2), dtype=complex)
        za[:, 0, 1] = z_abs * (1 + 1j) / np.sqrt(2)
        za[:, 1, 0] = -z_abs * (1 + 1j) / np.sqrt(2)
        zo = Z(z_array=za, freq=fr, z_err_array=err.copy())
        zs.append(zo)
    ax = plot_confidence(zs, view="1d")
    assert ax is not None
    plt.close("all")


# ─────────────────────────────────────────────────────────────────────────────
# emtools namespace: legacy functions reachable via top-level import
# ─────────────────────────────────────────────────────────────────────────────


def test_legacy_functions_in_emtools_namespace():
    import pycsamt.emtools as em

    for name in (
        "check_em_kind",
        "extract_z_list",
        "parse_tensor",
        "compute_qc",
        "full_freq",
        "tensor2d",
        "align_tensor",
        "export_edis",
        "plot_confidence",
        "plot_strike",
        "plot_tensors",
        "plot_station_tensors",
        "wrap_phase",
        "plot_lcurve",
    ):
        assert hasattr(em, name), f"pycsamt.emtools missing: {name!r}"
        assert callable(getattr(em, name))
