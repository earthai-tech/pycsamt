# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

"""Phase 8 tests for FCU-compatible EDI SPECTRA recovery."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from pycsamt.emtf import (
    EMTF,
    SpectraCovarianceWarning,
    SpectraRecoveryError,
    recover_spectra_transfer_functions,
    resolve_spectra_channels,
)
from pycsamt.exceptions import EdIDataError
from pycsamt.seg import EDIFile
from pycsamt.seg.spectra import (
    Spectra,
    SpectraIO,
    SpectraSECT,
    SpectraValidationWarning,
    _SpectraBlock,
    spectra_from_Z,
)


def _spd_complex(n: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    a = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    return a.conj().T @ a + np.eye(n) * (n + 2.0)


def _make_spectra(
    kinds: tuple[str, ...],
    matrices: list[np.ndarray],
    *,
    avgt: tuple[float, ...] | None = None,
    rotspec: tuple[float, ...] | None = None,
) -> Spectra:
    n_freq = len(matrices)
    sp = Spectra(name="SYN")
    sp._freq = np.asarray([100.0 / (10**i) for i in range(n_freq)])
    sp._S_fcu = np.stack(matrices)
    # Keep the historical pyCSAMT view deliberately conjugated, matching how
    # packed EDI SPECTRA were interpreted before the Phase-8 FCU view existed.
    sp._S = np.conjugate(sp._S_fcu)
    sp._missing_mask = np.zeros(sp._S.shape, dtype=bool)
    sp.chan_ids = [str(i + 1) for i in range(len(kinds))]
    sp.id_to_chtype = {
        str(i + 1): kind for i, kind in enumerate(kinds)
    }
    sp.declared_nfreq = n_freq
    sp.parsed_nfreq = n_freq
    sp._from_edi_spectra = True
    if avgt is None:
        avgt = tuple(8.0 for _ in range(n_freq))
    sp.avgt = np.asarray(avgt, dtype=float)
    sp.avgt_present = np.ones(n_freq, dtype=bool)
    if rotspec is None:
        rotspec = tuple(0.0 for _ in range(n_freq))
    sp.rotspec = np.asarray(rotspec, dtype=float)
    sp.rotspec_present = np.ones(n_freq, dtype=bool)
    sp.bw = np.ones(n_freq, dtype=float)
    sp.avgf = np.full(n_freq, np.nan)
    sp.segnum = np.zeros(n_freq, dtype=int)
    sp.band = [""] * n_freq
    return sp


def _manual_fcu(
    c: np.ndarray,
    *,
    h: tuple[int, int],
    r: tuple[int, int],
    outputs: tuple[int, ...],
    avgt: float,
):
    rh_h = c[np.ix_(r, h)]
    rh_r = c[np.ix_(r, r)]
    hh_h = c[np.ix_(h, h)]
    rh_e = c[np.ix_(r, outputs)]
    hh_e = c[np.ix_(h, outputs)]
    eh_e = c[np.ix_(outputs, outputs)]

    inv_rh_h = np.linalg.inv(rh_h)
    zh = inv_rh_h @ rh_e
    tf = zh.conj().T
    inv_sig = inv_rh_h @ rh_r @ np.linalg.inv(rh_h.conj().T)
    resid = (
        eh_e
        - zh.conj().T @ hh_e
        - hh_e.conj().T @ zh
        + zh.conj().T @ hh_h @ zh
    ) / avgt
    var = np.real(np.diag(resid)[:, None] * np.diag(inv_sig)[None, :])
    return tf, var, inv_sig, resid


def _pack_fcu(c: np.ndarray) -> list[float]:
    n = c.shape[0]
    matrix = np.zeros((n, n), dtype=float)
    for i in range(n):
        matrix[i, i] = float(c[i, i].real)
        for j in range(i + 1, n):
            matrix[j, i] = float(c[i, j].real)
            matrix[i, j] = float(-c[i, j].imag)
    return matrix.ravel().tolist()


def _from_packed(
    c: np.ndarray,
    *,
    declared_nfreq: int = 1,
    avgt: float | None = 8.0,
) -> Spectra:
    kinds = ("HX", "HY", "HZ", "EX", "EY")
    sect = SpectraSECT(
        sectid="PACKED",
        nchan=5,
        nfreq=declared_nfreq,
        meas_ids=["1", "2", "3", "4", "5"],
    )
    sect.id_to_chtype = {
        str(i + 1): kind for i, kind in enumerate(kinds)
    }
    io = SpectraIO()
    block = _SpectraBlock()
    block.freq = 10.0
    block.rotspec = 0
    block.avgt = avgt
    block.values = _pack_fcu(c)
    block.nvals_hint = 25
    io.blocks.append(block)
    return Spectra.from_io(sect, io)


def test_packed_spectra_exposes_exact_fcu_and_legacy_views():
    c = _spd_complex(5, seed=1)
    sp = _from_packed(c)

    np.testing.assert_allclose(sp.fcu_cross_spectra[0], c)
    np.testing.assert_allclose(sp.S[0], np.conjugate(c))


def test_channel_resolution_single_station_and_remote_reference():
    single = _make_spectra(
        ("HX", "HY", "HZ", "EX", "EY"),
        [_spd_complex(5, seed=2)],
    )
    cmap = resolve_spectra_channels(single)
    assert cmap.local_h == (0, 1)
    assert cmap.remote_h == (0, 1)
    assert cmap.outputs == (2, 3, 4)
    assert cmap.reference_type == "single_station"

    remote = _make_spectra(
        ("HX", "HY", "HZ", "EX", "EY", "HX", "HY"),
        [_spd_complex(7, seed=3)],
    )
    cmap = resolve_spectra_channels(remote)
    assert cmap.local_h == (0, 1)
    assert cmap.remote_h == (5, 6)
    assert cmap.outputs == (2, 3, 4)
    assert cmap.reference_type == "remote_reference"


def test_channel_resolution_rejects_incomplete_remote_pair():
    sp = _make_spectra(
        ("HX", "HY", "HZ", "EX", "EY", "HX"),
        [_spd_complex(6, seed=4)],
    )
    with pytest.raises(SpectraRecoveryError, match="incomplete remote"):
        resolve_spectra_channels(sp)


def test_single_station_recovery_matches_fcu_equations():
    matrices = [_spd_complex(5, 10), _spd_complex(5, 11)]
    sp = _make_spectra(
        ("HX", "HY", "HZ", "EX", "EY"),
        matrices,
        avgt=(8.0, 12.0),
    )
    result = recover_spectra_transfer_functions(sp)

    for k, (c, avgt) in enumerate(zip(matrices, (8.0, 12.0))):
        tf, var, inv_sig, resid = _manual_fcu(
            c,
            h=(0, 1),
            r=(0, 1),
            outputs=(2, 3, 4),
            avgt=avgt,
        )
        np.testing.assert_allclose(result.tipper.data[k], tf[0:1])
        np.testing.assert_allclose(result.impedance.data[k], tf[1:3])
        np.testing.assert_allclose(
            result.impedance.get_estimate("VAR").data[k],
            var[1:3],
        )
        np.testing.assert_allclose(
            result.impedance.get_estimate("INVSIGCOV").data[k],
            inv_sig,
        )
        np.testing.assert_allclose(
            result.impedance.get_estimate("RESIDCOV").data[k],
            resid[1:3, 1:3],
        )
        np.testing.assert_allclose(
            result.tipper.get_estimate("RESIDCOV").data[k],
            resid[0:1, 0:1],
        )
        np.testing.assert_allclose(
            result.combined_residual_covariance[k],
            resid,
        )


def test_remote_reference_recovery_matches_fcu_equations():
    matrices = [_spd_complex(7, 20), _spd_complex(7, 21)]
    sp = _make_spectra(
        ("HX", "HY", "HZ", "EX", "EY", "HX", "HY"),
        matrices,
        avgt=(6.0, 9.0),
    )
    result = recover_spectra_transfer_functions(sp)

    for k, (c, avgt) in enumerate(zip(matrices, (6.0, 9.0))):
        tf, var, inv_sig, resid = _manual_fcu(
            c,
            h=(0, 1),
            r=(5, 6),
            outputs=(2, 3, 4),
            avgt=avgt,
        )
        np.testing.assert_allclose(result.tipper.data[k], tf[0:1])
        np.testing.assert_allclose(result.impedance.data[k], tf[1:3])
        np.testing.assert_allclose(
            result.impedance.get_estimate("VAR").data[k],
            var[1:3],
        )
        np.testing.assert_allclose(
            result.impedance.get_estimate("INVSIGCOV").data[k],
            inv_sig,
        )
        np.testing.assert_allclose(
            result.combined_residual_covariance[k],
            resid,
        )


def test_nfreq_mismatch_is_explicit():
    c = _spd_complex(5, 30)
    sp = _from_packed(c, declared_nfreq=2)
    with pytest.raises(EdIDataError, match="NFREQ mismatch"):
        sp.validate_frequency_count()
    with pytest.warns(SpectraValidationWarning, match="NFREQ mismatch"):
        assert not sp.validate_frequency_count(policy="warn")


def test_missing_cross_spectra_raise_or_skip():
    matrices = [_spd_complex(5, 31), _spd_complex(5, 32)]
    sp = _make_spectra(
        ("HX", "HY", "HZ", "EX", "EY"),
        matrices,
    )
    sp._missing_mask[0, 0, 3] = True
    sp._missing_mask[0, 3, 0] = True

    with pytest.raises(SpectraRecoveryError, match="missing required"):
        recover_spectra_transfer_functions(sp)

    with pytest.warns(SpectraCovarianceWarning, match="skipping"):
        result = recover_spectra_transfer_functions(
            sp,
            missing_policy="skip",
        )
    assert result.used_indices.tolist() == [1]
    assert result.skipped_indices.tolist() == [0]


def test_missing_avgt_requires_explicit_policy():
    c = _spd_complex(5, 33)
    sp = _from_packed(c, avgt=None)
    with pytest.raises(SpectraRecoveryError, match="requires positive AVGT"):
        recover_spectra_transfer_functions(sp)

    with pytest.warns(SpectraCovarianceWarning, match="substituting AVGT=1"):
        result = recover_spectra_transfer_functions(
            sp,
            avgt_policy="warn",
        )
    assert result.avgt.tolist() == [1.0]


def _write_synthetic_edi(
    path: Path,
    matrices: list[np.ndarray],
    *,
    remote: bool,
    nf_override: int | None = None,
) -> Path:
    kinds = ["HX", "HY", "HZ", "EX", "EY"]
    if remote:
        kinds.extend(["HX", "HY"])
    nchan = len(kinds)
    nfreq = len(matrices) if nf_override is None else nf_override

    lines = [
        ">HEAD\n",
        ' DATAID="SPECTRA_TEST"\n',
        " LAT=5.0\n",
        " LONG=-3.0\n",
        " ELEV=100\n",
        " UNITS=M\n",
        ' COORDSYS="Geographic North"\n',
        " DATUM=WGS84\n",
        " EMPTY=1.0E+32\n\n",
        ">=DEFINEMEAS\n",
        f" MAXCHAN={nchan}\n",
        f" MAXMEAS={nchan}\n",
        " MAXRUN=1\n",
        " UNITS=M\n",
    ]
    for i, kind in enumerate(kinds, start=1):
        if kind.startswith("E"):
            lines.append(
                f">EMEAS ID={i} CHTYPE={kind} X=0 Y=0 "
                "X2=100 Y2=0\n"
            )
        else:
            lines.append(
                f">HMEAS ID={i} CHTYPE={kind} X=0 Y=0 Z=0 "
                "AZM=0\n"
            )
    lines.extend(
        [
            "\n>=SPECTRASECT\n",
            ' SECTID="SPECTRA_TEST"\n',
            f" NCHAN={nchan}\n",
            f" NFREQ={nfreq}\n",
            f" // {nchan}\n",
        ]
    )
    lines.extend(f" {i}\n" for i in range(1, nchan + 1))
    for k, matrix in enumerate(matrices):
        freq = 100.0 / (10**k)
        lines.append(
            f"\n>SPECTRA FREQ={freq:.12g} ROTSPEC=0 "
            f"BW=1 AVGT={8 + k} // {nchan * nchan}\n"
        )
        vals = _pack_fcu(matrix)
        for start in range(0, len(vals), 6):
            chunk = vals[start : start + 6]
            lines.append(" " + " ".join(f"{x:.17g}" for x in chunk) + "\n")
    lines.append("\n>END\n")
    path.write_text("".join(lines), encoding="utf-8")
    return path


def test_emtf_from_edi_prefers_spectra_and_recovers_full_covariance(
    tmp_path: Path,
):
    matrices = [_spd_complex(5, 41), _spd_complex(5, 42)]
    path = _write_synthetic_edi(
        tmp_path / "single_spectra.edi",
        matrices,
        remote=False,
    )
    edi = EDIFile(path)
    doc = EMTF.from_edi(edi)

    assert doc.attrs["source_format"] == "edi_spectra"
    assert doc.processing.remote_reference.reference_type == "Single Station"
    assert doc.impedance.get_estimate("INVSIGCOV") is not None
    assert doc.impedance.get_estimate("RESIDCOV") is not None
    assert doc.tipper_tf.get_estimate("RESIDCOV") is not None
    assert doc.metadata["edi_spectra"]["full_covariance_recovered"]
    assert "spectra_combined_residual_covariance" in doc.attrs


def test_spectra_conversion_uses_phase7_rotation_engine(tmp_path: Path):
    matrices = [_spd_complex(5, 43), _spd_complex(5, 44)]
    path = _write_synthetic_edi(
        tmp_path / "rot_spectra.edi",
        matrices,
        remote=False,
    )
    doc = EMTF.from_edi_spectra(path)
    expected = doc.rotate(30.0, source_angles={
        "impedance": np.zeros(doc.n_periods),
        "tipper": np.zeros(doc.n_periods),
    })
    actual = EMTF.from_edi_spectra(path, target_angle=30.0)

    np.testing.assert_allclose(actual.z, expected.z)
    np.testing.assert_allclose(actual.tipper, expected.tipper)
    np.testing.assert_allclose(
        actual.impedance.get_estimate("INVSIGCOV").data,
        expected.impedance.get_estimate("INVSIGCOV").data,
    )


def test_malformed_nfreq_is_rejected_from_edi(tmp_path: Path):
    matrices = [_spd_complex(5, 45), _spd_complex(5, 46)]
    path = _write_synthetic_edi(
        tmp_path / "bad_nfreq.edi",
        matrices,
        remote=False,
        nf_override=3,
    )
    with pytest.raises(EdIDataError, match="NFREQ mismatch"):
        EMTF.from_edi_spectra(path)


def test_remote_edi_spectra_resolves_second_h_pair(tmp_path: Path):
    matrices = [_spd_complex(7, 47), _spd_complex(7, 48)]
    path = _write_synthetic_edi(
        tmp_path / "remote_spectra.edi",
        matrices,
        remote=True,
    )
    doc = EMTF.from_edi_spectra(path)

    assert doc.processing.remote_reference.reference_type == "Remote Reference"
    info = doc.metadata["edi_spectra"]
    assert info["reference_type"] == "remote_reference"
    assert info["channel_types"] == [
        "HX",
        "HY",
        "HZ",
        "EX",
        "EY",
        "HX",
        "HY",
    ]


def test_prefer_spectra_false_keeps_legacy_fallback_explicit(tmp_path: Path):
    matrices = [_spd_complex(5, 49), _spd_complex(5, 50)]
    path = _write_synthetic_edi(
        tmp_path / "legacy_fallback.edi",
        matrices,
        remote=False,
    )
    rich = EMTF.from_edi(path)
    legacy = EMTF.from_edi(path, prefer_spectra=False)

    assert rich.attrs["source_format"] == "edi_spectra"
    assert legacy.attrs["source_format"] == "edi"
    assert rich.impedance.get_estimate("INVSIGCOV") is not None
    assert legacy.impedance.get_estimate("INVSIGCOV") is None


def test_from_edi_spectra_accepts_bare_spectra_container(tmp_path: Path):
    """``Spectra`` itself has no ``to_emtf`` convenience — EDI/spectra to
    EMTF conversion is entry-point-side only (mirrors ``EMTF.from_edi``),
    so ``pycsamt.seg`` stays independent of ``pycsamt.emtf``. This checks
    the entry point still accepts a bare ``Spectra`` (no ``EDIFile``) and
    matches the ``EDIFile``-backed conversion for the shared TF content.
    """
    matrices = [_spd_complex(5, 51), _spd_complex(5, 52)]
    path = _write_synthetic_edi(
        tmp_path / "spectra_convenience.edi",
        matrices,
        remote=False,
    )
    edi = EDIFile(path)
    sp = edi.get_section("spectra")

    bare = EMTF.from_edi_spectra(sp)
    expected = EMTF.from_edi_spectra(edi)
    np.testing.assert_allclose(bare.z, expected.z)
    np.testing.assert_allclose(bare.tipper, expected.tipper)


def test_legacy_synthetic_spectra_is_not_silently_reinterpreted():
    rng = np.random.default_rng(99)
    z_data = rng.normal(size=(2, 2, 2)) + 1j * rng.normal(
        size=(2, 2, 2)
    )

    class DummyZ:
        pass

    z_obj = DummyZ()
    z_obj.z = z_data
    z_obj.freq = np.asarray([10.0, 1.0])
    z_obj.name = "synthetic"
    z_obj.verbose = 0

    sp = spectra_from_Z(z_obj)
    assert sp._S_fcu is None
    np.testing.assert_allclose(sp.fcu_cross_spectra, np.conjugate(sp.S))

    with pytest.raises(
        SpectraRecoveryError,
        match="requires an EDI SPECTRA container",
    ):
        recover_spectra_transfer_functions(sp, avgt_policy="unit")

    # The historical synthetic path remains available and recovers the
    # original transfer function in the convention it was designed for.
    recovered, _ = sp.to_Z(estimate_error=False)
    np.testing.assert_allclose(recovered.z, z_data)
