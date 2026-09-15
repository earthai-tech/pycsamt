from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from pycsamt.exceptions import EdIDataError
from pycsamt.seg.spectra import Spectra, SpectraIO, SpectraSECT, _SpectraBlock


# ─────────────────────────────────────────────────────────────────────────
# SpectraSECT
# ─────────────────────────────────────────────────────────────────────────


def test_spectrasect_from_file_raises_when_section_missing(tmp_path: Path):
    p = tmp_path / "no_spectrasect.edi"
    p.write_text(
        "\n".join(
            [
                ">HEAD",
                "  DATAID=X",
                "",
                ">=OTHERSECT",
                "  SECTID=X",
                "",
                ">END",
            ]
        ),
        encoding="utf-8",
    )
    with pytest.raises(EdIDataError):
        SpectraSECT.from_file(p)


def test_spectrasect_write_skips_none_values_and_lists_meas_ids():
    sect = SpectraSECT(sectid="S1", nchan=None, nfreq=3, maxblks=None)
    sect.meas_ids = ["HX", "HY"]
    text = "".join(sect.write())
    assert "SECTID=S1" in text
    assert "NCHAN=" not in text
    assert "MAXBLKS=" not in text
    assert "// 2" in text
    assert "HX" in text and "HY" in text


def test_spectrasect_write_no_meas_ids_omits_block():
    sect = SpectraSECT(sectid="S1")
    text = "".join(sect.write())
    assert "//" not in text


# ─────────────────────────────────────────────────────────────────────────
# _SpectraBlock
# ─────────────────────────────────────────────────────────────────────────


def test_spectra_block_kwargs_override():
    blk = _SpectraBlock(freq=10.0, bw=2.0)
    assert blk.freq == 10.0
    assert blk.bw == 2.0


# ─────────────────────────────────────────────────────────────────────────
# SpectraIO
# ─────────────────────────────────────────────────────────────────────────


def test_spectraio_kwargs_override():
    blk = _SpectraBlock(freq=1.0)
    io = SpectraIO(blocks=[blk])
    assert io.blocks == [blk]


def test_spectraio_from_file_raises_when_no_blocks_found(tmp_path: Path):
    p = tmp_path / "no_spectra_blocks.edi"
    p.write_text(
        "\n".join(
            [
                ">HEAD",
                "  DATAID=X",
                "",
                ">=OTHERSECT",
                "  SECTID=X",
                "",
                ">END",
            ]
        ),
        encoding="utf-8",
    )
    with pytest.raises(EdIDataError):
        SpectraIO.from_file(p)


def test_spectraio_parse_block_handles_bad_comment_count_and_options(
    tmp_path: Path,
):
    p = tmp_path / "spectra_edge.edi"
    lines = [
        ">HEAD",
        "  DATAID=X",
        "",
        ">=SPECTRASECT",
        "  SECTID=X",
        "  NCHAN=1",
        "  NFREQ=1",
        "",
        "  HX",
        "",
        ">SPECTRA FREQ=10.0 BADOPT // not-a-number",
        "// a comment line to skip",
        "  1.0 bogus_token 2.0",
        "",
        ">END",
    ]
    p.write_text("\n".join(lines), encoding="utf-8")
    io = SpectraIO.from_file(p)
    assert len(io.blocks) == 1
    blk = io.blocks[0]
    assert blk.nvals_hint is None  # bad comment count -> exception fallback
    assert "badopt" not in blk.options  # no '=' -> skipped, not stored
    assert blk.values == [1.0, 2.0]  # bogus_token tolerated/skipped


def test_spectraio_parse_block_skips_bare_equals_free_token():
    # A token with no '=' inside the option list must be skipped (line 476).
    lines = [
        ">SPECTRA FREQ=10.0 JUSTAWORD // 1",
        "  1.0",
        ">END",
    ]
    blk, next_i = SpectraIO._parse_block(lines, 0)
    assert blk.freq == 10.0
    assert "justaword" not in blk.options


def test_spectraio_dunder_protocol():
    blk1, blk2 = _SpectraBlock(freq=1.0), _SpectraBlock(freq=2.0)
    io = SpectraIO(blocks=[blk1, blk2])
    assert len(io) == 2
    assert list(iter(io)) == [blk1, blk2]
    assert io[0] is blk1
    assert io[1] is blk2


# ─────────────────────────────────────────────────────────────────────────
# Spectra: simple properties
# ─────────────────────────────────────────────────────────────────────────


def test_fcu_cross_spectra_uses_stored_fcu_view_when_present():
    sp = Spectra()
    fake = np.ones((2, 2, 2), complex)
    sp._S_fcu = fake
    assert sp.fcu_cross_spectra is fake


def test_fcu_cross_spectra_falls_back_to_conjugate_of_s():
    sp = Spectra()
    sp._S = np.array([[[1 + 2j]]])
    out = sp.fcu_cross_spectra
    assert np.allclose(out, np.conjugate(sp._S))


def test_missing_mask_none_when_unset():
    sp = Spectra()
    assert sp.missing_mask is None


def test_missing_mask_returns_a_copy():
    sp = Spectra()
    mask = np.array([[True, False]])
    sp._missing_mask = mask
    out = sp.missing_mask
    assert out is not mask
    assert np.array_equal(out, mask)


# ─────────────────────────────────────────────────────────────────────────
# Spectra: _unpack / _unpack_fcu / _pack edge cases
# ─────────────────────────────────────────────────────────────────────────


def test_unpack_raises_when_payload_too_short():
    with pytest.raises(EdIDataError):
        Spectra._unpack(np.array([1.0, 2.0]), 3, empty=1e32)


def test_unpack_fcu_raises_when_payload_too_short():
    with pytest.raises(EdIDataError):
        Spectra._unpack_fcu(np.array([1.0, 2.0]), 3, empty=1e32)


def test_unpack_fcu_marks_missing_for_empty_and_nonfinite_values():
    n = 2
    empty = 1.0e32
    # M = [[empty, 0.5], [empty, 2.0]] (row-major from vals)
    vals = np.array([empty, 0.5, empty, 2.0])
    C, missing = Spectra._unpack_fcu(vals, n, empty=empty)
    assert missing[0, 0]  # diagonal itself is the empty sentinel
    assert missing[0, 1] and missing[1, 0]  # off-diagonal pair, M[1,0]=empty
    assert not missing[1, 1]


def test_pack_raises_for_non_square_matrix():
    with pytest.raises(ValueError):
        Spectra._pack(np.zeros((2, 3), complex))


def test_pack_roundtrips_with_unpack():
    H = np.array([[1.0 + 0j, 2.0 + 3.0j], [2.0 - 3.0j, 4.0 + 0j]])
    packed = Spectra._pack(H)
    rebuilt = Spectra._unpack(packed, 2, empty=1.0e32)
    assert np.allclose(rebuilt, H)


# ─────────────────────────────────────────────────────────────────────────
# Spectra.from_io: nchan-inference fallback branches
# ─────────────────────────────────────────────────────────────────────────


def test_from_io_infers_nchan_from_meas_ids_when_header_nchan_missing():
    sect = SpectraSECT(nchan=None, meas_ids=["HX", "HY"])
    blk = _SpectraBlock(freq=10.0, values=[1.0, 0.0, 0.0, 1.0])
    io = SpectraIO(blocks=[blk])
    sp = Spectra.from_io(sect, io)
    assert sp.n_chan == 2


def test_from_io_infers_nchan_from_block_nvals_hint_perfect_square():
    sect = SpectraSECT(nchan=None, meas_ids=[])
    blk = _SpectraBlock(freq=10.0, values=[1.0, 0.0, 0.0, 1.0])
    blk.nvals_hint = 4
    io = SpectraIO(blocks=[blk])
    sp = Spectra.from_io(sect, io)
    assert sp.n_chan == 2


def test_from_io_infers_nchan_from_values_length_when_no_hint():
    sect = SpectraSECT(nchan=None, meas_ids=[])
    blk = _SpectraBlock(freq=10.0, values=[1.0, 0.0, 0.0, 1.0])
    io = SpectraIO(blocks=[blk])
    sp = Spectra.from_io(sect, io)
    assert sp.n_chan == 2


def test_from_io_raises_when_nchan_cannot_be_inferred():
    sect = SpectraSECT(nchan=None, meas_ids=[])
    io = SpectraIO(blocks=[])
    with pytest.raises(EdIDataError):
        Spectra.from_io(sect, io)


def test_from_io_raises_when_no_blocks():
    sect = SpectraSECT(nchan=2, meas_ids=["HX", "HY"])
    io = SpectraIO(blocks=[])
    with pytest.raises(EdIDataError):
        Spectra.from_io(sect, io)


def test_from_io_bad_header_nchan_falls_back_to_meas_ids():
    sect = SpectraSECT(nchan="not-an-int", meas_ids=["HX", "HY"])
    blk = _SpectraBlock(freq=10.0, values=[1.0, 0.0, 0.0, 1.0])
    io = SpectraIO(blocks=[blk])
    sp = Spectra.from_io(sect, io)
    assert sp.n_chan == 2
