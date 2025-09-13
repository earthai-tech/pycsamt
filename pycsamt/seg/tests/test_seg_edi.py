# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

from __future__ import annotations

from pathlib import Path
import numpy as np

from pycsamt.seg.edi import EDIFile
from pycsamt.z.z import Z

# ---------------- helpers to craft minimal EDI ----------------

def _mk_head(station: str = "S1") -> str:
    return (
        ">HEAD\n"
        f"  DATAID={station}\n"
        "  EMPTY=1.0E+32\n"
    )


def _mk_info() -> str:
    return ">INFO\n  // minimal\n"


def _mk_dmeas() -> str:
    return ">DEFINEMEAS\n  // minimal\n"


def _mk_mtsect(sectid: str = "S1") -> str:
    return ">=MTSECT\n  SECTID={}\n".format(sectid)


def _mk_emapsect(sectid: str = "S1") -> str:
    return ">=EMAPSECT\n  SECTID={}\n".format(sectid)


def _mk_freq_block(freq: list[float]) -> str:
    body = "  " + "  ".join(f"{f: .6E}".upper() for f in freq)
    return (
        ">!****FREQUENCIES****!\n"
        f">FREQ  //{len(freq)}\n"
        f"{body}\n"
    )


def _mk_z_blocks(n: int, rot: str = "ROT=ZROT") -> str:
    # only ZXX* to keep tiny; parser fills zeros for rest
    zr = "  " + "  ".join(
        f"{v: .6E}".upper() for v in np.linspace(1e-3, 2e-3, n)
    )
    zi = "  " + "  ".join(
        f"{0.0: .6E}".upper() for _ in range(n)
    )
    zv = "  " + "  ".join(
        f"{1e-6: .6E}".upper() for _ in range(n)
    )
    return (
        ">!****IMPEDANCES****!\n"
        f">ZXXR {rot}  //{n}\n{zr}\n"
        f">ZXXI {rot}  //{n}\n{zi}\n"
        f">ZXX.VAR {rot}  //{n}\n{zv}\n"
    )


def _mk_tipper_blocks(n: int) -> str:
    txr = "  " + "  ".join(
        f"{v: .6E}".upper() for v in np.linspace(0.1, 0.2, n)
    )
    txi = "  " + "  ".join(
        f"{0.0: .6E}".upper() for _ in range(n)
    )
    tvr = "  " + "  ".join(
        f"{1e-4: .6E}".upper() for _ in range(n)
    )
    tyr = txr
    tyi = txi
    tyv = tvr
    trot = "  " + "  ".join(
        f"{0.0: .6E}".upper() for _ in range(n)
    )
    return (
        ">!****TIPPER ROTATION ANGLES****!\n"
        f">TROT  //{n}\n{trot}\n"
        ">!****TIPPER PARAMETERS****!\n"
        f">TXR.EXP ROT=TROT  //{n}\n{txr}\n"
        f">TXI.EXP ROT=TROT  //{n}\n{txi}\n"
        f">TXVAR.EXP ROT=TROT  //{n}\n{tvr}\n"
        f">TYR.EXP ROT=TROT  //{n}\n{tyr}\n"
        f">TYI.EXP ROT=TROT  //{n}\n{tyi}\n"
        f">TYVAR.EXP ROT=TROT  //{n}\n{tyv}\n"
    )


def _mk_spectra_section() -> str:
    # one spectra block with 3 values
    head = ">=SPECTRASECT\n  SECTID=SP1\n"
    blk = (
        ">SPECTRA FREQ=10 BW=1 ROTSPEC=0 // 3\n"
        "  1.0E+00  2.0E+00\n"
        "  3.0E+00\n"
    )
    return head + blk


def _mk_tseries_section() -> str:
    head = (
        ">=TSERIESSECT\n"
        "  SECTID=TSA\n"
        "  NCHAN=2\n"
        "  DT=0.5\n"
    )
    blk1 = (
        ">TSERIES ID=HX NPTS=3 DT=0.1 // 3\n"
        "  1.0E+00  2.0E+00\n"
        "  3.0E+00\n"
    )
    blk2 = (
        ">TSERIES ID=HY NPTS=2 DT=0.2 // 2\n"
        "  1.0E+00  0.0E+00\n"
    )
    return head + blk1 + blk2


def _mk_end() -> str:
    return ">END\n"


# ---------------------------- tests ---------------------------------


def test_read_mt_minimal(tmp_path: Path) -> None:
    # minimal MT with freq, ZXX*, end marker
    edi = (
        _mk_head("MT1")
        + _mk_info()
        + _mk_dmeas()
        + _mk_mtsect("MT1")
        + _mk_freq_block([1.0, 10.0])
        + _mk_z_blocks(2, "ROT=ZROT")
        + _mk_end()
    )
    fn = tmp_path / "mt.edi"
    fn.write_text(edi, encoding="utf-8")

    ed = EDIFile(fn)
    assert ed.Z.n_freq == 2
    assert np.allclose(ed.Z.freq, [10.0 , 1.0]) # reverse ascending 
    # ZXX real non-zero, imag zero as crafted
    assert np.allclose(ed.Z.z[:, 0, 0].imag, 0.0)
    assert np.all(ed.Z.z[:, 0, 0].real > 0.0)


def test_read_tipper_blocks(tmp_path: Path) -> None:
    edi = (
        _mk_head("MT2")
        + _mk_info()
        + _mk_dmeas()
        + _mk_mtsect("MT2")
        + _mk_freq_block([1.0, 2.0, 3.0])
        + _mk_z_blocks(3, "ROT=ZROT")
        + _mk_tipper_blocks(3)
        + _mk_end()
    )
    fn = tmp_path / "mt_tip.edi"
    fn.write_text(edi, encoding="utf-8")

    ed = EDIFile(fn)
    tip = ed.Tip.tipper
    assert tip is not None and tip.shape == (3, 1, 2)
    # magnitude should be non-zero for TX/|TY|
    assert np.all(np.abs(tip[:, 0, 0]) > 0.0)


def test_read_spectra_and_ts(tmp_path: Path) -> None:
    edi = (
        _mk_head("ALL1")
        + _mk_info()
        + _mk_dmeas()
        + _mk_mtsect("ALL1")
        + _mk_freq_block([5.0, 7.0])
        + _mk_z_blocks(2, "ROT=ZROT")
        + _mk_spectra_section()
        + _mk_tseries_section()
        + _mk_end()
    )
    fn = tmp_path / "all.edi"
    fn.write_text(edi, encoding="utf-8")

    ed = EDIFile(fn)

    # spectra was parsed to high-level container
    spec = ed.get_section("spectra")
    assert spec is not None
    assert hasattr(spec, "freq") or hasattr(spec, "to_io")

    # time series present; check channels
    ts = ed.get_section("timeseries")
    assert ts is not None
    ch = ts.channels()
    assert set(ch) == {"HX", "HY"}


def test_write_roundtrip_contains_blocks(tmp_path: Path) -> None:
    edi = (
        _mk_head("W1")
        + _mk_info()
        + _mk_dmeas()
        + _mk_mtsect("W1")
        + _mk_freq_block([1.0, 3.0, 9.0])
        + _mk_z_blocks(3, "ROT=ZROT")
        + _mk_spectra_section()
        + _mk_tseries_section()
        + _mk_end()
    )
    src = tmp_path / "src.edi"
    src.write_text(edi, encoding="utf-8")

    ed = EDIFile(src)
    out = ed.write(savepath=tmp_path)
    txt = Path(out).read_text(encoding="utf-8")

    assert ">FREQ" in txt
    assert ">ZXXR" in txt and ">ZXXI" in txt
    # not asserting headers for spectra/ts; blocks suffice
    assert ">SPECTRA" in txt
    assert ">TSERIES" in txt
    assert txt.strip().endswith(">END")


def test_interpolate_creates_new_grid() -> None:
    ed = EDIFile(path=None)
    f = np.array([1.0, 10.0, 100.0])
    z = np.zeros((3, 2, 2), complex)
    z[:, 0, 0] = (1e-3 + 0j)
    ed.Z._freq = f
    ed.Z._z = z
    ed.Z._z_err = np.zeros_like(z.real)

    new = ed.interpolate([2.0, 5.0, 20.0], kind="linear")
    assert new.n_freq == 3
    assert np.allclose(new.freq, [2.0, 5.0, 20.0])
    # component carried by interp
    assert np.allclose(new.z[:, 0, 0].real, 1e-3)


def test_write_new_edi_swaps(tmp_path: Path) -> None:
    # base file (2 freqs)
    edi = (
        _mk_head("S3")
        + _mk_info()
        + _mk_dmeas()
        + _mk_emapsect("S3")
        + _mk_freq_block([1.0, 10.0])
        + _mk_z_blocks(2, "ROT=NONE")
        + _mk_end()
    )
    src = tmp_path / "base.edi"
    src.write_text(edi, encoding="utf-8")

    ed = EDIFile(src)

    # build a new Z (3 freqs) to swap in
    f2 = np.array([1.0, 3.0, 9.0])
    z2 = np.zeros((3, 2, 2), complex)
    z2[:, 0, 0] = (2e-3 + 0j)
    newZ = Z(z_array=z2, z_err_array=np.zeros_like(z2.real),
             freq=f2)

    out = ed.write_new_edi(
        edi_fn=str(tmp_path / "swapped.edi"),
        Z=newZ,
    )
    txt = Path(out).read_text(encoding="utf-8")
    # new count present
    assert ">FREQ" in txt and "//3" in txt.split(">FREQ")[1]

