# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

from __future__ import annotations

import textwrap
from pathlib import Path

import pytest

from pycsamt.exceptions import EdIDataError
from pycsamt.seg.spectra import (
    SpectraIO,
    SpectraMixin,
    SpectraSECT,
)


def _write(tmp: Path, name: str, text: str) -> Path:
    p = tmp / name
    p.write_text(text, encoding="utf-8")
    return p


def test_spectrasect_from_file_parses_header_and_ids(
    tmp_path: Path,
):
    edi = textwrap.dedent(
        """
        >HEAD
          DATAID="A1"

        >=SPECTRASECT
          SECTID="SP"
          NCHAN=3
          NFREQ=2
          MAXBLKS=99
            // 3
            HX
            HY
            EX

        >SPECTRA FREQ=1.0 ROTSPEC=1 BW=0.25 AVGT=2.0 // 6
          1 2 3
          4 5 6

        >SPECTRA FREQ=2.0 // 4
          10 20
          30 40

        >END
        """
    ).strip("\n")
    p = _write(tmp_path, "sp_hdr.edi", edi)

    sect = SpectraSECT.from_file(str(p))
    assert sect.sectid == "SP"
    assert sect.nchan == 3
    assert sect.nfreq == 2
    assert sect.maxblks == 99
    assert sect.meas_ids == ["HX", "HY", "EX"]
    # start index should point at first >SPECTRA line
    assert isinstance(sect.start_data_lines_num, int)
    assert sect.start_data_lines_num > 0


def test_spectraio_from_file_parses_blocks(tmp_path: Path):
    edi = textwrap.dedent(
        """
        >HEAD
          DATAID="A2"

        >=SPECTRASECT
          SECTID="SP2"
          NCHAN=2
          NFREQ=2

        >SPECTRA FREQ=1.0 ROTSPEC=1 BW=0.25 AVGT=2.0 // 6
          1 2 3
          4 5 6

        >SPECTRA FREQ=2.0 EXTRA=YES // 4
          10 20
          30 40

        >END
        """
    ).strip("\n")
    p = _write(tmp_path, "sp_data.edi", edi)

    sect = SpectraSECT.from_file(str(p))
    io = SpectraIO.from_file(
        str(p), start_line=sect.start_data_lines_num
    )

    assert len(io.blocks) == 2

    b0, b1 = io.blocks
    assert b0.freq == 1.0
    assert b0.rotspec == 1
    assert b0.bw == 0.25
    assert b0.avgt == 2.0
    assert b0.nvals_hint == 6
    assert b0.values == [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]

    assert b1.freq == 2.0
    assert b1.rotspec is None
    assert b1.options.get("extra") == "YES"
    assert b1.nvals_hint == 4
    assert b1.values == [10.0, 20.0, 30.0, 40.0]


def test_spectraio_from_file_auto_finds_blocks(tmp_path: Path):
    edi = textwrap.dedent(
        """
        >HEAD
          DATAID="A2"

        >SPECTRA FREQ=3.0 // 3
          0.1 0.2 0.3

        >END
        """
    ).strip("\n")
    p = _write(tmp_path, "sp_auto.edi", edi)

    io = SpectraIO.from_file(str(p))
    assert len(io.blocks) == 1
    assert io.blocks[0].freq == 3.0
    assert io.blocks[0].values == [0.1, 0.2, 0.3]


def test_spectraio_write_roundtrip_formatting(tmp_path: Path):
    # Build blocks programmatically and serialize
    io = SpectraIO()
    b = io.blocks.append  # short alias

    # first block with known options and extras
    from pycsamt.seg.spectra import _SpectraBlock

    blk1 = _SpectraBlock()
    blk1.freq = 10.0
    blk1.rotspec = 1
    blk1.bw = 0.5
    blk1.avgt = 2.5
    blk1.options["foo"] = "bar"
    blk1.options["a"] = "1"
    blk1.values = [0.1, 0.2, 0.3, 0.4]
    b(blk1)

    # second block, minimal
    blk2 = _SpectraBlock()
    blk2.freq = 5.0
    blk2.values = [1.0]
    b(blk2)

    lines = io.write(per_line=3, float_fmt="{: .3E}")
    text = "".join(lines)

    # header line contains ordered known options first
    assert ">SPECTRA FREQ=10.0 ROTSPEC=1 BW=0.5 AVGT=2.5" in text
    # extras are sorted by key name and upper-cased
    assert " A=1 FOO=BAR" in text
    # values are formatted with the requested fmt and per_line
    # first block -> 0.1 0.2 0.3 on one line, 0.4 on next
    assert "  1.000E-01  2.000E-01  3.000E-01" in text
    assert "  4.000E-01" in text
    # second block single value
    assert ">SPECTRA FREQ=5.0" in text
    assert "  1.000E+00" in text


def test_spectraio_missing_blocks_raises(tmp_path: Path):
    edi = textwrap.dedent(
        """
        >HEAD
          DATAID="X"

        >END
        """
    ).strip("\n")
    p = _write(tmp_path, "sp_missing.edi", edi)

    with pytest.raises(EdIDataError):
        SpectraIO.from_file(str(p))


def test_spectramixin_helpers(tmp_path: Path):
    edi = textwrap.dedent(
        """

        >HEAD
          DATAID="A5"

        >=SPECTRASECT
          SECTID="MX"
          NCHAN=2
          NFREQ=1
            // 2
            HX
            EX

        >SPECTRA FREQ=1.0 // 2
          11 22

        >END
        """
    ).strip("\n")
    p = _write(tmp_path, "sp_mixin.edi", edi)

    class Host(SpectraMixin):
        pass

    sect = Host.from_file(str(p))
    assert isinstance(sect, SpectraSECT)

    io = Host.read_blocks(str(p))
    assert isinstance(io, SpectraIO)
    assert len(io.blocks) == 1
    assert io.blocks[0].values == [11.0, 22.0]

