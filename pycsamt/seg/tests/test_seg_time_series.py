# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
from __future__ import annotations

import textwrap
from pathlib import Path

import pytest

from pycsamt.exceptions import EdIDataError
from pycsamt.seg.time_series import (
    TSIO,
    TimeSeriesMixin,
    TSect,
)


def _write(tmp: Path, name: str, txt: str) -> Path:
    p = tmp / name
    p.write_text(txt, encoding="utf-8")
    return p


def test_tsect_from_file_parses_header_and_ids(tmp_path: Path):
    edi = textwrap.dedent(
        """
        >HEAD
          DATAID="TS"
        >=TSERIESSECT
          SECTID=SITE_A
          NCHAN=5
          NMEAS=3
          NPTS=1024
          MAXBLKS=2
          DT=0.125
          WINDOW=HANN

           // 3
            HX
            HY
            EX

        >TSERIES ID=HX NPTS=4 DT=0.125 // 4
          1 2
          3 4
        >END
        """
    ).strip("\n")
    p = _write(tmp_path, "ts_header.edi", edi)

    head = TSect.from_file(str(p))
    assert head.sectid == "SITE_A"
    assert head.nchan == 5
    assert head.nmeas == 3
    assert head.npts == 1024
    assert head.maxblks == 2
    assert head.dt == 0.125
    assert head.extra.get("window") == "HANN"
    assert head.meas_ids == ["HX", "HY", "EX"]
    # points to the first >TSERIES line
    with p.open("r", encoding="utf-8") as f:
        lines = f.readlines()
    assert lines[head.start_data_lines_num].lstrip().startswith(
        ">TSERIES"
    )

    # round-trip write includes ordered keys and ID list
    out = "".join(head.write())
    assert ">=TSERIESSECT" in out
    assert "  SECTID=SITE_A" in out
    assert "  NCHAN=5" in out
    assert "  NMEAS=3" in out
    assert "  NPTS=1024" in out
    assert "  MAXBLKS=2" in out
    assert "  DT=0.125" in out
    assert "  WINDOW=HANN" in out
    assert "// 3" in out
    assert "\n     HX\n" in out
    assert "\n     HY\n" in out
    assert "\n     EX\n" in out


def test_tsect_missing_section_raises(tmp_path: Path):
    edi = textwrap.dedent(
        """
        >HEAD
          DATAID="X"
        >END
        """
    ).strip("\n")
    p = _write(tmp_path, "no_ts.edi", edi)
    with pytest.raises(EdIDataError):
        TSect.from_file(str(p))


def test_tsio_from_file_parses_blocks(tmp_path: Path):
    edi = textwrap.dedent(
        """
        >HEAD
          DATAID="TS2"

        >=TSERIESSECT
          SECTID=S
          NCHAN=2

        >TSERIES ID=HX NPTS=4 DT=0.25 // 4
          1 2
          3 4

        >TSERIES ID=HY // 3
          0.1 0.2 0.3
        >END
        """
    ).strip("\n")
    p = _write(tmp_path, "ts_blocks.edi", edi)

    # parse header to get body start
    sect = TSect.from_file(str(p))
    io = TSIO.from_file(
        str(p), start_line=sect.start_data_lines_num
    )
    assert len(io.blocks) == 2

    b0, b1 = io.blocks
    # first block options
    assert b0.id == "HX"
    assert b0.npts == 4
    assert b0.dt == 0.25
    assert b0.nvals_hint == 4
    assert b0.values == [1.0, 2.0, 3.0, 4.0]

    # second block minimal header + values
    assert b1.id == "HY"
    assert b1.npts is None
    assert b1.dt is None
    assert b1.nvals_hint == 3
    assert b1.values == [0.1, 0.2, 0.3]


def test_tsio_parse_tolerates_bad_tokens(tmp_path: Path):
    edi = textwrap.dedent(
        """
        >HEAD
          DATAID="TS3"
        >=TSERIESSECT
          NCHAN=1

        >TSERIES ID=HX // 3
          bad 2.0 3.0
        >END
        """
    ).strip("\n")
    p = _write(tmp_path, "ts_bad.edi", edi)
    sect = TSect.from_file(str(p))
    io = TSIO.from_file(
        str(p), start_line=sect.start_data_lines_num
    )
    assert len(io.blocks) == 1
    # 'bad' ignored, numeric kept
    assert io.blocks[0].values == [2.0, 3.0]


def test_tsio_write_roundtrip_formatting(tmp_path: Path):
    # build blocks programmatically and serialize
    io = TSIO()

    from pycsamt.seg.time_series import _TSBlock

    b1 = _TSBlock()
    b1.options["id"] = "HX"
    b1.options["npts"] = 4
    b1.options["dt"] = 0.125
    b1.options["foo"] = "bar"
    b1.options["a"] = 1
    b1.values = [0.1, 0.2, 0.3, 0.4]
    io.blocks.append(b1)

    b2 = _TSBlock()
    b2.options["id"] = "HY"
    b2.values = [1.0]
    io.blocks.append(b2)

    lines = io.write(per_line=3, float_fmt="{: .3E}")
    text = "".join(lines)

    # canonical keys first (ID, NPTS, DT), then extras sorted
    assert ">TSERIES ID=HX NPTS=4 DT=0.125" in text
    assert " A=1 FOO=bar" in text
    # values chunked per_line and formatted
    assert "  1.000E-01  2.000E-01  3.000E-01" in text
    assert "  4.000E-01" in text
    # second block minimal
    assert ">TSERIES ID=HY" in text
    assert "  1.000E+00" in text


def test_timeseries_mixin_helpers(tmp_path: Path):
    edi = textwrap.dedent(
        """
        >HEAD
          DATAID="X2"
        >=TSERIESSECT
          SECTID=SITE_X
          NCHAN=1

        >TSERIES ID=HX // 2
          1 2
        >END
        """
    ).strip("\n")
    p = _write(tmp_path, "ts_mixin.edi", edi)

    class Host(TimeSeriesMixin):
        pass

    header = Host.read_tseries_header(str(p))
    blocks = Host.read_tseries_blocks(str(p))
    assert isinstance(header, TSect)
    assert isinstance(blocks, TSIO)
    assert len(blocks.blocks) == 1

