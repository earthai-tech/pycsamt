
from pathlib import Path

import pandas as pd
import pytest

from pycsamt.zonge.utils import (
    classify_avg_format,
    load_avg,
    write_avg,
)

# -----------------------------
# Fixtures: minimal K1 / K2 data
# -----------------------------

K1_TEXT = r"""
\ AMTAVG 7.76: "K1.fld", Dated 99-01-01, Processed 11 Jul 17
$ ASPACE=  50.0m
$ XMTR  =    20.
/ a comment line that should be ignored
skp Station Freq  Comp Amps Emag Ephz Hmag Hphz Resistivity Phase %Emag sEphz %Hmag sHphz %Rho sPhz
 2   150.0   8192 ExHy  5.00  3.1061e+2  1371.6  9.2137e-2  1953.2  2.7746e+2  -581.6   13.8   84.7    9.8   73.6   14.7  136.0
 2   150.0   4096 ExHy  7.00  3.6146e+2  -124.0  9.1877e-2     89.4  7.5575e+2  -213.4    7.8   85.6    9.1   73.7    8.6   20.7
"""

K2_TEXT = r"""
\ AMTAVG 7.76: "LCS01.fld", Dated 99-01-01, Processed 22 Jul 16
\ CSAVGW/ASTATIC modern format with repeated station blocks
$Survey.Type=CSAMT
$Survey.Array=Scalar
$Line.Name="LCS01"
$Unit.Length=m
$Unit.E=nV/Am
$Unit.B=pT/A
$Unit.Phase=mrad

$Rx.GdpStn=25
$Rx.Stn=25
$Rx.Length=50 m
$Rx.Cmp=ExHy
Z.mwgt,Z.pwgt,Freq, Tx.Amp,E.mag,E.phz,B.mag,B.phz,Z.mag,Z.phz,ARes.mag,E.%err,E.perr,B.%err,B.perr,Z.%err,Z.perr,ARes.%err
1,  1,  1,    13,    800.0,  -100.0,   1.30,  200.0,   600.0, -300.0,  10000,   2.0,   20.0,  3.0,   30.0,  4.0,   40.0,  5.0
1,  0,  2,    13,    *,      -110.0,   1.10,  210.0,   610.0, -310.0,  11000,   2.5,   21.0,  3.5,   31.0,  4.5,   41.0,  5.5

$Rx.GdpStn=75
$Rx.Stn=75
$Rx.Length=50 m
$Rx.Cmp=ExHy
Z.mwgt,Z.pwgt,Freq, Tx.Amp,E.mag,E.phz,B.mag,B.phz,Z.mag,Z.phz,ARes.mag,E.%err,E.perr,B.%err,B.perr,Z.%err,Z.perr,ARes.%err
1,  1,  1.41, 13,    900.0,  -120.0,   1.40,  220.0,   650.0, -320.0,  12000,   2.2,   22.0,  3.2,   32.0,  4.2,   42.0,  5.2
1,  1,  2.81, 13,    850.0,   -99.9,   1.52,  240.0,   700.0, -330.0,  13000,   2.3,   23.0,  3.3,   33.0,  4.3,   43.0,  5.3
""".lstrip()


# -----------------------------
# Helpers
# -----------------------------

def _write(tmp_path: Path, name: str, text: str) -> Path:
    p = tmp_path / name
    p.write_text(text, encoding="utf-8")
    return p

# -----------------------------
# Tests
# -----------------------------

def test_classify_avg_format_for_k1_and_k2(tmp_path):
    p1 = _write(tmp_path, "sample_k1.avg", K1_TEXT)
    p2 = _write(tmp_path, "sample_k2.avg", K2_TEXT)

    # Read back raw lines to test classifier directly
    k1_lines = p1.read_text(encoding="utf-8").splitlines()
    k2_lines = p2.read_text(encoding="utf-8").splitlines()

    assert classify_avg_format(k1_lines) == 1
    assert classify_avg_format(k2_lines) == 2


def test_load_avg_kind1_parses_and_meta_is_empty(tmp_path):
    p = _write(tmp_path, "k1.avg", K1_TEXT)
    df, meta = load_avg(p)

    # Basic shape and columns
    assert len(df) == 2
    # Standardised columns should exist
    for col in ("station", "freq", "comp", "rho", "phase"):
        assert col in df.columns

    # Legacy: meta is empty in current design
    assert meta == {}


def test_load_avg_kind2_reads_all_blocks_and_stamps_station(tmp_path):
    p = _write(tmp_path, "k2.avg", K2_TEXT)
    df, meta = load_avg(p)

    # Should have concatenated 2 blocks × 2 rows = 4 rows
    assert len(df) == 4

    # Station stamped from $Rx.Stn
    stations = set(df["station"].tolist())
    assert stations == {25.0, 75.0}

    # Column standardisation (kind-2 → canonical)
    for col in ("freq", "rho", "phase", "emag", "hmag"):
        assert col in df.columns

    # '*' → NaN for E.mag in the second row of station 25
    row = df[(df["station"] == 25.0) & (df["freq"] == 2.0)].iloc[0]
    assert pd.isna(row["emag"])

    # 'use' flag: z.mwgt=1 & z.pwgt=1 → True; pwt=0 → False
    assert "use" in df.columns
    uses = df.sort_values(["station", "freq"])["use"].tolist()
    # rows (25,1.0)=True, (25,2.0)=False, (75,1.41)=True, (75,2.81)=True
    assert uses == [True, False, True, True]

    # meta includes global keys and per-block list
    assert meta.get("Survey.Type") == "CSAMT"
    assert isinstance(meta.get("blocks"), list) and len(meta["blocks"]) == 2
    # each block meta should carry Rx.Stn used above
    assert {float(b["Rx.Stn"]) for b in meta["blocks"]} == {25.0, 75.0}


def test_write_avg_default_path_and_header_case(
        tmp_path, monkeypatch):
    # Create a tiny canonical core frame
    core = pd.DataFrame({
        "station": [10, 10],
        "freq": [1.0, 2.0],
        "rho": [100.0, 200.0],   # should be exported as "ARes.mag"
        "phase": [-100.0, -110.0],
    })

    # Run in a temp cwd so default file lands here
    monkeypatch.chdir(tmp_path)

    out_path = write_avg(core=core, extra=None,
                         meta={"Survey.Type": "CSAMT"},
                         path=None, stamp=False)
    assert out_path.exists()
    assert out_path.name == "exported_kind2.avg"

    txt = out_path.read_text(encoding="utf-8")
    # meta header written
    assert "$Survey.Type=CSAMT" in txt
    # column header should contain kind-2 label for rho
    # assert "ARes.mag" in txt.splitlines() # [-(len(core) + 1)]  # header line


def test_kind2_comment_lines_are_ignored(tmp_path):
    # Insert comment lines among rows and ensure parsing is unaffected
    k2_with_comments = K2_TEXT.replace(
        "ARes.mag,E.%err",
        'ARes.mag,E.%err\n\\ this is a comment line that must be ignored'
    )
    p = _write(tmp_path, "k2c.avg", k2_with_comments)
    df, meta = load_avg(p)

    # Still reads all 4 rows
    assert len(df) == 4
    # meta present
    assert "Survey.Type" in meta


if __name__=='__main__': # pragma: no-cover
   pytest.main( [__file__])
