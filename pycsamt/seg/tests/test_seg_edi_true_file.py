# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0-or-later

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from pycsamt.exceptions import EdIDataError
from pycsamt.seg.edi import EDIFile


def _has(ed: EDIFile, key: str) -> bool:
    return ed.get_section(key) is not None


def _non_empty(a: np.ndarray | None) -> bool:
    if a is None:
        return False
    try:
        return a.size > 0
    except Exception:
        return False


def test_read_any_edi(any_edi_file: Path) -> None:
    ed = EDIFile(any_edi_file)

    # basic sections
    assert _has(ed, "head")
    # mtsect may not exist for spectra-only EDI
    # (but your real files include it)
    mt = ed.has_section("mtsect")

    # Z container integrity when MT/EMAP present
    if mt:
        assert _non_empty(ed.Z.freq)
        assert ed.Z.z is not None
        assert ed.Z.z.shape[-2:] == (2, 2)

        # tipper arrays are optional; if present, keep shape
        tip = getattr(ed.Tip, "tipper", None)
        if tip is not None and tip.size:
            assert tip.shape[-2:] == (1, 2)

    # spectra high-level container if present
    if _has(ed, "spectra") or _has(ed, "spectra_sect"):
        spec = ed.get_section("spectra")
        # we expose a high-level Spectra object
        assert spec is not None
        # must be serializable back to IO
        assert hasattr(spec, "to_io")

    # time series high-level container if present
    if _has(ed, "timeseries") or _has(ed, "timeseries_sect"):
        ts = ed.get_section("timeseries")
        assert ts is not None
        # API surface used downstream
        assert hasattr(ts, "channels")
        assert hasattr(ts, "to_io")


def test_write_roundtrip_generic(any_edi_file: Path, tmp_path: Path) -> None:
    ed = EDIFile(any_edi_file)

    # write with auto-detected datatype, custom dir
    out = ed.write(savepath=tmp_path)
    outp = Path(out)
    assert outp.exists()
    # end sentinel present
    assert ">END" in outp.read_text(encoding="utf-8")

    # re-open written file to ensure readability
    ed2 = EDIFile(outp)
    assert _has(ed2, "head")
    # if MT/EMAP present, frequencies must be recovered
    if ed.has_section("mtsect"):
        assert _non_empty(ed2.Z.freq)


def test_interpolate_and_write_new(edi_imp_file: Path, tmp_path: Path) -> None:
    ed = EDIFile(edi_imp_file)
    if not ed.has_section("mtsect"):
        pytest.skip("no MT/EMAP in this sample")

    f = np.asarray(ed.Z.freq, float)
    assert f.size >= 4

    # new grid strictly inside the original span
    fnew = np.geomspace(f.min() * 1.05, f.max() * 0.95, 8)

    znew = ed.interpolate(fnew, kind="linear", bounds_error=True)

    # write a fresh file using the new Z
    out = ed.write_new_edi(
        edi_fn="interp.edi",
        Z=znew,
        savepath=tmp_path,
    )
    p = Path(out)
    assert p.exists()
    ed3 = EDIFile(p)
    got = np.asarray(ed3.Z.freq, float)
    assert np.allclose(
        np.sort(np.around(got, 2)),
        np.sort(np.around(fnew, 2)),
        atol=0.0,
    )

    # got = np.asarray(ed3.Z.freq, float)
    # exp = np.around(fnew, 2)[::-1]  # reader stores descending
    # assert np.allclose(got, exp, atol=0.0)


def test_repr_str_do_not_crash(any_edi_file: Path) -> None:
    ed = EDIFile(any_edi_file)
    # smoke test repr/str; no exceptions, non-empty str
    s = str(ed)
    r = repr(ed)
    assert isinstance(s, str) and len(s) > 0
    assert isinstance(r, str) and len(r) > 0


def test_spectra_roundtrip_if_present(edi_spe_file: Path, tmp_path: Path) -> None:
    ed = EDIFile(edi_spe_file)
    if not (_has(ed, "spectra") or _has(ed, "spectra_sect")):
        pytest.skip("no spectra section in sample")

    # force a spectra-only write path if mtsect is absent;
    # otherwise a full MT+spectra write will occur.
    out = ed.write(savepath=tmp_path)
    p = Path(out)
    txt = p.read_text(encoding="utf-8")

    # ensure spectra content is serialized back
    # check for either header or data block token
    assert (">=SPECTRASECT" in txt) or (">SPECTRA" in txt)


def test_timeseries_presence_if_any(edi_csamt_file: Path, tmp_path: Path) -> None:
    ed = EDIFile(edi_csamt_file)
    ts = ed.get_section("timeseries")
    # not all files will have TS; this is a soft check
    if ts is None:
        pytest.skip("no time series in sample")

    # ensure serialization calls succeed
    sect, io = ts.to_io()
    assert sect is not None and io is not None

    out = ed.write(savepath=tmp_path)
    p = Path(out)
    txt = p.read_text(encoding="utf-8")

    # look for either the header or a block tag
    assert (">=TSERIESSECT" in txt) or (">TSERIES" in txt)


def test_bounds_error_on_bad_interp(any_edi_file: Path) -> None:
    ed = EDIFile(any_edi_file)
    if not ed.has_section("mtsect"):
        pytest.skip("no MT/EMAP in this sample")

    f = np.asarray(ed.Z.freq, float)
    # request outside the original span → should raise
    bad = np.array([f.min() * 0.5, f.max() * 2.0])
    with pytest.raises((ValueError, EdIDataError)):
        ed.interpolate(bad, bounds_error=True)
