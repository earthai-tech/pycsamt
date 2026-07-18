from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from pycsamt.seg.edi import EDIFile
from pycsamt.site import edit as ed
from pycsamt.site.base import Site, Sites


def _load_edi(p: Path) -> EDIFile:
    return EDIFile(p)  # type: ignore


def _dup_edi(tmp_path: Path, src: Path, stem: str) -> Path:
    dst = tmp_path / f"{stem}.edi"
    dst.write_text(
        src.read_text(encoding="utf-8"),
        encoding="utf-8",
    )
    return dst


def _mk_two_edifiles(
    tmp_path: Path,
    simulated_edi: Path,
    s1: str = "S01",
    s2: str = "S02",
) -> tuple[EDIFile, EDIFile]:
    p1 = _dup_edi(tmp_path, simulated_edi, s1)
    p2 = _dup_edi(tmp_path, simulated_edi, s2)
    return _load_edi(p1), _load_edi(p2)


def test_rotate_inplace_and_copy(simulated_edi: Path) -> None:
    edf = _load_edi(simulated_edi)

    out = ed.rotate(edf, 45.0, inplace=False)
    assert out is not edf

    out2 = ed.rotate(edf, 30.0, inplace=True)
    assert out2 is edf  # mutated in place


def test_select_freq_by_keep_and_range(simulated_edi: Path) -> None:
    edf = _load_edi(simulated_edi)

    f0 = ed.get_freq(edf)
    if f0.size < 2:
        pytest.skip("simulated EDI has <2 freq")

    # keep first row by index
    out = ed.select_freq(edf, keep=[0], inplace=False)
    f1 = ed.get_freq(out)
    assert f1.size == 1

    # range selection (e.g., keep only top half of freq)
    fmin = float(np.median(f0))
    out2 = ed.select_freq(edf, fmin=fmin, inplace=False)
    f2 = ed.get_freq(out2)
    assert 1 <= f2.size <= f0.size
    # the retained rows must actually satisfy the requested range,
    # regardless of whether Z stores freq ascending or descending.
    assert np.all(f2 >= fmin)
    assert out2.Z.z.shape[0] == f2.size


def test_select_freq_keeps_z_aligned_with_descending_native_order(
    simulated_edi: Path,
) -> None:
    edf = _load_edi(simulated_edi)

    native_freq = np.asarray(edf.Z.freq, dtype=float)
    if native_freq.size < 2 or not np.all(np.diff(native_freq) < 0):
        pytest.skip("fixture Z.freq is not natively descending")

    fmin, fmax = 1.0, float(native_freq[len(native_freq) // 2])
    out = ed.select_freq(edf, fmin=fmin, fmax=fmax, inplace=False)

    # rows kept in Z's own (descending) order must match the same
    # band, not the first/last N rows by position.
    assert np.all(out.Z.freq >= fmin)
    assert np.all(out.Z.freq <= fmax)


def test_rename_explicit_and_policy(simulated_edi: Path) -> None:
    edf = _load_edi(simulated_edi)

    # explicit new name
    out = ed.rename(edf, name="NEW01", inplace=False)
    h = out.get_section("head")  # type: ignore
    assert str(h.dataid) == "NEW01"

    # policy-based rename
    def pol(n: str) -> str:
        return f"X_{n}"

    out2 = ed.rename(out, policy=pol, inplace=False)
    h2 = out2.get_section("head")  # type: ignore
    assert str(h2.dataid).startswith("X_")


def test_set_coords_single(simulated_edi: Path) -> None:
    edf = _load_edi(simulated_edi)

    out = ed.set_coords(
        edf,
        lat=11.25,
        lon=22.5,
        elev=333.0,
        inplace=False,
    )
    h = out.get_section("head")  # type: ignore
    got = (float(h.lat), float(h.lon), float(h.elev))
    assert got == (11.25, 22.5, 333.0)


def test_fill_missing_and_recompute(simulated_edi: Path) -> None:
    edf = _load_edi(simulated_edi)

    # Ensure Z arrays exist and are finite after filling
    out = ed.fill_missing(
        edf,
        how="zero",
        components=("Z",),
        inplace=False,
    )

    Z = getattr(out, "Z", None)
    if Z is None:
        pytest.skip("No Z section present after filling")

    zz = getattr(Z, "_z", None)
    assert zz is not None
    arr = np.asarray(zz)
    assert arr.shape[1:] == (2, 2)
    assert np.all(np.isfinite(arr))

    # Recompute rho/phi should not raise
    _ = ed.recompute_res_phase(out, inplace=True)
    # Accept best-effort: resistivity/phase may or may not exist
    # depending on downstream stack. We only assert no exception.


def test_rotate_all_and_select_freq_all(
    tmp_path: Path, simulated_edi: Path
) -> None:
    e1, e2 = _mk_two_edifiles(tmp_path, simulated_edi, "R01", "R02")
    src = [e1, e2]

    # rotate_all returns a Sites wrapper with two items
    ro = ed.rotate_all(src, 15.0, inplace=False)
    assert isinstance(ro, Sites)
    assert len(ro) == 2

    # Now select a subset of frequencies across all
    # (use median to keep at least one row)
    f = ed.get_freq(e1)
    if f.size < 2:
        pytest.skip("simulated EDI has <2 freq")
    fmin = float(np.median(f))

    sl = ed.select_freq_all(src, fmin=fmin, inplace=False)
    assert isinstance(sl, Sites)
    assert len(sl) == 2

    # Inspect first site frequency count decreased or stayed same
    s0 = sl.by_index(0)
    f_after = Site(s0.edi).freq
    assert f_after is not None
    assert 1 <= len(f_after) <= len(f)


def test_rename_all_and_set_coords_all(
    tmp_path: Path, simulated_edi: Path
) -> None:
    e1, e2 = _mk_two_edifiles(tmp_path, simulated_edi, "N01", "N02")
    src = [e1, e2]

    rn = ed.rename_all(
        src,
        name_fn=lambda edi: "X_{}".format(
            Path(getattr(edi, "path", getattr(edi, "file", ""))).stem
            or getattr(edi, "name", "site")
        ),
        inplace=False,
    )
    assert isinstance(rn, Sites)
    assert len(rn) == 2
    assert rn.by_index(0).name.startswith("X_")

    # Provide a mapping by current names to lat/lon/elev
    # Both stems carry same initial header name in the
    # simulated file; map the post-rename value seen above.
    n0 = rn.by_index(0).name
    n1 = rn.by_index(1).name
    coords = {
        n0: (10.0, 20.0, 100.0),
        n1: (11.0, 21.0, 110.0),
    }

    out = ed.set_coords_all(rn.as_list(), coords, inplace=False)
    assert isinstance(out, Sites)
    s0 = out.by_index(0)
    s1 = out.by_index(1)
    assert tuple(map(float, s0.coords)) == (10.0, 20.0, 100.0)
    assert tuple(map(float, s1.coords)) == (11.0, 21.0, 110.0)
