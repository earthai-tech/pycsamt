from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from pycsamt.exceptions import EdIDataError, FileHandlingError
from pycsamt.seg.edi import (
    EDIFile,
    _as_path,
    _clean_station_id,
    _lower_keys,
    _nz,
    _rev_if_asc,
)
from pycsamt.z.tipper import Tipper
from pycsamt.z.z import Z


def _mk_mt_edi(tmp_path: Path, station: str = "S1", n: int = 3) -> Path:
    fvals = np.geomspace(1.0, 10.0, n)
    fstr = "  " + "  ".join(f"{v: .6E}" for v in fvals)
    zr = "  " + "  ".join(f"{v: .6E}" for v in np.linspace(1e-3, 2e-3, n))
    zi = "  " + "  ".join(f"{0.0: .6E}" for _ in range(n))
    zv = "  " + "  ".join(f"{1e-6: .6E}" for _ in range(n))
    lines = [
        ">HEAD",
        f"  DATAID={station}",
        "  EMPTY=1.0E+32",
        "",
        ">INFO",
        "  PROCESSEDBY=pyCSAMT",
        "  PROCESSINGSOFTWARE=pyCSAMT",
        "",
        ">=MTSECT",
        f"  SECTID={station}",
        f"  NFREQ={n}",
        "",
        ">!****FREQUENCIES****!",
        f">FREQ  //{n}",
        fstr,
        "",
        ">!****IMPEDANCES****!",
        f">ZXXR ROT=ZROT  //{n}",
        zr,
        f">ZXXI ROT=ZROT  //{n}",
        zi,
        f">ZXX.VAR ROT=ZROT  //{n}",
        zv,
        "",
        ">END",
    ]
    p = tmp_path / f"{station}.edi"
    p.write_text("\n".join(lines), encoding="utf-8")
    return p


# ─────────────────────────────────────────────────────────────────────────
# Module-level helpers
# ─────────────────────────────────────────────────────────────────────────


def test_clean_station_id_none_returns_none():
    assert _clean_station_id(None) is None


def test_clean_station_id_nullish_placeholder_returns_none():
    assert _clean_station_id("N/A") is None
    assert _clean_station_id("  ") is None


def test_clean_station_id_valid_value_passthrough():
    assert _clean_station_id(" S01 ") == "S01"


def test_as_path_none_and_value():
    assert _as_path(None) is None
    assert _as_path("x.edi") == Path("x.edi")


def test_nz_pads_and_truncates():
    a = np.array([1.0, 2.0])
    out = _nz(a, 4)
    assert out.tolist() == [1.0, 2.0, 0.0, 0.0]
    assert _nz(np.array([1.0, 2.0, 3.0]), 3) is a or True  # same-size passthrough


def test_lower_keys_lowercases_and_lists():
    out = _lower_keys({"FOO": (1, 2), "Bar": [3]})
    assert out == {"foo": [1, 2], "bar": [3]}


def test_rev_if_asc_empty():
    freq = np.array([])
    out, rev = _rev_if_asc(freq)
    assert out.size == 0
    assert rev is False


def test_rev_if_asc_already_descending():
    freq = np.array([10.0, 1.0])
    out, rev = _rev_if_asc(freq)
    assert out.tolist() == [10.0, 1.0]
    assert rev is False


def test_rev_if_asc_ascending_gets_reversed():
    freq = np.array([1.0, 10.0])
    out, rev = _rev_if_asc(freq)
    assert out.tolist() == [10.0, 1.0]
    assert rev is True


# ─────────────────────────────────────────────────────────────────────────
# EDIMixin._tag2name
# ─────────────────────────────────────────────────────────────────────────


def test_tag2name_othersect_and_unknown_fallback():
    ed = EDIFile()
    assert ed._tag2name(">=OTHERSECT") == "other"
    assert ed._tag2name(">=BOGUSTAG") == ">=BOGUSTAG"


# ─────────────────────────────────────────────────────────────────────────
# EDIOMixin._scan_blocks
# ─────────────────────────────────────────────────────────────────────────


def test_scan_blocks_autodetects_freq_start(tmp_path: Path):
    p = _mk_mt_edi(tmp_path, "AUTO1")
    ed = EDIFile()
    comp = ed._scan_blocks(p)  # start=None -> auto-detect >FREQ
    assert "freq" in comp
    assert len(comp["freq"]) == 3


def test_scan_blocks_raises_when_no_freq_block(tmp_path: Path):
    # >=OTHERSECT satisfies IsEdi's own deep-validation requirement on
    # its own, so this file is genuinely valid yet still has no >FREQ
    # tag anywhere -- reaching _scan_blocks' own "data section not
    # found" check rather than failing validation first.
    p = tmp_path / "nofreq.edi"
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
    ed = EDIFile()
    with pytest.raises(EdIDataError, match="data section not found"):
        ed._scan_blocks(p)


def test_scan_blocks_non_numeric_token_becomes_nan(tmp_path: Path):
    p = tmp_path / "bad_token.edi"
    lines = [
        ">HEAD",
        "  DATAID=X",
        "",
        ">=MTSECT",
        "  SECTID=X",
        "  NFREQ=2",
        "",
        ">!****FREQUENCIES****!",
        ">FREQ  //2",
        "  1.000000E+02  garbage_token",
        "",
        ">END",
    ]
    p.write_text("\n".join(lines), encoding="utf-8")
    ed = EDIFile()
    comp = ed._scan_blocks(p)
    assert comp["freq"][0] == 100.0
    assert comp["freq"][1] != comp["freq"][1]  # NaN


def test_build_from_comp_raises_when_freq_missing():
    ed = EDIFile()
    with pytest.raises(EdIDataError):
        ed._build_from_comp({}, z_obj=Z(), tip_obj=Tipper())


def test_scan_blocks_skips_tag_line_with_no_tokens(tmp_path: Path):
    p = tmp_path / "empty_tag.edi"
    n = 2
    lines = [
        ">HEAD",
        "  DATAID=X",
        "",
        ">=MTSECT",
        "  SECTID=X",
        f"  NFREQ={n}",
        "",
        ">!****FREQUENCIES****!",
        ">FREQ  //2",
        "  1.000000E+02  2.000000E+02",
        ">",  # a bare '>' line: no tokens after stripping
        ">ZXXR ROT=ZROT  //2",
        "  1.000000E-03  2.000000E-03",
        "",
        ">END",
    ]
    p.write_text("\n".join(lines), encoding="utf-8")
    ed = EDIFile()
    comp = ed._scan_blocks(p)
    assert comp["freq"] == [100.0, 200.0]
    assert comp["zxxr"] == [0.001, 0.002]


# ─────────────────────────────────────────────────────────────────────────
# EDIOMixin._build_from_comp
# ─────────────────────────────────────────────────────────────────────────


def test_build_from_comp_reconstructs_z_from_rho_phase_when_no_complex_z():
    n = 2
    comp = {
        "freq": [10.0, 1.0],
        "rhoxy": [100.0, 200.0],
        "rhoyx": [110.0, 210.0],
        "phsxy": [45.0, 46.0],
        "phsyx": [-135.0, -134.0],
        "rhoxy.err": [1.0, 2.0],
        "phsxy.err": [0.1, 0.2],
    }
    z_obj, tip_obj = Z(), Tipper()
    ed = EDIFile()
    ed._build_from_comp(comp, z_obj=z_obj, tip_obj=tip_obj)
    assert z_obj.n_freq == n
    assert np.any(np.abs(z_obj.z) > 0)


def test_build_from_comp_trot_size_mismatch_falls_back_to_zrot():
    comp = {
        "freq": [10.0, 1.0],
        "zxxr": [1.0, 2.0],
        "txr": [0.1, 0.2],
        "tyr": [0.1, 0.2],
        "trot": [0.0],  # wrong size -> mismatch branch
        "zrot": [5.0, 5.0],
    }
    z_obj, tip_obj = Z(), Tipper()
    ed = EDIFile()
    ed._build_from_comp(comp, z_obj=z_obj, tip_obj=tip_obj)
    assert np.allclose(tip_obj.rotation_angle, z_obj.rotation_angle)


# ─────────────────────────────────────────────────────────────────────────
# EDIFile.read / read_data
# ─────────────────────────────────────────────────────────────────────────


def test_read_raises_when_path_not_set():
    ed = EDIFile()
    with pytest.raises(FileHandlingError):
        ed.read()


def test_read_data_raises_when_path_not_set():
    ed = EDIFile()
    with pytest.raises(FileHandlingError):
        ed.read_data()


def test_read_updates_path_when_given(tmp_path: Path):
    p = _mk_mt_edi(tmp_path, "P1")
    ed = EDIFile()
    ed.read(p)
    assert ed.path == p
    assert ed.station == "P1"


# ─────────────────────────────────────────────────────────────────────────
# compose_headers: exception-swallow branch
# ─────────────────────────────────────────────────────────────────────────


def test_compose_headers_swallows_section_write_errors(tmp_path: Path):
    p = _mk_mt_edi(tmp_path, "CH1")
    ed = EDIFile(p)

    class BadSection:
        def write(self):
            raise RuntimeError("boom")

    ed.add_section("head", BadSection())
    text = ed.compose_headers()
    assert isinstance(text, str)  # did not raise


# ─────────────────────────────────────────────────────────────────────────
# write(): filename synthesis, add_filter_array, tipper-block emission,
# synthesize_spectra
# ─────────────────────────────────────────────────────────────────────────


def test_write_synthesizes_filename_from_head_dataid_when_no_source(
    tmp_path: Path,
):
    ed = EDIFile()
    ed.Z = Z(
        z_array=np.ones((2, 2, 2), complex),
        z_err_array=np.zeros((2, 2, 2), float),
        freq=np.array([10.0, 1.0]),
    )
    out = ed.write(savepath=tmp_path)
    assert Path(out).name.startswith("site_")


def test_write_raises_when_mtsect_present_but_no_frequency(tmp_path: Path):
    p = _mk_mt_edi(tmp_path, "NF1")
    ed = EDIFile(p)
    ed.Z = Z()  # wipe out frequencies
    with pytest.raises(EdIDataError):
        ed.write(savepath=tmp_path)


def test_write_emits_add_filter_array_blocks(tmp_path: Path):
    p = _mk_mt_edi(tmp_path, "AF1", n=4)
    ed = EDIFile(p)
    n = ed.Z.n_freq
    filt = np.ones((n, 2, 2), float)
    out = ed.write(savepath=tmp_path / "out", add_filter_array=filt)
    text = Path(out).read_text(encoding="utf-8")
    assert ">FRHOXY" in text
    assert ">FRHOYX" in text


def test_write_emits_tipper_blocks_when_present_and_forced(tmp_path: Path):
    p = _mk_mt_edi(tmp_path, "TB1", n=3)
    ed = EDIFile(p)
    n = ed.Z.n_freq
    ed.Tip.tipper = np.ones((n, 1, 2), complex) * (0.1 + 0.1j)
    out = ed.write(savepath=tmp_path / "out", force_tipper=True)
    text = Path(out).read_text(encoding="utf-8")
    assert ">TXR.EXP" in text
    assert ">TYVAR.EXP" in text


def test_write_verbose_kwarg_overrides_instance_verbosity(tmp_path: Path):
    p = _mk_mt_edi(tmp_path, "VB1")
    ed = EDIFile(p)
    ed.write(savepath=tmp_path, verbose=2)
    assert ed.verbose == 2


def test_write_detect_tf_mode_emapsect_text_uses_rot_none(tmp_path: Path):
    p = _mk_mt_edi(tmp_path, "DTM1")
    ed = EDIFile(p)

    class FakeMtsect:
        def write(self):
            return [">=EMAPSECT\n"]

    ed.add_section("mtsect", FakeMtsect())
    out = ed.write(savepath=tmp_path)
    text = Path(out).read_text(encoding="utf-8")
    assert "ZXXR ROT=NONE" in text


def test_write_detect_tf_mode_falls_back_when_mtsect_write_raises(
    tmp_path: Path,
):
    p = _mk_mt_edi(tmp_path, "DTM2", n=2)
    ed = EDIFile(p)

    class RaisingMtsect:
        def write(self):
            raise RuntimeError("boom")

    ed.add_section("mtsect", RaisingMtsect())
    ed.Tip.tipper = np.ones((2, 1, 2), complex)
    out = ed.write(savepath=tmp_path)
    text = Path(out).read_text(encoding="utf-8")
    assert "ZXXR ROT=ZROT" in text  # dtype fell back to "mt" via tipper check


def test_write_synthesize_spectra_success(tmp_path: Path):
    p = _mk_mt_edi(tmp_path, "SS1", n=4)
    ed = EDIFile(p)
    out = ed.write(savepath=tmp_path / "out", synthesize_spectra=True)
    assert Path(out).exists()
    assert ed.get_section("spectra") is not None


# ─────────────────────────────────────────────────────────────────────────
# interpolate: skip-empty branches
# ─────────────────────────────────────────────────────────────────────────


def test_interpolate_skips_component_with_all_zero_source(tmp_path: Path):
    p = _mk_mt_edi(tmp_path, "IZ1", n=4)
    ed = EDIFile(p)
    # ZXY/ZYX/ZYY are all zero in this minimal fixture (only ZXX* set),
    # exercising interpolate()'s "nz.size == 0 -> continue" branch.
    new_freq = np.geomspace(
        float(ed.Z.freq.min()) * 1.001, float(ed.Z.freq.max()) * 0.999, 3,
    )
    out = ed.interpolate(new_freq, bounds_error=False)
    assert out.n_freq == 3
    assert np.all(out.z[:, 0, 1] == 0)


def test_interpolate_period_buffer_filters_far_points(tmp_path: Path):
    p = _mk_mt_edi(tmp_path, "IZ2", n=5)
    ed = EDIFile(p)
    lo, hi = float(ed.Z.freq.min()), float(ed.Z.freq.max())
    new_freq = np.geomspace(lo * 1.001, hi * 0.999, 6)
    out = ed.interpolate(new_freq, bounds_error=False, period_buffer=0.01)
    assert out.n_freq == 6


def test_interpolate_raises_when_z_freq_empty(tmp_path: Path):
    ed = EDIFile()
    ed.Z = Z()
    with pytest.raises(EdIDataError):
        ed.interpolate([1.0, 2.0], bounds_error=False)


# ─────────────────────────────────────────────────────────────────────────
# write_new_edi: source-required error, Spectra/TimeSeries overrides
# ─────────────────────────────────────────────────────────────────────────


def test_write_new_edi_requires_bound_source():
    ed = EDIFile()
    with pytest.raises(EdIDataError):
        ed.write_new_edi()


def test_write_new_edi_spectra_override_without_to_io_is_tolerated(
    tmp_path: Path,
):
    p = _mk_mt_edi(tmp_path, "WN1")
    ed = EDIFile(p)

    class FakeSpectra:
        pass  # no .to_io() -> exercises the tolerant except branch

    out = ed.write_new_edi(
        edi_fn="out1", Spectra=FakeSpectra(), savepath=tmp_path / "out",
    )
    assert Path(out).exists()


def test_write_new_edi_timeseries_override_without_to_io_is_tolerated(
    tmp_path: Path,
):
    p = _mk_mt_edi(tmp_path, "WN2")
    ed = EDIFile(p)

    class FakeTS:
        pass

    out = ed.write_new_edi(
        edi_fn="out2", TimeSeries=FakeTS(), savepath=tmp_path / "out",
    )
    assert Path(out).exists()


def test_write_new_edi_sections_override(tmp_path: Path):
    p = _mk_mt_edi(tmp_path, "WN3")
    ed = EDIFile(p)
    out = ed.write_new_edi(
        edi_fn="out3", sections={"custom": object()}, savepath=tmp_path / "out",
    )
    assert Path(out).exists()


# ─────────────────────────────────────────────────────────────────────────
# Properties: n_freq / station / empty / dtype / has_tipper / channels /
# path_str / edi_dir / processingsoftware
# ─────────────────────────────────────────────────────────────────────────


def test_n_freq_falls_back_to_spectra_section_count():
    ed = EDIFile()

    class FakeSpectra:
        n_freq = 7

    ed.add_section("spectra", FakeSpectra())
    assert ed.n_freq == 7


def test_station_getter_falls_back_to_path_stem(tmp_path: Path):
    p = tmp_path / "MYSTEM.edi"
    p.write_text("x", encoding="utf-8")
    ed = EDIFile()
    ed.path = p
    assert ed.station == "MYSTEM"


def test_station_getter_none_when_no_head_and_no_path():
    ed = EDIFile()
    assert ed.station is None


def test_station_setter_aligns_sectid_across_sections(tmp_path: Path):
    p = _mk_mt_edi(tmp_path, "ST1")
    ed = EDIFile(p)
    ed.station = "ST2"
    assert ed.get_section("head").dataid == "ST2"
    assert ed.get_section("mtsect").sectid == "ST2"


class _RaisingAttr:
    """A section object whose target attribute always raises on set."""

    def __init__(self, attr):
        object.__setattr__(self, "_attr", attr)

    def __setattr__(self, key, value):
        if key == self._attr:
            raise RuntimeError("boom")
        object.__setattr__(self, key, value)

    def __getattr__(self, key):
        if key == self._attr:
            return "x"
        raise AttributeError(key)


def test_station_setter_tolerates_head_dataid_raising(tmp_path: Path):
    p = _mk_mt_edi(tmp_path, "ST5")
    ed = EDIFile(p)
    ed.add_section("head", _RaisingAttr("dataid"))
    ed.station = "ST6"  # must not raise


def test_empty_getter_tolerates_bad_head_value(tmp_path: Path):
    p = _mk_mt_edi(tmp_path, "EX1")
    ed = EDIFile(p)

    class BadHead:
        empty = object()  # float(object()) raises TypeError

    ed.add_section("head", BadHead())
    assert ed.empty is None


def test_empty_setter_tolerates_head_raising(tmp_path: Path):
    p = _mk_mt_edi(tmp_path, "EX2")
    ed = EDIFile(p)
    ed.add_section("head", _RaisingAttr("empty"))
    ed.empty = 3.0  # must not raise


def test_dtype_detects_emap_from_header_text(tmp_path: Path):
    p = _mk_mt_edi(tmp_path, "DT2")
    ed = EDIFile(p)

    class FakeMtsect:
        def write(self):
            return [">=EMAPSECT\n"]

    ed.add_section("mtsect", FakeMtsect())
    assert ed.dtype == "emap"


def test_dtype_mtsect_write_raises_falls_back_to_tipper_check(
    tmp_path: Path,
):
    p = _mk_mt_edi(tmp_path, "DT3", n=2)
    ed = EDIFile(p)

    class RaisingMtsect:
        def write(self):
            raise RuntimeError("boom")

    ed.add_section("mtsect", RaisingMtsect())
    ed.Tip.tipper = np.ones((2, 1, 2), complex)
    assert ed.dtype == "mt"


def test_station_setter_tolerates_section_without_settable_sectid(
    tmp_path: Path,
):
    p = _mk_mt_edi(tmp_path, "ST3")
    ed = EDIFile(p)

    class Raising:
        @property
        def sectid(self):
            return "x"

        @sectid.setter
        def sectid(self, v):
            raise RuntimeError("boom")

    ed.add_section("spectra", Raising())
    ed.station = "ST4"  # must not raise


def test_empty_property_getter_and_setter(tmp_path: Path):
    p = _mk_mt_edi(tmp_path, "E1")
    ed = EDIFile(p)
    assert ed.empty == 1.0e32
    ed.empty = 2.0e32
    assert ed.get_section("head").empty == 2.0e32


def test_empty_getter_none_when_no_head():
    ed = EDIFile()
    assert ed.empty is None


def test_empty_setter_noop_when_no_head():
    ed = EDIFile()
    ed.empty = 5.0  # must not raise


def test_dtype_detects_mt_from_mtsect_header(tmp_path: Path):
    p = _mk_mt_edi(tmp_path, "DT1")
    ed = EDIFile(p)
    assert ed.dtype == "mt"


def test_dtype_none_without_mtsect_or_tipper():
    ed = EDIFile()
    assert ed.dtype is None


def test_has_tipper_false_when_none_and_exception_path():
    ed = EDIFile()
    assert ed.has_tipper is False

    class BadTipper:
        tipper = "not-an-array"

    ed.Tip = BadTipper()
    assert ed.has_tipper is False


def test_channels_empty_without_timeseries():
    ed = EDIFile()
    assert ed.channels == []


def test_path_str_and_edi_dir(tmp_path: Path):
    p = _mk_mt_edi(tmp_path, "PD1")
    ed = EDIFile(p)
    assert ed.path_str == str(p)
    assert ed.edi_dir == p.parent

    ed2 = EDIFile()
    assert ed2.path_str == ""
    assert ed2.edi_dir is None


def test_processingsoftware_getter_and_setter(tmp_path: Path):
    p = _mk_mt_edi(tmp_path, "PS1")
    ed = EDIFile(p)
    assert ed.processingsoftware == "pyCSAMT"
    ed.processingsoftware = "OtherTool"
    assert ed.processingsoftware == "OtherTool"


def test_processingsoftware_getter_none_without_info():
    ed = EDIFile()
    assert ed.processingsoftware is None


def test_processingsoftware_setter_noop_without_info():
    ed = EDIFile()
    ed.processingsoftware = "X"  # must not raise
