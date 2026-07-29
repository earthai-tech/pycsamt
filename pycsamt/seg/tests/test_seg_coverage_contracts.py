from __future__ import annotations

import io
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pytest

from pycsamt.exceptions import FileHandlingError
from pycsamt.seg import base as base_module
from pycsamt.seg import schema
from pycsamt.seg.base import Base, EDIComponentBase, EdiFileBase, SurveyBase
from pycsamt.seg.collection import EDICollection
from pycsamt.seg.config import (
    SEGConfig,
    _Facade,
    _import_first,
    discover_mixins,
)
from pycsamt.seg.edi import EDIFile
from pycsamt.seg.survey import EDIProfile, Stations, Topography
from pycsamt.seg.utils import (
    _dms_to_deg,
    _format_block_numbers,
    _format_kv,
    _haversine,
    _len,
    _quote,
    _safe_number,
    gather_measurement_key_value_with_str_parser,
    minimum_parser_to_write_edi,
    parse_kv_pairs,
    quick_edi_stats,
    show_edi_stats,
    sort_edis_by_location,
    strip_item,
)
from pycsamt.seg.validation import (
    IsEdi,
    _count_freq_values,
    _expected_nfreq_from_header,
    _extract_tag,
    _extract_tag_in,
    _has_any_data_block,
    _is_tag,
    _iter_blocks,
    _norm_str,
    _split_comment,
    _strip_norm,
    _to_float_or_none,
    _to_int_or_none,
)


@dataclass
class DemoComponent(EDIComponentBase):
    _section = "DEMO"
    name: str | None = None
    count: int = 0
    ratio: float = 0.0
    enabled: bool = False
    values: list[float] | None = None


class MemoryEdi(EdiFileBase):
    def __init__(self, text: str = ">END", **kwargs):
        self.text = text
        super().__init__(**kwargs)

    def read(self):
        return self

    def compose(self):
        return self.text


def test_schema_and_config_public_contracts() -> None:
    assert schema.normalize_tag("  zxxr ") == ">ZXXR"
    assert schema.normalize_tag(">=mtsect") == ">=MTSECT"
    assert schema.tag_family(">ZXYR") == "impedance"
    assert schema.tag_family("not-a-tag") is None
    assert ">HEAD" in schema.ALL_KEYWORDS
    assert "DATAID" in schema.HEAD_ALLOWED

    cfg = SEGConfig()
    assert cfg.lint_line == 62 and cfg.default_empty == 1e32
    mixins = discover_mixins()
    assert mixins and len(mixins) == len(set(mixins))
    assert _import_first("base", ["Missing", "Base"]).__name__ == "Base"
    assert _import_first("base", ["Missing"]).__name__ == "_NoOp"


def test_facade_delegation_and_fallbacks(tmp_path: Path) -> None:
    class Implements:
        @classmethod
        def from_file(cls, path, **kw):
            return (path, kw)

        @classmethod
        def loads(cls, text, **kw):
            return (text, kw)

        def to_file(self, path, **kw):
            Path(path).write_text(kw.get("text", ""), encoding="utf-8")
            return path

        def dumps(self, **kw):
            return kw.get("text", "dump")

        def asdict(self):
            return {"mixed": True}

    Facade = type("Facade", (Implements, _Facade), {"__mixins__": (Implements,)})
    assert Facade.from_file("x", strict=True) == ("x", {"strict": True})
    assert Facade.loads("raw", mode=1) == ("raw", {"mode": 1})
    obj = Facade()
    out = tmp_path / "facade.txt"
    assert obj.to_file(out, text="ok") == out
    assert out.read_text(encoding="utf-8") == "ok"
    assert obj.dumps(text="yes") == "yes"
    assert obj.asdict() == {"mixed": True}
    assert Facade.describe_mixins() == ("Implements",)

    class Bare(_Facade):
        pass

    with pytest.raises(NotImplementedError):
        Bare.from_file("x")
    with pytest.raises(NotImplementedError):
        Bare.loads("x")
    with pytest.raises(NotImplementedError):
        Bare().to_file(out)
    with pytest.raises(NotImplementedError):
        Bare().dumps()
    assert Bare().asdict() == {}


def test_component_serialization_formatting_and_clone() -> None:
    c = DemoComponent(
        name='two "words"',
        count=3,
        ratio=float("nan"),
        enabled=True,
        values=[1.5, 2.0, 3.0, 4.0, 5.0],
    )
    assert c.to_dict()["count"] == 3
    assert 'name="two words"' in c.to_text()
    lines = c.to_lines(only=["count", "enabled"], indent=4)
    assert lines == ["    count=3", "    enabled=1"]
    assert c._format_value(None) == "None"
    assert c._format_value(False) == "0"
    assert c._format_value(float("nan")) == "NaN"
    assert c._format_value((1, "a b")) == '1 "a b"'
    assert "…" in c._repr_value(c.values)

    clone = c.clone(count=9)
    assert clone.count == 9 and clone.name == c.name
    assert c.update(count=4) is c and c.count == 4
    rebuilt = DemoComponent.from_dict({"name": "rebuilt", "extra": 7})
    assert rebuilt.name == "rebuilt" and rebuilt.extra == 7


def test_base_cooperative_init_and_logger_fallback(monkeypatch) -> None:
    class RejectsInit:
        def __init__(self, *args, **kwargs):
            raise RuntimeError("not cooperative")

    class Tolerant(Base, RejectsInit):
        pass

    assert Tolerant("ignored", verbose=True).verbose == 1
    monkeypatch.setattr(base_module, "_get_logger", None)
    logger = Base._logger_factory("seg.coverage")
    assert logger.name == "seg.coverage"
    assert logger.handlers


def test_edifile_base_registry_io_and_helpers(tmp_path: Path, monkeypatch) -> None:
    ed = MemoryEdi("first\nsecond\n", strict_validate=False)
    head = DemoComponent(name="S01")
    head.dataid, head.lat, head.long, head.elev = "S01", 1.25, 2.5, 30.0
    ed.add_section(">HEAD", head)
    assert ed.get_section(" head ") is head
    assert (ed.dataid, ed.lat, ed.lon, ed.elev) == ("S01", 1.25, 2.5, 30.0)
    assert list(ed.iter_sections()) == [("head", head)]
    assert "dataid=S01" in ed.summary()
    assert ed.compose_to_string() == "first\nsecond\n"

    out = ed.write_file(tmp_path / "nested" / "x.edi")
    assert out.read_text(encoding="utf-8") == "first\nsecond\n"
    with pytest.raises(FileExistsError):
        ed.write_file(out, overwrite=False)
    monkeypatch.chdir(tmp_path)
    assert ed._default_out() == Path("S01.edi")

    assert ed.normalize_section_title(" >freq") == "FREQ"
    assert ed.format_block_header(">freq").startswith(">!****FREQ")
    assert ed.format_section_head("freq") == ">FREQ\n"
    assert ed.format_kv("id", "a b", quote=True, width=5) == '  ID="a b"\n'
    block = ed.format_data_block("freq", [1, 2, 3], per_line=2)
    assert "//3" in "".join(block)
    assert len(ed.format_data_block("empty", [], header_comment=False)) == 1
    assert ed.ensure_descending_frequency([]) == ([], None)
    assert ed.ensure_descending_frequency([3, 2]) == ([3, 2], None)
    assert ed.ensure_descending_frequency([1, 2, 3]) == ([3, 2, 1], [2, 1, 0])

    payload = ed.to_dict()
    copy = MemoryEdi.from_dict(payload)
    assert copy.get_section("head").name == "S01"
    ed.remove_section("head")
    assert ed.get_section("head") is None


def test_survey_base_summary_fallbacks_and_table() -> None:
    class Survey(SurveyBase):
        stations = ["A", "B", "C"]
        lat = [1.0, np.nan, 3.0]
        lon = [4.0, 5.0, 6.0]
        elev = [10.0, 20.0, 30.0]
        distance = [0.0, 5.0, 10.0]
        azimuth = 45.0

    s = Survey()
    summary = s.summary_dict()
    assert summary["n"] == 3
    assert summary["lat_rng"] == (1.0, 3.0)
    assert summary["step"] == 5.0
    assert SurveyBase._as_floats(["bad"]).size == 0
    assert SurveyBase._rng(np.array([])) is None
    assert SurveyBase._fmt_rng(None) == "-"
    assert "+2" in SurveyBase._fmt_list(range(5), maxn=4)
    assert s.format_table([]) == "<empty>"
    table = s.format_table(
        [{"station": "A", "x": 1}, {"station": "B", "x": 20}],
        cols=["station", "x"],
        max_rows=1,
    )
    assert "station" in table and "A" in table and "B" not in table


def test_utils_parsing_formatting_and_stats(capsys) -> None:
    assert _haversine(0, 0, 0, 1) == pytest.approx(111.19, rel=0.01)
    assert _safe_number(None) is None
    assert _safe_number("'12'") == 12
    assert _safe_number("1D2") == 100.0
    assert _safe_number("word") == "word"
    assert _dms_to_deg("10:30:00S") == -10.5
    assert _dms_to_deg("bad") is None
    parsed = parse_kv_pairs('ID=2 LAT=10:30:00S NAME="a b"')
    assert parsed == {"ID": 2, "LAT": -10.5, "NAME": "a b"}
    meas = gather_measurement_key_value_with_str_parser(
        ["ignored", ">HMEAS ID=1 CHTYPE=HX", ">emeas ID=2 CHTYPE=EX"]
    )
    assert [m["KIND"] for m in meas] == ["HMEAS", "EMEAS"]
    assert _format_block_numbers([1, 2, 3], per_line=2).count("\n") == 1
    assert _format_kv("LAT", "1D1") == "LAT=10.0"
    assert _format_kv("NAME", 'a "b"') == 'NAME="a b"'
    assert _quote(None) == '""' and _quote("a b") == '"a b"'
    assert _len(iter([1, 2])) == 1
    quick_edi_stats(total=4, ok=3, label="load", width=60)
    assert "75.00" in capsys.readouterr().out

    stream = io.StringIO()
    show_edi_stats([1, 2], [1], elapsed=0.25, stream=stream, width=40)
    text = stream.getvalue()
    assert "50.00" in text and "0.25" in text


def test_minimum_writer_covers_all_optional_blocks() -> None:
    text = minimum_parser_to_write_edi(
        {
            "head": {
                "DATAID": "S 1",
                "LAT": "1D1",
                "EMPTY": 1e32,
                "NOTE": None,
            },
            "info": ["one", "two"],
            "definemeas": {"MAXCHAN": 2, "UNITS": "M"},
            "measurements": [
                {"KIND": "HMEAS", "ID": 1, "CHTYPE": "HX"},
                {"KIND": "EMEAS", "ID": 2, "CHTYPE": "EX"},
            ],
            "mtsect": {"SECTID": "S 1", "NFREQ": 2},
            "freq": np.array([2.0, 1.0]),
            "zrot": [0.0, 0.0],
            "impedance": {
                key: [1.0, 2.0]
                for key in (
                    "ZXXR",
                    "ZXXI",
                    "ZXX.VAR",
                    "ZXYR",
                    "ZXYI",
                    "ZXY.VAR",
                    "ZYXR",
                    "ZYXI",
                    "ZYX.VAR",
                    "ZYYR",
                    "ZYYI",
                    "ZYY.VAR",
                )
            },
            "resistivity": {
                key: [1.0, 2.0]
                for key in (
                    "RHOROT",
                    "RHOXY",
                    "RHOXY.ERR",
                    "RHOYX",
                    "RHOYX.ERR",
                    "PHSXY",
                    "PHSXY.ERR",
                    "PHSYX",
                    "PHSYX.ERR",
                )
            },
            "coherence": [
                (),
                ("COH MEAS1=1 MEAS2=2", [0.5, 0.6]),
            ],
        }
    )
    assert text.startswith(">HEAD")
    assert ">ZYY.VAR" in text and ">PHSYX.ERR" in text
    assert text.rstrip().endswith(">END")

    multiline = minimum_parser_to_write_edi({"info": "a\nb"})
    assert "  a\n  b" in multiline


def test_strip_item_and_sort_location_edge_cases(tmp_path: Path) -> None:
    assert strip_item(None) is None
    assert strip_item("  x  ") == "x"
    assert strip_item(["//a//", "//b"], item="//") == ["a", "b"]
    arr = strip_item(np.array(["//a//", "//b//"]), item="//")
    assert arr.tolist() == ["a", "b"]
    with pytest.warns(RuntimeWarning):
        assert strip_item([" ", ""]) is None
    with pytest.raises(TypeError):
        strip_item("x", multi_space=0)
    with pytest.raises(TypeError):
        strip_item(123)
    with pytest.raises(FileHandlingError):
        sort_edis_by_location([])
    with pytest.raises(ValueError):
        sort_edis_by_location([object()], metric="bad")
    with pytest.raises(ValueError):
        sort_edis_by_location([object()], by="bad")


def _write_collection_edi(
    path: Path, station: str, lat: str, lon: str, elev: int = 100
) -> Path:
    text = minimum_parser_to_write_edi(
        {
            "head": {
                "DATAID": station,
                "LAT": lat,
                "LONG": lon,
                "ELEV": elev,
                "EMPTY": 1e32,
            },
            "info": ["coverage"],
            "mtsect": {"SECTID": station, "NFREQ": 2},
            "freq": [10.0, 1.0],
            "zrot": [0.0, 0.0],
            "impedance": {
                "ZXXR": [1, 2],
                "ZXXI": [2, 3],
                "ZXYR": [3, 4],
                "ZXYI": [4, 5],
                "ZYXR": [5, 6],
                "ZYXI": [6, 7],
                "ZYYR": [7, 8],
                "ZYYI": [8, 9],
            },
        }
    )
    path.write_text(text, encoding="utf-8")
    return path


def test_collection_get_set_adjust_fetch_export(tmp_path: Path) -> None:
    p1 = _write_collection_edi(tmp_path / "S10.edi", "S10", "10:00:00N", "020:00:00E")
    p2 = _write_collection_edi(
        tmp_path / "S20.edi", "S20", "11:00:00N", "021:00:00E", 200
    )
    col = EDICollection([EDIFile(p1), EDIFile(p2)])

    assert col._resolve("s10") is col[0]
    assert col._resolve("S20.edi") is col[1]
    with pytest.raises(Exception, match="site not found"):
        col._resolve("missing")

    for key in ("freq", "z", "zxx", "zxy", "zyx", "zyy"):
        assert col.get("S10", key) is not None
    assert col.get("S10", "station") == "S10"
    assert col.get("S10", "dataid") == "S10"
    assert col.get("S10", "lat") == pytest.approx(10.0)
    assert col.get("S10", "lon") == pytest.approx(20.0)
    assert col.get("S10", "elev") == pytest.approx(100)
    assert col.get("S10", "filename") == "S10.edi"
    assert str(col.get("S10", "path")).endswith("S10.edi")
    assert col.get("S10", "unknown", "fallback") == "fallback"
    assert col.get("S10", "spectra") is None
    assert col.get("S10", "timeseries") is None
    assert col.get("S10", "tip", "none") == "none"

    marker, ts = object(), object()
    z = np.ones((2, 2, 2), complex)
    tip = np.ones((2, 1, 2), complex)
    changed = col.set(
        "S10",
        update={
            "station": "S10",
            "dataid": "NEW10",
            "lat": 12,
            "lon": 22,
            "elev": 333,
            "freq": [20, 2],
            "z": z,
            "tip": tip,
            "spectra": marker,
            "ts": ts,
        },
    )
    changed_head = changed.get_section("head")
    assert changed_head.dataid == "NEW10"
    assert changed_head.lat == pytest.approx(12)
    assert changed_head.long == pytest.approx(22)
    assert np.array_equal(changed.Z.freq, [20, 2])
    assert col.get("S10", "spectra") is marker
    assert col.get("S10", "ts") is ts
    assert np.array_equal(col.get("S10", "tx"), np.ones(2))
    assert np.array_equal(col.get("S10", "ty"), np.ones(2))

    adjusted = col.adjust(
        "S10", lat=1, lon=2, elev=3, dlat=0.5, dlon=0.25, rename="S10"
    )
    adjusted_head = adjusted.get_section("head")
    assert adjusted_head.lat == pytest.approx(1.5)
    assert adjusted_head.long == pytest.approx(2.25)
    assert adjusted_head.elev == pytest.approx(3)

    assert col.fetch(site="s20", first=True) is col[1]
    assert col.fetch(lat=11, lon=21, first=True) is col[1]
    assert col.fetch(elev=200) == [col[1]]
    assert col.fetch(site="none", first=True) is None
    with pytest.raises(Exception, match="No station"):
        col.fetch(site="none", first=True, strict=True)
    with pytest.raises(Exception, match="No stations"):
        col.fetch(site="none", strict=True)

    result = col.export(tmp_path / "exports", file_pattern="copy_{station}.edi")
    assert len(result["successful"]) == 2
    assert result["failed"] == []
    assert all(Path(p).exists() for p in result["successful"])


def test_collection_merge_validation_and_replacement(tmp_path: Path) -> None:
    one = EDIFile(
        _write_collection_edi(tmp_path / "A1.edi", "A1", "1:00:00N", "001:00:00E")
    )
    two = EDIFile(
        _write_collection_edi(tmp_path / "A2.edi", "A2", "2:00:00N", "002:00:00E")
    )
    col = EDICollection([one])
    with pytest.raises(ValueError):
        col.merge(EDICollection([two]), on_dup="invalid")
    assert col.set("A1", edi=two) is two


def test_survey_views_dataframe_geometry_and_plots(tmp_path: Path) -> None:
    eds = [
        EDIFile(
            _write_collection_edi(
                tmp_path / f"P{i}.edi",
                f"P{i}",
                f"10:0{i}:00N",
                f"020:0{i}:00E",
                100 + i * 10,
            )
        )
        for i in range(3)
    ]
    stations = Stations(eds)
    assert stations.row("P1")["station"] == "P1"
    assert stations.row("none") is None
    assert stations.get("none") is None
    selected = stations.select(
        keys=["P1", "P2"],
        pattern="P*",
        regex=r"P[12]",
        pred=lambda row: row["elev"] > 100,
    )
    assert selected.names() == ["P1", "P2"]
    assert stations.sort(by="elev", reverse=True, inplace=False).names()[0] == "P2"
    assert stations.sort(by="station", inplace=True) is stations
    along, across = stations.offsets(
        origin=(stations._rows[0]["e"], stations._rows[0]["n"]), azimuth=45
    )
    assert along.size == across.size == 3
    stations.set_coords("none", lat=0)
    stations.set_coords("P0", lat=10.5, lon=20.5, elev=500)
    assert eds[0].lat == pytest.approx(10.5)

    frame = stations.to_dataframe(
        columns=["station", "lat", "lon", "elev", "missing"],
        index=None,
        coerce_numeric=False,
    )
    assert list(frame.columns) == ["station", "lat", "lon", "elev"]
    assert Stations([]).to_dataframe(columns=["station"]).empty

    profile = EDIProfile(eds)
    assert profile.get_bearing(method="linear") is not None
    assert profile.get_bearing(method="endpoints") is not None
    assert profile.get_step(method="median") > 0
    assert profile.get_step(as_array=True).size == 2
    px, py = profile.xy
    profile.adjust(origin=(px[0], py[0]), azimuth=30, spacing=25)
    profile.update(use_adjusted=False, update_elev=True)
    assert len(profile.as_table()) == 3
    ax1 = profile.plot_profile(annotate=True, title="raw")
    ax2 = profile.plot_profile(use_adjusted=True, title="adjusted")
    ax3 = profile.plot_track(title="track")
    ax4 = profile.plot_track(use_adjusted=True, title="adjusted track")
    assert all(ax is not None for ax in (ax1, ax2, ax3, ax4))

    topo = Topography(stations)
    topo.smooth(window=3, method="mean")
    topo.detrend().resample(step=max(1, topo.distance[-1] / 5))
    assert topo.gradient(as_degrees=True).size == topo.distance.size - 1
    assert set(topo.to_dict()) == {"distance", "elevation"}
    assert topo.plot(title="topo") is not None

    empty = object.__new__(Topography)
    SurveyBase.__init__(empty)
    empty._d = np.array([], float)
    empty._z = np.array([], float)
    empty._trend = None
    assert empty.smooth().resample(step=1).gradient().size == 0
    assert empty.detrend()._trend is None


def test_validation_helpers_and_lazy_exports(tmp_path: Path, monkeypatch) -> None:
    assert _strip_norm(" 'x' ") == "x"
    assert _to_int_or_none(None) is None
    assert _to_int_or_none("2.0") == 2
    assert _to_int_or_none("bad") is None
    assert _to_float_or_none("") is None
    assert _to_float_or_none("2.5") == 2.5
    assert _to_float_or_none("bad") is None
    assert _split_comment("x // note") == ("x ", "note")
    assert _split_comment("x") == ("x", None)
    assert _norm_str(" '' ") is None and _norm_str(None) is None
    assert _is_tag(" >HEAD x", ">head")
    assert not _is_tag("HEAD", ">HEAD")
    assert _extract_tag(">FREQ //2") == ">FREQ"
    assert _extract_tag("plain") is None
    assert _extract_tag_in(" >= mtsect") == ">=MTSECT"
    assert _extract_tag_in("plain") is None
    lines = [">FREQ //2", " 1.0 2D0", ">END"]
    assert _iter_blocks(lines) == [">FREQ", ">END"]
    assert _count_freq_values(lines, 0) == (2, 2)
    assert _expected_nfreq_from_header("NFREQ=12") == 12
    assert _expected_nfreq_from_header(">FREQ // 3") == 3
    assert _expected_nfreq_from_header(">FREQ") is None
    assert _has_any_data_block([">ZXYR"])
    assert not _has_any_data_block([">HEAD"])

    class Valid(IsEdi):
        @property
        def is_valid(self):
            return True

    class Invalid(IsEdi):
        @property
        def is_valid(self):
            return False

    assert IsEdi._assert_edi(Valid())
    with pytest.raises(Exception, match="invalid state"):
        IsEdi._assert_edi(Invalid())
    with pytest.raises(FileHandlingError):
        IsEdi._assert_edi(None)
    with pytest.raises(FileNotFoundError):
        IsEdi._assert_edi(tmp_path / "missing.edi")
    wrong = tmp_path / "wrong.txt"
    wrong.write_text(">HEAD\n>=OTHERSECT\n>END", encoding="utf-8")
    with pytest.raises(Exception, match="extension"):
        IsEdi._assert_edi(wrong, deep=False)

    import pycsamt.seg as seg

    assert seg.__getattr__("EDIFile") is EDIFile
    with pytest.raises(AttributeError):
        seg.__getattr__("does_not_exist")
    monkeypatch.setattr(seg, "import_module", lambda *a, **k: 1 / 0)
    with pytest.raises(AttributeError, match="unavailable"):
        seg.__getattr__("Spectra")
