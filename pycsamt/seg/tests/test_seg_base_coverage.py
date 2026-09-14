from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from pycsamt.seg.base import Base, EDIComponentBase, EdiFileBase, SurveyBase


# ─────────────────────────────────────────────────────────────────────────
# EDIComponentBase: __repr__ / __str__ without a section
# ─────────────────────────────────────────────────────────────────────────


class _NoSectionComponent(EDIComponentBase):
    """Non-dataclass component with plain class annotations."""

    name: str | None = None
    count: int = 0

    def __init__(self, **kw):
        super().__init__()
        for k, v in kw.items():
            setattr(self, k, v)


def test_repr_lists_key_value_pairs():
    c = _NoSectionComponent(name="S1", count=3)
    text = repr(c)
    assert text.startswith("_NoSectionComponent(")
    assert "name=S1" in text
    assert "count=3" in text


def test_repr_truncates_when_too_many_items():
    class ManyFields(EDIComponentBase):
        def __init__(self, **kw):
            super().__init__()
            self.__dict__.update(kw)

    c = ManyFields(**{f"f{i}": i for i in range(10)})
    text = repr(c)
    assert "…" in text


def test_repr_truncates_when_too_many_chars():
    c = _NoSectionComponent(name="x" * 200, count=1)
    text = repr(c)
    assert "…" in text


def test_str_without_section_uses_class_name_prefix():
    c = _NoSectionComponent(name="S1", count=2)
    assert c._section is None
    text = str(c)
    assert text.startswith("_NoSectionComponent ")
    assert "name=S1" in text


# ─────────────────────────────────────────────────────────────────────────
# _field_names: non-dataclass branches
# ─────────────────────────────────────────────────────────────────────────


def test_field_names_uses_annotations_when_not_dataclass():
    names = _NoSectionComponent._field_names()
    assert set(names) == {"name", "count"}


class _NoAnnotationsComponent(EDIComponentBase):
    # Shadow the inherited (non-empty) EDIComponentBase.__annotations__
    # so ``getattr(cls, "__annotations__", None)`` truly resolves empty.
    __annotations__ = {}


def test_field_names_empty_when_no_dataclass_and_no_annotations():
    assert _NoAnnotationsComponent._field_names() == []


# ─────────────────────────────────────────────────────────────────────────
# _iter_kv: annotated-but-unset attribute is skipped
# ─────────────────────────────────────────────────────────────────────────


class _PartiallySetComponent(EDIComponentBase):
    name: str | None = None
    missing_attr: str  # declared, never assigned, and not a class default

    def __init__(self, **kw):
        super().__init__()
        self.name = kw.get("name")


def test_iter_kv_skips_annotated_attribute_that_was_never_set():
    c = _PartiallySetComponent(name="S1")
    keys = dict(c._iter_kv(True))
    assert "name" in keys
    assert "missing_attr" not in keys


# ─────────────────────────────────────────────────────────────────────────
# _format_value / _repr_value: remaining branches
# ─────────────────────────────────────────────────────────────────────────


def test_format_value_plain_string_no_quotes_needed():
    c = _NoSectionComponent()
    assert c._format_value("plainvalue") == "plainvalue"


def test_repr_value_string_short_and_long():
    assert EDIComponentBase._repr_value("short") == "short"
    long_text = "x" * 30
    assert EDIComponentBase._repr_value(long_text) == f'"{long_text}"'
    assert EDIComponentBase._repr_value("has space") == '"has space"'


def test_repr_value_float():
    assert EDIComponentBase._repr_value(3.14159265) == "3.14159"


def test_repr_value_fallback_uses_repr():
    assert EDIComponentBase._repr_value(42) == repr(42)
    assert EDIComponentBase._repr_value(None) == repr(None)


# ─────────────────────────────────────────────────────────────────────────
# EdiFileBase: from_file / write_file list-compose / strict validation /
# location proxies fallback / format_kv / format_data_block
# ─────────────────────────────────────────────────────────────────────────


class _ListComposeEdi(EdiFileBase):
    def __init__(self, lines=None, **kw):
        self._lines = lines or ["a\n", "b\n"]
        super().__init__(**kw)

    def read(self):
        return self

    def compose(self):
        return self._lines


def test_write_file_joins_list_compose_output(tmp_path):
    ed = _ListComposeEdi(strict_validate=False)
    out = ed.write_file(tmp_path / "out.edi")
    assert out.read_text(encoding="utf-8") == "a\nb\n"


def test_from_file_reads_and_returns_instance(tmp_path):
    target = tmp_path / "mem.edi"
    target.write_text("dummy", encoding="utf-8")
    ed = _ListComposeEdi.from_file(target, strict_validate=False)
    assert isinstance(ed, _ListComposeEdi)
    assert ed.path == target


def test_strict_validate_true_calls_validate_path_on_real_edi(edi_imp_file):
    ed = _ListComposeEdi(path=None)
    # Constructing with a genuinely valid EDI path and strict_validate=True
    # exercises EdiFileBase._validate_path via IsEdi._assert_edi.
    ed2 = _ListComposeEdi(strict_validate=True, path=edi_imp_file)
    assert ed2.path == Path(edi_imp_file)


def test_location_proxies_fall_back_to_definemeasurement():
    class Head(EDIComponentBase):
        def __init__(self):
            super().__init__()

    class DM(EDIComponentBase):
        def __init__(self):
            super().__init__()
            self.reflat = 11.0
            self.reflong = 22.0
            self.refelev = 33.0

    ed = _ListComposeEdi(strict_validate=False)
    ed.add_section("head", Head())
    ed.add_section("definemeasurement", DM())
    assert ed.lat == 11.0
    assert ed.lon == 22.0
    assert ed.elev == 33.0


def test_format_kv_none_value_and_unquoted_string():
    ed = _ListComposeEdi(strict_validate=False)
    assert ed.format_kv("id", None) == "  ID=\n"
    assert ed.format_kv("name", "a b", quote=False) == "  NAME=a b\n"


def test_format_data_block_no_count_comment():
    ed = _ListComposeEdi(strict_validate=False)
    block = ed.format_data_block(
        "freq", [1, 2, 3], count_comment=False, header_comment=False,
    )
    assert block[0].strip() == ">FREQ"
    assert "//" not in block[0]


# ─────────────────────────────────────────────────────────────────────────
# SurveyBase: fallback branches with a minimally-attributed subclass
# ─────────────────────────────────────────────────────────────────────────


def test_fmt_rng_with_actual_range():
    assert SurveyBase._fmt_rng((1.0, 2.5)) == "1.000..2.500"


def test_fmt_list_short_list_no_truncation():
    assert SurveyBase._fmt_list(["a", "b"], maxn=4) == "a, b"
    assert SurveyBase._fmt_list([], maxn=4) == "-"


class _EmptySurvey(SurveyBase):
    """No lat/lon/elev/distance/stations/azimuth attributes at all."""


def test_stations_lat_lon_elev_distance_default_to_empty():
    s = _EmptySurvey()
    assert s._stations() == []
    assert s._lat().size == 0
    assert s._lon().size == 0
    assert s._elev().size == 0
    assert s._distance().size == 0


class _StationsAsTableSurvey(SurveyBase):
    def as_table(self):
        return [{"station": "A"}, {"station": "B"}]


def test_stations_falls_back_to_as_table():
    s = _StationsAsTableSurvey()
    assert s._stations() == ["A", "B"]


class _StationsCallableSurvey(SurveyBase):
    def stations(self):
        return ["X", "Y"]


def test_stations_callable_property():
    s = _StationsCallableSurvey()
    assert s._stations() == ["X", "Y"]


class _StationsRaisingSurvey(SurveyBase):
    # A callable ``stations`` exists fine (hasattr never executes it);
    # it only raises once actually *called*, which is what the
    # try/except around ``ss()`` in ``_stations`` is meant to catch.
    def stations(self):
        raise RuntimeError("boom")

    def as_table(self):
        return [{"station": "Z"}]


def test_stations_falls_back_to_as_table_when_stations_raises():
    s = _StationsRaisingSurvey()
    assert s._stations() == ["Z"]


class _AsTableRaisingSurvey(SurveyBase):
    def as_table(self):
        raise RuntimeError("boom")


def test_stations_returns_empty_when_as_table_raises_and_no_stations():
    assert _AsTableRaisingSurvey()._stations() == []


class _ElevationAliasSurvey(SurveyBase):
    elevation = [1.0, 2.0]


def test_elev_falls_back_to_elevation_alias():
    s = _ElevationAliasSurvey()
    assert s._elev().tolist() == [1.0, 2.0]


class _AzimuthRaisesSurvey(SurveyBase):
    azimuth = "not-a-number"


def test_azimuth_swallows_conversion_error():
    s = _AzimuthRaisesSurvey()
    assert s._azimuth() is None


class _BearingSurvey(SurveyBase):
    def get_bearing(self):
        return 12.5


def test_azimuth_falls_back_to_get_bearing():
    s = _BearingSurvey()
    assert s._azimuth() == 12.5


class _BearingNoneSurvey(SurveyBase):
    def get_bearing(self):
        return None


def test_azimuth_get_bearing_returns_none():
    s = _BearingNoneSurvey()
    assert s._azimuth() is None


class _BearingRaisesSurvey(SurveyBase):
    def get_bearing(self):
        raise RuntimeError("boom")


def test_azimuth_get_bearing_raises_is_swallowed():
    s = _BearingRaisesSurvey()
    assert s._azimuth() is None


def test_step_returns_none_when_not_enough_distance_points():
    s = _EmptySurvey()
    assert s._step() is None


class _GetStepScalarSurvey(SurveyBase):
    def get_step(self):
        return 5.0


def test_step_uses_get_step_scalar():
    s = _GetStepScalarSurvey()
    assert s._step() == 5.0


class _GetStepArraySurvey(SurveyBase):
    def get_step(self):
        return [1.0, 3.0, 6.0]


def test_step_uses_get_step_array_median_diff():
    s = _GetStepArraySurvey()
    assert s._step() == pytest.approx(2.5)


class _GetStepRaisesSurvey(SurveyBase):
    distance = [0.0, 4.0, 10.0]

    def get_step(self):
        raise RuntimeError("boom")


def test_step_falls_back_to_distance_when_get_step_raises():
    s = _GetStepRaisesSurvey()
    assert s._step() == pytest.approx(5.0)


def test_n_falls_back_through_lat_stations_len():
    class LatOnly(SurveyBase):
        lat = [1.0, 2.0, 3.0]

    assert LatOnly()._n() == 3

    class StationsOnly(SurveyBase):
        stations = ["A", "B"]

    assert StationsOnly()._n() == 2

    class LenOnly(SurveyBase):
        def __len__(self):
            return 7

    assert LenOnly()._n() == 7

    assert _EmptySurvey()._n() == 0


# ─────────────────────────────────────────────────────────────────────────
# format_table: default rows-from-as_table / default cols
# ─────────────────────────────────────────────────────────────────────────


class _TableSurvey(SurveyBase):
    def as_table(self):
        return [{"station": "A", "x": 1}, {"station": "B", "x": 22}]


def test_format_table_uses_as_table_and_default_columns():
    s = _TableSurvey()
    text = s.format_table()
    assert "station" in text
    assert "A" in text and "B" in text


def test_format_table_as_table_raises_returns_empty_marker():
    class Raising(SurveyBase):
        def as_table(self):
            raise RuntimeError("boom")

    assert Raising().format_table() == "<empty>"


def test_format_table_no_as_table_and_no_rows_returns_empty_marker():
    assert _EmptySurvey().format_table() == "<empty>"
