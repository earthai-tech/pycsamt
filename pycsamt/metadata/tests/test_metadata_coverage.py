from __future__ import annotations

from datetime import date
from types import SimpleNamespace

import numpy as np
import pytest

from pycsamt.metadata._property import (
    CopyrightInfo,
    MetadataError,
    Person,
    Reference,
    Software,
)
from pycsamt.metadata.geology import Formation, GeologyCatalog, geology_prior
from pycsamt.metadata.instrument import (
    InstrumentMeta,
    SensorSpec,
    known_system,
    list_presets,
)
from pycsamt.metadata.rocks import RockProperties
from pycsamt.metadata.survey import BBox, SurveyMeta


def test_reference_validation_serialization_bibtex_and_repr():
    with pytest.raises(MetadataError, match="Invalid year"):
        Reference("A", "T", year=-1)
    with pytest.raises(MetadataError, match="Invalid DOI"):
        Reference("A", "T", doi="invalid")

    ref = Reference.from_dict(
        {
            "author": "Doe, J.",
            "title": "Study",
            "journal": "Journal",
            "year": 2024,
            "volume": "2",
            "pages": "1--3",
            "doi": "10.1234/example",
            "publisher": "EarthAI",
        }
    )
    assert ref.to_dict()["publisher"] == "EarthAI"
    bib = ref.to_bibtex("Doe2024")
    assert "@article{Doe2024" in bib and "publisher = {EarthAI}" in bib
    assert "Study" in repr(ref)


def test_property_value_objects_flatten_extra_fields():
    ref = Reference("A", "T")
    copyright_info = CopyrightInfo("open", "cite", ref, {"license": "CC"})
    person = Person("Analyst", "a@example.org", extra={"orcid": "1"})
    software = Software("Tool", "1.0", author=person, extra={"build": 3})
    assert copyright_info.to_dict()["license"] == "CC"
    assert copyright_info.to_dict()["reference"]["title"] == "T"
    assert person.to_dict()["orcid"] == "1"
    assert software.to_dict()["author"]["name"] == "Analyst"
    assert software.to_dict()["build"] == 3
    assert "open" in repr(copyright_info)
    assert "Analyst" in repr(person)
    assert "Tool" in repr(software)


def test_rock_properties_queries_and_defensive_copies():
    rocks = RockProperties()
    ranges = rocks.resistivity_ranges
    patterns = rocks.hatch_patterns
    ranges.pop("shale")
    patterns.pop("shale")
    assert rocks.get_resistivity("shale") == [50.12, 32.0]
    assert rocks.get_pattern("shale")[0] == "="
    assert "metamorphic rock" in rocks.find_matching_rocks("MORPH")
    with pytest.raises(KeyError):
        rocks.get_resistivity("unobtainium")


def _formation(name="custom"):
    return Formation(
        name,
        resistivity_range=(1.0, 100.0),
        depth_range=(10.0, 1000.0),
        n_layers_range=(2, 6),
        description="test",
        rock_types=["shale"],
        extra={"region": "demo"},
    )


def test_formation_validation_computed_fields_and_roundtrips():
    with pytest.raises(ValueError, match="> 0"):
        Formation("bad", (0, 1), (0, 1), (1, 2))
    swapped = Formation("swap", (100, 1), (20, 10), (6, 2))
    assert swapped.resistivity_range == (1, 100)
    assert swapped.depth_range == (10, 20)
    assert swapped.n_layers_range == (2, 6)
    assert swapped.log_rho_range == (0.0, 2.0)
    assert swapped.rho_mid == 10.0
    assert swapped.depth_mid == 15.0
    assert swapped.n_layers_mid == 4

    formation = _formation()
    restored = Formation.from_dict(formation.to_dict())
    assert restored.extra["region"] == "demo"
    prior = formation.to_prior()
    assert Formation.from_prior("prior", prior).rho_mid == pytest.approx(10.0)
    assert "custom" in repr(formation)


def test_geology_catalog_mutation_queries_exports_and_compatibility():
    broad = _formation("broad")
    narrow = Formation("narrow", (5, 20), (100, 200), (2, 3), rock_types=["clay"])
    catalog = GeologyCatalog([broad, narrow])
    assert len(catalog) == 2 and "BROAD" in catalog
    assert catalog.names() == ["broad", "narrow"]
    assert list(catalog)[0] is broad
    assert catalog.get("BROAD") is broad
    with pytest.raises(KeyError, match="not found"):
        catalog.get("missing")
    assert catalog.lookup_by_resistivity(10, n=1) == [narrow]
    assert catalog.lookup_by_resistivity(10000, n=1)
    assert catalog.lookup_by_depth(150, n=1) == [narrow]
    assert catalog.lookup_by_depth(5000, n=1)
    assert catalog.lookup_by_rock_type("SHA") == [broad]
    assert catalog.all_scenarios() == {"broad": broad, "narrow": narrow}
    assert catalog.to_prior("broad")["n_layers"] == (2, 6)
    assert list(catalog.to_dataframe(api=False)["name"]) == ["broad", "narrow"]
    assert "2 formations" in repr(catalog)
    catalog.remove("narrow")
    catalog.add(narrow)
    catalog.reset()
    assert len(catalog) > 2
    assert geology_prior("sedimentary")["n_layers"]


def test_sensor_and_instrument_roundtrip_and_head_bridge(tmp_path):
    mag = SensorSpec("induction_coil", "MTC", (0.1, 100.0), "calibrated")
    elec = SensorSpec("electrode", "PbCl", None)
    assert mag.covers(1.0) and not mag.covers(1000.0)
    assert elec.covers(1e9)
    assert "MTC" in str(mag) and "full band" in str(elec)
    assert SensorSpec._from_dict(mag._to_dict()).frequency_range == (0.1, 100.0)

    inst = InstrumentMeta(
        system="System", serial="S1", magnetic_sensor=mag,
        electric_sensor=elec, software_version="2", notes="field",
    )
    assert inst.label == "System / S1"
    assert "Mag. sensor" in inst.summary() and "Notes" in inst.summary()
    assert "System" in repr(inst)
    assert inst.to_head_fields() == {"acqby": "System / S1", "progvers": "2"}
    from_head = InstrumentMeta.from_head(SimpleNamespace(acqby="Sys / 9", progvers="v"))
    assert from_head.system == "Sys" and from_head.serial == "9"
    assert InstrumentMeta.from_head(SimpleNamespace()).system == ""

    restored = InstrumentMeta.from_json(inst.to_json())
    assert restored.magnetic_sensor.model == "MTC"
    assert InstrumentMeta.from_yaml(inst.to_yaml()).electric_sensor.model == "PbCl"
    for suffix, fmt in ((".json", "json"), (".yaml", "yaml")):
        path = tmp_path / f"instrument{suffix}"
        inst.save(path, fmt=fmt)
        assert InstrumentMeta.load(path).system == "System"
    with pytest.raises(ValueError, match="Unknown format"):
        inst.save(tmp_path / "bad.txt", fmt="txt")


def test_instrument_presets_and_defaults():
    assert list_presets() == sorted(list_presets())
    assert known_system("Phoenix-V8").system == "Phoenix V8"
    with pytest.raises(KeyError, match="Unknown preset"):
        InstrumentMeta.from_preset("missing")
    fields = InstrumentMeta().to_head_fields()
    assert fields["progvers"].startswith("pyCSAMT")


def test_bbox_extra_membership_and_survey_serialization(tmp_path):
    bbox = BBox(0, 1, 2, 3)
    assert 0.5 in bbox
    assert (0.5, 2.5) in bbox
    assert [0.5] not in bbox
    assert "BBox" in repr(bbox)

    survey = SurveyMeta(
        name="Demo", project="P", operator="O", method="MT",
        bbox=bbox, n_stations=2, date_start=date(2024, 1, 1),
        date_end=date(2024, 1, 3), extra={"x": 1},
    )
    assert survey.to_dict()["bbox"] == bbox.to_dict()
    assert SurveyMeta.from_dict(survey.to_dict()).duration_days == 2
    yaml_path = tmp_path / "survey.yaml"
    survey.to_yaml(yaml_path)
    assert SurveyMeta.from_yaml(yaml_path).name == "Demo"
    assert "Survey" in survey.summary()
    assert "Demo" in repr(survey)


def test_survey_edge_bridges_invalid_sites_and_optional_fields(tmp_path):
    with pytest.raises(ValueError, match="lon_min"):
        BBox(0, 1, 3, 2)

    sites = [object(), SimpleNamespace(coords=("bad", 2, 0))]
    survey = SurveyMeta.from_sites(sites, name="edge")
    assert survey.n_stations == 2 and survey.bbox is None

    class Head:
        project = None
        acqby = None

    head = Head()
    described = SurveyMeta(
        name="edge",
        project="Project",
        operator="Operator",
        notes="Notes",
        date_start=date(2024, 1, 1),
        n_stations=0,
    )
    described.update_edi_head(head)
    assert head.project == "Project" and head.acqby == "Operator"
    summary = described.summary()
    assert "Notes" in summary and "?" in summary

    restored = SurveyMeta.from_dict(
        {"name": "dates", "date_start": date(2024, 1, 1), "date_end": None}
    )
    assert restored.date_start == date(2024, 1, 1)
    assert restored.date_end is None
    json_path = tmp_path / "minimal.json"
    described.to_json(json_path)
    assert SurveyMeta.from_json(json_path).operator == "Operator"
