"""Tests for the PCGS structural-geology mini-format."""

from __future__ import annotations

import pytest

from pycsamt.format.structure import (
    StructModel,
    StructModelValidationError,
    read_structure,
    structure_from_dict,
    structure_to_dict,
    write_structure,
)
from pycsamt.geology.structural import (
    FaultTrace,
    LinearMeasurement,
    StructuralMeasurement,
    StructuralModel,
)


def _model() -> StructuralModel:
    return StructuralModel(
        planar=[
            StructuralMeasurement(
                x=100.0, kind="bedding", strike_deg=45.0, dip_deg=30.0,
                dip_direction_deg=135.0, line="L1",
            ),
        ],
        linear=[
            LinearMeasurement(
                x=120.0, kind="fold_axis", trend_deg=210.0, plunge_deg=15.0,
                line="L1",
            ),
        ],
        faults=[
            FaultTrace(
                x=500.0, dip_deg=70.0, downthrown_side="right",
                throw_m=12.0, line="L1",
            ),
        ],
    )


def _doc() -> StructModel:
    return StructModel.from_structural_model(
        _model(), document_id="pcgs:test", title="Test structure"
    )


def test_json_round_trip(tmp_path):
    original = _doc()
    path = tmp_path / "structure.pcgs.json"
    write_structure(original, path)
    restored = read_structure(path)
    assert len(restored) == 3
    assert restored.faults[0].line == "L1"
    assert restored.planar[0].dip_deg == 30.0
    assert restored.linear[0].trend_deg == 210.0
    assert restored.pcgs_version == "0.1.0"
    assert restored.title == "Test structure"


def test_from_structural_model_default_is_valid():
    doc = StructModel.from_structural_model()
    assert len(doc) == 0
    assert doc.issues() == []


def test_csv_round_trip(tmp_path):
    planar_csv = tmp_path / "planar.csv"
    planar_csv.write_text(
        "x,kind,strike_deg,dip_deg,dip_direction_deg,line\n"
        "100,bedding,45,30,135,L1\n",
        encoding="utf-8",
    )
    faults_csv = tmp_path / "faults.csv"
    faults_csv.write_text(
        "x,dip_deg,downthrown_side,line\n500,70,right,L1\n",
        encoding="utf-8",
    )
    doc = StructModel.from_csv(
        planar_path=planar_csv, faults_path=faults_csv, title="from csv"
    )
    assert len(doc.planar) == 1 and len(doc.faults) == 1
    assert doc.faults[0].line == "L1"


def test_validation_rejects_wrong_version():
    doc = _doc()
    doc.pcgs_version = "9.0.0"
    with pytest.raises(StructModelValidationError):
        doc.validate()


def test_structure_from_dict_rejects_non_dict():
    with pytest.raises(StructModelValidationError):
        structure_from_dict([])  # type: ignore[arg-type]


def test_structure_to_dict_is_json_safe_and_stable_keys():
    payload = structure_to_dict(_doc())
    assert payload["pcgs_version"] == "0.1.0"
    assert payload["planar"][0]["kind"] == "bedding"
    assert payload["faults"][0]["downthrown_side"] == "right"
