# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Shared PCGS (structural geology) app-adapter tests."""

from __future__ import annotations

import base64
import json

from pycsamt.app._structure import (
    FAULT_COLUMNS,
    LINEAR_COLUMNS,
    PLANAR_COLUMNS,
    decode_structure_csv_upload,
    decode_structure_json_upload,
    model_from_table_rows,
    store_from_structure,
    structure_from_store,
    structure_scene_traces,
    structure_section_figure,
    table_rows_from_model,
)
from pycsamt.format.structure import StructModel, structure_to_dict
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
                z_top=50.0, line="L1",
            ),
        ],
    )


def _json_upload(model: StructuralModel) -> str:
    doc = StructModel.from_structural_model(model, document_id="pcgs:t")
    raw = json.dumps(structure_to_dict(doc)).encode()
    return "data:application/json;base64," + base64.b64encode(raw).decode()


def _csv_upload(text: str) -> str:
    return "data:text/csv;base64," + base64.b64encode(text.encode()).decode()


def test_decode_json_upload_round_trips():
    store = decode_structure_json_upload(_json_upload(_model()), "s.pcgs.json")
    assert store["n_planar"] == 1 and store["n_linear"] == 1
    assert store["n_faults"] == 1
    restored = structure_from_store(store)
    assert restored.faults[0].line == "L1"


def test_decode_planar_csv_upload():
    csv = (
        "x,kind,strike_deg,dip_deg,dip_direction_deg,line\n"
        "100,bedding,45,30,135,L1\n"
    )
    rows = decode_structure_csv_upload(
        _csv_upload(csv), "planar.csv", kind="planar"
    )
    assert len(rows) == 1
    assert set(rows[0]) == set(PLANAR_COLUMNS)
    assert rows[0]["kind"] == "bedding"


def test_decode_faults_csv_upload():
    csv = "x,dip_deg,downthrown_side,line\n500,70,right,L1\n"
    rows = decode_structure_csv_upload(
        _csv_upload(csv), "faults.csv", kind="faults"
    )
    assert len(rows) == 1 and set(rows[0]) == set(FAULT_COLUMNS)


def test_decode_linear_csv_upload():
    csv = "x,kind,trend_deg,plunge_deg,line\n120,fold_axis,210,15,L1\n"
    rows = decode_structure_csv_upload(
        _csv_upload(csv), "linear.csv", kind="linear"
    )
    assert len(rows) == 1 and set(rows[0]) == set(LINEAR_COLUMNS)


def test_table_round_trip_preserves_line():
    model = _model()
    rows = table_rows_from_model(model)
    restored = model_from_table_rows(
        rows["planar"], rows["linear"], rows["faults"]
    )
    assert restored.planar[0].line == "L1"
    assert restored.faults[0].line == "L1"


def test_model_from_table_rows_skips_incomplete_rows():
    planar_rows = [
        {"x": 100, "kind": "bedding", "strike_deg": 45, "dip_deg": 30,
         "dip_direction_deg": 135},
        {"x": None, "kind": "", "strike_deg": None, "dip_deg": None,
         "dip_direction_deg": None},
    ]
    model = model_from_table_rows(planar_rows, [], [])
    assert len(model.planar) == 1


def test_section_figure_empty_model_shows_placeholder():
    fig = structure_section_figure(None)
    assert fig.layout.annotations
    assert not fig.data


def test_section_figure_draws_one_trace_per_item():
    fig = structure_section_figure(_model())
    assert len(fig.data) == 3


def test_scene_traces_skip_items_on_unknown_line():
    model = _model()
    traces = structure_scene_traces(
        model, line_offsets={"L2": 5.0}, default_line=None,
    )
    assert traces == []


def test_scene_traces_place_items_using_line_offset():
    model = _model()
    traces = structure_scene_traces(
        model, line_offsets={"L1": 10.0}, default_line="L1",
    )
    kinds = [type(t).__name__ for t in traces]
    assert kinds.count("Mesh3d") == 1
    assert kinds.count("Scatter3d") == 2


def test_scene_traces_use_surface_for_elevation_datum():
    model = _model()
    traces = structure_scene_traces(
        model, line_offsets={"L1": 0.0}, default_line="L1",
        surface=lambda u: 200.0,
    )
    mesh = next(t for t in traces if type(t).__name__ == "Mesh3d")
    # z_top = surface(u) - fault.z_top = 200 - 50 = 150
    assert max(mesh.z) == 150.0


def test_store_from_structure_round_trips_document():
    doc = StructModel.from_structural_model(_model(), document_id="pcgs:x")
    store = store_from_structure(doc)
    restored = structure_from_store(store)
    assert len(restored.planar) == 1
