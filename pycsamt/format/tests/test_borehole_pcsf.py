"""PCBH embedding and PCSF block-alignment tests."""

from __future__ import annotations

import h5py
import numpy as np
import pytest

from pycsamt.app._borehole import add_embedded_pcbh_to_pcsf_figure
from pycsamt.format import (
    DerivedVolume,
    Grid2DGeometry,
    Grid3DGeometry,
    LineEntry,
    MultilineGeometry,
    PCSFModel,
    read_pcsf,
    write_pcsf,
)
from pycsamt.format.borehole import (
    Collar,
    CoordinateReferenceSystem,
    PCBHBorehole,
    PCBHDocument,
    align_pcbh_to_pcsf,
    embed_pcbh,
    extract_pcbh,
    pcbh_document_checksum,
    reference_pcbh,
)


def _document(x=500050.0, y=4500050.0) -> PCBHDocument:
    return PCBHDocument(
        document_id="test:pcsf-association",
        created_at="2026-08-28T00:00:00Z",
        created_by="pytest",
        crs=CoordinateReferenceSystem(
            "EPSG:32629",
            vertical="EPSG:5714",
        ),
        boreholes=[
            PCBHBorehole(
                id="BH-1",
                name="Borehole 1",
                kind="mining_exploration",
                status="completed",
                collar=Collar(x, y, 0.0),
                total_depth_md=80.0,
            )
        ],
    )


def _model() -> PCSFModel:
    geometry = Grid3DGeometry(
        x=np.array([25.0, 75.0]),
        y=np.array([25.0, 75.0]),
        z=np.array([20.0, 60.0, 100.0]),
        x_nodes=np.array([0.0, 50.0, 100.0]),
        y_nodes=np.array([0.0, 50.0, 100.0]),
        z_nodes=np.array([0.0, 40.0, 80.0, 120.0]),
        origin=np.array([500000.0, 4500000.0, 0.0]),
    )
    return PCSFModel(
        geometry=geometry,
        resistivity=np.full(geometry.resistivity_shape, 100.0),
        crs="EPSG:32629",
    )


def test_embedded_pcbh_round_trips_through_optional_hdf5_group(tmp_path):
    document = _document()
    model = embed_pcbh(_model(), document, uri="project/holes.pcbh.json")
    path = write_pcsf(model, tmp_path / "portable.pcsf")
    restored = read_pcsf(path)

    assert extract_pcbh(restored) == document
    assert restored.boreholes.reference.sha256 == pcbh_document_checksum(
        document
    )
    with h5py.File(path, "r") as handle:
        assert handle["boreholes"].attrs["kind"] == "pcbh"


def test_reference_only_round_trip_and_checksum_validation(tmp_path):
    checksum = pcbh_document_checksum(_document())
    model = reference_pcbh(_model(), "https://example.org/holes.json", checksum)
    restored = read_pcsf(write_pcsf(model, tmp_path / "reference.pcsf"))

    assert extract_pcbh(restored) is None
    assert restored.boreholes.reference.uri.startswith("https://")
    with pytest.raises(ValueError, match="64 hexadecimal"):
        reference_pcbh(_model(), "holes.json", "bad")


def test_alignment_requires_vertical_decision_and_reports_inside():
    with pytest.raises(ValueError, match="vertical datums"):
        align_pcbh_to_pcsf(_document(), _model())

    report = align_pcbh_to_pcsf(
        _document(),
        _model(),
        vertical_offset=0.0,
    )
    path = report.trajectories["BH-1"]
    assert (path[0].x, path[0].y, path[0].z) == (50.0, 50.0, 0.0)
    assert path[-1].z == pytest.approx(80.0)
    assert report.bounds_results[0].relation == "inside"
    assert not report.vertical_compatible


def test_alignment_reports_outside_and_supports_rotation():
    model = _model()
    model.geometry.rotation_deg = 90.0
    report = align_pcbh_to_pcsf(
        _document(x=500200.0, y=4500050.0),
        model,
        vertical_offset=0.0,
    )
    assert report.bounds_results[0].relation == "outside"
    assert report.trajectories["BH-1"][0].x == pytest.approx(50.0)
    assert report.trajectories["BH-1"][0].y == pytest.approx(-200.0)


def test_legacy_pcsf_without_borehole_group_still_reads(tmp_path):
    path = write_pcsf(_model(), tmp_path / "legacy.pcsf")
    restored = read_pcsf(path)
    assert restored.boreholes is None


def test_embedded_borehole_renders_in_native_block_coordinates():
    import plotly.graph_objects as go

    model = embed_pcbh(_model(), _document())
    figure = add_embedded_pcbh_to_pcsf_figure(
        go.Figure(),
        model,
        vertical_offset=0.0,
    )
    trace = figure.data[0]
    assert trace.x[0] == 50.0
    assert trace.y[0] == 50.0
    assert trace.z[-1] == pytest.approx(-80.0)
    assert trace.customdata[0][2] == "inside"


def test_synthesized_multiline_block_uses_derived_volume_frame():
    grid = _model().geometry
    line_geometry = Grid2DGeometry(
        x=np.array([0.0, 100.0]),
        z=np.array([20.0, 60.0]),
    )
    multiline = MultilineGeometry(
        lines=[
            LineEntry(
                "L1",
                line_geometry,
                np.full((2, 2), 100.0),
            )
        ],
        derived_volume=DerivedVolume(
            grid=grid,
            resistivity=np.full(grid.resistivity_shape, 100.0),
            derived_from=["L1"],
        ),
    )
    model = PCSFModel(
        geometry=multiline,
        crs="EPSG:32629",
        boreholes=embed_pcbh(_model(), _document()).boreholes,
    )
    report = align_pcbh_to_pcsf(
        _document(),
        model,
        vertical_offset=0.0,
    )
    assert report.bounds_results[0].relation == "inside"
