# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Round-trip tests for pycsamt.format.text — write_pcsm / read_pcsm.

Mirrors test_io.py's fixtures and coverage: PCSM must round-trip the
same four geometry kinds, both topography kinds, survey/metadata/
history, and enforce the same PCSM_VERSION policy as .pcsf.
"""

from __future__ import annotations

import numpy as np
import pytest

from pycsamt.format.io import read_pcsf, write_pcsf
from pycsamt.format.schema import (
    DerivedVolume,
    Grid2DGeometry,
    Grid3DGeometry,
    LineEntry,
    MultilineGeometry,
    PCSF_VERSION,
    PCSFModel,
    StationTable,
    TopographyPerStation,
    TopographyRaster,
    UnstructuredMeshGeometry,
)
from pycsamt.format.text import (
    _DEFAULT_ROW_WIDTH,
    pcsf_to_pcsm,
    pcsm_to_pcsf,
    read_pcsm,
    write_pcsm,
)


def _grid2d_model() -> PCSFModel:
    x = np.array([0.0, 100.0, 200.0, 300.0])
    z = np.array([10.0, 50.0, 150.0])
    rho_log10 = np.array(
        [
            [2.0, 2.1, 2.2, 2.3],
            [2.5, 2.6, 2.7, 2.8],
            [3.0, 3.1, 3.2, 3.3],
        ]
    )
    geometry = Grid2DGeometry(
        x=x, z=z, x_nodes=np.array([-50.0, 50.0, 150.0, 250.0, 350.0]),
        z_nodes=np.array([0.0, 30.0, 100.0, 200.0]),
        origin=np.array([500000.0, 4500000.0]), azimuth_deg=45.0,
    )
    return PCSFModel(
        geometry=geometry,
        resistivity=10.0**rho_log10,
        resistivity_native=rho_log10,
        resistivity_native_encoding="log10",
        uncertainty=np.full_like(rho_log10, 0.1),
        stations=StationTable(
            name=["S0", "S1", "S2", "S3"], x=x, y=np.zeros(4),
            z=np.array([10.0, 12.0, 9.0, 11.0]),
            lon=np.array([7.50, 7.51, 7.52, 7.53]),
            lat=np.array([45.10, 45.11, 45.12, 45.13]),
        ),
        topography=TopographyPerStation(
            station_id=["S0", "S1", "S2", "S3"],
            elevation=np.array([10.0, 12.0, 9.0, 11.0]),
        ),
        survey={"name": "demo survey", "n_stations": 4, "bbox": [0.0, 0.0, 300.0, 0.0]},
        history={"rms": np.array([3.2, 2.1, 1.4, 1.05]), "lambda": np.array([10.0, 5.0, 2.0, 1.0])},
        source_backend="occam2d",
        created_by="pytest",
        crs="EPSG:32633",
        description="synthetic round-trip fixture",
        metadata={"note": "unit test"},
    )


def _grid3d_model() -> PCSFModel:
    geometry = Grid3DGeometry(
        x=np.array([0.0, 100.0, 200.0]),
        y=np.array([0.0, 100.0]),
        z=np.array([10.0, 50.0, 150.0, 400.0]),
        origin=np.array([500000.0, 4500000.0, 0.0]),
        rotation_deg=12.5,
        n_air=3,
    )
    rho = np.full(geometry.resistivity_shape, 100.0) + np.arange(
        np.prod(geometry.resistivity_shape)
    ).reshape(geometry.resistivity_shape)
    return PCSFModel(geometry=geometry, resistivity=rho, source_backend="modem3d")


def _mesh_model() -> PCSFModel:
    geometry = UnstructuredMeshGeometry(
        nodes=np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]]),
        connectivity=np.array([[0, 1, 2], [1, 3, 2]]),
        region_ids=np.array([0, 1]),
        plane="xz",
    )
    return PCSFModel(
        geometry=geometry,
        resistivity=np.array([120.0, 340.0]),
        resistivity_by_region=np.array([120.0, 340.0]),
        resistivity_by_node=np.array([100.0, 150.0, 200.0, 250.0]),
        source_backend="mare2dem",
    )


def _multiline_model() -> PCSFModel:
    line_geo = Grid2DGeometry(x=np.array([0.0, 100.0]), z=np.array([10.0, 50.0]))
    lines = [
        LineEntry(
            line_id="L1", geometry=line_geo, resistivity=np.array([[100.0, 110.0], [50.0, 55.0]]),
            offset_y=0.0, offset_kind="real", azimuth_deg=0.0,
        ),
        LineEntry(
            line_id="L2", geometry=line_geo, resistivity=np.array([[200.0, 210.0], [90.0, 95.0]]),
            offset_y=250.0, offset_kind="real", azimuth_deg=0.0,
        ),
    ]
    derived = DerivedVolume(
        grid=Grid3DGeometry(x=np.array([0.0, 100.0]), y=np.array([0.0, 125.0, 250.0]), z=np.array([10.0, 50.0])),
        resistivity=np.full((2, 3, 2), 150.0),
        derivation_method="linear_interp",
        derived_from=["L1", "L2"],
        synthesized=True,
    )
    geometry = MultilineGeometry(lines=lines, derived_volume=derived)
    return PCSFModel(geometry=geometry, source_backend="generic")


def _raster_model() -> PCSFModel:
    model = _grid2d_model()
    x = np.array([0.0, 250.0, 500.0])
    y = np.array([0.0, 300.0])
    elevation = np.array([[100.0, 110.0, 120.0], [130.0, 140.0, 150.0]])
    model.topography = TopographyRaster(x=x, y=y, elevation=elevation)
    return model


class TestGrid2DRoundTrip:
    def test_round_trip(self, tmp_path):
        model = _grid2d_model()
        path = write_pcsm(model, tmp_path / "grid2d.pcsm")
        restored = read_pcsm(path)

        assert restored.kind == "grid2d"
        assert restored.source_backend == "occam2d"
        np.testing.assert_allclose(restored.geometry.x, model.geometry.x)
        np.testing.assert_allclose(restored.geometry.z, model.geometry.z)
        np.testing.assert_allclose(restored.geometry.x_nodes, model.geometry.x_nodes)
        np.testing.assert_allclose(restored.geometry.origin, model.geometry.origin)
        assert restored.geometry.azimuth_deg == pytest.approx(45.0)
        np.testing.assert_allclose(restored.resistivity, model.resistivity)
        np.testing.assert_allclose(restored.resistivity_native, model.resistivity_native)
        assert restored.resistivity_native_encoding == "log10"
        np.testing.assert_allclose(restored.uncertainty, model.uncertainty)
        assert restored.stations.name == model.stations.name
        np.testing.assert_allclose(restored.stations.z, model.stations.z)
        np.testing.assert_allclose(restored.stations.lon, model.stations.lon)
        np.testing.assert_allclose(restored.stations.lat, model.stations.lat)
        assert restored.topography.station_id == model.topography.station_id
        np.testing.assert_allclose(restored.topography.elevation, model.topography.elevation)
        assert restored.survey == model.survey
        np.testing.assert_allclose(restored.history["rms"], model.history["rms"])
        assert restored.crs == "EPSG:32633"
        assert restored.description == "synthetic round-trip fixture"
        assert restored.metadata == {"note": "unit test"}

    def test_bit_exact_float_round_trip(self, tmp_path):
        # repr()-precision floats: PCSM must not lose precision the way
        # ModEM's own ~5-sig-fig ASCII .rho format does.
        model = _grid2d_model()
        model.resistivity = np.array(
            [[123.456789012345, 1.0 / 3.0, 1e-7, 481684.923817],
             [2.0, 3.0, 4.0, 5.0],
             [6.0, 7.0, 8.0, 9.0]]
        )
        path = write_pcsm(model, tmp_path / "precise.pcsm")
        restored = read_pcsm(path)
        assert np.array_equal(restored.resistivity, model.resistivity)

    def test_stations_without_lon_lat_round_trip_as_none(self, tmp_path):
        model = _grid2d_model()
        model.stations = StationTable(
            name=model.stations.name,
            x=model.stations.x,
            y=model.stations.y,
            z=model.stations.z,
        )
        path = write_pcsm(model, tmp_path / "no_lonlat.pcsm")
        restored = read_pcsm(path)
        assert restored.stations.lon is None
        assert restored.stations.lat is None

    def test_reads_legacy_five_column_stations_block(self, tmp_path):
        """A .pcsm written before lon/lat columns existed (a 5-field
        ``name x y z line_id`` STATIONS row, no trailing lon/lat) must
        still load, with lon/lat coming back as None -- not a parse
        error."""
        model = _grid2d_model()
        model.stations = StationTable(
            name=model.stations.name,
            x=model.stations.x,
            y=model.stations.y,
            z=model.stations.z,
        )
        path = write_pcsm(model, tmp_path / "legacy.pcsm")
        text = path.read_text(encoding="utf-8")
        lines = text.splitlines()
        out: list[str] = []
        in_stations = False
        for ln in lines:
            if ln.startswith("STATIONS"):
                in_stations = True
                out.append("STATIONS  # name x y z line_id")
                continue
            if ln.strip() == "END_STATIONS":
                in_stations = False
                out.append(ln)
                continue
            if in_stations:
                out.append(" ".join(ln.split()[:5]))
            else:
                out.append(ln)
        path.write_text("\n".join(out) + "\n", encoding="utf-8")

        restored = read_pcsm(path)
        assert restored.stations.name == model.stations.name
        np.testing.assert_allclose(restored.stations.x, model.stations.x)
        assert restored.stations.lon is None
        assert restored.stations.lat is None


class TestGrid3DRoundTrip:
    def test_round_trip(self, tmp_path):
        model = _grid3d_model()
        path = write_pcsm(model, tmp_path / "grid3d.pcsm")
        restored = read_pcsm(path)

        assert restored.kind == "grid3d"
        assert restored.geometry.n_air == 3
        assert restored.geometry.rotation_deg == pytest.approx(12.5)
        np.testing.assert_allclose(restored.geometry.origin, model.geometry.origin)
        assert restored.resistivity.shape == (4, 2, 3)
        np.testing.assert_allclose(restored.resistivity, model.resistivity)


class TestMeshUnstructuredRoundTrip:
    def test_round_trip(self, tmp_path):
        model = _mesh_model()
        path = write_pcsm(model, tmp_path / "mesh.pcsm")
        restored = read_pcsm(path)

        assert restored.kind == "mesh_unstructured"
        np.testing.assert_allclose(restored.geometry.nodes, model.geometry.nodes)
        np.testing.assert_array_equal(restored.geometry.connectivity, model.geometry.connectivity)
        np.testing.assert_array_equal(restored.geometry.region_ids, model.geometry.region_ids)
        np.testing.assert_allclose(restored.resistivity, model.resistivity)
        np.testing.assert_allclose(restored.resistivity_by_region, model.resistivity_by_region)
        np.testing.assert_allclose(restored.resistivity_by_node, model.resistivity_by_node)


class TestMultilineRoundTrip:
    def test_round_trip_preserves_line_order_and_derived_volume(self, tmp_path):
        model = _multiline_model()
        path = write_pcsm(model, tmp_path / "multiline.pcsm")
        restored = read_pcsm(path)

        assert restored.kind == "multiline"
        assert [line.line_id for line in restored.geometry.lines] == ["L1", "L2"]
        np.testing.assert_allclose(
            restored.geometry.lines[1].resistivity, model.geometry.lines[1].resistivity
        )
        assert restored.geometry.lines[1].offset_y == pytest.approx(250.0)
        assert restored.geometry.lines[1].offset_kind == "real"
        assert restored.geometry.derived_volume is not None
        assert restored.geometry.derived_volume.synthesized is True
        assert restored.geometry.derived_volume.derived_from == ["L1", "L2"]
        np.testing.assert_allclose(
            restored.geometry.derived_volume.resistivity, model.geometry.derived_volume.resistivity
        )
        assert restored.resistivity is None


class TestRasterTopographyRoundTrip:
    def test_round_trip(self, tmp_path):
        model = _raster_model()
        path = write_pcsm(model, tmp_path / "raster.pcsm")
        restored = read_pcsm(path)

        assert isinstance(restored.topography, TopographyRaster)
        np.testing.assert_allclose(restored.topography.x, model.topography.x)
        np.testing.assert_allclose(restored.topography.y, model.topography.y)
        np.testing.assert_allclose(restored.topography.elevation, model.topography.elevation)


class TestComments:
    def test_hand_annotated_file_parses_identically(self, tmp_path):
        model = _grid2d_model()
        path = write_pcsm(model, tmp_path / "grid2d.pcsm")
        original = path.read_text(encoding="utf-8")

        annotated_lines = []
        for line in original.splitlines():
            annotated_lines.append(f"# note about: {line.split()[0] if line.split() else ''}")
            annotated_lines.append(line)
        annotated_lines.append("# trailing comment at end of file")
        annotated = "\n".join(annotated_lines) + "\n"

        annotated_path = tmp_path / "annotated.pcsm"
        annotated_path.write_text(annotated, encoding="utf-8")

        restored_plain = read_pcsm(path)
        restored_annotated = read_pcsm(annotated_path)
        np.testing.assert_allclose(restored_annotated.resistivity, restored_plain.resistivity)
        assert restored_annotated.stations.name == restored_plain.stations.name

    def test_inline_trailing_comment_on_data_line_is_stripped(self, tmp_path):
        model = PCSFModel(
            geometry=Grid2DGeometry(x=np.array([0.0, 1.0]), z=np.array([0.0, 1.0])),
            resistivity=np.array([[10.0, 20.0], [30.0, 40.0]]),
        )
        path = write_pcsm(model, tmp_path / "m.pcsm")
        text = path.read_text(encoding="utf-8")
        text = text.replace("NX 2", "NX 2  # number of stations along the profile")
        path.write_text(text, encoding="utf-8")
        restored = read_pcsm(path)
        assert restored.geometry.x.shape == (2,)


class TestResistivityLabelsAndAlignment:
    def test_resistivity_blocks_are_self_labelled(self, tmp_path):
        model = _grid2d_model()  # carries resistivity_native (log10) already
        path = write_pcsm(model, tmp_path / "labelled.pcsm")
        text = path.read_text(encoding="utf-8")
        assert "RESISTIVITY  # linear ohm.m (canonical" in text
        assert "RESISTIVITY_NATIVE  # source-native encoding: log10" in text

    def test_log10_view_off_by_default(self, tmp_path):
        model = _grid2d_model()
        path = write_pcsm(model, tmp_path / "no_log10.pcsm")
        assert "RESISTIVITY_LOG10" not in path.read_text(encoding="utf-8")

    def test_log10_view_written_and_labelled(self, tmp_path):
        model = _grid2d_model()
        path = write_pcsm(model, tmp_path / "with_log10.pcsm", log10_view=True)
        text = path.read_text(encoding="utf-8")
        assert "RESISTIVITY_LOG10  # log10(ohm.m)" in text

    def test_log10_view_values_are_correct_and_discarded_on_read(self, tmp_path):
        model = _grid2d_model()
        path = write_pcsm(model, tmp_path / "with_log10.pcsm", log10_view=True)
        text = path.read_text(encoding="utf-8")

        # Extract the RESISTIVITY_LOG10 block's numbers directly from
        # the file text and check them against a fresh np.log10 call --
        # this must survive independently of read_pcsm, which discards
        # the block entirely (verified below).
        lines = text.splitlines()
        start = lines.index("RESISTIVITY_LOG10  # log10(ohm.m) -- derived view for human inspection only; read_pcsm() discards this, it is never part of the model") + 1
        end = lines.index("END_RESISTIVITY_LOG10")
        values = np.asarray(
            [float(tok) for line in lines[start:end] for tok in line.split()]
        )
        np.testing.assert_allclose(
            values, np.log10(model.resistivity).reshape(-1)
        )

        restored = read_pcsm(path)
        np.testing.assert_allclose(restored.resistivity, model.resistivity)
        assert not hasattr(restored, "resistivity_log10")

    def test_log10_view_does_not_break_later_optional_blocks(self, tmp_path):
        # Regression guard: RESISTIVITY_LOG10 sits between RESISTIVITY
        # and UNCERTAINTY/SENSITIVITY in the grammar -- the reader must
        # still find those via its peek-based sequential parse.
        model = _grid2d_model()
        assert model.uncertainty is not None
        path = write_pcsm(model, tmp_path / "log10_then_uncertainty.pcsm", log10_view=True)
        restored = read_pcsm(path)
        np.testing.assert_allclose(restored.uncertainty, model.uncertainty)
        assert restored.stations.name == model.stations.name

    def test_log10_view_for_multiline_per_line_and_round_trips(self, tmp_path):
        model = _multiline_model()
        path = write_pcsm(model, tmp_path / "multiline_log10.pcsm", log10_view=True)
        text = path.read_text(encoding="utf-8")
        # Count header lines only ("END_RESISTIVITY_LOG10" also
        # contains the substring "RESISTIVITY_LOG10").
        assert text.count("RESISTIVITY_LOG10  #") == len(model.geometry.lines)

        restored = read_pcsm(path)
        for orig_line, restored_line in zip(model.geometry.lines, restored.geometry.lines):
            np.testing.assert_allclose(restored_line.resistivity, orig_line.resistivity)

    def test_pcsf_to_pcsm_passes_through_log10_view(self, tmp_path):
        model = _grid2d_model()
        pcsf_path = write_pcsf(model, tmp_path / "a.pcsf")
        pcsm_path = pcsf_to_pcsm(pcsf_path, tmp_path / "a.pcsm", log10_view=True)
        assert "RESISTIVITY_LOG10" in pcsm_path.read_text(encoding="utf-8")

    def test_coordinate_block_values_are_column_aligned(self, tmp_path):
        # Every token within one block shares a common right-justified
        # width, so full rows (all but a possibly-shorter trailing
        # partial row) come out the same length -- not ragged like
        # plain unpadded repr() would produce.
        model = PCSFModel(
            geometry=Grid2DGeometry(
                x=np.array([0.0, 1.5, 123456.789, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]),
                z=np.array([1.0, 2.0]),
            ),
            resistivity=np.ones((2, 9)) * 100.0,
        )
        path = write_pcsm(model, tmp_path / "aligned.pcsm")
        lines = path.read_text(encoding="utf-8").splitlines()
        start = lines.index("X_COORDS") + 1
        end = lines.index("END_X_COORDS")
        rows = lines[start:end]
        full_rows = [r for r in rows if len(r.split()) == _DEFAULT_ROW_WIDTH]
        assert len(full_rows) >= 1
        # Equal token count -> equal total length implies equal
        # per-token width (fixed single-space separators). ``.split()``
        # itself would swallow the left-padding, so length is the
        # right thing to compare here, not the split tokens.
        row_lengths = {len(line) for line in full_rows}
        assert len(row_lengths) == 1, f"full rows are not equal width: {full_rows}"


class TestGzipSupport:
    def test_gz_suffix_round_trips(self, tmp_path):
        model = _grid2d_model()
        path = write_pcsm(model, tmp_path / "grid2d.pcsm.gz")
        restored = read_pcsm(path)
        np.testing.assert_allclose(restored.resistivity, model.resistivity)
        assert restored.stations.name == model.stations.name

    def test_gz_file_is_smaller_than_plain(self, tmp_path):
        model = _grid3d_model()
        plain = write_pcsm(model, tmp_path / "grid3d.pcsm")
        gz = write_pcsm(model, tmp_path / "grid3d.pcsm.gz")
        assert gz.stat().st_size < plain.stat().st_size

    def test_gz_file_is_not_plain_text(self, tmp_path):
        model = _grid2d_model()
        path = write_pcsm(model, tmp_path / "grid2d.pcsm.gz")
        with pytest.raises(UnicodeDecodeError):
            path.read_text(encoding="utf-8")


class TestConversion:
    def test_pcsf_to_pcsm_and_back(self, tmp_path):
        model = _grid2d_model()
        pcsf_path = write_pcsf(model, tmp_path / "a.pcsf")
        pcsm_path = pcsf_to_pcsm(pcsf_path, tmp_path / "a.pcsm")
        restored_pcsm = read_pcsm(pcsm_path)
        np.testing.assert_allclose(restored_pcsm.resistivity, model.resistivity)

        pcsf_back_path = pcsm_to_pcsf(pcsm_path, tmp_path / "a_back.pcsf")
        restored_pcsf = read_pcsf(pcsf_back_path)
        np.testing.assert_allclose(restored_pcsf.resistivity, model.resistivity)
        assert restored_pcsf.stations.name == model.stations.name


class TestErrors:
    def test_write_rejects_non_pcsf_model(self, tmp_path):
        with pytest.raises(TypeError):
            write_pcsm(object(), tmp_path / "bad.pcsm")

    def test_write_rejects_invalid_model(self, tmp_path):
        model = PCSFModel(geometry=Grid2DGeometry(x=np.array([0.0]), z=np.array([0.0])))
        with pytest.raises(ValueError, match="resistivity is required"):
            write_pcsm(model, tmp_path / "invalid.pcsm")

    def test_read_rejects_missing_pcsm_version(self, tmp_path):
        path = tmp_path / "no_version.pcsm"
        path.write_text("SOURCE_BACKEND occam2d\n", encoding="utf-8")
        with pytest.raises(ValueError, match="PCSM_VERSION"):
            read_pcsm(path)

    def test_read_rejects_unrecognised_major_version(self, tmp_path):
        model = _grid2d_model()
        path = write_pcsm(model, tmp_path / "future.pcsm")
        text = path.read_text(encoding="utf-8").replace(
            f"PCSM_VERSION {PCSF_VERSION}", "PCSM_VERSION 99.0.0"
        )
        path.write_text(text, encoding="utf-8")
        with pytest.raises(ValueError, match="unsupported pcsf_version"):
            read_pcsm(path)

    def test_read_warns_on_newer_minor_version(self, tmp_path):
        current_major, current_minor, _ = (int(p) for p in PCSF_VERSION.split("."))
        model = _grid2d_model()
        path = write_pcsm(model, tmp_path / "newer_minor.pcsm")
        text = path.read_text(encoding="utf-8").replace(
            f"PCSM_VERSION {PCSF_VERSION}",
            f"PCSM_VERSION {current_major}.{current_minor + 1}.0",
        )
        path.write_text(text, encoding="utf-8")
        with pytest.warns(UserWarning, match="newer than this reader's"):
            restored = read_pcsm(path)
        assert restored.kind == "grid2d"

    def test_read_rejects_mismatched_array_length(self, tmp_path):
        model = _grid2d_model()
        path = write_pcsm(model, tmp_path / "bad_len.pcsm")
        text = path.read_text(encoding="utf-8").replace("NX 4", "NX 5")
        path.write_text(text, encoding="utf-8")
        with pytest.raises(ValueError, match="X_COORDS"):
            read_pcsm(path)

    def test_write_creates_missing_parent_directories(self, tmp_path):
        model = _grid2d_model()
        nested = tmp_path / "a" / "b" / "c.pcsm"
        path = write_pcsm(model, nested)
        assert path.exists()
