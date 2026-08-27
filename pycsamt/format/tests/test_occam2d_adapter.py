# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Phase 2 tests — the Occam2D adapter, against the real bundled
``data/occam2D`` dataset (47 sites, 17 frequencies, RMS ~ 0.998).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from pycsamt.format import read_pcsf, write_pcsf
from pycsamt.format.adapters.occam2d import occam2d_to_pcsf
from pycsamt.interp._base import ResistivityModel
from pycsamt.models.occam2d.results import InversionResult

_DATA_DIR = Path(__file__).parents[3] / "data" / "occam2D"
_SKIP = pytest.mark.skipif(
    not _DATA_DIR.exists(),
    reason=f"bundled occam2D data not found: {_DATA_DIR}",
)


@pytest.fixture(scope="module")
def result() -> InversionResult:
    return InversionResult(workdir=_DATA_DIR)


@_SKIP
class TestOccam2DAdapter:
    def test_geometry_kind_and_shape(self, result):
        model = occam2d_to_pcsf(result)
        assert model.kind == "grid2d"
        assert model.resistivity.shape == (
            result.mesh.n_zcells,
            result.mesh.n_xcells,
        )
        model.validate()

    def test_resistivity_matches_resistivity_model_cross_check(self, result):
        model = occam2d_to_pcsf(result)
        rm = ResistivityModel.from_occam2d(result)
        np.testing.assert_allclose(model.resistivity, 10.0**rm.rho_2d)
        np.testing.assert_allclose(model.resistivity_native, rm.rho_2d)
        np.testing.assert_allclose(model.geometry.x, rm.x_centers)
        np.testing.assert_allclose(model.geometry.z, rm.z_centers)
        assert model.resistivity_native_encoding == "log10"

    def test_x_nodes_bracket_x_centers_after_shift(self, result):
        # x_nodes must stay in the same real-chainage frame as x
        # (the cell-edge coordinates the ResistivityModel adapter does
        # not itself carry), so every centre sits strictly between its
        # neighbouring edges.
        model = occam2d_to_pcsf(result)
        x, x_nodes = model.geometry.x, model.geometry.x_nodes
        assert np.all(x >= x_nodes[:-1] - 1e-6)
        assert np.all(x <= x_nodes[1:] + 1e-6)

    def test_station_x_lands_inside_the_shifted_frame(self, result):
        # Regression guard for the exact bug ResistivityModel.from_occam2d
        # documents: an unshifted mesh-local x would put every station on
        # the same outermost padding column.
        model = occam2d_to_pcsf(result)
        sta_x = set(np.round(model.stations.x, 3))
        assert len(sta_x) == len(model.stations.name)

    def test_history_from_log(self, result):
        # A raw, lossless copy of OccamLog's own arrays: the log's
        # per-iteration RMS is not guaranteed to end on
        # InversionResult.final_rms (the loaded .iter file's own
        # embedded misfit_value) — Occam's log-iteration index
        # commonly describes the step that *produced* the next iter
        # file, so the two only need to agree at the iteration the log
        # itself reports as best, not at the array's last entry.
        model = occam2d_to_pcsf(result)
        assert set(model.history) == {
            "iteration", "rms", "roughness", "lagrange", "stepsize",
        }
        assert model.history["rms"].shape == (result.log.iterations.size,)
        best_ix = int(result.log.iterations.tolist().index(result.log.best_iteration))
        np.testing.assert_allclose(
            model.history["rms"][best_ix], result.final_rms
        )

    def test_station_elevations_populate_topography(self, result):
        names = list(InversionResult(workdir=_DATA_DIR).data.sites)
        elev = {names[0]: 1200.0, names[1]: 1195.5}
        model = occam2d_to_pcsf(result, station_elevations=elev)
        assert model.topography is not None
        assert model.topography.station_id[0] == names[0]
        assert model.topography.elevation[0] == pytest.approx(1200.0)
        # Stations without a known elevation are nan, not fabricated 0.
        assert np.isnan(model.topography.elevation[2])
        assert np.isnan(model.stations.z[2])

    def test_no_elevations_means_no_topography_group(self, result):
        model = occam2d_to_pcsf(result)
        assert model.topography is None
        assert model.stations is not None
        assert np.all(np.isnan(model.stations.z))

    def test_station_lonlat_populates_stations(self, result):
        names = list(InversionResult(workdir=_DATA_DIR).data.sites)
        lonlat = {names[0]: (7.5, 45.1), names[1]: (7.6, 45.2)}
        model = occam2d_to_pcsf(result, station_lonlat=lonlat)
        assert model.stations.lon is not None
        assert model.stations.lon[0] == pytest.approx(7.5)
        assert model.stations.lat[0] == pytest.approx(45.1)
        # Stations without a known lon/lat are nan, not fabricated 0.
        assert np.isnan(model.stations.lon[2])

    def test_no_lonlat_means_stations_lon_lat_stay_none(self, result):
        model = occam2d_to_pcsf(result)
        assert model.stations.lon is None
        assert model.stations.lat is None

    def test_station_lonlat_survives_pcsf_round_trip(self, result, tmp_path):
        names = list(InversionResult(workdir=_DATA_DIR).data.sites)
        lonlat = {names[0]: (7.5, 45.1)}
        model = occam2d_to_pcsf(result, station_lonlat=lonlat)
        path = write_pcsf(model, tmp_path / "occam_lonlat.pcsf")
        restored = read_pcsf(path)
        assert restored.stations.lon[0] == pytest.approx(7.5)
        assert restored.stations.lat[0] == pytest.approx(45.1)

    def test_topo_csv_overrides_station_lonlat_with_warning(self, result, tmp_path):
        names = list(InversionResult(workdir=_DATA_DIR).data.sites)
        csv_path = tmp_path / "topo.csv"
        csv_path.write_text(
            "station,lat,lon,elev\n"
            + "\n".join(
                f"{n},{45.0 + i * 0.001},{7.5 + i * 0.001},{1000.0 + i}"
                for i, n in enumerate(names)
            )
            + "\n",
            encoding="utf-8",
        )
        with pytest.warns(UserWarning, match="takes precedence"):
            model = occam2d_to_pcsf(
                result,
                topo=csv_path,
                station_lonlat={names[0]: (0.0, 0.0)},
                station_elevations={names[0]: -999.0},
            )
        assert model.stations.lon[0] == pytest.approx(7.5)
        assert model.stations.lat[0] == pytest.approx(45.0)
        assert model.stations.z[0] == pytest.approx(1000.0)
        assert model.topography is not None

    def test_topo_bln_positional_requires_exact_station_count(self, result, tmp_path):
        n = len(InversionResult(workdir=_DATA_DIR).data.sites)
        bln_path = tmp_path / "topo_wrong_count.bln"
        bln_path.write_text("2,1\n500000,4500000\n500010,4500000\n", encoding="utf-8")
        with pytest.raises(ValueError, match=f"{n} station"):
            occam2d_to_pcsf(result, topo=bln_path, epsg=32633)

    def test_topo_bln_positional_matches_survey_order_with_utm(self, result):
        names = list(InversionResult(workdir=_DATA_DIR).data.sites)
        n = len(names)
        rows = "\n".join(f"{500000 + 10 * i},4500000,{900.0 + i}" for i in range(n))
        import tempfile
        from pathlib import Path

        with tempfile.TemporaryDirectory() as td:
            bln_path = Path(td) / "topo.bln"
            bln_path.write_text(f"{n},1\n{rows}\n", encoding="utf-8")
            model = occam2d_to_pcsf(result, topo=bln_path, epsg=32633)
        assert model.stations.lon[0] == pytest.approx(15.0, abs=1e-3)
        assert model.stations.lat[0] == pytest.approx(40.650857, abs=1e-3)
        # station-count many points, each shifted 10 m east -> increasing lon
        assert np.all(np.diff(model.stations.lon) > 0)
        np.testing.assert_allclose(model.stations.z, 900.0 + np.arange(n))

    def test_no_topo_behaves_exactly_as_before(self, result):
        model = occam2d_to_pcsf(result)
        assert model.stations.lon is None
        assert model.stations.lat is None
        assert model.topography is None

    def test_survey_accepts_plain_mapping(self, result):
        model = occam2d_to_pcsf(result, survey={"name": "Tongkeng CSAMT"})
        assert model.survey == {"name": "Tongkeng CSAMT"}

    def test_survey_accepts_to_dict_object(self, result):
        class _Fake:
            def to_dict(self):
                return {"name": "from-object"}

        model = occam2d_to_pcsf(result, survey=_Fake())
        assert model.survey == {"name": "from-object"}

    def test_metadata_carries_provenance(self, result):
        model = occam2d_to_pcsf(result)
        assert model.metadata["n_iterations"] == result.n_iterations
        assert model.metadata["final_rms"] == pytest.approx(result.final_rms)
        assert model.metadata["converged"] is True

    def test_full_round_trip_through_pcsf_file(self, result, tmp_path):
        model = occam2d_to_pcsf(
            result,
            station_elevations={"S00": 1200.0},
            created_by="pytest",
            crs="local",
            description="Tongkeng CSAMT round trip",
        )
        path = write_pcsf(model, tmp_path / "occam2d_real.pcsf")
        restored = read_pcsf(path)

        assert restored.source_backend == "occam2d"
        np.testing.assert_allclose(restored.resistivity, model.resistivity)
        np.testing.assert_allclose(restored.geometry.x_nodes, model.geometry.x_nodes)
        assert restored.stations.name == model.stations.name
        assert restored.topography.station_id[0] == "S00"
        np.testing.assert_allclose(restored.history["rms"], model.history["rms"])
        assert restored.metadata["n_iterations"] == result.n_iterations

    def test_raises_without_rho_2d_or_mesh(self):
        empty = InversionResult.__new__(InversionResult)
        empty.rho_2d = None
        empty.mesh = None
        with pytest.raises(ValueError, match="rho_2d or mesh"):
            occam2d_to_pcsf(empty)
