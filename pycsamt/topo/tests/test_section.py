# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for pycsamt.topo.section — the one-call topo-embedded section API.

Covers the model adapters (raw arrays, ResistivityModel, backend-neutral
InversionResult, native Occam2D / ModEM duck-typed results, AI agent
results), the topography-source resolver (sites / array / model /
flat-fallback), depth cropping, and both plotting modes
(``pcolormesh`` and ``imshow``) — including an end-to-end run against
real WILLY AMT EDI data.
"""

from __future__ import annotations

import glob
import os
import warnings

import numpy as np
import pytest

matplotlib = pytest.importorskip("matplotlib")
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from pycsamt.agents._base import AgentResult
from pycsamt.interp import ResistivityModel
from pycsamt.inversion.results import InversionResult
from pycsamt.seg.collection import EDICollection
from pycsamt.seg.edi import EDIFile
from pycsamt.topo import build_topo_section, plot_topo_section, reset_topo
from pycsamt.topo.section import TopoSection, _cell_edges, _extract_grid, _infer_terrain_km

# ---------------------------------------------------------------------------
# Real WILLY EDI data (shared with the rest of pycsamt/topo/tests/)
# ---------------------------------------------------------------------------

_WILLY_ROOT = os.path.join(
    os.path.dirname(__file__), "..", "..", "..", "data", "AMT", "WILLY_DATA"
)


def _load(profile: str = "L18PLT") -> list[EDIFile]:
    edi_dir = os.path.join(_WILLY_ROOT, profile)
    paths = sorted(glob.glob(os.path.join(edi_dir, "*.edi")))
    if not paths:
        pytest.skip(f"WILLY {profile} EDI files not found at {edi_dir}")
    return [EDIFile(p) for p in paths]


@pytest.fixture(autouse=True)
def _reset_after():
    yield
    reset_topo()
    plt.close("all")


# ---------------------------------------------------------------------------
# Synthetic grid helpers
# ---------------------------------------------------------------------------


def _grid(n_x=10, n_z=20, seed=0, x_km=5.0, z_km=2.0):
    rng = np.random.default_rng(seed)
    x_c = np.linspace(0.0, x_km, n_x)
    z_c = np.linspace(z_km / n_z / 2, z_km, n_z)
    rho = rng.random((n_z, n_x)) * 3.0
    return x_c, z_c, rho


class _FakeMesh:
    def __init__(self, x_nodes, z_nodes):
        self.x_nodes = x_nodes
        self.z_nodes = z_nodes


class _FakeOccamData:
    def __init__(self, offsets, sites):
        self.offsets = offsets
        self.sites = sites


class _FakeOccam2DResult:
    """Minimal duck-type of pycsamt.models.occam2d.results.InversionResult."""

    def __init__(self, rho_2d, x_nodes, z_nodes, station_names, final_rms=1.1):
        self.rho_2d = rho_2d
        self.mesh = _FakeMesh(x_nodes, z_nodes)
        self.data = _FakeOccamData(
            0.5 * (x_nodes[:-1] + x_nodes[1:]), station_names
        )
        self.final_rms = final_rms


class _FakeModEmModel2D:
    def __init__(self, x_widths, z_widths, rho_loge):
        self.x_widths = x_widths
        self.z_widths = z_widths
        self.rho_loge = rho_loge


class _FakeModEmResult:
    """Minimal duck-type of pycsamt.models.modem.results.InversionResult."""

    def __init__(self, mode, model_final):
        self.mode = mode
        self.model_final = model_final
        self.model_initial = None


# ---------------------------------------------------------------------------
# _cell_edges
# ---------------------------------------------------------------------------


class TestCellEdges:
    def test_multi_point_edges_bracket_centres(self):
        centres = np.array([1.0, 2.0, 4.0, 8.0])
        edges = _cell_edges(centres)
        assert edges.shape == (5,)
        assert np.all(np.diff(edges) > 0)
        assert edges[0] < centres[0]
        assert edges[-1] > centres[-1]

    def test_single_point(self):
        edges = _cell_edges(np.array([3.0]))
        assert edges.shape == (2,)
        assert edges[0] < 3.0 < edges[1]

    def test_empty(self):
        edges = _cell_edges(np.array([]))
        assert edges.shape == (2,)


# ---------------------------------------------------------------------------
# _infer_terrain_km
# ---------------------------------------------------------------------------


class TestInferTerrainKm:
    def test_no_air_cells_returns_none(self):
        _, z_c, rho = _grid()
        assert _infer_terrain_km(z_c, rho, threshold=5.0) is None

    def test_uniform_earth_at_surface_returns_none(self):
        """Every column starting in 'earth' means no real air cap exists."""
        _, z_c, rho = _grid()
        rho[:, :] = 1.0  # well below threshold everywhere, including row 0
        assert _infer_terrain_km(z_c, rho, threshold=5.0) is None

    def test_varying_air_cap_gives_relief(self):
        _, z_c, rho = _grid(n_x=6, n_z=10)
        rho[:, :] = 1.0
        rho[0, :3] = 6.0  # shallow air only on the left half
        terrain = _infer_terrain_km(z_c, rho, threshold=5.0)
        assert terrain is not None
        assert terrain.shape == (6,)
        # left columns (air at row 0) sit lower than right columns (no air)
        assert terrain[0] < terrain[-1]


# ---------------------------------------------------------------------------
# _extract_grid adapters
# ---------------------------------------------------------------------------


class TestExtractGridAdapters:
    def test_raw_tuple(self):
        x_c, z_c, rho = _grid()
        info = _extract_grid((x_c, z_c, rho))
        np.testing.assert_allclose(info.x_centers, x_c)
        np.testing.assert_allclose(info.z_centers, z_c)
        assert info.rho_log10.shape == rho.shape
        assert info.method == "array"
        assert len(info.station_names) == x_c.size

    def test_raw_list_also_accepted(self):
        x_c, z_c, rho = _grid()
        info = _extract_grid([x_c, z_c, rho])
        assert info.method == "array"

    def test_bare_ndarray_raises(self):
        with pytest.raises(TypeError, match="bare 2-D array"):
            _extract_grid(np.zeros((5, 5)))

    def test_resistivity_model(self):
        x_c, z_c, rho = _grid()
        rm = ResistivityModel.from_array(
            rho, x_c, z_c, method="occam2d", rms=0.9
        )
        info = _extract_grid(rm)
        assert info.method == "occam2d"
        assert info.rms == pytest.approx(0.9)
        np.testing.assert_allclose(info.rho_log10, rho)

    def test_backend_neutral_inversion_result(self):
        x_c, z_c, rho = _grid()
        ir = InversionResult(
            method="mt",
            dimension="2d",
            backend="builtin",
            rms=0.42,
            model={"rho_2d": rho, "x_centers": x_c, "z_centers": z_c},
        )
        info = _extract_grid(ir)
        assert info.method == "builtin:mt"
        assert info.rms == pytest.approx(0.42)
        np.testing.assert_allclose(info.rho_log10, rho)

    def test_native_occam2d_duck_type(self):
        x_c, z_c, rho = _grid()
        x_nodes = _cell_edges(x_c)
        z_nodes = _cell_edges(z_c)
        names = [f"O{i:02d}" for i in range(x_c.size)]
        fake = _FakeOccam2DResult(rho, x_nodes, z_nodes, names, final_rms=1.7)
        info = _extract_grid(fake)
        assert info.method == "occam2d"
        assert info.rms == pytest.approx(1.7)
        assert info.rho_log10.shape == rho.shape

    def test_native_modem_2d_duck_type(self):
        x_c, z_c, rho = _grid()
        x_widths = np.diff(_cell_edges(x_c))
        z_widths = np.diff(_cell_edges(z_c))
        model = _FakeModEmModel2D(x_widths, z_widths, np.log(10.0**rho))
        fake = _FakeModEmResult("2d", model)
        info = _extract_grid(fake)
        assert info.method == "modem"
        np.testing.assert_allclose(info.rho_log10, rho, atol=1e-8)

    def test_native_modem_3d_raises_with_guidance(self):
        x_c, z_c, rho = _grid()
        x_widths = np.diff(_cell_edges(x_c))
        z_widths = np.diff(_cell_edges(z_c))
        model = _FakeModEmModel2D(x_widths, z_widths, np.log(10.0**rho))
        fake = _FakeModEmResult("3d", model)
        with pytest.raises(ValueError, match="station_curtain|PlotSection"):
            _extract_grid(fake)

    def test_ai_dict_result(self):
        x_c, z_c, rho = _grid(n_x=8, n_z=12)
        payload = {
            "pred_rho": rho.T,
            "depths_km": z_c,
            "station_names": [f"AI{i}" for i in range(8)],
            "rms_global": 0.33,
        }
        info = _extract_grid(payload)
        assert info.method == "ai"
        assert info.unit == "km"
        assert info.rms == pytest.approx(0.33)
        np.testing.assert_allclose(info.rho_log10, rho)

    def test_ai_agent_result_object(self):
        x_c, z_c, rho = _grid(n_x=5, n_z=6)
        result = AgentResult(
            status="success",
            summary="ok",
            data={
                "pred_rho": rho.T,
                "depths_km": z_c,
                "station_names": [f"S{i}" for i in range(5)],
                "station_coords": np.column_stack(
                    [x_c * 1000.0, np.zeros(5)]
                ),
            },
        )
        info = _extract_grid(result)
        assert info.method == "ai"
        assert info.x_centers.size == 5
        # cumulative-distance chainage should be increasing
        assert np.all(np.diff(info.x_centers) >= 0)

    def test_mare2dem_native_raises_not_implemented(self):
        class _FakeMare:
            pass

        _FakeMare.__module__ = "pycsamt.models.mare2dem.mesh"
        with pytest.raises(NotImplementedError, match="triangular mesh"):
            _extract_grid(_FakeMare())

    def test_unsupported_type_raises(self):
        with pytest.raises(TypeError, match="Unsupported model type"):
            _extract_grid(object())


# ---------------------------------------------------------------------------
# build_topo_section — topography source resolution
# ---------------------------------------------------------------------------


class TestTopoSourceResolution:
    def test_explicit_array_source(self):
        x_c, z_c, rho = _grid(n_x=6)
        elev = np.linspace(100.0, 160.0, 6)
        sec = build_topo_section(
            (x_c, z_c, rho), elevation=elev, chainage=x_c, topo_source="array"
        )
        assert sec.topo_source == "array"
        np.testing.assert_allclose(sec.elev_km, elev / 1000.0)

    def test_array_source_requires_elevation(self):
        x_c, z_c, rho = _grid(n_x=6)
        with pytest.raises(ValueError, match="topo_source='array'"):
            build_topo_section((x_c, z_c, rho), topo_source="array")

    def test_sites_source_requires_sites(self):
        x_c, z_c, rho = _grid(n_x=6)
        with pytest.raises(ValueError, match="topo_source='sites'"):
            build_topo_section((x_c, z_c, rho), topo_source="sites")

    def test_invalid_topo_source_raises(self):
        x_c, z_c, rho = _grid(n_x=6)
        with pytest.raises(ValueError, match="topo_source"):
            build_topo_section((x_c, z_c, rho), topo_source="bogus")

    def test_auto_prefers_sites_over_array(self):
        edis = _load("L18PLT")
        x_c, z_c, rho = _grid(n_x=len(edis))
        sec = build_topo_section(
            (x_c, z_c, rho),
            sites=edis,
            elevation=np.zeros(len(edis)),
        )
        assert sec.topo_source == "sites"

    def test_flat_fallback_warns(self):
        x_c, z_c, rho = _grid(n_x=6)
        with pytest.warns(UserWarning, match="flat datum"):
            sec = build_topo_section((x_c, z_c, rho))
        assert sec.topo_source == "flat"
        np.testing.assert_allclose(sec.elev_km, 0.0)

    def test_model_source_forced_without_air_cells_raises(self):
        x_c, z_c, rho = _grid(n_x=6)
        with pytest.raises(ValueError, match="topo_source='model'"):
            build_topo_section((x_c, z_c, rho), topo_source="model")


# ---------------------------------------------------------------------------
# build_topo_section — depth cropping, units, smoothing, log_rho
# ---------------------------------------------------------------------------


class TestBuildTopoSectionControls:
    def test_depth_crop_reduces_rows(self):
        x_c, z_c, rho = _grid(n_x=6, n_z=20, z_km=2.0)
        elev = np.linspace(100.0, 160.0, 6)
        full = build_topo_section((x_c, z_c, rho), elevation=elev, chainage=x_c)
        cropped = build_topo_section(
            (x_c, z_c, rho),
            elevation=elev,
            chainage=x_c,
            depth_max=1.0,
            model_unit="km",
        )
        assert cropped.values.shape[0] < full.values.shape[0]
        assert cropped.z_centers_km.max() <= 1.0 + 1e-9

    def test_depth_window_selecting_nothing_falls_back_with_warning(self):
        x_c, z_c, rho = _grid(n_x=6, n_z=20, z_km=2.0)
        elev = np.linspace(100.0, 160.0, 6)
        with pytest.warns(UserWarning, match="selects no layers"):
            sec = build_topo_section(
                (x_c, z_c, rho),
                elevation=elev,
                chainage=x_c,
                depth_min=5.0,
                depth_max=6.0,
                model_unit="km",
            )
        assert sec.values.shape[0] == rho.shape[0]

    def test_model_unit_metres_vs_km_equivalent(self):
        x_c, z_c, rho = _grid(n_x=6, n_z=10, x_km=5.0, z_km=2.0)
        elev = np.linspace(100.0, 160.0, 6)
        sec_km = build_topo_section(
            (x_c, z_c, rho), elevation=elev, chainage=x_c, model_unit="km"
        )
        sec_m = build_topo_section(
            (x_c * 1000.0, z_c * 1000.0, rho),
            elevation=elev,
            chainage=x_c,
            model_unit="m",
        )
        np.testing.assert_allclose(sec_km.x_centers_km, sec_m.x_centers_km)
        np.testing.assert_allclose(sec_km.z_centers_km, sec_m.z_centers_km)

    def test_log_rho_false_converts_to_linear(self):
        x_c, z_c, rho = _grid(n_x=6, n_z=10)
        elev = np.linspace(100.0, 160.0, 6)
        sec = build_topo_section(
            (x_c, z_c, rho), elevation=elev, chainage=x_c, log_rho=False
        )
        assert np.nanmin(sec.values) > 0.0  # linear resistivity is positive
        assert sec.log_rho is False

    def test_smooth_sigma_applies_without_error(self):
        x_c, z_c, rho = _grid(n_x=8, n_z=12)
        elev = np.linspace(100.0, 160.0, 8)
        sec = build_topo_section(
            (x_c, z_c, rho),
            elevation=elev,
            chainage=x_c,
            smooth_sigma=(1.0, 1.0),
        )
        assert sec.values.shape[1] == 8

    def test_exaggeration_widens_draped_relief(self):
        x_c, z_c, rho = _grid(n_x=6, n_z=8)
        elev = np.linspace(100.0, 300.0, 6)
        sec1 = build_topo_section(
            (x_c, z_c, rho), elevation=elev, chainage=x_c, exaggeration=1.0
        )
        sec2 = build_topo_section(
            (x_c, z_c, rho), elevation=elev, chainage=x_c, exaggeration=3.0
        )
        range1 = sec1.z_draped_km[0].max() - sec1.z_draped_km[-1].min()
        range2 = sec2.z_draped_km[0].max() - sec2.z_draped_km[-1].min()
        assert range2 > range1

    def test_returns_topo_section_dataclass(self):
        x_c, z_c, rho = _grid()
        sec = build_topo_section(
            (x_c, z_c, rho), elevation=np.full(x_c.size, 100.0), chainage=x_c
        )
        assert isinstance(sec, TopoSection)
        assert sec.z_draped_km.shape == (
            sec.z_nodes_km.size,
            sec.x_nodes_km.size,
        )


# ---------------------------------------------------------------------------
# plot_topo_section — synthetic data
# ---------------------------------------------------------------------------


class TestPlotTopoSectionSynthetic:
    def test_pcolormesh_default(self):
        x_c, z_c, rho = _grid(n_x=8, n_z=15)
        elev = np.linspace(80.0, 220.0, 8)
        ax = plot_topo_section((x_c, z_c, rho), elevation=elev, chainage=x_c)
        assert ax.get_xlabel() == "Profile distance (km)"
        assert ax.get_ylabel() == "Elevation (km)"
        ax.figure.canvas.draw()

    def test_imshow_mode(self):
        x_c, z_c, rho = _grid(n_x=8, n_z=15)
        elev = np.linspace(80.0, 220.0, 8)
        ax = plot_topo_section(
            (x_c, z_c, rho), elevation=elev, chainage=x_c, kind="imshow"
        )
        assert ax.get_ylabel() == "Depth (km)"
        ax.figure.canvas.draw()

    def test_invalid_kind_raises(self):
        x_c, z_c, rho = _grid()
        with pytest.raises(ValueError, match="kind must be"):
            plot_topo_section(
                (x_c, z_c, rho),
                elevation=np.zeros(x_c.size),
                chainage=x_c,
                kind="bogus",
            )

    def test_return_data_returns_section(self):
        x_c, z_c, rho = _grid()
        ax, data = plot_topo_section(
            (x_c, z_c, rho),
            elevation=np.full(x_c.size, 120.0),
            chainage=x_c,
            return_data=True,
        )
        assert isinstance(data, TopoSection)
        assert ax is not None

    def test_draws_into_existing_axes(self):
        x_c, z_c, rho = _grid()
        fig, ax = plt.subplots()
        out_ax = plot_topo_section(
            (x_c, z_c, rho),
            elevation=np.full(x_c.size, 120.0),
            chainage=x_c,
            ax=ax,
        )
        assert out_ax is ax

    def test_show_stations_false_skips_pins(self):
        x_c, z_c, rho = _grid(n_x=6)
        elev = np.linspace(80.0, 200.0, 6)
        ax = plot_topo_section(
            (x_c, z_c, rho),
            elevation=elev,
            chainage=x_c,
            show_stations=False,
        )
        ax.figure.canvas.draw()

    def test_explicit_vmin_vmax_used(self):
        x_c, z_c, rho = _grid()
        ax, data = plot_topo_section(
            (x_c, z_c, rho),
            elevation=np.full(x_c.size, 120.0),
            chainage=x_c,
            vmin=0.0,
            vmax=1.0,
            return_data=True,
        )
        mesh = ax.collections[0]
        assert mesh.get_clim() == (0.0, 1.0)

    def test_savepath_writes_file(self, tmp_path):
        x_c, z_c, rho = _grid()
        out = tmp_path / "topo_section"
        plot_topo_section(
            (x_c, z_c, rho),
            elevation=np.full(x_c.size, 120.0),
            chainage=x_c,
            savepath=str(out),
        )
        written = list(tmp_path.glob("topo_section*"))
        assert written, "expected at least one saved figure file"


# ---------------------------------------------------------------------------
# plot_topo_section's default styling: white air fill, white markers
# ---------------------------------------------------------------------------


class TestPlotTopoSectionDefaultStyling:
    def test_pcolormesh_air_fill_is_transparent_by_default(self):
        x_c, z_c, rho = _grid(n_x=8, n_z=15)
        elev = np.linspace(80.0, 220.0, 8)
        ax = plot_topo_section((x_c, z_c, rho), elevation=elev, chainage=x_c)
        alphas = [p.get_alpha() for p in ax.patches]
        assert alphas and all(a == 0.0 for a in alphas)

    def test_pcolormesh_markers_are_white_faced_by_default(self):
        x_c, z_c, rho = _grid(n_x=8, n_z=15)
        elev = np.linspace(80.0, 220.0, 8)
        ax = plot_topo_section((x_c, z_c, rho), elevation=elev, chainage=x_c)
        scatter = ax.collections[-1]
        assert tuple(scatter.get_facecolor()[0]) == (1.0, 1.0, 1.0, 1.0)
        assert tuple(scatter.get_edgecolor()[0]) == (0.0, 0.0, 0.0, 1.0)

    @staticmethod
    def _find_strip_axes(ax):
        for a in ax.figure.axes:
            if a is not ax and a.get_ylabel() == "Elev (m)":
                return a
        raise AssertionError("strip axes (ylabel='Elev (m)') not found")

    def test_imshow_markers_are_white_faced_by_default(self):
        x_c, z_c, rho = _grid(n_x=8, n_z=15)
        elev = np.linspace(80.0, 220.0, 8)
        ax = plot_topo_section(
            (x_c, z_c, rho), elevation=elev, chainage=x_c, kind="imshow"
        )
        strip_ax = self._find_strip_axes(ax)
        scatter = strip_ax.collections[-1]
        assert tuple(scatter.get_facecolor()[0]) == (1.0, 1.0, 1.0, 1.0)

    def test_imshow_strip_background_is_transparent_by_default(self):
        x_c, z_c, rho = _grid(n_x=8, n_z=15)
        elev = np.linspace(80.0, 220.0, 8)
        ax = plot_topo_section(
            (x_c, z_c, rho), elevation=elev, chainage=x_c, kind="imshow"
        )
        strip_ax = self._find_strip_axes(ax)
        r, g, b, a = strip_ax.get_facecolor()
        assert a == pytest.approx(0.0)

    def test_station_marker_override_applied_to_pcolormesh(self):
        from pycsamt.api.station import StationMarkerStyle

        x_c, z_c, rho = _grid(n_x=8, n_z=15)
        elev = np.linspace(80.0, 220.0, 8)
        custom = StationMarkerStyle(facecolor="#ff0000", edgecolor="#00ff00")
        ax = plot_topo_section(
            (x_c, z_c, rho), elevation=elev, chainage=x_c, station_marker=custom
        )
        scatter = ax.collections[-1]
        assert tuple(scatter.get_facecolor()[0]) == pytest.approx((1.0, 0.0, 0.0, 1.0))
        assert tuple(scatter.get_edgecolor()[0]) == pytest.approx((0.0, 1.0, 0.0, 1.0))

    def test_topo_cfg_override_disables_default_fill_alpha_zero(self):
        from pycsamt.topo import TopoConfig

        x_c, z_c, rho = _grid(n_x=8, n_z=15)
        elev = np.linspace(80.0, 220.0, 8)
        ax = plot_topo_section(
            (x_c, z_c, rho), elevation=elev, chainage=x_c,
            topo_cfg=TopoConfig(fill_alpha=0.4),
        )
        alphas = [p.get_alpha() for p in ax.patches]
        assert alphas and all(a == pytest.approx(0.4) for a in alphas)


# ---------------------------------------------------------------------------
# End-to-end with real WILLY EDI data
# ---------------------------------------------------------------------------


class TestPlotTopoSectionRealData:
    @pytest.mark.parametrize(
        "profile", ["L18PLT", "L22PLT"]
    )
    def test_pcolormesh_all_profiles(self, profile):
        edis = _load(profile)
        n = len(edis)
        x_c, z_c, rho = _grid(n_x=n, n_z=15)
        ax, data = plot_topo_section(
            (x_c, z_c, rho), sites=edis, return_data=True
        )
        assert data.topo_source == "sites"
        assert len(data.station_names) == n
        ax.figure.canvas.draw()

    def test_imshow_with_sites(self):
        edis = _load("L18PLT")
        n = len(edis)
        x_c, z_c, rho = _grid(n_x=n, n_z=15)
        ax = plot_topo_section((x_c, z_c, rho), sites=edis, kind="imshow")
        ax.figure.canvas.draw()

    def test_resistivity_model_with_sites_and_depth_crop(self):
        edis = _load("L18PLT")
        n = len(edis)
        x_c, z_c, rho = _grid(n_x=n, n_z=25, z_km=3.0)
        rm = ResistivityModel.from_array(
            rho, x_c * 1000.0, z_c * 1000.0, method="occam2d", rms=1.05
        )
        ax, data = plot_topo_section(
            rm, sites=edis, depth_max=1000.0, return_data=True
        )
        assert data.method == "occam2d"
        assert data.rms == pytest.approx(1.05)
        assert data.z_centers_km.max() <= 1.0 + 1e-6
        ax.figure.canvas.draw()

    def test_edi_collection_equivalent_to_list(self):
        edis = _load("L18PLT")
        coll = EDICollection(edis)
        n = len(edis)
        x_c, z_c, rho = _grid(n_x=n, n_z=10)
        sec_list = build_topo_section((x_c, z_c, rho), sites=edis)
        sec_coll = build_topo_section((x_c, z_c, rho), sites=coll)
        np.testing.assert_allclose(sec_list.elev_km, sec_coll.elev_km)
        np.testing.assert_allclose(sec_list.chainage_km, sec_coll.chainage_km)

    def test_ai_result_with_sites_topography(self):
        edis = _load("L18PLT")
        n = len(edis)
        n_layers = 12
        rng = np.random.default_rng(3)
        payload = {
            "pred_rho": rng.random((n, n_layers)) * 3.0,
            "depths_km": np.linspace(0.05, 1.5, n_layers),
            "station_names": [f"S{i:03d}" for i in range(n)],
        }
        ax, data = plot_topo_section(
            payload, sites=edis, model_unit="km", return_data=True
        )
        assert data.method == "ai"
        assert data.topo_source == "sites"
        ax.figure.canvas.draw()
