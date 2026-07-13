# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""End-to-end integration tests for the full topo pipeline.

These tests exercise the complete chain:
    EDI files → extract → drape → overlay → rendered figure

All five WILLY profiles (L18–L34, 128 stations) are tested.
"""

from __future__ import annotations

import contextlib
import glob
import os

import numpy as np
import pytest

matplotlib = pytest.importorskip("matplotlib")
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from pycsamt.seg.collection import EDICollection
from pycsamt.seg.edi import EDIFile
from pycsamt.topo.config import (
    PYCSAMT_TOPO,
    TopoConfig,
    configure_topo,
    reset_topo,
)
from pycsamt.topo.drape import (
    drape_section,
    interp_elev,
    station_surface_z,
)
from pycsamt.topo.extract import (
    extract_chainage,
    extract_elevation,
    extract_station_names,
)
from pycsamt.topo.overlay import (
    draw_topo_section,
    draw_topo_strip,
)

# ---------------------------------------------------------------------------
# Data helpers
# ---------------------------------------------------------------------------

_WILLY_ROOT = os.path.join(
    os.path.dirname(__file__),
    "..",
    "..",
    "..",
    "data",
    "AMT",
    "WILLY_DATA",
)

_PROFILES = ["L18PLT", "L22PLT", "L26PLT", "L30PLT", "L34PLT"]


def _load(profile: str) -> list[EDIFile]:
    edi_dir = os.path.join(_WILLY_ROOT, profile)
    paths = sorted(glob.glob(os.path.join(edi_dir, "*.edi")))
    if not paths:
        pytest.skip(f"WILLY {profile} data not found")
    return [EDIFile(p) for p in paths]


def _build_mock_inversion_grid(chain_km, n_depth=20):
    """Return (x_nodes, z_nodes, rho_2d) compatible with drape_section.

    Station positions are used directly as pcolormesh x-node boundaries.
    This gives nx = n_stations - 1 cells, and elev_at_centres should be
    the elevation interpolated to the midpoints between adjacent stations.
    """
    nx = len(chain_km) - 1  # number of cells
    x_nodes = np.array(chain_km, dtype=float)  # (nx+1,) ✓
    z_nodes = np.linspace(0.0, 2.0, n_depth + 1)  # (nz+1,)
    rho_2d = np.random.default_rng(7).random((n_depth, nx)) * 3.0  # (nz, nx)
    return x_nodes, z_nodes, rho_2d


@pytest.fixture(autouse=True)
def _reset_after():
    yield
    reset_topo()
    plt.close("all")


# ---------------------------------------------------------------------------
# Pipeline: extract → drape
# ---------------------------------------------------------------------------


class TestExtractToDrape:
    @pytest.mark.parametrize("profile", _PROFILES)
    def test_extract_then_drape(self, profile):
        edis = _load(profile)
        chain = extract_chainage(edis)
        elev = extract_elevation(edis) / 1000.0  # m → km

        x_nodes, z_nodes, rho = _build_mock_inversion_grid(chain)
        nx = len(chain) - 1
        elev_centres = interp_elev(
            chain, elev, (chain[:-1] + chain[1:]) / 2.0
        )

        xo, zd, do = drape_section(x_nodes, z_nodes, rho, elev_centres)

        assert xo.shape == (nx + 1,)
        assert zd.shape == (len(z_nodes), nx + 1)
        assert do.shape == rho.shape
        assert np.all(np.isfinite(zd))

    def test_l18_chainage_and_drape_shape(self):
        edis = _load("L18PLT")
        chain = extract_chainage(edis)
        elev = extract_elevation(edis) / 1000.0
        assert chain[0] == pytest.approx(0.0)
        assert chain[-1] > 1.0

        x_nodes, z_nodes, rho = _build_mock_inversion_grid(chain, n_depth=10)
        el_c = interp_elev(chain, elev, (chain[:-1] + chain[1:]) / 2.0)
        _, zd, _ = drape_section(x_nodes, z_nodes, rho, el_c)

        # top row of zd should span the elevation range
        assert zd[0].min() >= elev.min() - 0.01
        assert zd[0].max() <= elev.max() + 0.01

    def test_station_pins_z_at_surface(self):
        """Station markers must lie exactly on the terrain surface."""
        edis = _load("L18PLT")
        chain = extract_chainage(edis)
        elev = extract_elevation(edis) / 1000.0
        sz = station_surface_z(chain, elev, chain, exaggeration=1.0)
        np.testing.assert_allclose(sz, elev, rtol=1e-9)

    def test_exaggeration_doubles_relief(self):
        """Exaggeration=2 doubles the z-range of the draped grid.

        The top row z_draped[0] = elev - z_nodes[0]*exag = elev (z_nodes[0]=0
        always), so the surface is unchanged.  The depth extent of the grid
        is what doubles: range(zd) = z_nodes[-1] * exag.
        """
        edis = _load("L18PLT")
        chain = extract_chainage(edis)
        elev = extract_elevation(edis) / 1000.0

        x_nodes, z_nodes, rho = _build_mock_inversion_grid(chain, n_depth=5)
        el_c = interp_elev(chain, elev, (chain[:-1] + chain[1:]) / 2.0)

        _, zd1, _ = drape_section(
            x_nodes, z_nodes, rho, el_c, exaggeration=1.0
        )
        _, zd2, _ = drape_section(
            x_nodes, z_nodes, rho, el_c, exaggeration=2.0
        )

        # Surface row (k=0) is identical — exaggeration doesn't move the surface
        np.testing.assert_allclose(zd1[0], zd2[0], rtol=1e-9)

        # Bottom row z range should double: zd1_range = z_nodes[-1], zd2_range = 2×
        range1 = (zd1[0] - zd1[-1]).mean()
        range2 = (zd2[0] - zd2[-1]).mean()
        assert range2 / range1 == pytest.approx(2.0, rel=0.01)


# ---------------------------------------------------------------------------
# Pipeline: extract → overlay → figure
# ---------------------------------------------------------------------------


class TestExtractToOverlay:
    @pytest.mark.parametrize("profile", _PROFILES)
    def test_draw_section_all_profiles(self, profile):
        edis = _load(profile)
        chain = extract_chainage(edis)
        elev = extract_elevation(edis)
        names = extract_station_names(edis)

        fig, ax = plt.subplots()
        ax.set_xlim(chain[0], chain[-1])
        ax.set_ylim(0.02, 0.25)  # km scale
        draw_topo_section(ax, chain, elev, names)
        fig.canvas.draw()

    @pytest.mark.parametrize("profile", _PROFILES)
    def test_draw_strip_all_profiles(self, profile):
        edis = _load(profile)
        chain = extract_chainage(edis)
        elev = extract_elevation(edis)
        names = extract_station_names(edis)

        fig, ax = plt.subplots()
        ax.imshow(np.ones((10, len(edis))), aspect="auto")
        ax_strip = draw_topo_strip(fig, ax, chain, elev, names)
        assert ax_strip is not None
        fig.canvas.draw()

    def test_full_pipeline_l18_section_plot(self):
        edis = _load("L18PLT")
        chain = extract_chainage(edis)
        elev = extract_elevation(edis)
        names = extract_station_names(edis)
        elev_km = elev / 1000.0

        # Build a synthetic drape-ready section
        x_nodes, z_nodes, rho = _build_mock_inversion_grid(chain, n_depth=20)
        el_c = interp_elev(chain, elev_km, (chain[:-1] + chain[1:]) / 2.0)
        _, zd, rho_out = drape_section(x_nodes, z_nodes, rho, el_c)

        fig, ax = plt.subplots()
        ax.pcolormesh(
            x_nodes, zd, rho_out, shading="flat", cmap="jet", vmin=0, vmax=3
        )
        draw_topo_section(ax, chain, elev, names)
        ax.set_xlabel("Chainage (km)")
        ax.set_ylabel("Elevation (km)")
        fig.canvas.draw()  # must not raise


# ---------------------------------------------------------------------------
# configure_topo global singleton integration
# ---------------------------------------------------------------------------


class TestConfigureTopo:
    def test_configure_topo_affects_enabled_check(self):
        configure_topo(enabled=True, exaggeration=1.5)
        assert PYCSAMT_TOPO.enabled is True
        assert PYCSAMT_TOPO.exaggeration == pytest.approx(1.5)

    def test_context_manager_isolates_changes(self):
        configure_topo(enabled=True)
        with PYCSAMT_TOPO.context(exaggeration=5.0):
            assert PYCSAMT_TOPO.exaggeration == pytest.approx(5.0)
        assert PYCSAMT_TOPO.exaggeration == pytest.approx(1.0)
        assert PYCSAMT_TOPO.enabled is True  # unchanged by the context

    def test_reset_restores_disabled(self):
        configure_topo(enabled=True, exaggeration=3.0)
        reset_topo()
        assert PYCSAMT_TOPO.enabled is False
        assert PYCSAMT_TOPO.exaggeration == pytest.approx(1.0)

    def test_context_manager_survives_exception(self):
        configure_topo(exaggeration=2.0)
        with contextlib.suppress(ValueError):
            with PYCSAMT_TOPO.context(exaggeration=9.0):
                raise ValueError("test error")
        assert PYCSAMT_TOPO.exaggeration == pytest.approx(2.0)

    def test_cfg_passed_to_draw_section(self):
        edis = _load("L18PLT")
        chain = extract_chainage(edis)
        elev = extract_elevation(edis)

        cfg = TopoConfig()
        cfg.configure(exaggeration=3.0, fill_color="#0000ff", fill_alpha=0.5)

        fig, ax = plt.subplots()
        ax.set_xlim(chain[0], chain[-1])
        ax.set_ylim(0.02, 0.25)
        draw_topo_section(ax, chain, elev, cfg=cfg)
        fig.canvas.draw()


# ---------------------------------------------------------------------------
# All-profiles combined: consistency checks
# ---------------------------------------------------------------------------


class TestAllProfilesConsistency:
    def test_chainage_starts_at_zero_all_profiles(self):
        for p in _PROFILES:
            chain = extract_chainage(_load(p))
            assert chain[0] == pytest.approx(0.0), (
                f"{p} chain does not start at 0"
            )

    def test_elevation_in_expected_range_all_profiles(self):
        """WILLY elevations are between 20 m and 300 m."""
        for p in _PROFILES:
            elev = extract_elevation(_load(p))
            assert elev.min() > 20.0, f"{p}: min elev {elev.min()} < 20 m"
            assert elev.max() < 300.0, f"{p}: max elev {elev.max()} > 300 m"

    def test_station_count_all_profiles(self):
        expected = {
            "L18PLT": 28,
            "L22PLT": 25,
            "L26PLT": 25,
            "L30PLT": 25,
            "L34PLT": 25,
        }
        for p, n in expected.items():
            edis = _load(p)
            assert len(edis) == n, (
                f"{p}: expected {n} stations, got {len(edis)}"
            )

    def test_l34_has_highest_terrain(self):
        """L34PLT is documented as having the highest terrain (59–224 m)."""
        max_elevs = {}
        for p in _PROFILES:
            max_elevs[p] = extract_elevation(_load(p)).max()
        assert max_elevs["L34PLT"] == max(max_elevs.values())

    def test_edi_collection_same_as_list(self):
        for p in _PROFILES:
            edis = _load(p)
            coll = EDICollection(edis)
            np.testing.assert_array_equal(
                extract_elevation(edis),
                extract_elevation(coll),
                err_msg=f"{p}: list vs collection elevation mismatch",
            )
            np.testing.assert_allclose(
                extract_chainage(edis),
                extract_chainage(coll),
                rtol=1e-9,
                err_msg=f"{p}: list vs collection chainage mismatch",
            )

    def test_interp_elev_between_profiles(self):
        """Elevation should be positive and finite when interpolated to a fine grid."""
        for p in _PROFILES:
            edis = _load(p)
            chain = extract_chainage(edis)
            elev = extract_elevation(edis) / 1000.0
            fine = np.linspace(chain[0], chain[-1], 200)
            out = interp_elev(chain, elev, fine)
            assert np.all(np.isfinite(out)), (
                f"{p}: NaN in interpolated elevation"
            )
            assert out.min() > 0.0, (
                f"{p}: zero/negative interpolated elevation"
            )

    @pytest.mark.parametrize("profile", _PROFILES)
    def test_drape_zd_top_row_within_elev_bounds(self, profile):
        """Top row of z_draped grid must match the terrain elevation."""
        edis = _load(profile)
        chain = extract_chainage(edis)
        elev = extract_elevation(edis) / 1000.0

        x_nodes, z_nodes, rho = _build_mock_inversion_grid(chain, n_depth=10)
        el_c = interp_elev(chain, elev, (chain[:-1] + chain[1:]) / 2.0)
        _, zd, _ = drape_section(
            x_nodes, z_nodes, rho, el_c, exaggeration=1.0
        )

        # Inner columns (not boundary extrapolation) stay within elev bounds
        inner = zd[0, 1:-1]
        assert inner.min() >= elev.min() - 0.005
        assert inner.max() <= elev.max() + 0.005


# ---------------------------------------------------------------------------
# Robustness: topo pipeline never raises for degenerate inputs
# ---------------------------------------------------------------------------


class TestRobustness:
    def test_pipeline_with_two_stations(self):
        """Degenerate 2-station profile must not crash."""
        edis = _load("L18PLT")[:2]
        chain = extract_chainage(edis)
        elev = extract_elevation(edis)
        names = extract_station_names(edis)

        fig, ax = plt.subplots()
        ax.set_xlim(chain[0], chain[-1] + 0.1)
        ax.set_ylim(0.05, 0.2)
        draw_topo_section(ax, chain, elev, names)
        fig.canvas.draw()

    def test_pipeline_constant_elevation(self):
        """Perfectly flat terrain should not produce NaN in the drape grid."""
        chain = np.linspace(0.0, 3.0, 10)
        x_nodes, z_nodes, rho = _build_mock_inversion_grid(chain, n_depth=8)
        nx = len(chain) - 1
        el_c = np.full(nx, 0.15)  # flat terrain, one value per cell
        _, zd, _ = drape_section(x_nodes, z_nodes, rho, el_c)
        assert np.all(np.isfinite(zd))

    def test_draw_section_does_not_modify_axes_xlim(self):
        chain, elev, _ = (
            extract_chainage(_load("L18PLT")),
            extract_elevation(_load("L18PLT")),
            extract_station_names(_load("L18PLT")),
        )
        fig, ax = plt.subplots()
        ax.set_xlim(0.0, 5.0)
        xlim_before = ax.get_xlim()
        draw_topo_section(ax, chain, elev)
        assert ax.get_xlim() == pytest.approx(xlim_before)

    def test_context_manager_thread_safe_sequential(self):
        """Two sequential context uses must each restore state."""
        cfg = TopoConfig()
        for exag in (2.0, 4.0, 8.0):
            with cfg.context(exaggeration=exag):
                assert cfg.exaggeration == pytest.approx(exag)
            assert cfg.exaggeration == pytest.approx(1.0)
