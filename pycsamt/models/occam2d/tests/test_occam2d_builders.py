"""Tests for Phase-6 builder methods:
OccamData.from_edi, OccamMesh.from_data,
OccamModel.from_mesh, OccamStartup.from_model, InputBuilder.build
"""

import math

import numpy as np
import pytest

# ---------------------------------------------------------------------------
# Synthetic site fixture (no pyproj / full EDI stack required)
# ---------------------------------------------------------------------------


class _FakeSite:
    """Minimal duck-typed site object accepted by OccamData.from_edi."""

    def __init__(
        self,
        name,
        lat,
        lon,
        freqs,
        rho_xy,
        rho_yx,
        phs_xy,
        phs_yx,
        rho_xy_err=None,
        rho_yx_err=None,
        phs_xy_err=None,
        phs_yx_err=None,
    ):
        self.name = name
        self.coords = (lat, lon, 0.0)
        n = len(freqs)
        self.freq = np.asarray(freqs, dtype=float)
        # Build full 2×2 rho/phase arrays (only xy/yx used)
        rho = np.zeros((n, 2, 2), dtype=float)
        rho[:, 0, 1] = rho_xy
        rho[:, 1, 0] = rho_yx
        self.rho = rho

        phs = np.zeros((n, 2, 2), dtype=float)
        phs[:, 0, 1] = phs_xy
        phs[:, 1, 0] = phs_yx
        self.phase = phs

        if rho_xy_err is not None:
            re = np.zeros((n, 2, 2), dtype=float)
            re[:, 0, 1] = rho_xy_err
            re[:, 1, 0] = rho_yx_err or rho_xy_err
            self.rho_err = re
        else:
            self.rho_err = None

        if phs_xy_err is not None:
            pe = np.zeros((n, 2, 2), dtype=float)
            pe[:, 0, 1] = phs_xy_err
            pe[:, 1, 0] = phs_yx_err or phs_xy_err
            self.phase_err = pe
        else:
            self.phase_err = None


def _make_sites(n_sites=3, n_freqs=10):
    freqs = np.logspace(4, -1, n_freqs)  # 10 kHz → 0.1 Hz
    sites = []
    for i in range(n_sites):
        lon = 0.0 + i * 0.01  # spaced ~1 km east
        rho = np.full(n_freqs, 100.0)  # 100 Ω·m half-space
        phs = np.full(n_freqs, 45.0)  # 45° half-space (TE)
        phs_yx = np.full(n_freqs, -135.0)  # ≡ 45° after +180° normalisation
        sites.append(
            _FakeSite(
                name=f"S{i:02d}",
                lat=0.0,
                lon=lon,
                freqs=freqs,
                rho_xy=rho,
                rho_yx=rho,
                phs_xy=phs,
                phs_yx=phs_yx,
            )
        )
    return sites


# ---------------------------------------------------------------------------
# OccamData.from_edi
# ---------------------------------------------------------------------------


def test_from_edi_basic():
    from pycsamt.models.occam2d.data import OccamData

    sites = _make_sites(3, 8)
    d = OccamData.from_edi(sites, modes=["TE", "TM"])
    assert d.n_sites == 3
    assert d.n_frequencies == 8
    assert d.n_data > 0


def test_from_edi_site_names():
    from pycsamt.models.occam2d.data import OccamData

    sites = _make_sites(3)
    d = OccamData.from_edi(sites)
    assert d.sites[0] == "S00"


def test_from_edi_offsets_monotone():
    from pycsamt.models.occam2d.data import OccamData

    sites = _make_sites(4)
    d = OccamData.from_edi(sites)
    assert np.all(np.diff(d.offsets) >= 0), "offsets not monotone"


def test_from_edi_frequencies_descending():
    from pycsamt.models.occam2d.data import OccamData

    sites = _make_sites(2, 6)
    d = OccamData.from_edi(sites)
    assert np.all(np.diff(d.frequencies) < 0), "frequencies not descending"


def test_from_edi_type_codes():
    from pycsamt.models.occam2d.data import OccamData

    sites = _make_sites(2, 5)
    d = OccamData.from_edi(sites, modes=["TE", "TM"])
    codes = set(d.data_blocks[:, 2].astype(int))
    assert 1 in codes  # RhoTE
    assert 2 in codes  # PhsTE
    assert 5 in codes  # RhoTM
    assert 6 in codes  # PhsTM


def test_from_edi_te_only():
    from pycsamt.models.occam2d.data import OccamData

    sites = _make_sites(2, 5)
    d = OccamData.from_edi(sites, modes=["TE"])
    codes = set(d.data_blocks[:, 2].astype(int))
    assert 1 in codes and 2 in codes
    assert 5 not in codes and 6 not in codes


def test_from_edi_phs_tm_normalised():
    """PhsTM should be phase_yx + 180 = -135 + 180 = 45 degrees."""
    from pycsamt.models.occam2d.data import OccamData

    sites = _make_sites(2, 5)
    d = OccamData.from_edi(sites, modes=["TM"])
    phs = d.data_blocks[d.data_blocks[:, 2] == 6, 3]
    assert np.all(np.abs(phs - 45.0) < 0.1), (
        f"PhsTM not normalised: {phs[:3]}"
    )


def test_from_edi_rho_te_log10():
    """RhoTE datum should be log10(100) = 2.0 for half-space data."""
    from pycsamt.models.occam2d.data import OccamData

    sites = _make_sites(2, 5)
    d = OccamData.from_edi(sites, modes=["TE"])
    rho = d.data_blocks[d.data_blocks[:, 2] == 1, 3]
    assert np.all(np.abs(rho - 2.0) < 0.01), f"log10(rho) != 2.0: {rho[:3]}"


def test_from_edi_error_floor():
    from pycsamt.models.occam2d.config import OccamConfig
    from pycsamt.models.occam2d.data import OccamData

    sites = _make_sites(2, 4)
    cfg = OccamConfig(error_floor_rho=0.05, error_floor_phase=0.5)
    d = OccamData.from_edi(sites, config=cfg)
    rho_errs = d.data_blocks[
        np.isin(d.data_blocks[:, 2].astype(int), [1, 5]), 4
    ]
    assert np.all(rho_errs > 0)
    phs_errs = d.data_blocks[
        np.isin(d.data_blocks[:, 2].astype(int), [2, 6]), 4
    ]
    assert np.all(phs_errs >= cfg.error_floor_phase - 1e-9)


def test_from_edi_freq_band():
    from pycsamt.models.occam2d.config import OccamConfig
    from pycsamt.models.occam2d.data import OccamData

    sites = _make_sites(2, 10)
    cfg = OccamConfig(freq_min=1.0, freq_max=100.0)
    d = OccamData.from_edi(sites, config=cfg)
    assert d.frequencies.min() >= 1.0
    assert d.frequencies.max() <= 100.0


def test_from_edi_empty_source_raises():
    from pycsamt.models.occam2d.data import OccamData

    with pytest.raises(ValueError, match="no sites"):
        OccamData.from_edi([])


def test_from_edi_roundtrip(tmp_path):
    """from_edi → write → read should reproduce the data."""
    from pycsamt.models.occam2d.data import OccamData

    sites = _make_sites(3, 6)
    d = OccamData.from_edi(sites)
    out = tmp_path / "OccamDataFile.dat"
    d.write(out)
    back = OccamData.read(out)
    assert back.n_sites == d.n_sites
    assert back.n_frequencies == d.n_frequencies
    assert np.allclose(back.offsets, d.offsets, rtol=1e-4)
    assert np.allclose(back.frequencies, d.frequencies, rtol=1e-6)


# ---------------------------------------------------------------------------
# OccamMesh.from_data
# ---------------------------------------------------------------------------


def _simple_data(n=3):
    from pycsamt.models.occam2d.data import OccamData

    d = OccamData()
    d.sites = [f"S{i:02d}" for i in range(n)]
    d.offsets = np.linspace(0, (n - 1) * 500.0, n)
    d.frequencies = np.array([100.0, 10.0, 1.0])
    return d


def test_mesh_from_data_shape():
    from pycsamt.models.occam2d.config import OccamConfig
    from pycsamt.models.occam2d.mesh import OccamMesh

    d = _simple_data(5)
    cfg = OccamConfig(n_layers=10, n_airlayers=3, cell_size_horizontal=100.0)
    m = OccamMesh.from_data(d, config=cfg)
    assert m.n_xcells > 0
    assert m.n_zcells == cfg.n_layers + cfg.n_airlayers


def test_mesh_from_data_xcells_even():
    from pycsamt.models.occam2d.mesh import OccamMesh

    d = _simple_data(4)
    m = OccamMesh.from_data(d)
    assert m.n_xcells % 2 == 0, f"n_xcells={m.n_xcells} is not even"


def test_mesh_from_data_xwidths_positive():
    from pycsamt.models.occam2d.mesh import OccamMesh

    d = _simple_data(3)
    m = OccamMesh.from_data(d)
    assert np.all(m.x_widths > 0)


def test_mesh_from_data_zwidths_positive():
    from pycsamt.models.occam2d.mesh import OccamMesh

    d = _simple_data(3)
    m = OccamMesh.from_data(d)
    assert np.all(m.z_widths > 0)


def test_mesh_from_data_cell_rows_count():
    from pycsamt.models.occam2d.config import OccamConfig
    from pycsamt.models.occam2d.mesh import OccamMesh

    d = _simple_data(3)
    cfg = OccamConfig(n_layers=5, n_airlayers=2)
    m = OccamMesh.from_data(d, config=cfg)
    # 4 rows per z-cell
    assert len(m.cell_rows) == 4 * m.n_zcells


def test_mesh_from_data_air_rows_fixed():
    """Air rows should use '0'; active rows use '?'."""
    from pycsamt.models.occam2d.config import OccamConfig
    from pycsamt.models.occam2d.mesh import OccamMesh

    d = _simple_data(3)
    cfg = OccamConfig(n_layers=5, n_airlayers=2)
    m = OccamMesh.from_data(d, config=cfg)
    for row in m.cell_rows[:8]:  # first 2 z-cells × 4 rows = air
        assert "?" not in row
    for row in m.cell_rows[8:]:  # rest are active
        assert "0" not in row


def test_mesh_from_data_roundtrip(tmp_path):
    from pycsamt.models.occam2d.mesh import OccamMesh

    d = _simple_data(4)
    m = OccamMesh.from_data(d)
    out = tmp_path / "Occam2DMesh"
    m.write(out)
    back = OccamMesh.read(out)
    assert back.n_xcells == m.n_xcells
    assert back.n_zcells == m.n_zcells
    assert np.allclose(back.x_widths, m.x_widths, rtol=1e-4)
    assert np.allclose(back.z_widths, m.z_widths, rtol=1e-4)


def test_mesh_from_data_single_station():
    from pycsamt.models.occam2d.mesh import OccamMesh

    d = _simple_data(1)
    m = OccamMesh.from_data(d)
    assert m.n_xcells > 0


def test_mesh_from_data_no_stations_raises():
    from pycsamt.models.occam2d.data import OccamData
    from pycsamt.models.occam2d.mesh import OccamMesh

    d = OccamData()
    d.offsets = np.array([])
    with pytest.raises(ValueError, match="no station offsets"):
        OccamMesh.from_data(d)


# ---------------------------------------------------------------------------
# OccamModel.from_mesh
# ---------------------------------------------------------------------------


def _simple_mesh(n_xcells=50, n_zcells=8, n_air=2):
    from pycsamt.models.occam2d.mesh import OccamMesh

    # Build a mesh with exactly n_xcells x-cells (even, >= 14)
    # achieved by n_pad=7 each side + interior
    n_xcells - 14
    m = OccamMesh()
    m.x_widths = np.ones(n_xcells) * 100.0
    m.z_widths = np.ones(n_zcells) * 50.0
    m.x_nodes = np.concatenate([[0.0], np.cumsum(m.x_widths)])
    m.z_nodes = np.concatenate([[0.0], np.cumsum(m.z_widths)])
    m.n_airlayers = n_air
    m.cell_rows = ["?" * n_xcells] * (4 * n_zcells)
    return m


def test_model_from_mesh_n_layers():
    from pycsamt.models.occam2d.model import OccamModel

    m = _simple_mesh(50, 10, 2)
    mod = OccamModel.from_mesh(m)
    assert mod.n_layers == 8  # 10 - 2 air


def test_model_from_mesh_code_sum():
    """Sum of codes per layer must equal n_xcells."""
    from pycsamt.models.occam2d.model import OccamModel

    m = _simple_mesh(50, 8, 2)
    mod = OccamModel.from_mesh(m)
    for layer in mod.layers:
        assert sum(layer["params"]) == m.n_xcells


def test_model_from_mesh_odd_interior_width_absorbs_remainder():
    from pycsamt.models.occam2d.model import OccamModel

    m = _simple_mesh(17, 8, 2)
    mod = OccamModel.from_mesh(m)
    for layer in mod.layers:
        assert sum(layer["params"]) == m.n_xcells
        assert layer["params"].tolist() == [7, 3, 7]


def test_model_from_mesh_boundary_codes():
    """First and last code in each layer must be 7."""
    from pycsamt.models.occam2d.model import OccamModel

    m = _simple_mesh(50, 8, 2)
    mod = OccamModel.from_mesh(m)
    for layer in mod.layers:
        assert layer["params"][0] == 7
        assert layer["params"][-1] == 7


def test_model_from_mesh_interior_codes_even():
    """All interior codes (not first/last) must be even."""
    from pycsamt.models.occam2d.model import OccamModel

    m = _simple_mesh(50, 8, 2)
    mod = OccamModel.from_mesh(m)
    for layer in mod.layers:
        for code in layer["params"][1:-1]:
            assert code % 2 == 0, f"odd interior code: {code}"


def test_model_from_mesh_n_params_positive():
    from pycsamt.models.occam2d.model import OccamModel

    m = _simple_mesh(50, 8, 2)
    mod = OccamModel.from_mesh(m)
    assert mod.n_params > 0


def test_model_from_mesh_roundtrip(tmp_path):
    from pycsamt.models.occam2d.model import OccamModel

    m = _simple_mesh(50, 8, 2)
    mod = OccamModel.from_mesh(m)
    out = tmp_path / "Occam2DModel"
    mod.write(out)
    back = OccamModel.read(out)
    assert back.n_layers == mod.n_layers
    assert back.n_params == mod.n_params
    for bl, ol in zip(back.layers, mod.layers):
        assert np.array_equal(bl["params"], ol["params"])


def test_model_from_mesh_too_narrow_raises():
    from pycsamt.models.occam2d.model import OccamModel

    m = _simple_mesh(10, 5, 1)  # 10 < 14
    with pytest.raises(ValueError, match="too narrow"):
        OccamModel.from_mesh(m)


def test_model_from_mesh_no_active_raises():
    from pycsamt.models.occam2d.model import OccamModel

    m = _simple_mesh(50, 3, 3)  # all air
    with pytest.raises(ValueError, match="no active layers"):
        OccamModel.from_mesh(m)


# ---------------------------------------------------------------------------
# OccamStartup.from_model
# ---------------------------------------------------------------------------


def _simple_model(n_xcells=50, n_layers=5):
    from pycsamt.models.occam2d.model import OccamModel

    m = _simple_mesh(n_xcells, n_layers + 2, 2)
    return OccamModel.from_mesh(m)


def test_startup_from_model_n_params():
    from pycsamt.models.occam2d.startup import OccamStartup

    mod = _simple_model()
    su = OccamStartup.from_model(mod)
    assert su.n_params == mod.n_params


def test_startup_from_model_param_values():
    from pycsamt.models.occam2d.config import OccamConfig
    from pycsamt.models.occam2d.startup import OccamStartup

    mod = _simple_model()
    cfg = OccamConfig(initial_rho=100.0)
    su = OccamStartup.from_model(mod, config=cfg)
    expected = math.log10(100.0)
    assert np.allclose(su.param_values, expected, atol=1e-6)


def test_startup_from_model_iteration_zero():
    from pycsamt.models.occam2d.startup import OccamStartup

    mod = _simple_model()
    su = OccamStartup.from_model(mod)
    assert su.iteration == 0


def test_startup_from_model_param_count():
    from pycsamt.models.occam2d.startup import OccamStartup

    mod = _simple_model()
    su = OccamStartup.from_model(mod)
    assert len(su.param_values) == su.n_params


def test_startup_from_model_roundtrip(tmp_path):
    from pycsamt.models.occam2d.startup import OccamStartup

    mod = _simple_model()
    su = OccamStartup.from_model(mod)
    out = tmp_path / "Startup"
    su.write(out)
    back = OccamStartup.read(out)
    assert back.n_params == su.n_params
    assert back.iteration == 0
    assert np.allclose(back.param_values, su.param_values, atol=1e-4)


def test_startup_from_model_no_params_raises():
    from pycsamt.models.occam2d.model import OccamModel
    from pycsamt.models.occam2d.startup import OccamStartup

    mod = OccamModel()  # empty model, n_params = 0
    with pytest.raises(ValueError, match="no parameters"):
        OccamStartup.from_model(mod)


# ---------------------------------------------------------------------------
# InputBuilder (integration test — uses synthetic OccamData, skips from_edi)
# ---------------------------------------------------------------------------


def test_input_builder_integration(tmp_path):
    """Full pipeline: data → mesh → model → startup → files on disk."""
    from pycsamt.models.occam2d.config import OccamConfig
    from pycsamt.models.occam2d.data import OccamData
    from pycsamt.models.occam2d.mesh import OccamMesh
    from pycsamt.models.occam2d.model import OccamModel
    from pycsamt.models.occam2d.startup import OccamStartup
    from pycsamt.models.occam2d.validation import (
        is_data_file,
        is_mesh_file,
        is_model_file,
        is_startup_file,
    )

    # Build the four objects manually (bypassing from_edi)
    cfg = OccamConfig(
        n_layers=5, n_airlayers=2, cell_size_horizontal=200.0, depth_scale=1.3
    )

    data = OccamData.from_edi(_make_sites(3, 4), config=cfg)
    mesh = OccamMesh.from_data(data, config=cfg)
    mod = OccamModel.from_mesh(mesh, config=cfg)
    su = OccamStartup.from_model(mod, config=cfg)

    # Write to tmp dir
    data.write(tmp_path / cfg.data_file)
    mesh.write(tmp_path / cfg.mesh_file)
    mod.write(tmp_path / cfg.model_file)
    su.write(tmp_path / cfg.startup_file)

    # All files exist
    for fname in [
        cfg.data_file,
        cfg.mesh_file,
        cfg.model_file,
        cfg.startup_file,
    ]:
        assert (tmp_path / fname).exists(), f"missing {fname}"

    # Validation predicates recognise each file
    assert is_data_file(tmp_path / cfg.data_file)
    assert is_mesh_file(tmp_path / cfg.mesh_file)
    assert is_model_file(tmp_path / cfg.model_file)
    assert is_startup_file(tmp_path / cfg.startup_file)

    # Startup n_params == model n_params
    assert su.n_params == mod.n_params

    # Mesh and model are self-consistent
    assert sum(mod.layers[0]["params"]) == mesh.n_xcells


def test_input_builder_build_writes_all_required_files(tmp_path):
    from pycsamt.models.occam2d.builder import InputBuilder
    from pycsamt.models.occam2d.config import OccamConfig

    cfg = OccamConfig(
        n_layers=5, n_airlayers=2, cell_size_horizontal=200.0, depth_scale=1.3
    )
    builder = InputBuilder(_make_sites(3, 4), workdir=tmp_path, config=cfg)

    builder.build()

    for fname in [
        cfg.data_file,
        cfg.mesh_file,
        cfg.model_file,
        cfg.startup_file,
    ]:
        assert (tmp_path / fname).exists(), f"missing {fname}"
    assert builder.is_ready
    assert builder.startup.n_params == builder.model.n_params
    assert builder.data.path == tmp_path / cfg.data_file
    assert builder.mesh.path == tmp_path / cfg.mesh_file
    assert builder.model.path == tmp_path / cfg.model_file
    assert builder.startup.path == tmp_path / cfg.startup_file


def test_input_builder_summary(tmp_path):
    from pycsamt.models.occam2d.builder import InputBuilder
    from pycsamt.models.occam2d.config import OccamConfig
    from pycsamt.models.occam2d.data import OccamData
    from pycsamt.models.occam2d.mesh import OccamMesh
    from pycsamt.models.occam2d.model import OccamModel
    from pycsamt.models.occam2d.startup import OccamStartup

    cfg = OccamConfig(n_layers=4, n_airlayers=1)
    ib = InputBuilder(_make_sites(2, 3), workdir=tmp_path, config=cfg)
    ib.data = OccamData.from_edi(_make_sites(2, 3), config=cfg)
    ib.mesh = OccamMesh.from_data(ib.data, config=cfg)
    ib.model = OccamModel.from_mesh(ib.mesh, config=cfg)
    ib.startup = OccamStartup.from_model(ib.model, config=cfg)

    summary = ib.summary()
    assert "sites" in summary
    assert "mesh" in summary
    assert ib.is_ready
