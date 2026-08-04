"""Tests for the external MARE2DEM adapter (triangular-mesh MT backend).

A compiled ``MARE2DEM`` binary is *not* guaranteed to exist in any given
environment (same situation as ``modem3d.py`` -- see
``test_maxwell_modem3d.py``'s own module docstring), so most of these
tests validate everything that does not require genuine physics from a
real binary: input-file generation (round-tripped through the real
``pycsamt.models.mare2dem`` readers, not a second parser), the full
external-process pipeline against a scripted stand-in executable, the
region-uniformity/boundary-segments preflight checks, and the
missing-executable path.

A separate ``requires_real_mare2dem``-gated section near the bottom runs
a real half-space forward solve against a genuinely compiled MARE2DEM
binary, when one is found via the ``PYCSAMT_MARE2DEM_BINARY``
environment variable. Skipped, not failed, when no binary is present --
confirmed passing on 2026-07-31 against a real binary built from the
vendored source inside WSL2 Ubuntu with Intel's oneAPI toolchain (see
``mare2dem.py``'s own module docstring for the full build/validation
story and every bug that real run surfaced).
"""

from __future__ import annotations

import os
import shutil
import sys
import textwrap
from pathlib import Path

import numpy as np
import pytest

from pycsamt.forward.maxwell import (
    ExecutableNotFoundError,
    ExternalRunPolicy,
    IncompatibleProblemError,
    ReceiverSet,
    TriMesh,
    TriProblem,
    create_backend,
    list_backends,
)
from pycsamt.forward.maxwell.benchmarks import half_space_impedance
from pycsamt.forward.maxwell.mare2dem import (
    _DATA_CODES,
    Mare2DEMAdapter,
    register_mare2dem_backend,
)
from pycsamt.models.mare2dem.config import Mare2DEMConfig
from pycsamt.models.mare2dem.iotools.poly import read_poly


def _mesh(*, with_boundary=True):
    return TriMesh(
        [[0, 0], [1000, 0], [1000, 500], [0, 500]],
        [[0, 1, 2], [0, 2, 3]],
        region_ids=[1, 1],
        boundary_segments=(
            [[0, 1], [1, 2], [2, 3], [3, 0]] if with_boundary else None
        ),
    )


def _problem(*, conductivity=0.01, components=("zxy", "zyx"), mesh=None):
    mesh = mesh if mesh is not None else _mesh()
    return TriProblem(
        mesh,
        np.full(mesh.shape, conductivity),
        [10.0, 1.0],
        ReceiverSet([[200, 0], [800, 0]], ["S00", "S01"]),
        components,
    )


def _adapter():
    return Mare2DEMAdapter(run_policy=ExternalRunPolicy("does-not-exist-xyz"))


# ---------------------------------------------------------------------------
# assess()
# ---------------------------------------------------------------------------


def test_capabilities_declare_triangular_2d_scope():
    adapter = _adapter()
    assert adapter.capabilities.name == "mare2dem"
    assert adapter.capabilities.dimensions == (2,)
    assert adapter.capabilities.components == ("zxy", "zyx")


def test_assess_accepts_a_uniform_per_region_problem():
    report = _adapter().assess(_problem())
    assert report.compatible


def test_assess_rejects_missing_boundary_segments():
    mesh = _mesh(with_boundary=False)
    report = _adapter().assess(_problem(mesh=mesh))
    assert not report.compatible
    assert "boundary_segments" in report.errors[0]


def test_assess_rejects_conductivity_that_varies_within_one_region():
    mesh = _mesh()
    problem = TriProblem(
        mesh, np.array([0.01, 0.02]), [10.0], ReceiverSet([[500, 0]], ["S"])
    )
    report = _adapter().assess(problem)
    assert not report.compatible
    assert "not uniform within" in report.errors[-1]


def test_solve_raises_incompatible_before_writing_any_file(tmp_path):
    mesh = _mesh()
    problem = TriProblem(
        mesh, np.array([0.01, 0.02]), [10.0], ReceiverSet([[500, 0]], ["S"])
    )
    with pytest.raises(IncompatibleProblemError, match="not uniform"):
        _adapter().solve(problem)


# ---------------------------------------------------------------------------
# _prepare_inputs: real file round trips (no subprocess involved)
# ---------------------------------------------------------------------------


def test_prepare_inputs_writes_a_pslg_the_real_reader_round_trips(tmp_path):
    problem = _problem()
    adapter = _adapter()
    context = adapter._prepare_inputs(problem, tmp_path)

    pf = read_poly(tmp_path / "mare2dem.poly")
    np.testing.assert_allclose(pf.nodes, problem.mesh.nodes_m)
    # boundary_segments is 0-based; the .poly file is 1-based.
    np.testing.assert_array_equal(
        pf.segments, np.asarray(problem.mesh.boundary_segments) + 1
    )
    assert pf.n_regions == 1
    assert context["resistivity_path"].exists()
    assert context["data_path"].exists()
    assert context["settings_path"].exists()


def test_prepare_inputs_region_seed_is_interior_to_the_mesh(tmp_path):
    problem = _problem()
    adapter = _adapter()
    adapter._prepare_inputs(problem, tmp_path)

    pf = read_poly(tmp_path / "mare2dem.poly")
    seed_x, seed_z = pf.regions[0, 0], pf.regions[0, 1]
    x0, x1 = problem.mesh.nodes_m[:, 0].min(), problem.mesh.nodes_m[:, 0].max()
    z0, z1 = problem.mesh.nodes_m[:, 1].min(), problem.mesh.nodes_m[:, 1].max()
    assert x0 < seed_x < x1
    assert z0 < seed_z < z1


def test_prepare_inputs_resistivity_is_linear_not_log10(tmp_path):
    """Confirmed against a real compiled binary: the .resistivity file's
    region table is linear ohm-m (see mare2dem.py's module docstring for
    the full story -- writing log10(rho) here instead produced a
    self-consistent but wrong forward response, magnitude off by a
    constant factor at every frequency).
    """
    from pycsamt.models.mare2dem.iotools.resistivity import read_resistivity

    problem = _problem(conductivity=0.1)  # 10 ohm-m
    adapter = _adapter()
    context = adapter._prepare_inputs(problem, tmp_path)

    rf = read_resistivity(context["resistivity_path"])
    assert rf.num_regions == 1
    np.testing.assert_allclose(rf.resistivity[0, 0], 10.0)
    assert rf.max_iterations == 0


def test_prepare_inputs_resistivity_filename_has_iteration_number(tmp_path):
    """MARE2DEM's own CLI parser requires <stem>.<iteration>[.resistivity];
    a bare stem is a hard error ("no iteration number in resistivity
    file"), confirmed against a real compiled binary.
    """
    problem = _problem()
    adapter = _adapter()
    context = adapter._prepare_inputs(problem, tmp_path)

    assert context["run_stem"] == "mare2dem.0"
    assert context["resistivity_path"].name == "mare2dem.0.resistivity"


def test_prepare_inputs_writes_one_data_row_pair_per_freq_rx_component(
    tmp_path,
):
    from pycsamt.models.mare2dem.iotools.emdata import read_emdata

    problem = _problem()
    adapter = _adapter()
    context = adapter._prepare_inputs(problem, tmp_path)

    em = read_emdata(context["data_path"])
    n_freq, n_rx, n_comp = 2, 2, 2
    assert em.n_data == n_freq * n_rx * n_comp * 2  # real + imaginary rows
    codes = set(int(row[0]) for row in em.data)
    assert codes == {113, 114, 115, 116}


def test_prepare_inputs_maps_receiver_coordinates_to_mt_block(tmp_path):
    problem = _problem()
    adapter = _adapter()
    adapter._prepare_inputs(problem, tmp_path)
    from pycsamt.models.mare2dem.iotools.emdata import read_emdata

    em = read_emdata(tmp_path / "mare2dem.emdata")
    np.testing.assert_allclose(em.mt.receivers[:, 1], [200.0, 800.0])  # y
    np.testing.assert_allclose(em.mt.receivers[:, 2], [0.0, 0.0])  # z
    assert em.mt.receiver_name == ["S00", "S01"]


# ---------------------------------------------------------------------------
# solve(): scripted stand-in (no real MARE2DEM binary needed)
# ---------------------------------------------------------------------------

_FAKE_MARE2DEM_SCRIPT = textwrap.dedent(
    """
    import sys
    from pathlib import Path

    run_stem = sys.argv[-1]
    stem = run_stem.rsplit(".", 1)[0]
    workdir = Path.cwd()
    lines = (workdir / f"{stem}.emdata").read_text().splitlines()

    data_lines = []
    in_data = False
    for line in lines:
        s = line.strip()
        if not s or s.startswith("!"):
            continue
        if s.lower().startswith("# data"):
            in_data = True
            continue
        if in_data:
            data_lines.append(s)

    out = ["Format:  EMResp_2.3\\n", f"# Data:       {len(data_lines)}\\n"]
    for line in data_lines:
        parts = line.split()
        code, freq_idx, tx_idx, rx_idx = (int(v) for v in parts[:4])
        data_val, std_err = float(parts[4]), float(parts[5])
        predicted = code + freq_idx * 0.01 + rx_idx * 0.001
        out.append(
            f"{code:7d} {freq_idx:7d} {tx_idx:7d} {rx_idx:7d} "
            f"{data_val:22.15g} {std_err:22.15g} {predicted:20.15g} 0.0\\n"
        )

    (workdir / f"{run_stem}.resp").write_text("".join(out))
    """
)


def _write_fake_mare2dem(tmp_path: Path) -> Path:
    script = tmp_path / "fake_mare2dem.py"
    script.write_text(_FAKE_MARE2DEM_SCRIPT)
    return script


def _expected_value(code: int, freq_idx: int, rx_idx: int) -> float:
    return code + freq_idx * 0.01 + rx_idx * 0.001


def test_solve_end_to_end_against_a_scripted_stand_in(tmp_path):
    script = _write_fake_mare2dem(tmp_path)
    adapter = Mare2DEMAdapter(
        config=Mare2DEMConfig(use_mpi=False),
        run_policy=ExternalRunPolicy(sys.executable, workdir=str(tmp_path)),
    )
    adapter._build_command = lambda problem, workdir, executable, context: [
        str(executable),
        str(script),
        context["run_stem"],
    ]

    problem = _problem()
    result = adapter.solve(problem)

    assert result.shape == (2, 2, 2)
    assert result.success
    assert result.backend_name == "mare2dem"
    for fi in range(2):
        for ri in range(2):
            for ci, component in enumerate(problem.components):
                real_code, imag_code = _DATA_CODES[component]
                expected = complex(
                    _expected_value(real_code, fi + 1, ri + 1),
                    _expected_value(imag_code, fi + 1, ri + 1),
                )
                assert result.impedance_v_a[ri, fi, ci] == pytest.approx(
                    expected
                )
    assert result.diagnostics.runtime_s >= 0.0
    assert np.all(result.diagnostics.iterations == 0)


def test_solve_raises_when_no_response_file_is_produced(tmp_path):
    empty_script = tmp_path / "empty.py"
    empty_script.write_text("pass\n")
    adapter = Mare2DEMAdapter(
        config=Mare2DEMConfig(use_mpi=False),
        run_policy=ExternalRunPolicy(sys.executable, workdir=str(tmp_path)),
    )
    adapter._build_command = lambda problem, workdir, executable, context: [
        str(executable),
        str(empty_script),
        context["run_stem"],
    ]

    with pytest.raises(Exception, match="no MARE2DEM predicted-response"):
        adapter.solve(_problem())


def test_solve_raises_on_ambiguous_response_files(tmp_path):
    script = _write_fake_mare2dem(tmp_path)
    adapter = Mare2DEMAdapter(
        config=Mare2DEMConfig(use_mpi=False),
        run_policy=ExternalRunPolicy(sys.executable, workdir=str(tmp_path)),
    )

    def build_command(problem, workdir, executable, context):
        # a leftover/second .resp file makes auto-detection ambiguous
        (workdir / "stray.resp").write_text("Format:  EMResp_2.3\n")
        return [str(executable), str(script), context["run_stem"]]

    adapter._build_command = build_command
    with pytest.raises(Exception, match="multiple candidate"):
        adapter.solve(_problem())


def test_solve_raises_executable_not_found_without_a_real_binary():
    adapter = _adapter()
    with pytest.raises(ExecutableNotFoundError):
        adapter.solve(_problem())


# ---------------------------------------------------------------------------
# Registry
# ---------------------------------------------------------------------------


def test_register_mare2dem_backend_reports_unavailable_without_a_binary():
    register_mare2dem_backend(replace=True)
    assert "mare2dem" in list_backends()
    assert list_backends()["mare2dem"]["available"] is False
    with pytest.raises(RuntimeError, match="unavailable"):
        create_backend("mare2dem")


# ---------------------------------------------------------------------------
# requires_real_mare2dem: genuine physics against a real compiled binary
# ---------------------------------------------------------------------------

_MARE2DEM_BINARY = os.environ.get("PYCSAMT_MARE2DEM_BINARY")
_HAS_REAL_MARE2DEM = bool(
    _MARE2DEM_BINARY and Path(_MARE2DEM_BINARY).is_file()
)
_HAS_MPIRUN = shutil.which("mpirun") is not None
requires_real_mare2dem = pytest.mark.skipif(
    not (_HAS_REAL_MARE2DEM and _HAS_MPIRUN),
    reason=(
        "no compiled MARE2DEM binary configured; set "
        "PYCSAMT_MARE2DEM_BINARY to a real binary's path (and ensure "
        "mpirun is on PATH) to run this test -- see mare2dem.py's module "
        "docstring for how one was built for the 2026-07-31 validation"
    ),
)


def _half_space_mesh_and_receivers(resistivity_ohm_m: float, frequency_hz: float):
    """A domain padded to several skin depths, matching the safety-factor
    sizing convention already used elsewhere in this codebase (e.g.
    dataset2d.py's _solver_mesh_and_conductivity) -- an inadequately
    padded domain was mistaken for a units bug during initial manual
    validation before the domain was sized properly.
    """
    from pycsamt.forward.maxwell.mesh import skin_depth_m

    half_width = 8.0 * float(
        skin_depth_m(resistivity_ohm_m, frequency_hz)
    )
    depth = half_width
    mesh = TriMesh(
        nodes_m=[
            [-half_width, 0],
            [half_width, 0],
            [half_width, depth],
            [-half_width, depth],
        ],
        triangles=[[0, 1, 2], [0, 2, 3]],
        boundary_segments=[[0, 1], [1, 2], [2, 3], [3, 0]],
        region_ids=[0, 0],
    )
    receivers = ReceiverSet(
        [[-500.0, 0.0], [0.0, 0.0], [500.0, 0.0]], ["S00", "S01", "S02"]
    )
    return mesh, receivers


@requires_real_mare2dem
def test_real_mare2dem_passes_half_space_benchmark():
    resistivity_ohm_m = 100.0
    frequencies = [1.0, 0.1]
    mesh, receivers = _half_space_mesh_and_receivers(
        resistivity_ohm_m, min(frequencies)
    )
    problem = TriProblem(
        mesh,
        np.full(mesh.shape, 1.0 / resistivity_ohm_m),
        frequencies,
        receivers,
    )
    cfg = Mare2DEMConfig(binary=_MARE2DEM_BINARY, n_procs=2)
    adapter = Mare2DEMAdapter(config=cfg)
    result = adapter.solve(problem)

    assert result.success
    errors = []
    for fi, frequency in enumerate(frequencies):
        expected_zxy = half_space_impedance(resistivity_ohm_m, frequency)
        for ri in range(receivers.count):
            zxy, zyx = result.impedance_v_a[ri, fi, :]
            errors.append(abs(zxy - expected_zxy) / abs(expected_zxy))
            errors.append(abs(zyx - (-expected_zxy)) / abs(expected_zxy))
    assert max(errors) < 0.05, (
        f"max relative error {max(errors):.4f} exceeds 5% "
        f"(MARE2DEM's own adaptive refinement reports ~0.2-0.8% for "
        f"this domain, so this threshold has real margin)"
    )
