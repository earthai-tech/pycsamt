"""Tests for the external ModEM 3-D adapter.

No compiled ``Mod3DMT`` binary exists in this environment (see the
module docstring), so these tests validate everything that does not
require genuine Maxwell physics from a real binary:

* input-file generation, round-tripped through the real
  :class:`~pycsamt.models.modem.data.ModEmData`/
  :class:`~pycsamt.models.modem.model3d.ModEmModel3D` reader/writer
  (not a second, hand-rolled parser);
* the full external-process pipeline (write inputs -> run -> parse
  output -> build a ``ForwardResult``), exercised end-to-end against a
  small scripted stand-in "ModEM" executable that echoes back a known
  analytic half-space response in ModEM's own field units, so the
  unit-conversion and column-mapping logic is genuinely exercised;
* preflight rejections and the missing-executable path, which is the
  common case for anyone using this adapter without a compiled binary.
"""

from __future__ import annotations

import sys
import textwrap
from pathlib import Path

import numpy as np
import pytest

from pycsamt.forward.maxwell import (
    ExecutableNotFoundError,
    ExternalRunPolicy,
    IncompatibleProblemError,
    MaxwellMesh,
    MaxwellProblem,
    ReceiverSet,
    create_backend,
    list_backends,
)
from pycsamt.forward.maxwell.modem3d import (
    _MU0,
    ModEm3DAdapter,
    _unit_conversion_to_si,
    register_modem3d_backend,
)
from pycsamt.models.modem.config import ModEmConfig
from pycsamt.models.modem.data import ModEmData
from pycsamt.models.modem.model3d import ModEmModel3D


def _problem(**changes):
    mesh = MaxwellMesh([0, 1000, 2000], [0, 500, 1000], [0, 1000, 2000])
    values = dict(
        mesh=mesh,
        conductivity_s_m=np.full(mesh.shape, 0.01),
        frequencies_hz=[1.0, 10.0],
        receivers=ReceiverSet([[1000.0, 1000.0, 0.0]], ["S00"]),
        components=("zxx", "zxy", "zyx", "zyy"),
    )
    values.update(changes)
    return MaxwellProblem(**values)


# --------------------------------------------------------------------------- #
# Unit conversion
# --------------------------------------------------------------------------- #


def test_unit_conversion_recognizes_field_and_si_units():
    assert _unit_conversion_to_si("[mV/km]/[nT]") == pytest.approx(
        4.0e-4 * np.pi
    )
    assert _unit_conversion_to_si("[V/A]") == 1.0
    assert _unit_conversion_to_si("ohm") == 1.0


def test_unit_conversion_rejects_unknown_units():
    with pytest.raises(ValueError, match="unrecognized"):
        _unit_conversion_to_si("furlongs/fortnight")


# --------------------------------------------------------------------------- #
# Capabilities and assess()
# --------------------------------------------------------------------------- #


def test_capabilities_declare_3d_full_tensor_and_nonuniform_support():
    adapter = ModEm3DAdapter(
        run_policy=ExternalRunPolicy("does-not-exist-xyz")
    )
    cap = adapter.capabilities
    assert cap.name == "modem3d"
    assert cap.dimensions == (3,)
    assert set(cap.components) == {"zxx", "zxy", "zyx", "zyy"}
    assert cap.supports_nonuniform_mesh is True
    assert cap.verified_benchmarks == ()


def test_assess_rejects_buried_receivers():
    adapter = ModEm3DAdapter(
        run_policy=ExternalRunPolicy("does-not-exist-xyz")
    )
    problem = _problem(
        receivers=ReceiverSet([[1000.0, 1000.0, 50.0]], ["S00"])
    )
    report = adapter.assess(problem)
    assert not report.compatible
    assert any("surface" in message for message in report.errors)


def test_assess_rejects_receivers_outside_horizontal_bounds():
    adapter = ModEm3DAdapter(
        run_policy=ExternalRunPolicy("does-not-exist-xyz")
    )
    problem = _problem(
        receivers=ReceiverSet([[50_000.0, 1000.0, 0.0]], ["S00"])
    )
    report = adapter.assess(problem)
    assert not report.compatible
    assert any("mesh" in message for message in report.errors)


def test_assess_rejects_non_surface_aligned_mesh():
    adapter = ModEm3DAdapter(
        run_policy=ExternalRunPolicy("does-not-exist-xyz")
    )
    mesh = MaxwellMesh(
        [0, 1000, 2000], [100, 500, 1000], [0, 1000, 2000]
    )
    problem = _problem(
        mesh=mesh, conductivity_s_m=np.full(mesh.shape, 0.01)
    )
    report = adapter.assess(problem)
    assert not report.compatible
    assert any("top z edge" in message for message in report.errors)


def test_assess_rejects_non_vacuum_permeability():
    adapter = ModEm3DAdapter(
        run_policy=ExternalRunPolicy("does-not-exist-xyz")
    )
    problem = _problem(magnetic_permeability_h_m=_MU0 * 2.0)
    report = adapter.assess(problem)
    assert not report.compatible
    assert any("permeability" in message for message in report.errors)


def test_nonuniform_mesh_is_accepted_unlike_mt3d():
    # Regression guard: modem3d must NOT inherit mt3d's uniform-only
    # restriction, since ModEM itself has no such limitation.
    adapter = ModEm3DAdapter(
        run_policy=ExternalRunPolicy("does-not-exist-xyz")
    )
    mesh = MaxwellMesh([0, 500, 2000], [0, 500, 1000], [0, 800, 2000])
    problem = _problem(
        mesh=mesh, conductivity_s_m=np.full(mesh.shape, 0.01)
    )
    report = adapter.assess(problem)
    assert report.compatible


def test_incompatible_2d_problem_is_rejected_before_solving():
    mesh = MaxwellMesh([0, 1, 2], [0, 1, 2])
    problem = MaxwellProblem(
        mesh,
        np.full(mesh.shape, 0.01),
        [1.0],
        ReceiverSet([[0.5, 0.0]], ["S"]),
        ("zxy",),
    )
    adapter = ModEm3DAdapter(
        run_policy=ExternalRunPolicy("does-not-exist-xyz")
    )
    with pytest.raises(IncompatibleProblemError, match="2-D"):
        adapter.solve(problem)


# --------------------------------------------------------------------------- #
# Missing executable: the common case without a compiled binary
# --------------------------------------------------------------------------- #


def test_missing_binary_raises_executable_not_found():
    adapter = ModEm3DAdapter(
        run_policy=ExternalRunPolicy("does-not-exist-xyz")
    )
    with pytest.raises(ExecutableNotFoundError):
        adapter.solve(_problem())


# --------------------------------------------------------------------------- #
# _prepare_inputs: real ModEM file round trips (no subprocess involved)
# --------------------------------------------------------------------------- #


def test_prepare_inputs_writes_a_valid_modem_model_and_data_file(tmp_path):
    problem = _problem(
        conductivity_s_m=np.full(_problem().mesh.shape, 1.0 / 250.0)
    )
    adapter = ModEm3DAdapter(
        run_policy=ExternalRunPolicy("does-not-exist-xyz")
    )
    context = adapter._prepare_inputs(problem, tmp_path)

    model = ModEmModel3D.read(context["model_path"])
    assert model.shape == problem.mesh.shape
    # %.5E ASCII round-trip through log-space loses a little precision.
    np.testing.assert_allclose(model.rho_linear, 250.0, rtol=1e-5)
    assert model.n_air == 0

    data = ModEmData.read(context["data_path"])
    assert data.n_sites == 1
    assert data.site_names == ["S00"]
    # the mesh's own x/y minimum edge is (0, 0) here, so the shift is
    # a no-op and the receiver's coordinates pass through unchanged.
    x, y, z = data.site_coords["S00"]
    assert x == pytest.approx(1000.0)
    assert y == pytest.approx(1000.0)
    assert z == 0.0
    assert sorted(data.periods) == sorted(1.0 / np.array([1.0, 10.0]))
    assert data.blocks[0]["component_type"] == "Full_Impedance"
    assert len(data.blocks[0]["rows"]) == 2 * 1 * 4  # periods*sites*comps


def test_prepare_inputs_shifts_coordinates_to_the_mesh_origin(tmp_path):
    # MaxwellMesh(x_edges, z_edges, y_edges): x0=500, y0=100.
    mesh = MaxwellMesh(
        [500, 1500, 2500], [200, 700, 1200], [100, 1000, 2000]
    )
    problem = _problem(
        mesh=mesh,
        conductivity_s_m=np.full(mesh.shape, 0.01),
        receivers=ReceiverSet([[500.0, 100.0, 0.0]], ["A"]),
    )
    adapter = ModEm3DAdapter(
        run_policy=ExternalRunPolicy("does-not-exist-xyz")
    )
    context = adapter._prepare_inputs(problem, tmp_path)
    data = ModEmData.read(context["data_path"])
    x, y, _z = data.site_coords["A"]
    assert x == pytest.approx(0.0)
    assert y == pytest.approx(0.0)


# --------------------------------------------------------------------------- #
# _build_command
# --------------------------------------------------------------------------- #


def test_build_command_without_mpi():
    adapter = ModEm3DAdapter(
        run_policy=ExternalRunPolicy("does-not-exist-xyz")
    )
    context = {"model_path": Path("model.ws"), "data_path": Path("data.dat")}
    command = adapter._build_command(
        None, Path("."), Path("Mod3DMT"), context
    )
    assert command == ["Mod3DMT", "-F", "model.ws", "data.dat"]


def test_build_command_with_mpi_resolves_the_modem_binary(tmp_path):
    binary = tmp_path / "Mod3DMT"
    binary.write_text("x")
    cfg = ModEmConfig(mode="3d", use_mpi=True, n_procs=4, binary_3d="Mod3DMT")
    adapter = ModEm3DAdapter(
        config=cfg,
        run_policy=ExternalRunPolicy(
            "mpirun", search_paths=(str(tmp_path),)
        ),
    )
    context = {"model_path": Path("model.ws"), "data_path": Path("data.dat")}
    command = adapter._build_command(
        None, tmp_path, Path("mpirun"), context
    )
    assert command[:3] == ["mpirun", "-np", "4"]
    assert command[3] == str(binary)
    assert command[4:] == ["-F", "model.ws", "data.dat"]


# --------------------------------------------------------------------------- #
# End-to-end against a scripted stand-in "ModEM" executable
# --------------------------------------------------------------------------- #

_FAKE_MODEM_SCRIPT = textwrap.dedent(
    """
    # A dependency-free stand-in for Mod3DMT -F: parses the ModEM
    # data-file format directly (no pycsamt import -- importing the
    # package does an os.mkdir() relative to the current working
    # directory, which is this scripted run's temp workdir, not the
    # repo root) and echoes back a known analytic half-space response
    # in ModEM's own field units.
    import sys

    MU0 = 4.0e-7 * 3.141592653589793
    FIELD_PER_SI = 1.0 / (4.0e-4 * 3.141592653589793)
    RHO = {rho}

    _, _flag, model_name, data_name = sys.argv[:4]

    header_lines = []
    rows = []
    with open(data_name) as fh:
        for line in fh:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            if stripped.startswith(">"):
                header_lines.append(stripped)
                continue
            parts = stripped.split()
            rows.append(parts)

    out_lines = ["# predicted by fake_modem.py\\n"] + [
        h + "\\n" for h in header_lines
    ]
    for parts in rows:
        period = float(parts[0])
        code = parts[1]
        lat, lon = parts[2], parts[3]
        x, y, z = parts[4], parts[5], parts[6]
        comp = parts[7]
        err = parts[10]
        omega = 2.0 * 3.141592653589793 / period
        magnitude = (omega * MU0 * RHO) ** 0.5
        real_part = magnitude / 2.0**0.5
        imag_part = real_part
        if comp == "ZXY":
            value_real, value_imag = real_part, imag_part
        elif comp == "ZYX":
            value_real, value_imag = -real_part, -imag_part
        else:
            value_real, value_imag = 0.0, 0.0
        value_real *= FIELD_PER_SI
        value_imag *= FIELD_PER_SI
        out_lines.append(
            f"{{period:>15.6E}} {{code:>8}}  {{lat:>8}}  {{lon:>8}}  "
            f"{{x:>14}}  {{y:>14}}  {{z:>10}}  {{comp:<6}}  "
            f"{{value_real:>15.6E}}  {{value_imag:>15.6E}}  {{err:>15}}\\n"
        )

    with open("predicted.dat", "w") as fh:
        fh.writelines(out_lines)
    """
)


def _write_fake_modem(tmp_path: Path, rho: float) -> Path:
    script = tmp_path / "fake_modem.py"
    script.write_text(_FAKE_MODEM_SCRIPT.format(rho=rho))
    return script


def test_solve_end_to_end_against_a_scripted_stand_in(tmp_path):
    rho = 250.0
    script = _write_fake_modem(tmp_path, rho)
    adapter = ModEm3DAdapter(
        run_policy=ExternalRunPolicy(
            sys.executable, search_paths=(), workdir=None
        )
    )
    # Patch the command to run our script rather than a real binary.
    adapter._build_command = lambda problem, workdir, executable, context: [
        str(executable),
        str(script),
        "-F",
        context["model_path"].name,
        context["data_path"].name,
    ]

    mesh = MaxwellMesh([0, 1000, 2000], [0, 1000, 2000], [0, 1000, 2000])
    problem = MaxwellProblem(
        mesh,
        np.full(mesh.shape, 1.0 / rho),
        [1.0, 4.0],
        ReceiverSet([[1000.0, 1000.0, 0.0]], ["S00"]),
        ("zxy", "zyx", "zxx"),
    )
    result = adapter.solve(problem)

    assert result.shape == (1, 2, 3)
    assert result.success
    for fi, frequency in enumerate(problem.frequencies_hz):
        omega = 2.0 * np.pi * frequency
        expected = np.sqrt(1j * omega * _MU0 * rho)
        np.testing.assert_allclose(
            result.impedance_v_a[0, fi, 0], expected, rtol=1e-6
        )
        np.testing.assert_allclose(
            result.impedance_v_a[0, fi, 1], -expected, rtol=1e-6
        )
        assert result.impedance_v_a[0, fi, 2] == 0j
    assert result.diagnostics.runtime_s >= 0.0
    assert np.all(result.diagnostics.iterations == 0)


def test_solve_raises_on_ambiguous_predicted_output(tmp_path):
    script_body = textwrap.dedent(
        """
        import shutil
        import sys
        shutil.copy(sys.argv[3], "extra_one.dat")
        shutil.copy(sys.argv[3], "extra_two.dat")
        """
    )
    script = tmp_path / "ambiguous_modem.py"
    script.write_text(script_body)
    adapter = ModEm3DAdapter(run_policy=ExternalRunPolicy(sys.executable))
    adapter._build_command = lambda problem, workdir, executable, context: [
        str(executable),
        str(script),
        "-F",
        context["model_path"].name,
        context["data_path"].name,
    ]
    with pytest.raises(Exception, match="multiple candidate"):
        adapter.solve(_problem())


# --------------------------------------------------------------------------- #
# Registry integration
# --------------------------------------------------------------------------- #


def test_register_modem3d_backend_reports_unavailable_without_a_binary():
    register_modem3d_backend(
        config=ModEmConfig(mode="3d", binary_3d="does-not-exist-xyz"),
        replace=True,
    )
    assert "modem3d" in list_backends()
    assert list_backends()["modem3d"]["available"] is False
    with pytest.raises(RuntimeError, match="unavailable"):
        create_backend("modem3d")


def test_register_modem3d_backend_is_creatable_once_available(tmp_path):
    fake_binary = tmp_path / "Mod3DMT"
    fake_binary.write_text("placeholder")
    register_modem3d_backend(
        config=ModEmConfig(mode="3d", binary_3d=str(fake_binary)),
        replace=True,
    )
    assert list_backends()["modem3d"]["available"] is True
    backend = create_backend("modem3d")
    assert isinstance(backend, ModEm3DAdapter)
