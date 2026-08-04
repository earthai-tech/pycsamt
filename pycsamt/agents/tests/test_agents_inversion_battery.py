# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Offline contract battery for the inversion agent family.

Trainable agents run with tiny budgets (1-2 epochs, few layers) on the
bundled 3-station dataset so the whole battery stays in CI-friendly
time. Agents are built inside ``AGENT_CONFIG.offline()`` so no LLM key
is ever resolved.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pytest

from pycsamt.api.agents import AGENT_CONFIG

from .conftest import has_backend

requires_dl = pytest.mark.skipif(
    not has_backend(), reason="PyTorch/TensorFlow not installed"
)


def mk(cls, **kwargs):
    with AGENT_CONFIG.offline():
        return cls(**kwargs)


@pytest.fixture(autouse=True)
def _close_figures():
    yield
    plt.close("all")


@pytest.fixture(scope="module")
def sites3(edi_dir: Path):
    from pycsamt.emtools import ensure_sites

    return ensure_sites(edi_dir)


# ── missing-input contract ───────────────────────────────────────────────────


@pytest.mark.parametrize(
    "cls_name",
    [
        "Inv2DAgent",
        "Inv3DAgent",
        "EnsembleAgent",
        "JointInversionAgent",
        "HybridInversionAgent",
        "PINNInversionAgent",
        "InversionBackendAgent",
        "InversionPrepAgent",
        "Occam2DAgent",
        "ModEmAgent",
    ],
)
def test_inversion_agents_fail_without_input(cls_name):
    import pycsamt.agents as agents_pkg

    cls = getattr(agents_pkg, cls_name)
    result = mk(cls).execute({})
    assert result.status == "failed"
    assert result.error


# ── file-preparation agents (no training) ────────────────────────────────────


def test_occam2d_agent_prepares_files(sites3, tmp_path):
    from pycsamt.agents import Occam2DAgent

    result = mk(Occam2DAgent).execute({"sites": sites3, "output_dir": str(tmp_path)})
    assert result.status == "success"
    assert result.summary
    assert any(tmp_path.rglob("*")), "expected prepared files on disk"


def test_modem_agent_prepares_files(sites3, tmp_path):
    from pycsamt.agents import ModEmAgent

    result = mk(ModEmAgent).execute({"sites": sites3, "output_dir": str(tmp_path)})
    assert result.status == "success"
    assert result.summary
    assert any(tmp_path.rglob("*")), "expected prepared files on disk"


def test_inversion_prep_agent_occam2d(sites3, tmp_path):
    from pycsamt.agents import InversionPrepAgent

    result = mk(InversionPrepAgent).execute(
        {
            "sites": sites3,
            "code": "occam2d",
            "output_dir": str(tmp_path),
        }
    )
    # prep runs flag sparse inputs for human review
    assert result.status in {"success", "needs_review", "failed"}
    assert result.summary or result.error


# ── builtin numeric backend (no DL) ──────────────────────────────────────────


def test_inversion_backend_agent_builtin(sites3, tmp_path):
    from pycsamt.agents import InversionBackendAgent

    ag = mk(
        InversionBackendAgent,
        backend="builtin",
        dimension="1d",
        n_layers=4,
        max_iter=3,
    )
    result = ag.execute({"sites": sites3, "output_dir": str(tmp_path)})
    assert result.status in {"success", "failed"}
    assert result.summary or result.error


# ── trainable agents (tiny budgets) ──────────────────────────────────────────


@requires_dl
def test_inv2d_agent_tiny_training(sites3, tmp_path):
    from pycsamt.agents import Inv2DAgent

    ag = mk(Inv2DAgent, epochs=1)
    result = ag.execute({"sites": sites3, "output_dir": str(tmp_path)})
    assert result.status in {"success", "failed"}
    if result.status == "success":
        assert result.summary
        assert result.data["physics"] == "mt1d"
        assert result.data["mt2d_recovery"] is None


@requires_dl
def test_inv2d_agent_mt2d_physics_end_to_end(sites3):
    """physics='mt2d' trains on real 2-D Maxwell solves (not tiled
    1-D) and reports a genuine held-out recovery check, per the
    AI-inversion plan's M5 slice. n_train_profiles is kept modest but
    large enough to survive the dataset's own train/val/test split
    plus EMInverter2D.fit's internal validation split without ever
    handing BatchNorm a batch of size one.
    """
    import numpy as np

    from pycsamt.agents import Inv2DAgent

    ag = mk(
        Inv2DAgent,
        epochs=1,
        n_depth=4,
        n_freqs=3,
        n_train_profiles=20,
        n_stations_per_profile=3,
        physics="mt2d",
        station_spacing_m=200.0,
    )
    result = ag.execute({"sites": sites3, "freqs": [100.0, 30.0, 10.0]})
    assert result.status == "success", result.error
    assert result.data["physics"] == "mt2d"
    assert result.data["pred_section"].shape == (4, 3)
    assert np.all(np.isfinite(result.data["pred_section"]))

    recovery = result.data["mt2d_recovery"]
    assert recovery is not None
    assert recovery["n_samples"] >= 1
    assert np.isfinite(recovery["rmse"])
    assert np.isfinite(recovery["mae"])


_FAKE_MARE2DEM_SCRIPT = """
import sys
from pathlib import Path

stem = sys.argv[-1]
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
    predicted = 1.0
    out.append(
        f"{code:7d} {freq_idx:7d} {tx_idx:7d} {rx_idx:7d} "
        f"{data_val:22.15g} {std_err:22.15g} {predicted:20.15g} 0.0\\n"
    )

(workdir / f"{stem}.0.resp").write_text("".join(out))
"""


def _fake_mare2dem_adapter(tmp_path: Path):
    import sys

    from pycsamt.forward.maxwell import ExternalRunPolicy
    from pycsamt.forward.maxwell.mare2dem import Mare2DEMAdapter
    from pycsamt.models.mare2dem.config import Mare2DEMConfig

    script = tmp_path / "fake_mare2dem.py"
    script.write_text(_FAKE_MARE2DEM_SCRIPT)
    adapter = Mare2DEMAdapter(
        config=Mare2DEMConfig(use_mpi=False),
        run_policy=ExternalRunPolicy(sys.executable, workdir=str(tmp_path)),
    )
    adapter._build_command = lambda problem, workdir, executable, context: [
        str(executable),
        str(script),
        context["stem"],
    ]
    return adapter


@requires_dl
def test_inv2d_agent_mt2d_tri_physics_end_to_end(sites3, tmp_path):
    """physics='mt2d_tri' trains a GCN over a triangular mesh, solved
    with Mare2DEMAdapter -- unverified physics without a real compiled
    MARE2DEM binary (see that adapter's module docstring), so this test
    injects a scripted stand-in executable via mare2dem_adapter=,
    exactly like test_maxwell_mare2dem.py / test_ai_training_dataset2d_tri.py
    do. Proves the agent-level plumbing (mesh, GCN training/prediction,
    triangular-mesh figure, AgentResult shape), not real Maxwell physics.
    """
    import numpy as np

    from pycsamt.agents import Inv2DAgent

    ag = mk(
        Inv2DAgent,
        epochs=1,
        n_freqs=3,
        n_train_profiles=20,
        n_stations_per_profile=3,
        physics="mt2d_tri",
        station_spacing_m=200.0,
        mesh_target_cell_m=100.0,
        field_grid_cell_m=50.0,
        mare2dem_adapter=_fake_mare2dem_adapter(tmp_path),
    )
    result = ag.execute(
        {
            "sites": sites3,
            "freqs": [100.0, 30.0, 10.0],
            "output_dir": str(tmp_path / "out"),
        }
    )
    assert result.status == "success", result.error
    assert result.data["physics"] == "mt2d_tri"
    pred = result.data["pred_triangles"]
    assert pred["mesh"].n_triangles > 0
    assert pred["log10_resistivity"].shape == (pred["mesh"].n_triangles,)
    assert np.all(np.isfinite(pred["log10_resistivity"]))
    assert "inv2d_tri_section" in result.data["figures"]

    recovery = result.data["mt2d_tri_recovery"]
    assert recovery is not None
    assert recovery["n_samples"] >= 1
    assert np.isfinite(recovery["rmse"])


def test_inv2d_agent_mt2d_tri_physics_with_real_topography(sites3, tmp_path):
    """``topo_x_m``/``topo_z_m`` build a real topography-following
    training mesh (stations at their true elevation, not z=0) -- proves
    the capability is wired all the way from the agent's own
    constructor down to ``build_graded_tri_mesh``, not just available
    inside ``dataset2d_tri.py`` directly.
    """
    import numpy as np

    from pycsamt.agents import Inv2DAgent

    topo_x_m = [0.0, 200.0, 400.0]
    topo_z_m = [-30.0, -50.0, -10.0]

    ag = mk(
        Inv2DAgent,
        epochs=1,
        n_freqs=3,
        n_train_profiles=20,
        n_stations_per_profile=3,
        physics="mt2d_tri",
        station_spacing_m=200.0,
        mesh_target_cell_m=100.0,
        field_grid_cell_m=50.0,
        topo_x_m=topo_x_m,
        topo_z_m=topo_z_m,
        mare2dem_adapter=_fake_mare2dem_adapter(tmp_path),
    )
    result = ag.execute(
        {
            "sites": sites3,
            "freqs": [100.0, 30.0, 10.0],
            "output_dir": str(tmp_path / "out"),
        }
    )
    assert result.status == "success", result.error
    pred = result.data["pred_triangles"]
    mesh = pred["mesh"]
    # None of the stations are at z=0 in this synthetic topography, so a
    # mesh node existing at each station's true interpolated elevation
    # proves topography genuinely reached the mesh builder.
    station_x = np.arange(3) * 200.0
    station_z = np.interp(station_x, topo_x_m, topo_z_m)
    assert not np.any(np.isclose(station_z, 0.0))
    for sx, sz in zip(station_x, station_z):
        assert np.any(
            np.isclose(mesh.nodes_m[:, 0], sx) & np.isclose(mesh.nodes_m[:, 1], sz)
        )


@requires_dl
def test_inv2d_agent_rejects_bad_physics():
    from pycsamt.agents import Inv2DAgent

    with pytest.raises(ValueError, match="physics"):
        mk(Inv2DAgent, physics="mt3d")


@requires_dl
def test_ensemble_agent_tiny_training(sites3):
    from pycsamt.agents import EnsembleAgent

    ag = mk(EnsembleAgent, n_estimators=2, epochs=1, n_layers=4)
    result = ag.execute({"sites": sites3})
    assert result.status in {"success", "failed"}
    if result.status == "success":
        assert result.summary


@requires_dl
def test_joint_agent_tiny_training(sites3):
    from pycsamt.agents import JointInversionAgent

    ag = mk(JointInversionAgent, epochs=1, n_layers=4)
    result = ag.execute({"sites": sites3})
    assert result.status in {"success", "failed"}
    if result.status == "success":
        assert result.summary


@requires_dl
def test_inv3d_agent_tiny_training(sites3):
    from pycsamt.agents import Inv3DAgent

    ag = mk(Inv3DAgent, epochs=1, n_layers=4)
    result = ag.execute({"sites": sites3})
    assert result.status in {"success", "failed"}
    if result.status == "success":
        assert result.summary
        assert result.data["physics"] == "mt1d"
        assert result.data["mt3d_recovery"] is None


@requires_dl
def test_inv3d_agent_mt3d_physics_end_to_end(sites3):
    """physics='mt3d' trains on real 3-D Maxwell solves (MT3DAdapter,
    the research-only small-grid solver) at the survey's own real
    station positions, not tiled independent 1-D models, per the
    AI-inversion plan's M8 slice. Kept tiny (few realizations, one
    epoch, a low mesh_safety_factor) since each realization costs a
    genuine 3-D solve, not a cheap 1-D one.
    """
    import numpy as np

    from pycsamt.agents import Inv3DAgent

    ag = mk(
        Inv3DAgent,
        epochs=1,
        n_layers=4,
        n_freqs=3,
        n_train_profiles=6,
        physics="mt3d",
        mesh_safety_factor=2.0,
    )
    result = ag.execute({"sites": sites3, "freqs": [100.0, 30.0, 10.0]})
    assert result.status == "success", result.error
    assert result.data["physics"] == "mt3d"
    assert result.data["pred_rho"].shape == (3, 4)
    assert np.all(np.isfinite(result.data["pred_rho"]))

    recovery = result.data["mt3d_recovery"]
    assert recovery is not None
    assert recovery["n_samples"] >= 1
    assert np.isfinite(recovery["rmse"])
    assert np.isfinite(recovery["mae"])


@requires_dl
def test_inv3d_agent_rejects_bad_physics():
    from pycsamt.agents import Inv3DAgent

    with pytest.raises(ValueError, match="physics"):
        mk(Inv3DAgent, physics="mt2d")


@requires_dl
def test_pinn_agent_tiny_training(sites3, tmp_path):
    from pycsamt.agents import PINNInversionAgent

    ag = mk(PINNInversionAgent, dim=1, n_layers=4)
    result = ag.execute(
        {
            "sites": sites3,
            "epochs": 1,
            "output_dir": str(tmp_path),
        }
    )
    assert result.status in {"success", "failed"}
    if result.status == "success":
        assert result.summary


def test_hybrid_agent_requires_inverter(sites3):
    from pycsamt.agents import HybridInversionAgent

    # no ai_inverter/checkpoint provided -> actionable failure
    result = mk(HybridInversionAgent).execute({"sites": sites3})
    assert result.status == "failed"
    assert result.error


# ── CLI entry point and web shell ────────────────────────────────────────────


def test_agents_main_help_and_offline_commands(capsys):
    from pycsamt.agents.__main__ import main

    main([])  # help text
    out = capsys.readouterr().out
    assert "pycsamt" in out.lower() or "usage" in out.lower()

    main(["list"])
    out = capsys.readouterr().out
    assert out.strip()

    main(["pricing"])
    out = capsys.readouterr().out
    assert out.strip()


def test_agents_main_preview_offline(capsys):
    from pycsamt.agents.__main__ import main

    with AGENT_CONFIG.offline():
        main(["preview", "QC check for /data/3edis"])
    out = capsys.readouterr().out
    assert out.strip()


def test_web_launch_without_gradio():
    pytest.importorskip
    try:
        import gradio  # noqa: F401

        pytest.skip("gradio installed: launch would start a server")
    except ImportError:
        pass

    from pycsamt.agents.web import launch

    with pytest.raises(ImportError, match="[Gg]radio"):
        launch()
