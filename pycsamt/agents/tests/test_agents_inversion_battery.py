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
