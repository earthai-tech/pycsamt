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


@pytest.mark.parametrize("cls_name", [
    "Inv2DAgent", "Inv3DAgent", "EnsembleAgent", "JointInversionAgent",
    "HybridInversionAgent", "PINNInversionAgent", "InversionBackendAgent",
    "InversionPrepAgent", "Occam2DAgent", "ModEmAgent",
])
def test_inversion_agents_fail_without_input(cls_name):
    import pycsamt.agents as agents_pkg

    cls = getattr(agents_pkg, cls_name)
    result = mk(cls).execute({})
    assert result.status == "failed"
    assert result.error


# ── file-preparation agents (no training) ────────────────────────────────────


def test_occam2d_agent_prepares_files(sites3, tmp_path):
    from pycsamt.agents import Occam2DAgent

    result = mk(Occam2DAgent).execute(
        {"sites": sites3, "output_dir": str(tmp_path)}
    )
    assert result.status == "success"
    assert result.summary
    assert any(tmp_path.rglob("*")), "expected prepared files on disk"


def test_modem_agent_prepares_files(sites3, tmp_path):
    from pycsamt.agents import ModEmAgent

    result = mk(ModEmAgent).execute(
        {"sites": sites3, "output_dir": str(tmp_path)}
    )
    assert result.status == "success"
    assert result.summary
    assert any(tmp_path.rglob("*")), "expected prepared files on disk"


def test_inversion_prep_agent_occam2d(sites3, tmp_path):
    from pycsamt.agents import InversionPrepAgent

    result = mk(InversionPrepAgent).execute({
        "sites": sites3,
        "code": "occam2d",
        "output_dir": str(tmp_path),
    })
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
    result = ag.execute(
        {"sites": sites3, "output_dir": str(tmp_path)}
    )
    assert result.status in {"success", "failed"}
    assert result.summary or result.error


# ── trainable agents (tiny budgets) ──────────────────────────────────────────


@requires_dl
def test_inv2d_agent_tiny_training(sites3, tmp_path):
    from pycsamt.agents import Inv2DAgent

    ag = mk(Inv2DAgent, epochs=1)
    result = ag.execute(
        {"sites": sites3, "output_dir": str(tmp_path)}
    )
    assert result.status in {"success", "failed"}
    if result.status == "success":
        assert result.summary


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


@requires_dl
def test_pinn_agent_tiny_training(sites3, tmp_path):
    from pycsamt.agents import PINNInversionAgent

    ag = mk(PINNInversionAgent, dim=1, n_layers=4)
    result = ag.execute({
        "sites": sites3,
        "epochs": 1,
        "output_dir": str(tmp_path),
    })
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
