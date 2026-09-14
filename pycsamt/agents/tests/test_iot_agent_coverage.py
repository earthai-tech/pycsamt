# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Coverage-focused tests for :mod:`pycsamt.agents.iot_agent`.

Complements ``test_iot_agent.py`` (which already covers the packets happy
path, missing-source failure, session-object input, and manifest write to
``output_dir``) with: the ``assess()`` exception guard, the "ok" vs
"warn/critical" status-branch, deployment/sync/power optional tables, the
``session``-as-dict and ``packets``/``edis`` resolution branches (including
their exception paths), the manifest ``to_manifest``/sign/write exception
guards, the figure-``None``/exception/save-success branches, and the LLM
interpretation "no api_key" / query-exception paths.
"""

from __future__ import annotations

import types

import pytest

pytestmark = pytest.mark.filterwarnings("ignore::RuntimeWarning")


def _sim_packets(n_stations=2, n_samples=128, seed=3):
    from pycsamt.iot import simulate_iot_network

    return simulate_iot_network(
        n_stations=n_stations, n_samples=n_samples, seed=seed
    )


# ── status-assessment branches ───────────────────────────────────────────────


def test_assess_exception_is_recorded(monkeypatch):
    from pycsamt.agents import IoTFieldAgent
    from pycsamt.iot import FieldSession

    monkeypatch.setattr(
        FieldSession,
        "assess",
        lambda self, **k: (_ for _ in ()).throw(RuntimeError("assess boom")),
    )
    result = IoTFieldAgent().execute(
        {"packets": _sim_packets(), "figures": False}
    )
    assert any("assess:" in w for w in result.warnings)
    assert result["level"] == "unknown"
    assert result.status == "needs_review"


def test_level_ok_is_success(monkeypatch):
    from pycsamt.agents import IoTFieldAgent
    from pycsamt.iot import FieldSession

    fake_status = types.SimpleNamespace(level="ok", issues=[])
    monkeypatch.setattr(
        FieldSession, "assess", lambda self, **k: fake_status
    )
    result = IoTFieldAgent().execute(
        {"packets": _sim_packets(), "figures": False}
    )
    assert result["level"] == "ok"
    assert result.status == "success"


def test_level_critical_is_needs_review(monkeypatch):
    from pycsamt.agents import IoTFieldAgent
    from pycsamt.iot import FieldSession

    fake_status = types.SimpleNamespace(
        level="critical", issues=["battery_min_v below 11.0 V"]
    )
    monkeypatch.setattr(
        FieldSession, "assess", lambda self, **k: fake_status
    )
    result = IoTFieldAgent().execute(
        {"packets": _sim_packets(), "figures": False}
    )
    assert result["level"] == "critical"
    assert result.status == "needs_review"


# ── optional table branches ──────────────────────────────────────────────────


def test_deployment_table_is_built():
    from pycsamt.agents import IoTFieldAgent
    from pycsamt.iot import DeploymentConfig, DeviceConfig

    deployment = DeploymentConfig(
        survey_id="survey-a",
        devices=[DeviceConfig("node-1", station="S01", channels=["ex"])],
    )
    result = IoTFieldAgent().execute(
        {
            "packets": _sim_packets(),
            "deployment": deployment,
            "figures": False,
        }
    )
    assert "deployment_table" in result


def test_sync_table_is_built():
    from pycsamt.agents import IoTFieldAgent
    from pycsamt.iot import simulate_gps_drift

    sim = simulate_gps_drift(
        100, sample_interval_s=1.0, drift_ppm=2.0, jitter_ms=0.05, seed=1
    )
    result = IoTFieldAgent().execute(
        {
            "packets": _sim_packets(),
            "sync_references": {"node-1": (sim["local"], sim["reference"])},
            "figures": False,
        }
    )
    assert "sync_table" in result


def test_power_table_is_built():
    from pycsamt.agents import IoTFieldAgent
    from pycsamt.iot import EnergyConfig

    result = IoTFieldAgent().execute(
        {
            "packets": _sim_packets(),
            "energy_configs": [
                EnergyConfig(
                    battery_wh=80.0, active_power_w=1.0, duty_cycle=0.3
                )
            ],
            "figures": False,
        }
    )
    assert "power_table" in result


# ── session-resolution branches ──────────────────────────────────────────────


def test_session_given_as_dict():
    from pycsamt.agents import IoTFieldAgent
    from pycsamt.iot import FieldSession

    session = FieldSession("S1")
    session.add_packets(_sim_packets())
    as_dict = session.to_dict()
    result = IoTFieldAgent().execute(
        {"session": as_dict, "figures": False}
    )
    assert result.get("session") is not None


def test_session_dict_from_dict_exception_falls_through(monkeypatch):
    from pycsamt.agents import IoTFieldAgent
    from pycsamt.iot import FieldSession

    monkeypatch.setattr(
        FieldSession,
        "from_dict",
        classmethod(
            lambda cls, d: (_ for _ in ()).throw(ValueError("bad dict"))
        ),
    )
    result = IoTFieldAgent().execute(
        {"session": {"bogus": True}, "figures": False}
    )
    assert result.status == "failed"


def test_packets_add_exception_falls_through(monkeypatch):
    from pycsamt.agents import IoTFieldAgent
    from pycsamt.iot import FieldSession

    monkeypatch.setattr(
        FieldSession,
        "add_packets",
        lambda self, pkts: (_ for _ in ()).throw(RuntimeError("bad packets")),
    )
    result = IoTFieldAgent().execute(
        {"packets": _sim_packets(), "figures": False}
    )
    assert result.status == "failed"


def test_edis_branch_success(edi_dir):
    from pycsamt.agents import IoTFieldAgent

    result = IoTFieldAgent().execute(
        {"path": str(edi_dir), "figures": False}
    )
    assert result.get("session") is not None
    assert any("re-occupation session" in w for w in result.warnings)


def test_edis_branch_exception(monkeypatch, edi_dir):
    import pycsamt.iot as iot_pkg
    from pycsamt.agents import IoTFieldAgent

    # `_resolve_session` does `from ..iot import field_session_from_edis`
    # (the package re-export), so patch it there rather than on the
    # defining `pycsamt.iot.bridge` module (see the module-level-import
    # monkeypatch gotcha: the package binding is a separate reference).
    monkeypatch.setattr(
        iot_pkg,
        "field_session_from_edis",
        lambda *a, **k: (_ for _ in ()).throw(RuntimeError("edi boom")),
    )
    result = IoTFieldAgent().execute(
        {"path": str(edi_dir), "figures": False}
    )
    assert result.status == "failed"
    assert any("field_session_from_edis" in w for w in result.warnings)


# ── manifest branches ─────────────────────────────────────────────────────────


def test_to_manifest_exception_is_recorded(monkeypatch):
    from pycsamt.agents import IoTFieldAgent
    from pycsamt.iot import FieldSession

    monkeypatch.setattr(
        FieldSession,
        "to_manifest",
        lambda self: (_ for _ in ()).throw(RuntimeError("manifest boom")),
    )
    result = IoTFieldAgent().execute(
        {"packets": _sim_packets(), "manifest": True, "figures": False}
    )
    assert any("to_manifest:" in w for w in result.warnings)
    assert "manifest" not in result.data


def test_sign_key_writes_signed_manifest(tmp_output):
    from pycsamt.agents import IoTFieldAgent

    result = IoTFieldAgent().execute(
        {
            "packets": _sim_packets(),
            "manifest": True,
            "sign_key": "s3cr3t",
            "output_dir": str(tmp_output),
            "figures": False,
        }
    )
    assert "signature" in result
    assert result.get("manifest_path") is not None


def test_write_manifest_exception_is_recorded(monkeypatch):
    from pycsamt.agents import IoTFieldAgent
    from pycsamt.iot import FieldSession

    class _FakeManifest:
        def write(self, path):
            raise OSError("disk full")

    monkeypatch.setattr(
        FieldSession, "to_manifest", lambda self: _FakeManifest()
    )
    result = IoTFieldAgent().execute(
        {
            "packets": _sim_packets(),
            "manifest_path": "somewhere.json",
            "figures": False,
        }
    )
    assert any("write manifest:" in w for w in result.warnings)


# ── figure branches ───────────────────────────────────────────────────────────


def test_figure_returning_none_is_skipped(monkeypatch):
    import pycsamt.iot as iot_pkg
    from pycsamt.agents import IoTFieldAgent

    monkeypatch.setattr(iot_pkg, "plot_field_dashboard", lambda *a, **k: None)
    result = IoTFieldAgent().execute(
        {"packets": _sim_packets(), "figures": True}
    )
    assert "dashboard" not in (result.get("figures") or {})


def test_figure_exception_is_captured(monkeypatch):
    import pycsamt.iot as iot_pkg
    from pycsamt.agents import IoTFieldAgent

    monkeypatch.setattr(
        iot_pkg,
        "plot_field_dashboard",
        lambda *a, **k: (_ for _ in ()).throw(RuntimeError("plot boom")),
    )
    result = IoTFieldAgent().execute(
        {"packets": _sim_packets(), "figures": True}
    )
    assert any("plot dashboard:" in w for w in result.warnings)


def test_figures_saved_to_output_dir(tmp_output):
    from pycsamt.agents import IoTFieldAgent

    result = IoTFieldAgent().execute(
        {
            "packets": _sim_packets(),
            "figures": True,
            "output_dir": str(tmp_output),
        }
    )
    assert result.get("figure_paths")


# ── LLM interpretation branches ──────────────────────────────────────────────


def test_llm_interpretation_none_without_api_key():
    from pycsamt.agents import IoTFieldAgent

    result = IoTFieldAgent().execute(
        {"packets": _sim_packets(), "figures": False}
    )
    assert result.llm_interpretation is None


def test_llm_interpretation_exception_is_recorded():
    from pycsamt.agents import IoTFieldAgent

    agent = IoTFieldAgent(api_key="fake-key")
    agent.query_llm = lambda *a, **k: (_ for _ in ()).throw(
        RuntimeError("llm boom")
    )
    result = agent.execute({"packets": _sim_packets(), "figures": False})
    assert result.llm_interpretation is None
    assert any("llm interpretation:" in w for w in result.warnings)
