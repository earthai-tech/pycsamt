# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Local settings persist without credentials and surface local failures."""
from dash import no_update

from pycsamt.app.agent_master._ids import IDs
from pycsamt.app.agent_master._providers import is_llm, requires_api_key
from pycsamt.app.agent_master.callbacks import settings


def callback(app, input_id, output_id):
    for key, entry in app.callback_map.items():
        if output_id in key and entry["inputs"][0]["id"] == input_id:
            fn = entry["callback"]
            return getattr(fn, "__wrapped__", fn)
    raise AssertionError("Callback missing")


def test_local_provider_is_not_a_credential_provider():
    assert is_llm("ollama")
    assert not requires_api_key("ollama")
    assert requires_api_key("claude")


def test_local_settings_round_trip(agent_app, tmp_path, monkeypatch):
    monkeypatch.setattr(settings, "_CFG_DIR", tmp_path)
    monkeypatch.setattr(settings, "_CFG_FILE", tmp_path / "settings.json")
    save = callback(agent_app, IDs.BTN_SAVE_KEYS, IDs.STORE_SETTINGS)
    cfg, _, _ = save(1, "ollama", "ignored", None, "png", "", "", {},
                     "http://127.0.0.1:11434", "custom-coder:small", 60, 8192, 512, 0.1, 3)
    assert cfg["model_ollama"] == "custom-coder:small"
    assert "key_ollama" not in cfg
    assert settings._load_cfg() == cfg
    panel = callback(agent_app, IDs.ACTIVE_PROVIDER, IDs.LOCAL_PANEL)
    assert panel("ollama")[0] == {"display": "block"}
    assert panel("ollama")[2] == "custom-coder:small"
    assert panel("offline")[0] == {"display": "none"}
    bad, _, _ = save(1, "ollama", "", None, "png", "", "", {}, "https://remote.example")
    assert bad is no_update
    assert settings._load_cfg() == cfg


def test_local_check_error_is_visible(agent_app, monkeypatch):
    from pycsamt.agents import _local

    def unavailable(*args, **kwargs):
        raise _local.LocalModelError("Start ollama serve")

    monkeypatch.setattr(_local, "model_status", unavailable)
    check = callback(agent_app, IDs.LOCAL_CHECK, IDs.LOCAL_STATUS)
    assert "Start ollama serve" in check(1, "http://127.0.0.1:11434", "coder")


def test_cancelled_job_cannot_be_overwritten():
    from pycsamt.app.agent_master.callbacks import chat

    jid = chat._new_job()
    try:
        chat._update_job(jid, status="cancelled")
        chat._update_job(jid, status="done", result="late answer")
        assert chat._get_job(jid)["status"] == "cancelled"
        assert chat._get_job(jid).get("result") != "late answer"
    finally:
        with chat._JOBS_LOCK:
            chat._JOBS.pop(jid, None)
