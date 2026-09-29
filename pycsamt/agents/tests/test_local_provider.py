# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Local inference contract, privacy, budgets, and interrupted transport."""
import json
import threading
import time
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

import pytest

from pycsamt.agents import _local
from pycsamt.agents._base import BaseAgent
from pycsamt.api.agents import AGENT_CONFIG, _TLS


class Agent(BaseAgent):
    def execute(self, data):
        return self.query_llm(data["prompt"])


@pytest.fixture(autouse=True)
def inference_enabled(monkeypatch):
    monkeypatch.setattr(_TLS, "force_offline", False, raising=False)


@pytest.fixture
def server():
    state = {"requests": [], "delay": 0, "error": None, "remote": False, "length": False}

    class Handler(BaseHTTPRequestHandler):
        def log_message(self, *args):
            pass

        def do_POST(self):
            body = json.loads(self.rfile.read(int(self.headers["Content-Length"])))
            state["requests"].append((self.path, body))
            if state["error"]:
                self.send_response(state["error"])
                self.end_headers()
                self.wfile.write(b'{"error":"model not found"}')
                return
            if self.path == "/api/show":
                result = {"remote_host": "https://example.org"} if state["remote"] else {"details": {"quantization_level": "Q4_K_M"}}
                self.send_response(200)
                self.end_headers()
                self.wfile.write(json.dumps(result).encode())
                return
            time.sleep(state["delay"])
            try:
                self.send_response(200)
                self.end_headers()
                self.wfile.write(b'{"message":{"content":"grounded answer"},"done":false}\n')
                end = {"done": True, "prompt_eval_count": 15, "eval_count": 2,
                       "done_reason": "length" if state["length"] else "stop"}
                self.wfile.write(json.dumps(end).encode() + b"\n")
            except (BrokenPipeError, ConnectionResetError, ConnectionAbortedError):
                pass

    httpd = ThreadingHTTPServer(("127.0.0.1", 0), Handler)
    httpd.daemon_threads = True
    thread = threading.Thread(target=httpd.serve_forever, daemon=True)
    thread.start()
    yield f"http://127.0.0.1:{httpd.server_port}", state
    httpd.shutdown()
    httpd.server_close()


@pytest.mark.parametrize("endpoint", ["https://example.org", "http://10.0.0.1:11434",
                                     "http://localhost:11434/path", "http://user@localhost",
                                     "http://localhost?redirect=remote"])
def test_rejects_nonlocal_endpoints(endpoint):
    with pytest.raises(ValueError):
        _local.LocalSettings(endpoint=endpoint)


@pytest.mark.parametrize("kwargs", [{"timeout": 0}, {"temperature": float("nan")},
                                    {"model": "model:cloud"}, {"output_tokens": 8192}])
def test_settings_limits(kwargs):
    with pytest.raises(ValueError):
        _local.LocalSettings(**kwargs)


def test_real_http_key_free_usage_and_limits(server, monkeypatch):
    endpoint, state = server
    monkeypatch.setenv("OPENAI_API_KEY", "must-not-be-used")
    with _local.local_session(_local.LocalSettings(endpoint=endpoint, output_tokens=50, max_calls=1, temperature=0.7)) as req:
        agent = Agent("local", llm_provider="ollama")
        assert agent.api_key is None and agent.llm_available
        assert agent.query_llm("Explain loading", max_tokens=100) == "grounded answer"
        assert agent.last_usage["eval_count"] == 2
        assert agent._last_cost == 0
        assert req.calls == 1
        with pytest.raises(_local.LocalModelError, match="budget"):
            agent.query_llm("again")
    assert [p for p, _ in state["requests"]] == ["/api/show", "/api/chat"]
    assert state["requests"][-1][1]["options"]["num_predict"] == 50
    assert state["requests"][-1][1]["options"]["temperature"] == 0.7
    assert not _local.local_only()


def test_offline_disables_local_and_explicit_cloud_keys(monkeypatch):
    with AGENT_CONFIG.offline():
        local = Agent("local", llm_provider="ollama")
        cloud = Agent("cloud", api_key="explicit")
        assert local.query_llm("test") is None
        assert cloud.query_llm("test") is None
    assert not local.llm_available


def test_local_session_prevents_inherited_cloud_provider(server):
    endpoint, _ = server
    with _local.local_session(_local.LocalSettings(endpoint=endpoint)):
        child = Agent("nested", llm_provider="openai", api_key="cloud-key")
        assert child.llm_provider == "ollama"
        assert child.api_key is None


def test_cancel_during_model_load(server):
    endpoint, state = server
    state["delay"] = 3
    cancel = threading.Event()
    timer = threading.Timer(0.2, cancel.set)
    timer.start()
    start = time.monotonic()
    try:
        with _local.local_session(_local.LocalSettings(endpoint=endpoint), cancel.is_set):
            with pytest.raises(_local.LocalCancelled):
                Agent("local").query_llm("test")
    finally:
        timer.cancel()
    assert time.monotonic() - start < 1.5


def test_total_timeout_includes_model_load(server):
    endpoint, state = server
    state["delay"] = 3
    with _local.local_session(_local.LocalSettings(endpoint=endpoint, timeout=1)):
        with pytest.raises(_local.LocalModelError, match="timed out"):
            Agent("local").query_llm("test")


@pytest.mark.parametrize("flag,message", [("remote", "remotely"), ("length", "output limit")])
def test_remote_models_and_truncation_are_errors(server, flag, message):
    endpoint, state = server
    state[flag] = True
    with _local.local_session(_local.LocalSettings(endpoint=endpoint)):
        with pytest.raises(_local.LocalModelError, match=message):
            Agent("local").query_llm("test")


def test_missing_model_error(server):
    endpoint, state = server
    state["error"] = 404
    with _local.local_session(_local.LocalSettings(endpoint=endpoint)):
        with pytest.raises(_local.LocalModelError, match="not found"):
            Agent("local").query_llm("test")


def test_context_overflow_does_not_generate(server):
    endpoint, state = server
    with _local.local_session(_local.LocalSettings(endpoint=endpoint, context_tokens=512, output_tokens=10)):
        with pytest.raises(_local.LocalModelError, match="context budget"):
            Agent("local").query_llm("x" * 1000)
    assert all(p != "/api/chat" for p, _ in state["requests"])


def test_local_retrieval_cannot_reuse_dense_cache(monkeypatch):
    from pycsamt.assistant.rag import retriever as module
    from pycsamt.assistant.rag import index_store
    from pycsamt.assistant.rag.embeddings import resolve_embedding_backend

    module._CACHE.clear()
    monkeypatch.setattr(index_store, "load_index", lambda **kw: [object()])
    monkeypatch.setattr(index_store, "index_is_stale", lambda **kw: False)
    monkeypatch.setattr(module, "_resolve_dense", lambda *a, **k: pytest.fail("cloud embedding resolution"))
    monkeypatch.setattr(module, "Retriever", lambda corpus, **kw: kw)
    module._CACHE["fake-root"] = "cloud cached retriever"
    try:
        with _local.local_session(_local.LocalSettings()):
            result = module.build_retriever("fake-root", use_feedback=False, embed_api_key="cloud-key")
            assert result["embed_backend"] is None
            assert resolve_embedding_backend(api_key="cloud-key") is None
    finally:
        module._CACHE.clear()
