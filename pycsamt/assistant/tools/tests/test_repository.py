# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Developer scope must be current, bounded and independent of science RAG."""

import pytest

from pycsamt.assistant.tools.repository import (
    RepositoryTools,
    evidence_text,
    is_developer_question,
)


@pytest.fixture
def repository(tmp_path):
    app = tmp_path / "pycsamt/app"
    app.mkdir(parents=True)
    (app / "chat.py").write_text(
        'raise RuntimeError("never import")\ndef dispatch_code(request, *, output="results"):\n    """Preserve the original request."""\n    return request\n',
        encoding="utf-8",
    )
    tests = app / "tests"
    tests.mkdir()
    (tests / "test_chat.py").write_text(
        "def test_dispatch_code():\n    assert True\n", encoding="utf-8"
    )
    return RepositoryTools(tmp_path)


def test_static_search_signature_and_freshness(repository):
    result = repository.search("Where is dispatch_code implemented?")
    assert result["sources"][0]["path"] == "pycsamt/app/chat.py"
    card = repository.inspect_symbol("pycsamt/app/chat.py", "dispatch_code")
    assert "output='results'" in card["signature"]
    assert card["line"] == 2
    path = repository.root / card["path"]
    path.write_text(
        "def changed_api(new_parameter):\n    return new_parameter\n",
        encoding="utf-8",
    )
    refreshed = repository.search("changed_api")
    assert refreshed["sources"][0]["sha256"] != card["sha256"]
    assert not repository.search("dispatch_code")["sources"]
    assert "No matching developer evidence" in evidence_text(
        repository.search("nonexistent_api")
    )


def test_tests_are_opt_in_and_labelled(repository):
    assert all(
        "tests/" not in c["path"]
        for c in repository.search("dispatch_code")["sources"]
    )
    result = repository.search("test_dispatch_code", include_tests=True)
    assert "expected behavior" in result["sources"][0]["kind"]
    with pytest.raises(ValueError):
        repository.read_source("pycsamt/app/tests/test_chat.py")


@pytest.mark.parametrize(
    "path",
    [
        "../outside.py",
        ".env",
        "pycsamt/app/.env.py",
        "pycsamt/app/credentials.py",
        "pycsamt/app/__pycache__/cache.py",
        "data/private.py",
        "pycsamt/app/generated/output.py",
    ],
)
def test_excluded_paths(repository, path):
    with pytest.raises(ValueError):
        repository.read_source(path, include_tests=True)


def test_credentials_redacted_and_excerpt_bounded(repository):
    path = repository.root / "pycsamt/app/sample.py"
    path.write_text(
        'API_KEY = "sensitive-value"\n' + "print('fine')\n" * 300,
        encoding="utf-8",
    )
    card = repository.read_source("pycsamt/app/sample.py", lines=99999)
    assert "sensitive-value" not in card["excerpt"]
    assert len(card["excerpt"].splitlines()) <= 80
    assert len(card["excerpt"]) <= 4000


def test_provenance_warns_on_checkout_mismatch(repository):
    info = repository.provenance()
    assert not info["checkout_matches_loaded_package"]
    assert "may differ" in info["limitation"]
    assert "SHA-256" in info["freshness"]


def test_redaction_preserves_signature_and_line_numbers(repository):
    path = repository.root / "pycsamt/app/auth.py"
    path.write_text(
        'api_key = """secret\nsecond line"""\ndef configure(api_key=None, password="hidden"):\n    return api_key\n',
        encoding="utf-8",
    )
    card = repository.inspect_symbol("pycsamt/app/auth.py", "configure")
    assert card["line"] == 3
    assert "api_key=None" in card["signature"]
    assert "hidden" not in card["signature"]
    assert (
        "secret"
        not in repository.read_source("pycsamt/app/auth.py")["excerpt"]
    )


@pytest.mark.parametrize(
    "query",
    [
        "Where is _dispatch_code implemented?",
        "How does agent_master cancel a job?",
        "Explain pycsamt/app/agent_master/callbacks/chat.py",
        "Show the unit tests for dispatch_code",
    ],
)
def test_developer_routes(query):
    from pycsamt.agents.router import classify_intent_offline

    assert is_developer_question(query)
    assert classify_intent_offline(query)[0] == "question"


@pytest.mark.parametrize(
    "query",
    [
        "How do I load EDI files?",
        "What does StaticShiftAgent do?",
        "Run QC on data/3edis",
        "Write a script to load EDI files",
        "What is apparent resistivity?",
    ],
)
def test_science_not_redirected(query):
    assert not is_developer_question(query)


def test_agent_bypasses_science_rag_and_reports_missing(
    monkeypatch, repository
):
    from pycsamt.agents.package_qa import PackageQAAgent
    from pycsamt.api.agents import AGENT_CONFIG

    monkeypatch.setattr(
        "pycsamt.assistant.tools.repository.RepositoryTools",
        lambda: repository,
    )
    monkeypatch.setattr(
        PackageQAAgent,
        "_build_rag",
        lambda *a, **k: pytest.fail("science retrieval must stay separate"),
    )
    with AGENT_CONFIG.offline():
        agent = PackageQAAgent()
        result = agent.execute(
            {
                "question": "Where is dispatch_code implemented?",
                "scope": "developer",
            }
        )
        assert "pycsamt/app/chat.py:2" in result.data["answer"]
        assert result.data["source"] == "developer_offline"
        missing = agent.execute(
            {
                "question": "Where is nonexistent_api implemented?",
                "scope": "developer",
            }
        )
        assert "No matching developer evidence" in missing.data["answer"]


def test_source_instruction_is_data_in_model_prompt(monkeypatch, repository):
    from pycsamt.agents.package_qa import PackageQAAgent

    monkeypatch.setattr(
        "pycsamt.assistant.tools.repository.RepositoryTools",
        lambda: repository,
    )
    monkeypatch.setattr(
        PackageQAAgent, "llm_available", property(lambda self: True)
    )
    seen = {}

    def query(self, prompt, **kwargs):
        seen.update(kwargs)
        seen["prompt"] = prompt
        return "The request passes through dispatch_code [1]."

    monkeypatch.setattr(PackageQAAgent, "query_llm", query)
    result = PackageQAAgent(llm_provider="ollama").execute(
        {
            "question": "Where is dispatch_code implemented?",
            "scope": "developer",
        }
    )
    assert "untrusted reference data" in seen["system_message"]
    assert "Source references:" in result.data["answer"]
    assert "pycsamt/app/chat.py:2" in result.data["answer"]


def test_developer_read_question_bypasses_loaded_data_shortcut():
    from pycsamt.app.agent_master.callbacks.chat import _looks_like_data_read

    assert not _looks_like_data_read("How does agent_master read EDI files?")
    assert _looks_like_data_read("read the EDI data")
