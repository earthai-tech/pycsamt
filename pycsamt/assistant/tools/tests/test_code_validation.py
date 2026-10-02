# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Static validator contracts, skipped execution and scientific literals."""

import json
import time
from pathlib import Path
from types import SimpleNamespace

import pytest

from pycsamt.agents._code_validation import validate_and_repair
from pycsamt.agents._generation import GenerationInput
from pycsamt.agents._local import LocalSettings, local_session
from pycsamt.assistant.tools.validation_tools import validate_generated_code


@pytest.fixture
def source_root(tmp_path):
    package = tmp_path / "pycsamt"
    package.mkdir()
    (package / "__init__.py").write_text("", encoding="utf-8")
    (package / "fixture.py").write_text(
        'raise RuntimeError("Never import this module")\ndef function(required, /, *, option=1):\n    return required\nclass Solver:\n    def __init__(self, *, frequency):\n        pass\n    def run(self, model):\n        return model\n',
        encoding="utf-8",
    )
    (package / "dynamic.py").write_text(
        "def __getattr__(name):\n    return object()\n", encoding="utf-8"
    )
    (package / "broken.py").write_text("def broken(:\n", encoding="utf-8")
    return tmp_path


@pytest.mark.parametrize(
    "code",
    [
        "import pycsamt.no_such_module",
        "from pycsamt.fixture import absent",
        "from pycsamt.fixture import function as f\nf(required=1)",
        "from pycsamt.fixture import function as f\nf(1, bad=2)",
        "from pycsamt.fixture import Solver\ns = Solver(frequency=[1])\ns.run()",
    ],
)
def test_missing_modules_symbols_and_supported_arguments_fail(
    source_root, code
):
    result = validate_generated_code(code, root=source_root)
    assert not result["ok"]
    assert result["errors"]
    assert result["executed"] is False


def test_static_inspection_never_imports_target(source_root):
    result = validate_generated_code(
        "from pycsamt.fixture import function as f\nf(1, option=2)",
        root=source_root,
    )
    assert result["checks"]["imports"]["state"] == "passed"
    assert result["checks"]["arguments"]["state"] == "passed"
    assert result["api_evidence"][0]["path"] == "pycsamt/fixture.py"
    assert "sha256" in result["api_evidence"][0]


@pytest.mark.parametrize(
    "code",
    [
        "from pycsamt.broken import unknown",
        "from pycsamt.dynamic import unknown",
    ],
)
def test_unverifiable_import_not_reported_as_checked(source_root, code):
    result = validate_generated_code(code, root=source_root)
    assert result["checks"]["imports"]["state"] == "unverifiable"
    assert not result["checked"]
    assert result["warnings"]


def test_dynamic_kwargs_and_external_dependencies_explicitly_unverified(
    source_root,
):
    result = validate_generated_code(
        "import unavailable_external\nfrom pycsamt.fixture import function\nfunction(**settings)",
        root=source_root,
    )
    assert result["checks"]["arguments"]["state"] == "unverifiable"
    assert result["checks"]["imports"]["state"] == "unverifiable"
    assert result[
        "ok"
    ]  # compatibility: absence of detected static failures only


@pytest.mark.parametrize(
    "rho,thickness,valid",
    [
        ([100, 10, 1000], [500, 1000], True),
        ([100, 0], [500], False),
        ([100, 10], [500, 1000], False),
    ],
)
def test_layered_model_literal_invariants(rho, thickness, valid):
    code = f"from pycsamt.forward.synthetic import LayeredModel\nmodel = LayeredModel(resistivity={rho!r}, thickness={thickness!r})"
    result = validate_generated_code(code)
    assert result["ok"] is valid
    assert result["checks"]["scientific"]["state"] == (
        "passed" if valid else "failed"
    )


def test_literal_correction_factors_and_obsolete_constructor_keywords():
    result = validate_generated_code(
        "from pycsamt.emtools.ss import apply_ss_factors\napply_ss_factors(sites, {'S01': -1})"
    )
    assert result["checks"]["scientific"]["state"] == "failed"
    result = validate_generated_code(
        "from pycsamt.forward.synthetic import LayeredModel\nm = LayeredModel(resistivities=[100], thicknesses=[])"
    )
    assert result["checks"]["arguments"]["state"] == "failed"


def test_requested_execution_fails_closed(tmp_path):
    marker = tmp_path / "must_not_exist"
    code = f"open({str(marker)!r}, 'w').write('executed')"
    result = validate_generated_code(code, execute_fixture=True)
    assert result["checks"]["execution"]["state"] == "unverifiable"
    assert result["checks"]["artifacts"]["state"] == "not_checked"
    assert not marker.exists()


def test_repair_uses_errors_and_evidence_and_revalidates():
    prompts = []
    agent = SimpleNamespace(
        llm_provider="ollama",
        llm_available=True,
        query_llm=lambda prompt, **kw: (
            prompts.append((prompt, kw))
            or "from pycsamt.forward.synthetic import LayeredModel\nm = LayeredModel(resistivity=[100], thickness=[])"
        ),
    )
    bad = "from pycsamt.forward.synthetic import LayeredModel\nm = LayeredModel(resistivities=[100], thicknesses=[])"
    code, result = validate_and_repair(
        agent,
        bad,
        GenerationInput(original_request="Use one layer at 100 ohm-m"),
        started=time.monotonic(),
    )
    assert result["ok"]
    assert len(prompts) == 1
    assert "constructor fields" in prompts[0][0]
    assert "Use one layer" in prompts[0][0]
    assert prompts[0][1]["max_tokens"] == 1024
    assert result["repair"]["attempts"][0]["errors_after"] == []


def test_repair_hard_stop_and_no_cloud_calls():
    prompts = []

    def broken(prompt, **kwargs):
        prompts.append(prompt)
        return f"import pycsamt.no_such_module_{len(prompts)}"

    agent = SimpleNamespace(
        llm_provider="ollama", llm_available=True, query_llm=broken
    )
    _, result = validate_and_repair(
        agent,
        "import pycsamt.no_such_module",
        None,
        started=time.monotonic(),
        max_repairs=99,
    )
    assert len(prompts) == 2
    assert not result["ok"]
    assert result["repair"]["stop_reason"] == "repair attempt limit reached"
    agent.llm_provider = "claude"
    _, result = validate_and_repair(
        agent, "import pycsamt.no_such_module", None, started=time.monotonic()
    )
    assert len(prompts) == 2
    assert "shared request deadline" in result["repair"]["stop_reason"]


def test_expired_budget_and_cancellation_make_no_model_call():
    agent = SimpleNamespace(
        llm_provider="ollama",
        llm_available=True,
        query_llm=lambda *a, **k: pytest.fail("model must not be called"),
    )
    _, result = validate_and_repair(
        agent, "bad syntax :", None, started=time.monotonic() - 61
    )
    assert "deadline" in result["repair"]["stop_reason"]
    with local_session(LocalSettings(), cancelled=lambda: True):
        _, result = validate_and_repair(
            agent, "bad syntax :", None, started=time.monotonic()
        )
    assert "stopped" in result["repair"]["stop_reason"]


def test_publication_preserves_prior_script_and_saves_report(
    tmp_path, monkeypatch
):
    from pycsamt.agents.code_gen import CodeGenerationAgent
    from pycsamt.api.agents import _TLS

    monkeypatch.setattr(_TLS, "force_offline", False, raising=False)
    agent = CodeGenerationAgent(llm_provider="ollama")
    monkeypatch.setattr(
        CodeGenerationAgent, "llm_available", property(lambda self: True)
    )
    contract = GenerationInput(
        original_request="Write a custom script printing a message",
        workflow_config={"workflow": "custom"},
        retrieved_evidence="Python print() displays a message.",
    )
    prior = tmp_path / "workflow_script.py"
    prior.write_text("prior = 'keep'\n", encoding="utf-8")
    monkeypatch.setattr(
        agent, "query_llm", lambda *a, **k: "import pycsamt.no_such_module"
    )
    failed = agent.execute(
        {"generation_input": contract, "output_dir": str(tmp_path)}
    )
    assert failed.status == "needs_review"
    assert failed.data["script_path"] is None
    assert prior.read_text() == "prior = 'keep'\n"
    monkeypatch.setattr(
        agent, "query_llm", lambda *a, **k: "print('review this')"
    )
    success = agent.execute(
        {"generation_input": contract, "output_dir": str(tmp_path)}
    )
    assert Path(success.data["script_path"]).name == "workflow_script_2.py"
    saved = json.loads(Path(success.data["validation_path"]).read_text())
    assert saved["code_sha256"] == success.data["validation"]["code_sha256"]
    assert saved["executed"] is False


def test_shadowed_import_and_runtime_imports_are_not_certified(source_root):
    code = "from pycsamt.fixture import function\nfunction = lambda: 1\nfunction()"
    result = validate_generated_code(code, root=source_root)
    assert result["ok"]
    assert result["checks"]["arguments"]["state"] == "unverifiable"
    assert result["checks"]["runtime_imports"]["state"] == "not_checked"
