# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Request fidelity, clarification, and edits of actual prior artifacts."""
import ast
from pathlib import Path

import pytest

from pycsamt.agents._generation import GenerationInput, exact_artifact_edit
from pycsamt.agents.code_gen import CodeGenerationAgent
from pycsamt.agents.router import IntentRouter
from pycsamt.api.agents import AGENT_CONFIG, _TLS


@pytest.fixture(autouse=True)
def enable_inference(monkeypatch):
    monkeypatch.setattr(_TLS, "force_offline", False, raising=False)


def test_full_request_and_unusual_constraints_reach_model(tmp_path, monkeypatch):
    request = "Compare stations S01 and S03; three panels; phase in degrees; filter frequency < 10 Hz; save result.svg; keep my unusual blue dashed legend."
    contract = GenerationInput.from_chat(request, [{"role": "user", "content": "Keep original station order."}],
        workflow_config={"workflow": "qc", "data_path": "data/3edis"},
        project_context={"selected_lines": ["L18PLT", "L22PLT"]},
        retrieved_evidence="verified API evidence")
    captured = {}
    def model(self, prompt, **kwargs):
        captured["prompt"] = prompt
        return "fig.savefig('result.svg')\n"
    monkeypatch.setattr(CodeGenerationAgent, "query_llm", model)
    result = CodeGenerationAgent(llm_provider="ollama").execute({"generation_input": contract, "output_dir": str(tmp_path)})
    assert result.status == "success"
    for constraint in ["S01", "S03", "three panels", "phase in degrees", "frequency < 10 Hz", "result.svg", "blue dashed legend", "original station order", "L18PLT", "L22PLT"]:
        assert constraint in captured["prompt"]
    assert "build_qc_table" in captured["prompt"]
    assert result.data["generation"]["requirements_verified"] is False
    assert result.data["generation"]["executed"] is False


def test_custom_composition_does_not_receive_qc_template(tmp_path, monkeypatch):
    captured = []
    monkeypatch.setattr(CodeGenerationAgent, "query_llm", lambda self, prompt, **kw: captured.append(prompt) or "print('custom')")
    request = GenerationInput(original_request="Write code calling pycsamt.emtools._core.ensure_sites for two independent datasets.",
                              workflow_config={"workflow": "custom"}, retrieved_evidence="ensure_sites loads Sites")
    result = CodeGenerationAgent(llm_provider="ollama").execute({"generation_input": request, "output_dir": str(tmp_path)})
    assert result.status == "success"
    assert "build_qc_table" not in captured[0]
    assert "Optional template:" not in captured[0]


@pytest.mark.parametrize("text,fragment", [
    ("Write me a script", "What should"),
    ("Write code comparing L18PLT and L22PLT with three panels saved as line_comparison.png", "Which quantities"),
])
def test_missing_task_or_panel_science_is_clarified_before_writing(tmp_path, text, fragment):
    with AGENT_CONFIG.offline():
        result = CodeGenerationAgent().execute({"generation_input": GenerationInput(original_request=text), "output_dir": str(tmp_path)})
    assert result.status == "needs_review"
    assert fragment in result.data["clarification"]
    assert not list(tmp_path.iterdir())


def test_output_directory_edit_preserves_input_filename_and_previous_file(tmp_path):
    code = "import os\nsites = ensure_sites('data/3edis')\nos.makedirs('results/qc', exist_ok=True)\nqc_table.to_csv('results/qc/qc_results.csv')\n"
    prior = tmp_path / "workflow_script.py"
    prior.write_text(code)
    history = [{"role": "assistant", "code": code, "generation": {"output_dir": "results/qc"}}]
    contract = GenerationInput.from_chat("Change only the output directory to results/qc_review.", history)
    with AGENT_CONFIG.offline():
        assert IntentRouter().route(contract.original_request, history=history).intent == "code"
        result = CodeGenerationAgent().execute({"generation_input": contract, "output_dir": str(tmp_path)})
    assert result.status == "success"
    assert result.data["code"] == code.replace("'results/qc", "'results/qc_review")
    assert prior.read_text() == code
    assert Path(result.data["script_path"]) != prior
    assert result.data["generation"]["output_dir"] == "results/qc_review"


def test_dpi_edit_preserves_every_other_byte_including_unicode():
    code = "# Résistivité\nfig, ax = plt.subplots(3)\nfig.savefig('station.png', dpi=100)\nprint('unchanged')\n"
    contract = GenerationInput(original_request="Keep the same figure but save it at 300 dpi instead.", previous_code=code)
    edited, _ = exact_artifact_edit(contract)
    assert edited == code.replace("dpi=100", "dpi=300")
    ast.parse(edited)


def test_general_edit_passes_prior_code_instead_of_rebuilding_template(tmp_path, monkeypatch):
    previous = "factors = estimate_ss_ama(sites)\nfactors.to_csv('factors.csv')\ncorrected = correct_ss_ama(sites)\n"
    requested = "Remove the correction step but keep the factor-estimation table."
    contract = GenerationInput.from_chat(requested, [{"role": "assistant", "code": previous}])
    expected = "factors = estimate_ss_ama(sites)\nfactors.to_csv('factors.csv')\n"
    def model(self, prompt, **kwargs):
        assert "previous_code" in prompt and "corrected = correct_ss_ama" in prompt
        assert requested in prompt
        assert "Optional template:" not in prompt
        return expected
    monkeypatch.setattr(CodeGenerationAgent, "query_llm", model)
    result = CodeGenerationAgent(llm_provider="ollama").execute({"generation_input": contract, "output_dir": str(tmp_path)})
    assert result.data["code"] == expected.strip()


def test_explanation_of_prior_code_is_not_routed_to_edit():
    history = [{"role": "assistant", "code": "print(1)"}]
    with AGENT_CONFIG.offline():
        assert IntentRouter().route("Explain the third step without changing the code.", history=history).intent != "code"


def test_model_clarification_and_invalid_python_never_write(tmp_path, monkeypatch):
    agent = CodeGenerationAgent(llm_provider="ollama")
    contract = GenerationInput(original_request="Write a QC script", workflow_config={"workflow": "qc"})
    monkeypatch.setattr(agent, "query_llm", lambda *a, **k: "CLARIFY: Which stations should be selected?")
    result = agent.execute({"generation_input": contract, "output_dir": str(tmp_path)})
    assert result.data["clarification"] == "Which stations should be selected?"
    monkeypatch.setattr(agent, "query_llm", lambda *a, **k: "Here is your script: <broken>")
    assert agent.execute({"generation_input": contract, "output_dir": str(tmp_path)}).status == "failed"
    assert not list(tmp_path.iterdir())


def test_offline_detailed_template_does_not_claim_constraints_applied(tmp_path):
    request = GenerationInput(original_request="Write QC code selecting S03 only and save table.csv", workflow_config={"workflow": "qc"})
    with AGENT_CONFIG.offline():
        result = CodeGenerationAgent().execute({"generation_input": request, "output_dir": str(tmp_path)})
    assert result.status == "needs_review"
    assert result.data["code"]
    assert "placeholder" in result.data["code"]
    assert any("not been applied" in w for w in result.warnings)


def test_clarification_answer_keeps_original_requirements():
    original = "Compare L18PLT and L22PLT with three panels saved as comparison.png"
    history = [{"role": "assistant", "pending_request": original, "content": "Which quantities?"}]
    contract = GenerationInput.from_chat("Resistivity, phase and skew.", history)
    assert original in contract.task_text
    with AGENT_CONFIG.offline():
        assert IntentRouter().route(contract.original_request, history=history).intent == "code"


@pytest.mark.parametrize("code", [
    "fig.savefig(f'results/qc/{station}.png')",
    "open('results/qc/table.csv', 'w')",
    "os.makedirs('results/qc'); open('results/qc/table.csv', 'w')",
])
def test_unrecognized_output_expressions_require_general_edit(code):
    contract = GenerationInput(original_request="Change only the output directory to results/review.", previous_code=code, previous_output_dir="results/qc")
    assert exact_artifact_edit(contract) is None


@pytest.mark.parametrize("query,code", [
    ("Run QC and save qc_results.csv", "print('table')"),
    ("Estimate static shift without applying correction", "correct_ss_ama(sites)"),
])
def test_detectable_requirement_omissions_are_unsaved_drafts(tmp_path, monkeypatch, query, code):
    monkeypatch.setattr(CodeGenerationAgent, "query_llm", lambda *a, **k: code)
    contract = GenerationInput(original_request=query, workflow_config={"workflow": "qc"})
    result = CodeGenerationAgent(llm_provider="ollama").execute({"generation_input": contract, "output_dir": str(tmp_path)})
    assert result.status == "needs_review"
    assert "no script was saved" in result.data["review_reason"]
    assert not list(tmp_path.iterdir())


def test_explicit_example_routes_as_code_not_survey_metrics():
    from pycsamt.agents.router import classify_intent_offline

    assert classify_intent_offline("Write a 1-D forward-model example at 20 frequencies from 0.01 to 100 Hz.")[0] == "code"


def test_request_parser_does_not_treat_table_prose_as_directory():
    from pycsamt.agents.context import _regex_extract

    cfg = _regex_extract("Write a script to load data/3edis, run QC, and save the table to qc_results.csv in an output folder.")
    assert cfg["data_path"] == "data/3edis"
    assert "output_dir" not in cfg
    assert _regex_extract("run qc and save to results/qc")["output_dir"] == "results/qc"


def test_forward_values_and_units_survive_into_template(tmp_path):
    from pycsamt.agents._generation import forward_parameters

    text = "Write a 1-D forward-model example for resistivities 100, 10, 1000 ohm m and finite-layer thicknesses 500 and 1000 m, at 20 log-spaced frequencies from 0.01 to 100 Hz."
    assert forward_parameters(text) == {"resistivity": [100., 10., 1000.], "thickness": [500., 1000.], "frequency_grid": [0.01, 100., 20]}
    contract = GenerationInput(original_request=text, workflow_config={"workflow": "forward"})
    with AGENT_CONFIG.offline():
        result = CodeGenerationAgent().execute({"generation_input": contract, "output_dir": str(tmp_path)})
    code = result.data["code"]
    assert "resistivity=[100.0, 10.0, 1000.0]" in code
    assert "thickness=[500.0, 1000.0]" in code
    assert "np.geomspace(0.01, 100.0, 20)" in code
    assert "ensure_sites" not in code


def test_forward_parser_does_not_silently_treat_kilometres_as_metres():
    from pycsamt.agents._generation import forward_parameters

    parsed = forward_parameters("resistivities 100, 10 kohm m and thicknesses 1 and 2 km at 20 log-spaced frequencies from 1 to 10 kHz")
    assert parsed == {}


def test_local_compaction_retains_request_and_previous_script():
    contract = GenerationInput(original_request="Keep all 17 explicitly requested stations.", previous_code="print('prior script')", retrieved_evidence="reference " * 2000,
                               api_evidence=[{"symbol": str(i), "signature": "function(sites)", "doc": "docs" * 500} for i in range(4)])
    contract.compact_evidence()
    prompt = contract.prompt()
    assert contract.original_request in prompt
    assert "prior script" in prompt
    assert len(contract.retrieved_evidence) == 400
    assert len(contract.api_evidence) == 2
    assert contract.evidence_notes


def test_estimation_only_template_has_no_correction_call(tmp_path):
    contract = GenerationInput(original_request="Estimate static-shift factors and export before applying any correction", workflow_config={"workflow": "static_shift", "data_path": "data/3edis"})
    with AGENT_CONFIG.offline():
        result = CodeGenerationAgent().execute({"generation_input": contract, "output_dir": str(tmp_path)})
    code = result.data["code"]
    assert "estimate_ss_ama(sites)" in code
    assert "ss_table.empty" in code
    assert "ss_table.to_csv" in code
    assert "correct_ss_ama" not in code


def test_forward_literal_constraints():
    from pycsamt.agents._generation import GenerationInput, missing_request_constraints
    request = GenerationInput(workflow_config={"workflow": "forward", "frequency_grid": [0.01, 100, 20]})
    bad = "frequencies = np.logspace(0, 2, 20)\nmodel = LayeredModel(resistivities=[100], thicknesses=[])"
    assert len(missing_request_constraints(request, bad)) == 2
    assert not missing_request_constraints(request, "frequencies = np.logspace(-2, 2, 20)\nmodel = LayeredModel(resistivity=[100], thickness=[])")
