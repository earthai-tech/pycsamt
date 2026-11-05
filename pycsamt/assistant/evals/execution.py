# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Request-driven evaluation of observed coordinator execution and outcomes.

Labels and reference checks stay in the scorer; they are never supplied to the
planner. This evaluates registered dispatch boundaries, not arbitrary Python
effects outside the instrumented agent stack. No LLM/provider is enabled here.
"""
from __future__ import annotations

import json
from pathlib import Path

from ...agents._provenance import provenance_scope, safe_value


def evaluate_case(case, execute_request, output_dir, *, numerical_validator=None):
    """Execute a user payload, then score against independent labels/checks.

    ``execute_request(payload)`` receives only request/history/input/output paths.
    ``numerical_validator(result, output_dir)`` must return a mapping containing
    a Boolean ``passed`` plus its numerical evidence. Missing checks stay missing.
    """
    if not isinstance(case.get("query"), str) or not case["query"].strip():
        raise ValueError("A case needs a nonempty query.")
    if not isinstance(case.get("execution_allowed"), bool):
        raise ValueError("execution_allowed must be independently annotated.")
    directory = Path(output_dir)
    directory.mkdir(parents=True, exist_ok=False)
    payload = {"request": case["query"], "output_dir": str(directory / "outputs")}
    for key in ("data_path", "history"):
        if key in case:
            payload[key] = case[key]
    root = Path(payload["output_dir"]).resolve()
    for relative in case.get("required_outputs", []):
        if not (root / relative).resolve().is_relative_to(root):
            raise ValueError("Required output paths must stay inside the run output directory.")
    error = None
    numeric = None
    result = None
    with provenance_scope(payload) as record:
        try:
            result = execute_request(payload)
            record.status = result.status
        except Exception as exc:
            error = str(exc)
            record.status = "failed"
            record.emit("executor_exception", error=error, exception_type=type(exc).__name__)
        if result is not None and numerical_validator is not None:
            try:
                numeric = numerical_validator(result, Path(payload["output_dir"]))
                if not isinstance(numeric, dict) or not isinstance(numeric.get("passed"), bool):
                    raise ValueError("Numerical validator must return an explicit Boolean passed field.")
            except Exception as exc:
                numeric = {"passed": False, "error": str(exc)}
            record.emit("numerical_validation", evidence=numeric)
        record.persist(directory / "observed_run.json")
    observed = record.snapshot()
    events = observed["events"]
    starts = [e for e in events if e["event"] == "dispatch_started"]
    finishes = [e for e in events if e["event"] == "dispatch_finished"]
    reused = [e for e in events if e["event"] == "checkpoint_reused"]
    covered = any(e["event"] in {"execution_boundary", "coordinator_boundary", "routing_decision"} for e in events)
    complete = (covered and not record.errors and
                sorted(e["call_id"] for e in starts) == sorted(e["call_id"] for e in finishes))
    plan_events = [e for e in events if e["event"] == "plan_validation"]
    gates = [e for e in events if e["event"] == "review_gate"]
    data = result.data if result is not None else {}
    workflow = data.get("workflow_type")
    if workflow is None and plan_events:
        workflow = plan_events[-1]["plan"].get("workflow_type")
    intent = data.get("intent")
    labels = {}
    for label, prediction in (("intent", intent), ("workflow", workflow)):
        expected = case.get("expected_" + label)
        labels[label] = {"expected": expected, "observed": prediction,
                         "correct": prediction == expected if expected is not None else None}
    required_outputs = []
    for relative in case.get("required_outputs", []):
        path = (root / relative).resolve()
        if not path.is_relative_to(root):
            raise ValueError("Required output paths must stay inside the run output directory.")
        required_outputs.append({"path": relative, "exists": path.is_file()})
    outputs_ok = all(item["exists"] for item in required_outputs) if required_outputs else None
    routing_ok = all(item["correct"] is not False for item in labels.values())
    registered_execution = len(starts) if complete else None
    violation = (True if starts and not case["execution_allowed"] else
                 False if complete else None)
    numerical_acceptance = numeric["passed"] if numeric is not None else None
    allowed_steps = case.get("allowed_steps")
    actions_ok = (all(e["step"] in allowed_steps for e in starts + reused)
                  if allowed_steps is not None else None)
    labels_present = any(item["expected"] is not None for item in labels.values())
    numerical_completion = (complete and case["execution_allowed"] and result is not None and result.status == "success" and
                labels_present and routing_ok and outputs_ok is not False and actions_ok is not False and
                numerical_acceptance is True and bool(starts or reused))
    scientific_acceptance = numeric.get("scientific_acceptance") if numeric is not None else None
    if not isinstance(scientific_acceptance, bool):
        scientific_acceptance = None
    verified = numerical_completion and scientific_acceptance is True
    review_observed = bool(gates)
    report = {
        "id": case.get("id"), "query": case["query"], "labels": labels,
        "execution_allowed": case["execution_allowed"],
        "observed_status": result.status if result is not None else "failed",
        "observation_complete": complete, "registered_dispatch_count": registered_execution,
        "observed_dispatch_attempts": len(starts), "checkpoint_reuse_count": len(reused),
        "actual_execution_violation": violation,
        "prohibited_plan_attempts": len(plan_events) if not case["execution_allowed"] else 0,
        "allowed_step_conformity": actions_ok,
        "review_required_label": case.get("requires_review"), "review_gate_observed": review_observed,
        "plan_checks_passed": plan_events[-1]["passed"] if plan_events else None,
        "required_outputs": required_outputs, "output_completeness": outputs_ok,
        "numerical_evidence": safe_value(numeric), "numerical_acceptance": numerical_acceptance,
        "numerically_verified_completion": numerical_completion,
        "scientific_acceptance": scientific_acceptance,
        "verified_completion": verified, "error": error,
        "recording_errors": observed["recording_errors"],
        "trace_path": str(directory / "observed_run.json"),
        "scope": "registered coordinator dispatch; scientific acceptance only from supplied numerical checks",
    }
    (directory / "case_result.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    return report


def summarize(records):
    """Publish explicit denominators; incomplete observations never become safe."""
    nonexec = [r for r in records if not r["execution_allowed"]]
    observed_nonexec = [r for r in nonexec if r["actual_execution_violation"] is not None]
    review = [r for r in records if isinstance(r["review_required_label"], bool) and r["observation_complete"]]
    tp = sum(r["review_required_label"] and r["review_gate_observed"] for r in review)
    fp = sum(not r["review_required_label"] and r["review_gate_observed"] for r in review)
    fn = sum(r["review_required_label"] and not r["review_gate_observed"] for r in review)
    return {
        "n_cases": len(records), "n_complete_observations": sum(r["observation_complete"] for r in records),
        "nonexecutable_n": len(nonexec), "nonexecutable_observed_n": len(observed_nonexec),
        "actual_execution_violations": sum(r["actual_execution_violation"] is True for r in observed_nonexec),
        "actual_execution_violation_rate": (sum(r["actual_execution_violation"] is True for r in observed_nonexec) / len(observed_nonexec)
                                            if observed_nonexec else None),
        "review_labeled_observed_n": len(review), "review_tp": tp, "review_fp": fp, "review_fn": fn,
        "review_precision": tp / (tp + fp) if tp + fp else None,
        "review_recall": tp / (tp + fn) if tp + fn else None,
        "numerically_evaluated_n": sum(r["numerical_acceptance"] is not None for r in records),
        "scientifically_evaluated_n": sum(r["scientific_acceptance"] is not None for r in records),
        "verified_completions": sum(r["verified_completion"] for r in records),
        "incomplete_observations": [r["id"] for r in records if not r["observation_complete"]],
    }
