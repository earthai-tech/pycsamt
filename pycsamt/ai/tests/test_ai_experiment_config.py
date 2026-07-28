"""Contracts for :mod:`pycsamt.ai.experiments.config`."""

from __future__ import annotations

import json

import pytest

from pycsamt.ai.experiments import (
    AcceptanceCriterion,
    DatasetReference,
    ExperimentConfig,
    GateEvaluation,
    SeedPlan,
)


def _dataset() -> DatasetReference:
    return DatasetReference(
        "willy-v1",
        "a" * 64,
        "b" * 64,
        normalization_hash="c" * 64,
        manifest_uri="artifacts/willy/manifest.json",
    )


def _config() -> ExperimentConfig:
    return ExperimentConfig(
        experiment_id="m0-baseline",
        stage="baseline",
        dataset=_dataset(),
        seeds=SeedPlan(42, "willy"),
        model={"architecture": "unet", "channels": [16, 32]},
        training={"epochs": 100, "learning_rate": 1e-3},
        physics={"backend": "mt1d", "frequency_order": "descending"},
        acceptance=(
            AcceptanceCriterion("test.nrms", "<=", 2.0),
            AcceptanceCriterion("test.coverage", ">=", 0.9),
        ),
        description="Frozen baseline",
        tags=("willy", "baseline"),
        created_utc="2026-07-26T02:00:00+02:00",
    )


def test_seed_plan_is_stable_order_independent_and_serializable():
    plan = SeedPlan(42, "willy")
    assert plan.derive("data") == plan.derive("data")
    assert plan.derive("data") != plan.derive("network")
    assert plan.derive_many(["data", "network"]) == {
        "data": plan.derive("data"),
        "network": plan.derive("network"),
    }
    assert SeedPlan.from_dict(plan.to_dict()) == plan
    with pytest.raises(ValueError, match="unique"):
        plan.derive_many(["same", "same"])


def test_dataset_reference_validates_and_roundtrips():
    reference = _dataset()
    assert DatasetReference.from_dict(reference.to_dict()) == reference
    assert reference.normalization_hash == "c" * 64
    with pytest.raises(ValueError, match="manifest_hash"):
        DatasetReference("dataset", "bad", "b" * 64)
    with pytest.raises(ValueError, match="portable"):
        DatasetReference("bad id", "a" * 64, "b" * 64)


def test_acceptance_criterion_and_gate_results():
    criterion = AcceptanceCriterion("nrms", "<=", 2.0, "Response fit")
    assert criterion.evaluate(2.0)
    assert not criterion.evaluate(2.1)
    assert AcceptanceCriterion.from_dict(criterion.to_dict()) == criterion
    with pytest.raises(ValueError, match="finite"):
        criterion.evaluate(float("nan"))

    passed = GateEvaluation(True, {"nrms": True}, {"nrms": 1.0})
    assert passed.passed
    assert passed.failed_metrics == ()
    failed = GateEvaluation(False, {"nrms": False}, {"nrms": 3.0})
    assert failed.failed_metrics == ("nrms",)
    with pytest.raises(ValueError, match="inconsistent"):
        GateEvaluation(True, {"nrms": False}, {"nrms": 3.0})


def test_experiment_config_is_deeply_immutable_and_hash_stable():
    config = _config()
    assert config.created_utc == "2026-07-26T00:00:00Z"
    assert config.model["channels"] == (16, 32)
    with pytest.raises(TypeError):
        config.model["architecture"] = "gcn"
    with pytest.raises(TypeError):
        config.training["nested"] = {}
    restored = ExperimentConfig.from_dict(json.loads(json.dumps(config.to_dict())))
    assert restored.config_hash == config.config_hash
    assert restored.child_seed("network") == config.child_seed("network")


def test_experiment_gate_pass_fail_missing_and_partial_status():
    config = _config()
    passed = config.evaluate_gate({"test.nrms": 1.5, "test.coverage": 0.95})
    assert passed.passed and passed.complete
    failed = config.evaluate_gate({"test.nrms": 3.0, "test.coverage": 0.95})
    assert not failed.passed
    assert failed.failed_metrics == ("test.nrms",)
    missing = config.evaluate_gate({"test.nrms": 1.5})
    assert not missing.passed and not missing.complete
    assert missing.missing == ("test.coverage",)
    partial = config.evaluate_gate({"test.nrms": 1.5}, require_all=False)
    assert not partial.passed and not partial.complete
    assert partial.missing == ()


def test_experiment_json_roundtrip_and_overwrite_guard(tmp_path):
    config = _config()
    path = config.write_json(tmp_path / "experiment.json")
    restored = ExperimentConfig.read_json(path)
    assert restored.to_dict() == config.to_dict()
    with pytest.raises(FileExistsError):
        config.write_json(path, overwrite=False)


def test_experiment_rejects_invalid_stage_duplicate_metrics_and_nonfinite_sections():
    with pytest.raises(ValueError, match="stage must"):
        ExperimentConfig("e", "unknown", _dataset(), SeedPlan(0), {}, {})
    with pytest.raises(ValueError, match="metric names"):
        ExperimentConfig(
            "e",
            "baseline",
            _dataset(),
            SeedPlan(0),
            {},
            {},
            acceptance=(
                AcceptanceCriterion("same", "<", 1),
                AcceptanceCriterion("same", ">", 0),
            ),
        )
    with pytest.raises(ValueError, match="finite JSON"):
        ExperimentConfig(
            "e",
            "baseline",
            _dataset(),
            SeedPlan(0),
            {"bad": float("nan")},
            {},
        )
