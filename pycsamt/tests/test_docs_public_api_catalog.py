"""Regression tests for the generated public API long table."""

from __future__ import annotations

import importlib
import importlib.util
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]


def _load(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _catalog_module():
    return _load(
        REPO_ROOT / "docs/source/_ext/public_api_catalog.py",
        "test_public_api_catalog_extension",
    )


def _validator_module():
    return _load(
        REPO_ROOT / "docs/scripts/validate_api_catalog.py",
        "test_validate_api_catalog",
    )


def test_reexported_objects_keep_public_facade_paths():
    members = _catalog_module()._members("pycsamt.ai.data")

    assert "pycsamt.ai.data.SurveyData" in members["Classes"]
    assert "pycsamt.ai.data.merge_surveys" in members["Functions"]
    assert "pycsamt.ai.data.canonical_hash" in members["Functions"]
    assert "pycsamt.ai.data.contracts.merge_surveys" not in members["Functions"]
    assert "pycsamt.ai.data.manifest.canonical_hash" not in members["Functions"]


def test_cross_package_class_aliases_use_the_documented_owner():
    members = _catalog_module()._members("pycsamt.core")

    assert "pycsamt.transformers.jedi.AVGtoEDI" in members["Classes"]
    assert "pycsamt.transformers.jedi.JtoEDI" in members["Classes"]
    assert "pycsamt.core.AVGtoEDI" not in members["Classes"]

    inversion_members = _catalog_module()._members("pycsamt.inversion")
    map_members = _catalog_module()._members("pycsamt.map")
    assert "pycsamt.inversion.model.StartingModel" in inversion_members["Classes"]
    assert "pycsamt.map.volume.Map3D" in map_members["Classes"]


def test_undocumented_exceptions_and_tuple_factories_are_not_class_rows():
    core_members = _catalog_module()._members("pycsamt.core")
    maxwell_members = _catalog_module()._members("pycsamt.forward.maxwell")

    assert "pycsamt.core.RegistryError" not in core_members["Classes"]
    assert "pycsamt.core.Packer" not in core_members["Classes"]
    assert not any(
        name.endswith("Error") for name in maxwell_members["Classes"]
    )


def test_lazy_agent_exports_are_visible_and_use_facade_paths():
    agents = importlib.import_module("pycsamt.agents")
    assert "AnomalyDetectionAgent" in dir(agents)
    assert "BatchSurveyAgent" in dir(agents)

    members = _catalog_module()._members("pycsamt.agents")
    assert "pycsamt.agents.AnomalyDetectionAgent" in members["Classes"]
    assert "pycsamt.agents.BatchSurveyAgent" in members["Classes"]
    assert not any(
        name.startswith("pycsamt.agents.anomaly_agent.")
        for name in members["Classes"]
    )


def test_validator_reads_the_catalogue_package_list():
    validator = _validator_module()
    packages = validator.catalogue_packages()

    assert packages[0] == "pycsamt.api"
    assert "pycsamt.ai.data" in packages
    assert "pycsamt.agents" in packages
    assert packages[-1] == "pycsamt.app"


def test_validator_reports_only_unindexed_objects(monkeypatch):
    validator = _validator_module()
    monkeypatch.setattr(
        validator, "catalogue_packages", lambda: ["pycsamt.example"]
    )
    monkeypatch.setattr(
        validator,
        "catalogue_objects",
        lambda _packages: ["pycsamt.example.ClassA", "pycsamt.example.func_b"],
    )
    monkeypatch.setattr(
        validator,
        "inventory_objects",
        lambda _path: {"pycsamt.example.ClassA"},
    )

    assert validator.unresolved_objects(Path("objects.inv")) == [
        "pycsamt.example.func_b"
    ]
