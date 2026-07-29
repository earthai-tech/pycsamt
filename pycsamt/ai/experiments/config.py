# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Immutable source-of-truth configuration for AI inversion experiments.

This module describes *what* an experiment is and which evidence must pass.
It intentionally does not start training jobs, import machine-learning
frameworks, or construct Maxwell solvers.  Runtime orchestration and result
records belong to later modules in :mod:`pycsamt.ai.experiments`.
"""

from __future__ import annotations

import hashlib
import json
import operator
import re
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from types import MappingProxyType
from typing import Any

import numpy as np

from ..data.manifest import canonical_hash

__all__ = [
    "SeedPlan",
    "DatasetReference",
    "AcceptanceCriterion",
    "GateEvaluation",
    "ExperimentConfig",
]

_HEX_64 = re.compile(r"^[0-9a-f]{64}$")
_PORTABLE_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")
_STAGES = {
    "baseline",
    "data_audit",
    "geology",
    "domain_gap",
    "forward_2d",
    "learning_2d",
    "feasibility_3d",
    "forward_3d",
    "learning_3d",
    "hybrid",
    "field_evaluation",
}
_OPERATORS = {
    "<": operator.lt,
    "<=": operator.le,
    ">": operator.gt,
    ">=": operator.ge,
    "==": operator.eq,
}


def _freeze(value: Any) -> Any:
    if isinstance(value, Mapping):
        return MappingProxyType(
            {str(key): _freeze(item) for key, item in value.items()}
        )
    if isinstance(value, (list, tuple)):
        return tuple(_freeze(item) for item in value)
    return value


def _thaw(value: Any) -> Any:
    if isinstance(value, Mapping):
        return {str(key): _thaw(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_thaw(item) for item in value]
    return value


def _json_section(value: Mapping[str, Any], name: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise TypeError(f"{name} must be a mapping.")
    canonical_hash(value)
    plain = _thaw(value)
    return _freeze(json.loads(json.dumps(plain, allow_nan=False)))


def _digest(value: str | None, name: str, *, required: bool = False) -> str | None:
    if value is None:
        if required:
            raise ValueError(f"{name} is required.")
        return None
    result = str(value).lower().strip()
    if result.startswith("sha256:"):
        result = result[7:]
    if not _HEX_64.fullmatch(result):
        raise ValueError(f"{name} must be a 64-character SHA-256 digest.")
    return result


def _utc(value: str | None) -> str | None:
    if value is None:
        return None
    text = str(value).strip()
    try:
        parsed = datetime.fromisoformat(
            text[:-1] + "+00:00" if text.endswith("Z") else text
        )
    except ValueError as exc:
        raise ValueError("created_utc must be an ISO-8601 timestamp.") from exc
    if parsed.tzinfo is None or parsed.utcoffset() is None:
        raise ValueError("created_utc must include a timezone.")
    return (
        parsed.astimezone(timezone.utc)
        .isoformat(timespec="seconds")
        .replace("+00:00", "Z")
    )


@dataclass(frozen=True)
class SeedPlan:
    """Derive stable, labeled child seeds from one experiment root seed.

    Parameters
    ----------
    root_seed : int
        Non-negative root seed smaller than ``2**64``.
    namespace : str, default="pycsamt.ai"
        Domain separator preventing identical labels in unrelated projects
        from producing the same child sequence.

    Examples
    --------
    Child seeds are label-stable and request-order independent:

    >>> plan = SeedPlan(42, namespace="willy")
    >>> plan.derive("dataset") == plan.derive("dataset")
    True
    >>> plan.derive("dataset") != plan.derive("network")
    True
    """

    root_seed: int
    namespace: str = "pycsamt.ai"

    def __post_init__(self) -> None:
        if not isinstance(self.root_seed, (int, np.integer)) or isinstance(
            self.root_seed, bool
        ):
            raise TypeError("root_seed must be an integer.")
        seed = int(self.root_seed)
        if seed < 0 or seed >= 2**64:
            raise ValueError("root_seed must be in [0, 2**64).")
        namespace = str(self.namespace).strip()
        if not namespace:
            raise ValueError("namespace cannot be empty.")
        object.__setattr__(self, "root_seed", seed)
        object.__setattr__(self, "namespace", namespace)

    def derive(self, label: str) -> int:
        """Derive one unsigned 32-bit seed for a named subsystem.

        Parameters
        ----------
        label : str
            Stable descriptive label such as ``"geology"``, ``"noise"``, or
            ``"network/seed-0"``.

        Returns
        -------
        int
            Deterministic value in ``[0, 2**32)``.

        Raises
        ------
        ValueError
            If ``label`` is empty.

        Examples
        --------
        >>> child = SeedPlan(1).derive("training")
        >>> 0 <= child < 2**32
        True
        """
        name = str(label).strip()
        if not name:
            raise ValueError("label cannot be empty.")
        payload = f"{self.namespace}\0{self.root_seed}\0{name}".encode()
        return int.from_bytes(hashlib.sha256(payload).digest()[:4], "big")

    def derive_many(self, labels: Sequence[str]) -> dict[str, int]:
        """Derive seeds for multiple unique labels.

        Parameters
        ----------
        labels : sequence of str
            Unique non-empty subsystem labels.

        Returns
        -------
        dict
            Labels in caller order mapped to stable child seeds.

        Examples
        --------
        >>> sorted(SeedPlan(2).derive_many(["data", "model"]))
        ['data', 'model']
        """
        names = tuple(str(label).strip() for label in labels)
        if any(not name for name in names) or len(set(names)) != len(names):
            raise ValueError("labels must be non-empty and unique.")
        return {name: self.derive(name) for name in names}

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-serializable seed plan.

        Returns
        -------
        dict
            Schema version, root seed, and namespace.

        Examples
        --------
        >>> SeedPlan(3).to_dict()["root_seed"]
        3
        """
        return {
            "schema_version": 1,
            "root_seed": self.root_seed,
            "namespace": self.namespace,
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> SeedPlan:
        """Restore a validated seed plan.

        Parameters
        ----------
        data : mapping
            State returned by :meth:`to_dict`.

        Returns
        -------
        SeedPlan
            Immutable seed plan.

        Examples
        --------
        >>> SeedPlan.from_dict(SeedPlan(4).to_dict()) == SeedPlan(4)
        True
        """
        if data.get("schema_version", 1) != 1:
            raise ValueError("unsupported SeedPlan schema version.")
        return cls(data["root_seed"], data.get("namespace", "pycsamt.ai"))


@dataclass(frozen=True)
class DatasetReference:
    """Pin an experiment to exact dataset preparation artifacts.

    Parameters
    ----------
    dataset_id : str
        Portable identifier from the dataset manifest.
    manifest_hash, split_hash : str
        SHA-256 digests of the complete dataset manifest and realization split.
    normalization_hash : str or None, optional
        Digest of fitted normalization state. It may be absent for an audit or
        forward-only experiment that has not normalized data.
    manifest_uri : str or None, optional
        Informational local path or remote URI. Integrity relies on the hash,
        not on this location.

    Examples
    --------
    >>> reference = DatasetReference("willy-v1", "a" * 64, "b" * 64)
    >>> reference.dataset_id
    'willy-v1'
    """

    dataset_id: str
    manifest_hash: str
    split_hash: str
    normalization_hash: str | None = None
    manifest_uri: str | None = None

    def __post_init__(self) -> None:
        dataset_id = str(self.dataset_id).strip()
        if not _PORTABLE_ID.fullmatch(dataset_id):
            raise ValueError("dataset_id must be a portable identifier.")
        uri = None if self.manifest_uri is None else str(self.manifest_uri).strip()
        if self.manifest_uri is not None and not uri:
            raise ValueError("manifest_uri cannot be empty.")
        object.__setattr__(self, "dataset_id", dataset_id)
        object.__setattr__(
            self,
            "manifest_hash",
            _digest(self.manifest_hash, "manifest_hash", required=True),
        )
        object.__setattr__(
            self,
            "split_hash",
            _digest(self.split_hash, "split_hash", required=True),
        )
        object.__setattr__(
            self,
            "normalization_hash",
            _digest(self.normalization_hash, "normalization_hash"),
        )
        object.__setattr__(self, "manifest_uri", uri)

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-serializable artifact reference.

        Returns
        -------
        dict
            Schema version and pinned artifact identifiers.

        Examples
        --------
        >>> DatasetReference("d", "a" * 64, "b" * 64).to_dict()[
        ...     "schema_version"
        ... ]
        1
        """
        return {
            "schema_version": 1,
            "dataset_id": self.dataset_id,
            "manifest_hash": self.manifest_hash,
            "split_hash": self.split_hash,
            "normalization_hash": self.normalization_hash,
            "manifest_uri": self.manifest_uri,
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> DatasetReference:
        """Restore a validated dataset reference.

        Parameters
        ----------
        data : mapping
            State returned by :meth:`to_dict`.

        Returns
        -------
        DatasetReference
            Immutable pinned reference.

        Examples
        --------
        >>> ref = DatasetReference("d", "a" * 64, "b" * 64)
        >>> DatasetReference.from_dict(ref.to_dict()) == ref
        True
        """
        if data.get("schema_version", 1) != 1:
            raise ValueError("unsupported DatasetReference schema version.")
        return cls(
            data["dataset_id"],
            data["manifest_hash"],
            data["split_hash"],
            data.get("normalization_hash"),
            data.get("manifest_uri"),
        )


@dataclass(frozen=True)
class AcceptanceCriterion:
    """Predeclare one numerical condition required for an experiment gate.

    Parameters
    ----------
    metric : str
        Exact metric key, preferably namespaced, for example
        ``"test.impedance_nrms"``.
    operator : {"<", "<=", ">", ">=", "=="}
        Comparison applied as ``observed operator threshold``.
    threshold : float
        Finite value fixed before results are inspected.
    description : str, optional
        Human-readable scientific justification.

    Examples
    --------
    >>> criterion = AcceptanceCriterion("test.impedance_nrms", "<=", 2.0)
    >>> criterion.evaluate(1.7)
    True
    >>> criterion.evaluate(2.5)
    False
    """

    metric: str
    operator: str
    threshold: float
    description: str = ""

    def __post_init__(self) -> None:
        metric = str(self.metric).strip()
        if not metric:
            raise ValueError("metric cannot be empty.")
        if self.operator not in _OPERATORS:
            raise ValueError(f"operator must be one of {tuple(_OPERATORS)}.")
        threshold = float(self.threshold)
        if not np.isfinite(threshold):
            raise ValueError("threshold must be finite.")
        object.__setattr__(self, "metric", metric)
        object.__setattr__(self, "threshold", threshold)
        object.__setattr__(self, "description", str(self.description).strip())

    def evaluate(self, observed: float) -> bool:
        """Evaluate an observed metric against the frozen threshold.

        Parameters
        ----------
        observed : float
            Finite measured value.

        Returns
        -------
        bool
            Result of ``observed operator threshold``.

        Raises
        ------
        ValueError
            If ``observed`` is NaN or infinite.

        Examples
        --------
        >>> AcceptanceCriterion("coverage", ">=", 0.9).evaluate(0.95)
        True
        """
        value = float(observed)
        if not np.isfinite(value):
            raise ValueError("observed metric must be finite.")
        return bool(_OPERATORS[self.operator](value, self.threshold))

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-serializable criterion.

        Returns
        -------
        dict
            Metric, comparison, threshold, and description.

        Examples
        --------
        >>> AcceptanceCriterion("x", "<", 1).to_dict()["operator"]
        '<'
        """
        return {
            "schema_version": 1,
            "metric": self.metric,
            "operator": self.operator,
            "threshold": self.threshold,
            "description": self.description,
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> AcceptanceCriterion:
        """Restore a validated acceptance criterion.

        Parameters
        ----------
        data : mapping
            State returned by :meth:`to_dict`.

        Returns
        -------
        AcceptanceCriterion
            Immutable numerical gate condition.

        Examples
        --------
        >>> c = AcceptanceCriterion("x", ">", 0)
        >>> AcceptanceCriterion.from_dict(c.to_dict()) == c
        True
        """
        if data.get("schema_version", 1) != 1:
            raise ValueError("unsupported AcceptanceCriterion schema version.")
        return cls(
            data["metric"],
            data["operator"],
            data["threshold"],
            data.get("description", ""),
        )


@dataclass(frozen=True)
class GateEvaluation:
    """Immutable result of evaluating configured acceptance criteria.

    Parameters
    ----------
    passed : bool
        Whether every criterion passed and no required metric was missing.
    criteria : mapping of str to bool
        Per-metric pass/fail results.
    observed : mapping of str to float
        Finite observed values for metrics that were available.
    missing : sequence of str
        Required metric keys absent from the supplied result mapping.
    complete : bool, default=True
        Whether every configured criterion was evaluated. A partial status
        report can never pass the final gate.

    Examples
    --------
    >>> result = GateEvaluation(True, {"nrms": True}, {"nrms": 1.2}, ())
    >>> result.failed_metrics
    ()
    """

    passed: bool
    criteria: Mapping[str, bool]
    observed: Mapping[str, float]
    missing: tuple[str, ...] = ()
    complete: bool = True

    def __post_init__(self) -> None:
        criteria = {str(key): bool(value) for key, value in self.criteria.items()}
        observed = {str(key): float(value) for key, value in self.observed.items()}
        if any(not np.isfinite(value) for value in observed.values()):
            raise ValueError("observed gate metrics must be finite.")
        missing = tuple(str(value) for value in self.missing)
        if len(set(missing)) != len(missing) or any(not value for value in missing):
            raise ValueError("missing metric names must be non-empty and unique.")
        complete = bool(self.complete)
        expected_passed = complete and not missing and all(criteria.values())
        if bool(self.passed) != expected_passed:
            raise ValueError(
                "passed is inconsistent with criteria and missing metrics."
            )
        object.__setattr__(self, "passed", bool(self.passed))
        object.__setattr__(self, "criteria", MappingProxyType(criteria))
        object.__setattr__(self, "observed", MappingProxyType(observed))
        object.__setattr__(self, "missing", missing)
        object.__setattr__(self, "complete", complete)

    @property
    def failed_metrics(self) -> tuple[str, ...]:
        """Return metrics that failed their configured comparisons.

        Returns
        -------
        tuple of str
            Failed available metrics; missing metrics are reported separately.

        Examples
        --------
        >>> GateEvaluation(False, {"x": False}, {"x": 2.0}).failed_metrics
        ('x',)
        """
        return tuple(key for key, passed in self.criteria.items() if not passed)

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-serializable gate result.

        Returns
        -------
        dict
            Overall and per-metric outcomes.

        Examples
        --------
        >>> GateEvaluation(True, {"x": True}, {"x": 0.5}).to_dict()["passed"]
        True
        """
        return {
            "schema_version": 1,
            "passed": self.passed,
            "criteria": dict(self.criteria),
            "observed": dict(self.observed),
            "missing": list(self.missing),
            "complete": self.complete,
        }


@dataclass(frozen=True)
class ExperimentConfig:
    """Immutable source of truth for one reproducible inversion experiment.

    Parameters
    ----------
    experiment_id : str
        Portable unique experiment identifier.
    stage : str
        Roadmap stage such as ``"baseline"``, ``"forward_2d"``, or
        ``"field_evaluation"``.
    dataset : DatasetReference
        Exact dataset, split, and optional normalization artifacts.
    seeds : SeedPlan
        Root and namespace used for labeled child seeds.
    model, training, physics : mapping
        Finite JSON-compatible configuration sections. Their schemas are owned
        by later model/trainer/solver adapters; this class freezes and hashes
        them without importing optional dependencies.
    acceptance : sequence of AcceptanceCriterion, optional
        Criteria fixed before results are viewed. Metric names must be unique.
    description : str, optional
        Human-readable experiment objective.
    tags : sequence of str, optional
        Unique searchable labels.
    created_utc : str or None, optional
        Timezone-aware ISO-8601 timestamp normalized to UTC.
    schema_version : int, default=1
        Configuration schema version.

    Examples
    --------
    >>> dataset = DatasetReference("willy-v1", "a" * 64, "b" * 64)
    >>> config = ExperimentConfig(
    ...     "m0-baseline",
    ...     "baseline",
    ...     dataset,
    ...     SeedPlan(42, "willy"),
    ...     model={"architecture": "unet"},
    ...     training={"epochs": 100},
    ...     acceptance=[AcceptanceCriterion("test.nrms", "<=", 2.0)],
    ... )
    >>> len(config.config_hash)
    64
    >>> config.child_seed("network") == config.seeds.derive("network")
    True
    """

    experiment_id: str
    stage: str
    dataset: DatasetReference
    seeds: SeedPlan
    model: Mapping[str, Any]
    training: Mapping[str, Any]
    physics: Mapping[str, Any] = field(default_factory=dict)
    acceptance: tuple[AcceptanceCriterion, ...] = ()
    description: str = ""
    tags: tuple[str, ...] = ()
    created_utc: str | None = None
    schema_version: int = 1

    def __post_init__(self) -> None:
        if self.schema_version != 1:
            raise ValueError("unsupported ExperimentConfig schema version.")
        experiment_id = str(self.experiment_id).strip()
        if not _PORTABLE_ID.fullmatch(experiment_id):
            raise ValueError("experiment_id must be a portable identifier.")
        stage = str(self.stage).strip()
        if stage not in _STAGES:
            raise ValueError(f"stage must be one of {tuple(sorted(_STAGES))}.")
        if not isinstance(self.dataset, DatasetReference):
            raise TypeError("dataset must be a DatasetReference.")
        if not isinstance(self.seeds, SeedPlan):
            raise TypeError("seeds must be a SeedPlan.")
        acceptance = tuple(self.acceptance)
        if any(not isinstance(item, AcceptanceCriterion) for item in acceptance):
            raise TypeError("acceptance entries must be AcceptanceCriterion objects.")
        metrics = [item.metric for item in acceptance]
        if len(set(metrics)) != len(metrics):
            raise ValueError("acceptance metric names must be unique.")
        tags = tuple(str(tag).strip() for tag in self.tags)
        if any(not tag for tag in tags) or len(set(tags)) != len(tags):
            raise ValueError("tags must be non-empty and unique.")
        object.__setattr__(self, "experiment_id", experiment_id)
        object.__setattr__(self, "stage", stage)
        object.__setattr__(self, "model", _json_section(self.model, "model"))
        object.__setattr__(self, "training", _json_section(self.training, "training"))
        object.__setattr__(self, "physics", _json_section(self.physics, "physics"))
        object.__setattr__(self, "acceptance", acceptance)
        object.__setattr__(self, "description", str(self.description).strip())
        object.__setattr__(self, "tags", tags)
        object.__setattr__(self, "created_utc", _utc(self.created_utc))

    @property
    def config_hash(self) -> str:
        """Return the canonical digest of the complete configuration.

        Returns
        -------
        str
            SHA-256 digest covering every serialized field.

        Examples
        --------
        >>> d = DatasetReference("d", "a" * 64, "b" * 64)
        >>> len(
        ...     ExperimentConfig(
        ...         "e", "baseline", d, SeedPlan(0), {}, {}
        ...     ).config_hash
        ... )
        64
        """
        return canonical_hash(self.to_dict())

    def child_seed(self, label: str) -> int:
        """Derive a stable subsystem seed from this experiment.

        Parameters
        ----------
        label : str
            Subsystem label passed to :meth:`SeedPlan.derive`.

        Returns
        -------
        int
            Deterministic unsigned 32-bit seed.

        Examples
        --------
        >>> d = DatasetReference("d", "a" * 64, "b" * 64)
        >>> c = ExperimentConfig("e", "baseline", d, SeedPlan(0), {}, {})
        >>> c.child_seed("data") == c.child_seed("data")
        True
        """
        return self.seeds.derive(label)

    def evaluate_gate(
        self, metrics: Mapping[str, float], *, require_all: bool = True
    ) -> GateEvaluation:
        """Evaluate observed metrics against predeclared criteria.

        Parameters
        ----------
        metrics : mapping of str to float
            Observed metric values.
        require_all : bool, default=True
            Treat absent configured metrics as missing failures. When false,
            absent metrics are omitted; this is useful only for partial status
            reports and cannot prove the final gate passed unless all criteria
            were supplied.

        Returns
        -------
        GateEvaluation
            Immutable overall, per-metric, observed, and missing results.

        Examples
        --------
        >>> d = DatasetReference("d", "a" * 64, "b" * 64)
        >>> c = ExperimentConfig(
        ...     "e",
        ...     "baseline",
        ...     d,
        ...     SeedPlan(0),
        ...     {},
        ...     {},
        ...     acceptance=[AcceptanceCriterion("nrms", "<=", 2)],
        ... )
        >>> c.evaluate_gate({"nrms": 1.5}).passed
        True
        """
        supplied = {str(key): float(value) for key, value in metrics.items()}
        results = {}
        observed = {}
        missing = []
        for criterion in self.acceptance:
            if criterion.metric not in supplied:
                if require_all:
                    missing.append(criterion.metric)
                continue
            value = supplied[criterion.metric]
            results[criterion.metric] = criterion.evaluate(value)
            observed[criterion.metric] = value
        all_evaluated = len(results) == len(self.acceptance)
        passed = all_evaluated and not missing and all(results.values())
        return GateEvaluation(
            passed,
            results,
            observed,
            tuple(missing),
            complete=all_evaluated,
        )

    def to_dict(self) -> dict[str, Any]:
        """Return the complete JSON-compatible configuration.

        Returns
        -------
        dict
            Mutable schema-1 representation.

        Examples
        --------
        >>> d = DatasetReference("d", "a" * 64, "b" * 64)
        >>> ExperimentConfig(
        ...     "e", "baseline", d, SeedPlan(0), {}, {}
        ... ).to_dict()["stage"]
        'baseline'
        """
        return {
            "schema_version": self.schema_version,
            "experiment_id": self.experiment_id,
            "stage": self.stage,
            "dataset": self.dataset.to_dict(),
            "seeds": self.seeds.to_dict(),
            "model": _thaw(self.model),
            "training": _thaw(self.training),
            "physics": _thaw(self.physics),
            "acceptance": [item.to_dict() for item in self.acceptance],
            "description": self.description,
            "tags": list(self.tags),
            "created_utc": self.created_utc,
        }

    def write_json(self, path: str | Path, *, overwrite: bool = True) -> Path:
        """Write a deterministic UTF-8 JSON configuration file.

        Parameters
        ----------
        path : str or pathlib.Path
            Destination file.
        overwrite : bool, default=True
            Permit replacement of an existing file.

        Returns
        -------
        pathlib.Path
            Destination path.

        Examples
        --------
        >>> from tempfile import TemporaryDirectory
        >>> d = DatasetReference("d", "a" * 64, "b" * 64)
        >>> c = ExperimentConfig("e", "baseline", d, SeedPlan(0), {}, {})
        >>> with TemporaryDirectory() as directory:
        ...     path = c.write_json(Path(directory) / "experiment.json")
        ...     loaded = ExperimentConfig.read_json(path)
        >>> loaded.config_hash == c.config_hash
        True
        """
        target = Path(path)
        if target.exists() and not overwrite:
            raise FileExistsError(f"experiment configuration already exists: {target}")
        target.write_text(
            json.dumps(self.to_dict(), indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        return target

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> ExperimentConfig:
        """Restore and validate a serialized configuration.

        Parameters
        ----------
        data : mapping
            State returned by :meth:`to_dict`.

        Returns
        -------
        ExperimentConfig
            Immutable source-of-truth configuration.

        Examples
        --------
        >>> d = DatasetReference("d", "a" * 64, "b" * 64)
        >>> c = ExperimentConfig("e", "baseline", d, SeedPlan(0), {}, {})
        >>> ExperimentConfig.from_dict(
        ...     c.to_dict()
        ... ).config_hash == c.config_hash
        True
        """
        if data.get("schema_version", 1) != 1:
            raise ValueError("unsupported ExperimentConfig schema version.")
        return cls(
            experiment_id=data["experiment_id"],
            stage=data["stage"],
            dataset=DatasetReference.from_dict(data["dataset"]),
            seeds=SeedPlan.from_dict(data["seeds"]),
            model=data["model"],
            training=data["training"],
            physics=data.get("physics", {}),
            acceptance=tuple(
                AcceptanceCriterion.from_dict(item)
                for item in data.get("acceptance", [])
            ),
            description=data.get("description", ""),
            tags=tuple(data.get("tags", [])),
            created_utc=data.get("created_utc"),
            schema_version=data.get("schema_version", 1),
        )

    @classmethod
    def read_json(cls, path: str | Path) -> ExperimentConfig:
        """Read and validate a JSON experiment configuration.

        Parameters
        ----------
        path : str or pathlib.Path
            Existing UTF-8 JSON file.

        Returns
        -------
        ExperimentConfig
            Validated immutable configuration.

        Examples
        --------
        >>> from tempfile import TemporaryDirectory
        >>> d = DatasetReference("d", "a" * 64, "b" * 64)
        >>> c = ExperimentConfig("e", "baseline", d, SeedPlan(0), {}, {})
        >>> with TemporaryDirectory() as directory:
        ...     path = c.write_json(Path(directory) / "config.json")
        ...     loaded = ExperimentConfig.read_json(path)
        >>> loaded.experiment_id
        'e'
        """
        return cls.from_dict(json.loads(Path(path).read_text(encoding="utf-8")))
