# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Deterministic, lineage-aware geological-realization dataset splits.

Stations, frequency patches, noise variants, and augmented copies derived from
one geological model must never be divided across training and evaluation
partitions.  This module therefore splits unique realization IDs and can bind
related IDs through explicit lineage labels.
"""

from __future__ import annotations

import hashlib
import json
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field, replace
from types import MappingProxyType
from typing import Any

import numpy as np

__all__ = [
    "RealizationSplit",
    "split_realizations",
    "realization_folds",
]

_PARTITIONS = ("train", "validation", "test")


def _ids(values: Sequence[str], label: str) -> tuple[str, ...]:
    result = tuple(str(value).strip() for value in values)
    if any(not value for value in result):
        raise ValueError(f"{label} cannot contain empty IDs.")
    if len(set(result)) != len(result):
        raise ValueError(f"{label} contains duplicate IDs.")
    return result


def _seed(value: int | None) -> int | None:
    if value is None:
        return None
    if not isinstance(value, (int, np.integer)) or isinstance(value, bool):
        raise TypeError("seed must be an integer or None.")
    if value < 0:
        raise ValueError("seed must be non-negative or None.")
    return int(value)


def _fraction(value: float, name: str) -> float:
    result = float(value)
    if not np.isfinite(result) or result < 0 or result >= 1:
        raise ValueError(f"{name} must be finite and in [0, 1).")
    return result


def _lineage_map(
    ids: Sequence[str], lineage: Mapping[str, str] | None
) -> Mapping[str, str]:
    if lineage is None:
        return MappingProxyType({})
    normalized = {
        str(realization_id).strip(): str(parent).strip()
        for realization_id, parent in lineage.items()
    }
    if set(normalized) != set(ids):
        missing = sorted(set(ids) - set(normalized))
        extra = sorted(set(normalized) - set(ids))
        raise ValueError(
            f"lineage keys must exactly match split IDs; missing={missing}, extra={extra}."
        )
    if any(not parent for parent in normalized.values()):
        raise ValueError("lineage labels cannot be empty.")
    return MappingProxyType(normalized)


def _groups(
    ids: Sequence[str], lineage: Mapping[str, str] | None
) -> dict[str, tuple[str, ...]]:
    if lineage is None:
        return {realization_id: (realization_id,) for realization_id in ids}
    result: dict[str, list[str]] = {}
    for realization_id in sorted(ids):
        result.setdefault(str(lineage[realization_id]), []).append(
            realization_id
        )
    return {name: tuple(members) for name, members in result.items()}


def _target_count(size: int, fraction: float, available: int) -> int:
    target = int(np.floor(size * fraction))
    if fraction > 0 and target == 0 and available >= 3:
        target = 1
    return target


def _take_groups(
    available: list[str], sizes: Mapping[str, int], target: int, reserve: int
) -> tuple[list[str], list[str]]:
    if target <= 0:
        return [], available
    chosen: list[str] = []
    total = 0
    remaining = list(available)
    while len(remaining) > reserve:
        candidates = remaining
        best = min(
            candidates,
            key=lambda name: abs(target - (total + sizes[name])),
        )
        previous_distance = abs(target - total)
        next_distance = abs(target - (total + sizes[best]))
        if chosen and next_distance >= previous_distance:
            break
        chosen.append(best)
        total += sizes[best]
        remaining.remove(best)
    if not chosen and len(remaining) > reserve:
        chosen.append(remaining.pop(0))
    return chosen, remaining


@dataclass(frozen=True)
class RealizationSplit:
    """Immutable train/validation/test assignment of geological realizations.

    Parameters
    ----------
    train, validation, test : sequence of str
        Unique realization identifiers. The three partitions must be disjoint,
        and training cannot be empty.
    seed : int or None, optional
        Random seed used to create the assignment. ``None`` records that the
        seed is unknown or intentionally nondeterministic.
    lineage : mapping, optional
        Complete mapping from every realization ID to its parent geological
        lineage. IDs sharing a lineage must occupy the same partition.
    strategy : str, default="random"
        Human-readable splitting strategy identifier.

    Examples
    --------
    >>> split = RealizationSplit(("r1", "r2"), ("r3",), ("r4",), seed=7)
    >>> split.sizes
    {'train': 2, 'validation': 1, 'test': 1}
    >>> split.partition("r3")
    'validation'
    """

    train: tuple[str, ...]
    validation: tuple[str, ...]
    test: tuple[str, ...]
    seed: int | None = None
    lineage: Mapping[str, str] = field(default_factory=dict)
    strategy: str = "random"

    def __post_init__(self) -> None:
        train = _ids(self.train, "train")
        validation = _ids(self.validation, "validation")
        test = _ids(self.test, "test")
        if (
            set(train) & set(validation)
            or set(train) & set(test)
            or set(validation) & set(test)
        ):
            raise ValueError(
                "train, validation, and test realization IDs must be disjoint."
            )
        if not train:
            raise ValueError("the training split cannot be empty.")
        seed = _seed(self.seed)
        strategy = str(self.strategy).strip()
        if not strategy:
            raise ValueError("strategy cannot be empty.")
        all_ids = train + validation + test
        lineage = _lineage_map(all_ids, self.lineage if self.lineage else None)
        object.__setattr__(self, "train", train)
        object.__setattr__(self, "validation", validation)
        object.__setattr__(self, "test", test)
        object.__setattr__(self, "seed", seed)
        object.__setattr__(self, "lineage", lineage)
        object.__setattr__(self, "strategy", strategy)
        self.assert_no_lineage_leakage()

    @property
    def all_ids(self) -> tuple[str, ...]:
        """Return all realization IDs in partition order.

        Returns
        -------
        tuple of str
            Training IDs followed by validation IDs and test IDs.

        Examples
        --------
        >>> RealizationSplit(("a",), ("b",), ("c",)).all_ids
        ('a', 'b', 'c')
        """
        return self.train + self.validation + self.test

    @property
    def sizes(self) -> dict[str, int]:
        """Return the number of realizations in each partition.

        Returns
        -------
        dict
            ``train``, ``validation``, and ``test`` counts.

        Examples
        --------
        >>> RealizationSplit(("a", "b"), (), ("c",)).sizes["train"]
        2
        """
        return {name: len(getattr(self, name)) for name in _PARTITIONS}

    @property
    def fractions(self) -> dict[str, float]:
        """Return realized partition fractions.

        Returns
        -------
        dict
            Counts divided by total realization count. These may differ from
            requested targets when lineages contain multiple realizations.

        Examples
        --------
        >>> RealizationSplit(("a", "b"), ("c",), ("d",)).fractions
        {'train': 0.5, 'validation': 0.25, 'test': 0.25}
        """
        total = len(self.all_ids)
        return {name: len(getattr(self, name)) / total for name in _PARTITIONS}

    @property
    def split_hash(self) -> str:
        """Return a deterministic SHA-256 digest of the assignment.

        Returns
        -------
        str
            Digest covering partitions, seed, lineage, and strategy.

        Examples
        --------
        >>> len(RealizationSplit(("a",), (), ()).split_hash)
        64
        """
        payload = json.dumps(
            self.to_dict(), sort_keys=True, separators=(",", ":")
        )
        return hashlib.sha256(payload.encode("utf-8")).hexdigest()

    def partition(self, realization_id: str) -> str:
        """Return the partition containing a realization.

        Parameters
        ----------
        realization_id : str
            Exact realization identifier.

        Returns
        -------
        {"train", "validation", "test"}
            Partition name.

        Raises
        ------
        KeyError
            If the identifier is unknown.

        Examples
        --------
        >>> RealizationSplit(("a",), ("b",), ()).partition("b")
        'validation'
        """
        key = str(realization_id)
        for name in _PARTITIONS:
            if key in getattr(self, name):
                return name
        raise KeyError(f"unknown realization ID {key!r}.")

    def ids_for(self, partition: str) -> tuple[str, ...]:
        """Return IDs assigned to a named partition.

        Parameters
        ----------
        partition : {"train", "validation", "test"}
            Partition to retrieve.

        Returns
        -------
        tuple of str
            Immutable realization IDs.

        Raises
        ------
        ValueError
            If ``partition`` is unsupported.

        Examples
        --------
        >>> RealizationSplit(("a",), ("b",), ()).ids_for("train")
        ('a',)
        """
        if partition not in _PARTITIONS:
            raise ValueError(f"partition must be one of {_PARTITIONS}.")
        return getattr(self, partition)

    def mask(
        self,
        realization_ids: Sequence[str],
        partition: str,
        *,
        unknown: str = "raise",
    ) -> np.ndarray:
        """Build a Boolean sample mask from realization IDs.

        Parameters
        ----------
        realization_ids : sequence of str
            Per-sample realization IDs; duplicates are allowed because many
            stations or noise variants may belong to one realization.
        partition : {"train", "validation", "test"}
            Partition selected as ``True``.
        unknown : {"raise", "false"}, default="raise"
            Raise for IDs absent from the split or mark them false.

        Returns
        -------
        ndarray of bool
            One value per supplied sample ID.

        Examples
        --------
        >>> split = RealizationSplit(("a",), (), ("b",))
        >>> split.mask(["a", "b", "a"], "train").tolist()
        [True, False, True]
        """
        selected = set(self.ids_for(partition))
        if unknown not in {"raise", "false"}:
            raise ValueError("unknown must be 'raise' or 'false'.")
        values = tuple(str(value) for value in realization_ids)
        if unknown == "raise":
            missing = sorted(set(values) - set(self.all_ids))
            if missing:
                raise KeyError(f"unknown realization IDs: {missing}.")
        result = np.asarray(
            [value in selected for value in values], dtype=bool
        )
        result.setflags(write=False)
        return result

    def assert_complete(self, expected_ids: Sequence[str]) -> None:
        """Assert exact coverage of an expected realization collection.

        Parameters
        ----------
        expected_ids : sequence of str
            Unique IDs expected across all partitions.

        Returns
        -------
        None
            Successful return proves no expected ID is missing or unexpected.

        Raises
        ------
        ValueError
            If expected IDs are duplicated or coverage differs.

        Examples
        --------
        >>> split = RealizationSplit(("a",), (), ("b",))
        >>> split.assert_complete(["b", "a"]) is None
        True
        """
        expected = _ids(expected_ids, "expected_ids")
        missing = sorted(set(expected) - set(self.all_ids))
        extra = sorted(set(self.all_ids) - set(expected))
        if missing or extra:
            raise ValueError(
                f"split coverage differs; missing={missing}, extra={extra}."
            )

    def assert_no_lineage_leakage(
        self, lineage: Mapping[str, str] | None = None
    ) -> None:
        """Assert that every parent lineage occurs in one partition only.

        Parameters
        ----------
        lineage : mapping or None, optional
            Complete external ID-to-lineage mapping. When omitted, use the
            mapping persisted on the split. If neither exists, each realization
            is already an independent unit and the check succeeds.

        Returns
        -------
        None
            Successful return proves no supplied lineage crosses partitions.

        Raises
        ------
        ValueError
            If mapping coverage is incomplete or a lineage leaks.

        Examples
        --------
        >>> split = RealizationSplit(("a1", "a2"), (), ("b1",))
        >>> split.assert_no_lineage_leakage(
        ...     {"a1": "a", "a2": "a", "b1": "b"}
        ... ) is None
        True
        """
        mapping = (
            self.lineage
            if lineage is None
            else _lineage_map(self.all_ids, lineage)
        )
        if not mapping:
            return
        owners: dict[str, str] = {}
        for realization_id in self.all_ids:
            parent = mapping[realization_id]
            partition = self.partition(realization_id)
            previous = owners.setdefault(parent, partition)
            if previous != partition:
                raise ValueError(
                    f"lineage {parent!r} crosses partitions {previous!r} and {partition!r}."
                )

    def reassign(
        self, realization_ids: Sequence[str], partition: str
    ) -> RealizationSplit:
        """Return a copy with complete realizations moved to one partition.

        Parameters
        ----------
        realization_ids : sequence of str
            Known IDs to move. When lineage is recorded, all members of each
            affected lineage must be supplied together.
        partition : {"train", "validation", "test"}
            Destination partition.

        Returns
        -------
        RealizationSplit
            New validated assignment. Existing relative order is retained;
            moved IDs are appended in the requested order.

        Raises
        ------
        KeyError
            If an ID is unknown.
        ValueError
            If a lineage is moved partially or training would become empty.

        Examples
        --------
        >>> split = RealizationSplit(("a", "b"), (), ("c",))
        >>> split.reassign(["b"], "validation").validation
        ('b',)
        """
        target = self.ids_for(partition)
        moving = _ids(realization_ids, "realization_ids")
        unknown = sorted(set(moving) - set(self.all_ids))
        if unknown:
            raise KeyError(f"unknown realization IDs: {unknown}.")
        if self.lineage:
            parents = {self.lineage[item] for item in moving}
            required = {
                item for item in self.all_ids if self.lineage[item] in parents
            }
            if set(moving) != required:
                raise ValueError(
                    "all members of an affected lineage must be reassigned together."
                )
        moved = set(moving)
        values = {
            name: tuple(
                item for item in getattr(self, name) if item not in moved
            )
            for name in _PARTITIONS
        }
        values[partition] = (
            tuple(item for item in target if item not in moved) + moving
        )
        return replace(
            self,
            train=values["train"],
            validation=values["validation"],
            test=values["test"],
            strategy=f"{self.strategy}+manual",
        )

    def to_dict(self) -> dict[str, Any]:
        """Return the complete schema-2 JSON representation.

        Returns
        -------
        dict
            Mutable copy of partitions, seed, lineage, and strategy.

        Examples
        --------
        >>> RealizationSplit(("a",), (), ()).to_dict()["schema_version"]
        2
        """
        return {
            "schema_version": 2,
            "train": list(self.train),
            "validation": list(self.validation),
            "test": list(self.test),
            "seed": self.seed,
            "lineage": dict(self.lineage),
            "strategy": self.strategy,
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> RealizationSplit:
        """Restore schema-1 or schema-2 split state.

        Parameters
        ----------
        data : mapping
            Serialized split dictionary.

        Returns
        -------
        RealizationSplit
            Validated immutable runtime split.

        Raises
        ------
        ValueError
            If the schema is unsupported or assignments leak/overlap.

        Examples
        --------
        >>> split = RealizationSplit(("a",), (), ("b",), seed=1)
        >>> RealizationSplit.from_dict(split.to_dict()) == split
        True
        """
        version = data.get("schema_version", 1)
        if version not in {1, 2}:
            raise ValueError("unsupported RealizationSplit schema version.")
        return cls(
            tuple(data["train"]),
            tuple(data["validation"]),
            tuple(data["test"]),
            data.get("seed"),
            lineage={} if version == 1 else data.get("lineage", {}),
            strategy="random"
            if version == 1
            else data.get("strategy", "random"),
        )


def split_realizations(
    realization_ids: Sequence[str],
    *,
    validation_fraction: float = 0.1,
    test_fraction: float = 0.1,
    seed: int | None = 0,
    lineage: Mapping[str, str] | None = None,
) -> RealizationSplit:
    """Create a deterministic realization- or lineage-level random split.

    Parameters
    ----------
    realization_ids : sequence of str
        Unique geological realization IDs. Input order does not affect output.
    validation_fraction, test_fraction : float, default=0.1
        Target fractions in ``[0, 1)`` whose sum is less than one.
    seed : int or None, default=0
        NumPy random-generator seed.
    lineage : mapping, optional
        Complete ID-to-parent mapping. Whole lineages are assigned together,
        so realized fractions can differ from targets.

    Returns
    -------
    RealizationSplit
        Immutable leakage-checked assignment.

    Raises
    ------
    ValueError
        If IDs, fractions, or lineage coverage are invalid or no training group
        can remain.

    Examples
    --------
    >>> ids = [f"r{i}" for i in range(10)]
    >>> first = split_realizations(
    ...     ids, validation_fraction=0.2, test_fraction=0.2, seed=4
    ... )
    >>> second = split_realizations(
    ...     list(reversed(ids)),
    ...     validation_fraction=0.2,
    ...     test_fraction=0.2,
    ...     seed=4,
    ... )
    >>> first == second
    True

    Keep multiple noise variants of one parent in a single partition:

    >>> lineage = {"a-clean": "a", "a-noisy": "a", "b": "b", "c": "c"}
    >>> split = split_realizations(
    ...     list(lineage),
    ...     validation_fraction=0.25,
    ...     test_fraction=0.25,
    ...     lineage=lineage,
    ... )
    >>> split.partition("a-clean") == split.partition("a-noisy")
    True
    """
    ids = _ids(realization_ids, "realization_ids")
    if not ids:
        raise ValueError("realization_ids cannot be empty.")
    validation_fraction = _fraction(validation_fraction, "validation_fraction")
    test_fraction = _fraction(test_fraction, "test_fraction")
    if validation_fraction + test_fraction >= 1:
        raise ValueError(
            "validation_fraction + test_fraction must be less than 1."
        )
    seed = _seed(seed)
    lineage_state = (
        _lineage_map(ids, lineage)
        if lineage is not None
        else MappingProxyType({})
    )
    grouped = _groups(ids, lineage_state if lineage_state else None)
    required_partitions = (
        1 + int(validation_fraction > 0) + int(test_fraction > 0)
    )
    if len(grouped) < required_partitions:
        raise ValueError(
            "too few independent lineages for requested non-empty partitions."
        )
    group_names = np.asarray(sorted(grouped), dtype=object)
    order = np.random.default_rng(seed).permutation(len(group_names))
    available = [str(group_names[index]) for index in order]
    sizes = {name: len(grouped[name]) for name in grouped}
    n_validation = _target_count(len(ids), validation_fraction, len(ids))
    n_test = _target_count(len(ids), test_fraction, len(ids) - n_validation)
    test_groups, available = _take_groups(
        available,
        sizes,
        n_test,
        reserve=1 + int(validation_fraction > 0),
    )
    validation_groups, train_groups = _take_groups(
        available, sizes, n_validation, reserve=1
    )

    def _members(names: Sequence[str]) -> tuple[str, ...]:
        return tuple(item for name in names for item in grouped[name])

    strategy = "lineage_random" if lineage_state else "random"
    return RealizationSplit(
        train=_members(train_groups),
        validation=_members(validation_groups),
        test=_members(test_groups),
        seed=seed,
        lineage=lineage_state,
        strategy=strategy,
    )


def realization_folds(
    realization_ids: Sequence[str],
    *,
    n_splits: int = 5,
    seed: int | None = 0,
    lineage: Mapping[str, str] | None = None,
) -> tuple[RealizationSplit, ...]:
    """Build deterministic group-safe cross-validation test folds.

    Parameters
    ----------
    realization_ids : sequence of str
        Unique realization IDs.
    n_splits : int, default=5
        Number of folds. It cannot exceed independent lineage count.
    seed : int or None, default=0
        Seed controlling shuffled group order.
    lineage : mapping, optional
        Complete ID-to-parent mapping. Each lineage appears in exactly one test
        fold and never crosses train/test within a fold.

    Returns
    -------
    tuple of RealizationSplit
        Splits with empty validation partitions. Across the tuple, every input
        realization appears in test exactly once.

    Raises
    ------
    ValueError
        If fewer than two folds are requested or independent groups are
        insufficient.

    Examples
    --------
    >>> folds = realization_folds(["a", "b", "c", "d"], n_splits=2, seed=1)
    >>> len(folds)
    2
    >>> sorted(item for fold in folds for item in fold.test)
    ['a', 'b', 'c', 'd']
    """
    ids = _ids(realization_ids, "realization_ids")
    if not ids:
        raise ValueError("realization_ids cannot be empty.")
    if (
        not isinstance(n_splits, int)
        or isinstance(n_splits, bool)
        or n_splits < 2
    ):
        raise ValueError("n_splits must be an integer of at least two.")
    seed = _seed(seed)
    lineage_state = (
        _lineage_map(ids, lineage)
        if lineage is not None
        else MappingProxyType({})
    )
    grouped = _groups(ids, lineage_state if lineage_state else None)
    if n_splits > len(grouped):
        raise ValueError("n_splits cannot exceed independent lineage count.")
    names = np.asarray(sorted(grouped), dtype=object)
    order = np.random.default_rng(seed).permutation(len(names))
    shuffled = [str(names[index]) for index in order]
    buckets: list[list[str]] = [[] for _ in range(n_splits)]
    loads = [0] * n_splits
    for name in sorted(
        shuffled, key=lambda item: len(grouped[item]), reverse=True
    ):
        index = min(range(n_splits), key=lambda fold: (loads[fold], fold))
        buckets[index].append(name)
        loads[index] += len(grouped[name])
    folds = []
    for bucket in buckets:
        test_groups = set(bucket)
        test = tuple(item for name in bucket for item in grouped[name])
        train = tuple(
            item
            for name in shuffled
            if name not in test_groups
            for item in grouped[name]
        )
        folds.append(
            RealizationSplit(
                train,
                (),
                test,
                seed=seed,
                lineage=lineage_state,
                strategy=f"lineage_kfold_{n_splits}"
                if lineage_state
                else f"kfold_{n_splits}",
            )
        )
    return tuple(folds)
