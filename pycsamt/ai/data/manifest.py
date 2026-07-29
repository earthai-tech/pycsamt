# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Immutable, versioned provenance manifests for generated EM datasets.

A manifest identifies how a dataset was generated, which geological
realizations belong to each split, and which external artifacts belong to the
dataset.  It contains no training arrays itself.  All hashes use SHA-256 and
all serialized content uses deterministic canonical JSON.
"""

from __future__ import annotations

import hashlib
import json
import re
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field, replace
from datetime import datetime, timezone
from pathlib import Path, PurePosixPath
from types import MappingProxyType
from typing import Any

from .splits import RealizationSplit

__all__ = [
    "ArtifactRecord",
    "DatasetManifest",
    "canonical_hash",
    "sha256_file",
]

_HEX_64 = re.compile(r"^[0-9a-f]{64}$")
_DATASET_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")


def _json_value(value: Any) -> Any:
    if isinstance(value, Mapping):
        result = {}
        for key, item in value.items():
            normalized = str(key)
            if normalized in result:
                raise ValueError(f"mapping contains colliding JSON key {normalized!r}.")
            result[normalized] = _json_value(item)
        return result
    if isinstance(value, (list, tuple)):
        return [_json_value(item) for item in value]
    return value


def _freeze(value: Any) -> Any:
    if isinstance(value, Mapping):
        return MappingProxyType(
            {str(key): _freeze(item) for key, item in value.items()}
        )
    if isinstance(value, (list, tuple)):
        return tuple(_freeze(item) for item in value)
    return value


def _canonical_json(value: Any) -> str:
    normalized = _json_value(value)
    try:
        return json.dumps(
            normalized,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=False,
            allow_nan=False,
        )
    except (TypeError, ValueError) as exc:
        raise ValueError(
            "value must contain finite JSON-serializable data with string-compatible keys."
        ) from exc


def _timestamp(value: str | None) -> str | None:
    if value is None:
        return None
    text = str(value).strip()
    if not text:
        raise ValueError("created_utc cannot be empty.")
    try:
        parsed = datetime.fromisoformat(
            text[:-1] + "+00:00" if text.endswith("Z") else text
        )
    except ValueError as exc:
        raise ValueError("created_utc must be an ISO-8601 timestamp.") from exc
    if parsed.tzinfo is None or parsed.utcoffset() is None:
        raise ValueError("created_utc must include a timezone.")
    utc = parsed.astimezone(timezone.utc)
    return utc.isoformat(timespec="seconds").replace("+00:00", "Z")


def _artifact_path(value: str) -> str:
    text = str(value).strip().replace("\\", "/")
    path = PurePosixPath(text)
    if (
        not text
        or ":" in text
        or any(ord(character) < 32 for character in text)
        or path.is_absolute()
        or ".." in path.parts
        or "." in path.parts
    ):
        raise ValueError(
            "artifact paths must be normalized relative paths without '..'."
        )
    if text.endswith("/"):
        raise ValueError("artifact paths must identify files, not directories.")
    return path.as_posix()


def _resolved_artifact(root: Path, relative_path: str) -> Path:
    base = root.resolve()
    candidate = base.joinpath(*PurePosixPath(relative_path).parts).resolve()
    try:
        candidate.relative_to(base)
    except ValueError as exc:
        raise ValueError(
            f"artifact path resolves outside verification root: {relative_path}"
        ) from exc
    return candidate


def canonical_hash(value: Any) -> str:
    """Return the SHA-256 digest of deterministic canonical JSON.

    Parameters
    ----------
    value : Any
        Finite JSON-serializable value. Mapping keys are converted to strings,
        mappings are sorted recursively, and insignificant whitespace is
        removed before hashing.

    Returns
    -------
    str
        Lowercase 64-character hexadecimal SHA-256 digest.

    Raises
    ------
    ValueError
        If ``value`` contains NaN, infinity, bytes, arrays, or another object
        that cannot be represented safely as JSON.

    Examples
    --------
    Mapping insertion order does not affect the digest:

    >>> canonical_hash({"a": 1, "b": 2}) == canonical_hash({"b": 2, "a": 1})
    True
    >>> len(canonical_hash({"frequencies_hz": [100.0, 10.0]}))
    64
    """
    return hashlib.sha256(_canonical_json(value).encode("utf-8")).hexdigest()


def sha256_file(path: str | Path, *, chunk_size: int = 1024 * 1024) -> str:
    """Hash a file without loading the complete artifact into memory.

    Parameters
    ----------
    path : str or pathlib.Path
        Existing regular file to hash.
    chunk_size : int, default=1048576
        Positive number of bytes read per iteration.

    Returns
    -------
    str
        Lowercase hexadecimal SHA-256 digest.

    Raises
    ------
    ValueError
        If ``chunk_size`` is not a positive integer.
    OSError
        If the file cannot be opened or read.

    Examples
    --------
    >>> from tempfile import TemporaryDirectory
    >>> with TemporaryDirectory() as directory:
    ...     path = Path(directory) / "artifact.bin"
    ...     _ = path.write_bytes(b"pycsamt")
    ...     digest = sha256_file(path)
    >>> len(digest)
    64
    """
    if (
        not isinstance(chunk_size, int)
        or isinstance(chunk_size, bool)
        or chunk_size <= 0
    ):
        raise ValueError("chunk_size must be a positive integer.")
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        while True:
            block = stream.read(chunk_size)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


@dataclass(frozen=True)
class ArtifactRecord:
    """Integrity metadata for one external dataset artifact.

    Parameters
    ----------
    sha256 : str
        Lowercase 64-character SHA-256 digest of the complete file.
    size_bytes : int or None, optional
        Exact file size. When present it is checked before hashing.
    media_type : str or None, optional
        MIME type such as ``"application/x-npz"``.
    role : str or None, optional
        Human-readable role such as ``"responses"`` or ``"models"``.

    Examples
    --------
    >>> record = ArtifactRecord("0" * 64, size_bytes=1024, role="models")
    >>> record.algorithm
    'sha256'
    """

    sha256: str
    size_bytes: int | None = None
    media_type: str | None = None
    role: str | None = None

    def __post_init__(self) -> None:
        digest = str(self.sha256).lower()
        if digest.startswith("sha256:"):
            digest = digest[7:]
        if not _HEX_64.fullmatch(digest):
            raise ValueError("sha256 must contain exactly 64 hexadecimal characters.")
        if self.size_bytes is not None and (
            not isinstance(self.size_bytes, int)
            or isinstance(self.size_bytes, bool)
            or self.size_bytes < 0
        ):
            raise ValueError("size_bytes must be a non-negative integer or None.")
        media_type = None if self.media_type is None else str(self.media_type).strip()
        role = None if self.role is None else str(self.role).strip()
        if self.media_type is not None and not media_type:
            raise ValueError("media_type cannot be empty.")
        if self.role is not None and not role:
            raise ValueError("role cannot be empty.")
        object.__setattr__(self, "sha256", digest)
        object.__setattr__(self, "media_type", media_type)
        object.__setattr__(self, "role", role)

    @property
    def algorithm(self) -> str:
        """Return the checksum algorithm identifier.

        Returns
        -------
        str
            Always ``"sha256"`` for the current artifact schema.

        Examples
        --------
        >>> ArtifactRecord("a" * 64).algorithm
        'sha256'
        """
        return "sha256"

    @classmethod
    def from_file(
        cls,
        path: str | Path,
        *,
        media_type: str | None = None,
        role: str | None = None,
        chunk_size: int = 1024 * 1024,
    ) -> ArtifactRecord:
        """Create an integrity record from an existing file.

        Parameters
        ----------
        path : str or pathlib.Path
            Existing regular file.
        media_type, role : str, optional
            Optional descriptive metadata.
        chunk_size : int, default=1048576
            Bytes read per hashing iteration.

        Returns
        -------
        ArtifactRecord
            Digest and exact file size captured from disk.

        Examples
        --------
        >>> from tempfile import TemporaryDirectory
        >>> with TemporaryDirectory() as directory:
        ...     path = Path(directory) / "models.npz"
        ...     _ = path.write_bytes(b"model-data")
        ...     record = ArtifactRecord.from_file(path, role="models")
        >>> record.size_bytes
        10
        """
        source = Path(path)
        if not source.is_file():
            raise FileNotFoundError(f"artifact is not a regular file: {source}")
        return cls(
            sha256_file(source, chunk_size=chunk_size),
            size_bytes=source.stat().st_size,
            media_type=media_type,
            role=role,
        )

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-serializable artifact record.

        Returns
        -------
        dict
            Versioned checksum, size, media type, and role fields.

        Examples
        --------
        >>> ArtifactRecord("f" * 64, size_bytes=2).to_dict()["size_bytes"]
        2
        """
        return {
            "schema_version": 1,
            "algorithm": self.algorithm,
            "sha256": self.sha256,
            "size_bytes": self.size_bytes,
            "media_type": self.media_type,
            "role": self.role,
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> ArtifactRecord:
        """Restore a validated artifact record.

        Parameters
        ----------
        data : mapping
            Versioned state returned by :meth:`to_dict`.

        Returns
        -------
        ArtifactRecord
            Immutable integrity record.

        Raises
        ------
        ValueError
            If the schema or checksum algorithm is unsupported.

        Examples
        --------
        >>> state = ArtifactRecord("1" * 64).to_dict()
        >>> ArtifactRecord.from_dict(state).sha256 == "1" * 64
        True
        """
        if data.get("schema_version", 1) != 1:
            raise ValueError("unsupported ArtifactRecord schema version.")
        if data.get("algorithm", "sha256") != "sha256":
            raise ValueError("unsupported artifact checksum algorithm.")
        return cls(
            sha256=data["sha256"],
            size_bytes=data.get("size_bytes"),
            media_type=data.get("media_type"),
            role=data.get("role"),
        )


@dataclass(frozen=True)
class DatasetManifest:
    """Identify a generated dataset and its complete reproducibility state.

    Parameters
    ----------
    dataset_id : str
        Portable identifier containing letters, digits, dots, underscores, or
        hyphens. It must start with a letter or digit.
    generator, generator_version : str
        Fully qualified generator name and its version or source revision.
    configuration : mapping
        Finite JSON-compatible generator configuration. It is recursively
        copied and frozen.
    split : RealizationSplit
        Disjoint realization-level train/validation/test assignment.
    sample_count : int
        Number of samples represented by the dataset.
    created_utc : str or None, optional
        Timezone-aware ISO-8601 creation time. It is normalized to UTC.
    artifacts : mapping, optional
        Normalized relative paths mapped to :class:`ArtifactRecord` objects or
        their serialized dictionaries.
    schema_version : int, default=2
        Manifest format version. New manifests use version 2.

    Examples
    --------
    >>> split = RealizationSplit(("r1", "r2"), ("r3",), ("r4",), seed=7)
    >>> manifest = DatasetManifest(
    ...     dataset_id="willy-2d-v1",
    ...     generator="pycsamt.ai.geology.correlated2d",
    ...     generator_version="0.1.0",
    ...     configuration={"seed": 7, "correlation_m": [1000, 100]},
    ...     split=split,
    ...     sample_count=4,
    ... )
    >>> len(manifest.configuration_hash)
    64
    """

    dataset_id: str
    generator: str
    generator_version: str
    configuration: Mapping[str, Any]
    split: RealizationSplit
    sample_count: int
    created_utc: str | None = None
    artifacts: Mapping[str, ArtifactRecord | Mapping[str, Any] | str] = field(
        default_factory=dict
    )
    schema_version: int = 2

    def __post_init__(self) -> None:
        if self.schema_version != 2:
            raise ValueError("new DatasetManifest objects require schema version 2.")
        dataset_id = str(self.dataset_id).strip()
        if not _DATASET_ID.fullmatch(dataset_id):
            raise ValueError(
                "dataset_id must be portable and contain only letters, digits, '.', '_', or '-'."
            )
        generator = str(self.generator).strip()
        generator_version = str(self.generator_version).strip()
        if not generator or not generator_version:
            raise ValueError("generator and generator_version cannot be empty.")
        if not isinstance(self.split, RealizationSplit):
            raise TypeError("split must be a RealizationSplit.")
        if (
            not isinstance(self.sample_count, int)
            or isinstance(self.sample_count, bool)
            or self.sample_count < 0
        ):
            raise ValueError("sample_count must be a non-negative integer.")
        configuration = json.loads(_canonical_json(self.configuration))
        artifacts: dict[str, ArtifactRecord] = {}
        for path, value in self.artifacts.items():
            key = _artifact_path(path)
            if key in artifacts:
                raise ValueError(f"duplicate normalized artifact path {key!r}.")
            if isinstance(value, ArtifactRecord):
                record = value
            elif isinstance(value, Mapping):
                record = ArtifactRecord.from_dict(value)
            elif isinstance(value, str):
                record = ArtifactRecord(value)
            else:
                raise TypeError(
                    "artifact values must be records, mappings, or SHA-256 strings."
                )
            artifacts[key] = record
        object.__setattr__(self, "dataset_id", dataset_id)
        object.__setattr__(self, "generator", generator)
        object.__setattr__(self, "generator_version", generator_version)
        object.__setattr__(self, "configuration", _freeze(configuration))
        object.__setattr__(self, "created_utc", _timestamp(self.created_utc))
        object.__setattr__(self, "artifacts", MappingProxyType(artifacts))

    @property
    def configuration_hash(self) -> str:
        """Return the canonical configuration digest.

        Returns
        -------
        str
            SHA-256 digest of generator configuration only.

        Examples
        --------
        >>> split = RealizationSplit(("r1",), (), ())
        >>> m = DatasetManifest("d", "g", "1", {"seed": 0}, split, 1)
        >>> m.configuration_hash == canonical_hash({"seed": 0})
        True
        """
        return canonical_hash(self.configuration)

    @property
    def manifest_hash(self) -> str:
        """Return a digest of the complete serialized manifest.

        Returns
        -------
        str
            SHA-256 digest covering configuration, split, timestamps, and all
            artifact records.

        Examples
        --------
        >>> split = RealizationSplit(("r1",), (), ())
        >>> m = DatasetManifest("d", "g", "1", {}, split, 1)
        >>> len(m.manifest_hash)
        64
        """
        return canonical_hash(self.to_dict())

    @property
    def realization_count(self) -> int:
        """Return the total number of split realizations.

        Returns
        -------
        int
            Length of the combined train, validation, and test ID sets.

        Examples
        --------
        >>> split = RealizationSplit(("a", "b"), ("c",), ())
        >>> DatasetManifest("d", "g", "1", {}, split, 3).realization_count
        3
        """
        return len(self.split.all_ids)

    def with_artifact(
        self, path: str, record: ArtifactRecord | Mapping[str, Any] | str
    ) -> DatasetManifest:
        """Return a copy containing or replacing one artifact record.

        Parameters
        ----------
        path : str
            Portable relative artifact path.
        record : ArtifactRecord, mapping, or str
            Integrity record, serialized record, or SHA-256 digest.

        Returns
        -------
        DatasetManifest
            New immutable manifest; the original is unchanged.

        Examples
        --------
        >>> split = RealizationSplit(("r1",), (), ())
        >>> m = DatasetManifest("d", "g", "1", {}, split, 1)
        >>> updated = m.with_artifact("data/models.npz", "a" * 64)
        >>> list(updated.artifacts)
        ['data/models.npz']
        """
        artifacts = dict(self.artifacts)
        artifacts[path] = record
        return replace(self, artifacts=artifacts)

    def verify_artifacts(
        self,
        root: str | Path,
        *,
        paths: Sequence[str] | None = None,
        raise_on_error: bool = False,
    ) -> dict[str, bool]:
        """Verify recorded artifact sizes and SHA-256 digests on disk.

        Parameters
        ----------
        root : str or pathlib.Path
            Directory against which relative artifact paths are resolved.
        paths : sequence of str, optional
            Subset of recorded paths. By default all artifacts are checked.
        raise_on_error : bool, default=False
            Raise on the first missing, size-mismatched, or hash-mismatched
            artifact instead of returning ``False`` for it.

        Returns
        -------
        dict
            Normalized artifact paths mapped to verification results.

        Raises
        ------
        KeyError
            If a requested path is not recorded.
        ValueError
            If ``root`` is not a directory or verification fails while
            ``raise_on_error`` is true.

        Examples
        --------
        >>> from tempfile import TemporaryDirectory
        >>> split = RealizationSplit(("r1",), (), ())
        >>> with TemporaryDirectory() as directory:
        ...     root = Path(directory)
        ...     file = root / "data.bin"
        ...     _ = file.write_bytes(b"data")
        ...     record = ArtifactRecord.from_file(file)
        ...     manifest = DatasetManifest(
        ...         "d", "g", "1", {}, split, 1, artifacts={"data.bin": record}
        ...     )
        ...     result = manifest.verify_artifacts(root)
        >>> result
        {'data.bin': True}
        """
        base = Path(root)
        if not base.is_dir():
            raise ValueError("root must be an existing directory.")
        selected = (
            list(self.artifacts)
            if paths is None
            else [_artifact_path(path) for path in paths]
        )
        results: dict[str, bool] = {}
        for path in selected:
            if path not in self.artifacts:
                raise KeyError(f"unrecorded artifact path {path!r}.")
            record = self.artifacts[path]
            candidate = _resolved_artifact(base, path)
            ok = candidate.is_file()
            if ok and record.size_bytes is not None:
                ok = candidate.stat().st_size == record.size_bytes
            if ok:
                ok = sha256_file(candidate) == record.sha256
            results[path] = ok
            if not ok and raise_on_error:
                raise ValueError(f"artifact integrity verification failed: {path}")
        return results

    def to_dict(self) -> dict[str, Any]:
        """Return the complete schema-2 JSON representation.

        Returns
        -------
        dict
            Mutable JSON-compatible copy including the configuration digest.

        Examples
        --------
        >>> split = RealizationSplit(("r1",), (), ())
        >>> m = DatasetManifest("d", "g", "1", {}, split, 1)
        >>> m.to_dict()["schema_version"]
        2
        """
        return {
            "schema_version": self.schema_version,
            "dataset_id": self.dataset_id,
            "generator": self.generator,
            "generator_version": self.generator_version,
            "configuration": _json_value(self.configuration),
            "configuration_hash": self.configuration_hash,
            "split": self.split.to_dict(),
            "sample_count": self.sample_count,
            "created_utc": self.created_utc,
            "artifacts": {
                path: record.to_dict() for path, record in self.artifacts.items()
            },
        }

    def write_json(self, path: str | Path, *, overwrite: bool = True) -> Path:
        """Write a deterministic, human-readable manifest file.

        Parameters
        ----------
        path : str or pathlib.Path
            Destination JSON file.
        overwrite : bool, default=True
            Permit replacement of an existing file.

        Returns
        -------
        pathlib.Path
            Destination path.

        Raises
        ------
        FileExistsError
            If the destination exists and ``overwrite`` is false.

        Examples
        --------
        >>> from tempfile import TemporaryDirectory
        >>> split = RealizationSplit(("r1",), (), ())
        >>> m = DatasetManifest("d", "g", "1", {}, split, 1)
        >>> with TemporaryDirectory() as directory:
        ...     path = m.write_json(Path(directory) / "manifest.json")
        ...     loaded = DatasetManifest.read_json(path)
        >>> loaded.manifest_hash == m.manifest_hash
        True
        """
        target = Path(path)
        if target.exists() and not overwrite:
            raise FileExistsError(f"manifest already exists: {target}")
        target.write_text(
            json.dumps(self.to_dict(), indent=2, sort_keys=True, ensure_ascii=False)
            + "\n",
            encoding="utf-8",
        )
        return target

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> DatasetManifest:
        """Restore a manifest and verify its recorded configuration digest.

        Parameters
        ----------
        data : mapping
            Schema-1 or schema-2 serialized manifest.

        Returns
        -------
        DatasetManifest
            Validated schema-2 runtime object.

        Raises
        ------
        ValueError
            If the schema is unsupported, required content is invalid, or the
            recorded configuration hash does not match its configuration.

        Examples
        --------
        >>> split = RealizationSplit(("r1",), (), ())
        >>> original = DatasetManifest("d", "g", "1", {"seed": 2}, split, 1)
        >>> restored = DatasetManifest.from_dict(original.to_dict())
        >>> restored.configuration_hash == original.configuration_hash
        True
        """
        version = data.get("schema_version", 1)
        if version not in {1, 2}:
            raise ValueError(f"unsupported DatasetManifest schema version {version!r}.")
        configuration = dict(data["configuration"])
        recorded_hash = data.get("configuration_hash")
        if recorded_hash is not None and recorded_hash != canonical_hash(configuration):
            raise ValueError("configuration_hash does not match configuration.")
        raw_artifacts = data.get("artifacts", {})
        artifacts = {}
        for path, record in raw_artifacts.items():
            if version == 1 and isinstance(record, str):
                artifacts[path] = record
            else:
                artifacts[path] = ArtifactRecord.from_dict(record)
        return cls(
            dataset_id=data["dataset_id"],
            generator=data["generator"],
            generator_version=data["generator_version"],
            configuration=configuration,
            split=RealizationSplit.from_dict(data["split"]),
            sample_count=data["sample_count"],
            created_utc=data.get("created_utc"),
            artifacts=artifacts,
            schema_version=2,
        )

    @classmethod
    def read_json(cls, path: str | Path) -> DatasetManifest:
        """Read and validate a UTF-8 JSON manifest.

        Parameters
        ----------
        path : str or pathlib.Path
            Existing manifest file.

        Returns
        -------
        DatasetManifest
            Validated immutable manifest.

        Raises
        ------
        OSError
            If the file cannot be read.
        json.JSONDecodeError
            If its content is not valid JSON.
        ValueError
            If decoded content violates the manifest contract.

        Examples
        --------
        >>> from tempfile import TemporaryDirectory
        >>> split = RealizationSplit(("r1",), (), ())
        >>> source = DatasetManifest("d", "g", "1", {}, split, 1)
        >>> with TemporaryDirectory() as directory:
        ...     path = source.write_json(Path(directory) / "manifest.json")
        ...     loaded = DatasetManifest.read_json(path)
        >>> loaded.dataset_id
        'd'
        """
        return cls.from_dict(json.loads(Path(path).read_text(encoding="utf-8")))
