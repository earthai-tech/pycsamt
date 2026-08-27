# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Machine-learning-model provenance for a PCSF/PCSM AI-inversion result.

This is deliberately distinct from :class:`pycsamt.metadata.ProvenanceMeta`,
which describes *survey/data* provenance (who collected/submitted an EM
dataset). :class:`ModelProvenance` describes the *model* that produced a
resistivity result -- its architecture, training framework, checkpoint, and
hyperparameters -- so a shared ``.pcsf``/``.pcsm`` file is independently
checkable by anyone re-running the model, not just re-readable.

Every field is optional free text: this is meant to be usable by literally
any AI/DL inversion tool, including third-party ones with no dependency on
pycsamt at all, so nothing here is enforced beyond basic typing.
"""

from __future__ import annotations

import hashlib
from dataclasses import dataclass, field
from os import PathLike
from pathlib import Path
from typing import Any

from ..api.property import PyCSAMTObject

__all__ = ["ModelProvenance", "compute_checkpoint_hash"]


@dataclass(repr=False)
class ModelProvenance(PyCSAMTObject):
    """Describe the AI/DL model that produced a resistivity result.

    Parameters
    ----------
    architecture : str, optional
        Free-text model family/architecture, e.g. ``"UNet"``, ``"GCN"``,
        ``"ResNet18"``, or any third-party name.
    framework, framework_version : str, optional
        e.g. ``"pytorch"``/``"2.3.0"``, ``"tensorflow"``/``"2.16.1"``.
    checkpoint : str, optional
        Path or identifier of the trained weights used to produce the
        result (e.g. a filename, a model-hub id, a DOI).
    checkpoint_sha256 : str, optional
        SHA-256 hex digest of the checkpoint file, so a reader can
        verify they are re-running the exact weights this result
        claims -- see :func:`compute_checkpoint_hash`.
    training_data : str, optional
        Free-text reference to the training dataset (name, DOI, path).
    hyperparameters : dict, default {}
        Free-form training/model hyperparameters.
    random_seed : int, optional
        Seed used for training and/or inference, when reproducibility
        depends on it.
    git_commit : str, optional
        Commit hash of the code that produced this result.
    authors : list of str, default []
        Model authors/maintainers.
    contact : str, optional
        Contact e-mail or URL for questions about this result.
    notes : str, optional
        Free-text notes not covered by the fields above.
    extra : dict, default {}
        Unmodelled provenance fields, retained losslessly.

    Examples
    --------
    >>> from pycsamt.format.provenance import ModelProvenance
    >>> prov = ModelProvenance(
    ...     architecture="UNet",
    ...     framework="pytorch",
    ...     framework_version="2.3.0",
    ...     checkpoint="unet_v3.pt",
    ...     random_seed=42,
    ... )
    >>> prov.to_dict()["architecture"]
    'UNet'
    """

    architecture: str = ""
    framework: str = ""
    framework_version: str = ""
    checkpoint: str = ""
    checkpoint_sha256: str = ""
    training_data: str = ""
    hyperparameters: dict[str, Any] = field(default_factory=dict)
    random_seed: int | None = None
    git_commit: str = ""
    authors: list[str] = field(default_factory=list)
    contact: str = ""
    notes: str = ""
    extra: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.hyperparameters = dict(self.hyperparameters or {})
        self.authors = list(self.authors or [])
        self.extra = dict(self.extra or {})

    def validate(self) -> None:
        if self.random_seed is not None and not isinstance(self.random_seed, int):
            raise TypeError("ModelProvenance.random_seed must be an int or None")

    def to_dict(self) -> dict[str, Any]:
        """Plain-dict form, suitable for ``PCSFModel.metadata['model_provenance']``."""
        return {
            "architecture": self.architecture,
            "framework": self.framework,
            "framework_version": self.framework_version,
            "checkpoint": self.checkpoint,
            "checkpoint_sha256": self.checkpoint_sha256,
            "training_data": self.training_data,
            "hyperparameters": dict(self.hyperparameters),
            "random_seed": self.random_seed,
            "git_commit": self.git_commit,
            "authors": list(self.authors),
            "contact": self.contact,
            "notes": self.notes,
            "extra": dict(self.extra),
        }


def compute_checkpoint_hash(path: str | PathLike, *, chunk_size: int = 1 << 20) -> str:
    """SHA-256 hex digest of a model-checkpoint file, streamed in chunks.

    Parameters
    ----------
    path : path-like
        The checkpoint file to hash (e.g. a ``.pt``/``.h5``/``.onnx`` file).
    chunk_size : int, default 1 MiB
        Read block size; large checkpoints are hashed without loading
        the whole file into memory.

    Returns
    -------
    str
        Lowercase hex digest, directly comparable to
        :attr:`ModelProvenance.checkpoint_sha256`.

    Examples
    --------
    >>> from pycsamt.format.provenance import compute_checkpoint_hash
    >>> compute_checkpoint_hash("unet_v3.pt")  # doctest: +SKIP
    '3b1c...'
    """
    digest = hashlib.sha256()
    with open(Path(path), "rb") as fh:
        for chunk in iter(lambda: fh.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()
