# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

"""Inversion-model format detection and reader registration.

Sibling to :mod:`pycsamt.io.formats` (transfer-function formats): this
module knows how to identify registered inversion-*model* formats and
obtain their reader/writer callables, but carries no geophysical
science logic of its own. PCSF (:mod:`pycsamt.format`) is registered
here rather than in :mod:`pycsamt.io.formats` because it is not a
transfer function — conflating the two registries would make
``list_tf_formats()`` misreport what it lists.
"""

from __future__ import annotations

import gzip
from dataclasses import dataclass
from os import PathLike
from pathlib import Path
from typing import Any, Callable

from ..core._named_registry import (
    NamedRegistry,
    normalize_extension as _normalize_extension,
    normalize_key as _normalize_name,
)

__all__ = [
    "ModelFormat",
    "ModelFormatError",
    "register_model_format",
    "get_model_format",
    "list_model_formats",
    "detect_model_format",
    "get_model_format_for_target",
]


ModelReader = Callable[..., Any]
ModelWriter = Callable[..., Any]
ModelDetector = Callable[[Any], bool]


class ModelFormatError(ValueError):
    """Raised when an inversion-model format cannot be resolved."""


@dataclass(frozen=True)
class ModelFormat:
    """Definition of one registered inversion-model serialization.

    Parameters
    ----------
    name : str
        Canonical, case-insensitive format name.
    reader : callable
        Callable accepting the source as its first argument.
    writer : callable, optional
        Callable accepting ``(object, target)`` for serialization.
    extensions : tuple of str, optional
        Filename extensions used only as hints. Detection remains
        content based whenever a detector is available.
    detector : callable, optional
        Callable returning ``True`` when the source belongs to this
        format.
    aliases : tuple of str, optional
        Alternative names accepted by :func:`get_model_format`.
    description : str, optional
        Human-readable summary.
    """

    name: str
    reader: ModelReader
    writer: ModelWriter | None = None
    extensions: tuple[str, ...] = ()
    detector: ModelDetector | None = None
    aliases: tuple[str, ...] = ()
    description: str = ""


_registry: NamedRegistry[ModelFormat] = NamedRegistry(
    kind="inversion-model format",
    error_cls=ModelFormatError,
)


def register_model_format(
    name: str,
    *,
    reader: ModelReader,
    writer: ModelWriter | None = None,
    extensions: tuple[str, ...] | list[str] = (),
    detector: ModelDetector | None = None,
    aliases: tuple[str, ...] | list[str] = (),
    description: str = "",
    replace: bool = False,
) -> ModelFormat:
    """Register one inversion-model reader definition."""
    canonical = _normalize_name(name)
    if not callable(reader):
        raise TypeError("reader must be callable")
    if writer is not None and not callable(writer):
        raise TypeError("writer must be callable or None")
    if detector is not None and not callable(detector):
        raise TypeError("detector must be callable or None")

    ext_values = tuple(
        ext
        for ext in (_normalize_extension(item) for item in extensions)
        if ext
    )
    alias_values = tuple(_normalize_name(item) for item in aliases)

    definition = ModelFormat(
        name=canonical,
        reader=reader,
        writer=writer,
        extensions=ext_values,
        detector=detector,
        aliases=alias_values,
        description=str(description or ""),
    )
    _registry.register(
        canonical, definition, aliases=alias_values, replace=replace
    )
    return definition


def get_model_format(name: str) -> ModelFormat:
    """Return a registered format by canonical name or alias."""
    result = _registry.get(name)
    if result is None:
        available = ", ".join(_registry.names()) or "<none>"
        raise ModelFormatError(
            f"unknown inversion-model format {name!r}; "
            f"available: {available}"
        )
    return result


def list_model_formats() -> dict[str, dict[str, Any]]:
    """Return stable public information about registered formats."""
    return {
        spec.name: {
            "extensions": spec.extensions,
            "aliases": spec.aliases,
            "description": spec.description,
            "readable": callable(spec.reader),
            "writable": callable(spec.writer),
        }
        for spec in sorted(_registry.all(), key=lambda item: item.name)
    }


def get_model_format_for_target(target: Any) -> ModelFormat:
    """Resolve a writable format from a target filename extension.

    Content detection cannot be used for a file that does not yet
    exist. Therefore writing uses explicit format selection or a
    unique registered extension.
    """
    if not isinstance(target, (str, PathLike)):
        raise ModelFormatError(
            "a writable stream requires an explicit inversion-model format"
        )
    try:
        suffix = Path(target).suffix.lower()
    except (OSError, ValueError) as exc:
        raise ModelFormatError(
            "cannot infer output format from target"
        ) from exc
    if not suffix:
        raise ModelFormatError(
            "output target has no extension; specify format explicitly"
        )
    candidates = [
        spec
        for spec in _registry.all()
        if spec.writer is not None and suffix in spec.extensions
    ]
    if len(candidates) == 1:
        return candidates[0]
    if not candidates:
        raise ModelFormatError(
            f"no writable inversion-model format registered for {suffix!r}"
        )
    raise ModelFormatError(
        f"ambiguous output extension {suffix!r}; specify format explicitly"
    )


def detect_model_format(source: Any) -> str:
    """Detect the canonical inversion-model format from source content."""
    matches = _registry.match_by_detector(source)

    if len(matches) == 1:
        return matches[0]
    if len(matches) > 1:
        raise ModelFormatError(
            "ambiguous inversion-model format: " + ", ".join(matches)
        )

    hint = ""
    if isinstance(source, (str, PathLike)):
        try:
            suffix = Path(source).suffix.lower()
        except (OSError, ValueError):
            suffix = ""
        if suffix:
            candidates = _registry.match_by_extension(suffix)
            if candidates:
                hint = (
                    f" Extension {suffix!r} suggests "
                    f"{', '.join(candidates)}, but the content did not "
                    "validate as that format."
                )
    raise ModelFormatError(
        "unable to detect a supported inversion-model format from "
        f"content.{hint}"
    )


def _is_pcsf(source: Any) -> bool:
    """Detect a PCSF file by its HDF5 signature + root attribute."""
    try:
        if isinstance(source, (str, PathLike)):
            path = Path(source)
            if not path.is_file():
                return False
            with path.open("rb") as stream:
                signature = stream.read(8)
            if signature != b"\x89HDF\r\n\x1a\n":
                return False
            import h5py

            with h5py.File(path, "r") as fh:
                return "pcsf_version" in fh.attrs
    except (OSError, ValueError, ImportError):
        return False
    return False


def _read_pcsf(source: Any, **kwargs: Any) -> Any:
    # Lazy import keeps ``import pycsamt.io`` independent of h5py/PCSF
    # internals, matching pycsamt.io.formats's lazy-reader convention.
    from ..format import read_pcsf

    return read_pcsf(source, **kwargs)


def _write_pcsf(obj: Any, target: Any, **kwargs: Any) -> Any:
    from ..format import write_pcsf

    return write_pcsf(obj, target, **kwargs)


def _pcsm_head(path: Path, limit: int = 4096) -> str:
    """First *limit* decoded chars, transparently ungzipping ``.gz``
    content (by trying, not by trusting the extension — this is a
    content detector).
    """
    try:
        with gzip.open(path, "rt", encoding="utf-8", errors="ignore") as fh:
            return fh.read(limit)
    except OSError:
        with path.open("rt", encoding="utf-8", errors="ignore") as fh:
            return fh.read(limit)


def _is_pcsm(source: Any) -> bool:
    """Detect a PCSM file by a ``PCSM_VERSION`` header line among its
    first few lines (plain or gzip-compressed).
    """
    try:
        if not isinstance(source, (str, PathLike)):
            return False
        path = Path(source)
        if not path.is_file():
            return False
        for line in _pcsm_head(path).splitlines()[:10]:
            if line.split(None, 1)[:1] == ["PCSM_VERSION"]:
                return True
        return False
    except (OSError, ValueError, UnicodeDecodeError):
        return False


def _read_pcsm(source: Any, **kwargs: Any) -> Any:
    from ..format import read_pcsm

    return read_pcsm(source, **kwargs)


def _write_pcsm(obj: Any, target: Any, **kwargs: Any) -> Any:
    from ..format import write_pcsm

    return write_pcsm(obj, target, **kwargs)


# Built-in format definitions. Reader/writer imports remain lazy by design.
register_model_format(
    "pcsf",
    reader=_read_pcsf,
    writer=_write_pcsf,
    extensions=(".pcsf",),
    detector=_is_pcsf,
    description=(
        "pyCSAMT Common Subsurface Format — backend-neutral inversion "
        "result container (see PYCSAMT-PCSF-INVERSION-FORMAT-PLAN.md)"
    ),
)
register_model_format(
    "pcsm",
    reader=_read_pcsm,
    writer=_write_pcsm,
    # ".gz" is a hint too: Path("model.pcsm.gz").suffix is ".gz" alone
    # (only the last suffix), and get_model_format_for_target resolves
    # a writer by suffix when the target doesn't exist yet to sniff.
    extensions=(".pcsm", ".gz"),
    detector=_is_pcsm,
    description=(
        "PCSM — pyCSAMT Common Subsurface Markup, the hand-editable "
        "ASCII sibling of PCSF (see pycsamt/format/SPEC.md section 9)"
    ),
)
