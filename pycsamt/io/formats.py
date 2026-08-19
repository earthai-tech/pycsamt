# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

"""Transfer-function format detection and reader registration.

This module is intentionally serialization-focused. It knows how to identify
registered transfer-function formats and how to obtain their reader callables,
but it does not contain electromagnetic science logic.
"""

from __future__ import annotations

from dataclasses import dataclass
from os import PathLike
from pathlib import Path
import re
from typing import Any, Callable

from ..core._named_registry import (
    NamedRegistry,
    normalize_extension as _normalize_extension,
    normalize_key as _normalize_name,
)

__all__ = [
    "TransferFunctionFormat",
    "TransferFunctionFormatError",
    "register_tf_format",
    "get_tf_format",
    "list_tf_formats",
    "detect_tf_format",
    "get_tf_format_for_target",
]


TFReader = Callable[..., Any]
TFWriter = Callable[..., Any]
TFDetector = Callable[[Any], bool]


class TransferFunctionFormatError(ValueError):
    """Raised when a transfer-function format cannot be resolved."""


@dataclass(frozen=True)
class TransferFunctionFormat:
    """Definition of one registered transfer-function serialization.

    Parameters
    ----------
    name : str
        Canonical, case-insensitive format name.
    reader : callable
        Callable accepting the source as its first argument.
    writer : callable, optional
        Callable accepting ``(object, target)`` for serialization.
    extensions : tuple of str, optional
        Filename extensions used only as hints. Detection remains content
        based whenever a detector is available.
    detector : callable, optional
        Callable returning ``True`` when the source belongs to this format.
    aliases : tuple of str, optional
        Alternative names accepted by :func:`get_tf_format`.
    description : str, optional
        Human-readable summary.
    """

    name: str
    reader: TFReader
    writer: TFWriter | None = None
    extensions: tuple[str, ...] = ()
    detector: TFDetector | None = None
    aliases: tuple[str, ...] = ()
    description: str = ""


_registry: NamedRegistry[TransferFunctionFormat] = NamedRegistry(
    kind="transfer-function format",
    error_cls=ValueError,
)


def register_tf_format(
    name: str,
    *,
    reader: TFReader,
    writer: TFWriter | None = None,
    extensions: tuple[str, ...] | list[str] = (),
    detector: TFDetector | None = None,
    aliases: tuple[str, ...] | list[str] = (),
    description: str = "",
    replace: bool = False,
) -> TransferFunctionFormat:
    """Register one transfer-function reader definition.

    Registration is deliberately small and additive. A format may be
    readable without being writable; write support is advertised explicitly
    through the optional ``writer`` callable.
    """
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

    definition = TransferFunctionFormat(
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


def get_tf_format(name: str) -> TransferFunctionFormat:
    """Return a registered format by canonical name or alias."""
    result = _registry.get(name)
    if result is None:
        available = ", ".join(_registry.names()) or "<none>"
        raise TransferFunctionFormatError(
            f"unknown transfer-function format {name!r}; "
            f"available: {available}"
        )
    return result


def list_tf_formats() -> dict[str, dict[str, Any]]:
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


def _source_probe(source: Any, *, limit: int = 65536) -> str:
    """Return a text probe without consuming seekable streams."""
    if isinstance(source, bytes):
        return source[:limit].decode("utf-8-sig", errors="replace")

    if isinstance(source, str):
        # XML text and other inline payloads should never be interpreted as
        # enormous filenames. Newlines or markup strongly indicate content.
        if "\n" in source or "\r" in source or "<" in source:
            return source[:limit]
        try:
            path = Path(source)
            if path.exists() and path.is_file():
                with path.open("rb") as stream:
                    return stream.read(limit).decode(
                        "utf-8-sig", errors="replace"
                    )
        except (OSError, ValueError):
            pass
        return source[:limit]

    if isinstance(source, PathLike):
        path = Path(source)
        with path.open("rb") as stream:
            return stream.read(limit).decode("utf-8-sig", errors="replace")

    read = getattr(source, "read", None)
    if callable(read):
        tell = getattr(source, "tell", None)
        seek = getattr(source, "seek", None)
        position = None
        if callable(tell) and callable(seek):
            try:
                position = tell()
            except Exception:
                position = None
        raw = read(limit)
        if position is not None:
            try:
                seek(position)
            except Exception:
                pass
        if isinstance(raw, bytes):
            return raw.decode("utf-8-sig", errors="replace")
        return str(raw)

    raise TypeError(
        "source must be a path, inline text/bytes, or readable file object"
    )


def _is_emtf_xml(source: Any) -> bool:
    probe = _source_probe(source).lstrip("\ufeff \t\r\n")
    # Allow an XML declaration, comments, and a namespace prefix on EM_TF.
    return bool(
        re.search(
            r"<(?:[A-Za-z_][\w.\-]*:)?EM_TF(?:\s|>)",
            probe,
            flags=re.IGNORECASE,
        )
    )


def _is_edi(source: Any) -> bool:
    probe = _source_probe(source)
    first_tag = None
    for line in probe.splitlines():
        stripped = line.lstrip("\ufeff \t")
        if not stripped or stripped.startswith("#"):
            continue
        if stripped.startswith(">"):
            first_tag = stripped.upper()
            break
    if first_tag is None or not first_tag.startswith(">HEAD"):
        return False
    # A second structural signal reduces accidental classification of a text
    # file whose first line happens to be '>HEAD'.
    upper = probe.upper()
    structural = (
        ">=MTSECT",
        ">=EMAPSECT",
        ">=SPECTRASECT",
        ">=TSERIESSECT",
        ">=DEFINEMEAS",
        ">FREQ",
    )
    return any(token in upper for token in structural)


def _read_emtf_xml(source: Any, **kwargs: Any) -> Any:
    # Lazy import keeps ``import pycsamt.io`` independent of EMTF internals.
    from ..emtf.xml import EMTFXMLReader

    return EMTFXMLReader(**kwargs).read(source)


def _read_edi(source: Any, **kwargs: Any) -> Any:
    # EDIFile continues to own path-oriented historical parsing. Conversion
    # to the neutral EMTF model is explicit through pycsamt.emtf.from_edi().
    if not isinstance(source, (str, PathLike)):
        raise TypeError(
            "EDI reading currently requires a filesystem path; "
            "in-memory EDI parsing is not provided by the historical backend"
        )
    if isinstance(source, str) and (
        "\n" in source or "\r" in source or source.lstrip().startswith(">")
    ):
        raise TypeError(
            "EDI reading currently requires a filesystem path; "
            "in-memory EDI parsing is not provided by the historical backend"
        )
    from ..seg.edi import EDIFile

    return EDIFile(source, **kwargs)


def _write_edi(obj: Any, target: Any, **kwargs: Any) -> Any:
    # Keep the registry serialization-only and import the adapter lazily.
    from ..emtf.converters.edi import write_edi

    return write_edi(obj, target, **kwargs)


def _write_emtf_xml(obj: Any, target: Any, **kwargs: Any) -> Any:
    # Lazy import keeps ``import pycsamt.io`` independent of XML internals.
    from ..emtf.xml import EMTFXMLWriter

    return EMTFXMLWriter(
        strict=kwargs.pop("strict", True),
        precision=kwargs.pop("precision", 17),
        pretty=kwargs.pop("pretty", True),
    ).write(obj, target, **kwargs)


def get_tf_format_for_target(target: Any) -> TransferFunctionFormat:
    """Resolve a writable format from a target filename extension.

    Content detection cannot be used for a file that does not yet exist.
    Therefore writing uses explicit format selection or a unique registered
    extension. Streams require an explicit format at the higher-level API.
    """
    if not isinstance(target, (str, PathLike)):
        raise TransferFunctionFormatError(
            "a writable stream requires an explicit transfer-function format"
        )
    try:
        suffix = Path(target).suffix.lower()
    except (OSError, ValueError) as exc:
        raise TransferFunctionFormatError(
            "cannot infer output format from target"
        ) from exc
    if not suffix:
        raise TransferFunctionFormatError(
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
        raise TransferFunctionFormatError(
            f"no writable transfer-function format registered for {suffix!r}"
        )
    raise TransferFunctionFormatError(
        f"ambiguous output extension {suffix!r}; specify format explicitly"
    )


def detect_tf_format(source: Any) -> str:
    """Detect the canonical transfer-function format from source content.

    Detection is content-first. Filename extensions are intentionally not
    accepted as proof of a format; they are metadata hints stored in the
    registry for interfaces and future writer dispatch.
    """
    matches = _registry.match_by_detector(source)

    if len(matches) == 1:
        return matches[0]
    if len(matches) > 1:
        raise TransferFunctionFormatError(
            "ambiguous transfer-function format: " + ", ".join(matches)
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
    raise TransferFunctionFormatError(
        "unable to detect a supported transfer-function format from "
        f"content.{hint}"
    )


# Built-in format definitions. Reader imports remain lazy by design.
register_tf_format(
    "edi",
    reader=_read_edi,
    writer=_write_edi,
    extensions=(".edi",),
    detector=_is_edi,
    aliases=("seg_edi",),
    description="SEG Electrical Data Interchange transfer functions",
)
register_tf_format(
    "emtf_xml",
    reader=_read_emtf_xml,
    writer=_write_emtf_xml,
    extensions=(".xml",),
    detector=_is_emtf_xml,
    aliases=("xml", "emtf"),
    description="EMTF XML electromagnetic transfer functions",
)
