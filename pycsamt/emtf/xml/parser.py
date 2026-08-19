# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

"""Low-level, namespace-safe helpers for EMTF XML parsing."""

from __future__ import annotations

from dataclasses import dataclass
from os import PathLike
from pathlib import Path
from typing import Any, IO
import warnings
import xml.etree.ElementTree as ET

__all__ = [
    "EMTFXMLParseError",
    "EMTFXMLWarning",
    "XMLDataTypeSpec",
    "child",
    "children",
    "element_to_mapping",
    "local_name",
    "load_xml_root",
    "parse_numeric_text",
    "text",
]


class EMTFXMLParseError(ValueError):
    """Raised when an EMTF XML document cannot be interpreted safely."""


class EMTFXMLWarning(UserWarning):
    """Warning emitted for recoverable EMTF XML inconsistencies."""


@dataclass(frozen=True)
class XMLDataTypeSpec:
    """Per-document data-type declaration used while parsing XML.

    This is deliberately separate from the global scientific registry. EMTF
    XML files are self-describing and may contain data types unknown to the
    installed pyCSAMT version; reading one file must not mutate global state.
    """

    name: str
    tag: str
    data_kind: str
    input_kind: str | None = None
    output_kind: str | None = None
    units: str | None = None
    intention: str = "primary"
    description: str = ""
    derived_from: str | None = None
    see_also: tuple[str, ...] = ()
    external_url: str | None = None

    @property
    def is_scalar(self) -> bool:
        return not self.input_kind and not self.output_kind


def local_name(tag: str) -> str:
    """Return an XML local name, discarding namespace or prefix."""
    raw = str(tag)
    if "}" in raw:
        raw = raw.rsplit("}", 1)[-1]
    if ":" in raw:
        raw = raw.rsplit(":", 1)[-1]
    return raw


def children(node: ET.Element | None, name: str) -> list[ET.Element]:
    """Return direct child elements matching *name* by local name."""
    if node is None:
        return []
    return [item for item in list(node) if local_name(item.tag) == name]


def child(node: ET.Element | None, name: str) -> ET.Element | None:
    """Return the first direct child matching *name*."""
    matches = children(node, name)
    return matches[0] if matches else None


def text(node: ET.Element | None, name: str | None = None) -> str | None:
    """Return stripped text from *node* or one named direct child."""
    target = child(node, name) if name is not None else node
    if target is None or target.text is None:
        return None
    value = target.text.strip()
    return value or None


def _coerce_source(
    source: str | bytes | PathLike[str] | IO[str] | IO[bytes],
) -> tuple[ET.Element, str]:
    if hasattr(source, "read"):
        try:
            tree = ET.parse(source)  # type: ignore[arg-type]
        except (ET.ParseError, OSError, ValueError) as exc:
            raise EMTFXMLParseError(f"invalid XML stream: {exc}") from exc
        return tree.getroot(), getattr(source, "name", "<stream>")

    if isinstance(source, bytes):
        try:
            return ET.fromstring(source), "<bytes>"
        except ET.ParseError as exc:
            raise EMTFXMLParseError(f"invalid XML bytes: {exc}") from exc

    raw = str(source)
    if raw.lstrip().startswith("<"):
        try:
            return ET.fromstring(raw), "<string>"
        except ET.ParseError as exc:
            raise EMTFXMLParseError(f"invalid XML string: {exc}") from exc

    path = Path(raw).expanduser()
    try:
        tree = ET.parse(path)
    except (ET.ParseError, OSError, ValueError) as exc:
        raise EMTFXMLParseError(f"cannot parse {path}: {exc}") from exc
    return tree.getroot(), str(path)


def load_xml_root(
    source: str | bytes | PathLike[str] | IO[str] | IO[bytes],
) -> tuple[ET.Element, str]:
    """Parse *source* and return ``(root, source_label)``."""
    return _coerce_source(source)


def parse_numeric_text(
    value: str | None,
    *,
    complex_: bool,
) -> complex | float:
    """Parse an FCU-style real or ``real imag`` numeric value."""
    if value is None:
        raise EMTFXMLParseError("numeric XML value is empty")
    normalized = value.replace("D", "E").replace("d", "e").strip()
    parts = normalized.replace(",", " ").split()
    if complex_:
        if len(parts) == 2:
            try:
                return complex(float(parts[0]), float(parts[1]))
            except ValueError as exc:
                raise EMTFXMLParseError(
                    f"invalid complex value: {value!r}"
                ) from exc
        if len(parts) == 1:
            token = parts[0].replace("i", "j")
            try:
                return complex(token)
            except ValueError as exc:
                raise EMTFXMLParseError(
                    f"invalid complex value: {value!r}"
                ) from exc
        raise EMTFXMLParseError(f"invalid complex value: {value!r}")

    if len(parts) != 1:
        raise EMTFXMLParseError(f"invalid real value: {value!r}")
    try:
        return float(parts[0])
    except ValueError as exc:
        raise EMTFXMLParseError(f"invalid real value: {value!r}") from exc


def element_to_mapping(node: ET.Element) -> dict[str, Any]:
    """Convert a small metadata subtree to a loss-preserving Python mapping."""
    result: dict[str, Any] = {}
    if node.attrib:
        result["@attributes"] = dict(node.attrib)
    raw_text = (node.text or "").strip()
    if raw_text:
        result["#text"] = raw_text
    for item in list(node):
        key = local_name(item.tag)
        value: Any
        if list(item) or item.attrib:
            value = element_to_mapping(item)
        else:
            value = (item.text or "").strip()
        if key in result:
            current = result[key]
            if not isinstance(current, list):
                current = [current]
            current.append(value)
            result[key] = current
        else:
            result[key] = value
    return result


def warn(message: str) -> None:
    """Emit a reader warning with a stable warning category."""
    warnings.warn(message, EMTFXMLWarning, stacklevel=3)
