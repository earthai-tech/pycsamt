# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

"""EMTF scientific objects -> deterministic XML element trees.

The serializer is intentionally separate from filesystem I/O.  It translates
format-neutral :class:`pycsamt.emtf.EMTF` objects into the EMTF XML vocabulary
without routing through EDI or exposing ``ElementTree`` nodes to the scientific
model.
"""

from __future__ import annotations

from datetime import datetime
import re
from typing import Any, Iterable
import warnings
import xml.etree.ElementTree as ET

import numpy as np

from ...api.property import PyCSAMTObject
from ...metadata import ChannelMeta, Person
from ..document import EMTF
from ..estimates import StatisticalEstimate
from ..transfer import TransferFunction
from .constants import EMTF_ROOT, ESTIMATE_CODES

__all__ = [
    "EMTFXMLSerializationError",
    "EMTFXMLWriteWarning",
    "EMTFXMLSerializer",
]


_XML_NAME = re.compile(r"^[A-Za-z_][A-Za-z0-9_.-]*$")

_ESTIMATE_DEFAULTS = {
    "VAR": {
        "name": "VAR",
        "type": "real",
        "description": "Variance",
        "intention": "error estimate",
        "tag": "variance",
        "external_url": None,
    },
    "INVSIGCOV": {
        "name": "INVSIGCOV",
        "type": "complex",
        "description": "Inverse Coherent Signal Power Matrix (S)",
        "intention": "signal power estimate",
        "tag": "inverse_signal_covariance",
        "external_url": None,
    },
    "RESIDCOV": {
        "name": "RESIDCOV",
        "type": "complex",
        "description": "Residual Covariance (N)",
        "intention": "error estimate",
        "tag": "residual_covariance",
        "external_url": None,
    },
}


class EMTFXMLSerializationError(ValueError):
    """Raised when an EMTF object cannot be serialized safely."""


class EMTFXMLWriteWarning(UserWarning):
    """Warning emitted when permissive XML writing omits unsafe content."""


def _warn(message: str) -> None:
    warnings.warn(message, EMTFXMLWriteWarning, stacklevel=3)


def _value_text(value: Any) -> str:
    if isinstance(value, datetime):
        return value.isoformat()
    return str(value)


def _append_text(
    parent: ET.Element,
    tag: str,
    value: Any,
    *,
    allow_empty: bool = False,
    attrs: dict[str, Any] | None = None,
) -> ET.Element | None:
    if value is None:
        return None
    raw = _value_text(value)
    if not raw and not allow_empty:
        return None
    node = ET.SubElement(parent, tag)
    for key, attr_value in (attrs or {}).items():
        if attr_value is not None:
            node.set(str(key), _value_text(attr_value))
    if raw:
        node.text = raw
    return node


def _apply_mapping(node: ET.Element, value: Any) -> None:
    """Reverse the reader's ``element_to_mapping`` representation."""
    if isinstance(value, dict):
        attributes = value.get("@attributes", {})
        if isinstance(attributes, dict):
            for key, attr_value in attributes.items():
                if attr_value is not None:
                    node.set(str(key), _value_text(attr_value))
        if "#text" in value and value["#text"] is not None:
            node.text = _value_text(value["#text"])
        for key, child_value in value.items():
            if key in {"@attributes", "#text"}:
                continue
            values = (
                child_value
                if isinstance(child_value, list)
                else [child_value]
            )
            for item in values:
                child_node = ET.SubElement(node, str(key))
                _apply_mapping(child_node, item)
        return
    if value is not None:
        node.text = _value_text(value)


def _append_mapping(parent: ET.Element, tag: str, value: Any) -> ET.Element:
    node = ET.SubElement(parent, tag)
    _apply_mapping(node, value)
    return node


class EMTFXMLSerializer(PyCSAMTObject):
    """Build a deterministic EMTF XML element tree from an :class:`EMTF`.

    Parameters
    ----------
    strict : bool, default=True
        Reject unsupported or scientifically ambiguous content.  Permissive
        mode warns and omits only the problematic component.
    precision : int, default=17
        Significant decimal digits used for floating-point response data.
        Seventeen digits are sufficient to round-trip IEEE float64 values.
    """

    def __init__(self, *, strict: bool = True, precision: int = 17) -> None:
        self.strict = bool(strict)
        self.precision = int(precision)
        if not 6 <= self.precision <= 20:
            raise ValueError("precision must be between 6 and 20")

    def to_element(self, document: EMTF) -> ET.Element:
        """Return the root ``<EM_TF>`` element for *document*."""
        if not isinstance(document, EMTF):
            raise TypeError("document must be a pycsamt.emtf.EMTF")

        root = ET.Element(EMTF_ROOT)
        self._append_header(root, document)
        self._append_provenance(root, document)
        self._append_copyright(root, document)
        self._append_site(root, document)
        self._append_field_notes(root, document)
        self._append_processing(root, document)

        data_specs = self._data_type_specs(document)
        estimate_specs = self._estimate_specs(document)
        self._append_estimate_declarations(root, estimate_specs)
        self._append_data_type_declarations(root, data_specs)
        self._append_site_layout(root, document)
        self._append_data(root, document, data_specs)
        self._append_period_range(root, document)
        return root

    # ------------------------------------------------------------------
    # Header / metadata
    # ------------------------------------------------------------------

    def _append_header(self, root: ET.Element, document: EMTF) -> None:
        _append_text(root, "Description", document.description)
        _append_text(root, "ProductId", document.product_id)
        _append_text(root, "SubType", document.subtype)

        metadata = document.metadata or {}
        if "notes" in metadata:
            _append_mapping(root, "Notes", metadata["notes"])

        if document.tags:
            _append_text(root, "Tags", ",".join(document.tags))

        for key, tag in (
            ("externalurl", "ExternalUrl"),
            ("primarydata", "PrimaryData"),
            ("attachment", "Attachment"),
            ("gridorigin", "GridOrigin"),
        ):
            if key not in metadata:
                continue
            raw = metadata[key]
            values = raw if isinstance(raw, list) else [raw]
            for value in values:
                _append_mapping(root, tag, value)

    def _append_provenance(self, root: ET.Element, document: EMTF) -> None:
        provenance = document.provenance
        if provenance is None:
            return
        node = ET.SubElement(root, "Provenance")
        _append_text(node, "CreateTime", provenance.create_time)
        _append_text(
            node,
            "CreatingApplication",
            provenance.creating_application,
        )
        self._append_person(node, "Creator", provenance.creator)
        self._append_person(node, "Submitter", provenance.submitter)

    @staticmethod
    def _append_person(
        parent: ET.Element,
        tag: str,
        person: Person | None,
    ) -> None:
        if person is None:
            return
        if not any(
            (
                person.name,
                person.email,
                person.organization,
                person.organization_url,
            )
        ):
            return
        node = ET.SubElement(parent, tag)
        _append_text(node, "Name", person.name)
        _append_text(node, "Email", person.email)
        _append_text(node, "Org", person.organization)
        _append_text(node, "OrgUrl", person.organization_url)

    def _append_copyright(self, root: ET.Element, document: EMTF) -> None:
        info = document.copyright
        orientation = document.orientation
        rotation_info = (
            orientation.rotation_info if orientation is not None else None
        )
        if info is None and rotation_info is None:
            return
        if info is None:
            # Rotation history is scientifically useful but must not force the
            # creation of an otherwise fabricated Copyright metadata block.
            _append_text(root, "RotationInfo", rotation_info)
            return

        node = ET.SubElement(root, "Copyright")
        if info is not None:
            reference = info.reference
            citation = ET.SubElement(node, "Citation")
            _append_text(citation, "Title", reference.title)
            _append_text(citation, "Authors", reference.author)
            _append_text(citation, "Year", reference.year)
            survey_doi = reference.extra.get("survey_doi")
            _append_text(citation, "SurveyDOI", survey_doi)
            doi = reference.doi or reference.extra.get("doi")
            _append_text(citation, "DOI", doi)

            extras = reference.extra
            _append_text(
                node,
                "SelectedPublications",
                extras.get("selectedpublications"),
            )
            _append_text(
                node,
                "Acknowledgement",
                extras.get("acknowledgement"),
            )
            _append_text(node, "ReleaseStatus", info.release_status)
            _append_text(node, "ConditionsOfUse", info.conditions_of_use)
        if rotation_info:
            _append_text(node, "RotationInfo", rotation_info)
        if info is not None:
            _append_text(
                node,
                "AdditionalInfo",
                info.reference.extra.get("additionalinfo"),
            )

    def _append_site(self, root: ET.Element, document: EMTF) -> None:
        site = document.site
        orientation = document.orientation
        quality = document.quality
        if site is None and orientation is None and quality is None:
            return

        node = ET.SubElement(root, "Site")
        if site is not None:
            _append_text(node, "Project", site.project)
            _append_text(node, "Survey", site.survey)
            _append_text(node, "YearCollected", site.year_collected)
            _append_text(node, "Country", site.country)
            _append_text(node, "Id", site.site_id)
            _append_text(node, "Name", site.name)
            self._append_location(node, site.location)

        if orientation is not None and orientation.mode is not None:
            attrs: dict[str, Any] = {}
            if orientation.is_orthogonal:
                attrs["angle_to_geographic_north"] = self._format_float(
                    orientation.angle_to_geographic_north
                )
            _append_text(
                node,
                "Orientation",
                orientation.mode,
                attrs=attrs,
            )

        if site is not None:
            _append_text(node, "AcquiredBy", site.acquired_by)
            _append_text(node, "Start", site.start)
            _append_text(node, "End", site.end)
            run_list = site.extra.get("run_list")
            if isinstance(run_list, (tuple, list)):
                run_list = " ".join(str(item) for item in run_list)
            _append_text(node, "RunList", run_list)

        self._append_quality(node, quality)
        if site is not None:
            comments = site.extra.get("comments", [])
            if not isinstance(comments, list):
                comments = [comments]
            for comment in comments:
                _append_text(node, "Comments", comment)

    def _append_location(self, parent: ET.Element, location: Any) -> None:
        if location is None:
            return
        attrs: dict[str, str] = {}
        if location.datum:
            attrs["datum"] = str(location.datum)
        node = ET.SubElement(parent, "Location", attrs)
        _append_text(
            node,
            "Latitude",
            self._format_float(location.latitude),
        )
        _append_text(
            node,
            "Longitude",
            self._format_float(location.longitude),
        )
        if location.elevation is not None:
            _append_text(
                node,
                "Elevation",
                self._format_float(location.elevation),
                attrs={"units": location.elevation_units},
            )
        if location.declination is not None:
            attrs = {}
            if location.declination_epoch is not None:
                attrs["epoch"] = self._format_float(
                    location.declination_epoch
                )
            _append_text(
                node,
                "Declination",
                self._format_float(location.declination),
                attrs=attrs,
            )

    def _append_quality(self, parent: ET.Element, quality: Any) -> None:
        if quality is None:
            return
        if any(
            value is not None
            for value in (
                quality.rating,
                quality.good_from_period,
                quality.good_to_period,
            )
        ) or quality.comments:
            node = ET.SubElement(parent, "DataQualityNotes")
            _append_text(node, "Rating", quality.rating)
            _append_text(
                node,
                "GoodFromPeriod",
                self._format_float(quality.good_from_period),
            )
            _append_text(
                node,
                "GoodToPeriod",
                self._format_float(quality.good_to_period),
            )
            for comment in quality.comments:
                attrs = {"author": comment.author} if comment.author else {}
                _append_text(node, "Comments", comment.text, attrs=attrs)

        if quality.warning_flag is not None or quality.warnings:
            node = ET.SubElement(parent, "DataQualityWarnings")
            _append_text(node, "Flag", quality.warning_flag)
            for comment in quality.warnings:
                attrs = {"author": comment.author} if comment.author else {}
                _append_text(node, "Comments", comment.text, attrs=attrs)

    def _append_field_notes(self, root: ET.Element, document: EMTF) -> None:
        for key, raw in document.field_notes.items():
            values = raw if isinstance(raw, list) else [raw]
            for value in values:
                node = ET.SubElement(root, "FieldNotes")
                _apply_mapping(node, value)
                if "run" not in node.attrib and key and not key.startswith(
                    "run_"
                ):
                    node.set("run", str(key))

    def _append_processing(self, root: ET.Element, document: EMTF) -> None:
        processing = document.processing
        if processing is None:
            return
        node = ET.SubElement(root, "ProcessingInfo")
        _append_text(
            node,
            "SignConvention",
            self._xml_sign_convention(processing.sign_convention),
        )
        remote = processing.remote_reference
        if remote is not None and remote.reference_type is not None:
            ET.SubElement(
                node,
                "RemoteRef",
                {"type": str(remote.reference_type)},
            )

        remote_info = processing.extra.get("remote_info")
        if remote_info is not None:
            _append_mapping(node, "RemoteInfo", remote_info)
        elif remote is not None and remote.site is not None:
            info = ET.SubElement(node, "RemoteInfo")
            site = ET.SubElement(info, "Site")
            _append_text(site, "Id", remote.site)

        _append_text(node, "ProcessedBy", processing.processed_by)
        _append_text(
            node,
            "ProcessDate",
            processing.extra.get("process_date"),
        )
        if processing.software is not None:
            software = ET.SubElement(node, "ProcessingSoftware")
            _append_text(software, "Name", processing.software.name)
            _append_text(software, "LastMod", processing.software.release)
            author = processing.software.author
            _append_text(
                software,
                "Author",
                author.name if author is not None else None,
            )
        _append_text(node, "ProcessingTag", processing.processing_tag)

    @staticmethod
    def _xml_sign_convention(value: str | None) -> str | None:
        if value == "exp(+i ω t)":
            return r"exp(+ i\omega t)"
        if value == "exp(-i ω t)":
            return r"exp(- i\omega t)"
        return value

    # ------------------------------------------------------------------
    # Declarations
    # ------------------------------------------------------------------

    def _data_type_specs(self, document: EMTF) -> list[dict[str, Any]]:
        stored = document.metadata.get("xml_data_types", [])
        stored_specs = [
            dict(item) for item in stored if isinstance(item, dict)
        ]
        by_tag = {
            str(item.get("tag", "")).lower(): item
            for item in stored_specs
            if item.get("tag")
        }
        by_name = {
            str(item.get("name", "")).upper(): item
            for item in stored_specs
            if item.get("name")
        }

        result = list(stored_specs)
        known_keys = {
            (
                str(item.get("name", "")).upper(),
                str(item.get("tag", "")).lower(),
            )
            for item in result
        }
        for tf in document.transfer_functions.values():
            xml_name = str(tf.attrs.get("xml_name") or "").upper()
            source = by_tag.get(tf.name.lower()) or by_name.get(xml_name)
            spec = (
                dict(source) if source is not None else self._spec_from_tf(tf)
            )
            key = (
                str(spec.get("name", "")).upper(),
                str(spec.get("tag", "")).lower(),
            )
            if key not in known_keys:
                result.append(spec)
                known_keys.add(key)
        return result

    def _spec_from_tf(self, tf: TransferFunction) -> dict[str, Any]:
        definition = tf.definition
        xml_name = tf.attrs.get("xml_name")
        if not xml_name and definition is not None:
            xml_name = definition.name
        if not xml_name:
            candidate = tf.name.upper()
            if _XML_NAME.match(candidate):
                xml_name = candidate
        if not xml_name or not _XML_NAME.match(str(xml_name)):
            self._problem(
                f"transfer function {tf.name!r} has no safe EMTF XML code"
            )
            xml_name = "UNKNOWN"

        data_kind = tf.attrs.get("xml_data_kind")
        if not data_kind:
            if definition is not None:
                data_kind = definition.data_kind
            else:
                data_kind = "complex" if np.iscomplexobj(tf.data) else "real"

        return {
            "name": str(xml_name).upper(),
            "tag": tf.name,
            "data_kind": str(data_kind).lower(),
            "input_kind": tf.attrs.get("xml_input_kind")
            or (definition.input_kind if definition is not None else None),
            "output_kind": tf.attrs.get("xml_output_kind")
            or (definition.output_kind if definition is not None else None),
            "units": tf.units
            if tf.units is not None
            else (definition.units if definition is not None else None),
            "intention": tf.attrs.get("xml_intention")
            or (definition.intention if definition is not None else "primary"),
            "description": definition.description if definition else "",
            "derived_from": definition.derived_from if definition else None,
            "see_also": definition.see_also if definition else (),
            "external_url": None,
        }

    def _estimate_specs(self, document: EMTF) -> list[dict[str, Any]]:
        stored = document.metadata.get("xml_statistical_estimates", [])
        result = [dict(item) for item in stored if isinstance(item, dict)]
        known = {
            str(item.get("name", "")).upper()
            for item in result
            if item.get("name")
        }
        for tf in document.transfer_functions.values():
            for estimate in tf.estimates.values():
                code = estimate.name.upper()
                if code not in ESTIMATE_CODES:
                    self._problem(
                        "EMTF XML writer does not yet define a safe mapping "
                        f"for statistical estimate {code!r} on {tf.name!r}"
                    )
                    continue
                if code not in known:
                    result.append(dict(_ESTIMATE_DEFAULTS[code]))
                    known.add(code)
        return result

    def _append_estimate_declarations(
        self,
        root: ET.Element,
        specs: list[dict[str, Any]],
    ) -> None:
        if not specs:
            return
        parent = ET.SubElement(root, "StatisticalEstimates")
        for spec in specs:
            name = spec.get("name")
            if not name:
                continue
            attrs = {"name": str(name)}
            kind = spec.get("type") or spec.get("data_kind")
            if kind:
                attrs["type"] = str(kind)
            node = ET.SubElement(parent, "Estimate", attrs)
            _append_text(node, "Description", spec.get("description"))
            _append_text(node, "Intention", spec.get("intention"))
            _append_text(node, "Tag", spec.get("tag"))
            _append_text(node, "ExternalUrl", spec.get("external_url"))

    def _append_data_type_declarations(
        self,
        root: ET.Element,
        specs: list[dict[str, Any]],
    ) -> None:
        if not specs:
            return
        parent = ET.SubElement(root, "DataTypes")
        for spec in specs:
            name = spec.get("name")
            tag = spec.get("tag")
            if not name or not tag:
                self._problem("EMTF DataType declaration needs name and tag")
                continue
            attrs: dict[str, str] = {"name": str(name)}
            data_kind = spec.get("data_kind") or spec.get("type")
            if data_kind:
                attrs["type"] = str(data_kind)
            if spec.get("output_kind"):
                attrs["output"] = str(spec["output_kind"])
            elif spec.get("output"):
                attrs["output"] = str(spec["output"])
            if spec.get("input_kind"):
                attrs["input"] = str(spec["input_kind"])
            elif spec.get("input"):
                attrs["input"] = str(spec["input"])
            if spec.get("units"):
                attrs["units"] = str(spec["units"])
            node = ET.SubElement(parent, "DataType", attrs)
            _append_text(node, "Description", spec.get("description"))
            intention = str(spec.get("intention") or "primary")
            if intention in {"primary", "derived"}:
                intention = f"{intention} data type"
            _append_text(node, "Intention", intention)
            _append_text(node, "Tag", tag)
            _append_text(node, "DerivedFrom", spec.get("derived_from"))
            see_also = spec.get("see_also")
            if isinstance(see_also, (tuple, list)):
                see_also = ",".join(str(item) for item in see_also if item)
            _append_text(node, "SeeAlso", see_also)
            _append_text(node, "ExternalUrl", spec.get("external_url"))

    # ------------------------------------------------------------------
    # Site layout
    # ------------------------------------------------------------------

    def _append_site_layout(self, root: ET.Element, document: EMTF) -> None:
        layout = document.site_layout
        if layout is None:
            return
        parent = ET.SubElement(root, "SiteLayout")
        self._append_channel_group(
            parent,
            "InputChannels",
            layout.input_channels,
            units=layout.input_units,
            reference=layout.input_reference,
        )
        self._append_channel_group(
            parent,
            "OutputChannels",
            layout.output_channels,
            units=layout.output_units,
            reference=layout.output_reference,
        )

    def _append_channel_group(
        self,
        parent: ET.Element,
        tag: str,
        channels: Iterable[ChannelMeta],
        *,
        units: str | None,
        reference: str | None,
    ) -> None:
        channel_list = list(channels)
        if not channel_list:
            return
        attrs: dict[str, str] = {}
        if reference:
            attrs["ref"] = reference
        if units:
            attrs["units"] = units
        node = ET.SubElement(parent, tag, attrs)
        for channel in channel_list:
            tag_name = (
                "Electric"
                if channel.is_electric
                else "Magnetic"
                if channel.is_magnetic
                else "Channel"
            )
            values: dict[str, Any] = {
                "name": channel.name,
                "orientation": channel.orientation,
                "tilt": channel.tilt,
                "x": channel.x,
                "y": channel.y,
                "z": channel.z,
                "x2": channel.x2,
                "y2": channel.y2,
                "z2": channel.z2,
                "units": channel.units,
                "ref": channel.reference,
                "id": channel.sensor_id,
            }
            values.update(channel.extra)
            attrs = {
                key: self._format_float(value)
                if isinstance(value, (float, int, np.floating, np.integer))
                else str(value)
                for key, value in values.items()
                if value is not None
            }
            ET.SubElement(node, tag_name, attrs)

    # ------------------------------------------------------------------
    # Data matrices
    # ------------------------------------------------------------------

    def _append_data(
        self,
        root: ET.Element,
        document: EMTF,
        specs: list[dict[str, Any]],
    ) -> None:
        if document.is_empty():
            return
        periods = document.periods
        if periods is None:
            self._problem(
                "cannot serialize transfer-function data without periods"
            )
            return
        periods = np.asarray(periods, dtype=float)
        if periods.size != document.n_periods:
            self._problem("document period count is inconsistent")
            return

        spec_by_tag = {
            str(item.get("tag", "")).lower(): item
            for item in specs
            if item.get("tag")
        }
        parent = ET.SubElement(root, "Data", {"count": str(periods.size)})
        for iper, period in enumerate(periods):
            pnode = ET.SubElement(
                parent,
                "Period",
                {"value": self._format_period_number(period), "units": "secs"},
            )
            for tf in document.transfer_functions.values():
                spec = spec_by_tag.get(tf.name.lower())
                if spec is None:
                    spec = self._spec_from_tf(tf)
                self._append_tf_period(pnode, tf, spec, iper)

    def _append_tf_period(
        self,
        parent: ET.Element,
        tf: TransferFunction,
        spec: dict[str, Any],
        iper: int,
    ) -> None:
        code = str(spec.get("name") or "").upper()
        if not code or not _XML_NAME.match(code):
            self._problem(f"invalid XML data-type name {code!r}")
            return

        data_kind = str(spec.get("data_kind") or "complex").lower()
        self._append_matrix(
            parent,
            code,
            tf.data[iper],
            rows=tf.output_channels,
            cols=tf.input_channels,
            complex_=data_kind == "complex",
            units=tf.units or spec.get("units"),
            component_code=code,
            scalar=tf.n_output == 1
            and tf.n_input == 1
            and not tf.output_channels
            and not tf.input_channels,
            declared_size=(tf.n_output, tf.n_input),
        )

        for estimate in tf.estimates.values():
            est_code = estimate.name.upper()
            if est_code not in ESTIMATE_CODES:
                # The declaration path already issued a diagnostic.  Avoid a
                # second warning in permissive mode.
                continue
            matrix = self._estimate_matrix(tf, estimate, est_code, iper)
            if matrix is None:
                continue
            if est_code == "VAR":
                rows = tf.output_channels
                cols = tf.input_channels
                complex_ = False
                scalar = not rows and not cols
            elif est_code == "INVSIGCOV":
                rows = tf.input_channels
                cols = tf.input_channels
                complex_ = True
                scalar = False
            else:
                rows = tf.output_channels
                cols = tf.output_channels
                complex_ = True
                scalar = False
            self._append_matrix(
                parent,
                f"{code}.{est_code}",
                matrix,
                rows=rows,
                cols=cols,
                complex_=complex_,
                units=estimate.units,
                component_code=code,
                scalar=scalar,
                # FCU v4.1 writes the parent TF size for covariance blocks.
                # The value channel attributes remain the source of truth.
                declared_size=(tf.n_output, tf.n_input),
            )

    def _estimate_matrix(
        self,
        tf: TransferFunction,
        estimate: StatisticalEstimate,
        code: str,
        iper: int,
    ) -> np.ndarray | None:
        data = np.asarray(estimate.data)
        if data.ndim != 3:
            self._problem(
                f"estimate {code} on {tf.name!r} must be a 3-D array"
            )
            return None
        if code == "VAR":
            expected = (tf.n_periods, tf.n_output, tf.n_input)
        elif code == "INVSIGCOV":
            expected = (tf.n_periods, tf.n_input, tf.n_input)
        else:
            expected = (tf.n_periods, tf.n_output, tf.n_output)
        if data.shape != expected:
            self._problem(
                f"estimate {code} on {tf.name!r} has shape {data.shape}; "
                f"expected {expected}"
            )
            return None
        if code == "VAR" and np.iscomplexobj(data):
            imag = np.imag(data)
            if np.any(np.isfinite(imag) & (imag != 0.0)):
                self._problem(
                    f"variance estimate on {tf.name!r} contains non-zero "
                    "imaginary values"
                )
                return None
            data = np.real(data)
        return np.asarray(data[iper])

    def _append_matrix(
        self,
        parent: ET.Element,
        tag: str,
        matrix: np.ndarray,
        *,
        rows: tuple[str, ...],
        cols: tuple[str, ...],
        complex_: bool,
        units: str | None,
        component_code: str,
        scalar: bool,
        declared_size: tuple[int, int],
    ) -> None:
        matrix = np.asarray(matrix)
        values: list[tuple[int, int, complex | float]] = []
        for row in range(matrix.shape[0]):
            for col in range(matrix.shape[1]):
                raw = matrix[row, col]
                if not self._is_finite(raw, complex_=complex_):
                    if self._is_missing(raw):
                        continue
                    self._problem(
                        f"{tag} contains a non-finite value at "
                        f"[{row}, {col}]"
                    )
                    continue
                values.append((row, col, raw))
        if not values:
            return

        attrs = {
            "type": "complex" if complex_ else "real",
            "size": f"{declared_size[0]} {declared_size[1]}",
        }
        if units:
            attrs["units"] = str(units)
        node = ET.SubElement(parent, tag, attrs)
        for row, col, raw in values:
            value = ET.SubElement(node, "value")
            if scalar:
                value.set("name", component_code)
            else:
                out_name = rows[row]
                in_name = cols[col]
                if ".INVSIGCOV" not in tag and ".RESIDCOV" not in tag:
                    value.set(
                        "name",
                        self._component_name(
                            component_code,
                            output=out_name,
                            input_=in_name,
                        ),
                    )
                value.set("output", out_name)
                value.set("input", in_name)
            value.text = self._format_numeric(raw, complex_=complex_)

    @staticmethod
    def _component_name(code: str, *, output: str, input_: str) -> str:
        def axis(name: str) -> str:
            lowered = str(name).strip().lower()
            for candidate in ("x", "y", "z"):
                if lowered.endswith(candidate):
                    return candidate
            return lowered

        upper = code.upper()
        if upper in {"T", "TI"} and output.lower().endswith("z"):
            return f"{upper}{axis(input_)}"
        return f"{upper}{axis(output)}{axis(input_)}"

    # ------------------------------------------------------------------
    # Period range / formatting
    # ------------------------------------------------------------------

    def _append_period_range(self, root: ET.Element, document: EMTF) -> None:
        if document.periods is not None and len(document.periods):
            periods = np.asarray(document.periods, dtype=float)
            ET.SubElement(
                root,
                "PeriodRange",
                {
                    "min": self._format_period_number(
                        float(np.min(periods))
                    ),
                    "max": self._format_period_number(
                        float(np.max(periods))
                    ),
                },
            )
            return
        stored = document.metadata.get("period_range")
        if isinstance(stored, dict) and stored:
            ET.SubElement(
                root,
                "PeriodRange",
                {str(key): str(value) for key, value in stored.items()},
            )

    def _format_data_number(self, value: float) -> str:
        return format(float(value), f".{self.precision}g")

    @staticmethod
    def _format_period_number(value: float) -> str:
        # Periods are positive coordinate values rather than noisy matrix
        # coefficients. Fifteen significant digits avoid binary display
        # artefacts while preserving practical float64 period precision.
        return format(float(value), ".15g")

    @staticmethod
    def _format_float(value: Any | None) -> str | None:
        if value is None:
            return None
        return format(float(value), ".15g")

    def _format_numeric(
        self,
        value: complex | float,
        *,
        complex_: bool,
    ) -> str:
        if complex_:
            number = complex(value)
            return (
                f"{self._format_data_number(number.real)} "
                f"{self._format_data_number(number.imag)}"
            )
        return self._format_data_number(float(np.real(value)))

    @staticmethod
    def _is_missing(value: Any) -> bool:
        try:
            if np.iscomplexobj(value):
                return bool(
                    np.isnan(np.real(value)) or np.isnan(np.imag(value))
                )
            return bool(np.isnan(value))
        except TypeError:
            return False

    @staticmethod
    def _is_finite(value: Any, *, complex_: bool) -> bool:
        if complex_:
            number = complex(value)
            return bool(np.isfinite(number.real) and np.isfinite(number.imag))
        try:
            return bool(np.isfinite(float(np.real(value))))
        except (TypeError, ValueError):
            return False

    def _problem(self, message: str) -> None:
        if self.strict:
            raise EMTFXMLSerializationError(message)
        _warn(message)
