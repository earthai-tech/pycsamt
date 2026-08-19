# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

"""EMTF XML -> format-neutral :class:`pycsamt.emtf.EMTF` reader."""

from __future__ import annotations

from dataclasses import asdict
from os import PathLike
from typing import Any, IO
import xml.etree.ElementTree as ET

import numpy as np

from ...api.property import PyCSAMTObject
from ...metadata import (
    ChannelMeta,
    CopyrightInfo,
    LocationMeta,
    OrientationMeta,
    Person,
    ProcessingMeta,
    ProvenanceMeta,
    QualityComment,
    Reference,
    RemoteReferenceMeta,
    SiteLayout,
    SiteMeta,
    Software,
    TransferFunctionQuality,
)
from ..datatypes import get_emtf_datatype
from ..document import EMTF
from ..estimates import StatisticalEstimate
from ..transfer import TransferFunction
from .constants import (
    EMTF_ROOT,
    ESTIMATE_CODES,
    FREQUENCY_UNIT_ALIASES,
    PERIOD_UNIT_ALIASES,
)
from .parser import (
    EMTFXMLParseError,
    XMLDataTypeSpec,
    child,
    children,
    element_to_mapping,
    load_xml_root,
    local_name,
    parse_numeric_text,
    text,
    warn,
)

__all__ = ["EMTFXMLReader", "read_emtf_xml"]


class EMTFXMLReader(PyCSAMTObject):
    """Read EMTF XML without routing the document through EDI.

    Parameters
    ----------
    strict : bool, default=True
        In strict mode malformed scientific content raises
        :class:`EMTFXMLParseError`. In permissive mode recoverable problems are
        warned about and the reader preserves as much content as possible.
    """

    def __init__(self, *, strict: bool = True) -> None:
        self.strict = bool(strict)
        self.source_label = "<unknown>"

    def read(
        self,
        source: str | bytes | PathLike[str] | IO[str] | IO[bytes],
    ) -> EMTF:
        """Read one EMTF XML source into the scientific EMTF model."""
        root, self.source_label = load_xml_root(source)
        if local_name(root.tag) != EMTF_ROOT:
            self._problem(
                f"root element must be <{EMTF_ROOT}>, got "
                f"<{local_name(root.tag)}>"
            )

        tags = self._parse_tags(text(root, "Tags"))
        provenance = self._parse_provenance(child(root, "Provenance"))
        copyright_info = self._parse_copyright(child(root, "Copyright"))
        site_node = child(root, "Site")
        site = self._parse_site(site_node)
        orientation = self._parse_orientation(site_node, root)
        quality = self._parse_quality(site_node)
        processing = self._parse_processing(child(root, "ProcessingInfo"))
        layout = self._parse_site_layout(root)
        field_notes = self._parse_field_notes(root)

        metadata = self._parse_auxiliary_metadata(root)
        specs = self._parse_data_type_specs(root)
        metadata["xml_data_types"] = [asdict(spec) for spec in specs.values()]
        metadata["xml_statistical_estimates"] = (
            self._parse_estimate_declarations(root)
        )
        metadata["source_format"] = "emtf_xml"
        metadata["source"] = self.source_label

        period_nodes = self._period_nodes(root)
        periods, period_nodes = self._parse_periods(root, period_nodes)

        document = EMTF(
            product_id=text(root, "ProductId"),
            description=text(root, "Description"),
            subtype=text(root, "SubType"),
            tags=tags,
            periods=periods,
            provenance=provenance,
            copyright=copyright_info,
            site=site,
            site_layout=layout,
            orientation=orientation,
            processing=processing,
            quality=quality,
            field_notes=field_notes,
            metadata=metadata,
        )

        if period_nodes:
            if not specs:
                specs = self._infer_data_type_specs(period_nodes)
                document.metadata["xml_data_types_inferred"] = True
                document.metadata["xml_data_types"] = [
                    asdict(spec) for spec in specs.values()
                ]
            for spec in specs.values():
                tf = self._parse_transfer_function(
                    spec,
                    period_nodes=period_nodes,
                    periods=periods,
                    layout=layout,
                )
                if tf is not None:
                    document.add_transfer_function(tf, replace=True)

        return document

    # ------------------------------------------------------------------
    # Document/header metadata
    # ------------------------------------------------------------------

    @staticmethod
    def _parse_tags(value: str | None) -> tuple[str, ...]:
        if not value:
            return ()
        raw = value.replace(";", ",")
        return tuple(
            item.strip().lower() for item in raw.split(",") if item.strip()
        )

    def _parse_provenance(
        self, node: ET.Element | None
    ) -> ProvenanceMeta | None:
        if node is None:
            return None
        return ProvenanceMeta(
            create_time=text(node, "CreateTime"),
            creating_application=text(node, "CreatingApplication"),
            creator=self._parse_person(child(node, "Creator")),
            submitter=self._parse_person(child(node, "Submitter")),
        )

    @staticmethod
    def _parse_person(node: ET.Element | None) -> Person | None:
        if node is None:
            return None
        values = {
            "name": text(node, "Name"),
            "email": text(node, "Email"),
            "organization": text(node, "Org"),
            "organization_url": text(node, "OrgUrl"),
        }
        if not any(values.values()):
            return None
        return Person(**values)

    def _parse_copyright(
        self, node: ET.Element | None
    ) -> CopyrightInfo | None:
        if node is None:
            return None
        citation = child(node, "Citation")
        title = text(citation, "Title") or ""
        authors = text(citation, "Authors") or ""
        year_text = text(citation, "Year")
        year: int | None = None
        extra: dict[str, Any] = {}
        if year_text:
            try:
                year = int(year_text)
            except ValueError:
                extra["year_text"] = year_text

        doi_raw = text(citation, "DOI")
        doi = None
        if doi_raw:
            candidate = doi_raw.strip()
            if candidate.lower().startswith("doi:"):
                candidate = candidate[4:].strip()
            if Reference.DOI_PATTERN.match(candidate):
                doi = candidate
            else:
                extra["doi"] = doi_raw

        survey_doi = text(citation, "SurveyDOI")
        if survey_doi:
            extra["survey_doi"] = survey_doi
        for key in (
            "SelectedPublications",
            "Acknowledgement",
            "AdditionalInfo",
        ):
            value = text(node, key)
            if value:
                extra[key.lower()] = value

        reference = Reference(
            author=authors,
            title=title,
            year=year,
            doi=doi,
            extra=extra,
        )
        return CopyrightInfo(
            release_status=text(node, "ReleaseStatus") or "",
            conditions_of_use=text(node, "ConditionsOfUse") or "",
            reference=reference,
        )

    def _parse_site(self, node: ET.Element | None) -> SiteMeta | None:
        if node is None:
            return None
        location = self._parse_location(child(node, "Location"))
        year = self._integer(text(node, "YearCollected"), "YearCollected")
        extra: dict[str, Any] = {}
        run_list = text(node, "RunList")
        if run_list:
            extra["run_list"] = run_list
        comments = [
            (item.text or "").strip()
            for item in children(node, "Comments")
            if (item.text or "").strip()
        ]
        if comments:
            extra["comments"] = comments
        try:
            return SiteMeta(
                project=text(node, "Project"),
                survey=text(node, "Survey"),
                year_collected=year,
                country=text(node, "Country"),
                site_id=text(node, "Id"),
                name=text(node, "Name"),
                location=location,
                acquired_by=text(node, "AcquiredBy"),
                start=text(node, "Start"),
                end=text(node, "End"),
                extra=extra,
            )
        except (TypeError, ValueError) as exc:
            self._problem(f"invalid <Site> metadata: {exc}")
            return None

    def _parse_location(
        self, node: ET.Element | None
    ) -> LocationMeta | None:
        if node is None:
            return None
        elevation = child(node, "Elevation")
        declination = child(node, "Declination")
        try:
            return LocationMeta(
                latitude=self._float(text(node, "Latitude"), "Latitude"),
                longitude=self._float(text(node, "Longitude"), "Longitude"),
                elevation=self._float(
                    text(elevation), "Elevation", required=False
                ),
                datum=node.attrib.get("datum") or "WGS84",
                elevation_units=(
                    elevation.attrib.get("units", "meters")
                    if elevation is not None
                    else "meters"
                ),
                declination=self._float(
                    text(declination), "Declination", required=False
                ),
                declination_epoch=self._float(
                    declination.attrib.get("epoch")
                    if declination is not None
                    else None,
                    "Declination epoch",
                    required=False,
                ),
            )
        except (TypeError, ValueError) as exc:
            self._problem(f"invalid <Location> metadata: {exc}")
            return None

    def _parse_orientation(
        self,
        site: ET.Element | None,
        root: ET.Element,
    ) -> OrientationMeta | None:
        node = child(site, "Orientation")
        rotation_info = None
        copyright_node = child(root, "Copyright")
        for parent in (copyright_node, root):
            candidate = text(parent, "RotationInfo")
            if candidate:
                rotation_info = candidate
                break
        if node is None and rotation_info is None:
            return None
        mode = text(node)
        angle = self._float(
            node.attrib.get("angle_to_geographic_north")
            if node is not None
            else None,
            "Orientation angle_to_geographic_north",
            required=False,
        )
        try:
            return OrientationMeta(
                mode=mode,
                angle_to_geographic_north=angle,
                rotation_info=rotation_info,
            )
        except (TypeError, ValueError) as exc:
            self._problem(f"invalid orientation metadata: {exc}")
            return None

    def _parse_quality(
        self, site: ET.Element | None
    ) -> TransferFunctionQuality | None:
        if site is None:
            return None
        notes = child(site, "DataQualityNotes")
        warnings_node = child(site, "DataQualityWarnings")
        if notes is None and warnings_node is None:
            return None

        comments = self._quality_comments(notes)
        warning_comments = self._quality_comments(warnings_node)
        rating = self._integer(text(notes, "Rating"), "quality Rating")
        warning_flag = self._integer(
            text(warnings_node, "Flag"), "quality Flag"
        )
        try:
            return TransferFunctionQuality(
                rating=rating,
                good_from_period=self._float(
                    text(notes, "GoodFromPeriod"),
                    "GoodFromPeriod",
                    required=False,
                ),
                good_to_period=self._float(
                    text(notes, "GoodToPeriod"),
                    "GoodToPeriod",
                    required=False,
                ),
                comments=comments,
                warning_flag=warning_flag,
                warnings=warning_comments,
            )
        except (TypeError, ValueError) as exc:
            self._problem(f"invalid data-quality metadata: {exc}")
            return None

    @staticmethod
    def _quality_comments(
        node: ET.Element | None,
    ) -> list[QualityComment]:
        if node is None:
            return []
        result: list[QualityComment] = []
        for item in children(node, "Comments"):
            value = (item.text or "").strip()
            if not value:
                continue
            result.append(
                QualityComment(text=value, author=item.attrib.get("author"))
            )
        return result

    def _parse_processing(
        self, node: ET.Element | None
    ) -> ProcessingMeta | None:
        if node is None:
            return None
        software_node = child(node, "ProcessingSoftware")
        software = None
        if software_node is not None and any(
            text(software_node, name)
            for name in ("Name", "LastMod", "Author")
        ):
            software = Software(
                name=text(software_node, "Name") or "unknown",
                release=text(software_node, "LastMod"),
                author=Person(name=text(software_node, "Author")),
            )

        remote_node = child(node, "RemoteRef")
        remote_site = None
        remote_info = child(node, "RemoteInfo")
        if remote_info is not None:
            remote_site_node = child(remote_info, "Site")
            remote_site = text(remote_site_node, "Id")
        remote = None
        if remote_node is not None or remote_site is not None:
            remote = RemoteReferenceMeta(
                reference_type=(
                    remote_node.attrib.get("type")
                    if remote_node is not None
                    else None
                ),
                site=remote_site,
            )

        extra: dict[str, Any] = {}
        process_date = text(node, "ProcessDate")
        if process_date:
            extra["process_date"] = process_date
        if remote_info is not None:
            extra["remote_info"] = element_to_mapping(remote_info)

        return ProcessingMeta(
            sign_convention=text(node, "SignConvention"),
            processed_by=text(node, "ProcessedBy"),
            software=software,
            remote_reference=remote,
            processing_tag=text(node, "ProcessingTag"),
            extra=extra,
        )

    def _parse_site_layout(self, root: ET.Element) -> SiteLayout | None:
        layout_node = child(root, "SiteLayout")
        legacy = False
        if layout_node is None:
            if child(root, "InputChannels") is None and child(
                root, "OutputChannels"
            ) is None:
                return None
            layout_node = root
            legacy = True
            warn(
                f"{self.source_label}: reading legacy EMTF XML without "
                "<SiteLayout>; data orientation must be interpreted with "
                "channel geometry"
            )

        input_node = child(layout_node, "InputChannels")
        output_node = child(layout_node, "OutputChannels")
        input_channels = self._parse_channel_group(input_node)
        output_channels = self._parse_channel_group(output_node)
        try:
            return SiteLayout(
                input_channels=input_channels,
                output_channels=output_channels,
                input_units=(
                    input_node.attrib.get("units")
                    if input_node is not None
                    else None
                ),
                output_units=(
                    output_node.attrib.get("units")
                    if output_node is not None
                    else None
                ),
                input_reference=(
                    input_node.attrib.get("ref")
                    if input_node is not None
                    else None
                ),
                output_reference=(
                    output_node.attrib.get("ref")
                    if output_node is not None
                    else None
                ),
                extra={"legacy_without_sitelayout": True} if legacy else {},
            )
        except (TypeError, ValueError) as exc:
            self._problem(f"invalid <SiteLayout>: {exc}")
            return None

    def _parse_channel_group(
        self, node: ET.Element | None
    ) -> list[ChannelMeta]:
        if node is None:
            return []
        result: list[ChannelMeta] = []
        for item in list(node):
            kind = local_name(item.tag)
            if kind not in {"Magnetic", "Electric"}:
                continue
            known = {
                "name",
                "orientation",
                "tilt",
                "x",
                "y",
                "z",
                "x2",
                "y2",
                "z2",
                "units",
                "ref",
                "id",
            }
            extra = {k: v for k, v in item.attrib.items() if k not in known}
            name = item.attrib.get("name")
            if not name:
                self._problem(f"{kind} channel is missing required name")
                continue
            try:
                result.append(
                    ChannelMeta(
                        name=name,
                        field_type=kind,
                        orientation=self._float(
                            item.attrib.get("orientation"),
                            f"channel {name} orientation",
                            required=False,
                        ),
                        tilt=self._float(
                            item.attrib.get("tilt"),
                            f"channel {name} tilt",
                            required=False,
                        ),
                        x=self._float(
                            item.attrib.get("x"),
                            f"channel {name} x",
                            required=False,
                        ),
                        y=self._float(
                            item.attrib.get("y"),
                            f"channel {name} y",
                            required=False,
                        ),
                        z=self._float(
                            item.attrib.get("z"),
                            f"channel {name} z",
                            required=False,
                        ),
                        x2=self._float(
                            item.attrib.get("x2"),
                            f"channel {name} x2",
                            required=False,
                        ),
                        y2=self._float(
                            item.attrib.get("y2"),
                            f"channel {name} y2",
                            required=False,
                        ),
                        z2=self._float(
                            item.attrib.get("z2"),
                            f"channel {name} z2",
                            required=False,
                        ),
                        units=item.attrib.get("units"),
                        reference=item.attrib.get("ref"),
                        sensor_id=item.attrib.get("id"),
                        extra=extra,
                    )
                )
            except (TypeError, ValueError) as exc:
                self._problem(f"invalid channel {name!r}: {exc}")
        return result

    def _parse_field_notes(self, root: ET.Element) -> dict[str, Any]:
        result: dict[str, Any] = {}
        for index, node in enumerate(children(root, "FieldNotes"), start=1):
            key = node.attrib.get("run") or f"run_{index}"
            value = element_to_mapping(node)
            if key in result:
                current = result[key]
                if not isinstance(current, list):
                    current = [current]
                current.append(value)
                result[key] = current
            else:
                result[key] = value
        return result

    @staticmethod
    def _parse_auxiliary_metadata(root: ET.Element) -> dict[str, Any]:
        result: dict[str, Any] = {}
        notes = child(root, "Notes")
        if notes is not None:
            result["notes"] = element_to_mapping(notes)
        for key in ("ExternalUrl", "PrimaryData", "Attachment", "GridOrigin"):
            values = [element_to_mapping(item) for item in children(root, key)]
            if values:
                result[key.lower()] = values if len(values) > 1 else values[0]
        period_range = child(root, "PeriodRange")
        if period_range is not None:
            result["period_range"] = dict(period_range.attrib)
        return result

    # ------------------------------------------------------------------
    # Data-type declarations and periods
    # ------------------------------------------------------------------

    def _parse_data_type_specs(
        self, root: ET.Element
    ) -> dict[str, XMLDataTypeSpec]:
        parent = child(root, "DataTypes")
        if parent is None:
            return {}
        result: dict[str, XMLDataTypeSpec] = {}
        for node in children(parent, "DataType"):
            raw_name = (node.attrib.get("name") or "").strip().upper()
            raw_tag = (text(node, "Tag") or "").strip().lower()
            registered = get_emtf_datatype(raw_tag or raw_name)
            if not raw_name and registered is not None:
                raw_name = registered.name
            if not raw_tag and registered is not None:
                raw_tag = registered.tag
            if not raw_name or not raw_tag:
                self._problem("<DataType> requires a name and semantic Tag")
                continue

            raw_kind = (node.attrib.get("type") or "").strip().lower()
            if raw_kind not in {"real", "complex"}:
                raw_kind = registered.data_kind if registered else "complex"
                self._problem(
                    f"DataType {raw_name} has invalid/missing type; "
                    f"using {raw_kind!r}",
                    recoverable=True,
                )

            intention_text = (text(node, "Intention") or "").lower()
            intention = (
                "derived" if "derived" in intention_text else "primary"
            )
            see_also = tuple(
                item.strip().lower()
                for item in (text(node, "SeeAlso") or "").split(",")
                if item.strip()
            )
            spec = XMLDataTypeSpec(
                name=raw_name,
                tag=raw_tag,
                data_kind=raw_kind,
                input_kind=(node.attrib.get("input") or None),
                output_kind=(node.attrib.get("output") or None),
                units=node.attrib.get("units"),
                intention=intention,
                description=text(node, "Description") or "",
                derived_from=text(node, "DerivedFrom"),
                see_also=see_also,
                external_url=text(node, "ExternalUrl"),
            )
            if raw_name in result:
                self._problem(f"duplicate DataType declaration: {raw_name}")
                continue
            result[raw_name] = spec
        return result

    @staticmethod
    def _parse_estimate_declarations(root: ET.Element) -> list[dict[str, Any]]:
        parent = child(root, "StatisticalEstimates")
        if parent is None:
            return []
        out: list[dict[str, Any]] = []
        for node in children(parent, "Estimate"):
            out.append(
                {
                    "name": node.attrib.get("name"),
                    "type": node.attrib.get("type"),
                    "description": text(node, "Description"),
                    "intention": text(node, "Intention"),
                    "tag": text(node, "Tag"),
                    "external_url": text(node, "ExternalUrl"),
                }
            )
        return out

    @staticmethod
    def _period_nodes(root: ET.Element) -> list[ET.Element]:
        data = child(root, "Data")
        if data is None:
            return []
        return children(data, "Period")

    def _parse_periods(
        self,
        root: ET.Element,
        period_nodes: list[ET.Element],
    ) -> tuple[np.ndarray | None, list[ET.Element]]:
        """Return ``(periods, kept_period_nodes)``.

        In strict mode an invalid period raises immediately, as before. In
        permissive mode an invalid period is dropped -- together with its
        ``<Period>`` node -- rather than substituted with ``NaN``, since a
        ``NaN`` period cannot survive :class:`~pycsamt.emtf.document.EMTF`'s
        own validation and would otherwise defeat permissive mode's purpose
        by crashing anyway with an unrelated, uncaught ``ValueError``.
        """
        if not period_nodes:
            return None, []
        data = child(root, "Data")
        declared_count = self._integer(
            data.attrib.get("count") if data is not None else None,
            "Data count",
        )
        if declared_count is not None and declared_count != len(period_nodes):
            self._problem(
                "<Data count> does not match the number of <Period> "
                f"elements: {declared_count} != {len(period_nodes)}"
            )

        values: list[float] = []
        kept_nodes: list[ET.Element] = []
        for index, node in enumerate(period_nodes):
            raw = node.attrib.get("value")
            numeric = self._float(raw, f"Period[{index}] value", required=True)
            if numeric is None or numeric <= 0.0:
                self._problem(f"Period[{index}] must be positive")
                continue
            units = (node.attrib.get("units") or "secs").strip().lower()
            if units in FREQUENCY_UNIT_ALIASES:
                numeric = 1.0 / numeric
            elif units not in PERIOD_UNIT_ALIASES:
                self._problem(
                    f"Period[{index}] has unsupported units {units!r}; "
                    "assuming seconds",
                    recoverable=True,
                )
            values.append(float(numeric))
            kept_nodes.append(node)
        if not values:
            self._problem("no usable <Period> elements remain")
            return None, []
        return np.asarray(values, dtype=float), kept_nodes

    def _infer_data_type_specs(
        self, period_nodes: list[ET.Element]
    ) -> dict[str, XMLDataTypeSpec]:
        names: list[str] = []
        for period in period_nodes:
            for node in list(period):
                name = local_name(node.tag)
                if "." in name:
                    continue
                code = name.upper()
                if code not in names:
                    names.append(code)
        result: dict[str, XMLDataTypeSpec] = {}
        for code in names:
            definition = get_emtf_datatype(code)
            if definition is None:
                self._problem(
                    f"cannot infer undeclared data type {code!r}; skipping",
                    recoverable=True,
                )
                continue
            result[code] = XMLDataTypeSpec(
                name=code,
                tag=definition.tag,
                data_kind=definition.data_kind,
                input_kind=definition.input_kind,
                output_kind=definition.output_kind,
                units=definition.units,
                intention=definition.intention,
                description=definition.description,
                derived_from=definition.derived_from,
                see_also=definition.see_also,
            )
        return result

    # ------------------------------------------------------------------
    # Matrix data and statistical estimates
    # ------------------------------------------------------------------

    def _parse_transfer_function(
        self,
        spec: XMLDataTypeSpec,
        *,
        period_nodes: list[ET.Element],
        periods: np.ndarray | None,
        layout: SiteLayout | None,
    ) -> TransferFunction | None:
        input_names, output_names = self._channel_names_for_spec(spec, layout)
        if not spec.is_scalar and (not input_names or not output_names):
            inferred_in, inferred_out = self._infer_component_channels(
                spec, period_nodes
            )
            input_names = input_names or inferred_in
            output_names = output_names or inferred_out

        if not spec.is_scalar and (not input_names or not output_names):
            self._problem(
                f"cannot determine input/output channels for {spec.name}"
            )
            return None

        nper = len(period_nodes)
        nout = max(1, len(output_names))
        nin = max(1, len(input_names))
        dtype = complex if spec.data_kind == "complex" else float
        fill = np.nan + 1j * np.nan if dtype is complex else np.nan
        data = np.full((nper, nout, nin), fill, dtype=dtype)
        seen_data = False
        for iper, period_node in enumerate(period_nodes):
            data_node = self._direct_named_child(period_node, spec.name)
            if data_node is None:
                continue
            seen_data = True
            self._fill_matrix_values(
                data[iper],
                data_node,
                input_names=input_names,
                output_names=output_names,
                complex_=spec.data_kind == "complex",
                label=f"{spec.name} period {iper}",
                scalar=spec.is_scalar,
            )

        if not seen_data:
            return None

        tf = TransferFunction(
            name=spec.tag,
            data=data,
            input_channels=input_names,
            output_channels=output_names,
            units=spec.units,
            periods=periods,
            attrs={
                "xml_name": spec.name,
                "xml_data_kind": spec.data_kind,
                "xml_input_kind": spec.input_kind,
                "xml_output_kind": spec.output_kind,
                "xml_intention": spec.intention,
            },
        )

        for code, kind in ESTIMATE_CODES.items():
            estimate = self._parse_estimate(
                spec,
                code=code,
                kind=kind,
                period_nodes=period_nodes,
                input_names=input_names,
                output_names=output_names,
            )
            if estimate is not None:
                tf.add_estimate(estimate)
        return tf

    @staticmethod
    def _direct_named_child(
        node: ET.Element,
        name: str,
    ) -> ET.Element | None:
        wanted = name.upper()
        for item in list(node):
            if local_name(item.tag).upper() == wanted:
                return item
        return None

    def _channel_names_for_spec(
        self,
        spec: XMLDataTypeSpec,
        layout: SiteLayout | None,
    ) -> tuple[tuple[str, ...], tuple[str, ...]]:
        if spec.is_scalar:
            return (), ()
        if layout is None:
            return (), ()

        def family(
            channels: list[ChannelMeta],
            kind: str | None,
        ) -> tuple[str, ...]:
            if not kind:
                return ()
            normalized = kind.strip().upper()
            if normalized == "H":
                return tuple(ch.name for ch in channels if ch.is_magnetic)
            if normalized == "E":
                return tuple(ch.name for ch in channels if ch.is_electric)
            return tuple(ch.name for ch in channels)

        return (
            family(layout.input_channels, spec.input_kind),
            family(layout.output_channels, spec.output_kind),
        )

    def _infer_component_channels(
        self,
        spec: XMLDataTypeSpec,
        period_nodes: list[ET.Element],
    ) -> tuple[tuple[str, ...], tuple[str, ...]]:
        inputs: list[str] = []
        outputs: list[str] = []
        for period in period_nodes:
            node = self._direct_named_child(period, spec.name)
            if node is None:
                continue
            for value in children(node, "value"):
                raw_input = value.attrib.get("input")
                raw_output = value.attrib.get("output")
                if raw_input and raw_input not in inputs:
                    inputs.append(raw_input)
                if raw_output and raw_output not in outputs:
                    outputs.append(raw_output)
        return tuple(inputs), tuple(outputs)

    def _fill_matrix_values(
        self,
        target: np.ndarray,
        node: ET.Element,
        *,
        input_names: tuple[str, ...],
        output_names: tuple[str, ...],
        complex_: bool,
        label: str,
        scalar: bool,
    ) -> None:
        if scalar:
            values = children(node, "value")
            if not values:
                return
            if len(values) > 1:
                self._problem(f"{label} contains multiple scalar values")
            try:
                target[0, 0] = parse_numeric_text(
                    text(values[-1]), complex_=complex_
                )
            except EMTFXMLParseError as exc:
                self._problem(f"{label}: {exc}")
            return

        input_map = {name.lower(): i for i, name in enumerate(input_names)}
        output_map = {name.lower(): i for i, name in enumerate(output_names)}
        assigned: set[tuple[int, int]] = set()
        for value in children(node, "value"):
            in_name = (value.attrib.get("input") or "").strip()
            out_name = (value.attrib.get("output") or "").strip()
            if not in_name or not out_name:
                self._problem(
                    f"{label} value requires input and output "
                    "channel attributes"
                )
                continue
            i = input_map.get(in_name.lower())
            o = output_map.get(out_name.lower())
            if i is None or o is None:
                self._problem(
                    f"{label} references unknown channel "
                    f"output={out_name!r}, input={in_name!r}"
                )
                continue
            key = (o, i)
            if key in assigned:
                self._problem(
                    f"{label} has duplicate component "
                    f"output={out_name!r}, input={in_name!r}"
                )
            try:
                target[o, i] = parse_numeric_text(
                    text(value), complex_=complex_
                )
                assigned.add(key)
            except EMTFXMLParseError as exc:
                self._problem(f"{label}: {exc}")

    def _parse_estimate(
        self,
        spec: XMLDataTypeSpec,
        *,
        code: str,
        kind: str,
        period_nodes: list[ET.Element],
        input_names: tuple[str, ...],
        output_names: tuple[str, ...],
    ) -> StatisticalEstimate | None:
        if spec.is_scalar and code in {"INVSIGCOV", "RESIDCOV"}:
            return None

        if code == "VAR":
            rows = output_names
            cols = input_names
            complex_ = False
        elif code == "INVSIGCOV":
            rows = input_names
            cols = input_names
            complex_ = True
        else:
            rows = output_names
            cols = output_names
            complex_ = True

        if spec.is_scalar:
            rows = cols = ()
        nrow = max(1, len(rows))
        ncol = max(1, len(cols))
        dtype = complex if complex_ else float
        fill = np.nan + 1j * np.nan if complex_ else np.nan
        data = np.full((len(period_nodes), nrow, ncol), fill, dtype=dtype)
        seen = False
        for iper, period in enumerate(period_nodes):
            node = self._direct_named_child(period, f"{spec.name}.{code}")
            if node is None:
                continue
            seen = True
            self._fill_matrix_values(
                data[iper],
                node,
                input_names=tuple(cols),
                output_names=tuple(rows),
                complex_=complex_,
                label=f"{spec.name}.{code} period {iper}",
                scalar=spec.is_scalar,
            )
        if not seen:
            return None
        return StatisticalEstimate(
            name=code,
            kind=kind,
            data=data,
            attrs={"xml_parent": spec.name},
        )

    # ------------------------------------------------------------------
    # Conversion/diagnostic helpers
    # ------------------------------------------------------------------

    def _float(
        self,
        value: str | None,
        label: str,
        *,
        required: bool = False,
    ) -> float | None:
        if value is None or not str(value).strip():
            if required:
                self._problem(f"missing required numeric field {label}")
            return None
        try:
            return float(str(value).replace("D", "E").replace("d", "e"))
        except ValueError:
            self._problem(f"invalid numeric field {label}: {value!r}")
            return None

    def _integer(self, value: str | None, label: str) -> int | None:
        if value is None or not str(value).strip():
            return None
        try:
            return int(str(value).strip())
        except ValueError:
            self._problem(f"invalid integer field {label}: {value!r}")
            return None

    def _problem(
        self,
        message: str,
        *,
        recoverable: bool = False,
    ) -> None:
        full = f"{self.source_label}: {message}"
        if self.strict and not recoverable:
            raise EMTFXMLParseError(full)
        warn(full)


def read_emtf_xml(
    source: str | bytes | PathLike[str] | IO[str] | IO[bytes],
    *,
    strict: bool = True,
) -> EMTF:
    """Convenience wrapper around :class:`EMTFXMLReader`."""
    return EMTFXMLReader(strict=strict).read(source)
