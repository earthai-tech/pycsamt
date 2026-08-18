# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

"""Shared structural QC and inspection for airborne EM datasets.

These utilities deliberately assess representation quality, completeness, and
metadata consistency. They do not invent universal geophysical thresholds for
signal-to-noise ratio, frequency range, or anomaly quality. Technology-specific
scientific QC can later be registered on top of this common structural layer.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any

import numpy as np

from ..core.base import CoreObject

if TYPE_CHECKING:
    from ..emtf.document import EMTF
from .base import AirborneEMDataset, AirborneEMLine, AirborneEMRecord
from .registry import (
    get_airborne_technology,
    identify_airborne_technologies,
)

def _emtf_class():
    from ..emtf.document import EMTF

    return EMTF


__all__ = [
    "AirborneQCIssue",
    "AirborneInspection",
    "AirborneQCReport",
    "inspect_airborne",
    "assess_airborne_qc",
]


@dataclass(frozen=True)
class AirborneQCIssue:
    """One structural/metadata QC finding."""

    code: str
    severity: str
    message: str
    line_id: str | None = None
    sample_id: str | None = None

    def __post_init__(self) -> None:
        code = str(self.code).strip().lower()
        severity = str(self.severity).strip().lower()
        if not code:
            raise ValueError("QC issue code must be non-empty")
        if severity not in {"info", "warning", "error"}:
            raise ValueError("QC severity must be info, warning, or error")
        object.__setattr__(self, "code", code)
        object.__setattr__(self, "severity", severity)
        object.__setattr__(self, "message", str(self.message).strip())


@dataclass(repr=False)
class AirborneInspection(CoreObject):
    """Compact scientific inventory of an airborne object."""

    object_type: str
    technologies: tuple[str, ...] = field(default_factory=tuple)
    n_lines: int = 0
    n_samples: int = 0
    n_records: int = 0
    transfer_function_names: tuple[str, ...] = field(default_factory=tuple)
    bbox: Any | None = None
    attrs: dict[str, Any] = field(default_factory=dict)


@dataclass(repr=False)
class AirborneQCReport(CoreObject):
    """Common structural QC report for an airborne dataset."""

    technologies: tuple[str, ...]
    metrics: dict[str, Any]
    line_metrics: dict[str, dict[str, Any]]
    issues: tuple[AirborneQCIssue, ...] = field(default_factory=tuple)

    @property
    def status(self) -> str:
        """Return ``error``, ``warning``, or ``pass``."""
        severities = {issue.severity for issue in self.issues}
        if "error" in severities:
            return "error"
        if "warning" in severities:
            return "warning"
        return "pass"

    @property
    def errors(self) -> tuple[AirborneQCIssue, ...]:
        return tuple(
            issue for issue in self.issues if issue.severity == "error"
        )

    @property
    def warnings(self) -> tuple[AirborneQCIssue, ...]:
        return tuple(
            issue for issue in self.issues if issue.severity == "warning"
        )


def inspect_airborne(obj: Any) -> AirborneInspection:
    """Return a compact inventory for dataset, line, record, or EMTF object."""
    technologies = identify_airborne_technologies(obj)
    if isinstance(obj, AirborneEMDataset):
        return AirborneInspection(
            object_type="dataset",
            technologies=technologies,
            n_lines=obj.n_lines,
            n_samples=obj.n_samples,
            n_records=obj.n_records,
            transfer_function_names=obj.transfer_function_names,
            bbox=obj.bbox,
            attrs={"name": obj.name, "method": obj.method},
        )
    if isinstance(obj, AirborneEMLine):
        return AirborneInspection(
            object_type="line",
            technologies=technologies,
            n_lines=1,
            n_samples=obj.n_samples,
            n_records=obj.n_records,
            transfer_function_names=obj.transfer_function_names,
            bbox=obj.bbox,
            attrs={"line_id": obj.line_id},
        )
    if isinstance(obj, AirborneEMRecord):
        return AirborneInspection(
            object_type="record",
            technologies=technologies,
            n_records=1,
            transfer_function_names=obj.transfer_function_names,
            attrs={"sample_id": obj.sample_id},
        )
    if isinstance(obj, _emtf_class()):
        return AirborneInspection(
            object_type="emtf",
            technologies=technologies,
            n_records=1,
            transfer_function_names=tuple(obj.transfer_functions),
            attrs={"product_id": obj.product_id, "subtype": obj.subtype},
        )
    raise TypeError(
        "inspect_airborne expects AirborneEMDataset, AirborneEMLine, "
        "AirborneEMRecord, or EMTF"
    )


def _finite_complex_fraction(value: Any) -> tuple[int, int]:
    arr = np.asarray(value)
    if arr.size == 0:
        return 0, 0
    if np.iscomplexobj(arr):
        finite = np.isfinite(arr.real) & np.isfinite(arr.imag)
    else:
        finite = np.isfinite(arr)
    return int(np.count_nonzero(finite)), int(arr.size)


def _navigation_location_count(line: AirborneEMLine) -> int:
    nav = line.navigation
    mask = np.zeros(nav.n_samples, dtype=bool)
    if nav.latitude is not None and nav.longitude is not None:
        mask |= np.isfinite(nav.latitude) & np.isfinite(nav.longitude)
    if nav.easting is not None and nav.northing is not None:
        mask |= np.isfinite(nav.easting) & np.isfinite(nav.northing)
    return int(np.count_nonzero(mask))


def _clearance_count(line: AirborneEMLine) -> int:
    values = line.navigation.clearance_values
    if values is None:
        return 0
    return int(np.count_nonzero(np.isfinite(values)))


def _reference_present(doc: "EMTF", technology: str) -> bool:
    if technology == "mobilemt":
        meta = doc.attrs.get("mobilemt")
        return bool(
            isinstance(meta, dict)
            and isinstance(meta.get("reference_station"), dict)
            and meta["reference_station"].get("station_id")
        )
    if technology in {"ztem", "airmt"}:
        processing = doc.processing
        remote = None if processing is None else processing.remote_reference
        return bool(remote is not None and remote.site)
    return True


def _line_metrics(line: AirborneEMLine) -> dict[str, Any]:
    n_samples = line.n_samples
    n_records = line.n_records
    return {
        "n_samples": n_samples,
        "n_records": n_records,
        "record_coverage_fraction": (
            float(n_records / n_samples) if n_samples else 0.0
        ),
        "location_fraction": (
            float(_navigation_location_count(line) / n_samples)
            if n_samples else 0.0
        ),
        "clearance_fraction": (
            float(_clearance_count(line) / n_samples)
            if n_samples else 0.0
        ),
        "transfer_function_names": line.transfer_function_names,
    }


def assess_airborne_qc(dataset: AirborneEMDataset) -> AirborneQCReport:
    """Assess common structural completeness and metadata consistency.

    The report is intentionally descriptive. A missing record, missing
    covariance, or absent coordinates can be important without being a
    universal processing failure, so only internally inconsistent scientific
    states are classified as errors.
    """
    if not isinstance(dataset, AirborneEMDataset):
        raise TypeError("dataset must be an AirborneEMDataset")

    technologies = identify_airborne_technologies(dataset)
    issues: list[AirborneQCIssue] = []
    line_metrics = {
        line.line_id: _line_metrics(line)
        for line in dataset.iter_lines()
    }

    total_samples = dataset.n_samples
    total_records = dataset.n_records
    emtf_records = 0
    valid_frequency_records = 0
    finite_count = 0
    value_count = 0
    tf_count = 0
    primary_tf_count = 0
    derived_tf_count = 0
    variance_count = 0
    invsig_count = 0
    resid_count = 0
    required_reference_records = 0
    present_reference_records = 0

    for line in dataset.iter_lines():
        if line.n_records < line.n_samples:
            issues.append(
                AirborneQCIssue(
                    code="missing_em_records",
                    severity="info",
                    message=(
                        f"{line.n_samples - line.n_records} navigation "
                        "samples have no attached EM record"
                    ),
                    line_id=line.line_id,
                )
            )
        if _navigation_location_count(line) == 0:
            issues.append(
                AirborneQCIssue(
                    code="missing_navigation_coordinates",
                    severity="warning",
                    message="line has no finite geographic/projected position",
                    line_id=line.line_id,
                )
            )

        explicit_line_tech = identify_airborne_technologies(line)
        if len(explicit_line_tech) > 1:
            issues.append(
                AirborneQCIssue(
                    code="mixed_line_technology",
                    severity="warning",
                    message=(
                        "line contains multiple technologies: "
                        + ", ".join(explicit_line_tech)
                    ),
                    line_id=line.line_id,
                )
            )

        for record in line.iter_records():
            doc = record.emtf
            if doc is None:
                issues.append(
                    AirborneQCIssue(
                        code="record_without_emtf",
                        severity="info",
                        message="record has no EMTF transfer-function payload",
                        line_id=line.line_id,
                        sample_id=record.sample_id,
                    )
                )
                continue
            emtf_records += 1
            freq = doc.frequency
            if (
                freq is not None
                and np.asarray(freq).size > 0
                and np.all(np.isfinite(freq))
                and np.all(np.asarray(freq) > 0.0)
            ):
                valid_frequency_records += 1
            else:
                issues.append(
                    AirborneQCIssue(
                        code="invalid_frequency_axis",
                        severity="error",
                        message=(
                            "EMTF record has no valid positive frequency axis"
                        ),
                        line_id=line.line_id,
                        sample_id=record.sample_id,
                    )
                )

            record_technologies = identify_airborne_technologies(record)
            for technology in record_technologies:
                definition = get_airborne_technology(technology)
                if definition is None or not definition.reference_required:
                    continue
                required_reference_records += 1
                if _reference_present(doc, technology):
                    present_reference_records += 1
                else:
                    issues.append(
                        AirborneQCIssue(
                            code="missing_reference_station",
                            severity="warning",
                            message=(
                                f"{definition.label} record lacks explicit "
                                "reference-station metadata"
                            ),
                            line_id=line.line_id,
                            sample_id=record.sample_id,
                        )
                    )

            for tf in doc.transfer_functions.values():
                tf_count += 1
                finite, total = _finite_complex_fraction(tf.data)
                finite_count += finite
                value_count += total
                definition = tf.definition
                is_primary = (
                    definition is None or definition.intention == "primary"
                )
                if not is_primary:
                    derived_tf_count += 1
                    continue
                primary_tf_count += 1
                if tf.get_estimate("VAR") is not None:
                    variance_count += 1
                if tf.get_estimate("INVSIGCOV") is not None:
                    invsig_count += 1
                if tf.get_estimate("RESIDCOV") is not None:
                    resid_count += 1

    metrics: dict[str, Any] = {
        "n_lines": dataset.n_lines,
        "n_samples": total_samples,
        "n_records": total_records,
        "n_emtf_records": emtf_records,
        "record_coverage_fraction": (
            float(total_records / total_samples) if total_samples else 0.0
        ),
        "emtf_coverage_fraction": (
            float(emtf_records / total_samples) if total_samples else 0.0
        ),
        "valid_frequency_record_fraction": (
            float(valid_frequency_records / emtf_records)
            if emtf_records else 0.0
        ),
        "finite_response_fraction": (
            float(finite_count / value_count) if value_count else 0.0
        ),
        "n_transfer_functions": tf_count,
        "n_primary_transfer_functions": primary_tf_count,
        "n_derived_transfer_functions": derived_tf_count,
        "variance_tf_fraction": (
            float(variance_count / primary_tf_count)
            if primary_tf_count else 0.0
        ),
        "inverse_signal_covariance_tf_fraction": (
            float(invsig_count / primary_tf_count)
            if primary_tf_count else 0.0
        ),
        "residual_covariance_tf_fraction": (
            float(resid_count / primary_tf_count)
            if primary_tf_count else 0.0
        ),
        "reference_metadata_fraction": (
            float(present_reference_records / required_reference_records)
            if required_reference_records else 1.0
        ),
    }

    if len(technologies) > 1:
        issues.append(
            AirborneQCIssue(
                code="mixed_dataset_technology",
                severity="info",
                message=(
                    "dataset contains multiple airborne technologies: "
                    + ", ".join(technologies)
                ),
            )
        )

    return AirborneQCReport(
        technologies=technologies,
        metrics=metrics,
        line_metrics=line_metrics,
        issues=tuple(issues),
    )
