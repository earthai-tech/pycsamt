# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

"""Array-to-scientific-object adapter contract for MobileMT products.

This module deliberately does not define a native vendor file reader. It
accepts already decoded arrays and maps them into the common pyCSAMT airborne
and EMTF model. A genuine delivery sample can later define the parser without
changing these scientific contracts.
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

import numpy as np

from ...emtf.document import EMTF
from ...emtf.estimates import StatisticalEstimate
from ...emtf.transfer import TransferFunction
from ...metadata import SurveyMeta
from ..base import AirborneEMDataset, AirborneEMLine, AirborneEMRecord
from ..navigation import NavigationTrack
from .base import MobileMTReferenceStation, MobileMTSystemSpec
from .constants import (
    MOBILEMT_ADMITTANCE_TAG,
    MOBILEMT_APPARENT_CONDUCTIVITY_FIELD,
    MOBILEMT_APPARENT_CONDUCTIVITY_UNITS,
    MOBILEMT_INPUT_CHANNELS,
    MOBILEMT_OUTPUT_CHANNELS,
)
from .datatypes import register_mobilemt_datatypes

__all__ = [
    "MobileMTValidationError",
    "validate_mobilemt_transfer_function",
    "build_mobilemt_emtf",
    "build_mobilemt_record",
    "build_mobilemt_line",
    "build_mobilemt_dataset",
]


class MobileMTValidationError(ValueError):
    """Raised when decoded MobileMT scientific arrays are inconsistent."""


def _as_frequency(value: Any, *, name: str = "frequency") -> np.ndarray:
    arr = np.asarray(value, dtype=float)
    if arr.ndim == 0:
        arr = arr.reshape(1)
    if arr.ndim != 1:
        raise MobileMTValidationError(f"{name} must be a 1-D array")
    if arr.size == 0:
        raise MobileMTValidationError(f"{name} must be non-empty")
    if not np.all(np.isfinite(arr)) or np.any(arr <= 0.0):
        raise MobileMTValidationError(
            f"{name} must contain finite positive values"
        )
    return arr


def _normalize_estimate(
    value: Any,
    *,
    n_frequency: int,
    tail: tuple[int, int],
    name: str,
) -> np.ndarray:
    arr = np.asarray(value)
    if arr.ndim == 2 and arr.shape == tail:
        arr = arr[None, ...]
    if arr.ndim != 3 or arr.shape != (n_frequency, *tail):
        raise MobileMTValidationError(
            f"{name} must have shape {(n_frequency, *tail)}, "
            f"got {arr.shape}"
        )
    if arr.dtype.kind not in "biufc":
        raise MobileMTValidationError(f"{name} must be numeric")
    return arr


def _reference_mapping(
    reference_station: MobileMTReferenceStation | None,
) -> dict[str, Any] | None:
    if reference_station is None:
        return None
    out: dict[str, Any] = {
        "station_id": reference_station.preferred_id,
        "electric_channels": list(reference_station.electric_channels),
    }
    if reference_station.site is not None:
        site = reference_station.site
        out["site"] = site.to_dict(max_depth=2)
    if reference_station.attrs:
        out["attrs"] = dict(reference_station.attrs)
    return out




def _xml_notes_mapping(
    spec: MobileMTSystemSpec,
    reference_station: MobileMTReferenceStation | None,
) -> dict[str, Any]:
    low, high = spec.nominal_frequency_range_hz
    mobile: dict[str, Any] = {
        "NominalFrequencyMinHz": low,
        "NominalFrequencyMaxHz": high,
        "NominalMaxFrequencyWindows": spec.nominal_max_frequency_windows,
        "NominalSamplingRateHz": spec.nominal_sampling_rate_hz,
        "InputChannels": ",".join(spec.input_channels),
        "OutputChannels": ",".join(spec.output_channels),
    }
    if reference_station is not None:
        if reference_station.preferred_id is not None:
            mobile["ReferenceStationId"] = reference_station.preferred_id
        site = reference_station.site
        location = site.location if site is not None else None
        if location is not None:
            if location.latitude is not None:
                mobile["ReferenceLatitude"] = location.latitude
            if location.longitude is not None:
                mobile["ReferenceLongitude"] = location.longitude
            if location.elevation is not None:
                mobile["ReferenceElevation"] = location.elevation
            if location.datum is not None:
                mobile["ReferenceDatum"] = location.datum
    return {"MobileMT": mobile}


def validate_mobilemt_transfer_function(
    tf: TransferFunction,
) -> TransferFunction:
    """Validate and return a MobileMT 3x2 admittance transfer function."""
    if not isinstance(tf, TransferFunction):
        raise TypeError("tf must be a TransferFunction")
    if tf.name != MOBILEMT_ADMITTANCE_TAG:
        raise MobileMTValidationError(
            "MobileMT transfer function must use the "
            f"{MOBILEMT_ADMITTANCE_TAG!r} datatype"
        )
    if tuple(tf.input_channels) != MOBILEMT_INPUT_CHANNELS:
        raise MobileMTValidationError(
            "MobileMT admittance inputs must be ('Ex', 'Ey')"
        )
    if tuple(tf.output_channels) != MOBILEMT_OUTPUT_CHANNELS:
        raise MobileMTValidationError(
            "MobileMT admittance outputs must be ('Hx', 'Hy', 'Hz')"
        )
    if tf.data.shape[1:] != (3, 2):
        raise MobileMTValidationError(
            "MobileMT admittance must have matrix shape (3, 2)"
        )
    return tf


def build_mobilemt_emtf(
    admittance: Any,
    *,
    frequency: Any | None = None,
    periods: Any | None = None,
    units: str | None = None,
    variance: Any | None = None,
    inverse_signal_covariance: Any | None = None,
    residual_covariance: Any | None = None,
    product_id: str | None = None,
    description: str | None = None,
    reference_station: MobileMTReferenceStation | None = None,
    system_spec: MobileMTSystemSpec | None = None,
    attrs: Mapping[str, Any] | None = None,
) -> EMTF:
    """Build one sample-level :class:`EMTF` MobileMT response.

    Parameters
    ----------
    admittance : array-like
        Complex admittance with shape ``(nf, 3, 2)`` or one ``(3, 2)``
        matrix. Rows are ``Hx, Hy, Hz`` and columns are ``Ex, Ey``.
    frequency, periods : array-like, optional
        Exactly one frequency or period grid must be supplied.
    units : str, optional
        Units of the supplied admittance. No raw-unit convention is invented
        by this adapter because processed deliveries may differ.
    variance : array-like, optional
        Component variance, shape ``(nf, 3, 2)``.
    inverse_signal_covariance : array-like, optional
        Input covariance factor ``S``, shape ``(nf, 2, 2)``.
    residual_covariance : array-like, optional
        Output residual covariance ``N``, shape ``(nf, 3, 3)``.
    reference_station : MobileMTReferenceStation, optional
        Ground electric reference metadata retained as MobileMT metadata.
    system_spec : MobileMTSystemSpec, optional
        System-description metadata; defaults to published nominal values.

    Notes
    -----
    Apparent conductivity is intentionally not computed here. The published
    system exposes it as a processed output, but a verified delivery schema is
    required before pyCSAMT should codify its exact exported representation.
    """
    register_mobilemt_datatypes()
    if (frequency is None) == (periods is None):
        raise MobileMTValidationError(
            "exactly one of frequency or periods must be supplied"
        )

    if frequency is not None:
        freq = _as_frequency(frequency)
        period_arr = 1.0 / freq
    else:
        period_arr = _as_frequency(periods, name="periods")
        freq = 1.0 / period_arr

    data = np.asarray(admittance)
    if data.ndim == 2 and data.shape == (3, 2):
        data = data[None, ...]
    if data.ndim != 3 or data.shape[1:] != (3, 2):
        raise MobileMTValidationError(
            "admittance must have shape (nf, 3, 2) or (3, 2)"
        )
    if data.shape[0] != freq.size:
        raise MobileMTValidationError(
            "frequency count does not match admittance: "
            f"{freq.size} != {data.shape[0]}"
        )
    if data.dtype.kind not in "biufc":
        raise MobileMTValidationError("admittance must be numeric")

    tf = TransferFunction(
        name=MOBILEMT_ADMITTANCE_TAG,
        data=data,
        input_channels=MOBILEMT_INPUT_CHANNELS,
        output_channels=MOBILEMT_OUTPUT_CHANNELS,
        units=units,
        periods=period_arr,
        attrs={"technology": "MobileMT"},
    )

    if variance is not None:
        tf.add_estimate(
            StatisticalEstimate(
                name="VAR",
                kind="variance",
                data=_normalize_estimate(
                    variance,
                    n_frequency=freq.size,
                    tail=(3, 2),
                    name="variance",
                ),
            )
        )
    if inverse_signal_covariance is not None:
        tf.add_estimate(
            StatisticalEstimate(
                name="INVSIGCOV",
                kind="inverse_signal_covariance",
                data=_normalize_estimate(
                    inverse_signal_covariance,
                    n_frequency=freq.size,
                    tail=(2, 2),
                    name="inverse_signal_covariance",
                ),
            )
        )
    if residual_covariance is not None:
        tf.add_estimate(
            StatisticalEstimate(
                name="RESIDCOV",
                kind="residual_covariance",
                data=_normalize_estimate(
                    residual_covariance,
                    n_frequency=freq.size,
                    tail=(3, 3),
                    name="residual_covariance",
                ),
            )
        )

    spec = system_spec or MobileMTSystemSpec()
    if not isinstance(spec, MobileMTSystemSpec):
        raise TypeError("system_spec must be a MobileMTSystemSpec or None")
    if reference_station is not None and not isinstance(
        reference_station,
        MobileMTReferenceStation,
    ):
        raise TypeError(
            "reference_station must be MobileMTReferenceStation or None"
        )

    mobilemt_meta: dict[str, Any] = {
        "system": spec.to_dict(max_depth=2),
        "reference_station": _reference_mapping(reference_station),
    }
    document_attrs = dict(attrs or {})
    document_attrs["mobilemt"] = mobilemt_meta

    doc = EMTF(
        product_id=product_id,
        description=(
            description
            or "MobileMT three-component magnetic admittance response"
        ),
        subtype="mobilemt",
        tags=("mobilemt", MOBILEMT_ADMITTANCE_TAG),
        periods=period_arr,
        metadata={
            "notes": _xml_notes_mapping(spec, reference_station),
        },
        attrs=document_attrs,
    )
    doc.add_transfer_function(tf)
    validate_mobilemt_transfer_function(tf)
    return doc


def build_mobilemt_record(
    sample_id: str,
    admittance: Any,
    *,
    frequency: Any | None = None,
    periods: Any | None = None,
    apparent_conductivity: Any | None = None,
    apparent_conductivity_units: str = (
        MOBILEMT_APPARENT_CONDUCTIVITY_UNITS
    ),
    quality: Mapping[str, Any] | None = None,
    record_attrs: Mapping[str, Any] | None = None,
    **emtf_kwargs: Any,
) -> AirborneEMRecord:
    """Build one airborne record from decoded MobileMT arrays."""
    doc = build_mobilemt_emtf(
        admittance,
        frequency=frequency,
        periods=periods,
        product_id=emtf_kwargs.pop("product_id", str(sample_id)),
        **emtf_kwargs,
    )
    fields: dict[str, Any] = {}
    attrs = dict(record_attrs or {})
    if apparent_conductivity is not None:
        sigma = np.asarray(apparent_conductivity, dtype=float)
        if sigma.ndim == 0:
            sigma = sigma.reshape(1)
        if sigma.ndim != 1 or sigma.size != doc.n_periods:
            raise MobileMTValidationError(
                "apparent_conductivity must have one value per frequency"
            )
        if np.any(np.isinf(sigma)):
            raise MobileMTValidationError(
                "apparent_conductivity must not contain infinite values"
            )
        fields[MOBILEMT_APPARENT_CONDUCTIVITY_FIELD] = sigma
        attrs["apparent_conductivity_units"] = str(
            apparent_conductivity_units
        ).strip()

    return AirborneEMRecord(
        sample_id=str(sample_id),
        emtf=doc,
        fields=fields,
        quality=dict(quality or {}),
        attrs=attrs,
    )


def _sample_axis_array(
    value: Any,
    *,
    name: str,
    n_samples: int,
    tail_ndim: int,
) -> np.ndarray:
    arr = np.asarray(value)
    if n_samples == 1 and arr.ndim == tail_ndim:
        arr = arr[None, ...]
    if arr.ndim != tail_ndim + 1 or arr.shape[0] != n_samples:
        raise MobileMTValidationError(
            f"{name} must have a leading sample axis of length {n_samples}"
        )
    return arr


def build_mobilemt_line(
    line_id: str,
    navigation: NavigationTrack,
    admittance: Any,
    *,
    frequency: Any,
    record_mask: Any | None = None,
    apparent_conductivity: Any | None = None,
    variance: Any | None = None,
    inverse_signal_covariance: Any | None = None,
    residual_covariance: Any | None = None,
    units: str | None = None,
    reference_station: MobileMTReferenceStation | None = None,
    system_spec: MobileMTSystemSpec | None = None,
    attrs: Mapping[str, Any] | None = None,
) -> AirborneEMLine:
    """Build one MobileMT flight line from decoded sample-aligned arrays.

    ``frequency`` may be a common ``(nf,)`` vector or a sample-specific
    ``(n_samples, nf)`` matrix. ``record_mask`` controls sparse EM coverage;
    navigation points are never deleted when their EM record is absent.
    """
    if not isinstance(navigation, NavigationTrack):
        raise TypeError("navigation must be a NavigationTrack")
    n_samples = navigation.n_samples
    data = _sample_axis_array(
        admittance,
        name="admittance",
        n_samples=n_samples,
        tail_ndim=3,
    )
    if data.shape[2:] != (3, 2):
        raise MobileMTValidationError(
            "line admittance must have shape (samples, nf, 3, 2)"
        )
    n_frequency = data.shape[1]

    freq = np.asarray(frequency, dtype=float)
    if freq.ndim == 1:
        common_frequency = _as_frequency(freq)
        if common_frequency.size != n_frequency:
            raise MobileMTValidationError(
                "shared frequency length does not match admittance"
            )
        frequency_rows = None
    elif freq.ndim == 2 and freq.shape == (n_samples, n_frequency):
        common_frequency = None
        frequency_rows = freq
    else:
        raise MobileMTValidationError(
            "frequency must have shape (nf,) or (n_samples, nf)"
        )

    if record_mask is None:
        mask = np.ones(n_samples, dtype=bool)
    else:
        mask = np.asarray(record_mask, dtype=bool)
        if mask.ndim != 1 or mask.size != n_samples:
            raise MobileMTValidationError(
                "record_mask must have one boolean per navigation sample"
            )

    sigma_all = None
    if apparent_conductivity is not None:
        sigma_all = np.asarray(apparent_conductivity, dtype=float)
        if n_samples == 1 and sigma_all.ndim == 1:
            sigma_all = sigma_all[None, ...]
        if sigma_all.shape != (n_samples, n_frequency):
            raise MobileMTValidationError(
                "apparent_conductivity must have shape (samples, nf)"
            )

    def _optional_sample(value: Any, tail: tuple[int, int], name: str):
        if value is None:
            return None
        arr = _sample_axis_array(
            value,
            name=name,
            n_samples=n_samples,
            tail_ndim=3,
        )
        expected = (n_samples, n_frequency, *tail)
        if arr.shape != expected:
            raise MobileMTValidationError(
                f"{name} must have shape {expected}, got {arr.shape}"
            )
        return arr

    var_all = _optional_sample(variance, (3, 2), "variance")
    inv_all = _optional_sample(
        inverse_signal_covariance,
        (2, 2),
        "inverse_signal_covariance",
    )
    res_all = _optional_sample(
        residual_covariance,
        (3, 3),
        "residual_covariance",
    )

    line_attrs = dict(attrs or {})
    line_attrs.setdefault("technology", "MobileMT")
    line_attrs.setdefault(
        "reference_station",
        _reference_mapping(reference_station),
    )
    line = AirborneEMLine(
        line_id=str(line_id),
        navigation=navigation,
        attrs=line_attrs,
    )

    for index, sample_id in enumerate(navigation.sample_ids):
        if not mask[index]:
            continue
        sample_frequency = (
            common_frequency
            if common_frequency is not None
            else _as_frequency(
                frequency_rows[index],
                name=f"frequency[{sample_id}]",
            )
        )
        record = build_mobilemt_record(
            sample_id,
            data[index],
            frequency=sample_frequency,
            apparent_conductivity=(
                None if sigma_all is None else sigma_all[index]
            ),
            variance=None if var_all is None else var_all[index],
            inverse_signal_covariance=(
                None if inv_all is None else inv_all[index]
            ),
            residual_covariance=(
                None if res_all is None else res_all[index]
            ),
            units=units,
            reference_station=reference_station,
            system_spec=system_spec,
            product_id=f"{line.line_id}.{sample_id}",
            record_attrs={"navigation_index": index},
        )
        line.add_record(record)
    return line


def build_mobilemt_dataset(
    name: str,
    lines: Any,
    *,
    survey: SurveyMeta | None = None,
    system_spec: MobileMTSystemSpec | None = None,
    instrument_serial: str | None = None,
    software_version: str = "",
    attrs: Mapping[str, Any] | None = None,
) -> AirborneEMDataset:
    """Build a common airborne dataset from already constructed MobileMT lines.

    This helper does not parse files. It only supplies technology metadata and
    reuses :class:`AirborneEMDataset` as the survey-level container.
    """
    if survey is not None and not isinstance(survey, SurveyMeta):
        raise TypeError("survey must be a SurveyMeta or None")
    spec = system_spec or MobileMTSystemSpec()
    if not isinstance(spec, MobileMTSystemSpec):
        raise TypeError("system_spec must be a MobileMTSystemSpec or None")

    dataset_attrs = dict(attrs or {})
    dataset_attrs.setdefault("technology", "MobileMT")
    dataset_attrs.setdefault("mobilemt_system", spec.to_dict(max_depth=2))
    dataset = AirborneEMDataset(
        name=str(name),
        survey=survey,
        instrument=spec.to_instrument_meta(
            serial=instrument_serial,
            software_version=software_version,
        ),
        method="AEM",
        attrs=dataset_attrs,
    )

    iterable = lines.values() if isinstance(lines, Mapping) else lines
    for line in iterable:
        if not isinstance(line, AirborneEMLine):
            raise TypeError("lines must contain AirborneEMLine instances")
        technology = str(line.attrs.get("technology", "")).strip().lower()
        if technology and technology != "mobilemt":
            raise MobileMTValidationError(
                f"line {line.line_id!r} is tagged as {technology!r}"
            )
        dataset.add_line(line)
    return dataset
