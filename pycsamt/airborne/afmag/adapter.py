# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

"""Scientific adapter contracts for historical and tensor AFMAG products.

No native vendor/archive parser is defined here.  The module accepts decoded
arrays and maps two distinct AFMAG generations into the common pyCSAMT model:

* original comparator AFMAG -> scalar line-direction tilt response;
* tensor AFMAG / AirMt -> 3 x 2 interstation magnetic transfer function plus
  the optional rotationally invariant amplification parameter.
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

import numpy as np

from ...emtf.document import EMTF
from ...emtf.estimates import StatisticalEstimate
from ...emtf.transfer import TransferFunction
from ...metadata import (
    OrientationMeta,
    ProcessingMeta,
    RemoteReferenceMeta,
    SiteMeta,
    SurveyMeta,
)
from ..base import AirborneEMDataset, AirborneEMLine, AirborneEMRecord
from ..navigation import NavigationTrack
from .base import (
    AFMAGReferenceStation,
    AirMtSystemSpec,
    OriginalAFMAGSystemSpec,
)
from .constants import (
    AFMAG_AP_TAG,
    AFMAG_AP_UNITS,
    AFMAG_TENSOR_COMPONENTS,
    AFMAG_TENSOR_INPUT_CHANNELS,
    AFMAG_TENSOR_OUTPUT_CHANNELS,
    AFMAG_TENSOR_TAG,
    AFMAG_TENSOR_UNITS,
    AFMAG_TILT_DEFAULT_UNITS,
    AFMAG_TILT_TAG,
)
from .datatypes import register_afmag_datatypes

__all__ = [
    "AFMAGValidationError",
    "compute_airmt_amplification_parameter",
    "validate_airmt_transfer_function",
    "validate_original_afmag_tilt",
    "build_airmt_emtf",
    "build_airmt_record",
    "build_airmt_line",
    "build_airmt_dataset",
    "build_original_afmag_emtf",
    "build_original_afmag_record",
    "build_original_afmag_line",
    "build_original_afmag_dataset",
]


class AFMAGValidationError(ValueError):
    """Raised when decoded AFMAG scientific arrays are inconsistent."""


def _as_frequency(value: Any, *, name: str = "frequency") -> np.ndarray:
    arr = np.asarray(value, dtype=float)
    if arr.ndim == 0:
        arr = arr.reshape(1)
    if arr.ndim != 1:
        raise AFMAGValidationError(f"{name} must be a 1-D array")
    if arr.size == 0:
        raise AFMAGValidationError(f"{name} must be non-empty")
    if not np.all(np.isfinite(arr)) or np.any(arr <= 0.0):
        raise AFMAGValidationError(
            f"{name} must contain finite positive values"
        )
    return arr


def _period_axis(
    *,
    frequency: Any | None,
    periods: Any | None,
) -> tuple[np.ndarray, np.ndarray]:
    if (frequency is None) == (periods is None):
        raise AFMAGValidationError(
            "exactly one of frequency or periods must be supplied"
        )
    if frequency is not None:
        freq = _as_frequency(frequency)
        return freq, 1.0 / freq
    period_arr = _as_frequency(periods, name="periods")
    return 1.0 / period_arr, period_arr


def _normalize_tensor(value: Any) -> np.ndarray:
    arr = np.asarray(value)
    if arr.ndim == 2 and arr.shape == (3, 2):
        arr = arr[None, ...]
    if arr.ndim != 3 or arr.shape[1:] != (3, 2):
        raise AFMAGValidationError(
            "AirMt tensor must have shape (nf, 3, 2) or (3, 2)"
        )
    if arr.dtype.kind not in "biufc":
        raise AFMAGValidationError("AirMt tensor must be numeric")
    return arr


def _normalize_tilt(value: Any) -> np.ndarray:
    arr = np.asarray(value)
    if arr.ndim == 0:
        arr = arr.reshape(1)
    if arr.ndim != 1:
        raise AFMAGValidationError(
            "original AFMAG tilt must be scalar or a 1-D array"
        )
    if arr.dtype.kind not in "biufc":
        raise AFMAGValidationError("original AFMAG tilt must be numeric")
    if np.iscomplexobj(arr) and np.any(np.imag(arr) != 0.0):
        raise AFMAGValidationError("original AFMAG tilt must be real-valued")
    return np.asarray(np.real(arr), dtype=float)


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
        raise AFMAGValidationError(
            f"{name} must have shape {(n_frequency, *tail)}, "
            f"got {arr.shape}"
        )
    if arr.dtype.kind not in "biufc":
        raise AFMAGValidationError(f"{name} must be numeric")
    return arr


def _normalize_scalar_variance(
    value: Any,
    *,
    n_frequency: int,
) -> np.ndarray:
    arr = np.asarray(value)
    if arr.ndim == 0:
        arr = arr.reshape(1)
    if arr.ndim == 1:
        arr = arr[:, None, None]
    elif arr.ndim == 2 and arr.shape == (1, 1):
        arr = arr[None, ...]
    if arr.shape != (n_frequency, 1, 1):
        raise AFMAGValidationError(
            "tilt variance must have shape "
            f"{(n_frequency, 1, 1)} or ({n_frequency},)"
        )
    if np.iscomplexobj(arr) and np.any(np.imag(arr) != 0.0):
        raise AFMAGValidationError("tilt variance must be real-valued")
    return np.asarray(np.real(arr), dtype=float)


def compute_airmt_amplification_parameter(
    tensor: Any,
    *,
    zero_policy: str = "nan",
) -> np.ndarray:
    """Compute the AirMt rotation-invariant complex amplification parameter.

    For the two column vectors ``T1`` and ``T2`` of a 3 x 2 magnetic transfer
    function, the implemented patent formulation is::

        Kvec = T1 x T2
        AP = Kvec . Re(Kvec) / |Re(Kvec)|

    Parameters
    ----------
    tensor : array-like
        One ``(3, 2)`` or ``(3, 3)`` matrix, or a family with those
        trailing dimensions.  The patent parameter uses the first two
        column vectors in either case.
    zero_policy : {"nan", "raise"}
        Behavior when ``|Re(Kvec)|`` is zero.  Such a sample has no defined
        projection direction for this parameter.

    Returns
    -------
    numpy.ndarray
        Complex AP with shape equal to the leading tensor dimensions.
    """
    arr = np.asarray(tensor)
    if (
        arr.ndim < 2
        or arr.shape[-2] != 3
        or arr.shape[-1] not in {2, 3}
    ):
        raise AFMAGValidationError(
            "AirMt amplification parameter requires (..., 3, 2) or "
            "(..., 3, 3) data"
        )
    if arr.dtype.kind not in "biufc":
        raise AFMAGValidationError("AirMt tensor must be numeric")

    policy = str(zero_policy).strip().lower()
    if policy not in {"nan", "raise"}:
        raise ValueError("zero_policy must be 'nan' or 'raise'")

    first = arr[..., :, 0]
    second = arr[..., :, 1]
    cross = np.cross(first, second, axis=-1)
    real_cross = np.real(cross)
    norm = np.linalg.norm(real_cross, axis=-1)
    numerator = np.sum(cross * real_cross, axis=-1)

    bad = np.isfinite(norm) & (norm <= np.finfo(float).eps)
    if policy == "raise" and np.any(bad):
        raise AFMAGValidationError(
            "AirMt amplification parameter is undefined where "
            "|Re(T1 x T2)| is zero"
        )

    result = np.full(norm.shape, np.nan + 1j * np.nan, dtype=complex)
    good = np.isfinite(norm) & (norm > np.finfo(float).eps)
    result[good] = numerator[good] / norm[good]
    return result


def validate_airmt_transfer_function(
    tf: TransferFunction,
) -> TransferFunction:
    """Validate and return an AirMt 3 x 2 interstation magnetic TF."""
    if not isinstance(tf, TransferFunction):
        raise TypeError("tf must be a TransferFunction")
    if tf.name != AFMAG_TENSOR_TAG:
        raise AFMAGValidationError(
            "AirMt must use the standard EMTF interstation magnetic TF"
        )
    if tuple(tf.input_channels) != AFMAG_TENSOR_INPUT_CHANNELS:
        raise AFMAGValidationError("AirMt transfer inputs must be Hx/Hy")
    if tuple(tf.output_channels) != AFMAG_TENSOR_OUTPUT_CHANNELS:
        raise AFMAGValidationError(
            "AirMt airborne outputs must be Hx/Hy/Hz"
        )
    if tf.data.shape[1:] != (3, 2):
        raise AFMAGValidationError(
            "AirMt transfer function must have matrix shape (3, 2)"
        )
    return tf


def validate_original_afmag_tilt(
    tf: TransferFunction,
) -> TransferFunction:
    """Validate and return an original AFMAG scalar tilt response."""
    if not isinstance(tf, TransferFunction):
        raise TypeError("tf must be a TransferFunction")
    if tf.name != AFMAG_TILT_TAG:
        raise AFMAGValidationError(
            "original AFMAG response must use the afmag_tilt datatype"
        )
    if tf.input_channels or tf.output_channels:
        raise AFMAGValidationError("original AFMAG tilt must be scalar")
    if tf.data.shape[1:] != (1, 1):
        raise AFMAGValidationError("original AFMAG tilt must be scalar")
    if np.iscomplexobj(tf.data) and np.any(np.imag(tf.data) != 0.0):
        raise AFMAGValidationError("original AFMAG tilt must be real")
    return tf


def _reference_mapping(
    reference_station: AFMAGReferenceStation | None,
) -> dict[str, Any] | None:
    if reference_station is None:
        return None
    result: dict[str, Any] = {
        "station_id": reference_station.preferred_id,
        "measured_channels": list(reference_station.measured_channels),
        "transfer_input_channels": list(
            reference_station.transfer_input_channels
        ),
    }
    if reference_station.site is not None:
        result["site"] = reference_station.site.to_dict(max_depth=2)
    if reference_station.attrs:
        result["attrs"] = dict(reference_station.attrs)
    return result


def _airmt_notes(
    spec: AirMtSystemSpec,
    reference_station: AFMAGReferenceStation | None,
) -> dict[str, Any]:
    low, high = spec.practical_frequency_range_hz
    values: dict[str, Any] = {
        "Family": "tensor_afmag_airmt",
        "PracticalFrequencyMinHz": low,
        "PracticalFrequencyMaxHz": high,
        "TypicalFrequencyCountMin": spec.typical_frequency_count[0],
        "TypicalFrequencyCountMax": spec.typical_frequency_count[1],
        "TimeSeriesSamplingRateHz": spec.time_series_sampling_rate_hz,
        "TransferInputChannels": ",".join(spec.input_channels),
        "AirborneOutputChannels": ",".join(spec.output_channels),
        "ReferenceMeasuredChannels": ",".join(spec.reference_channels),
        "TransferFunction": (
            "[Hx,Hy,Hz](r) = TI * [Hx,Hy](r0)"
        ),
    }
    if reference_station is not None:
        if reference_station.preferred_id is not None:
            values["ReferenceStationId"] = reference_station.preferred_id
        site = reference_station.site
        location = site.location if site is not None else None
        if location is not None:
            if location.latitude is not None:
                values["ReferenceLatitude"] = location.latitude
            if location.longitude is not None:
                values["ReferenceLongitude"] = location.longitude
            if location.elevation is not None:
                values["ReferenceElevation"] = location.elevation
            if location.datum is not None:
                values["ReferenceDatum"] = location.datum
    return {"AFMAG": values}


def _original_notes(
    spec: OriginalAFMAGSystemSpec,
    *,
    response_kind: str,
) -> dict[str, Any]:
    low, high = spec.historical_frequency_band_hz
    values: dict[str, Any] = {
        "Family": "original_comparator_afmag",
        "HistoricalFrequencyMinHz": low,
        "HistoricalFrequencyMaxHz": high,
        "TypicalFrequenciesHz": ",".join(
            f"{value:g}" for value in spec.typical_frequencies_hz
        ),
        "CoilCount": spec.coil_count,
        "CoilTiltDeg": spec.coil_tilt_deg,
        "CoilSeparationDeg": spec.coil_separation_deg,
        "DigitalRecording": spec.digital_recording,
        "ResponseKind": response_kind,
    }
    return {"AFMAG": values}


def _processing_for_reference(
    reference_station: AFMAGReferenceStation | None,
    processing: ProcessingMeta | None,
) -> ProcessingMeta | None:
    if processing is not None and not isinstance(processing, ProcessingMeta):
        raise TypeError("processing must be a ProcessingMeta or None")
    if reference_station is None:
        return processing

    remote = RemoteReferenceMeta(
        reference_type="fixed_ground_magnetic",
        site=reference_station.preferred_id,
        extra={
            "measured_channels": list(reference_station.measured_channels),
            "transfer_input_channels": list(
                reference_station.transfer_input_channels
            ),
            "technology": "AirMt",
        },
    )
    if processing is None:
        return ProcessingMeta(remote_reference=remote)
    if processing.remote_reference is None:
        return ProcessingMeta(
            sign_convention=processing.sign_convention,
            processed_by=processing.processed_by,
            software=processing.software,
            remote_reference=remote,
            processing_tag=processing.processing_tag,
            run_list=(
                None
                if processing.run_list is None
                else list(processing.run_list)
            ),
            extra=dict(processing.extra),
        )

    existing = processing.remote_reference
    reference_id = reference_station.preferred_id
    if existing.site is not None and reference_id is not None:
        if str(existing.site) != str(reference_id):
            raise AFMAGValidationError(
                "processing remote-reference site conflicts with "
                "AFMAGReferenceStation"
            )
    return processing


def build_airmt_emtf(
    tensor: Any,
    *,
    frequency: Any | None = None,
    periods: Any | None = None,
    units: str | None = AFMAG_TENSOR_UNITS,
    variance: Any | None = None,
    inverse_signal_covariance: Any | None = None,
    residual_covariance: Any | None = None,
    include_amplification_parameter: bool = True,
    amplification_zero_policy: str = "nan",
    product_id: str | None = None,
    description: str | None = None,
    reference_station: AFMAGReferenceStation | None = None,
    system_spec: AirMtSystemSpec | None = None,
    site: SiteMeta | None = None,
    orientation: OrientationMeta | None = None,
    processing: ProcessingMeta | None = None,
    attrs: Mapping[str, Any] | None = None,
) -> EMTF:
    """Build one sample-level tensor AFMAG/AirMt :class:`EMTF` response."""
    register_afmag_datatypes()
    freq, period_arr = _period_axis(frequency=frequency, periods=periods)
    data = _normalize_tensor(tensor)
    if data.shape[0] != freq.size:
        raise AFMAGValidationError(
            "frequency count does not match AirMt tensor: "
            f"{freq.size} != {data.shape[0]}"
        )

    tf = TransferFunction(
        name=AFMAG_TENSOR_TAG,
        data=data,
        input_channels=AFMAG_TENSOR_INPUT_CHANNELS,
        output_channels=AFMAG_TENSOR_OUTPUT_CHANNELS,
        units=units,
        periods=period_arr,
        attrs={
            "technology": "AirMt",
            "afmag_family": "tensor",
            "reference_geometry": (
                "airborne_Hx_Hy_Hz_to_fixed_ground_Hx_Hy"
            ),
            "components": list(AFMAG_TENSOR_COMPONENTS),
            "axis_convention": "as_supplied",
        },
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

    spec = system_spec or AirMtSystemSpec()
    if not isinstance(spec, AirMtSystemSpec):
        raise TypeError("system_spec must be an AirMtSystemSpec or None")
    if reference_station is not None and not isinstance(
        reference_station,
        AFMAGReferenceStation,
    ):
        raise TypeError(
            "reference_station must be AFMAGReferenceStation or None"
        )
    if site is not None and not isinstance(site, SiteMeta):
        raise TypeError("site must be a SiteMeta or None")
    if orientation is not None and not isinstance(
        orientation,
        OrientationMeta,
    ):
        raise TypeError("orientation must be an OrientationMeta or None")

    document_attrs = dict(attrs or {})
    document_attrs["afmag"] = {
        "family": "tensor_airmt",
        "system": spec.to_dict(max_depth=2),
        "reference_station": _reference_mapping(reference_station),
    }

    doc = EMTF(
        product_id=product_id,
        description=(
            description
            or "Tensor AFMAG/AirMt interstation magnetic response"
        ),
        subtype="afmag_airmt",
        tags=("afmag", "airmt", AFMAG_TENSOR_TAG),
        periods=period_arr,
        site=site,
        orientation=orientation,
        processing=_processing_for_reference(reference_station, processing),
        metadata={"notes": _airmt_notes(spec, reference_station)},
        attrs=document_attrs,
    )
    doc.add_transfer_function(tf)
    validate_airmt_transfer_function(tf)

    if include_amplification_parameter:
        ap = compute_airmt_amplification_parameter(
            data,
            zero_policy=amplification_zero_policy,
        )
        ap_tf = TransferFunction(
            name=AFMAG_AP_TAG,
            data=ap,
            input_channels=(),
            output_channels=(),
            units=AFMAG_AP_UNITS,
            periods=period_arr,
            attrs={
                "technology": "AirMt",
                "derived_from": AFMAG_TENSOR_TAG,
                "rotation_invariant": True,
                "formula": (
                    "K=T1xT2; AP=K.Re(K)/|Re(K)|"
                ),
            },
        )
        doc.add_transfer_function(ap_tf)
    return doc


def build_original_afmag_emtf(
    tilt: Any,
    *,
    frequency: Any | None = None,
    periods: Any | None = None,
    response_kind: str = "tilt_angle",
    units: str | None = None,
    variance: Any | None = None,
    product_id: str | None = None,
    description: str | None = None,
    system_spec: OriginalAFMAGSystemSpec | None = None,
    site: SiteMeta | None = None,
    orientation: OrientationMeta | None = None,
    attrs: Mapping[str, Any] | None = None,
) -> EMTF:
    """Build one historical comparator-AFMAG scalar tilt response.

    ``response_kind='tilt_angle'`` represents a calibrated tilt angle and
    defaults to degrees. ``'comparator_deflection'`` preserves an archival
    instrument deflection without pretending it is an angle; units should then
    be supplied when known.
    """
    register_afmag_datatypes()
    freq, period_arr = _period_axis(frequency=frequency, periods=periods)
    data = _normalize_tilt(tilt)
    if data.size != freq.size:
        raise AFMAGValidationError(
            "frequency count does not match original AFMAG tilt: "
            f"{freq.size} != {data.size}"
        )

    kind = str(response_kind).strip().lower().replace("-", "_")
    aliases = {
        "tilt": "tilt_angle",
        "angle": "tilt_angle",
        "deflection": "comparator_deflection",
    }
    kind = aliases.get(kind, kind)
    if kind not in {"tilt_angle", "comparator_deflection"}:
        raise ValueError(
            "response_kind must be 'tilt_angle' or "
            "'comparator_deflection'"
        )
    if units is None and kind == "tilt_angle":
        units = AFMAG_TILT_DEFAULT_UNITS

    tf = TransferFunction(
        name=AFMAG_TILT_TAG,
        data=data,
        input_channels=(),
        output_channels=(),
        units=units,
        periods=period_arr,
        attrs={
            "technology": "AFMAG",
            "afmag_family": "original_comparator",
            "response_kind": kind,
            "line_direction_response": True,
            "calibrated_angle": kind == "tilt_angle",
        },
    )
    if variance is not None:
        tf.add_estimate(
            StatisticalEstimate(
                name="VAR",
                kind="variance",
                data=_normalize_scalar_variance(
                    variance,
                    n_frequency=freq.size,
                ),
                units=(None if units is None else f"{units}^2"),
            )
        )

    spec = system_spec or OriginalAFMAGSystemSpec()
    if not isinstance(spec, OriginalAFMAGSystemSpec):
        raise TypeError(
            "system_spec must be OriginalAFMAGSystemSpec or None"
        )
    if site is not None and not isinstance(site, SiteMeta):
        raise TypeError("site must be a SiteMeta or None")
    if orientation is not None and not isinstance(
        orientation,
        OrientationMeta,
    ):
        raise TypeError("orientation must be an OrientationMeta or None")

    document_attrs = dict(attrs or {})
    document_attrs["afmag"] = {
        "family": "original_comparator",
        "system": spec.to_dict(max_depth=2),
    }
    doc = EMTF(
        product_id=product_id,
        description=(
            description or "Original comparator AFMAG tilt response"
        ),
        subtype="afmag_original",
        tags=("afmag", AFMAG_TILT_TAG),
        periods=period_arr,
        site=site,
        orientation=orientation,
        metadata={
            "notes": _original_notes(spec, response_kind=kind),
        },
        attrs=document_attrs,
    )
    doc.add_transfer_function(tf)
    validate_original_afmag_tilt(tf)
    return doc


def build_airmt_record(
    sample_id: str,
    tensor: Any,
    *,
    frequency: Any | None = None,
    periods: Any | None = None,
    fields: Mapping[str, Any] | None = None,
    quality: Mapping[str, Any] | None = None,
    record_attrs: Mapping[str, Any] | None = None,
    **emtf_kwargs: Any,
) -> AirborneEMRecord:
    """Build one airborne record from decoded tensor AFMAG/AirMt data."""
    doc = build_airmt_emtf(
        tensor,
        frequency=frequency,
        periods=periods,
        product_id=emtf_kwargs.pop("product_id", str(sample_id)),
        **emtf_kwargs,
    )
    return AirborneEMRecord(
        sample_id=str(sample_id),
        emtf=doc,
        fields=dict(fields or {}),
        quality=dict(quality or {}),
        attrs=dict(record_attrs or {}),
    )


def build_original_afmag_record(
    sample_id: str,
    tilt: Any,
    *,
    frequency: Any | None = None,
    periods: Any | None = None,
    fields: Mapping[str, Any] | None = None,
    quality: Mapping[str, Any] | None = None,
    record_attrs: Mapping[str, Any] | None = None,
    **emtf_kwargs: Any,
) -> AirborneEMRecord:
    """Build one record from a historical AFMAG tilt response."""
    doc = build_original_afmag_emtf(
        tilt,
        frequency=frequency,
        periods=periods,
        product_id=emtf_kwargs.pop("product_id", str(sample_id)),
        **emtf_kwargs,
    )
    return AirborneEMRecord(
        sample_id=str(sample_id),
        emtf=doc,
        fields=dict(fields or {}),
        quality=dict(quality or {}),
        attrs=dict(record_attrs or {}),
    )


def _line_frequency_rows(
    frequency: Any,
    *,
    n_samples: int,
    n_frequency: int,
) -> tuple[np.ndarray | None, np.ndarray | None]:
    freq = np.asarray(frequency, dtype=float)
    if freq.ndim == 1:
        common = _as_frequency(freq)
        if common.size != n_frequency:
            raise AFMAGValidationError(
                "shared frequency length does not match response"
            )
        return common, None
    if freq.ndim == 2 and freq.shape == (n_samples, n_frequency):
        if not np.all(np.isfinite(freq)) or np.any(freq <= 0.0):
            raise AFMAGValidationError(
                "per-sample frequency grids must be finite and positive"
            )
        return None, freq
    raise AFMAGValidationError(
        "frequency must have shape (nf,) or (n_samples, nf)"
    )


def _record_mask(value: Any | None, *, n_samples: int) -> np.ndarray:
    if value is None:
        return np.ones(n_samples, dtype=bool)
    mask = np.asarray(value, dtype=bool)
    if mask.ndim != 1 or mask.size != n_samples:
        raise AFMAGValidationError(
            "record_mask must have one boolean per navigation sample"
        )
    return mask


def _optional_sample_array(
    value: Any | None,
    *,
    name: str,
    expected: tuple[int, ...],
    n_samples: int,
) -> np.ndarray | None:
    if value is None:
        return None
    arr = np.asarray(value)
    if n_samples == 1 and arr.shape == expected[1:]:
        arr = arr[None, ...]
    if arr.shape != expected:
        raise AFMAGValidationError(
            f"{name} must have shape {expected}, got {arr.shape}"
        )
    return arr


def build_airmt_line(
    line_id: str,
    navigation: NavigationTrack,
    tensor: Any,
    *,
    frequency: Any,
    record_mask: Any | None = None,
    variance: Any | None = None,
    inverse_signal_covariance: Any | None = None,
    residual_covariance: Any | None = None,
    units: str | None = AFMAG_TENSOR_UNITS,
    include_amplification_parameter: bool = True,
    reference_station: AFMAGReferenceStation | None = None,
    system_spec: AirMtSystemSpec | None = None,
    orientation: OrientationMeta | None = None,
    attrs: Mapping[str, Any] | None = None,
) -> AirborneEMLine:
    """Build one tensor-AFMAG/AirMt flight line from decoded arrays."""
    if not isinstance(navigation, NavigationTrack):
        raise TypeError("navigation must be a NavigationTrack")
    n_samples = navigation.n_samples
    data = np.asarray(tensor)
    if n_samples == 1 and data.shape[-2:] == (3, 2):
        if data.ndim == 2:
            data = data[None, None, ...]
        elif data.ndim == 3:
            data = data[None, ...]
    if data.ndim != 4 or data.shape[0] != n_samples:
        raise AFMAGValidationError(
            "AirMt line tensor must have shape (samples, nf, 3, 2)"
        )
    if data.shape[2:] != (3, 2):
        raise AFMAGValidationError(
            "AirMt line tensor must have matrix shape (3, 2)"
        )
    if data.dtype.kind not in "biufc":
        raise AFMAGValidationError("AirMt line tensor must be numeric")

    n_frequency = data.shape[1]
    common_freq, freq_rows = _line_frequency_rows(
        frequency,
        n_samples=n_samples,
        n_frequency=n_frequency,
    )
    mask = _record_mask(record_mask, n_samples=n_samples)

    var = _optional_sample_array(
        variance,
        name="variance",
        expected=(n_samples, n_frequency, 3, 2),
        n_samples=n_samples,
    )
    inv = _optional_sample_array(
        inverse_signal_covariance,
        name="inverse_signal_covariance",
        expected=(n_samples, n_frequency, 2, 2),
        n_samples=n_samples,
    )
    res = _optional_sample_array(
        residual_covariance,
        name="residual_covariance",
        expected=(n_samples, n_frequency, 3, 3),
        n_samples=n_samples,
    )

    line = AirborneEMLine(
        line_id=str(line_id),
        navigation=navigation,
        attrs={"technology": "AirMt", **dict(attrs or {})},
    )
    for index, sample_id in enumerate(navigation.sample_ids):
        if not mask[index]:
            continue
        freq_row = common_freq if common_freq is not None else freq_rows[index]
        record = build_airmt_record(
            sample_id,
            data[index],
            frequency=freq_row,
            variance=None if var is None else var[index],
            inverse_signal_covariance=None if inv is None else inv[index],
            residual_covariance=None if res is None else res[index],
            units=units,
            include_amplification_parameter=include_amplification_parameter,
            reference_station=reference_station,
            system_spec=system_spec,
            orientation=orientation,
            product_id=f"{line_id}.{sample_id}",
            record_attrs={"navigation_index": index},
        )
        line.add_record(record)
    return line


def build_original_afmag_line(
    line_id: str,
    navigation: NavigationTrack,
    tilt: Any,
    *,
    frequency: Any,
    record_mask: Any | None = None,
    response_kind: str = "tilt_angle",
    units: str | None = None,
    variance: Any | None = None,
    system_spec: OriginalAFMAGSystemSpec | None = None,
    orientation: OrientationMeta | None = None,
    attrs: Mapping[str, Any] | None = None,
) -> AirborneEMLine:
    """Build one historical comparator-AFMAG flight line."""
    if not isinstance(navigation, NavigationTrack):
        raise TypeError("navigation must be a NavigationTrack")
    n_samples = navigation.n_samples
    data = np.asarray(tilt)
    if n_samples == 1:
        if data.ndim == 0:
            data = data.reshape(1, 1)
        elif data.ndim == 1:
            data = data[None, ...]
    if data.ndim != 2 or data.shape[0] != n_samples:
        raise AFMAGValidationError(
            "original AFMAG line tilt must have shape (samples, nf)"
        )
    if data.dtype.kind not in "biufc":
        raise AFMAGValidationError("original AFMAG line tilt must be numeric")

    n_frequency = data.shape[1]
    common_freq, freq_rows = _line_frequency_rows(
        frequency,
        n_samples=n_samples,
        n_frequency=n_frequency,
    )
    mask = _record_mask(record_mask, n_samples=n_samples)

    var = None
    if variance is not None:
        var = np.asarray(variance)
        if n_samples == 1 and var.ndim == 1:
            var = var[None, ...]
        if var.shape != (n_samples, n_frequency):
            raise AFMAGValidationError(
                "tilt variance must have shape (samples, nf)"
            )

    line = AirborneEMLine(
        line_id=str(line_id),
        navigation=navigation,
        attrs={"technology": "AFMAG", **dict(attrs or {})},
    )
    for index, sample_id in enumerate(navigation.sample_ids):
        if not mask[index]:
            continue
        freq_row = common_freq if common_freq is not None else freq_rows[index]
        record = build_original_afmag_record(
            sample_id,
            data[index],
            frequency=freq_row,
            response_kind=response_kind,
            units=units,
            variance=None if var is None else var[index],
            system_spec=system_spec,
            orientation=orientation,
            product_id=f"{line_id}.{sample_id}",
            record_attrs={"navigation_index": index},
        )
        line.add_record(record)
    return line


def build_airmt_dataset(
    name: str,
    lines: Any,
    *,
    survey: SurveyMeta | None = None,
    system_spec: AirMtSystemSpec | None = None,
    instrument_serial: str | None = None,
    attrs: Mapping[str, Any] | None = None,
) -> AirborneEMDataset:
    """Build a common airborne dataset for tensor AFMAG/AirMt lines."""
    spec = system_spec or AirMtSystemSpec()
    if not isinstance(spec, AirMtSystemSpec):
        raise TypeError("system_spec must be an AirMtSystemSpec or None")
    line_map = {
        line.line_id: line
        for line in lines
    }
    meta = survey or SurveyMeta(name=str(name), method="AEM")
    return AirborneEMDataset(
        name=str(name),
        lines=line_map,
        survey=meta,
        instrument=spec.to_instrument_meta(serial=instrument_serial),
        attrs={"technology": "AirMt", **dict(attrs or {})},
    )


def build_original_afmag_dataset(
    name: str,
    lines: Any,
    *,
    survey: SurveyMeta | None = None,
    system_spec: OriginalAFMAGSystemSpec | None = None,
    instrument_serial: str | None = None,
    attrs: Mapping[str, Any] | None = None,
) -> AirborneEMDataset:
    """Build a common airborne dataset for original AFMAG lines."""
    spec = system_spec or OriginalAFMAGSystemSpec()
    if not isinstance(spec, OriginalAFMAGSystemSpec):
        raise TypeError(
            "system_spec must be OriginalAFMAGSystemSpec or None"
        )
    line_map = {
        line.line_id: line
        for line in lines
    }
    meta = survey or SurveyMeta(name=str(name), method="AEM")
    return AirborneEMDataset(
        name=str(name),
        lines=line_map,
        survey=meta,
        instrument=spec.to_instrument_meta(serial=instrument_serial),
        attrs={"technology": "AFMAG", **dict(attrs or {})},
    )
