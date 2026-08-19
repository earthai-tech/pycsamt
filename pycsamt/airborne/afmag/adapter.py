# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

"""Scientific adapter contracts for historical and tensor AFMAG products.

No native vendor/archive parser is defined here.  The module accepts decoded
arrays and maps two distinct AFMAG generations into the common pyCSAMT model:

* original comparator AFMAG -> scalar line-direction tilt response;
* tensor AFMAG / AirMt -> 3 x 2 interstation magnetic transfer function plus
  the optional rotationally invariant amplification parameter.

Shape/frequency/mask/metadata normalization that is not specific to
either AFMAG generation is delegated to
:mod:`pycsamt.airborne.validation` (shared with
:mod:`pycsamt.airborne.mobilemt` and :mod:`pycsamt.airborne.ztem`) via
its ``error_cls`` parameter, so every ``AFMAGValidationError`` raised
here still comes from this module even though the check itself is not
reimplemented per technology. ``build_original_afmag_*`` and
``build_airmt_*`` remain genuinely separate below them, because the
two AFMAG generations are scientifically distinct responses (a real
scalar tilt versus a complex 3x2 interstation tensor), not two
configurations of one response.
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
    SiteMeta,
    SurveyMeta,
)
from ..base import AirborneEMDataset, AirborneEMLine, AirborneEMRecord
from ..navigation import NavigationTrack
from ..validation import (
    merge_remote_reference_processing,
    normalize_estimate_array,
    normalize_frequency,
    normalize_record_mask,
    normalize_sample_axis_array,
    reference_station_mapping,
    resolve_frequency_or_periods,
)
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


def _normalize_tilt(value: Any) -> np.ndarray:
    """Return one sample-batched real-valued ``(nf,)`` tilt array.

    Original comparator AFMAG reports a real angle or deflection, not
    a complex EM transfer function; a genuinely complex input (with
    non-zero imaginary part) is rejected rather than silently
    truncated to its real part.

    Raises
    ------
    AFMAGValidationError
        If *value* is not scalar/1-D, is not numeric, or is complex
        with a non-zero imaginary component.
    """
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


def _normalize_scalar_variance(
    value: Any,
    *,
    n_frequency: int,
) -> np.ndarray:
    """Return one sample-batched real-valued ``(nf, 1, 1)`` variance array.

    Distinct from :func:`~pycsamt.airborne.validation.normalize_estimate_array`
    in two ways beyond shape: a bare ``(nf,)`` vector is also accepted
    (there is no separate output/input channel axis to distinguish it
    from), and the result is coerced real, matching
    :func:`_normalize_tilt`'s real-valued response.

    Raises
    ------
    AFMAGValidationError
        If *value* does not match one of the accepted shapes, or is
        complex with a non-zero imaginary component.
    """
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
    """Validate and return an AirMt 3 x 2 interstation magnetic TF.

    Parameters
    ----------
    tf : TransferFunction
        Transfer function to validate in place.

    Returns
    -------
    TransferFunction
        *tf*, unchanged, for convenient chaining after
        :meth:`~pycsamt.emtf.EMTF.add_transfer_function`.

    Raises
    ------
    TypeError
        If *tf* is not a :class:`~pycsamt.emtf.TransferFunction`.
    AFMAGValidationError
        If *tf* does not use the standard EMTF interstation
        ``TI``/``interstation_transfer_functions`` datatype with
        ``Hx``/``Hy`` inputs, ``Hx``/``Hy``/``Hz`` outputs, and matrix
        shape ``(3, 2)``.
    """
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
    """Validate and return an original AFMAG scalar tilt response.

    Parameters
    ----------
    tf : TransferFunction
        Transfer function to validate in place.

    Returns
    -------
    TransferFunction
        *tf*, unchanged, for convenient chaining after
        :meth:`~pycsamt.emtf.EMTF.add_transfer_function`.

    Raises
    ------
    TypeError
        If *tf* is not a :class:`~pycsamt.emtf.TransferFunction`.
    AFMAGValidationError
        If *tf* does not use the ``afmag_tilt`` datatype, is not
        scalar (no input/output channels, matrix shape ``(1, 1)``),
        or carries a non-zero imaginary component.
    """
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


def _airmt_notes(
    spec: AirMtSystemSpec,
    reference_station: AFMAGReferenceStation | None,
) -> dict[str, Any]:
    """Return a ``{"AFMAG": {...}}`` block for ``EMTF.metadata["notes"]``.

    Presentation metadata only, mirroring
    :func:`~pycsamt.airborne.ztem.adapter._xml_notes_mapping`: it does
    not feed :attr:`~pycsamt.emtf.EMTF.processing` (see
    :func:`_processing_for_reference` for that) and is not read back
    by any adapter.
    """
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
    """Return a ``{"AFMAG": {...}}`` block for ``EMTF.metadata["notes"]``.

    Presentation metadata only; see :func:`_airmt_notes`. The original
    generation has no reference-station geometry to record, unlike
    AirMt, so this block is purely descriptive of the historical
    instrument.
    """
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
    """Fold *reference_station* into *processing*'s remote reference.

    Thin AirMt-specific wrapper around
    :func:`~pycsamt.airborne.validation.merge_remote_reference_processing`;
    see that function's docstring for the merge/conflict semantics
    (shared verbatim with
    :func:`~pycsamt.airborne.ztem.adapter._processing_for_reference`).
    """
    extra = {}
    if reference_station is not None:
        extra = {
            "measured_channels": list(reference_station.measured_channels),
            "transfer_input_channels": list(
                reference_station.transfer_input_channels
            ),
        }
    return merge_remote_reference_processing(
        reference_station,
        processing,
        reference_type="fixed_ground_magnetic",
        technology="AirMt",
        extra=extra,
        error_cls=AFMAGValidationError,
    )


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
    """Build one sample-level tensor AFMAG/AirMt :class:`EMTF` response.

    Parameters
    ----------
    tensor : array-like
        Complex interstation magnetic transfer function with shape
        ``(nf, 3, 2)``, or one ``(3, 2)`` matrix for a single
        frequency.
    frequency, periods : array-like, optional
        Exactly one positive frequency or period vector must be
        supplied.
    variance : array-like, optional
        Component variance with shape ``(nf, 3, 2)``.
    inverse_signal_covariance : array-like, optional
        Input covariance factor ``S`` with shape ``(nf, 2, 2)``.
    residual_covariance : array-like, optional
        Output residual covariance ``N`` with shape ``(nf, 3, 3)``.
    include_amplification_parameter : bool, default True
        Whether to also attach the derived
        ``airmt_amplification_parameter`` transfer function; see
        :func:`compute_airmt_amplification_parameter`.
    amplification_zero_policy : {"nan", "raise"}, default "nan"
        Forwarded to :func:`compute_airmt_amplification_parameter`.
    reference_station : AFMAGReferenceStation, optional
        Fixed ground reference metadata; folded into
        ``EMTF.processing`` and ``EMTF.attrs["afmag"]`` (see
        :func:`_processing_for_reference`).
    system_spec : AirMtSystemSpec, optional
        System-description metadata; defaults to published nominal
        values.

    Returns
    -------
    EMTF
        Document carrying the ``interstation_transfer_functions``
        (``TI``) transfer function and, unless disabled, the derived
        amplification parameter.

    Raises
    ------
    AFMAGValidationError
        If *tensor*, *frequency*/*periods*, or any statistical
        estimate does not match its expected shape.
    TypeError
        If *system_spec*, *reference_station*, *site*, or
        *orientation* has the wrong type.
    """
    register_afmag_datatypes()
    freq, period_arr = resolve_frequency_or_periods(
        frequency=frequency,
        periods=periods,
        error_cls=AFMAGValidationError,
    )
    data = normalize_estimate_array(
        tensor,
        n_frequency=freq.size,
        tail=(3, 2),
        name="tensor",
        error_cls=AFMAGValidationError,
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
                data=normalize_estimate_array(
                    variance,
                    n_frequency=freq.size,
                    tail=(3, 2),
                    name="variance",
                    error_cls=AFMAGValidationError,
                ),
            )
        )
    if inverse_signal_covariance is not None:
        tf.add_estimate(
            StatisticalEstimate(
                name="INVSIGCOV",
                kind="inverse_signal_covariance",
                data=normalize_estimate_array(
                    inverse_signal_covariance,
                    n_frequency=freq.size,
                    tail=(2, 2),
                    name="inverse_signal_covariance",
                    error_cls=AFMAGValidationError,
                ),
            )
        )
    if residual_covariance is not None:
        tf.add_estimate(
            StatisticalEstimate(
                name="RESIDCOV",
                kind="residual_covariance",
                data=normalize_estimate_array(
                    residual_covariance,
                    n_frequency=freq.size,
                    tail=(3, 3),
                    name="residual_covariance",
                    error_cls=AFMAGValidationError,
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
        "reference_station": reference_station_mapping(
            reference_station,
            channel_fields=(
                "measured_channels",
                "transfer_input_channels",
            ),
        ),
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

    Parameters
    ----------
    tilt : array-like
        Real-valued tilt/deflection response, scalar or shape ``(nf,)``.
    frequency, periods : array-like, optional
        Exactly one positive frequency or period vector must be
        supplied.
    response_kind : str, default "tilt_angle"
        Either ``"tilt_angle"`` (a calibrated angle, defaulting to
        degrees) or ``"comparator_deflection"`` (a raw archival
        instrument deflection with no implied angular unit).
    variance : array-like, optional
        Real-valued response variance, scalar or shape ``(nf,)`` /
        ``(nf, 1, 1)``.
    system_spec : OriginalAFMAGSystemSpec, optional
        System-description metadata; defaults to published nominal
        values.

    Returns
    -------
    EMTF
        Document carrying the ``afmag_tilt`` transfer function.

    Raises
    ------
    AFMAGValidationError
        If *tilt*, *frequency*/*periods*, or *variance* does not
        match its expected shape.
    ValueError
        If *response_kind* is not a recognized value.
    TypeError
        If *system_spec*, *site*, or *orientation* has the wrong type.
    """
    register_afmag_datatypes()
    freq, period_arr = resolve_frequency_or_periods(
        frequency=frequency,
        periods=periods,
        error_cls=AFMAGValidationError,
    )
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
    """Build one airborne record from decoded tensor AFMAG/AirMt data.

    Parameters
    ----------
    sample_id : str
        Navigation sample identifier for the new record.
    tensor : array-like
        Forwarded to :func:`build_airmt_emtf`.
    frequency, periods : array-like, optional
        Exactly one must be supplied; forwarded to
        :func:`build_airmt_emtf`.
    fields, quality, record_attrs : dict, optional
        Forwarded to :class:`~pycsamt.airborne.base.AirborneEMRecord`.
    **emtf_kwargs
        Forwarded to :func:`build_airmt_emtf`.

    Returns
    -------
    AirborneEMRecord
        The record, with its EMTF ``product_id`` defaulted to
        ``str(sample_id)`` unless overridden in ``emtf_kwargs``.
    """
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
    """Build one record from a historical AFMAG tilt response.

    Parameters
    ----------
    sample_id : str
        Navigation sample identifier for the new record.
    tilt : array-like
        Forwarded to :func:`build_original_afmag_emtf`.
    frequency, periods : array-like, optional
        Exactly one must be supplied; forwarded to
        :func:`build_original_afmag_emtf`.
    fields, quality, record_attrs : dict, optional
        Forwarded to :class:`~pycsamt.airborne.base.AirborneEMRecord`.
    **emtf_kwargs
        Forwarded to :func:`build_original_afmag_emtf`.

    Returns
    -------
    AirborneEMRecord
        The record, with its EMTF ``product_id`` defaulted to
        ``str(sample_id)`` unless overridden in ``emtf_kwargs``.
    """
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
    """Return ``(common_frequency, frequency_rows)`` for one flight line.

    Deliberately not
    :func:`~pycsamt.airborne.validation.resolve_line_frequency_grid`:
    that shared helper only shape-checks a per-sample ``(n_samples,
    nf)`` grid, deferring each row's own finiteness/positivity to a
    later per-sample call so a row belonging to a ``record_mask``-
    excluded sample never needs to be valid. AFMAG validates the
    entire per-sample grid eagerly instead, catching a bad frequency
    row at line-construction time even for samples that end up masked
    out. That is a real, intentional behavioral difference from
    MobileMT/ZTEM, not an oversight -- do not "fix" it into matching
    the shared helper without confirming the stricter behavior is no
    longer wanted.
    """
    freq = np.asarray(frequency, dtype=float)
    if freq.ndim == 1:
        common = normalize_frequency(freq, error_cls=AFMAGValidationError)
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


def _optional_sample_array(
    value: Any | None,
    *,
    name: str,
    expected: tuple[int, ...],
    n_samples: int,
) -> np.ndarray | None:
    """Return one line's optional per-sample estimate array, or ``None``.

    Thin ``None``-tolerant wrapper around
    :func:`~pycsamt.airborne.validation.normalize_sample_axis_array`.
    """
    if value is None:
        return None
    return normalize_sample_axis_array(
        value,
        name=name,
        n_samples=n_samples,
        expected=expected,
        error_cls=AFMAGValidationError,
    )


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
    """Build one tensor-AFMAG/AirMt flight line from decoded arrays.

    Parameters
    ----------
    line_id : str
        Flight-line identifier.
    navigation : NavigationTrack
        Sample-aligned navigation defining the line's sample axis.
    tensor : array-like
        Line-batched interstation tensor. For ``n_samples == 1``, a
        single ``(3, 2)`` matrix or one ``(nf, 3, 2)`` stack is also
        accepted (``n_frequency`` is not known in advance, unlike
        :func:`~pycsamt.airborne.validation.normalize_sample_axis_array`'s
        contract); for ``n_samples > 1``, the canonical
        ``(samples, nf, 3, 2)`` shape is required directly.
    frequency : array-like
        Either one shared ``(nf,)`` vector or a per-sample
        ``(n_samples, nf)`` grid; see :func:`_line_frequency_rows`.
    record_mask : array-like of bool, optional
        Marks which navigation samples get an attached record;
        ``None`` means every sample does.
    variance, inverse_signal_covariance, residual_covariance :
    array-like, optional
        Per-sample statistical estimates, each shaped
        ``(n_samples, nf, *tail)``; forwarded per sample to
        :func:`build_airmt_record`.
    include_amplification_parameter : bool, default True
        Forwarded to :func:`build_airmt_emtf` for every sample.
    units, reference_station, system_spec, orientation : optional
        Forwarded to :func:`build_airmt_emtf` for every sample.
    attrs : dict, optional
        Line-level extension metadata; ``"technology"`` is set to
        ``"AirMt"`` here.

    Returns
    -------
    AirborneEMLine
        The line, with one record per sample where ``record_mask``
        (or its default) is ``True``.

    Raises
    ------
    TypeError
        If *navigation* is not a :class:`NavigationTrack`.
    AFMAGValidationError
        If *tensor*, *frequency*, *record_mask*, or any statistical
        estimate does not match its expected shape.
    """
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
    mask = normalize_record_mask(
        record_mask,
        n_samples=n_samples,
        error_cls=AFMAGValidationError,
    )

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
    """Build one historical comparator-AFMAG flight line.

    Parameters
    ----------
    line_id : str
        Flight-line identifier.
    navigation : NavigationTrack
        Sample-aligned navigation defining the line's sample axis.
    tilt : array-like
        Line-batched real tilt/deflection response, shape
        ``(n_samples, nf)`` (or unbatched when ``n_samples == 1``).
    frequency : array-like
        Either one shared ``(nf,)`` vector or a per-sample
        ``(n_samples, nf)`` grid; see :func:`_line_frequency_rows`.
    record_mask : array-like of bool, optional
        Marks which navigation samples get an attached record;
        ``None`` means every sample does.
    response_kind, units : optional
        Forwarded to :func:`build_original_afmag_emtf` for every
        sample.
    variance : array-like, optional
        Per-sample response variance, shape ``(n_samples, nf)``.
    system_spec, orientation : optional
        Forwarded to :func:`build_original_afmag_emtf` for every
        sample.
    attrs : dict, optional
        Line-level extension metadata; ``"technology"`` is set to
        ``"AFMAG"`` here.

    Returns
    -------
    AirborneEMLine
        The line, with one record per sample where ``record_mask``
        (or its default) is ``True``.

    Raises
    ------
    TypeError
        If *navigation* is not a :class:`NavigationTrack`.
    AFMAGValidationError
        If *tilt*, *frequency*, *record_mask*, or *variance* does not
        match its expected shape.
    """
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
    mask = normalize_record_mask(
        record_mask,
        n_samples=n_samples,
        error_cls=AFMAGValidationError,
    )

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
    """Build a common airborne dataset for tensor AFMAG/AirMt lines.

    Parameters
    ----------
    name : str
        Dataset/survey name.
    lines : iterable of AirborneEMLine
        Lines to attach, typically previously built by
        :func:`build_airmt_line`.
    survey : SurveyMeta, optional
        Survey-level metadata; defaults to
        ``SurveyMeta(name=name, method="AEM")``.
    system_spec : AirMtSystemSpec, optional
        Used to build the dataset's ``instrument`` metadata; defaults
        to published nominal values.
    instrument_serial : str, optional
        Forwarded to :meth:`AirMtSystemSpec.to_instrument_meta`.
    attrs : dict, optional
        Dataset-level extension metadata; ``"technology"`` is set to
        ``"AirMt"`` here.

    Returns
    -------
    AirborneEMDataset
        The dataset, with every line attached.

    Raises
    ------
    TypeError
        If *system_spec* has the wrong type, or an entry of *lines*
        is not an :class:`~pycsamt.airborne.base.AirborneEMLine` (this
        is enforced by :class:`AirborneEMDataset` construction, unlike
        :func:`~pycsamt.airborne.ztem.build_ztem_dataset`, which
        checks explicitly beforehand and additionally rejects a line
        tagged with a conflicting technology; this function does not
        perform that second, stricter check).
    """
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
    """Build a common airborne dataset for original AFMAG lines.

    Parameters
    ----------
    name : str
        Dataset/survey name.
    lines : iterable of AirborneEMLine
        Lines to attach, typically previously built by
        :func:`build_original_afmag_line`.
    survey : SurveyMeta, optional
        Survey-level metadata; defaults to
        ``SurveyMeta(name=name, method="AEM")``.
    system_spec : OriginalAFMAGSystemSpec, optional
        Used to build the dataset's ``instrument`` metadata; defaults
        to published nominal values.
    instrument_serial : str, optional
        Forwarded to
        :meth:`OriginalAFMAGSystemSpec.to_instrument_meta`.
    attrs : dict, optional
        Dataset-level extension metadata; ``"technology"`` is set to
        ``"AFMAG"`` here.

    Returns
    -------
    AirborneEMDataset
        The dataset, with every line attached.

    Raises
    ------
    TypeError
        If *system_spec* has the wrong type, or an entry of *lines* is
        not an :class:`~pycsamt.airborne.base.AirborneEMLine`; see
        :func:`build_airmt_dataset` for how this differs from the
        ZTEM/MobileMT dataset builders.
    """
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
