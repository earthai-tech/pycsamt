# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

"""Array-to-scientific-object adapter contract for ZTEM products.

ZTEM is represented as the standard vertical magnetic transfer function
(tipper), while acquisition geometry records that the airborne Hz output is
referenced to fixed ground Hx/Hy measurements. No proprietary native file
schema is assumed here.

Shape/frequency/mask normalization that is not specific to the ZTEM
tipper layout is delegated to :mod:`pycsamt.airborne.validation`
(shared with :mod:`pycsamt.airborne.mobilemt` and, once it moves onto
the same helpers, :mod:`pycsamt.airborne.afmag`) via its ``error_cls``
parameter, so every ``ZTEMValidationError`` raised here still comes
from this module even though the check itself is not reimplemented
per technology.
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
from ..validation import (
    normalize_estimate_array,
    normalize_frequency,
    normalize_record_mask,
    normalize_sample_axis_array,
    resolve_frequency_or_periods,
    resolve_line_frequency_grid,
)
from .base import ZTEMReferenceStation, ZTEMSystemSpec
from .constants import (
    ZTEM_INPUT_CHANNELS,
    ZTEM_OUTPUT_CHANNELS,
    ZTEM_TIPPER_COMPONENTS,
    ZTEM_TIPPER_TAG,
    ZTEM_TIPPER_UNITS,
)

__all__ = [
    "ZTEMValidationError",
    "validate_ztem_transfer_function",
    "build_ztem_emtf",
    "build_ztem_record",
    "build_ztem_line",
    "build_ztem_dataset",
]


class ZTEMValidationError(ValueError):
    """Raised when decoded ZTEM scientific arrays are inconsistent."""


def _normalize_tipper(value: Any) -> np.ndarray:
    """Return one sample-batched ``(nf, 1, 2)`` complex tipper array.

    Accepts, and promotes to that canonical shape, a single
    ``(Tzx, Tzy)`` pair of shape ``(2,)``, a stack of shape
    ``(nf, 2)``, or the canonical EMTF matrix shape ``(nf, 1, 2)``
    already.

    Raises
    ------
    ZTEMValidationError
        If *value* matches none of the accepted shapes, or is not
        numeric.
    """
    arr = np.asarray(value)
    if arr.ndim == 1 and arr.shape == (2,):
        arr = arr[None, None, :]
    elif arr.ndim == 2 and arr.shape[1:] == (2,):
        arr = arr[:, None, :]
    elif arr.ndim != 3 or arr.shape[1:] != (1, 2):
        raise ZTEMValidationError(
            "tipper must have shape (2,), (nf, 2), or (nf, 1, 2)"
        )
    if arr.dtype.kind not in "biufc":
        raise ZTEMValidationError("tipper must be numeric")
    return arr


def _reference_mapping(
    reference_station: ZTEMReferenceStation | None,
) -> dict[str, Any] | None:
    """Return *reference_station* as a plain ``EMTF.attrs`` mapping."""
    if reference_station is None:
        return None
    result: dict[str, Any] = {
        "station_id": reference_station.preferred_id,
        "magnetic_channels": list(reference_station.magnetic_channels),
    }
    if reference_station.site is not None:
        result["site"] = reference_station.site.to_dict(max_depth=2)
    if reference_station.attrs:
        result["attrs"] = dict(reference_station.attrs)
    return result


def _xml_notes_mapping(
    spec: ZTEMSystemSpec,
    reference_station: ZTEMReferenceStation | None,
) -> dict[str, Any]:
    """Return a ``{"ZTEM": {...}}`` block for ``EMTF.metadata["notes"]``.

    This is presentation metadata only: it does not feed
    :attr:`~pycsamt.emtf.EMTF.processing` (see
    :func:`_processing_for_reference` for that) and is not read back
    by any adapter, so its keys are free to describe the system
    however is most useful for a human or archive reading the
    serialized EMTF XML.
    """
    low, high = spec.practical_frequency_range_hz
    ztem: dict[str, Any] = {
        "PracticalFrequencyMinHz": low,
        "PracticalFrequencyMaxHz": high,
        "TypicalFrequencyCountMin": spec.typical_frequency_count[0],
        "TypicalFrequencyCountMax": spec.typical_frequency_count[1],
        "TimeSeriesSamplingRateHz": spec.time_series_sampling_rate_hz,
        "NominalOutputRateHz": spec.nominal_output_rate_hz,
        "InputChannels": ",".join(spec.input_channels),
        "OutputChannels": ",".join(spec.output_channels),
        "TransferFunction": "Hz = Tzx*Hx + Tzy*Hy",
    }
    if reference_station is not None:
        if reference_station.preferred_id is not None:
            ztem["ReferenceStationId"] = reference_station.preferred_id
        site = reference_station.site
        location = site.location if site is not None else None
        if location is not None:
            if location.latitude is not None:
                ztem["ReferenceLatitude"] = location.latitude
            if location.longitude is not None:
                ztem["ReferenceLongitude"] = location.longitude
            if location.elevation is not None:
                ztem["ReferenceElevation"] = location.elevation
            if location.datum is not None:
                ztem["ReferenceDatum"] = location.datum
    return {"ZTEM": ztem}


def _processing_for_reference(
    reference_station: ZTEMReferenceStation | None,
    processing: ProcessingMeta | None,
) -> ProcessingMeta | None:
    """Fold *reference_station* into *processing*'s remote reference.

    ``reference_station`` is the ZTEM-specific, scientifically typed
    way to supply the fixed ground magnetic reference; a caller-
    supplied *processing* is the general EMTF way. When both are
    given, they must describe the same reference site -- this raises
    :class:`ZTEMValidationError` on conflict rather than silently
    preferring one over the other. When only *reference_station* is
    given, this synthesizes a :class:`~pycsamt.metadata.ProcessingMeta`
    around it so :attr:`~pycsamt.emtf.EMTF.processing` is always the
    one place downstream code looks, in
    :func:`~pycsamt.airborne.qc.assess_airborne_qc` and elsewhere.
    """
    if processing is not None and not isinstance(processing, ProcessingMeta):
        raise TypeError("processing must be a ProcessingMeta or None")
    if reference_station is None:
        return processing

    remote = RemoteReferenceMeta(
        reference_type="fixed_ground_horizontal_magnetic",
        site=reference_station.preferred_id,
        extra={
            "channels": list(reference_station.magnetic_channels),
            "technology": "ZTEM",
        },
    )
    if processing is None:
        return ProcessingMeta(remote_reference=remote)
    if processing.remote_reference is None:
        processing = ProcessingMeta(
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
        return processing

    existing = processing.remote_reference
    reference_id = reference_station.preferred_id
    if existing.site is not None and reference_id is not None:
        if str(existing.site) != str(reference_id):
            raise ZTEMValidationError(
                "processing remote-reference site conflicts with "
                "ZTEMReferenceStation"
            )
    return processing


def validate_ztem_transfer_function(
    tf: TransferFunction,
) -> TransferFunction:
    """Validate and return a ZTEM 1x2 vertical magnetic tipper.

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
    ZTEMValidationError
        If *tf* does not use the standard EMTF ``tipper`` datatype
        with ``Hx``/``Hy`` inputs, ``Hz`` output, and matrix shape
        ``(1, 2)``.
    """
    if not isinstance(tf, TransferFunction):
        raise TypeError("tf must be a TransferFunction")
    if tf.name != ZTEM_TIPPER_TAG:
        raise ZTEMValidationError(
            "ZTEM must use the standard EMTF 'tipper' datatype"
        )
    if tuple(tf.input_channels) != ZTEM_INPUT_CHANNELS:
        raise ZTEMValidationError("ZTEM tipper inputs must be ('Hx', 'Hy')")
    if tuple(tf.output_channels) != ZTEM_OUTPUT_CHANNELS:
        raise ZTEMValidationError("ZTEM tipper output must be ('Hz',)")
    if tf.data.shape[1:] != (1, 2):
        raise ZTEMValidationError(
            "ZTEM tipper must have matrix shape (1, 2)"
        )
    return tf


def build_ztem_emtf(
    tipper: Any,
    *,
    frequency: Any | None = None,
    periods: Any | None = None,
    units: str | None = ZTEM_TIPPER_UNITS,
    variance: Any | None = None,
    inverse_signal_covariance: Any | None = None,
    residual_covariance: Any | None = None,
    product_id: str | None = None,
    description: str | None = None,
    reference_station: ZTEMReferenceStation | None = None,
    system_spec: ZTEMSystemSpec | None = None,
    site: SiteMeta | None = None,
    orientation: OrientationMeta | None = None,
    processing: ProcessingMeta | None = None,
    attrs: Mapping[str, Any] | None = None,
) -> EMTF:
    """Build one sample-level :class:`EMTF` ZTEM response.

    Parameters
    ----------
    tipper : array-like
        Complex ``(Tzx, Tzy)`` values with shape ``(nf, 2)`` or canonical
        EMTF shape ``(nf, 1, 2)``. A single ``(2,)`` vector is accepted.
    frequency, periods : array-like, optional
        Exactly one positive frequency or period vector must be supplied.
    variance : array-like, optional
        Component variance with shape ``(nf, 1, 2)``.
    inverse_signal_covariance : array-like, optional
        Input covariance factor ``S`` with shape ``(nf, 2, 2)``.
    residual_covariance : array-like, optional
        Output residual covariance ``N`` with shape ``(nf, 1, 1)``.

    Notes
    -----
    ZTEM reuses the standard EMTF tipper datatype. The technology-specific
    distinction is acquisition geometry: airborne ``Hz`` is related to fixed
    ground-reference ``Hx`` and ``Hy``. No line-axis orientation is inferred.
    """
    freq, period_arr = resolve_frequency_or_periods(
        frequency=frequency,
        periods=periods,
        error_cls=ZTEMValidationError,
    )
    data = _normalize_tipper(tipper)
    if data.shape[0] != freq.size:
        raise ZTEMValidationError(
            "frequency count does not match tipper: "
            f"{freq.size} != {data.shape[0]}"
        )

    tf = TransferFunction(
        name=ZTEM_TIPPER_TAG,
        data=data,
        input_channels=ZTEM_INPUT_CHANNELS,
        output_channels=ZTEM_OUTPUT_CHANNELS,
        units=units,
        periods=period_arr,
        attrs={
            "technology": "ZTEM",
            "reference_geometry": (
                "airborne_Hz_to_fixed_ground_Hx_Hy"
            ),
            "components": list(ZTEM_TIPPER_COMPONENTS),
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
                    tail=(1, 2),
                    name="variance",
                    error_cls=ZTEMValidationError,
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
                    error_cls=ZTEMValidationError,
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
                    tail=(1, 1),
                    name="residual_covariance",
                    error_cls=ZTEMValidationError,
                ),
            )
        )

    spec = system_spec or ZTEMSystemSpec()
    if not isinstance(spec, ZTEMSystemSpec):
        raise TypeError("system_spec must be a ZTEMSystemSpec or None")
    if reference_station is not None and not isinstance(
        reference_station,
        ZTEMReferenceStation,
    ):
        raise TypeError(
            "reference_station must be ZTEMReferenceStation or None"
        )
    if site is not None and not isinstance(site, SiteMeta):
        raise TypeError("site must be a SiteMeta or None")
    if orientation is not None and not isinstance(
        orientation,
        OrientationMeta,
    ):
        raise TypeError("orientation must be an OrientationMeta or None")

    ztem_meta: dict[str, Any] = {
        "system": spec.to_dict(max_depth=2),
        "reference_station": _reference_mapping(reference_station),
    }
    document_attrs = dict(attrs or {})
    document_attrs["ztem"] = ztem_meta

    doc = EMTF(
        product_id=product_id,
        description=(
            description
            or "ZTEM airborne vertical magnetic transfer-function response"
        ),
        subtype="ztem",
        tags=("ztem", ZTEM_TIPPER_TAG),
        periods=period_arr,
        site=site,
        orientation=orientation,
        processing=_processing_for_reference(reference_station, processing),
        metadata={
            "notes": _xml_notes_mapping(spec, reference_station),
        },
        attrs=document_attrs,
    )
    doc.add_transfer_function(tf)
    validate_ztem_transfer_function(tf)
    return doc


def build_ztem_record(
    sample_id: str,
    tipper: Any,
    *,
    frequency: Any | None = None,
    periods: Any | None = None,
    fields: Mapping[str, Any] | None = None,
    quality: Mapping[str, Any] | None = None,
    record_attrs: Mapping[str, Any] | None = None,
    **emtf_kwargs: Any,
) -> AirborneEMRecord:
    """Build one airborne record from decoded ZTEM tipper values.

    Parameters
    ----------
    sample_id : str
        Navigation sample identifier for the new record.
    tipper : array-like
        Forwarded to :func:`build_ztem_emtf`.
    frequency, periods : array-like, optional
        Exactly one must be supplied; forwarded to
        :func:`build_ztem_emtf`.
    fields, quality, record_attrs : dict, optional
        Forwarded to :class:`~pycsamt.airborne.base.AirborneEMRecord`.
    **emtf_kwargs
        Forwarded to :func:`build_ztem_emtf`.

    Returns
    -------
    AirborneEMRecord
        The record, with its EMTF ``product_id`` defaulted to
        ``str(sample_id)`` unless overridden in ``emtf_kwargs``.
    """
    doc = build_ztem_emtf(
        tipper,
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


def _sample_axis_tipper(
    value: Any,
    *,
    n_samples: int,
) -> np.ndarray:
    """Return one line's tipper data as ``(n_samples, nf, 1, 2)``.

    Combines two independent promotions: the leading sample axis may
    be omitted when ``n_samples == 1`` (see
    :func:`~pycsamt.airborne.validation.normalize_sample_axis_array`),
    and the tipper tail axis may be supplied as ``(..., 2)`` instead
    of the canonical ``(..., 1, 2)`` (see :func:`_normalize_tipper`).
    Both are handled together here because unlike the generic
    per-estimate arrays (variance, covariances), a stacked ``(nf, 2)``
    per-sample tipper is ambiguous between "one sample, ``nf`` rows"
    and "``nf`` samples, one row" without also knowing ``n_samples``.

    Raises
    ------
    ZTEMValidationError
        If *value* does not match one of the accepted shapes for the
        given *n_samples*, or is not numeric.
    """
    arr = np.asarray(value)
    if arr.ndim == 1 and arr.shape == (2,) and n_samples == 1:
        arr = arr[None, None, None, :]
    elif arr.ndim == 2 and arr.shape[-1] == 2 and n_samples == 1:
        arr = arr[None, :, None, :]
    elif (
        arr.ndim == 3
        and arr.shape[-2:] == (1, 2)
        and n_samples == 1
    ):
        arr = arr[None, ...]
    elif arr.ndim == 3 and arr.shape[0] == n_samples:
        if arr.shape[-1] != 2:
            raise ZTEMValidationError(
                "line tipper component axis must contain Tzx/Tzy"
            )
        arr = arr[:, :, None, :]
    if arr.ndim != 4 or arr.shape[0] != n_samples:
        raise ZTEMValidationError(
            "tipper must have shape (samples, nf, 2) or "
            "(samples, nf, 1, 2)"
        )
    if arr.shape[2:] != (1, 2):
        raise ZTEMValidationError(
            "line tipper must have matrix shape (samples, nf, 1, 2)"
        )
    if arr.dtype.kind not in "biufc":
        raise ZTEMValidationError("tipper must be numeric")
    return arr


def build_ztem_line(
    line_id: str,
    navigation: NavigationTrack,
    tipper: Any,
    *,
    frequency: Any,
    record_mask: Any | None = None,
    variance: Any | None = None,
    inverse_signal_covariance: Any | None = None,
    residual_covariance: Any | None = None,
    units: str | None = ZTEM_TIPPER_UNITS,
    reference_station: ZTEMReferenceStation | None = None,
    system_spec: ZTEMSystemSpec | None = None,
    orientation: OrientationMeta | None = None,
    attrs: Mapping[str, Any] | None = None,
) -> AirborneEMLine:
    """Build one ZTEM flight line from decoded sample-aligned arrays.

    Parameters
    ----------
    line_id : str
        Flight-line identifier.
    navigation : NavigationTrack
        Sample-aligned navigation defining the line's sample axis.
    tipper : array-like
        Forwarded to :func:`_sample_axis_tipper`; accepts the
        per-sample analogues of every shape :func:`_normalize_tipper`
        accepts for one sample.
    frequency : array-like
        Either one shared ``(nf,)`` vector or a per-sample
        ``(n_samples, nf)`` grid; see
        :func:`~pycsamt.airborne.validation.resolve_line_frequency_grid`.
    record_mask : array-like of bool, optional
        Marks which navigation samples get an attached record;
        ``None`` means every sample does. Samples excluded here never
        need a valid ``tipper``/``frequency`` row -- navigation points
        are never deleted to represent a rejected EM sample.
    variance, inverse_signal_covariance, residual_covariance :
    array-like, optional
        Per-sample statistical estimates, each shaped
        ``(n_samples, nf, *tail)`` (or unbatched when
        ``n_samples == 1``); forwarded per sample to
        :func:`build_ztem_record`.
    units : str, optional
        Forwarded to :func:`build_ztem_emtf` for every sample.
    reference_station, system_spec, orientation : optional
        Forwarded to :func:`build_ztem_emtf` for every sample.
    attrs : dict, optional
        Line-level extension metadata; ``"technology"`` and
        ``"reference_station"`` are set here unless already present.

    Returns
    -------
    AirborneEMLine
        The line, with one record per sample where ``record_mask``
        (or its default) is ``True``.

    Raises
    ------
    TypeError
        If *navigation* is not a :class:`NavigationTrack`.
    ZTEMValidationError
        If *tipper*, *frequency*, *record_mask*, or any statistical
        estimate does not match its expected shape.
    """
    if not isinstance(navigation, NavigationTrack):
        raise TypeError("navigation must be a NavigationTrack")
    n_samples = navigation.n_samples
    data = _sample_axis_tipper(tipper, n_samples=n_samples)
    n_frequency = data.shape[1]

    common_frequency, frequency_rows = resolve_line_frequency_grid(
        frequency,
        n_samples=n_samples,
        n_frequency=n_frequency,
        error_cls=ZTEMValidationError,
    )
    mask = normalize_record_mask(
        record_mask,
        n_samples=n_samples,
        error_cls=ZTEMValidationError,
    )

    def _optional_sample(
        value: Any | None,
        *,
        name: str,
        tail: tuple[int, int],
    ) -> np.ndarray | None:
        if value is None:
            return None
        expected = (n_samples, n_frequency, *tail)
        return normalize_sample_axis_array(
            value,
            name=name,
            n_samples=n_samples,
            expected=expected,
            error_cls=ZTEMValidationError,
        )

    var_all = _optional_sample(
        variance,
        name="variance",
        tail=(1, 2),
    )
    inv_all = _optional_sample(
        inverse_signal_covariance,
        name="inverse_signal_covariance",
        tail=(2, 2),
    )
    res_all = _optional_sample(
        residual_covariance,
        name="residual_covariance",
        tail=(1, 1),
    )

    line_attrs = dict(attrs or {})
    line_attrs.setdefault("technology", "ZTEM")
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
            else normalize_frequency(
                frequency_rows[index],
                name=f"frequency[{sample_id}]",
                error_cls=ZTEMValidationError,
            )
        )
        record = build_ztem_record(
            sample_id,
            data[index],
            frequency=sample_frequency,
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
            orientation=orientation,
            product_id=f"{line.line_id}.{sample_id}",
            record_attrs={"navigation_index": index},
        )
        line.add_record(record)
    return line


def build_ztem_dataset(
    name: str,
    lines: Any,
    *,
    survey: SurveyMeta | None = None,
    system_spec: ZTEMSystemSpec | None = None,
    instrument_serial: str | None = None,
    software_version: str = "",
    attrs: Mapping[str, Any] | None = None,
) -> AirborneEMDataset:
    """Build a common airborne dataset from constructed ZTEM lines.

    Parameters
    ----------
    name : str
        Dataset/survey name.
    lines : iterable of AirborneEMLine, or mapping of str to AirborneEMLine
        Lines to attach, typically previously built by
        :func:`build_ztem_line`.
    survey : SurveyMeta, optional
        Survey-level metadata.
    system_spec : ZTEMSystemSpec, optional
        Used to build the dataset's ``instrument`` metadata; defaults
        to published nominal values.
    instrument_serial : str, optional
        Forwarded to :meth:`ZTEMSystemSpec.to_instrument_meta`.
    software_version : str, optional
        Forwarded to :meth:`ZTEMSystemSpec.to_instrument_meta`.
    attrs : dict, optional
        Dataset-level extension metadata; ``"technology"`` and
        ``"ztem_system"`` are set here unless already present.

    Returns
    -------
    AirborneEMDataset
        The dataset, with every line attached via
        :meth:`~pycsamt.airborne.base.AirborneEMDataset.add_line`.

    Raises
    ------
    TypeError
        If *survey*, *system_spec*, or an entry of *lines* has the
        wrong type.
    ZTEMValidationError
        If a line is explicitly tagged with a different technology.
    """
    if survey is not None and not isinstance(survey, SurveyMeta):
        raise TypeError("survey must be a SurveyMeta or None")
    spec = system_spec or ZTEMSystemSpec()
    if not isinstance(spec, ZTEMSystemSpec):
        raise TypeError("system_spec must be a ZTEMSystemSpec or None")

    dataset_attrs = dict(attrs or {})
    dataset_attrs.setdefault("technology", "ZTEM")
    dataset_attrs.setdefault("ztem_system", spec.to_dict(max_depth=2))
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
        if technology and technology != "ztem":
            raise ZTEMValidationError(
                f"line {line.line_id!r} is tagged as {technology!r}"
            )
        dataset.add_line(line)
    return dataset
