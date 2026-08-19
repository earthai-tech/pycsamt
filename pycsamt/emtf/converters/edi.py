# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

"""Explicit conversion between SEG EDI and the format-neutral EMTF model.

This module is deliberately an adapter.  It does not make :class:`EDIFile`
the scientific core and it does not move rotation mathematics into parsing.
The bridge preserves the historical pyCSAMT ``Z.z_err`` convention while
representing EDI ``*.VAR`` blocks explicitly as EMTF complex variances.
"""

from __future__ import annotations

from datetime import datetime
from math import atan2, degrees
from os import PathLike
from pathlib import Path
from typing import Any, Iterable
import warnings

import numpy as np

from ...metadata import (
    ChannelMeta,
    LocationMeta,
    OrientationMeta,
    Person,
    ProcessingMeta,
    ProvenanceMeta,
    RemoteReferenceMeta,
    SiteLayout,
    SiteMeta,
    Software,
)
from ..base import (
    IMPEDANCE_INPUT_CHANNELS,
    IMPEDANCE_OUTPUT_CHANNELS,
    LEGACY_STANDARD_ERROR_KIND,
    TIPPER_INPUT_CHANNELS,
    TIPPER_OUTPUT_CHANNELS,
)
from ..document import EMTF
from ..estimates import StatisticalEstimate
from ..transfer import TransferFunction

__all__ = [
    "DataLossWarning",
    "EMTFEDIConversionError",
    "edi_to_emtf",
    "emtf_to_edi",
    "write_edi",
]


class DataLossWarning(UserWarning):
    """Warning emitted when an EDI conversion cannot preserve information."""


class EMTFEDIConversionError(ValueError):
    """Raised when an EDI conversion would be scientifically invalid."""


_VALID_LOSS_POLICIES = frozenset({"warn", "raise", "ignore"})


def _loss(policy: str, message: str) -> None:
    policy = str(policy).strip().lower()
    if policy not in _VALID_LOSS_POLICIES:
        raise ValueError(
            "on_loss must be one of " f"{sorted(_VALID_LOSS_POLICIES)}"
        )
    if policy == "raise":
        raise EMTFEDIConversionError(message)
    if policy == "warn":
        warnings.warn(message, DataLossWarning, stacklevel=3)


def _nonempty(value: Any) -> bool:
    return value not in (None, "", "None")


def _first(*values: Any) -> Any:
    for value in values:
        if _nonempty(value):
            return value
    return None


def _explicit_keys(lines: Iterable[str] | None) -> set[str]:
    keys: set[str] = set()
    for raw in lines or ():
        text = str(raw).strip()
        if "=" not in text:
            continue
        keys.add(text.split("=", 1)[0].strip().lower())
    return keys


def _year_from_date(value: Any) -> int | None:
    if value is None:
        return None
    if isinstance(value, datetime):
        return int(value.year)
    text = str(value).strip()
    if not text:
        return None
    for fmt in (
        "%m/%d/%y",
        "%m/%d/%Y",
        "%Y/%m/%d",
        "%Y-%m-%d",
        "%Y-%m-%dT%H:%M:%S",
    ):
        try:
            return int(datetime.strptime(text[:19], fmt).year)
        except ValueError:
            pass
    for token in text.replace("-", "/").split("/"):
        if token.isdigit() and len(token) == 4:
            year = int(token)
            if 1800 <= year <= 2500:
                return year
    return None


def _as_edi_date(value: Any) -> str | None:
    """Return a conservative EDI date string without fabricating a date."""
    if value is None:
        return None
    if isinstance(value, datetime):
        return value.strftime("%m/%d/%y")
    text = str(value).strip()
    if not text:
        return None
    for fmt in (
        "%Y-%m-%dT%H:%M:%S",
        "%Y-%m-%d",
        "%Y/%m/%d",
        "%m/%d/%Y",
        "%m/%d/%y",
    ):
        try:
            return datetime.strptime(text[:19], fmt).strftime("%m/%d/%y")
        except ValueError:
            pass
    return text


def _elevation_units_to_emtf(value: Any) -> str:
    text = str(value or "m").strip().lower()
    if text in {"m", "meter", "meters", "metre", "metres"}:
        return "meters"
    if text in {"ft", "foot", "feet"}:
        return "feet"
    return text or "meters"


def _elevation_units_to_edi(value: Any) -> str:
    text = str(value or "meters").strip().lower()
    if text in {"m", "meter", "meters", "metre", "metres"}:
        return "m"
    if text in {"ft", "foot", "feet"}:
        return "ft"
    return text or "m"


def _array_or_none(value: Any) -> np.ndarray | None:
    if value is None:
        return None
    arr = np.asarray(value)
    if arr.size == 0:
        return None
    return arr


def _rotation_vector(value: Any, n_freq: int) -> np.ndarray | None:
    arr = _array_or_none(value)
    if arr is None:
        return None
    arr = np.asarray(arr, dtype=float).ravel()
    if arr.size == 1:
        return np.full(n_freq, float(arr[0]))
    if arr.size == n_freq:
        return arr
    return None


def _raw_file_has_block(edi: Any, tag: str) -> bool:
    path = getattr(edi, "path", None)
    if path is None:
        return False
    try:
        text = Path(path).read_text(
            encoding="utf-8-sig", errors="replace"
        ).upper()
    except (OSError, TypeError, ValueError):
        return False
    wanted = f">{str(tag).strip().upper()}"
    return any(line.lstrip().startswith(wanted) for line in text.splitlines())


def _measurement_angle(measurement: Any) -> float | None:
    azm = getattr(measurement, "azm", None)
    if azm is not None:
        try:
            return float(azm)
        except (TypeError, ValueError):
            pass
    x = getattr(measurement, "x", None)
    y = getattr(measurement, "y", None)
    x2 = getattr(measurement, "x2", None)
    y2 = getattr(measurement, "y2", None)
    if None in (x, y, x2, y2):
        return None
    dx = float(x2) - float(x)
    dy = float(y2) - float(y)
    if dx == 0.0 and dy == 0.0:
        return None
    # EDI site-layout convention: +X is 0 degrees and +Y is 90 degrees.
    return float(degrees(atan2(dy, dx)) % 360.0)


def _canonical_channel_name(value: Any) -> str | None:
    if not _nonempty(value):
        return None
    raw = str(value).strip().upper()
    mapping = {
        "HX": "Hx",
        "HY": "Hy",
        "HZ": "Hz",
        "EX": "Ex",
        "EY": "Ey",
        "RHX": "RHx",
        "RHY": "RHy",
    }
    return mapping.get(raw, str(value).strip())


def _channel_from_measurement(
    measurement: Any,
    *,
    fallback_name: str,
    units: str | None,
    reference: str | None = "site",
) -> ChannelMeta:
    name = _canonical_channel_name(
        getattr(measurement, "chtype", None)
    ) or fallback_name
    upper = name.upper()
    field_type = "electric" if upper.startswith("E") else "magnetic"
    extra: dict[str, Any] = {}
    for source, target in (
        ("id", "edi_id"),
        ("acqchan", "acquisition_channel"),
        ("filter", "filter"),
        ("gain", "gain"),
        ("measdate", "measurement_date"),
    ):
        value = getattr(measurement, source, None)
        if _nonempty(value):
            extra[target] = value
    return ChannelMeta(
        name=name,
        field_type=field_type,
        orientation=_measurement_angle(measurement),
        tilt=getattr(measurement, "dip", None),
        x=getattr(measurement, "x", None),
        y=getattr(measurement, "y", None),
        z=getattr(measurement, "z", None),
        x2=getattr(measurement, "x2", None),
        y2=getattr(measurement, "y2", None),
        z2=getattr(measurement, "z2", None),
        units=units,
        reference=reference,
        sensor_id=getattr(measurement, "sensor", None),
        extra=extra,
    )


def _edi_site_layout(edi: Any) -> SiteLayout | None:
    define = edi.get_section("definemeas")
    mtsect = edi.get_section("mtsect")
    if define is None:
        return None

    measurements = list(getattr(define, "all_meas", lambda: [])())
    by_id = {
        str(getattr(item, "id", "")).strip(): item
        for item in measurements
        if _nonempty(getattr(item, "id", None))
    }
    by_type: dict[str, Any] = {}
    for item in measurements:
        chtype = getattr(item, "chtype", None)
        if _nonempty(chtype):
            by_type.setdefault(str(chtype).strip().upper(), item)

    units = getattr(define, "units", None)
    input_channels: list[ChannelMeta] = []
    output_channels: list[ChannelMeta] = []

    def resolve(role: str, name: str) -> Any | None:
        ref = getattr(mtsect, role, None) if mtsect is not None else None
        if _nonempty(ref) and str(ref).strip() in by_id:
            return by_id[str(ref).strip()]
        return by_type.get(name.upper())

    for role, name in (("hx", "Hx"), ("hy", "Hy")):
        item = resolve(role, name)
        if item is not None:
            input_channels.append(
                _channel_from_measurement(
                    item,
                    fallback_name=name,
                    units=units,
                )
            )

    for role, name in (("hz", "Hz"), ("ex", "Ex"), ("ey", "Ey")):
        item = resolve(role, name)
        if item is not None:
            output_channels.append(
                _channel_from_measurement(
                    item,
                    fallback_name=name,
                    units=units,
                )
            )

    if not input_channels and not output_channels:
        return None

    extra: dict[str, Any] = {}
    if mtsect is not None:
        for key in ("rx", "ry"):
            value = getattr(mtsect, key, None)
            if _nonempty(value):
                extra[f"edi_{key}"] = str(value)
    for key in ("reflat", "reflong", "refelev", "reftype"):
        value = getattr(define, key, None)
        if _nonempty(value):
            extra[f"edi_{key}"] = value

    return SiteLayout(
        input_channels=input_channels,
        output_channels=output_channels,
        input_units=units,
        output_units=units,
        input_reference="site",
        output_reference="site",
        extra=extra,
    )


def _edi_site(edi: Any) -> SiteMeta | None:
    head = edi.get_section("head")
    info = edi.get_section("info")
    mtsect = edi.get_section("mtsect")
    if head is None and info is None and mtsect is None:
        return None

    head_keys = _explicit_keys(getattr(head, "edi_header", None))
    info_keys = _explicit_keys(getattr(info, "ediinfo", None))
    source = getattr(info, "Source", None) if info is not None else None

    head_project = getattr(head, "project", None) if head else None
    head_survey = getattr(head, "survey", None) if head else None
    info_project = (
        getattr(source, "project", None) if "project" in info_keys else None
    )
    info_survey = (
        getattr(source, "survey", None) if "survey" in info_keys else None
    )
    site_name = (
        getattr(source, "sitename", None)
        if "sitename" in info_keys
        else None
    )
    if site_name is None and head is not None:
        site_name = _first(
            getattr(head, "loc", None) if "loc" in head_keys else None,
            getattr(head, "prospect", None)
            if "prospect" in head_keys
            else None,
        )

    site_id = _first(
        getattr(head, "dataid", None) if head is not None else None,
        getattr(mtsect, "sectid", None) if mtsect is not None else None,
    )

    location = None
    if head is not None:
        lat = getattr(head, "lat", None) if "lat" in head_keys else None
        lon = getattr(head, "long", None) if "long" in head_keys else None
        elev = getattr(head, "elev", None) if "elev" in head_keys else None
        decl = (
            getattr(head, "declination", None)
            if "declination" in head_keys
            else None
        )
        datum = getattr(head, "datum", None) if "datum" in head_keys else None
        if any(value is not None for value in (lat, lon, elev, decl, datum)):
            location = LocationMeta(
                latitude=lat,
                longitude=lon,
                elevation=elev,
                datum=datum,
                elevation_units=_elevation_units_to_emtf(
                    getattr(head, "units", None)
                ),
                declination=decl,
            )

    start = (
        getattr(head, "acqdate", None)
        if head is not None and "acqdate" in head_keys
        else None
    )
    end = (
        getattr(head, "enddate", None)
        if head is not None and "enddate" in head_keys
        else None
    )
    extra: dict[str, Any] = {}
    if head is not None:
        for key in (
            "state",
            "county",
            "prospect",
            "loc",
            "chainage",
            "coordsys",
            "stdvers",
        ):
            value = getattr(head, key, None)
            if _nonempty(value):
                extra[f"edi_{key}"] = value
    if (
        head_project
        and info_project
        and str(head_project) != str(info_project)
    ):
        extra["edi_head_project"] = head_project
    if head_survey and info_survey and str(head_survey) != str(info_survey):
        extra["edi_head_survey"] = head_survey

    year = _year_from_date(start)
    country = (
        getattr(head, "country", None)
        if head is not None and "country" in head_keys
        else None
    )
    acquired_by = (
        getattr(head, "acqby", None)
        if head is not None and "acqby" in head_keys
        else None
    )

    if not any(
        value is not None
        for value in (
            site_id,
            site_name,
            location,
            start,
            end,
            country,
            acquired_by,
            head_project,
            info_project,
            head_survey,
            info_survey,
        )
    ):
        return None

    return SiteMeta(
        project=_first(info_project, head_project),
        survey=_first(info_survey, head_survey),
        year_collected=year,
        country=country,
        site_id=site_id,
        name=site_name,
        location=location,
        acquired_by=acquired_by,
        start=start,
        end=end,
        extra=extra,
    )


def _edi_provenance(edi: Any) -> ProvenanceMeta | None:
    head = edi.get_section("head")
    info = edi.get_section("info")
    head_keys = _explicit_keys(getattr(head, "edi_header", None))
    info_keys = _explicit_keys(getattr(info, "ediinfo", None))
    source = getattr(info, "Source", None) if info is not None else None

    create_time = None
    if head is not None and "filedate" in head_keys:
        create_time = getattr(head, "filedate", None)
    elif source is not None and "creationdate" in info_keys:
        create_time = getattr(source, "creationdate", None)

    creating_application = None
    if head is not None and "progvers" in head_keys:
        creating_application = getattr(head, "progvers", None)

    creator = None
    if head is not None and "fileby" in head_keys:
        name = getattr(head, "fileby", None)
        if _nonempty(name):
            creator = Person(name=str(name))

    if (
        create_time is None
        and creating_application is None
        and creator is None
    ):
        return None
    return ProvenanceMeta(
        create_time=create_time,
        creating_application=creating_application,
        creator=creator,
    )


def _edi_processing(edi: Any) -> ProcessingMeta | None:
    info = edi.get_section("info")
    mtsect = edi.get_section("mtsect")
    if info is None and mtsect is None:
        return None

    info_keys = _explicit_keys(getattr(info, "ediinfo", None))
    old = getattr(info, "Processing", None) if info is not None else None

    sign = (
        getattr(old, "signconvention", None)
        if old is not None and "signconvention" in info_keys
        else None
    )
    processed_by = (
        getattr(old, "processedby", None)
        if old is not None and "processedby" in info_keys
        else None
    )
    processing_tag = (
        getattr(old, "processingtag", None)
        if old is not None and "processingtag" in info_keys
        else None
    )

    software = None
    if old is not None and "processingsoftware" in info_keys:
        old_sw = getattr(old, "ProcessingSoftware", None)
        sw_name = getattr(old_sw, "name", None)
        if _nonempty(sw_name):
            software = Software(name=str(sw_name))

    remote_type = (
        getattr(old, "remoteref", None)
        if old is not None and "remoteref" in info_keys
        else None
    )
    remote_site = (
        getattr(old, "remotesite", None)
        if old is not None and "remotesite" in info_keys
        else None
    )
    remote_extra: dict[str, Any] = {}
    if mtsect is not None:
        for key in ("rx", "ry"):
            value = getattr(mtsect, key, None)
            if _nonempty(value):
                remote_extra[f"edi_{key}"] = str(value)
    remote = None
    if remote_type or remote_site or remote_extra:
        if remote_type is None and remote_extra:
            # RX/RY are explicit reference-channel identifiers in MTSECT.
            # Preserve that scientific fact in the format-neutral metadata
            # even when the historical INFO block did not name the method.
            remote_type = "Remote Reference"
        remote = RemoteReferenceMeta(
            reference_type=remote_type,
            site=remote_site,
            extra=remote_extra,
        )

    run_list = None
    if old is not None and "runlist" in info_keys:
        raw = getattr(old, "runlist", None)
        if isinstance(raw, list):
            run_list = [str(item) for item in raw]
        elif _nonempty(raw):
            text = str(raw).replace(",", " ")
            run_list = [item for item in text.split() if item]

    if not any(
        value is not None
        for value in (
            sign,
            processed_by,
            software,
            remote,
            processing_tag,
            run_list,
        )
    ):
        return None
    return ProcessingMeta(
        sign_convention=sign,
        processed_by=processed_by,
        software=software,
        remote_reference=remote,
        processing_tag=processing_tag,
        run_list=run_list,
    )


def _edi_orientation(edi: Any, n_freq: int) -> OrientationMeta | None:
    head = edi.get_section("head")
    coordsys = getattr(head, "coordsys", None) if head is not None else None
    zrot = _rotation_vector(getattr(edi.Z, "rotation_angle", None), n_freq)
    trot = _rotation_vector(getattr(edi.Tip, "rotation_angle", None), n_freq)
    had_zrot = _raw_file_has_block(edi, "ZROT")
    had_trot = _raw_file_has_block(edi, "TROT")

    if zrot is None and trot is None:
        return None

    extra: dict[str, Any] = {
        "edi_coordsys": coordsys,
        "edi_had_zrot": had_zrot,
        "edi_had_trot": had_trot,
    }
    if zrot is not None:
        extra["edi_zrot"] = zrot.tolist()
    if trot is not None:
        extra["edi_trot"] = trot.tolist()

    mode = None
    angle = None
    rotation_info = None
    if zrot is not None and zrot.size and np.allclose(zrot, zrot[0]):
        coordinate_text = str(coordsys or "").lower()
        if "geographic" in coordinate_text:
            mode = "orthogonal"
            angle = float(zrot[0])
        else:
            rotation_info = (
                "EDI ZROT is constant at "
                f"{float(zrot[0]):g} deg relative to "
                f"HEAD.COORDSYS={coordsys!r}; it was not promoted to an "
                "angle_to_geographic_north because historical EDI rotation "
                "metadata can be ambiguous."
            )
    elif zrot is not None:
        rotation_info = (
            "EDI contains frequency-dependent ZROT values. They are retained "
            "in OrientationMeta.extra for in-memory EDI round-tripping; "
            "Phase 7 rotation is required before archival normalization."
        )

    return OrientationMeta(
        mode=mode,
        angle_to_geographic_north=angle,
        rotation_info=rotation_info,
        extra=extra,
    )


def _edi_field_notes(edi: Any) -> dict[str, Any]:
    info = edi.get_section("info")
    lines = list(getattr(info, "info_text", []) or []) if info else []
    values = [line.rstrip("\n") for line in lines if str(line).strip()]
    return {"edi_info": values} if values else {}


def _edi_auxiliary_metadata(edi: Any) -> dict[str, Any]:
    head = edi.get_section("head")
    info = edi.get_section("info")
    mtsect = edi.get_section("mtsect")
    define = edi.get_section("definemeas")
    result: dict[str, Any] = {
        "source_format": "edi",
        "edi_dtype": getattr(edi, "dtype", None),
    }
    if getattr(edi, "path", None) is not None:
        result["source_path"] = str(edi.path)
    if head is not None:
        result["edi_head"] = getattr(head, "as_dict", lambda: {})()
        result["edi_head_explicit"] = list(
            getattr(head, "edi_header", []) or []
        )
    if info is not None:
        result["edi_info"] = getattr(info, "as_dict", lambda: {})()
        result["edi_info_explicit"] = list(
            getattr(info, "ediinfo", []) or []
        )
    if mtsect is not None:
        result["edi_mtsect"] = {
            key: getattr(mtsect, key, None)
            for key in (
                "sectid",
                "nfreq",
                "maxblks",
                "hx",
                "hy",
                "hz",
                "ex",
                "ey",
                "rx",
                "ry",
            )
        }
    if define is not None:
        result["edi_definemeas"] = {
            key: getattr(define, key, None)
            for key in (
                "maxchan",
                "maxrun",
                "maxmeas",
                "units",
                "reftype",
                "reflat",
                "reflong",
                "refelev",
            )
        }
        measurements: list[dict[str, Any]] = []
        for item in list(getattr(define, "hmeas", []) or []):
            payload = getattr(item, "to_dict", lambda: {})()
            measurements.append({"kind": "H", **dict(payload)})
        for item in list(getattr(define, "emeas", []) or []):
            payload = getattr(item, "to_dict", lambda: {})()
            measurements.append({"kind": "E", **dict(payload)})
        if measurements:
            # Keep the complete historical measurement declarations as an
            # adapter-level escape hatch. SiteLayout holds the scientifically
            # important local channels; these raw records additionally retain
            # remote-reference HMEAS entries needed for an EDI round-trip.
            result["edi_measurements"] = measurements
    return result


def _ensure_edi_file(source: Any) -> Any:
    from ...seg.edi import EDIFile

    if isinstance(source, EDIFile):
        return source
    if isinstance(source, (str, PathLike)):
        return EDIFile(source)
    raise TypeError("EDI source must be an EDIFile or filesystem path")


def edi_to_emtf(
    source: Any,
    *,
    prefer_spectra: bool = True,
    spectra_nfreq_policy: str = "raise",
    spectra_missing_policy: str = "raise",
    spectra_avgt_policy: str = "raise",
) -> EMTF:
    """Convert a historical SEG EDI object/path into :class:`EMTF`.

    EDI SPECTRA are preferred when present because they retain enough
    cross-power information to recover the full inverse-signal and residual
    covariance matrices. Set ``prefer_spectra=False`` to force the historical
    impedance/tipper blocks.

    Notes
    -----
    pyCSAMT historically stores ``Z.z_err = sqrt(EDI complex variance)``.
    The non-SPECTRA path therefore reconstructs EMTF ``VAR`` as ``z_err**2``
    exactly, preserving the EDI file's complex-variance convention without
    redefining the public ``Z.z_err`` API. The standard error of a real or
    imaginary component, when needed statistically, is ``sqrt(VAR / 2)``.
    """
    edi = _ensure_edi_file(source)

    spectra = edi.get_section("spectra")
    if (
        prefer_spectra
        and spectra is not None
        and int(getattr(spectra, "n_freq", 0)) > 0
    ):
        from .spectra import spectra_to_emtf

        return spectra_to_emtf(
            edi,
            spectra=spectra,
            nfreq_policy=spectra_nfreq_policy,
            missing_policy=spectra_missing_policy,
            avgt_policy=spectra_avgt_policy,
        )

    freq = _array_or_none(getattr(edi.Z, "freq", None))
    if freq is None:
        if spectra is not None and int(getattr(spectra, "n_freq", 0)) > 0:
            raise EMTFEDIConversionError(
                "EDI contains SPECTRA but no parsed impedance blocks; use "
                "prefer_spectra=True for FCU-compatible recovery"
            )
        raise EMTFEDIConversionError(
            "EDI conversion requires parsed impedance/tipper or usable "
            "SPECTRA content"
        )
    freq = np.asarray(freq, dtype=float)
    periods = 1.0 / freq

    document = EMTF(
        subtype="MT_TF" if getattr(edi, "dtype", None) == "mt" else None,
        periods=periods,
        site=_edi_site(edi),
        site_layout=_edi_site_layout(edi),
        orientation=_edi_orientation(edi, freq.size),
        processing=_edi_processing(edi),
        provenance=_edi_provenance(edi),
        field_notes=_edi_field_notes(edi),
        metadata=_edi_auxiliary_metadata(edi),
        attrs={"source_format": "edi"},
    )

    z = _array_or_none(getattr(edi.Z, "z", None))
    if z is not None:
        ztf = TransferFunction(
            name="impedance",
            data=np.asarray(z),
            input_channels=IMPEDANCE_INPUT_CHANNELS,
            output_channels=IMPEDANCE_OUTPUT_CHANNELS,
            periods=periods,
            attrs={
                "source_format": "edi",
                "edi_error_carrier": "sqrt_complex_variance",
            },
        )
        zerr = _array_or_none(getattr(edi.Z, "z_err", None))
        if zerr is not None:
            ztf.add_estimate(
                StatisticalEstimate(
                    name="VAR",
                    kind="variance",
                    data=np.square(np.asarray(zerr, dtype=float)),
                    attrs={
                        "source_format": "edi",
                        "definition": "complex_transfer_function_variance",
                    },
                )
            )
        document.add_transfer_function(ztf, replace=True)

    tip = _array_or_none(getattr(edi.Tip, "tipper", None))
    if tip is not None:
        tip = np.asarray(tip)
        if tip.ndim == 3 and tip.shape[0] == freq.size:
            ttf = TransferFunction(
                name="tipper",
                data=tip,
                input_channels=TIPPER_INPUT_CHANNELS,
                output_channels=TIPPER_OUTPUT_CHANNELS,
                periods=periods,
                attrs={
                    "source_format": "edi",
                    "edi_error_carrier": "sqrt_complex_variance",
                },
            )
            terr = _array_or_none(getattr(edi.Tip, "tipper_err", None))
            if terr is not None:
                ttf.add_estimate(
                    StatisticalEstimate(
                        name="VAR",
                        kind="variance",
                        data=np.square(np.asarray(terr, dtype=float)),
                        attrs={
                            "source_format": "edi",
                            "definition": (
                                "complex_transfer_function_variance"
                            ),
                        },
                    )
                )
            document.add_transfer_function(ttf, replace=True)

    if document.impedance is None:
        raise EMTFEDIConversionError(
            "EDI did not expose a standard impedance tensor; conversion of "
            "spectra-only or nonstandard EDI data is deferred to Phase 8"
        )
    document.tags = tuple(document.transfer_functions)
    return document


def _extract_edi_info_lines(document: EMTF) -> list[str]:
    raw = document.field_notes.get("edi_info")
    values = (
        raw if isinstance(raw, list) else ([raw] if raw is not None else [])
    )
    lines: list[str] = []
    for value in values:
        if isinstance(value, str):
            lines.append(value)
            continue
        if isinstance(value, dict):
            text = value.get("#text")
            if text is not None:
                lines.append(str(text))
    return lines


def _preferred_variance(tf: TransferFunction) -> np.ndarray | None:
    estimate = tf.get_estimate("variance") or tf.get_estimate("VAR")
    if estimate is None:
        return None
    data = np.asarray(estimate.data)
    if data.shape != tf.data.shape:
        raise EMTFEDIConversionError(
            f"{tf.name} VAR shape {data.shape!r} does not match TF shape "
            f"{tf.data.shape!r}"
        )
    return np.asarray(data, dtype=float)


def _legacy_error(tf: TransferFunction) -> np.ndarray | None:
    estimate = tf.get_estimate(LEGACY_STANDARD_ERROR_KIND)
    if estimate is None:
        return None
    data = np.asarray(estimate.data, dtype=float)
    if data.shape != tf.data.shape:
        raise EMTFEDIConversionError(
            f"legacy standard-error shape {data.shape!r} does not match "
            f"{tf.data.shape!r}"
        )
    semantics = str(estimate.attrs.get("semantics", "")).strip().lower()
    if semantics != "pycsamt_legacy_z_err":
        raise EMTFEDIConversionError(
            "a generic standard error cannot be mapped to EDI *.VAR without "
            "knowing its statistical convention; expected explicit "
            "semantics='pycsamt_legacy_z_err'"
        )
    return data


def _tf_to_edi_error(tf: TransferFunction) -> np.ndarray | None:
    variance = _preferred_variance(tf)
    if variance is not None:
        # EDIFile.write squares its historical z_err/tipper_err carrier.
        # sqrt(VAR) therefore preserves the complex EDI variance exactly.
        return np.sqrt(np.abs(variance))
    return _legacy_error(tf)


def _collect_unrepresentable(document: EMTF) -> list[str]:
    losses: list[str] = []
    representable = {"impedance", "tipper"}
    extras = [
        key for key in document.transfer_functions if key not in representable
    ]
    if extras:
        losses.append(
            "EDI writer will omit transfer-function/data types: "
            + ", ".join(sorted(extras))
        )

    for key in ("impedance", "tipper"):
        tf = document.get_transfer_function(key)
        if tf is None:
            continue
        covariance = [
            estimate.name
            for estimate in tf.estimates.values()
            if estimate.kind
            in {
                "inverse_signal_covariance",
                "residual_covariance",
                "covariance",
            }
        ]
        if covariance:
            losses.append(
                f"standard EDI cannot preserve {key} full covariance "
                f"estimate(s): {', '.join(sorted(set(covariance)))}"
            )
        supported_names: set[str] = set()
        variance = tf.get_estimate("variance") or tf.get_estimate("VAR")
        if variance is not None:
            supported_names.add(variance.name)
        legacy = tf.get_estimate(LEGACY_STANDARD_ERROR_KIND)
        if variance is None and legacy is not None:
            semantics = str(
                legacy.attrs.get("semantics", "")
            ).strip().lower()
            if semantics == "pycsamt_legacy_z_err":
                supported_names.add(legacy.name)
        other_estimates = [
            estimate.name
            for estimate in tf.estimates.values()
            if estimate.name not in supported_names
            and estimate.kind
            not in {
                "inverse_signal_covariance",
                "residual_covariance",
                "covariance",
            }
        ]
        if other_estimates:
            losses.append(
                f"standard EDI conversion will omit unsupported {key} "
                "statistical estimate(s): "
                + ", ".join(sorted(set(other_estimates)))
            )

    if document.copyright is not None:
        losses.append(
            "current SEG EDI writer has no lossless mapping for EMTF "
            "Copyright/Citation metadata"
        )
    if document.quality is not None:
        losses.append(
            "current SEG EDI writer has no lossless mapping for EMTF "
            "DataQualityNotes/DataQualityWarnings"
        )
    provenance = document.provenance
    if provenance is not None:
        if provenance.submitter is not None:
            losses.append("EDI cannot preserve EMTF Submitter provenance")
        creator = provenance.creator
        if creator is not None and any(
            (
                creator.email,
                creator.organization,
                creator.organization_url,
            )
        ):
            losses.append(
                "EDI FILEBY can preserve only the creator name, not full "
                "creator contact metadata"
            )
    software = document.processing.software if document.processing else None
    if software is not None and any(
        (
            software.version,
            software.release,
            software.author.name if software.author else None,
            software.author.email if software.author else None,
        )
    ):
        losses.append(
            "EDI INFO PROCESSINGSOFTWARE preserves only the software name"
        )
    return losses


def _channel_lookup(document: EMTF) -> dict[str, ChannelMeta]:
    layout = document.site_layout
    if layout is None:
        return {}
    return {
        channel.name.lower(): channel
        for channel in (*layout.input_channels, *layout.output_channels)
    }


def _measurement_id(
    channel: ChannelMeta | None,
    *,
    fallback: str,
) -> str:
    if channel is not None:
        value = channel.extra.get("edi_id")
        if _nonempty(value):
            return str(value)
    return fallback


def _build_measurements(document: EMTF, n_freq: int) -> tuple[Any, Any]:
    from ...seg.meas import DefineMeas, Emeasurement, Hmeasurement
    from ...seg.mtemap import MTEMAP

    channels = _channel_lookup(document)
    units = None
    if document.site_layout is not None:
        units = _first(
            document.site_layout.input_units,
            document.site_layout.output_units,
        )
    units = str(units or "m")

    ids = {
        "hx": _measurement_id(channels.get("hx"), fallback="1001.001"),
        "hy": _measurement_id(channels.get("hy"), fallback="1002.001"),
        "hz": _measurement_id(channels.get("hz"), fallback="1003.001"),
        "ex": _measurement_id(channels.get("ex"), fallback="1004.001"),
        "ey": _measurement_id(channels.get("ey"), fallback="1005.001"),
    }

    def magnetic(name: str, key: str) -> Hmeasurement:
        channel = channels.get(name.lower())
        kwargs: dict[str, Any] = {
            "id": ids[key],
            "chtype": name.upper(),
            "x": None,
            "y": None,
            "z": None,
            "azm": None,
            "dip": None,
            "acqchan": None,
            "filter": None,
            "sensor": None,
            "gain": None,
            "measdate": None,
        }
        if channel is not None:
            kwargs.update(
                x=channel.x,
                y=channel.y,
                z=channel.z,
                azm=channel.orientation,
                dip=channel.tilt,
                sensor=channel.sensor_id,
                acqchan=channel.extra.get("acquisition_channel"),
                filter=channel.extra.get("filter"),
                gain=channel.extra.get("gain"),
                measdate=channel.extra.get("measurement_date"),
            )
        return Hmeasurement(**kwargs)

    def electric(name: str, key: str) -> Emeasurement:
        channel = channels.get(name.lower())
        kwargs: dict[str, Any] = {
            "id": ids[key],
            "chtype": name.upper(),
            "x": None,
            "y": None,
            "z": None,
            "x2": None,
            "y2": None,
            "z2": None,
            "acqchan": None,
            "filter": None,
            "sensor": None,
            "gain": None,
            "measdate": None,
        }
        if channel is not None:
            kwargs.update(
                x=channel.x,
                y=channel.y,
                z=channel.z,
                x2=channel.x2,
                y2=channel.y2,
                z2=channel.z2,
                sensor=channel.sensor_id,
                acqchan=channel.extra.get("acquisition_channel"),
                filter=channel.extra.get("filter"),
                gain=channel.extra.get("gain"),
                measdate=channel.extra.get("measurement_date"),
            )
        return Emeasurement(**kwargs)

    original_define = document.metadata.get("edi_definemeas", {})
    if not isinstance(original_define, dict):
        original_define = {}
    define = DefineMeas(
        maxchan=5,
        maxrun=original_define.get("maxrun") or 999,
        maxmeas=5,
        units=units,
        reftype=original_define.get("reftype") or "CART",
        reflat=original_define.get("reflat"),
        reflong=original_define.get("reflong"),
        refelev=original_define.get("refelev"),
    )
    define.hmeas = [
        magnetic("Hx", "hx"),
        magnetic("Hy", "hy"),
        magnetic("Hz", "hz"),
    ]
    define.emeas = [electric("Ex", "ex"), electric("Ey", "ey")]

    site_id = (
        document.site.site_id
        if document.site is not None
        else document.station
    )
    mtsect = MTEMAP(
        sectid=str(site_id) if site_id is not None else None,
        nfreq=n_freq,
        hx=ids["hx"],
        hy=ids["hy"],
        hz=ids["hz"],
        ex=ids["ex"],
        ey=ids["ey"],
    )

    original = document.metadata.get("edi_mtsect", {})
    if not isinstance(original, dict):
        original = {}
    remote = (
        document.processing.remote_reference
        if document.processing is not None
        else None
    )

    requested_remote: dict[str, str] = {}
    for key in ("rx", "ry"):
        value = original.get(key)
        if remote is not None and _nonempty(remote.extra.get(f"edi_{key}")):
            value = remote.extra.get(f"edi_{key}")
        if _nonempty(value):
            requested_remote[key] = str(value)

    raw_records = document.metadata.get("edi_measurements", [])
    if not isinstance(raw_records, list):
        raw_records = []
    by_id: dict[str, dict[str, Any]] = {}
    for record in raw_records:
        if not isinstance(record, dict):
            continue
        ident = record.get("id")
        if _nonempty(ident):
            by_id[str(ident).strip()] = record

    existing_ids = {
        str(item.id).strip()
        for item in define.all_meas()
        if _nonempty(getattr(item, "id", None))
    }

    # An EDI-derived EMTF keeps all original HMEAS/EMEAS records in adapter
    # metadata. Restore records that are not already represented by the
    # canonical local SiteLayout so a direct EDI -> EMTF -> EDI conversion
    # does not discard auxiliary/remote measurement declarations.
    for record in raw_records:
        if not isinstance(record, dict):
            continue
        ident = record.get("id")
        if not _nonempty(ident) or str(ident).strip() in existing_ids:
            continue
        payload = {
            field: value
            for field, value in record.items()
            if field != "kind"
        }
        kind = str(record.get("kind", "H")).upper()
        item = (
            Hmeasurement(**payload)
            if kind == "H"
            else Emeasurement(**payload)
        )
        if kind == "H":
            define.hmeas.append(item)
        else:
            define.emeas.append(item)
        existing_ids.add(str(ident).strip())

    restored_remote = 0
    missing_remote: list[str] = []
    for key, ident in requested_remote.items():
        record = by_id.get(ident)
        if record is None:
            missing_remote.append(f"{key.upper()}={ident}")
            continue
        if ident not in existing_ids:
            payload = {
                field: value
                for field, value in record.items()
                if field != "kind"
            }
            kind = str(record.get("kind", "H")).upper()
            if kind != "H":
                missing_remote.append(f"{key.upper()}={ident}")
                continue
            define.hmeas.append(Hmeasurement(**payload))
            existing_ids.add(ident)
            restored_remote += 1
        setattr(mtsect, key, ident)

    # Remote IDs without their HMEAS declarations would create structurally
    # inconsistent EDI. Do not write dangling references. The caller turns
    # this marker into an explicit DataLossWarning.
    if missing_remote:
        setattr(mtsect, "_pycsamt_missing_remote_measurements", missing_remote)

    define.maxchan = len(define.all_meas())
    define.maxmeas = len(define.all_meas())
    if restored_remote:
        define.maxchan = max(define.maxchan, 5 + restored_remote)

    return define, mtsect


def _build_head_info(document: EMTF) -> tuple[Any, Any]:
    from ...seg.heads import Head, Info

    site = document.site
    location = site.location if site is not None else None
    station = _first(
        site.site_id if site is not None else None,
        document.station,
        document.product_id,
    )
    head = Head()
    # ``Head`` has convenience defaults intended for newly authored EDI
    # files. A format converter must not silently claim geomagnetic
    # coordinates or invent file/program provenance, so clear those defaults
    # and repopulate only mapped values. Structural defaults such as STDVERS,
    # datum/units, and EMPTY remain available when needed by the EDI writer.
    head.filedate = None
    head.progvers = None
    head.progdate = None
    head.coordsys = None
    if station is not None:
        head.dataid = str(station)
    if site is not None:
        head.acqby = site.acquired_by
        head.acqdate = _as_edi_date(site.start)
        head.enddate = _as_edi_date(site.end)
        head.country = site.country
        head.project = site.project
        head.survey = site.survey
        head.loc = site.name
        for source, target in (
            ("edi_state", "state"),
            ("edi_county", "county"),
            ("edi_prospect", "prospect"),
            ("edi_chainage", "chainage"),
        ):
            value = site.extra.get(source)
            if _nonempty(value):
                setattr(head, target, value)
    if location is not None:
        head.lat = location.latitude
        head.long = location.longitude
        head.elev = location.elevation
        head.units = _elevation_units_to_edi(location.elevation_units)
        head.datum = location.datum or "WGS84"
        head.declination = location.declination
    if document.provenance is not None:
        creator = document.provenance.creator
        if creator is not None:
            head.fileby = creator.name
        head.filedate = _as_edi_date(document.provenance.create_time)
        head.progvers = document.provenance.creating_application

    orientation = document.orientation
    if orientation is not None:
        source_coordsys = orientation.extra.get("edi_coordsys")
        if _nonempty(source_coordsys):
            head.coordsys = str(source_coordsys)
        elif orientation.is_orthogonal:
            head.coordsys = "Geographic North"

    # Restore EDI-specific header details only when this document actually
    # originated from EDI and the key was explicit in the source. Scientific
    # metadata above remain authoritative for fields shared with EMTF.
    raw_head = document.metadata.get("edi_head", {})
    explicit_head = _explicit_keys(
        document.metadata.get("edi_head_explicit", [])
    )
    if isinstance(raw_head, dict):
        for key in (
            "filedate",
            "stdvers",
            "progdate",
            "maxsect",
            "bindata",
            "empty",
        ):
            if key in explicit_head and _nonempty(raw_head.get(key)):
                setattr(head, key, raw_head[key])
        if "progvers" in explicit_head and document.provenance is None:
            head.progvers = raw_head.get("progvers")
        if "coordsys" in explicit_head and orientation is None:
            head.coordsys = raw_head.get("coordsys")

    info = Info()
    # Avoid Info's convenience processing defaults becoming fabricated
    # scientific provenance during conversion.
    info.Processing.processedby = None
    info.Processing.ProcessingSoftware.name = None
    if site is not None:
        info.Source.project = site.project
        info.Source.survey = site.survey
        info.Source.sitename = site.name
    if document.provenance is not None:
        if document.provenance.create_time is not None:
            info.Source.creationdate = str(document.provenance.create_time)
    if document.processing is not None:
        processing = document.processing
        info.Processing.processedby = processing.processed_by
        info.Processing.processingtag = processing.processing_tag
        info.Processing.signconvention = processing.sign_convention
        info.Processing.runlist = processing.run_list
        if processing.software is not None:
            info.Processing.ProcessingSoftware.name = processing.software.name
        if processing.remote_reference is not None:
            info.Processing.remoteref = (
                processing.remote_reference.reference_type
            )
            info.Processing.remotesite = processing.remote_reference.site
    info_lines = _extract_edi_info_lines(document)
    raw_info = document.metadata.get("edi_info", {})
    explicit_info = list(document.metadata.get("edi_info_explicit", []) or [])
    explicit_keys = _explicit_keys(explicit_info)
    if isinstance(raw_info, dict) and "maxinfo" in explicit_keys:
        value = raw_info.get("maxinfo")
        if value is not None:
            try:
                info.maxinfo = int(value)
            except (TypeError, ValueError):
                pass

    known_info = {
        "maxinfo",
        "project",
        "survey",
        "creationdate",
        "processedby",
        "processingsoftware",
        "processingtag",
        "sitename",
        "runlist",
        "remoteref",
        "remotesite",
        "signconvention",
    }
    for raw in explicit_info:
        text = str(raw).strip()
        if "=" not in text:
            continue
        key = text.split("=", 1)[0].strip().lower()
        if key not in known_info:
            # INFO supports free text, so retaining an unknown historical KV
            # line here is preferable to silently dropping it.
            info_lines.append(text)
    info.info_text = info_lines
    return head, info


def _restore_rotation(document: EMTF, edi: Any, n_freq: int) -> None:
    orientation = document.orientation
    zrot = None
    trot = None
    if orientation is not None:
        raw_zrot = orientation.extra.get("edi_zrot")
        raw_trot = orientation.extra.get("edi_trot")
        zrot = _rotation_vector(raw_zrot, n_freq)
        trot = _rotation_vector(raw_trot, n_freq)
        if zrot is None and orientation.is_orthogonal:
            angle = orientation.angle_to_geographic_north
            if angle is not None:
                zrot = np.full(n_freq, float(angle))
    if zrot is None:
        zrot = np.zeros(n_freq, dtype=float)
    if trot is None:
        trot = np.array(zrot, copy=True)
    edi.Z.rotation_angle = zrot
    if getattr(edi.Tip, "tipper", None) is not None:
        edi.Tip.rotation_angle = trot


def emtf_to_edi(
    document: EMTF,
    *,
    on_loss: str = "warn",
) -> Any:
    """Convert :class:`EMTF` into an in-memory historical ``EDIFile``.

    The standard EDI representation used here stores impedance/tipper and
    component variances.  Rich covariance matrices, arbitrary TF types, and
    several EMTF metadata blocks have no lossless destination; ``on_loss``
    controls whether those reductions warn, raise, or are explicitly ignored.
    """
    if not isinstance(document, EMTF):
        raise TypeError("document must be a pycsamt.emtf.EMTF")
    policy = str(on_loss).strip().lower()
    if policy not in _VALID_LOSS_POLICIES:
        raise ValueError(
            "on_loss must be one of " f"{sorted(_VALID_LOSS_POLICIES)}"
        )

    ztf = document.impedance
    if ztf is None:
        raise EMTFEDIConversionError(
            "standard EDI writing requires an impedance transfer function"
        )
    if ztf.data.shape[1:] != (2, 2):
        raise EMTFEDIConversionError(
            "EDI impedance must have matrix shape (n_period, 2, 2)"
        )
    freq = document.frequency
    if freq is None or len(freq) != ztf.n_periods:
        raise EMTFEDIConversionError(
            "EDI writing requires a frequency vector matching impedance"
        )
    if tuple(name.lower() for name in ztf.input_channels) != ("hx", "hy"):
        raise EMTFEDIConversionError(
            "EDI impedance input channels must be Hx/Hy"
        )
    if tuple(name.lower() for name in ztf.output_channels) != ("ex", "ey"):
        raise EMTFEDIConversionError(
            "EDI impedance output channels must be Ex/Ey"
        )

    for message in _collect_unrepresentable(document):
        _loss(policy, message)

    station = _first(
        document.site.site_id if document.site is not None else None,
        document.station,
        document.product_id,
    )
    if station is None:
        _loss(
            policy,
            "EMTF has no site identifier; EDI requires DATAID and would need "
            "a synthesized identifier",
        )
        station = "site"

    from ...seg.edi import EDIFile
    from ...z.tipper import Tipper
    from ...z.z import Z

    edi = EDIFile()
    head, info = _build_head_info(document)
    if head.dataid is None:
        head.dataid = str(station)
    define, mtsect = _build_measurements(document, len(freq))
    if mtsect.sectid is None:
        mtsect.sectid = str(station)
    missing_remote = getattr(
        mtsect, "_pycsamt_missing_remote_measurements", []
    )
    if missing_remote:
        _loss(
            policy,
            "EDI remote-reference IDs could not be restored because their "
            "HMEAS declarations are unavailable: "
            + ", ".join(missing_remote)
            + ". The dangling RX/RY references were omitted.",
        )

    edi.add_section("head", head)
    edi.add_section("info", info)
    edi.add_section("definemeas", define)
    edi.add_section("mtsect", mtsect)

    zerr = _tf_to_edi_error(ztf)
    # Build Z without uncertainty first. The compatibility Z class computes
    # rho/phase during ``freq`` assignment and requires finite errors, while
    # EMTF legitimately allows missing variance components. The EDI serializer
    # can represent those missing uncertainties using EMPTY, so attach the
    # legacy error carrier after Z initialization.
    edi.Z = Z(
        z_array=np.array(ztf.data, copy=True),
        freq=np.asarray(freq, dtype=float),
        name=str(station),
    )
    if zerr is not None:
        edi.Z._z_err = np.array(zerr, copy=True)

    ttf = document.tipper_tf
    if ttf is not None:
        if ttf.data.shape[1:] != (1, 2):
            raise EMTFEDIConversionError(
                "EDI tipper must have matrix shape (n_period, 1, 2)"
            )
        if tuple(name.lower() for name in ttf.input_channels) != ("hx", "hy"):
            raise EMTFEDIConversionError("EDI tipper inputs must be Hx/Hy")
        if tuple(name.lower() for name in ttf.output_channels) != ("hz",):
            raise EMTFEDIConversionError("EDI tipper output must be Hz")
        terr = _tf_to_edi_error(ttf)
        edi.Tip = Tipper(
            tipper_array=np.array(ttf.data, copy=True),
            tipper_err_array=(
                None if terr is None else np.array(terr, copy=True)
            ),
            freq=np.asarray(freq, dtype=float),
            name=str(station),
        )
        edi._force_tipper_on_write = True

    _restore_rotation(document, edi, len(freq))
    return edi


def write_edi(
    obj: Any,
    target: str | PathLike[str],
    *,
    on_loss: str = "warn",
    preserve_zero: bool = True,
    stamp_headers: bool = False,
    force_tipper: bool | None = None,
    **kwargs: Any,
) -> Path:
    """Write an ``EMTF`` or ``EDIFile`` to an exact filesystem target.

    Existing :meth:`EDIFile.write` behavior remains unchanged by default.
    This neutral adapter opts into ``preserve_zero=True`` so finite physical
    zeros are not turned into the EDI missing-value sentinel during format
    conversion; NaN/Inf values are still serialized as the sentinel.
    """
    if not isinstance(target, (str, PathLike)):
        raise TypeError("EDI writing currently requires a filesystem target")
    target_path = Path(target).expanduser()
    parent = target_path.parent if str(target_path.parent) else Path.cwd()
    parent.mkdir(parents=True, exist_ok=True)

    from ...seg.edi import EDIFile

    if isinstance(obj, EMTF):
        edi = emtf_to_edi(obj, on_loss=on_loss)
        if force_tipper is None:
            force_tipper = obj.tipper_tf is not None
    elif isinstance(obj, EDIFile):
        edi = obj
        if force_tipper is None:
            force_tipper = False
    else:
        raise TypeError("EDI writer accepts EMTF or EDIFile objects")

    result = Path(
        edi.write(
            new_edifn=target_path.name,
            savepath=parent,
            preserve_zero=preserve_zero,
            stamp_headers=stamp_headers,
            force_tipper=bool(force_tipper),
            **kwargs,
        )
    )
    resolved_target = target_path.resolve()
    if result.resolve() != resolved_target:
        result.replace(resolved_target)
    return resolved_target
