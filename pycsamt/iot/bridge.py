"""Bridge between IoT field telemetry and the pyCSAMT data model.

The rest of the :mod:`pycsamt.iot` subpackage is deliberately lean: it
processes raw arrays and telemetry with numpy only and never imports the
heavier pyCSAMT data model. That keeps ``import pycsamt.iot`` cheap on a
field gateway. This module is the one place where the two worlds meet, so
every data-model import (:mod:`pycsamt.seg`, :mod:`pycsamt.site`) is done
lazily inside the functions that need it.

Two directions are supported.

**Forward — acquisition to data model.** Turn edge impedance estimates
into first-class objects that the processing/inversion flow understands:

* :func:`impedance_to_z` — per-window edge impedance to a :class:`~pycsamt.z.z.Z`,
* :func:`build_edifile` / :func:`z_to_edi` — a :class:`~pycsamt.z.z.Z` to an
  in-memory ``EDIFile`` or a written ``.edi`` file,
* :func:`field_session_to_edifiles` — a :class:`~pycsamt.iot.session.FieldSession`
  plus impedance to per-station EDIs enriched with station geometry,
* :func:`z_to_site` — a :class:`~pycsamt.z.z.Z` to an EDI-backed
  :class:`~pycsamt.site.base.Site` (needs the optional geospatial stack).

**Reverse — data model to acquisition planning.** Seed an IoT deployment
from a previous survey so known station geometry can be re-occupied:

* :func:`read_edi_survey` / :func:`edi_survey_table` — summarise an EDI survey,
* :func:`field_session_from_edis` — a ready-to-occupy ``FieldSession``,
* :func:`deployment_from_edis` — a :class:`~pycsamt.iot.core.DeploymentConfig`.

Example
-------
>>> import numpy as np
>>> from pycsamt.iot.bridge import impedance_to_z, z_to_edi
>>> freq = np.logspace(4, 0, 12)
>>> zxy = (1 + 1j) * np.sqrt(freq)          # one scalar sounding
>>> z = impedance_to_z(zxy, freq, station="S01")
>>> path = z_to_edi(z, station="S01", lat=6.5, lon=3.4,   # doctest: +SKIP
...                 savepath="edi_out")
"""

from __future__ import annotations

import glob
import os
from collections.abc import Iterable, Mapping, Sequence
from typing import Any

import numpy as np
import pandas as pd

from ..api.view import maybe_wrap_frame
from . import _common as _c

__all__ = [
    "impedance_to_z",
    "build_edifile",
    "z_to_edi",
    "field_session_to_edifiles",
    "z_to_site",
    "emtools_qc",
    "read_edi_survey",
    "edi_survey_table",
    "field_session_from_edis",
    "deployment_from_edis",
]

_DEFAULT_CHANNELS = ("ex", "ey", "hx", "hy")


# ---------------------------------------------------------------------------
# forward: impedance -> Z
# ---------------------------------------------------------------------------
def _aggregate_windows(
    arr: np.ndarray, how: str
) -> tuple[np.ndarray, np.ndarray]:
    """Collapse a window axis (axis 0), returning ``(value, std_error)``.

    The standard deviation across windows is a natural absolute error for
    the aggregated impedance and is exactly what an edge node can report
    cheaply alongside a running mean.
    """
    how = str(how).strip().lower()
    if how not in {"mean", "median"}:
        raise ValueError("aggregate must be 'mean' or 'median'.")
    n_windows = arr.shape[0]
    if how == "median":
        # Median of complex data: take it per component on real/imag.
        value = np.median(arr.real, axis=0) + 1j * np.median(arr.imag, axis=0)
    else:
        value = np.mean(arr, axis=0)
    if n_windows >= 2:
        # Magnitude of the complex standard deviation, a real error scale.
        err = np.std(arr, axis=0)
        err = np.abs(err).astype(float)
    else:
        err = np.zeros(value.shape, dtype=float)
    return value, err


def _scalar_to_tensor(
    scalar: np.ndarray,
    err: np.ndarray | None,
    component: str,
) -> tuple[np.ndarray, np.ndarray | None]:
    """Place a per-frequency scalar impedance into a ``(n, 2, 2)`` stack."""
    component = str(component).strip().lower()
    n = scalar.shape[0]
    z = np.zeros((n, 2, 2), dtype=complex)
    z_err = None if err is None else np.zeros((n, 2, 2), dtype=float)
    if component in {"offdiag", "off-diagonal", "xy-yx"}:
        z[:, 0, 1] = scalar
        z[:, 1, 0] = -scalar
        if z_err is not None:
            z_err[:, 0, 1] = err
            z_err[:, 1, 0] = err
    elif component == "xy":
        z[:, 0, 1] = scalar
        if z_err is not None:
            z_err[:, 0, 1] = err
    elif component == "yx":
        z[:, 1, 0] = scalar
        if z_err is not None:
            z_err[:, 1, 0] = err
    else:
        raise ValueError(
            "component must be one of: 'offdiag', 'xy', 'yx'."
        )
    return z, z_err


def impedance_to_z(
    impedance: Any,
    freq: Any,
    *,
    impedance_err: Any = None,
    station: str | None = None,
    method: str | None = None,
    component: str = "offdiag",
    aggregate: str = "mean",
) -> Any:
    r"""Build a :class:`~pycsamt.z.z.Z` from edge impedance estimates.

    This is the natural continuation of the edge diagnostics: the field
    node already forms per-window impedance estimates (see
    :func:`~pycsamt.iot.edge_amt.assess_impedance_stability`); this turns
    them into the impedance-tensor container the processing and inversion
    flow consumes.

    Parameters
    ----------
    impedance : array-like of complex
        Accepted shapes, with ``n = len(freq)``:

        * ``(n,)`` — one scalar impedance per frequency, placed according
          to *component* (an off-diagonal antisymmetric tensor by default);
        * ``(n, 2, 2)`` — a full impedance tensor per frequency;
        * ``(n_windows, n)`` — per-window scalar estimates, aggregated
          across windows (error taken from the window spread);
        * ``(n_windows, n, 2, 2)`` — per-window tensors, aggregated;
        * ``(2, 2)`` — a single tensor (only when ``n == 1``).
    freq : array-like of float
        Strictly positive frequency vector in Hz.
    impedance_err : array-like, optional
        Explicit absolute errors matching the *aggregated* impedance
        shape. Overrides any error derived from window spread.
    station : str, optional
        Name attached to the returned :class:`~pycsamt.z.z.Z`.
    method : str, optional
        Acquisition method (``"amt"``, ``"mt"``, ``"csamt"``) recorded in
        the object metadata.
    component : {'offdiag', 'xy', 'yx'}
        Placement for scalar inputs.
    aggregate : {'mean', 'median'}
        Reduction across a window axis when one is present.

    Returns
    -------
    pycsamt.z.z.Z
    """
    from ..z.z import Z

    f = np.asarray(freq, dtype=float).ravel()
    if f.size == 0:
        raise ValueError("freq cannot be empty.")
    if np.any(~np.isfinite(f)) or np.any(f <= 0):
        raise ValueError("freq must be finite and strictly positive.")
    n = f.size

    arr = np.asarray(impedance, dtype=complex)
    z_err: np.ndarray | None = None

    if arr.ndim == 4:
        if arr.shape[1:] != (n, 2, 2):
            raise ValueError(
                f"4-D impedance must have shape (n_windows, {n}, 2, 2)."
            )
        z_stack, z_err = _aggregate_windows(arr, aggregate)
    elif arr.ndim == 3:
        if arr.shape != (n, 2, 2):
            raise ValueError(f"3-D impedance must have shape ({n}, 2, 2).")
        z_stack = arr
    elif arr.ndim == 2:
        if arr.shape == (2, 2) and n == 1:
            z_stack = arr[None, ...]
        elif arr.shape[1] == n:
            scalar, serr = _aggregate_windows(arr, aggregate)
            z_stack, z_err = _scalar_to_tensor(scalar, serr, component)
        else:
            raise ValueError(
                "2-D impedance must be (n_windows, n_freq) scalar windows "
                "or a single (2, 2) tensor when n_freq == 1."
            )
    elif arr.ndim == 1:
        if arr.shape[0] != n:
            raise ValueError(
                f"1-D impedance must match len(freq) = {n}."
            )
        z_stack, _ = _scalar_to_tensor(arr, None, component)
    else:
        raise ValueError("impedance has an unsupported number of dimensions.")

    if impedance_err is not None:
        supplied = np.asarray(impedance_err, dtype=float)
        if supplied.shape == z_stack.shape:
            z_err = supplied
        elif supplied.ndim == 1 and supplied.shape[0] == n:
            z_err, _ = _scalar_to_tensor(
                supplied.astype(complex), None, component
            )
            z_err = np.abs(z_err).astype(float)
        else:
            raise ValueError(
                "impedance_err shape is incompatible with the impedance."
            )

    meta: dict[str, Any] = {}
    if method is not None:
        meta["method"] = _c.as_optional_str(method, "method")
    z_obj = Z(
        z_array=z_stack,
        z_err_array=z_err,
        freq=f,
        name=_c.as_optional_str(station, "station"),
        meta=meta,
    )
    return z_obj


# ---------------------------------------------------------------------------
# forward: Z -> EDI
# ---------------------------------------------------------------------------
def build_edifile(
    z: Any,
    *,
    station: str,
    lat: float | None = None,
    lon: float | None = None,
    elevation: float | None = None,
    method: str | None = None,
    acqby: str = "pycsamt.iot",
    metadata: Mapping[str, Any] | None = None,
) -> Any:
    """Assemble an in-memory ``EDIFile`` from an impedance tensor.

    The container carries the minimum valid SEG-EDI structure — a
    ``>HEAD`` block (station id and geometry), a ``>=MTSECT`` measurement
    section, and the impedance transfer function — so it round-trips
    through :meth:`pycsamt.seg.edi.EDIFile.write`.

    Parameters
    ----------
    z : pycsamt.z.z.Z
        Impedance tensor for one station.
    station : str
        Station / DATAID label.
    lat, lon, elevation : float, optional
        Station location written to the EDI head.
    method : str, optional
        Acquisition method recorded as an EDI info note.
    acqby : str
        ``ACQBY`` value; defaults to a pyCSAMT-IoT provenance tag.
    metadata : mapping, optional
        Extra head attributes (e.g. ``county``, ``prospect``).

    Returns
    -------
    pycsamt.seg.edi.EDIFile
    """
    from ..seg.edi import EDIFile
    from ..seg.heads import Head
    from ..seg.mtemap import MTEMAP
    from ..z.z import Z

    if not isinstance(z, Z):
        raise TypeError("z must be a pycsamt.z.z.Z instance.")
    station = _c.as_nonempty_str(station, "station")
    n_freq = int(getattr(z, "n_freq", 0) or 0)
    if n_freq <= 0:
        raise ValueError("z has no frequencies; cannot build an EDI.")

    head = Head()
    head.dataid = station
    lat = _c.as_latitude(lat)
    lon = _c.as_longitude(lon)
    elevation = _c.as_elevation(elevation)
    if lat is not None:
        head.lat = lat
    if lon is not None:
        head.lon = lon
    if elevation is not None:
        head.elev = elevation
    if acqby:
        head.acqby = str(acqby)
    for key, value in dict(metadata or {}).items():
        try:
            setattr(head, str(key), value)
        except Exception:  # pragma: no cover - defensive on unknown attrs
            pass

    ed = EDIFile()
    ed.add_section("head", head)
    m = MTEMAP()
    m.sectid = station
    m.nfreq = n_freq
    ed.add_section("mtsect", m)
    ed.Z = z
    if method is not None:
        info = ed.get_section("info")
        note = f"method={_c.as_optional_str(method, 'method')}"
        try:
            if info is not None and hasattr(info, "notes"):
                info.notes = note
        except Exception:  # pragma: no cover - info section is optional
            pass
    return ed


def z_to_edi(
    z: Any,
    *,
    station: str,
    savepath: str | None = None,
    filename: str | None = None,
    lat: float | None = None,
    lon: float | None = None,
    elevation: float | None = None,
    method: str | None = None,
    datatype: str | None = None,
    acqby: str = "pycsamt.iot",
    metadata: Mapping[str, Any] | None = None,
) -> str:
    """Write an impedance tensor to a ``.edi`` file and return its path.

    See :func:`build_edifile` for the assembled structure. ``datatype``
    is forwarded to :meth:`pycsamt.seg.edi.EDIFile.write` and is
    auto-detected (MT) when left as ``None``.
    """
    ed = build_edifile(
        z,
        station=station,
        lat=lat,
        lon=lon,
        elevation=elevation,
        method=method,
        acqby=acqby,
        metadata=metadata,
    )
    return ed.write(
        new_edifn=filename or station,
        savepath=savepath,
        datatype=datatype,
    )


def _normalise_impedance_map(
    impedance: Mapping[str, Any],
    freq: Any,
) -> dict[str, Any]:
    """Return ``{station: Z}`` from a mapping of arrays or ``Z`` objects."""
    from ..z.z import Z

    out: dict[str, Any] = {}
    for station_id, value in dict(impedance).items():
        key = str(station_id)
        if isinstance(value, Z):
            out[key] = value
        elif isinstance(value, tuple) and len(value) == 2:
            arr, f = value
            out[key] = impedance_to_z(arr, f, station=key)
        else:
            if freq is None:
                raise ValueError(
                    f"freq is required to build Z for station {key!r}."
                )
            out[key] = impedance_to_z(value, freq, station=key)
    return out


def field_session_to_edifiles(
    session: Any,
    impedance: Mapping[str, Any],
    freq: Any = None,
    *,
    savepath: str | None = None,
    method: str | None = None,
    write: bool = False,
    acqby: str = "pycsamt.iot",
) -> dict[str, Any]:
    """Build per-station EDIs from a session plus impedance estimates.

    Station geometry (lat/lon/elevation) is taken from the session's
    registered stations, so the resulting EDIs carry the field-recorded
    location without any manual bookkeeping.

    Parameters
    ----------
    session : pycsamt.iot.session.FieldSession
        Source of station geometry and the observed method.
    impedance : mapping
        ``{station_id: value}`` where *value* is a :class:`~pycsamt.z.z.Z`,
        a ``(impedance_array, freq)`` pair, or an impedance array paired
        with the shared *freq* argument.
    freq : array-like, optional
        Shared frequency vector used for array values that do not carry
        their own frequencies.
    savepath : str, optional
        Directory for written files (used only when ``write=True``).
    method : str, optional
        Overrides the session method recorded in each EDI.
    write : bool
        When ``True`` the EDIs are written and the returned mapping holds
        file paths; otherwise it holds in-memory ``EDIFile`` objects.
    acqby : str
        ``ACQBY`` provenance tag.

    Returns
    -------
    dict
        ``{station_id: EDIFile}`` or ``{station_id: path}``.
    """
    z_by_station = _normalise_impedance_map(impedance, freq)
    stations = {s.station_id: s for s in session.to_sites()}
    survey_method = method or getattr(session, "method", None)

    out: dict[str, Any] = {}
    for station_id, z_obj in z_by_station.items():
        station = stations.get(station_id)
        lat = getattr(station, "lat", None) if station else None
        lon = getattr(station, "lon", None) if station else None
        elevation = getattr(station, "elevation", None) if station else None
        ed = build_edifile(
            z_obj,
            station=station_id,
            lat=lat,
            lon=lon,
            elevation=elevation,
            method=survey_method,
            acqby=acqby,
        )
        if write:
            out[station_id] = ed.write(
                new_edifn=station_id,
                savepath=savepath,
            )
        else:
            out[station_id] = ed
    return out


def z_to_site(
    z: Any,
    *,
    station: str,
    lat: float | None = None,
    lon: float | None = None,
    elevation: float | None = None,
    method: str | None = None,
) -> Any:
    """Return an EDI-backed :class:`~pycsamt.site.base.Site` for *z*.

    Unlike :meth:`FieldSession.to_sites`, which returns acquisition
    descriptors, this yields a real ``Site`` built on an ``EDIFile`` and
    therefore requires the optional geospatial stack (``pyproj``).
    """
    ed = build_edifile(
        z,
        station=station,
        lat=lat,
        lon=lon,
        elevation=elevation,
        method=method,
    )
    try:
        from ..site.base import Site
    except ImportError as exc:  # pragma: no cover - optional dependency
        raise ImportError(
            "z_to_site requires the pyCSAMT geospatial stack "
            "(e.g. pyproj). Install it, or use build_edifile/z_to_edi "
            "which have no geospatial dependency."
        ) from exc
    return Site(ed)


def _resolve_sites(
    source: Any,
    impedance: Mapping[str, Any] | None,
    freq: Any,
    *,
    method: str | None,
) -> Any:
    """Return a ``site.Sites`` collection from a session, sites, or EDIs."""
    try:
        from ..site.base import Sites
    except ImportError as exc:  # pragma: no cover - optional dependency
        raise ImportError(
            "emtools_qc requires the pyCSAMT geospatial stack (e.g. pyproj)."
        ) from exc

    if isinstance(source, Sites):
        return source
    if hasattr(source, "to_sites_collection"):
        if impedance is None:
            raise ValueError(
                "impedance is required when source is a FieldSession."
            )
        return source.to_sites_collection(impedance, freq, method=method)
    # Otherwise treat *source* as an EDI source (paths, dir, EDIFile, ...).
    edifiles = list(_iter_edifiles(source))
    if not edifiles:
        raise ValueError(
            "emtools_qc could not resolve any sites from the source; pass a "
            "FieldSession with impedance, a Sites collection, or EDI files."
        )
    return Sites(edifiles)


def emtools_qc(
    source: Any,
    impedance: Mapping[str, Any] | None = None,
    freq: Any = None,
    *,
    method: str | None = None,
    flags: bool = False,
    api: bool | None = None,
    **kwargs: Any,
) -> Any:
    """Route IoT-acquired impedance through :mod:`pycsamt.emtools` QC.

    Running the *same* coherence/skew/SNR quality control that downstream
    processing uses keeps field-side and post-processing QC consistent,
    rather than having two independent notions of a "good" station.

    Parameters
    ----------
    source : FieldSession, site.Sites, or EDI source
        A :class:`~pycsamt.iot.session.FieldSession` (with *impedance*), an
        already-built ``Sites`` collection, or any EDI source accepted by
        :func:`read_edi_survey` (directory, glob, file, ``EDIFile``, ...).
    impedance : mapping, optional
        Per-station impedance, required when *source* is a session (see
        :func:`field_session_to_edifiles` for the accepted forms).
    freq : array-like, optional
        Shared frequency vector for impedance arrays without their own.
    method : str, optional
        Acquisition method recorded on the built EDIs.
    flags : bool
        When ``True`` return :func:`pycsamt.emtools.qc.qc_flags` (per-station
        pass/flag verdicts); otherwise :func:`pycsamt.emtools.qc.build_qc_table`
        (the full QC metrics table).
    api : bool, optional
        Forwarded to the emtools helper's table wrapper.
    **kwargs
        Forwarded to the emtools helper (thresholds such as ``min_frac_ok``,
        ``min_snr_med``, ``max_skew_med``).

    Returns
    -------
    A pyCSAMT QC table (or flags table).

    Notes
    -----
    Requires the optional geospatial stack (``pyproj``), since emtools QC
    operates on EDI-backed sites.
    """
    from ..emtools.qc import build_qc_table, qc_flags

    sites = _resolve_sites(source, impedance, freq, method=method)
    if flags:
        return qc_flags(sites, **kwargs)
    return build_qc_table(sites, api=api, **kwargs)


# ---------------------------------------------------------------------------
# reverse: EDI survey -> acquisition planning
# ---------------------------------------------------------------------------
def _iter_edifiles(sources: Any) -> Iterable[Any]:
    """Yield ``EDIFile`` objects from paths, directories, globs, or objects."""
    from ..seg.edi import EDIFile

    if sources is None:
        return
    if isinstance(sources, (str, os.PathLike)):
        candidates: list[Any] = [sources]
    elif isinstance(sources, EDIFile):
        candidates = [sources]
    elif isinstance(sources, (Sequence, Iterable)):
        candidates = list(sources)
    else:
        candidates = [sources]

    for item in candidates:
        if isinstance(item, EDIFile):
            yield item
            continue
        path = os.fspath(item)
        if os.path.isdir(path):
            for edi_path in sorted(glob.glob(os.path.join(path, "*.edi"))):
                yield EDIFile(edi_path)
        elif any(ch in path for ch in "*?["):
            for edi_path in sorted(glob.glob(path)):
                yield EDIFile(edi_path)
        else:
            yield EDIFile(path)


def _clean_coord(value: Any, kind: str) -> float | None:
    """Return a finite coordinate, mapping NaN/None to ``None``."""
    if value is None:
        return None
    try:
        out = float(value)
    except (TypeError, ValueError):
        return None
    if not np.isfinite(out):
        return None
    if kind == "lat":
        return _c.as_latitude(out)
    if kind == "lon":
        return _c.as_longitude(out)
    return out


def _channels_from_edi(ed: Any) -> list[str]:
    """Return the acquisition channels implied by a station's EDI.

    A tensor EDI records the four-component AMT/MT layout, which is the
    right default for re-occupation planning. The ``ed`` argument is kept
    so this can specialise per station if a future EDI carries an explicit
    channel list.
    """
    return list(_DEFAULT_CHANNELS)


def _edi_record(ed: Any) -> dict[str, Any]:
    """Summarise one ``EDIFile`` into a plain dictionary."""
    head = ed.get_section("head")
    station = ed.station or getattr(head, "dataid", None) or "site"
    lat = _clean_coord(getattr(head, "lat", None), "lat") if head else None
    lon = _clean_coord(getattr(head, "lon", None), "lon") if head else None
    elev = _clean_coord(getattr(head, "elev", None), "elev") if head else None

    z = getattr(ed, "Z", None)
    freq = getattr(z, "freq", None) if z is not None else None
    if freq is not None and getattr(freq, "size", 0):
        f = np.asarray(freq, dtype=float)
        f = f[np.isfinite(f) & (f > 0)]
        f_min = float(f.min()) if f.size else None
        f_max = float(f.max()) if f.size else None
        n_freq = int(f.size)
    else:
        f_min = f_max = None
        n_freq = 0
    return dict(
        station=str(station),
        lat=lat,
        lon=lon,
        elevation=elev,
        n_freq=n_freq,
        f_min_hz=f_min,
        f_max_hz=f_max,
        channels=_channels_from_edi(ed),
    )


def read_edi_survey(sources: Any) -> list[dict[str, Any]]:
    """Read an EDI survey into per-station summary records.

    Parameters
    ----------
    sources : str, path, EDIFile, or iterable thereof
        A directory of ``.edi`` files, a glob pattern, a single file or
        ``EDIFile``, or any iterable mixing these.

    Returns
    -------
    list of dict
        One record per station with ``station``, ``lat``, ``lon``,
        ``elevation``, ``n_freq``, ``f_min_hz``, ``f_max_hz``, and
        ``channels``.
    """
    return [_edi_record(ed) for ed in _iter_edifiles(sources)]


def edi_survey_table(sources: Any, *, api: bool | None = None) -> Any:
    """Return an EDI-survey summary as a pyCSAMT table."""
    records = read_edi_survey(sources)
    rows = []
    for rec in records:
        row = dict(rec)
        row["channels"] = ";".join(rec["channels"])
        row["n_channels"] = len(rec["channels"])
        rows.append(row)
    df = pd.DataFrame.from_records(rows)
    return maybe_wrap_frame(
        df,
        api=api,
        name="iot_edi_survey_table",
        kind="iot.bridge.survey",
        source=sources,
        description="Station geometry and frequency coverage from an EDI survey.",
    )


def _sample_rate_hint(f_max: float | None) -> float | None:
    """Suggest an acquisition sample rate from the highest EDI frequency.

    A conventional acquisition margin of about five times the highest
    recovered frequency keeps that band comfortably below Nyquist.
    """
    if f_max is None or not np.isfinite(f_max) or f_max <= 0:
        return None
    return float(5.0 * f_max)


def field_session_from_edis(
    sources: Any,
    *,
    survey_id: str,
    protocol: str = "mqtt",
    role: str = "sensor_node",
    operator: str | None = None,
    method: str | None = None,
    device_suffix: str = "-node",
) -> Any:
    """Seed a :class:`~pycsamt.iot.session.FieldSession` from an EDI survey.

    Every station in *sources* becomes a registered station (with its
    recorded geometry and channels) and a sensor node, ready to be
    re-occupied in a follow-up IoT-enabled campaign.

    Parameters
    ----------
    sources : str, path, EDIFile, or iterable thereof
        The EDI survey to import (see :func:`read_edi_survey`).
    survey_id : str
        Identifier for the new session.
    protocol : str
        Telemetry protocol assigned to each seeded device.
    role : str
        Device role for each seeded node.
    operator : str, optional
        Operator recorded on the session.
    method : str, optional
        Acquisition method recorded on the session.
    device_suffix : str
        Suffix appended to a station id to form its device id.

    Returns
    -------
    pycsamt.iot.session.FieldSession
    """
    from .core import DeviceConfig
    from .session import FieldSession
    from .station import StationConfig

    records = read_edi_survey(sources)
    session = FieldSession(
        survey_id, operator=operator, method=method
    )
    for rec in records:
        station_id = rec["station"]
        device_id = f"{station_id}{device_suffix}"
        session.add_station(
            StationConfig(
                station_id=station_id,
                lat=rec["lat"],
                lon=rec["lon"],
                elevation=rec["elevation"],
                channels=list(rec["channels"]),
            )
        )
        session.add_device(
            DeviceConfig(
                device_id=device_id,
                station=station_id,
                protocol=protocol,
                sample_rate_hz=_sample_rate_hint(rec["f_max_hz"]),
                channels=list(rec["channels"]),
                role=role,
            )
        )
    return session


def deployment_from_edis(
    sources: Any,
    *,
    survey_id: str,
    protocol: str = "mqtt",
    role: str = "sensor_node",
    capabilities: Iterable[str] | None = None,
    device_suffix: str = "-node",
) -> Any:
    """Seed a :class:`~pycsamt.iot.core.DeploymentConfig` from an EDI survey.

    A lighter-weight counterpart to :func:`field_session_from_edis` that
    returns just the deployment inventory (one device per station), useful
    for planning device provisioning from a prior survey's station list.
    """
    from .core import DeploymentConfig, DeviceConfig

    records = read_edi_survey(sources)
    devices = [
        DeviceConfig(
            device_id=f"{rec['station']}{device_suffix}",
            station=rec["station"],
            protocol=protocol,
            sample_rate_hz=_sample_rate_hint(rec["f_max_hz"]),
            channels=list(rec["channels"]),
            role=role,
        )
        for rec in records
    ]
    return DeploymentConfig(
        survey_id=survey_id,
        devices=devices,
        capabilities=list(capabilities or []),
    )
