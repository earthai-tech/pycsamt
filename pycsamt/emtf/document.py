# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

"""Document-level electromagnetic transfer-function scientific container."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np

from ..core.base import MTBase, TFBundle
from ..metadata import (
    CopyrightInfo,
    LocationMeta,
    OrientationMeta,
    ProcessingMeta,
    ProvenanceMeta,
    SiteLayout,
    SiteMeta,
    TransferFunctionQuality,
)
from .base import (
    IMPEDANCE_INPUT_CHANNELS,
    IMPEDANCE_OUTPUT_CHANNELS,
    LEGACY_STANDARD_ERROR_KIND,
    LEGACY_STANDARD_ERROR_NAME,
    TIPPER_INPUT_CHANNELS,
    TIPPER_OUTPUT_CHANNELS,
)
from .estimates import StatisticalEstimate
from .transfer import TransferFunction
from .validation import normalize_periods

__all__ = ["EMTF"]


def _as_frequency(value: Any | None) -> np.ndarray | None:
    if value is None:
        return None
    arr = np.asarray(value, dtype=float)
    if arr.ndim == 0:
        arr = arr.reshape(1)
    if arr.ndim != 1:
        raise ValueError("frequency must be a 1-D array")
    if arr.size and (not np.all(np.isfinite(arr)) or np.any(arr <= 0.0)):
        raise ValueError("frequency must contain finite positive values")
    return arr


def _tipper_to_matrix(value: Any) -> np.ndarray:
    arr = np.asarray(value)
    if arr.ndim == 1 and arr.shape == (2,):
        return arr[None, None, :]
    if arr.ndim == 2 and arr.shape[1] == 2:
        return arr[:, None, :]
    if arr.ndim == 3 and arr.shape[1:] == (1, 2):
        return arr
    raise ValueError(
        "legacy tipper must have shape (2,), (n, 2), or (n, 1, 2)"
    )


def _legacy_scalar_or_tensor(
    name: str,
    value: Any,
    *,
    periods: np.ndarray | None,
) -> TransferFunction:
    arr = np.asarray(value)
    if arr.ndim == 3 and arr.shape[1:] == (2, 2):
        return TransferFunction(
            name=name,
            data=arr,
            input_channels=IMPEDANCE_INPUT_CHANNELS,
            output_channels=IMPEDANCE_OUTPUT_CHANNELS,
            periods=periods,
        )
    if arr.ndim == 2 and arr.shape == (2, 2):
        return TransferFunction(
            name=name,
            data=arr,
            input_channels=IMPEDANCE_INPUT_CHANNELS,
            output_channels=IMPEDANCE_OUTPUT_CHANNELS,
            periods=periods,
        )
    if arr.ndim <= 1:
        return TransferFunction(name=name, data=arr, periods=periods)
    raise ValueError(
        f"legacy {name} data must be scalar/frequency-vector or (n, 2, 2)"
    )


@dataclass(repr=False)
class EMTF(MTBase):
    """Format-neutral electromagnetic transfer-function document.

    ``EMTF`` is the scientific object that EDI and EMTF XML adapters will
    eventually populate. Phase 2 adds reusable metadata objects while keeping
    the model independent of XML parsing and :mod:`pycsamt.seg`.
    """

    product_id: str | None = None
    description: str | None = None
    subtype: str | None = None
    tags: tuple[str, ...] = field(default_factory=tuple)
    periods: Any | None = None
    transfer_functions: dict[str, TransferFunction] = field(
        default_factory=dict
    )

    # Format-neutral document metadata. These objects are independent of
    # EDI and XML syntax and can therefore be reused by future formats.
    provenance: ProvenanceMeta | None = None
    copyright: CopyrightInfo | None = None
    site: SiteMeta | None = None
    site_layout: SiteLayout | None = None
    orientation: OrientationMeta | None = None
    processing: ProcessingMeta | None = None
    quality: TransferFunctionQuality | None = None
    field_notes: dict[str, Any] = field(default_factory=dict)

    # Legacy TFBundle bridge. These fields remain available so existing
    # pyCSAMT callers do not need to understand the richer metadata model.
    station: str | None = None
    station_id: str | int | None = None
    lat: float | None = None
    lon: float | None = None
    elev: float | None = None
    azimuth: float | None = None

    # Unmodelled document metadata and application-specific attributes.
    metadata: dict[str, Any] = field(default_factory=dict)
    attrs: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.validate()

    def validate(self) -> None:
        """Normalize document state and validate attached TF period grids."""
        self.periods = normalize_periods(self.periods)
        self.tags = tuple(
            str(tag).strip().lower()
            for tag in self.tags
            if str(tag).strip()
        )
        self.metadata = dict(self.metadata or {})
        self.attrs = dict(self.attrs or {})
        self.field_notes = dict(self.field_notes or {})
        self._validate_metadata_types()
        self._synchronize_site_bridge()

        incoming = dict(self.transfer_functions or {})
        self.transfer_functions = {}
        for key, tf in incoming.items():
            self.add_transfer_function(tf, key=key, replace=True)

    def _validate_metadata_types(self) -> None:
        expected = (
            ("provenance", self.provenance, ProvenanceMeta),
            ("copyright", self.copyright, CopyrightInfo),
            ("site", self.site, SiteMeta),
            ("site_layout", self.site_layout, SiteLayout),
            ("orientation", self.orientation, OrientationMeta),
            ("processing", self.processing, ProcessingMeta),
            ("quality", self.quality, TransferFunctionQuality),
        )
        for name, value, cls in expected:
            if value is not None and not isinstance(value, cls):
                raise TypeError(f"{name} must be {cls.__name__} or None")

    def _synchronize_site_bridge(self) -> None:
        """Expose explicit ``SiteMeta`` values through legacy aliases.

        Phase 2 intentionally does not infer a ``SiteMeta`` object from the
        legacy TFBundle fields. Their historical ``station``/``station_id``
        semantics are not identical to EMTF ``Site/Id`` and ``Site/Name``.
        The later EDI interoperability layer will perform that mapping with
        format-specific knowledge.
        """
        if self.site is None:
            return

        station_key = (
            self.site.site_id
            if self.site.site_id is not None
            else self.site.name
        )
        self._merge_bridge_value("station", station_key)
        self._merge_bridge_value("station_id", self.site.site_id)
        location = self.site.location
        if location is not None:
            self._merge_bridge_value("lat", location.latitude)
            self._merge_bridge_value("lon", location.longitude)
            self._merge_bridge_value("elev", location.elevation)

    def _merge_bridge_value(self, name: str, metadata_value: Any) -> None:
        if metadata_value is None:
            return
        legacy_value = getattr(self, name)
        if legacy_value is None:
            setattr(self, name, metadata_value)
            return
        if name in {"lat", "lon", "elev"}:
            if not np.isclose(float(legacy_value), float(metadata_value)):
                raise ValueError(
                    f"legacy {name} conflicts with SiteMeta: "
                    f"{legacy_value!r} != {metadata_value!r}"
                )
            return
        if str(legacy_value) != str(metadata_value):
            raise ValueError(
                f"legacy {name} conflicts with SiteMeta: "
                f"{legacy_value!r} != {metadata_value!r}"
            )

    @property
    def frequency(self) -> np.ndarray | None:
        """Return the document frequency vector in Hz when available."""
        if self.periods is not None:
            return 1.0 / np.asarray(self.periods, dtype=float)
        for tf in self.transfer_functions.values():
            if tf.periods is not None:
                return tf.frequency
        return None

    @property
    def n_periods(self) -> int:
        """Return the common number of period samples when known."""
        if self.periods is not None:
            return int(np.asarray(self.periods).size)
        for tf in self.transfer_functions.values():
            return tf.n_periods
        return 0

    def add_transfer_function(
        self,
        tf: TransferFunction,
        *,
        key: str | None = None,
        replace: bool = False,
    ) -> "EMTF":
        """Attach one scientific transfer function to the document."""
        if not isinstance(tf, TransferFunction):
            raise TypeError("tf must be a TransferFunction")
        tf_key = str(key or tf.name).strip().lower()
        if not tf_key:
            raise ValueError("transfer-function key must be non-empty")
        if tf_key in self.transfer_functions and not replace:
            raise ValueError(f"transfer function already exists: {tf_key}")

        if self.periods is not None and tf.periods is not None:
            doc_periods = np.asarray(self.periods, dtype=float)
            tf_periods = np.asarray(tf.periods, dtype=float)
            if doc_periods.shape != tf_periods.shape or not np.allclose(
                doc_periods,
                tf_periods,
                rtol=1.0e-10,
                atol=0.0,
            ):
                raise ValueError(
                    f"period grid for {tf_key!r} differs from EMTF document"
                )
        elif self.periods is None and tf.periods is not None:
            self.periods = np.array(tf.periods, copy=True)

        self.transfer_functions[tf_key] = tf
        if tf_key not in self.tags:
            self.tags = (*self.tags, tf_key)
        return self

    def get_transfer_function(self, key: str) -> TransferFunction | None:
        """Return a TF by semantic key or registered short code."""
        raw = str(key).strip()
        if not raw:
            return None
        direct = self.transfer_functions.get(raw.lower())
        if direct is not None:
            return direct
        for tf in self.transfer_functions.values():
            definition = tf.definition
            if definition is not None:
                names = {definition.name, *definition.aliases}
                if raw.upper() in names:
                    return tf
        return None

    @property
    def impedance(self) -> TransferFunction | None:
        """Return the matrix-oriented impedance TF."""
        return self.get_transfer_function("impedance")

    @property
    def tipper_tf(self) -> TransferFunction | None:
        """Return the matrix-oriented tipper TF."""
        return self.get_transfer_function("tipper")

    @property
    def z(self) -> np.ndarray | None:
        """Return the impedance array for legacy-style access."""
        tf = self.impedance
        return None if tf is None else tf.data

    @property
    def z_err(self) -> np.ndarray | None:
        """Return legacy standard errors when explicitly available."""
        tf = self.impedance
        if tf is None:
            return None
        estimate = tf.get_estimate(LEGACY_STANDARD_ERROR_KIND)
        return None if estimate is None else estimate.data

    @property
    def tipper(self) -> np.ndarray | None:
        """Return tipper data in the pyCSAMT ``(n, 1, 2)`` shape."""
        tf = self.tipper_tf
        return None if tf is None else tf.data

    @property
    def tipper_err(self) -> np.ndarray | None:
        """Return legacy tipper standard errors when available."""
        tf = self.tipper_tf
        if tf is None:
            return None
        estimate = tf.get_estimate(LEGACY_STANDARD_ERROR_KIND)
        return None if estimate is None else estimate.data

    @property
    def rho(self) -> np.ndarray | None:
        """Return stored apparent resistivity, or derive it from ``Z``."""
        tf = self.get_transfer_function("apparent_resistivity")
        if tf is not None:
            data = tf.data
            if data.shape[1:] == (1, 1):
                return data[:, 0, 0]
            return data
        z_obj = self.Z
        return None if z_obj is None else z_obj.resistivity

    @property
    def phase(self) -> np.ndarray | None:
        """Return stored impedance phase, or derive it from ``Z``."""
        tf = self.get_transfer_function("impedance_phase")
        if tf is not None:
            data = tf.data
            if data.shape[1:] == (1, 1):
                return data[:, 0, 0]
            return data
        z_obj = self.Z
        return None if z_obj is None else z_obj.phase

    @property
    def Z(self):
        """Build the existing :class:`pycsamt.z.Z` compatibility object."""
        if self.z is None:
            return None
        if self.z.shape[1:] != (2, 2):
            return None
        from ..z.z import Z

        return Z(
            z_array=self.z,
            z_err_array=self.z_err,
            freq=self.frequency,
            name=self.station,
        )

    @property
    def Tip(self):
        """Build the existing :class:`pycsamt.z.tipper.Tipper` object."""
        if self.tipper is None:
            return None
        from ..z.tipper import Tipper

        return Tipper(
            tipper_array=self.tipper,
            tipper_err_array=self.tipper_err,
            freq=self.frequency,
            name=self.station,
        )

    @classmethod
    def from_xml(cls, source, *, strict: bool = True) -> "EMTF":
        """Read an EMTF XML document through the XML adapter.

        The import is intentionally local so constructing the scientific core
        never requires XML serialization support unless this method is called.
        """
        from .xml import EMTFXMLReader

        return EMTFXMLReader(strict=strict).read(source)

    @classmethod
    def from_edi(cls, source, **kwargs) -> "EMTF":
        """Convert a historical SEG EDI object or path into ``EMTF``.

        EDI SPECTRA are preferred by default when present because they retain
        the full information needed to recover transfer-function covariance.
        Pass ``prefer_spectra=False`` to force the traditional impedance/tipper
        blocks instead. The adapter import remains local so the scientific core
        is independent of :mod:`pycsamt.seg` until requested.
        """
        from .converters.edi import edi_to_emtf

        return edi_to_emtf(source, **kwargs)

    @classmethod
    def from_edi_spectra(cls, source, **kwargs) -> "EMTF":
        """Recover EMTFs and full covariance directly from EDI SPECTRA.

        This is the explicit Phase-8 entry point for single-station and
        remote-reference cross-power spectra. Any requested rotation is
        delegated to the format-neutral Phase-7 rotation engine.
        """
        from .converters.spectra import spectra_to_emtf

        return spectra_to_emtf(source, **kwargs)

    def to_edi(self, *, on_loss: str = "warn"):
        """Return an in-memory :class:`pycsamt.seg.EDIFile` representation.

        Parameters
        ----------
        on_loss : {"warn", "raise", "ignore"}
            Policy for EMTF content that standard EDI cannot preserve.
        """
        from .converters.edi import emtf_to_edi

        return emtf_to_edi(self, on_loss=on_loss)

    def to_xml(
        self,
        *,
        strict: bool = True,
        precision: int = 17,
        pretty: bool = True,
        xml_declaration: bool = True,
    ) -> str:
        """Serialize this scientific document to an EMTF XML string."""
        from .xml import EMTFXMLWriter

        return EMTFXMLWriter(
            strict=strict,
            precision=precision,
            pretty=pretty,
        ).dumps(self, xml_declaration=xml_declaration)

    def write_xml(
        self,
        target,
        *,
        strict: bool = True,
        precision: int = 17,
        pretty: bool = True,
        xml_declaration: bool = True,
        encoding: str = "utf-8",
    ):
        """Write this document in EMTF XML format."""
        from .xml import write_emtf_xml

        return write_emtf_xml(
            self,
            target,
            strict=strict,
            precision=precision,
            pretty=pretty,
            xml_declaration=xml_declaration,
            encoding=encoding,
        )

    def write(self, target, *, format: str = "emtf_xml", **kwargs):
        """Write this document using a supported serialization format.

        Phase 6 supports EMTF XML and historical SEG EDI. EDI conversion is
        explicit and may emit :class:`DataLossWarning` for content that the
        historical format cannot represent.
        """
        normalized = str(format).strip().lower().replace("-", "_")
        if normalized in {"emtf_xml", "xml", "emtf"}:
            return self.write_xml(target, **kwargs)
        if normalized in {"edi", "seg_edi"}:
            from .converters.edi import write_edi

            return write_edi(self, target, **kwargs)
        raise ValueError("unsupported EMTF output format: " f"{format!r}")

    def rotate(
        self,
        angle: float | None = 0.0,
        *,
        target: str = "orthogonal",
        inplace: bool = False,
        source_angles=None,
        use_legacy_edi_rotation: bool = False,
        variance_policy: str = "drop",
        unsupported_estimates: str = "drop",
        derived_policy: str = "drop",
    ) -> "EMTF":
        """Rotate transfer functions with the format-neutral EMTF engine.

        The original :class:`~pycsamt.metadata.SiteLayout` remains physical
        acquisition metadata and is never rotated. Full covariance factors
        are transformed consistently when present.
        """
        from .orientation import rotate_emtf

        return rotate_emtf(
            self,
            angle,
            target=target,
            inplace=inplace,
            source_angles=source_angles,
            use_legacy_edi_rotation=use_legacy_edi_rotation,
            variance_policy=variance_policy,
            unsupported_estimates=unsupported_estimates,
            derived_policy=derived_policy,
        )

    def is_empty(self) -> bool:
        """Return ``True`` when no transfer-function content is attached."""
        return not bool(self.transfer_functions)

    def to_bundle(self) -> TFBundle:
        """Return a backward-compatible :class:`TFBundle` view."""
        tip = self.tipper
        tip_err = self.tipper_err

        def _stored_derived(tag: str):
            tf = self.get_transfer_function(tag)
            if tf is None:
                return None
            data = tf.data
            if data.shape[1:] == (1, 1):
                return data[:, 0, 0]
            return data

        rho = _stored_derived("apparent_resistivity")
        phase = _stored_derived("impedance_phase")
        if tip is not None and tip.ndim == 3 and tip.shape[1:] == (1, 2):
            tip = tip[:, 0, :]
        if (
            tip_err is not None
            and tip_err.ndim == 3
            and tip_err.shape[1:] == (1, 2)
        ):
            tip_err = tip_err[:, 0, :]

        return TFBundle(
            freq=self.frequency,
            z=self.z,
            z_err=self.z_err,
            tipper=tip,
            tipper_err=tip_err,
            rho=rho,
            phase=phase,
            station=self.station,
            station_id=self.station_id,
            lat=self.lat,
            lon=self.lon,
            elev=self.elev,
            azimuth=self.azimuth,
            attrs=dict(self.attrs),
            transfer_functions=dict(self.transfer_functions),
            estimates={
                key: dict(tf.estimates)
                for key, tf in self.transfer_functions.items()
                if tf.estimates
            },
        )

    @classmethod
    def from_bundle(cls, bundle: TFBundle) -> "EMTF":
        """Build an EMTF document from a legacy/generalized bundle."""
        if not isinstance(bundle, TFBundle):
            raise TypeError("bundle must be a TFBundle")

        freq = _as_frequency(bundle.freq)
        periods = None if freq is None else 1.0 / freq
        out = cls(
            periods=periods,
            station=bundle.station,
            station_id=bundle.station_id,
            lat=bundle.lat,
            lon=bundle.lon,
            elev=bundle.elev,
            azimuth=bundle.azimuth,
            attrs=dict(bundle.attrs or {}),
        )

        for key, value in dict(bundle.transfer_functions or {}).items():
            if isinstance(value, TransferFunction):
                out.add_transfer_function(value.copy(), key=key, replace=True)

        if bundle.z is not None and out.impedance is None:
            ztf = TransferFunction(
                name="impedance",
                data=np.asarray(bundle.z),
                input_channels=IMPEDANCE_INPUT_CHANNELS,
                output_channels=IMPEDANCE_OUTPUT_CHANNELS,
                periods=periods,
            )
            if bundle.z_err is not None:
                ztf.add_estimate(
                    StatisticalEstimate(
                        name=LEGACY_STANDARD_ERROR_NAME,
                        kind=LEGACY_STANDARD_ERROR_KIND,
                        data=np.asarray(bundle.z_err),
                        attrs={
                            "semantics": "pycsamt_legacy_z_err",
                            "source": "TFBundle",
                        },
                    )
                )
            out.add_transfer_function(ztf, replace=True)

        if bundle.tipper is not None and out.tipper_tf is None:
            tip_data = _tipper_to_matrix(bundle.tipper)
            ttf = TransferFunction(
                name="tipper",
                data=tip_data,
                input_channels=TIPPER_INPUT_CHANNELS,
                output_channels=TIPPER_OUTPUT_CHANNELS,
                periods=periods,
            )
            if bundle.tipper_err is not None:
                tip_err = _tipper_to_matrix(bundle.tipper_err)
                ttf.add_estimate(
                    StatisticalEstimate(
                        name=LEGACY_STANDARD_ERROR_NAME,
                        kind=LEGACY_STANDARD_ERROR_KIND,
                        data=tip_err,
                        attrs={
                            "semantics": "pycsamt_legacy_z_err",
                            "source": "TFBundle",
                        },
                    )
                )
            out.add_transfer_function(ttf, replace=True)

        if (
            bundle.rho is not None
            and out.get_transfer_function("apparent_resistivity") is None
        ):
            out.add_transfer_function(
                _legacy_scalar_or_tensor(
                    "apparent_resistivity",
                    bundle.rho,
                    periods=periods,
                ),
                replace=True,
            )
        if (
            bundle.phase is not None
            and out.get_transfer_function("impedance_phase") is None
        ):
            out.add_transfer_function(
                _legacy_scalar_or_tensor(
                    "impedance_phase",
                    bundle.phase,
                    periods=periods,
                ),
                replace=True,
            )

        return out
