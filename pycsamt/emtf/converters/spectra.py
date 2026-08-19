# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

"""FCU-compatible recovery of EMTFs from historical EDI SPECTRA.

The SEG SPECTRA section stores a Hermitian cross-power matrix at each
frequency. EMTF FCU uses those matrices to recover the transfer functions and
the two covariance factors required for statistically correct rotation. This
module implements the same single-station and remote-reference equations while
keeping parsing in :mod:`pycsamt.seg.spectra` and rotation in
:mod:`pycsamt.emtf.orientation`.
"""

from __future__ import annotations

import warnings
from copy import deepcopy
from dataclasses import dataclass
from os import PathLike
from typing import TYPE_CHECKING, Any

import numpy as np

from ...exceptions import EdIDataError
from ...metadata import OrientationMeta, ProcessingMeta, RemoteReferenceMeta
from ..base import (
    IMPEDANCE_INPUT_CHANNELS,
    IMPEDANCE_OUTPUT_CHANNELS,
    TIPPER_INPUT_CHANNELS,
    TIPPER_OUTPUT_CHANNELS,
)
from ..document import EMTF
from ..estimates import StatisticalEstimate
from ..transfer import TransferFunction

if TYPE_CHECKING:
    from ...seg.spectra import Spectra

__all__ = [
    "SpectraCovarianceWarning",
    "SpectraRecoveryError",
    "SpectraChannelMap",
    "SpectraRecoveryResult",
    "resolve_spectra_channels",
    "recover_spectra_transfer_functions",
    "spectra_to_emtf",
]


class SpectraCovarianceWarning(UserWarning):
    """Warning emitted for recoverable SPECTRA covariance issues."""


class SpectraRecoveryError(EdIDataError):
    """Raised when an exact EDI SPECTRA recovery cannot be performed."""


@dataclass(frozen=True)
class SpectraChannelMap:
    """Resolved local, remote, and output channel indices."""

    hx: int
    hy: int
    ex: int
    ey: int
    hz: int | None = None
    rx: int | None = None
    ry: int | None = None
    channel_types: tuple[str, ...] = ()

    @property
    def local_h(self) -> tuple[int, int]:
        return self.hx, self.hy

    @property
    def remote_h(self) -> tuple[int, int]:
        return (
            self.hx if self.rx is None else self.rx,
            self.hy if self.ry is None else self.ry,
        )

    @property
    def outputs(self) -> tuple[int, ...]:
        if self.hz is None:
            return self.ex, self.ey
        return self.hz, self.ex, self.ey

    @property
    def reference_type(self) -> str:
        return (
            "single_station"
            if self.remote_h == self.local_h
            else "remote_reference"
        )


@dataclass(repr=False)
class SpectraRecoveryResult:
    """Numerical result of one FCU-compatible SPECTRA recovery."""

    frequency: np.ndarray
    periods: np.ndarray
    transfer_functions: dict[str, TransferFunction]
    channel_map: SpectraChannelMap
    used_indices: np.ndarray
    skipped_indices: np.ndarray
    avgt: np.ndarray
    rotspec: np.ndarray
    combined_residual_covariance: np.ndarray | None = None
    combined_output_channels: tuple[str, ...] = ()

    @property
    def impedance(self) -> TransferFunction:
        return self.transfer_functions["impedance"]

    @property
    def tipper(self) -> TransferFunction | None:
        return self.transfer_functions.get("tipper")


_VALID_NFREQ_POLICIES = frozenset({"raise", "warn", "ignore"})
_VALID_MISSING_POLICIES = frozenset({"raise", "skip"})
_VALID_AVGT_POLICIES = frozenset({"raise", "warn", "unit"})


def _norm_channel(value: Any) -> str:
    return "".join(ch for ch in str(value).upper() if ch.isalpha())


def _channel_types(spectra: Any) -> tuple[str, ...]:
    ids = list(getattr(spectra, "chan_ids", []) or [])
    mapping = dict(getattr(spectra, "id_to_chtype", {}) or {})
    if not ids:
        raise SpectraRecoveryError(
            "SPECTRA channel order is unavailable; SPECTRASECT measurement "
            "IDs are required for FCU-compatible recovery"
        )
    raw = [mapping.get(str(mid), mid) for mid in ids]
    kinds = tuple(_norm_channel(value) for value in raw)
    if len(kinds) != int(getattr(spectra, "n_chan", 0)):
        raise SpectraRecoveryError(
            "SPECTRA channel count does not match the measurement-ID list"
        )
    return kinds


def _indices(kinds: tuple[str, ...], labels: set[str]) -> list[int]:
    return [i for i, kind in enumerate(kinds) if kind in labels]


def resolve_spectra_channels(spectra: Any) -> SpectraChannelMap:
    """Resolve FCU local/remote channel roles from SPECTRA channel order.

    Historical files commonly encode remote H channels as a second ``HX/HY``
    pair. Explicit ``RX/RY`` or ``RHX/RHY`` labels are also accepted. When no
    remote pair is present, the local ``HX/HY`` pair is reused exactly as FCU
    does for a single-station analysis.
    """
    kinds = _channel_types(spectra)

    hx_all = _indices(kinds, {"HX"})
    hy_all = _indices(kinds, {"HY"})
    ex_all = _indices(kinds, {"EX"})
    ey_all = _indices(kinds, {"EY"})
    hz_all = _indices(kinds, {"HZ"})
    rx_explicit = _indices(kinds, {"RX", "RHX", "REMOTEHX"})
    ry_explicit = _indices(kinds, {"RY", "RHY", "REMOTEHY"})

    if not hx_all or not hy_all or not ex_all or not ey_all:
        raise SpectraRecoveryError(
            "SPECTRA recovery requires local HX, HY, EX, and EY channels; "
            f"resolved channel types are {kinds!r}"
        )

    hx = hx_all[0]
    hy = hy_all[0]
    ex = ex_all[0]
    ey = ey_all[0]
    hz = hz_all[0] if hz_all else None

    rx = rx_explicit[0] if rx_explicit else None
    ry = ry_explicit[0] if ry_explicit else None
    if rx is None and len(hx_all) > 1:
        rx = hx_all[1]
    if ry is None and len(hy_all) > 1:
        ry = hy_all[1]

    if (rx is None) != (ry is None):
        raise SpectraRecoveryError(
            "incomplete remote-reference pair: both RX/RHX and RY/RHY are "
            "required when either remote channel is present"
        )

    return SpectraChannelMap(
        hx=hx,
        hy=hy,
        ex=ex,
        ey=ey,
        hz=hz,
        rx=rx,
        ry=ry,
        channel_types=kinds,
    )


def _validate_policy(value: str, allowed: frozenset[str], name: str) -> str:
    value = str(value).strip().lower()
    if value not in allowed:
        raise ValueError(f"{name} must be one of {sorted(allowed)}")
    return value


def _resolve_avgt(
    spectra: Any,
    indices: np.ndarray,
    *,
    policy: str,
) -> np.ndarray:
    policy = _validate_policy(policy, _VALID_AVGT_POLICIES, "avgt_policy")
    values = np.asarray(getattr(spectra, "avgt", []), dtype=float)
    present = np.asarray(
        getattr(spectra, "avgt_present", np.ones(values.size, dtype=bool)),
        dtype=bool,
    )
    if values.size != int(getattr(spectra, "n_freq", 0)):
        raise SpectraRecoveryError("SPECTRA AVGT vector has invalid length")
    if present.size != values.size:
        present = np.ones(values.size, dtype=bool)

    out = values[indices].copy()
    supplied = present[indices]
    bad = (~supplied) | (~np.isfinite(out)) | (out <= 0.0)
    if not np.any(bad):
        return out

    where = indices[bad].tolist()
    message = (
        "FCU covariance recovery requires positive AVGT for every frequency; "
        f"missing/invalid AVGT at original index/indices {where}"
    )
    if policy == "raise":
        raise SpectraRecoveryError(message)
    if policy == "warn":
        warnings.warn(
            message + "; substituting AVGT=1 explicitly",
            SpectraCovarianceWarning,
            stacklevel=3,
        )
    out[bad] = 1.0
    return out


def _relevant_missing(
    missing: np.ndarray | None,
    channel_map: SpectraChannelMap,
) -> np.ndarray:
    if missing is None:
        return np.zeros(0, dtype=bool)
    relevant = sorted(
        set(channel_map.local_h)
        | set(channel_map.remote_h)
        | set(channel_map.outputs)
    )
    sub = missing[:, relevant][:, :, relevant]
    return np.any(sub, axis=(1, 2))


def _inv2(matrix: np.ndarray, *, index: int, label: str) -> np.ndarray:
    try:
        return np.linalg.inv(matrix)
    except np.linalg.LinAlgError as exc:
        raise SpectraRecoveryError(
            f"{label} is singular at SPECTRA frequency index {index}"
        ) from exc


def _fcu_period_recovery(
    cavg: np.ndarray,
    channel_map: SpectraChannelMap,
    *,
    avgt: float,
    index: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Apply equations implemented by FCU ``read_edi_spectra``."""
    h = channel_map.local_h
    r = channel_map.remote_h
    out = channel_map.outputs

    rh_h = cavg[np.ix_(r, h)]
    rh_r = cavg[np.ix_(r, r)]
    hh_h = cavg[np.ix_(h, h)]
    rh_e = cavg[np.ix_(r, out)]
    hh_e = cavg[np.ix_(h, out)]
    eh_e = cavg[np.ix_(out, out)]

    inv_rh_h = _inv2(rh_h, index=index, label="remote/local H block")
    zh = inv_rh_h @ rh_e
    tf = np.conjugate(zh.T)

    inv_rh_h_h = _inv2(
        np.conjugate(rh_h.T),
        index=index,
        label="conjugate remote/local H block",
    )
    inv_sig_cov = inv_rh_h @ rh_r @ inv_rh_h_h

    resid_cov = (
        eh_e
        - np.conjugate(zh.T) @ hh_e
        - np.conjugate(hh_e.T) @ zh
        + np.conjugate(zh.T) @ hh_h @ zh
    ) / float(avgt)

    variance = np.real(
        np.diag(resid_cov)[:, None] * np.diag(inv_sig_cov)[None, :]
    )
    return tf, variance, inv_sig_cov, resid_cov


def recover_spectra_transfer_functions(
    spectra: Any,
    *,
    nfreq_policy: str = "raise",
    missing_policy: str = "raise",
    avgt_policy: str = "raise",
) -> SpectraRecoveryResult:
    """Recover Z/T and full covariance from an EDI SPECTRA container.

    The implementation mirrors EMTF FCU's ``read_edi_spectra`` equations. It
    intentionally does not rotate the spectra or the recovered TFs; Phase 7's
    EMTF rotation engine performs any requested coordinate transformation.
    """
    nfreq_policy = _validate_policy(
        nfreq_policy,
        _VALID_NFREQ_POLICIES,
        "nfreq_policy",
    )
    missing_policy = _validate_policy(
        missing_policy,
        _VALID_MISSING_POLICIES,
        "missing_policy",
    )

    # Exact EMTF-FCU recovery depends on the packed SEG SPECTRA sign
    # convention. Historical ``Spectra.S`` uses pyCSAMT's conventional
    # cross-spectrum view and is also used by synthetic ``Spectra.from_Z``
    # objects. Do not silently reinterpret such objects as an EDI SPECTRA
    # archive: that can conjugate the recovered transfer function.
    if getattr(spectra, "_S_fcu", None) is None:
        raise SpectraRecoveryError(
            "exact FCU covariance recovery requires an EDI SPECTRA "
            "container parsed with its packed SEG convention; synthetic "
            "or legacy Spectra objects should continue to use to_Z()"
        )

    if hasattr(spectra, "validate_frequency_count"):
        spectra.validate_frequency_count(policy=nfreq_policy)

    n_freq = int(getattr(spectra, "n_freq", 0))
    if n_freq <= 0:
        raise SpectraRecoveryError("EDI SPECTRA contains no usable frequency")

    freq = np.asarray(getattr(spectra, "freq", []), dtype=float)
    if freq.size != n_freq or np.any(~np.isfinite(freq)) or np.any(freq <= 0):
        raise SpectraRecoveryError("SPECTRA frequency vector is invalid")

    channel_map = resolve_spectra_channels(spectra)
    cube = np.asarray(
        getattr(spectra, "fcu_cross_spectra", getattr(spectra, "S", None)),
        dtype=complex,
    )
    if cube.shape != (n_freq, int(spectra.n_chan), int(spectra.n_chan)):
        raise SpectraRecoveryError(
            "SPECTRA cross-power cube has invalid shape"
        )

    missing = getattr(spectra, "missing_mask", None)
    missing_freq = _relevant_missing(missing, channel_map)
    all_indices = np.arange(n_freq, dtype=int)
    if missing_freq.size == 0:
        missing_freq = np.zeros(n_freq, dtype=bool)
    if np.any(missing_freq) and missing_policy == "raise":
        bad = all_indices[missing_freq].tolist()
        raise SpectraRecoveryError(
            "exact covariance recovery cannot use frequencies with missing "
            f"required cross-spectra; bad index/indices: {bad}"
        )
    used = all_indices[~missing_freq]
    skipped = all_indices[missing_freq]
    if used.size == 0:
        raise SpectraRecoveryError("all SPECTRA frequencies are incomplete")
    if skipped.size:
        warnings.warn(
            "skipping incomplete SPECTRA frequencies at original indices "
            f"{skipped.tolist()}",
            SpectraCovarianceWarning,
            stacklevel=2,
        )

    avgt = _resolve_avgt(spectra, used, policy=avgt_policy)
    n_used = used.size
    n_out = 3 if channel_map.hz is not None else 2
    tf_all = np.empty((n_used, n_out, 2), dtype=complex)
    var_all = np.empty((n_used, n_out, 2), dtype=float)
    inv_all = np.empty((n_used, 2, 2), dtype=complex)
    resid_all = np.empty((n_used, n_out, n_out), dtype=complex)

    for out_index, source_index in enumerate(used):
        tf, var, inv_cov, resid_cov = _fcu_period_recovery(
            cube[source_index],
            channel_map,
            avgt=float(avgt[out_index]),
            index=int(source_index),
        )
        tf_all[out_index] = tf
        var_all[out_index] = var
        inv_all[out_index] = inv_cov
        resid_all[out_index] = resid_cov

    frequency = freq[used]
    periods = 1.0 / frequency
    rots = np.asarray(
        getattr(spectra, "rotspec", np.full(n_freq, np.nan)),
        dtype=float,
    )
    rots = rots[used] if rots.size == n_freq else np.full(n_used, np.nan)

    estimates_common = {
        "inverse_signal_covariance": StatisticalEstimate(
            name="INVSIGCOV",
            kind="inverse_signal_covariance",
            data=inv_all,
            attrs={"source_format": "edi_spectra", "full_covariance": True},
        )
    }

    if channel_map.hz is None:
        z_data = tf_all
        z_var = var_all
        z_resid = resid_all
        tipper = None
    else:
        z_data = tf_all[:, 1:3, :]
        z_var = var_all[:, 1:3, :]
        z_resid = resid_all[:, 1:3, 1:3]
        t_data = tf_all[:, 0:1, :]
        t_var = var_all[:, 0:1, :]
        t_resid = resid_all[:, 0:1, 0:1]
        tipper = TransferFunction(
            name="tipper",
            data=t_data,
            input_channels=TIPPER_INPUT_CHANNELS,
            output_channels=TIPPER_OUTPUT_CHANNELS,
            periods=periods,
            estimates={
                "variance": StatisticalEstimate(
                    name="VAR",
                    kind="variance",
                    data=t_var,
                    attrs={"source_format": "edi_spectra"},
                ),
                "inverse_signal_covariance": estimates_common[
                    "inverse_signal_covariance"
                ].copy(),
                "residual_covariance": StatisticalEstimate(
                    name="RESIDCOV",
                    kind="residual_covariance",
                    data=t_resid,
                    attrs={
                        "source_format": "edi_spectra",
                        "full_covariance": True,
                    },
                ),
            },
            attrs={
                "source_format": "edi_spectra",
                "reference_type": channel_map.reference_type,
                "edi_rotspec": rots.tolist(),
            },
        )

    impedance = TransferFunction(
        name="impedance",
        data=z_data,
        input_channels=IMPEDANCE_INPUT_CHANNELS,
        output_channels=IMPEDANCE_OUTPUT_CHANNELS,
        periods=periods,
        estimates={
            "variance": StatisticalEstimate(
                name="VAR",
                kind="variance",
                data=z_var,
                attrs={"source_format": "edi_spectra"},
            ),
            "inverse_signal_covariance": estimates_common[
                "inverse_signal_covariance"
            ],
            "residual_covariance": StatisticalEstimate(
                name="RESIDCOV",
                kind="residual_covariance",
                data=z_resid,
                attrs={
                    "source_format": "edi_spectra",
                    "full_covariance": True,
                },
            ),
        },
        attrs={
            "source_format": "edi_spectra",
            "reference_type": channel_map.reference_type,
            "edi_rotspec": rots.tolist(),
        },
    )

    transfer_functions = {"impedance": impedance}
    if tipper is not None:
        transfer_functions["tipper"] = tipper

    return SpectraRecoveryResult(
        frequency=frequency,
        periods=periods,
        transfer_functions=transfer_functions,
        channel_map=channel_map,
        used_indices=used,
        skipped_indices=skipped,
        avgt=avgt,
        rotspec=rots,
        combined_residual_covariance=resid_all.copy(),
        combined_output_channels=(
            ("Hz", "Ex", "Ey")
            if channel_map.hz is not None
            else ("Ex", "Ey")
        ),
    )


def _spectra_orientation(
    spectra: Any,
    result: SpectraRecoveryResult,
    edi: Any | None,
) -> OrientationMeta:
    present = np.asarray(
        getattr(
            spectra,
            "rotspec_present",
            np.isfinite(np.asarray(getattr(spectra, "rotspec", []))),
        ),
        dtype=bool,
    )
    if present.size == int(getattr(spectra, "n_freq", 0)):
        present = present[result.used_indices]
    else:
        present = np.isfinite(result.rotspec)

    extra: dict[str, Any] = {
        "edi_rotspec": result.rotspec.tolist(),
        "edi_zrot": result.rotspec.tolist(),
    }
    if result.tipper is not None:
        extra["edi_trot"] = result.rotspec.tolist()

    if present.size and np.all(present):
        coordsys = None
        if edi is not None:
            head = edi.get_section("head")
            coordsys = getattr(head, "coordsys", None) if head else None
        extra["edi_coordsys"] = coordsys
        if result.rotspec.size and np.allclose(
            result.rotspec,
            result.rotspec[0],
            equal_nan=False,
        ):
            angle = float(result.rotspec[0])
            if "geographic" in str(coordsys or "").lower():
                return OrientationMeta(
                    mode="orthogonal",
                    angle_to_geographic_north=angle,
                    rotation_info=(
                        "EDI SPECTRA ROTSPEC is present and constant; "
                        "HEAD.COORDSYS identifies a geographic reference."
                    ),
                    extra=extra,
                )
            return OrientationMeta(
                rotation_info=(
                    "EDI SPECTRA ROTSPEC is present and FCU treats the data "
                    "as orthogonal, but its angle was not promoted to a "
                    "geographic azimuth because HEAD.COORDSYS is ambiguous."
                ),
                extra=extra,
            )
        return OrientationMeta(
            rotation_info=(
                "EDI SPECTRA contains frequency-dependent ROTSPEC values. "
                "They are retained for Phase-7 opt-in rotation."
            ),
            extra=extra,
        )

    if not np.any(present):
        return OrientationMeta(
            mode="sitelayout",
            rotation_info=(
                "No ROTSPEC value was supplied; FCU-compatible interpretation "
                "therefore keeps the transfer functions in site layout."
            ),
            extra=extra,
        )
    return OrientationMeta(
        rotation_info=(
            "ROTSPEC is present for only part of the SPECTRA frequency set; "
            "orientation is retained as ambiguous."
        ),
        extra=extra,
    )


def _ensure_sources(
    source: Any,
    spectra: Any | None,
) -> tuple[Any | None, Any]:
    from ...seg.edi import EDIFile
    from ...seg.spectra import Spectra

    if spectra is not None and not isinstance(spectra, Spectra):
        raise TypeError("spectra must be a pycsamt.seg.spectra.Spectra")
    if isinstance(source, Spectra):
        return None, source
    if isinstance(source, EDIFile):
        edi = source
    elif isinstance(source, (str, PathLike)):
        edi = EDIFile(source)
    else:
        raise TypeError(
            "SPECTRA source must be a Spectra, EDIFile, or filesystem path"
        )
    spec = spectra if spectra is not None else edi.get_section("spectra")
    if spec is None or int(getattr(spec, "n_freq", 0)) <= 0:
        raise SpectraRecoveryError("EDI file does not contain usable SPECTRA")
    return edi, spec


def spectra_to_emtf(
    source: Any,
    *,
    spectra: Any | None = None,
    nfreq_policy: str = "raise",
    missing_policy: str = "raise",
    avgt_policy: str = "raise",
    target_angle: float | None = None,
) -> EMTF:
    """Convert EDI SPECTRA to the format-neutral EMTF scientific model.

    If ``target_angle`` is supplied, the recovered full-covariance TFs are
    passed to the Phase-7 rotation engine; no spectra-specific rotation code is
    used here.
    """
    edi, spec = _ensure_sources(source, spectra)
    result = recover_spectra_transfer_functions(
        spec,
        nfreq_policy=nfreq_policy,
        missing_policy=missing_policy,
        avgt_policy=avgt_policy,
    )

    if edi is None:
        document = EMTF(
            subtype="MT_TF",
            periods=result.periods,
            orientation=_spectra_orientation(spec, result, None),
            station=getattr(spec, "name", None),
            attrs={"source_format": "edi_spectra"},
        )
    else:
        # Reuse the mature Phase-6 EDI metadata mappings without routing the
        # scientific data through EDI impedance blocks.
        from .edi import (
            _edi_auxiliary_metadata,
            _edi_field_notes,
            _edi_processing,
            _edi_provenance,
            _edi_site,
            _edi_site_layout,
        )

        processing = deepcopy(_edi_processing(edi))
        if processing is None:
            processing = ProcessingMeta()
        reference = result.channel_map.reference_type
        processing.remote_reference = RemoteReferenceMeta(
            reference_type=(
                "Single Station"
                if reference == "single_station"
                else "Remote Reference"
            ),
            extra={
                "edi_rx_index": result.channel_map.remote_h[0],
                "edi_ry_index": result.channel_map.remote_h[1],
            },
        )
        metadata = _edi_auxiliary_metadata(edi)
        metadata["source_format"] = "edi_spectra"
        metadata["edi_spectra"] = {
            "declared_nfreq": getattr(spec, "declared_nfreq", None),
            "parsed_nfreq": getattr(spec, "parsed_nfreq", None),
            "used_nfreq": int(result.frequency.size),
            "skipped_indices": result.skipped_indices.tolist(),
            "channel_types": list(result.channel_map.channel_types),
            "reference_type": reference,
            "avgt": result.avgt.tolist(),
            "rotspec": result.rotspec.tolist(),
            "full_covariance_recovered": True,
            "covariance_algorithm": "EMTF_FCU_v4.1",
            "combined_output_channels": list(
                result.combined_output_channels
            ),
            "cross_data_type_residual_covariance_in_emtf": False,
        }
        document = EMTF(
            subtype="MT_TF",
            periods=result.periods,
            site=_edi_site(edi),
            site_layout=_edi_site_layout(edi),
            orientation=_spectra_orientation(spec, result, edi),
            processing=processing,
            provenance=_edi_provenance(edi),
            field_notes=_edi_field_notes(edi),
            metadata=metadata,
            attrs={"source_format": "edi_spectra"},
        )

    if result.tipper is not None:
        # EMTF XML deliberately separates tipper and impedance residual
        # covariance and therefore cannot archive the Hz/E cross terms. Keep
        # the original combined matrix on the in-memory document so Phase 8
        # recovery itself remains lossless until a serialization boundary.
        document.attrs[
            "spectra_combined_residual_covariance"
        ] = result.combined_residual_covariance.copy()
        document.attrs[
            "spectra_combined_output_channels"
        ] = result.combined_output_channels

    for key, tf in result.transfer_functions.items():
        document.add_transfer_function(tf, key=key, replace=True)
    document.tags = tuple(document.transfer_functions)

    if target_angle is not None:
        orientation = document.orientation
        rots = result.rotspec
        if orientation is not None and orientation.mode == "sitelayout":
            return document.rotate(float(target_angle))
        source_angles = {
            "impedance": rots,
            "tipper": rots,
        }
        return document.rotate(
            float(target_angle),
            source_angles=source_angles,
        )
    return document
