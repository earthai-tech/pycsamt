# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

"""Format-neutral EMTF rotation and covariance transformations.

The implementation follows the general matrix formulation used by EMTF FCU:

``TF' = V @ TF @ U.T``

where ``U`` transforms input-channel coordinates and ``V`` transforms
output-channel coordinates.  When full error covariance factors are present,
``S`` (inverse signal covariance) and ``N`` (residual covariance) transform as
``S' = U @ S @ U.T`` and ``N' = V @ N @ V.T``.

The physical :class:`~pycsamt.metadata.SiteLayout` is never rotated here.  It
continues to describe the original field geometry while
:class:`~pycsamt.metadata.OrientationMeta` records the coordinate frame of the
transfer-function data.
"""

from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass
from math import cos, radians, sin
from typing import Any, Mapping
import warnings

import numpy as np

from ..metadata import OrientationMeta, SiteLayout
from .document import EMTF
from .transfer import TransferFunction

__all__ = [
    "ApproximateVarianceRotationWarning",
    "DerivedDataRotationWarning",
    "EMTFRotationError",
    "EMTFRotationWarning",
    "LegacyRotationAssumptionWarning",
    "RotationMatrices",
    "UnsupportedEstimateRotationWarning",
    "horizontal_rotation_matrix",
    "horizontal_inverse_rotation_matrix",
    "rotate_covariance",
    "rotate_emtf",
    "rotate_transfer_function",
]


class EMTFRotationError(ValueError):
    """Raised when an EMTF rotation cannot be defined unambiguously."""


class EMTFRotationWarning(UserWarning):
    """Base warning for scientifically incomplete EMTF rotations."""


class ApproximateVarianceRotationWarning(EMTFRotationWarning):
    """Warn that variance was handled without complete covariance."""


class UnsupportedEstimateRotationWarning(EMTFRotationWarning):
    """Warn that a statistical estimate cannot be rotated safely."""


class DerivedDataRotationWarning(EMTFRotationWarning):
    """Warn that derived products were dropped or retained stale."""


class LegacyRotationAssumptionWarning(EMTFRotationWarning):
    """Warn that historical EDI rotation metadata are being interpreted."""


@dataclass(frozen=True)
class RotationMatrices:
    """Per-period input and output matrices used in an EMTF rotation.

    Parameters
    ----------
    input_matrix : ndarray
        Array of shape ``(n_period, n_input, n_input)`` containing ``U``.
    output_matrix : ndarray
        Array of shape ``(n_period, n_output, n_output)`` containing ``V``.
    source_mode, target_mode : str
        Coordinate-frame modes used to construct the matrices.
    source_angles : ndarray or None
        Orthogonal source azimuth(s), in degrees clockwise from geographic
        north.  ``None`` when the source follows the physical site layout.
    target_angle : float or None
        Orthogonal target azimuth.  ``None`` for site-layout targets.
    """

    input_matrix: np.ndarray
    output_matrix: np.ndarray
    source_mode: str
    target_mode: str
    source_angles: np.ndarray | None
    target_angle: float | None

    @property
    def n_periods(self) -> int:
        return int(self.input_matrix.shape[0])

    @property
    def is_identity(self) -> bool:
        """Return whether both matrix families are identities."""
        ni = self.input_matrix.shape[1]
        no = self.output_matrix.shape[1]
        return bool(
            np.allclose(
                self.input_matrix,
                np.eye(ni)[None, :, :],
                rtol=0.0,
                atol=1.0e-12,
            )
            and np.allclose(
                self.output_matrix,
                np.eye(no)[None, :, :],
                rtol=0.0,
                atol=1.0e-12,
            )
        )


def _angle(value: Any, *, name: str) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise EMTFRotationError(f"{name} must be a finite angle") from exc
    if not np.isfinite(number):
        raise EMTFRotationError(f"{name} must be a finite angle")
    return number


def _angle_vector(value: Any, n: int, *, name: str) -> np.ndarray:
    arr = np.asarray(value, dtype=float).reshape(-1)
    if arr.size == 1:
        arr = np.full(n, float(arr[0]), dtype=float)
    if arr.size != n:
        raise EMTFRotationError(
            f"{name} must be scalar or length {n}; got {arr.size}"
        )
    if not np.all(np.isfinite(arr)):
        raise EMTFRotationError(f"{name} must contain finite angles")
    return arr


def horizontal_rotation_matrix(
    theta1: float,
    theta2: float,
    target_angle: float,
) -> np.ndarray:
    """Return the FCU ``Q`` transform from two channels to an orthogonal frame.

    Angles are degrees clockwise from geographic north.  The source channels
    need not be orthogonal.  The matrix maps source vector components into a
    right-handed orthogonal coordinate frame whose x-axis has azimuth
    ``target_angle``.
    """
    target_angle = _angle(target_angle, name="target_angle")
    t1 = radians(_angle(theta1, name="theta1") - target_angle)
    t2 = radians(_angle(theta2, name="theta2") - target_angle)
    return np.array(
        [
            [cos(t1), cos(t2)],
            [sin(t1), sin(t2)],
        ],
        dtype=float,
    )


def horizontal_inverse_rotation_matrix(
    theta1: float,
    theta2: float,
    target_angle: float,
) -> np.ndarray:
    """Return the inverse of :func:`horizontal_rotation_matrix`.

    The explicit formula mirrors EMTF FCU ``rot2inv`` and permits a
    non-orthogonal original site layout.  Parallel channels are singular and
    therefore rejected.
    """
    theta1 = _angle(theta1, name="theta1")
    theta2 = _angle(theta2, name="theta2")
    target_angle = _angle(target_angle, name="target_angle")
    det = sin(radians(theta2 - theta1))
    if abs(det) <= 1.0e-12:
        raise EMTFRotationError(
            "horizontal channel orientations are parallel or nearly parallel"
        )
    return np.array(
        [
            [
                sin(radians(theta2 - target_angle)) / det,
                -cos(radians(theta2 - target_angle)) / det,
            ],
            [
                -sin(radians(theta1 - target_angle)) / det,
                cos(radians(theta1 - target_angle)) / det,
            ],
        ],
        dtype=float,
    )


def _is_missing(data: np.ndarray) -> np.ndarray:
    if np.iscomplexobj(data):
        return ~np.isfinite(data.real) | ~np.isfinite(data.imag)
    return ~np.isfinite(data)


def _masked_bilinear(
    left: np.ndarray,
    data: np.ndarray,
    right: np.ndarray,
    *,
    coefficient_tol: float = 1.0e-13,
) -> np.ndarray:
    """Compute ``left @ data @ right.T`` while respecting missing values.

    A missing source entry contaminates only output entries that actually
    depend on it.  This avoids the usual ``0 * NaN -> NaN`` propagation during
    identity or right-angle rotations while still refusing to manufacture a
    rotated value when a required component is absent.
    """
    arr = np.asarray(data)
    missing = _is_missing(arr)
    filled = np.where(missing, 0.0, arr)
    result = left @ filled @ right.T
    if not np.any(missing):
        return result

    left_dep = np.abs(left) > coefficient_tol
    right_dep = np.abs(right) > coefficient_tol
    dependencies = np.einsum(
        "ia,ab,jb->ij",
        left_dep.astype(int),
        missing.astype(int),
        right_dep.astype(int),
    )
    required_missing = dependencies > 0
    if np.iscomplexobj(result):
        result = np.asarray(result, dtype=complex)
        result[required_missing] = np.nan + 1j * np.nan
    else:
        result = np.asarray(result, dtype=float)
        result[required_missing] = np.nan
    return result


def _channel_pair_indices(names: tuple[str, ...]) -> list[tuple[int, int]]:
    """Identify horizontal x/y output pairs while leaving z channels fixed."""
    if len(names) == 1:
        return []
    if len(names) == 2:
        return [(0, 1)]

    groups: dict[str, dict[str, int]] = {}
    for index, raw_name in enumerate(names):
        name = str(raw_name).strip().lower()
        if not name or name[-1:] not in {"x", "y"}:
            continue
        groups.setdefault(name[:-1], {})[name[-1]] = index

    pairs: list[tuple[int, int]] = []
    paired: set[int] = set()
    for group in groups.values():
        if "x" in group and "y" in group:
            ix, iy = group["x"], group["y"]
            pairs.append((ix, iy))
            paired.update((ix, iy))

    unpaired_horizontal = []
    for index, raw_name in enumerate(names):
        name = str(raw_name).strip().lower()
        if name[-1:] in {"x", "y"} and index not in paired:
            unpaired_horizontal.append(raw_name)
    if unpaired_horizontal:
        raise EMTFRotationError(
            "unable to identify complete x/y output pairs for channels: "
            + ", ".join(map(str, unpaired_horizontal))
        )
    return pairs


def _layout_channels(
    layout: SiteLayout,
    names: tuple[str, ...],
    *,
    role: str,
):
    channels = []
    for name in names:
        channel = layout.get_channel(name, role=role)
        if channel is None:
            raise EMTFRotationError(
                f"site layout has no {role} channel matching {name!r}"
            )
        channels.append(channel)
    return channels


def _require_horizontal_pair(channels, *, label: str) -> tuple[float, float]:
    if len(channels) != 2:
        raise EMTFRotationError(f"{label} must contain exactly two channels")
    values = []
    for channel in channels:
        if channel.orientation is None:
            raise EMTFRotationError(
                f"{label} channel {channel.name!r} has no orientation"
            )
        if channel.tilt is not None and abs(float(channel.tilt)) > 1.0e-8:
            raise EMTFRotationError(
                f"tilted {label} channel {channel.name!r} is not supported"
            )
        values.append(float(channel.orientation))
    return values[0], values[1]


def _site_input_transform(
    tf: TransferFunction,
    layout: SiteLayout,
    target_angle: float,
    *,
    inverse: bool,
) -> np.ndarray:
    channels = _layout_channels(
        layout,
        tf.input_channels,
        role="input",
    )
    theta1, theta2 = _require_horizontal_pair(
        channels,
        label="input",
    )
    if inverse:
        return horizontal_rotation_matrix(theta1, theta2, target_angle).T
    return horizontal_inverse_rotation_matrix(
        theta1,
        theta2,
        target_angle,
    ).T


def _site_output_transform(
    tf: TransferFunction,
    layout: SiteLayout,
    target_angle: float,
    *,
    inverse: bool,
) -> np.ndarray:
    matrix = np.eye(tf.n_output, dtype=float)
    channels = _layout_channels(
        layout,
        tf.output_channels,
        role="output",
    )
    for ix, iy in _channel_pair_indices(tf.output_channels):
        pair = [channels[ix], channels[iy]]
        theta1, theta2 = _require_horizontal_pair(pair, label="output")
        if inverse:
            block = horizontal_inverse_rotation_matrix(
                theta1,
                theta2,
                target_angle,
            )
        else:
            block = horizontal_rotation_matrix(
                theta1,
                theta2,
                target_angle,
            )
        indices = np.ix_([ix, iy], [ix, iy])
        matrix[indices] = block
    return matrix


def _orthogonal_input_transform(
    old_angle: float,
    new_angle: float,
) -> np.ndarray:
    # For an orthogonal source frame, FCU U == Q because Q^-1 == Q.T.
    return horizontal_rotation_matrix(
        old_angle,
        old_angle + 90.0,
        new_angle,
    )


def _orthogonal_output_transform(
    tf: TransferFunction,
    old_angle: float,
    new_angle: float,
) -> np.ndarray:
    matrix = np.eye(tf.n_output, dtype=float)
    block = horizontal_rotation_matrix(
        old_angle,
        old_angle + 90.0,
        new_angle,
    )
    for ix, iy in _channel_pair_indices(tf.output_channels):
        matrix[np.ix_([ix, iy], [ix, iy])] = block
    return matrix


def _rotation_matrices(
    tf: TransferFunction,
    *,
    source_mode: str,
    target_mode: str,
    site_layout: SiteLayout | None,
    source_angles: Any | None,
    target_angle: float | None,
) -> RotationMatrices:
    source_mode = str(source_mode).strip().lower()
    target_mode = str(target_mode).strip().lower()
    if source_mode not in {"orthogonal", "sitelayout"}:
        raise EMTFRotationError(
            "source_mode must be 'orthogonal' or 'sitelayout'"
        )
    if target_mode not in {"orthogonal", "sitelayout"}:
        raise EMTFRotationError(
            "target_mode must be 'orthogonal' or 'sitelayout'"
        )
    if tf.n_input != 2:
        raise EMTFRotationError(
            "Phase 7 EMTF rotation currently requires exactly two input "
            "channels"
        )

    n = tf.n_periods
    source_vector = None
    if source_mode == "orthogonal":
        if source_angles is None:
            raise EMTFRotationError(
                "orthogonal source data require source angle metadata"
            )
        source_vector = _angle_vector(
            source_angles,
            n,
            name="source_angles",
        )

    if target_mode == "orthogonal":
        if target_angle is None:
            raise EMTFRotationError(
                "orthogonal target rotation requires target_angle"
            )
        target_angle = _angle(target_angle, name="target_angle")
    else:
        target_angle = None

    if source_mode == "sitelayout" or target_mode == "sitelayout":
        if site_layout is None:
            raise EMTFRotationError(
                "site-layout rotation requires original SiteLayout metadata"
            )

    inputs = np.empty((n, tf.n_input, tf.n_input), dtype=float)
    outputs = np.empty((n, tf.n_output, tf.n_output), dtype=float)

    for index in range(n):
        if source_mode == target_mode == "sitelayout":
            inputs[index] = np.eye(tf.n_input)
            outputs[index] = np.eye(tf.n_output)
            continue

        if source_mode == "sitelayout" and target_mode == "orthogonal":
            inputs[index] = _site_input_transform(
                tf,
                site_layout,
                target_angle,
                inverse=False,
            )
            outputs[index] = _site_output_transform(
                tf,
                site_layout,
                target_angle,
                inverse=False,
            )
            continue

        old_angle = float(source_vector[index])
        if source_mode == "orthogonal" and target_mode == "orthogonal":
            inputs[index] = _orthogonal_input_transform(
                old_angle,
                target_angle,
            )
            outputs[index] = _orthogonal_output_transform(
                tf,
                old_angle,
                target_angle,
            )
            continue

        # Orthogonal -> original site layout.  FCU applies the inverse site
        # transforms at the old orthogonal azimuth.
        inputs[index] = _site_input_transform(
            tf,
            site_layout,
            old_angle,
            inverse=True,
        )
        outputs[index] = _site_output_transform(
            tf,
            site_layout,
            old_angle,
            inverse=True,
        )

    return RotationMatrices(
        input_matrix=inputs,
        output_matrix=outputs,
        source_mode=source_mode,
        target_mode=target_mode,
        source_angles=source_vector,
        target_angle=target_angle,
    )


def rotate_covariance(
    covariance: Any,
    matrices: Any,
    *,
    side: str,
) -> np.ndarray:
    """Rotate an inverse-signal or residual covariance matrix family.

    Parameters
    ----------
    covariance : array-like
        Shape ``(n_period, n, n)``.
    matrices : array-like
        The corresponding ``U`` or ``V`` matrices for every period.
    side : {"input", "output"}
        Descriptive validation label.  Both covariance kinds use
        ``R @ C @ R.T``; the side identifies which matrix family is supplied.
    """
    normalized = str(side).strip().lower()
    if normalized not in {"input", "output"}:
        raise ValueError("side must be 'input' or 'output'")
    cov = np.asarray(covariance)
    transforms = np.asarray(matrices, dtype=float)
    if cov.ndim != 3 or cov.shape[1] != cov.shape[2]:
        raise EMTFRotationError(
            f"{normalized} covariance must have shape (n, m, m)"
        )
    if transforms.shape != cov.shape:
        raise EMTFRotationError(
            f"{normalized} covariance/rotation shape mismatch: "
            f"{cov.shape} vs {transforms.shape}"
        )
    out = np.empty_like(cov, dtype=np.result_type(cov, complex))
    for index in range(cov.shape[0]):
        out[index] = _masked_bilinear(
            transforms[index],
            cov[index],
            transforms[index],
        )
    if not np.iscomplexobj(cov):
        out = np.real(out)
    return out


def _variance_from_covariance(
    inverse_signal: np.ndarray,
    residual: np.ndarray,
) -> np.ndarray:
    sdiag = np.diagonal(inverse_signal, axis1=1, axis2=2)
    ndiag = np.diagonal(residual, axis1=1, axis2=2)
    product = ndiag[:, :, None] * sdiag[:, None, :]
    if np.iscomplexobj(product):
        scale = np.maximum(1.0, np.abs(product.real))
        bad_imag = np.abs(product.imag) > 1.0e-7 * scale
        if np.any(bad_imag):
            warnings.warn(
                "rotated covariance diagonals yield variance with a "
                "non-negligible imaginary part; the real part is retained",
                EMTFRotationWarning,
                stacklevel=3,
            )
        product = product.real
    tiny_negative = (product < 0.0) & (product > -1.0e-12)
    if np.any(tiny_negative):
        product = np.array(product, copy=True)
        product[tiny_negative] = 0.0
    if np.any(product < 0.0):
        warnings.warn(
            "rotated covariance factors produced negative variance entries",
            EMTFRotationWarning,
            stacklevel=3,
        )
    return np.asarray(product, dtype=float)


def _rotate_variance_independent(
    variance: np.ndarray,
    matrices: RotationMatrices,
) -> np.ndarray:
    """Propagate component variances assuming zero cross-covariance."""
    var = np.asarray(variance, dtype=float)
    expected = (
        matrices.n_periods,
        matrices.output_matrix.shape[1],
        matrices.input_matrix.shape[1],
    )
    if var.shape != expected:
        raise EMTFRotationError(
            f"variance shape must be {expected}; got {var.shape}"
        )
    out = np.empty_like(var, dtype=float)
    for index in range(var.shape[0]):
        left = np.square(matrices.output_matrix[index])
        right = np.square(matrices.input_matrix[index])
        out[index] = _masked_bilinear(left, var[index], right)
    return out


def _rotate_variance_fcu(
    variance: np.ndarray,
    matrices: RotationMatrices,
) -> np.ndarray:
    """Replicate FCU's legacy direct variance-matrix rotation."""
    var = np.asarray(variance, dtype=float)
    expected = (
        matrices.n_periods,
        matrices.output_matrix.shape[1],
        matrices.input_matrix.shape[1],
    )
    if var.shape != expected:
        raise EMTFRotationError(
            f"variance shape must be {expected}; got {var.shape}"
        )
    out = np.empty_like(var, dtype=float)
    for index in range(var.shape[0]):
        out[index] = _masked_bilinear(
            matrices.output_matrix[index],
            var[index],
            matrices.input_matrix[index],
        )
    return out


def _handle_variance_without_full_covariance(
    variance: np.ndarray,
    matrices: RotationMatrices,
    *,
    policy: str,
) -> np.ndarray | None:
    normalized = str(policy).strip().lower().replace("-", "_")
    aliases = {
        "independent_components": "independent",
        "legacy": "fcu",
        "legacy_fcu": "fcu",
    }
    normalized = aliases.get(normalized, normalized)
    allowed = {"drop", "raise", "independent", "fcu"}
    if normalized not in allowed:
        raise ValueError(f"variance_policy must be one of {sorted(allowed)}")
    if normalized == "raise":
        raise EMTFRotationError(
            "VAR cannot be rotated exactly without both INVSIGCOV and "
            "RESIDCOV"
        )
    if normalized == "drop":
        warnings.warn(
            "dropping VAR because full covariance is unavailable for an "
            "exact rotation",
            ApproximateVarianceRotationWarning,
            stacklevel=3,
        )
        return None
    if normalized == "independent":
        warnings.warn(
            "rotating VAR under an independent-component assumption because "
            "full covariance is unavailable",
            ApproximateVarianceRotationWarning,
            stacklevel=3,
        )
        return _rotate_variance_independent(variance, matrices)

    warnings.warn(
        "using EMTF FCU legacy direct variance-matrix rotation; error bars "
        "are not statistically valid without full covariance",
        ApproximateVarianceRotationWarning,
        stacklevel=3,
    )
    return _rotate_variance_fcu(variance, matrices)


def _copy_estimate_with_data(estimate, data):
    result = estimate.copy()
    result.data = np.asarray(data)
    return result


def rotate_transfer_function(
    tf: TransferFunction,
    *,
    source_mode: str,
    target_mode: str = "orthogonal",
    target_angle: float | None = 0.0,
    source_angles: Any | None = None,
    site_layout: SiteLayout | None = None,
    variance_policy: str = "drop",
    unsupported_estimates: str = "drop",
) -> TransferFunction:
    """Rotate one matrix-valued transfer function and supported estimates.

    Full ``INVSIGCOV`` + ``RESIDCOV`` factors are rotated exactly and ``VAR``
    is recomputed from their diagonal products.  If full covariance is absent,
    ``variance_policy`` controls whether ``VAR`` is dropped, rejected, or
    transformed under an explicit approximation.
    """
    if not isinstance(tf, TransferFunction):
        raise TypeError("tf must be a TransferFunction")
    matrices = _rotation_matrices(
        tf,
        source_mode=source_mode,
        target_mode=target_mode,
        site_layout=site_layout,
        source_angles=source_angles,
        target_angle=target_angle,
    )
    if matrices.is_identity:
        return tf.copy()

    out = tf.copy()
    rotated = np.empty_like(tf.data, dtype=np.result_type(tf.data, complex))
    for index in range(tf.n_periods):
        rotated[index] = _masked_bilinear(
            matrices.output_matrix[index],
            tf.data[index],
            matrices.input_matrix[index],
        )
    if not np.iscomplexobj(tf.data):
        rotated = np.real(rotated)
    out.data = rotated
    out.estimates = {}

    invsig = tf.get_estimate("inverse_signal_covariance")
    resid = tf.get_estimate("residual_covariance")
    variance = tf.get_estimate("variance")

    rotated_s = None
    rotated_n = None
    if invsig is not None:
        expected = (tf.n_periods, tf.n_input, tf.n_input)
        if invsig.data.shape != expected:
            raise EMTFRotationError(
                "INVSIGCOV shape mismatch: expected "
                f"{expected}, got {invsig.data.shape}"
            )
        rotated_s = rotate_covariance(
            invsig.data,
            matrices.input_matrix,
            side="input",
        )
        out.add_estimate(
            _copy_estimate_with_data(invsig, rotated_s),
            key="inverse_signal_covariance",
            replace=True,
        )

    if resid is not None:
        expected = (tf.n_periods, tf.n_output, tf.n_output)
        if resid.data.shape != expected:
            raise EMTFRotationError(
                "RESIDCOV shape mismatch: expected "
                f"{expected}, got {resid.data.shape}"
            )
        rotated_n = rotate_covariance(
            resid.data,
            matrices.output_matrix,
            side="output",
        )
        out.add_estimate(
            _copy_estimate_with_data(resid, rotated_n),
            key="residual_covariance",
            replace=True,
        )

    if rotated_s is not None and rotated_n is not None:
        new_var = _variance_from_covariance(rotated_s, rotated_n)
        if variance is None:
            from .estimates import StatisticalEstimate

            variance = StatisticalEstimate(
                name="VAR",
                kind="variance",
                data=new_var,
            )
        out.add_estimate(
            _copy_estimate_with_data(variance, new_var),
            key="variance",
            replace=True,
        )
    elif variance is not None:
        new_var = _handle_variance_without_full_covariance(
            variance.data,
            matrices,
            policy=variance_policy,
        )
        if new_var is not None:
            out.add_estimate(
                _copy_estimate_with_data(variance, new_var),
                key="variance",
                replace=True,
            )

    handled = {
        id(item)
        for item in (invsig, resid, variance)
        if item is not None
    }
    policy = str(unsupported_estimates).strip().lower()
    if policy not in {"drop", "keep", "raise"}:
        raise ValueError(
            "unsupported_estimates must be 'drop', 'keep', or 'raise'"
        )
    for key, estimate in tf.estimates.items():
        if id(estimate) in handled:
            continue
        message = (
            f"statistical estimate {estimate.name!r} has no defined Phase 7 "
            "rotation rule"
        )
        if policy == "raise":
            raise EMTFRotationError(message)
        warnings.warn(
            message
            + (
                "; keeping it unchanged"
                if policy == "keep"
                else "; dropping it"
            ),
            UnsupportedEstimateRotationWarning,
            stacklevel=2,
        )
        if policy == "keep":
            preserved = estimate.copy()
            preserved.attrs["stale_after_rotation"] = True
            out.add_estimate(preserved, key=key, replace=True)

    out.attrs = dict(out.attrs)
    out.attrs["rotation"] = {
        "source_mode": matrices.source_mode,
        "target_mode": matrices.target_mode,
        "source_angles": (
            None
            if matrices.source_angles is None
            else matrices.source_angles.tolist()
        ),
        "target_angle": matrices.target_angle,
    }
    return out


def _source_angles_for_document(
    document: EMTF,
    tf: TransferFunction,
    explicit: Any | None,
    *,
    use_legacy_edi_rotation: bool,
) -> tuple[str, Any | None]:
    orientation = document.orientation
    if explicit is not None:
        return "orthogonal", explicit
    if orientation is not None and orientation.mode == "sitelayout":
        return "sitelayout", None
    if orientation is not None and orientation.mode == "orthogonal":
        angle = orientation.angle_to_geographic_north
        if angle is None:
            raise EMTFRotationError(
                "orthogonal document orientation has no angle metadata"
            )
        return "orthogonal", angle

    if not use_legacy_edi_rotation:
        raise EMTFRotationError(
            "transfer-function orientation is ambiguous; provide "
            "source_angles or enable use_legacy_edi_rotation"
        )
    if orientation is None:
        raise EMTFRotationError(
            "no orientation metadata are available for legacy EDI rotation"
        )

    extra = orientation.extra or {}
    tag = tf.name.lower()
    key = "edi_trot" if tag == "tipper" else "edi_zrot"
    raw = extra.get(key)
    if raw is None and tag == "tipper":
        raw = extra.get("edi_zrot")
        if raw is not None:
            key = "edi_zrot"
    if raw is None:
        raise EMTFRotationError(
            f"no historical EDI rotation vector is available for {tf.name}"
        )
    warnings.warn(
        f"interpreting historical {key.upper()} as an orthogonal source "
        "orientation; this opt-in assumption is needed for principal-axis "
        "or otherwise frequency-dependent EDI rotations",
        LegacyRotationAssumptionWarning,
        stacklevel=3,
    )
    return "orthogonal", raw


def _source_override_for_tf(
    source_angles: Any | Mapping[str, Any] | None,
    tf: TransferFunction,
) -> Any | None:
    if not isinstance(source_angles, Mapping):
        return source_angles
    candidates = [tf.name]
    definition = tf.definition
    if definition is not None:
        candidates.extend(
            [definition.tag, definition.name, *definition.aliases]
        )
    for candidate in candidates:
        if candidate in source_angles:
            return source_angles[candidate]
        low = str(candidate).lower()
        for key, value in source_angles.items():
            if str(key).lower() == low:
                return value
    return None


def _handle_derived(
    tf: TransferFunction,
    *,
    policy: str,
) -> TransferFunction | None:
    normalized = str(policy).strip().lower()
    if normalized not in {"drop", "keep", "raise"}:
        raise ValueError("derived_policy must be 'drop', 'keep', or 'raise'")
    message = (
        f"derived data type {tf.name!r} should be recomputed after rotating "
        "its primary transfer function"
    )
    if normalized == "raise":
        raise EMTFRotationError(message)
    warnings.warn(
        message
        + (
            "; keeping it marked stale"
            if normalized == "keep"
            else "; dropping it"
        ),
        DerivedDataRotationWarning,
        stacklevel=3,
    )
    if normalized == "drop":
        return None
    result = tf.copy()
    result.attrs["stale_after_rotation"] = True
    return result


def _rotation_history(
    orientation: OrientationMeta | None,
    *,
    target_mode: str,
    target_angle: float | None,
) -> OrientationMeta:
    old = orientation
    extra = dict(old.extra if old is not None else {})
    history = list(extra.get("rotation_history", []) or [])
    history.append(
        {
            "operation": "emtf_rotation",
            "source_mode": None if old is None else old.mode,
            "source_angle": (
                None if old is None else old.angle_to_geographic_north
            ),
            "target_mode": target_mode,
            "target_angle": target_angle,
        }
    )
    extra["rotation_history"] = history
    previous = None if old is None else old.rotation_info
    statement = (
        "pyCSAMT rotated transfer functions to original site layout"
        if target_mode == "sitelayout"
        else (
            "pyCSAMT rotated transfer functions to an orthogonal frame at "
            f"{float(target_angle):g} degrees clockwise from geographic north"
        )
    )
    rotation_info = statement if not previous else f"{previous}\n{statement}"
    return OrientationMeta(
        mode=target_mode,
        angle_to_geographic_north=(
            float(target_angle) if target_mode == "orthogonal" else None
        ),
        rotation_info=rotation_info,
        extra=extra,
    )


def rotate_emtf(
    document: EMTF,
    angle: float | None = 0.0,
    *,
    target: str = "orthogonal",
    inplace: bool = False,
    source_angles: Any | Mapping[str, Any] | None = None,
    use_legacy_edi_rotation: bool = False,
    variance_policy: str = "drop",
    unsupported_estimates: str = "drop",
    derived_policy: str = "drop",
) -> EMTF:
    """Rotate all primary transfer functions in an :class:`EMTF` document.

    Parameters
    ----------
    document : EMTF
        Source scientific document.
    angle : float, optional
        Target azimuth in degrees clockwise from geographic north. Ignored for
        ``target='sitelayout'``.
    target : {"orthogonal", "sitelayout"}
        Coordinate frame to produce.
    inplace : bool, default=False
        If ``True``, mutate and return ``document``.  The original physical
        ``SiteLayout`` object is never modified.
    source_angles : scalar, array-like, mapping, optional
        Explicit orthogonal source angle(s). A mapping may provide different
        vectors by TF tag/code, useful for historical principal-axis EDI.
    use_legacy_edi_rotation : bool, default=False
        Opt in to interpreting retained ``edi_zrot`` / ``edi_trot`` vectors as
        orthogonal source frames when document-level orientation is ambiguous.
    variance_policy : {"drop", "raise", "independent", "fcu"}
        Behavior for ``VAR`` when full covariance factors are unavailable.
    unsupported_estimates : {"drop", "keep", "raise"}
        Policy for estimates without a defined rotation law.
    derived_policy : {"drop", "keep", "raise"}
        Derived products are not rotated directly; normally they should be
        recomputed after rotating their primary transfer function.
    """
    if not isinstance(document, EMTF):
        raise TypeError("document must be an EMTF")
    target_mode = str(target).strip().lower().replace("_", "")
    if target_mode in {"site", "layout", "sitelayout"}:
        target_mode = "sitelayout"
        target_angle = None
    elif target_mode in {"orthogonal", "ortho"}:
        target_mode = "orthogonal"
        target_angle = _angle(angle, name="angle")
    else:
        raise ValueError("target must be 'orthogonal' or 'sitelayout'")

    result = document if inplace else deepcopy(document)
    original_layout = document.site_layout
    rotated_tfs: dict[str, TransferFunction] = {}

    for key, original_tf in document.transfer_functions.items():
        definition = original_tf.definition
        if definition is not None and definition.is_derived:
            derived = _handle_derived(original_tf, policy=derived_policy)
            if derived is not None:
                rotated_tfs[key] = derived
            continue

        explicit = _source_override_for_tf(source_angles, original_tf)
        source_mode, angles = _source_angles_for_document(
            document,
            original_tf,
            explicit,
            use_legacy_edi_rotation=use_legacy_edi_rotation,
        )
        rotated_tfs[key] = rotate_transfer_function(
            original_tf,
            source_mode=source_mode,
            target_mode=target_mode,
            target_angle=target_angle,
            source_angles=angles,
            site_layout=original_layout,
            variance_policy=variance_policy,
            unsupported_estimates=unsupported_estimates,
        )

    result.transfer_functions = rotated_tfs
    result.tags = tuple(rotated_tfs)
    result.orientation = _rotation_history(
        document.orientation,
        target_mode=target_mode,
        target_angle=target_angle,
    )
    # Historical EDI rotation vectors describe the old numerical frame.
    # Never carry them forward unchanged after rotating the data.
    rotation_extra = dict(result.orientation.extra)
    if target_mode == "orthogonal":
        normalized_angle = float(target_angle)
        count = 0 if result.periods is None else len(result.periods)
        rotation_extra["edi_zrot"] = [normalized_angle] * count
        rotation_extra["edi_had_zrot"] = True
        if result.tipper_tf is not None:
            rotation_extra["edi_trot"] = [normalized_angle] * count
            rotation_extra["edi_had_trot"] = True
        else:
            rotation_extra.pop("edi_trot", None)
            rotation_extra["edi_had_trot"] = False
        rotation_extra["edi_rotation_metadata_normalized"] = True
    else:
        rotation_extra.pop("edi_zrot", None)
        rotation_extra.pop("edi_trot", None)
        rotation_extra["edi_had_zrot"] = False
        rotation_extra["edi_had_trot"] = False
        rotation_extra["edi_rotation_metadata_normalized"] = True
    result.orientation.extra = rotation_extra
    # Never rotate/replace the physical field geometry.  In inplace mode the
    # exact original object remains attached; in copy mode deepcopy created a
    # detached but numerically identical layout.
    if inplace:
        result.site_layout = original_layout
    result.validate()
    return result
