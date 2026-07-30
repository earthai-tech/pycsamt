# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Multi-method EM fusion — merge shallow and deep resistivity models.

Different EM methods sample different depth windows:

+----------+--------------+--------------------------------------------+
| Method   | Depth range  | Hydrogeological target                     |
+==========+==============+============================================+
| TDEM     | 10 – 500 m   | Water table, shallow aquifer               |
+----------+--------------+--------------------------------------------+
| AMT      | 100 – 5 000 m| Regional aquifer, fractured basement        |
+----------+--------------+--------------------------------------------+
| MT       | 500 m – 50 km| Deep basins, geothermal                    |
+----------+--------------+--------------------------------------------+
| EMAP     | 10 – 1 000 m | Lateral conductors, saline intrusion       |
+----------+--------------+--------------------------------------------+

:class:`MultiMethodEMModel` merges two :class:`~pycsamt.interp.ResistivityModel`
instances onto a unified depth grid.  In the *overlap zone* where both methods
have coverage, the two models are blended by one of three strategies:

* ``'linear'``       — linear weight ramp across the overlap (default)
* ``'sigmoid'``      — smooth S-curve transition (no depth kinks)
* ``'rms_weighted'`` — constant weights derived from each model's RMS misfit

The fused :class:`~pycsamt.interp.ResistivityModel` feeds directly into
:class:`~pycsamt.interp.hydromodel.EMHydroModel`.

Typical use (TDEM + AMT)
------------------------
>>> from pycsamt.interp.fusion import MultiMethodEMModel
>>> from pycsamt.interp.hydromodel import EMHydroModel, PetrophysicalConfig
>>>
>>> fused = MultiMethodEMModel(
...     primary=tdem_model,
...     secondary=amt_model,
...     blend="sigmoid",
... ).merge()
>>>
>>> result = EMHydroModel(
...     fused, PetrophysicalConfig(), method_tag="TDEM+AMT"
... ).fit()
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..api.property import PyCSAMTObject
from ._base import ResistivityModel

__all__ = [
    "MultiMethodEMModel",
    "FusionDiagnostics",
]

_BLEND_MODES = ("linear", "sigmoid", "rms_weighted")


# ─────────────────────────────────────────────────────────────────────────────
# Diagnostics dataclass
# ─────────────────────────────────────────────────────────────────────────────


@dataclass
class FusionDiagnostics(PyCSAMTObject):
    """Metadata about a completed fusion operation.

    Attributes
    ----------
    z_overlap_start : float
        Top of the depth zone where both methods contribute (m).
    z_overlap_end : float
        Bottom of the overlap zone (m).
    has_overlap : bool
        ``False`` if the two depth ranges are disjoint (simple concatenation).
    blend_mode : str
    primary_method : str
    secondary_method : str
    primary_rms : float
    secondary_rms : float
    n_z_fused : int
        Total number of depth cells in the fused model.
    blend_weights : ndarray (n_z,)
        Primary-model weight at each depth cell (1 = all primary,
        0 = all secondary).
    """

    z_overlap_start: float
    z_overlap_end: float
    has_overlap: bool
    blend_mode: str
    primary_method: str
    secondary_method: str
    primary_rms: float
    secondary_rms: float
    n_z_fused: int
    blend_weights: np.ndarray


# ─────────────────────────────────────────────────────────────────────────────
# MultiMethodEMModel
# ─────────────────────────────────────────────────────────────────────────────


class MultiMethodEMModel(PyCSAMTObject):
    """Fuse two EM resistivity models onto a single depth grid.

    Parameters
    ----------
    primary : ResistivityModel
        The model trusted at **shallow** depths (e.g. TDEM, EMAP).
    secondary : ResistivityModel
        The model trusted at **deeper** depths (e.g. AMT, MT).
    primary_max_depth : float, optional
        Override the primary model's maximum contributing depth (m).
        Defaults to ``primary.z_centers[-1]``.
    secondary_min_depth : float, optional
        Override the secondary model's minimum contributing depth (m).
        Defaults to ``secondary.z_centers[0]``.
    blend : str
        Blend strategy in the overlap zone:

        ``'linear'`` (default)
            Linear weight ramp from primary (top of overlap) to secondary
            (bottom of overlap).
        ``'sigmoid'``
            Smooth S-curve parameterised by *sigmoid_k*.  Avoids the slope
            discontinuity at the ends of the transition zone.
        ``'rms_weighted'``
            Constant weights throughout: primary weight = rms_secondary /
            (rms_primary + rms_secondary).  Requires both models to carry a
            valid ``rms`` value.  Falls back to ``'linear'`` if either RMS
            is ``nan``.

    blend_overlap : float, optional
        Restrict the blend transition to a window of this width (m) centred
        on the mid-point of the natural overlap zone.  ``None`` uses the
        full overlap.
    z_grid : ndarray, optional
        Explicit output depth-cell centres (m).  Overrides the automatic
        union grid.  Both models are interpolated onto this grid.
    sigmoid_k : float
        Shape parameter for the sigmoid blend (m⁻¹; default 0.02, giving a
        smooth ~100 m transition for a 500 m overlap zone).

    Attributes
    ----------
    diagnostics_ : FusionDiagnostics or None
        Set after :meth:`merge` is called.

    Examples
    --------
    >>> fused_model = MultiMethodEMModel(
    ...     tdem_model, amt_model, blend="sigmoid", sigmoid_k=0.03
    ... ).merge()
    >>> fused_model.method
    'TDEM+AMT'
    """

    def __init__(
        self,
        primary: ResistivityModel,
        secondary: ResistivityModel,
        *,
        primary_max_depth: float | None = None,
        secondary_min_depth: float | None = None,
        blend: str = "linear",
        blend_overlap: float | None = None,
        z_grid: np.ndarray | None = None,
        sigmoid_k: float = 0.02,
    ) -> None:
        if blend not in _BLEND_MODES:
            raise ValueError(
                f"blend must be one of {_BLEND_MODES}, got {blend!r}."
            )
        self.primary = primary
        self.secondary = secondary
        self.primary_max_depth = primary_max_depth
        self.secondary_min_depth = secondary_min_depth
        self.blend = blend
        self.blend_overlap = blend_overlap
        self.z_grid = z_grid
        self.sigmoid_k = float(sigmoid_k)
        self.diagnostics_: FusionDiagnostics | None = None

    # ── public ─────────────────────────────────────────────────────────────

    def merge(self) -> ResistivityModel:
        """Produce the fused :class:`~pycsamt.interp.ResistivityModel`.

        Returns
        -------
        ResistivityModel
            Unified model on the output depth grid.  ``method`` is set to
            ``'<primary_method>+<secondary_method>'`` (e.g. ``'TDEM+AMT'``).
        """
        p, s = self.primary, self.secondary

        # ── build unified depth grid ────────────────────────────────────────
        z_out = self._build_z_grid()

        # ── interpolate both models to unified z and primary x ─────────────
        x_out = p.x_centers.copy()
        rho_p = _interp_model_to_grid(p, x_out, z_out)
        rho_s = _interp_model_to_grid(s, x_out, z_out)

        # ── effective contributing depth ranges ─────────────────────────────
        z_p_max = float(self.primary_max_depth or p.z_centers[-1])
        z_s_min = float(self.secondary_min_depth or s.z_centers[0])

        float(p.z_centers[0])
        float(s.z_centers[-1])

        has_overlap = z_s_min < z_p_max

        # overlap window (may be narrowed by blend_overlap)
        z_ov_start = z_s_min
        z_ov_end = z_p_max
        if self.blend_overlap is not None and has_overlap:
            z_mid = 0.5 * (z_ov_start + z_ov_end)
            half = 0.5 * float(self.blend_overlap)
            z_ov_start = z_mid - half
            z_ov_end = z_mid + half

        # ── compute blend weights per depth cell ───────────────────────────
        w = self._blend_weights(z_out, z_ov_start, z_ov_end, has_overlap)

        # ── secondary coverage mask (where secondary model has data) ────────
        s_has_data = (z_out >= s.z_centers[0]) & (z_out <= s.z_centers[-1])
        p_has_data = (z_out >= p.z_centers[0]) & (z_out <= p.z_centers[-1])

        # ── fuse ────────────────────────────────────────────────────────────
        w3d = w[:, np.newaxis]  # broadcast over x
        rho_fuse = np.where(
            p_has_data[:, np.newaxis] & s_has_data[:, np.newaxis],
            w3d * rho_p + (1.0 - w3d) * rho_s,
            np.where(p_has_data[:, np.newaxis], rho_p, rho_s),
        )

        # ── station metadata from primary ────────────────────────────────────
        sta_x = p.station_x if len(p.station_x) else x_out
        sta_names = p.station_names

        method_tag = f"{p.method}+{s.method}"

        self.diagnostics_ = FusionDiagnostics(
            z_overlap_start=float(z_ov_start),
            z_overlap_end=float(z_ov_end),
            has_overlap=bool(has_overlap),
            blend_mode=self.blend,
            primary_method=p.method,
            secondary_method=s.method,
            primary_rms=float(p.rms),
            secondary_rms=float(s.rms),
            n_z_fused=len(z_out),
            blend_weights=w.copy(),
        )

        return ResistivityModel.from_array(
            rho_fuse,
            x_out,
            z_out,
            station_x=sta_x,
            station_names=sta_names,
            method=method_tag,
            rms=float("nan"),
        )

    # ── private ────────────────────────────────────────────────────────────

    def _build_z_grid(self) -> np.ndarray:
        """Unified depth grid spanning both models."""
        if self.z_grid is not None:
            return np.asarray(self.z_grid, dtype=float)

        p_z = self.primary.z_centers
        s_z = self.secondary.z_centers

        z_p_max = float(self.primary_max_depth or p_z[-1])

        # primary grid + secondary cells strictly deeper than primary's max
        p_part = p_z[p_z <= z_p_max]
        s_part = s_z[s_z > z_p_max]

        # also include secondary cells in the overlap zone
        z_s_min = float(self.secondary_min_depth or s_z[0])
        s_overlap = s_z[(s_z >= z_s_min) & (s_z <= z_p_max)]

        all_z = np.unique(np.concatenate([p_part, s_overlap, s_part]))
        return np.sort(all_z)

    def _blend_weights(
        self,
        z: np.ndarray,
        z_lo: float,
        z_hi: float,
        has_overlap: bool,
    ) -> np.ndarray:
        """Primary-model weight array (1 at shallow, 0 at deep)."""
        w = np.ones(len(z))

        if not has_overlap:
            # hard boundary at the midpoint between primary max and secondary min
            w[z > z_lo] = 0.0
            return w

        if self.blend == "rms_weighted":
            rms_p = float(self.primary.rms)
            rms_s = float(self.secondary.rms)
            if (
                np.isfinite(rms_p)
                and np.isfinite(rms_s)
                and (rms_p + rms_s) > 0
            ):
                # lower RMS → higher weight (better fit)
                w_const = rms_s / (rms_p + rms_s)
                w[:] = 1.0
                in_overlap = (z > z_lo) & (z < z_hi)
                w[in_overlap] = w_const
                w[z >= z_hi] = 0.0
                return w
            # fall through to linear if RMS is not available

        # linear or sigmoid
        span = max(z_hi - z_lo, 1e-6)
        t = np.clip((z - z_lo) / span, 0.0, 1.0)  # 0 at z_lo, 1 at z_hi

        if self.blend == "sigmoid":
            # sigmoid centred at t=0.5 with sharpness proportional to sigmoid_k * span
            k = self.sigmoid_k * span
            w_blend = 1.0 - _sigmoid(t, k=k)
        else:  # linear
            w_blend = 1.0 - t

        # outside overlap: pure primary or pure secondary
        w[:] = np.where(z <= z_lo, 1.0, np.where(z >= z_hi, 0.0, w_blend))
        return w

    def __repr__(self) -> str:
        return (
            f"MultiMethodEMModel("
            f"primary='{self.primary.method}', "
            f"secondary='{self.secondary.method}', "
            f"blend='{self.blend}')"
        )


# ─────────────────────────────────────────────────────────────────────────────
# Internal helpers
# ─────────────────────────────────────────────────────────────────────────────


def _interp_model_to_grid(
    model: ResistivityModel,
    x_out: np.ndarray,
    z_out: np.ndarray,
) -> np.ndarray:
    """Bilinear interpolation of model.rho_2d onto (z_out, x_out) grid.

    * z-axis: ``numpy.interp`` per column (linear, extrapolates boundary values)
    * x-axis: ``numpy.interp`` per row after z-resampling
    """
    # Step 1: resample in z for each original x column
    rho_z = np.empty((len(z_out), model.n_x))
    for ix in range(model.n_x):
        rho_z[:, ix] = np.interp(z_out, model.z_centers, model.rho_2d[:, ix])

    # Step 2: resample in x for each output z row
    rho_zx = np.empty((len(z_out), len(x_out)))
    for iz in range(len(z_out)):
        rho_zx[iz, :] = np.interp(x_out, model.x_centers, rho_z[iz, :])

    return rho_zx


def _sigmoid(t: np.ndarray, *, k: float = 5.0) -> np.ndarray:
    """Sigmoid from 0 to 1 as t goes from 0 to 1 (centred at t=0.5)."""
    x = k * (t - 0.5)
    return 1.0 / (1.0 + np.exp(-x))
