# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Time-lapse EM monitoring for hydrogeological change detection.

Compares a sequence of EM resistivity models (acquired at different times)
to detect and quantify:

* Resistivity change Δlog₁₀(ρ) — raw geophysical signal
* Water-saturation change ΔSw — via Archie inverse
* Volumetric water-content change Δθ = φ · ΔSw
* Water-table displacement (rise/fall in metres)

Primary applications in an EM-hydro context:

* **TDEM time-lapse** — monitoring recharge, dewatering, or saline intrusion
  at 10–500 m depth.  High sensitivity to near-surface saturation changes.
* **AMT time-lapse** — seasonal or induced changes in regional aquifers and
  fractured basement (100–2 000 m).

All outputs are referenced to a *baseline* survey (default: index 0).  The
grids of all surveys must match (same ``x_centers`` and ``z_centers``).  Use
:func:`assert_compatible_grids` to check before computing.

Typical use
-----------
>>> from pycsamt.interp import ResistivityModel
>>> from pycsamt.interp.petrophysics import ArchieModel
>>> from pycsamt.interp.timelapse import TimeLapseEM
>>>
>>> tl = TimeLapseEM(
...     surveys=[model_dry, model_wet, model_recharge],
...     labels=["dry", "wet", "recharge"],
... )
>>> delta_rho = tl.resistivity_change()    # list of (n_z, n_x) arrays
>>> delta_Sw  = tl.saturation_change(ArchieModel(), rho_w=0.025)
>>> delta_wt  = tl.water_table_displacement(ArchieModel())  # (n_surveys-1, n_x)
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import Union

import numpy as np

from ..api.property import PyCSAMTObject
from ._base import ResistivityModel
from .petrophysics import (
    ArchieModel,
    WaxmanSmitsModel,
    water_table_from_profile,
)

__all__ = [
    "TimeLapseEM",
    "assert_compatible_grids",
]

_PetroModel = Union[ArchieModel, WaxmanSmitsModel]


# ─────────────────────────────────────────────────────────────────────────────
# Grid compatibility check
# ─────────────────────────────────────────────────────────────────────────────

def assert_compatible_grids(
    surveys: Sequence[ResistivityModel],
    *,
    rtol: float = 1e-4,
) -> None:
    """Raise ``ValueError`` if surveys have incompatible grids.

    Parameters
    ----------
    surveys : sequence of ResistivityModel
    rtol : float
        Relative tolerance for coordinate comparison (default 1e-4).
    """
    ref = surveys[0]
    for i, s in enumerate(surveys[1:], start=1):
        if s.rho_2d.shape != ref.rho_2d.shape:
            raise ValueError(
                f"Survey {i} has shape {s.rho_2d.shape} "
                f"but survey 0 has shape {ref.rho_2d.shape}. "
                "All surveys must share the same model grid."
            )
        if not np.allclose(s.x_centers, ref.x_centers, rtol=rtol):
            raise ValueError(
                f"Survey {i} x_centers differ from survey 0 by more than rtol={rtol}."
            )
        if not np.allclose(s.z_centers, ref.z_centers, rtol=rtol):
            raise ValueError(
                f"Survey {i} z_centers differ from survey 0 by more than rtol={rtol}."
            )


# ─────────────────────────────────────────────────────────────────────────────
# TimeLapseEM
# ─────────────────────────────────────────────────────────────────────────────

class TimeLapseEM(PyCSAMTObject):
    """Time-lapse EM analysis for hydrogeological change detection.

    Parameters
    ----------
    surveys : sequence of ResistivityModel
        Time-ordered EM inversion results.  All must share the same grid
        (checked on construction via :func:`assert_compatible_grids`).
    times : sequence of float, optional
        Time stamps for each survey (any consistent unit: days, months, etc.).
        Used for labelling only — no numerical time-derivative is computed.
    labels : sequence of str, optional
        Human-readable survey labels (e.g. ``['dry2022', 'wet2023']``).

    Attributes
    ----------
    n_surveys : int
    n_x, n_z : int
    """

    def __init__(
        self,
        surveys: Sequence[ResistivityModel],
        *,
        times: Sequence[float] | None = None,
        labels: Sequence[str] | None = None,
    ) -> None:
        if len(surveys) < 2:
            raise ValueError("At least two surveys are required.")
        assert_compatible_grids(surveys)
        self.surveys = list(surveys)
        self.times   = list(times) if times is not None else list(range(len(surveys)))
        self.labels  = list(labels) if labels is not None else [
            f"T{i:02d}" for i in range(len(surveys))
        ]
        if len(self.times) != len(surveys):
            raise ValueError("len(times) must equal len(surveys).")
        if len(self.labels) != len(surveys):
            raise ValueError("len(labels) must equal len(surveys).")

    # ── properties ────────────────────────────────────────────────────────

    @property
    def n_surveys(self) -> int:
        return len(self.surveys)

    @property
    def n_x(self) -> int:
        return int(self.surveys[0].rho_2d.shape[1])

    @property
    def n_z(self) -> int:
        return int(self.surveys[0].rho_2d.shape[0])

    # ── resistivity change ─────────────────────────────────────────────────

    def resistivity_change(
        self,
        baseline_idx: int = 0,
    ) -> list[np.ndarray]:
        r"""Resistivity change relative to the baseline survey.

        .. math::

            \Delta\log_{10}\rho_i = \log_{10}\rho_i - \log_{10}\rho_\text{baseline}

        Parameters
        ----------
        baseline_idx : int
            Index of the baseline survey (default 0).

        Returns
        -------
        list of ndarray (n_z, n_x)
            One array per non-baseline survey, in time order.
            Positive values indicate resistivity increase (drying/desaturation).
            Negative values indicate resistivity decrease (wetting/salinisation).
        """
        base = self.surveys[baseline_idx].rho_2d
        return [
            s.rho_2d - base
            for i, s in enumerate(self.surveys)
            if i != baseline_idx
        ]

    # ── saturation change ──────────────────────────────────────────────────

    def saturation_change(
        self,
        petro: _PetroModel,
        *,
        phi: float | np.ndarray = 0.25,
        rho_w: float = 0.025,
        baseline_idx: int = 0,
    ) -> list[np.ndarray]:
        r"""Water-saturation change Δ*Sw* via Archie/WS inverse.

        .. math::

            \Delta S_{w,i} = S_w(\rho_i) - S_w(\rho_\text{baseline})

        Parameters
        ----------
        petro : ArchieModel or WaxmanSmitsModel
            Petrophysical model for ρ → Sw inversion.
        phi : float or ndarray
            Porosity used in the inversion.  Scalar for a uniform profile;
            2-D array (n_z, n_x) for a spatially varying prior.
        rho_w : float
            Pore-water resistivity (Ω·m).
        baseline_idx : int
            Index of the baseline survey (default 0).

        Returns
        -------
        list of ndarray (n_z, n_x)
            ΔSw arrays, one per non-baseline survey.
            Positive = wetting (saturation increase).
            Negative = drying (saturation decrease).
        """
        archie = _to_archie(petro)
        phi_arr = np.broadcast_to(np.asarray(phi, dtype=float),
                                   self.surveys[0].rho_2d.shape).copy()

        def _sw_map(survey: ResistivityModel) -> np.ndarray:
            rho = 10.0 ** survey.rho_2d
            return np.clip(archie.saturation(rho, phi_arr, rho_w), 0.0, 1.0)

        sw_base = _sw_map(self.surveys[baseline_idx])
        return [
            _sw_map(s) - sw_base
            for i, s in enumerate(self.surveys)
            if i != baseline_idx
        ]

    # ── water-content change ───────────────────────────────────────────────

    def water_content_change(
        self,
        petro: _PetroModel,
        *,
        phi: float | np.ndarray = 0.25,
        rho_w: float = 0.025,
        baseline_idx: int = 0,
    ) -> list[np.ndarray]:
        r"""Volumetric water-content change Δθ = φ · ΔSw.

        Parameters
        ----------
        petro, phi, rho_w, baseline_idx
            Same as :meth:`saturation_change`.

        Returns
        -------
        list of ndarray (n_z, n_x)
            Δθ arrays (dimensionless, range ≈ −φ to +φ).
        """
        phi_arr = np.broadcast_to(
            np.asarray(phi, dtype=float), self.surveys[0].rho_2d.shape
        ).copy()
        delta_sw_list = self.saturation_change(
            petro, phi=phi, rho_w=rho_w, baseline_idx=baseline_idx
        )
        return [delta_sw * phi_arr for delta_sw in delta_sw_list]

    # ── water-table displacement ───────────────────────────────────────────

    def water_table_displacement(
        self,
        petro: _PetroModel,
        *,
        rho_w: float = 0.025,
        Sw_threshold: float = 0.85,
        min_depth: float = 0.5,
        baseline_idx: int = 0,
    ) -> np.ndarray:
        r"""Water-table depth change relative to the baseline survey.

        .. math::

            \Delta z_{wt} = z_{wt}(t_i) - z_{wt}(t_\text{baseline})

        Positive values → water table dropped (deeper, e.g. dry season).
        Negative values → water table rose (shallower, e.g. recharge).

        ``nan`` is returned where the water table cannot be detected in
        either the baseline or the comparison survey.

        Parameters
        ----------
        petro : ArchieModel or WaxmanSmitsModel
        rho_w : float
        Sw_threshold : float
            Saturation level that defines the water table (default 0.85).
        min_depth : float
            Minimum search depth in metres (default 0.5).
        baseline_idx : int

        Returns
        -------
        ndarray (n_surveys−1, n_x)
            Water-table displacement in metres, one row per non-baseline survey.
        """
        archie = _to_archie(petro)
        base_survey = self.surveys[baseline_idx]
        base_wt = self._water_table_map(
            base_survey, archie, rho_w, Sw_threshold, min_depth
        )

        others = [s for i, s in enumerate(self.surveys) if i != baseline_idx]
        rows = []
        for s in others:
            wt_i = self._water_table_map(s, archie, rho_w, Sw_threshold, min_depth)
            rows.append(wt_i - base_wt)

        return np.vstack(rows) if len(rows) > 1 else rows[0]

    # ── specific-survey water table ────────────────────────────────────────

    def water_table_map(
        self,
        petro: _PetroModel,
        *,
        rho_w: float = 0.025,
        Sw_threshold: float = 0.85,
        min_depth: float = 0.5,
    ) -> np.ndarray:
        """Water-table depth (m) for every survey and every column.

        Returns
        -------
        ndarray (n_surveys, n_x)
            ``nan`` where the water table could not be detected.
        """
        archie = _to_archie(petro)
        rows = [
            self._water_table_map(s, archie, rho_w, Sw_threshold, min_depth)
            for s in self.surveys
        ]
        return np.vstack(rows)

    # ── statistics across time ─────────────────────────────────────────────

    def resistivity_stats(
        self,
        baseline_idx: int = 0,
    ) -> dict:
        """Summary statistics of resistivity change across all surveys.

        Returns
        -------
        dict with keys:
            ``mean_delta``, ``std_delta``, ``max_increase``, ``max_decrease``
            — each an ndarray (n_z, n_x).
        """
        deltas = self.resistivity_change(baseline_idx=baseline_idx)
        stack  = np.stack(deltas, axis=0)   # (n_surveys-1, n_z, n_x)
        return {
            "mean_delta":    np.nanmean(stack, axis=0),
            "std_delta":     np.nanstd(stack,  axis=0),
            "max_increase":  np.nanmax(stack,  axis=0),
            "max_decrease":  np.nanmin(stack,  axis=0),
        }

    # ── private ────────────────────────────────────────────────────────────

    def _water_table_map(
        self,
        survey: ResistivityModel,
        archie: ArchieModel,
        rho_w: float,
        Sw_threshold: float,
        min_depth: float,
    ) -> np.ndarray:
        """Water-table depth per column for one survey."""
        wt = np.full(survey.n_x, np.nan)
        for ix in range(survey.n_x):
            depth = water_table_from_profile(
                survey.rho_2d[:, ix],
                survey.z_centers,
                archie,
                rho_w=rho_w,
                Sw_threshold=Sw_threshold,
                min_depth=min_depth,
            )
            if depth is not None:
                wt[ix] = depth
        return wt

    def __repr__(self) -> str:
        return (
            f"TimeLapseEM(n_surveys={self.n_surveys}, "
            f"n_z={self.n_z}, n_x={self.n_x}, "
            f"labels={self.labels})"
        )


# ─────────────────────────────────────────────────────────────────────────────
# Internal helpers
# ─────────────────────────────────────────────────────────────────────────────

def _to_archie(petro: _PetroModel) -> ArchieModel:
    """Return an ArchieModel; if WaxmanSmitsModel, approximate as Archie."""
    if isinstance(petro, ArchieModel):
        return petro
    if isinstance(petro, WaxmanSmitsModel):
        return ArchieModel(m=petro.m, n=petro.n, a=petro.a)
    raise TypeError(f"petro must be ArchieModel or WaxmanSmitsModel, got {type(petro)}")
