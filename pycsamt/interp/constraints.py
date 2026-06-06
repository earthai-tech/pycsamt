# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Field-measurement constraints for petrophysical calibration.

Hydrogeological field measurements — piezometer levels, pumping-test
transmissivity, slug-test K, and EC logs — can be used to anchor the
petrophysical parameters (ρ_w, Archie m/n, porosity prior) that drive the
quantitative :class:`~pycsamt.interp.hydromodel.EMHydroModel`.

:class:`ConstrainedCalibrator` accepts a list of constraint objects and
optimises one or more parameters of a :class:`~pycsamt.interp.hydromodel.PetrophysicalConfig`
to minimise the weighted least-squares misfit between model-predicted
hydrogeological quantities and the observations.

EM method relevance
-------------------
* **TDEM** — water-level constraints from piezometers are the most direct.
  ``rho_w`` is the primary free parameter (calibrated from well EC or TDS).
* **AMT** — pumping-test T constraints calibrate porosity_prior + rho_w
  jointly; fracture aperture is not fitted here but can be set manually.
* **MT** — EC constraints from deep formation waters calibrate rho_w for
  basin-scale aquifer studies.
* **EMAP** — EC logs along the profile calibrate rho_w laterally.

Typical use
-----------
>>> from pycsamt.interp.constraints import (
...     ConstrainedCalibrator,
...     WaterLevelConstraint,
...     PumpingTestConstraint,
... )
>>> from pycsamt.interp.hydromodel import EMHydroModel, PetrophysicalConfig
>>>
>>> cal = ConstrainedCalibrator(
...     constraints=[
...         WaterLevelConstraint(x=500.0,  depth_m=18.5),
...         WaterLevelConstraint(x=1200.0, depth_m=22.0),
...         PumpingTestConstraint(x=800.0, T_m2s=1.2e-3),
...     ],
...     calibrate_rho_w=True,
...     calibrate_phi_prior=True,
... )
>>> model   = EMHydroModel(rm, PetrophysicalConfig())
>>> result  = cal.fit(model)      # returns calibrated EMHydroResult
>>> print(cal.calibrated_config_) # fitted PetrophysicalConfig
>>> print(cal.misfit_history_)    # optimiser convergence

References
----------
.. [1] Schön, J. H. (2015). Physical Properties of Rocks, 2nd ed.
   Elsevier.  (Chapter 6: calibration of Archie parameters)
.. [2] Bussian, A. E. (1983). Geophysics, 48, 1258–1268.
"""

from __future__ import annotations

import dataclasses
import warnings
from dataclasses import dataclass, field
from typing import Any, List, Optional, Sequence, Tuple, Union

import numpy as np

try:
    from scipy.optimize import minimize, OptimizeResult
    _SCIPY_OK = True
except ImportError:
    _SCIPY_OK = False

from ..api.property import PyCSAMTObject
from ._base import ResistivityModel
from .petrophysics import ec_mscm_to_rho
from .hydromodel import EMHydroModel, EMHydroResult, PetrophysicalConfig

__all__ = [
    "WaterLevelConstraint",
    "PumpingTestConstraint",
    "SlugTestConstraint",
    "ECConstraint",
    "ConstrainedCalibrator",
]


# ─────────────────────────────────────────────────────────────────────────────
# Constraint dataclasses
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class WaterLevelConstraint(PyCSAMTObject):
    """Piezometer or well water-level measurement.

    The calibrator matches the EM-derived water-table depth at the nearest
    model column to the observed depth.

    Parameters
    ----------
    x : float
        Profile position of the measurement (m, same reference as the model).
    depth_m : float
        Observed water-table depth below surface (m, positive downward).
    uncertainty_m : float
        Measurement + representativeness uncertainty (m; default 1.0 m).
    station : str
        Optional label for reporting.
    """

    x: float
    depth_m: float
    uncertainty_m: float = 1.0
    station: str = ""

    def __post_init__(self) -> None:
        if self.depth_m < 0:
            raise ValueError("depth_m must be ≥ 0.")
        if self.uncertainty_m <= 0:
            raise ValueError("uncertainty_m must be positive.")


@dataclass
class PumpingTestConstraint(PyCSAMTObject):
    """Transmissivity (and optionally storativity) from a pumping test.

    Misfit is computed in log₁₀ space because T spans orders of magnitude.

    Parameters
    ----------
    x : float
        Profile position of the pumping well (m).
    T_m2s : float
        Observed transmissivity (m²/s).
    S : float, optional
        Observed storativity (dimensionless).  Currently used for reporting
        only; storativity calibration will be added in a future release.
    uncertainty_factor : float
        Multiplicative uncertainty on T (default 3 → ±0.5 in log₁₀).
    station : str
        Optional label.
    """

    x: float
    T_m2s: float
    S: Optional[float] = None
    uncertainty_factor: float = 3.0
    station: str = ""

    def __post_init__(self) -> None:
        if self.T_m2s <= 0:
            raise ValueError("T_m2s must be positive.")
        if self.uncertainty_factor < 1:
            raise ValueError("uncertainty_factor must be ≥ 1.")


@dataclass
class SlugTestConstraint(PyCSAMTObject):
    """Hydraulic conductivity K from a slug test or bail test.

    Matches the model K at the nearest cell to ``(x, depth_m)``.

    Parameters
    ----------
    x : float
        Profile position (m).
    K_ms : float
        Observed K (m/s).
    depth_m : float
        Depth of the tested interval mid-point (m).
    uncertainty_factor : float
        Multiplicative uncertainty (default 5 → ~0.7 in log₁₀).
    station : str
        Optional label.
    """

    x: float
    K_ms: float
    depth_m: float
    uncertainty_factor: float = 5.0
    station: str = ""

    def __post_init__(self) -> None:
        if self.K_ms <= 0:
            raise ValueError("K_ms must be positive.")
        if self.depth_m < 0:
            raise ValueError("depth_m must be ≥ 0.")


@dataclass
class ECConstraint(PyCSAMTObject):
    """Electrical conductivity of formation water (e.g., from an EC log or sample).

    Directly constrains ``rho_w`` via the EC ↔ ρ_w conversion.

    Parameters
    ----------
    x : float
        Profile position (m).
    ec_mscm : float
        Observed pore-water EC (mS/cm; ≈ 1 mS/cm for fresh water at 25 °C).
    uncertainty_mscm : float
        Measurement uncertainty (mS/cm; default 2.0).
    station : str
        Optional label.
    """

    x: float
    ec_mscm: float
    uncertainty_mscm: float = 2.0
    station: str = ""

    def __post_init__(self) -> None:
        if self.ec_mscm <= 0:
            raise ValueError("ec_mscm must be positive.")
        if self.uncertainty_mscm <= 0:
            raise ValueError("uncertainty_mscm must be positive.")


# ─────────────────────────────────────────────────────────────────────────────
# Calibrator
# ─────────────────────────────────────────────────────────────────────────────

_AnyConstraint = Union[
    WaterLevelConstraint,
    PumpingTestConstraint,
    SlugTestConstraint,
    ECConstraint,
]


class ConstrainedCalibrator(PyCSAMTObject):
    """Calibrate petrophysical parameters to match hydrogeological field data.

    Uses ``scipy.optimize.minimize`` (L-BFGS-B) to find the parameter set
    that minimises the weighted sum of squared residuals between the
    EM-derived hydrogeological quantities and the supplied observations.

    The optimisation is deterministic (single run from the current parameter
    values).  If convergence is poor, set ``n_restarts > 1`` to try multiple
    random starting points within the parameter bounds.

    Parameters
    ----------
    constraints : sequence of constraint objects
        Any mix of :class:`WaterLevelConstraint`, :class:`PumpingTestConstraint`,
        :class:`SlugTestConstraint`, :class:`ECConstraint`.
    calibrate_rho_w : bool
        Fit pore-water resistivity ρ_w (default True).  Almost always needed.
    calibrate_m : bool
        Fit Archie cementation exponent m (default False).
    calibrate_phi_prior : bool
        Fit prior porosity (default False).
    rho_w_bounds : (float, float)
        Search bounds for ρ_w in Ω·m (default 0.003 – 2.0).
    m_bounds : (float, float)
        Search bounds for Archie m (default 1.2 – 2.8).
    phi_bounds : (float, float)
        Search bounds for porosity prior (default 0.03 – 0.60).
    n_restarts : int
        Number of optimisation restarts with random initialisations (default 1).
    verbose : bool
        Print iteration count and final misfit (default False).

    Attributes (set after :meth:`fit`)
    -----------------------------------
    calibrated_config_ : PetrophysicalConfig
        The fitted configuration.
    misfit_history_ : list of float
        Final misfit value per restart.
    opt_result_ : OptimizeResult or None
        Raw scipy result from the best restart.
    """

    def __init__(
        self,
        constraints: Sequence[_AnyConstraint],
        *,
        calibrate_rho_w: bool = True,
        calibrate_m: bool = False,
        calibrate_phi_prior: bool = False,
        rho_w_bounds: Tuple[float, float] = (0.003, 2.0),
        m_bounds: Tuple[float, float] = (1.2, 2.8),
        phi_bounds: Tuple[float, float] = (0.03, 0.60),
        n_restarts: int = 1,
        verbose: bool = False,
    ) -> None:
        if not constraints:
            raise ValueError("At least one constraint is required.")
        if not (calibrate_rho_w or calibrate_m or calibrate_phi_prior):
            raise ValueError(
                "At least one of calibrate_rho_w, calibrate_m, "
                "calibrate_phi_prior must be True."
            )
        self.constraints     = list(constraints)
        self.calibrate_rho_w = calibrate_rho_w
        self.calibrate_m     = calibrate_m
        self.calibrate_phi_prior = calibrate_phi_prior
        self.rho_w_bounds    = rho_w_bounds
        self.m_bounds        = m_bounds
        self.phi_bounds      = phi_bounds
        self.n_restarts      = int(n_restarts)
        self.verbose         = verbose

        # fitted attributes
        self.calibrated_config_: Optional[PetrophysicalConfig] = None
        self.misfit_history_: List[float] = []
        self.opt_result_: Any = None

    # ── public ─────────────────────────────────────────────────────────────

    def fit(self, em_hydro_model: EMHydroModel) -> EMHydroResult:
        """Calibrate and return the best :class:`~pycsamt.interp.hydromodel.EMHydroResult`.

        Parameters
        ----------
        em_hydro_model : EMHydroModel
            The model to calibrate.  Its ``config`` provides the starting
            point for the optimisation.

        Returns
        -------
        EMHydroResult
            Model result evaluated at the calibrated parameters.
        """
        if not _SCIPY_OK:
            raise ImportError(
                "scipy is required for ConstrainedCalibrator.fit()."
            )

        x0, bounds = self._pack_params(em_hydro_model.config)
        rng = np.random.default_rng(42)

        best_x: Optional[np.ndarray] = None
        best_f: float = float("inf")
        best_opt: Any = None

        for restart in range(self.n_restarts):
            if restart == 0:
                x_start = x0.copy()
            else:
                lo = np.array([b[0] for b in bounds])
                hi = np.array([b[1] for b in bounds])
                x_start = rng.uniform(lo, hi)

            opt = minimize(
                fun=lambda x: self._objective(x, em_hydro_model),
                x0=x_start,
                method="L-BFGS-B",
                bounds=bounds,
                options={"maxiter": 200, "ftol": 1e-10, "gtol": 1e-7},
            )
            self.misfit_history_.append(float(opt.fun))
            if opt.fun < best_f:
                best_f = float(opt.fun)
                best_x = opt.x.copy()
                best_opt = opt

        self.opt_result_ = best_opt
        fitted_cfg = self._unpack_params(best_x, em_hydro_model.config)
        self.calibrated_config_ = fitted_cfg

        if self.verbose:
            print(
                f"ConstrainedCalibrator: final misfit={best_f:.4f}  "
                f"restarts={self.n_restarts}"
            )
            for k, v in dataclasses.asdict(fitted_cfg).items():
                if k in ("petro",):
                    print(f"  {k}: {fitted_cfg.petro!r}")
                elif k not in (
                    "kozeny_C", "kozeny_tortuosity", "fracture_depth_m",
                    "fracture_rho_matrix", "specific_storage", "min_wt_search_depth",
                ):
                    print(f"  {k}: {v}")

        calibrated_model = EMHydroModel(
            em_hydro_model.resistivity_model,
            fitted_cfg,
            method_tag=em_hydro_model.method_tag,
        )
        return calibrated_model.fit()

    def constraint_residuals(
        self, result: EMHydroResult
    ) -> List[dict]:
        """Return per-constraint residuals for the given result (for diagnostics)."""
        rows = []
        for c in self.constraints:
            ix = _nearest_x_ix(result.resistivity_model, c.x)
            row = {"type": type(c).__name__, "x_m": c.x, "station": getattr(c, "station", "")}

            if isinstance(c, WaterLevelConstraint):
                obs  = c.depth_m
                pred = float(result.water_table[ix])
                row.update({"observed": obs, "predicted": pred,
                            "residual_m": pred - obs,
                            "normalized": (pred - obs) / c.uncertainty_m})

            elif isinstance(c, PumpingTestConstraint):
                obs  = np.log10(c.T_m2s)
                pred = np.log10(max(result.transmissivity[ix], 1e-20))
                row.update({"observed_T_m2s": c.T_m2s,
                            "predicted_T_m2s": result.transmissivity[ix],
                            "residual_log10": pred - obs,
                            "normalized": (pred - obs) / np.log10(c.uncertainty_factor)})

            elif isinstance(c, SlugTestConstraint):
                iz   = _nearest_z_ix(result.resistivity_model, c.depth_m)
                obs  = np.log10(c.K_ms)
                pred = np.log10(max(float(result.hydraulic_K[iz, ix]), 1e-20))
                row.update({"observed_K_ms": c.K_ms,
                            "predicted_K_ms": float(result.hydraulic_K[iz, ix]),
                            "residual_log10": pred - obs,
                            "normalized": (pred - obs) / np.log10(c.uncertainty_factor)})

            elif isinstance(c, ECConstraint):
                rho_w_from_ec = float(ec_mscm_to_rho(c.ec_mscm))
                rho_w_model   = result.config.rho_w
                row.update({"observed_ec_mscm": c.ec_mscm,
                            "rho_w_from_ec": rho_w_from_ec,
                            "rho_w_model": rho_w_model,
                            "residual_rho_w": rho_w_model - rho_w_from_ec})
            rows.append(row)
        return rows

    # ── private optimisation helpers ───────────────────────────────────────

    def _pack_params(
        self, cfg: PetrophysicalConfig
    ) -> Tuple[np.ndarray, list]:
        """Extract free parameters as a 1-D array with bounds."""
        x0_list: List[float] = []
        bounds_list: List[Tuple[float, float]] = []

        if self.calibrate_rho_w:
            x0_list.append(cfg.rho_w)
            bounds_list.append(self.rho_w_bounds)

        if self.calibrate_m:
            m = cfg.petro.m if hasattr(cfg.petro, "m") else 1.8
            x0_list.append(m)
            bounds_list.append(self.m_bounds)

        if self.calibrate_phi_prior:
            x0_list.append(cfg.porosity_prior)
            bounds_list.append(self.phi_bounds)

        return np.array(x0_list, dtype=float), bounds_list

    def _unpack_params(
        self, x: np.ndarray, cfg: PetrophysicalConfig
    ) -> PetrophysicalConfig:
        """Reconstruct a PetrophysicalConfig from the optimised vector."""
        idx = 0
        updates: dict = {}
        if self.calibrate_rho_w:
            updates["rho_w"] = float(np.clip(x[idx], *self.rho_w_bounds))
            idx += 1
        if self.calibrate_m:
            petro = cfg.petro
            updates["petro"] = dataclasses.replace(
                petro, m=float(np.clip(x[idx], *self.m_bounds))
            )
            idx += 1
        if self.calibrate_phi_prior:
            updates["porosity_prior"] = float(
                np.clip(x[idx], *self.phi_bounds)
            )
        return dataclasses.replace(cfg, **updates)

    def _objective(
        self, x: np.ndarray, base_model: EMHydroModel
    ) -> float:
        """Weighted sum of squared normalised residuals."""
        cfg = self._unpack_params(x, base_model.config)
        try:
            result = EMHydroModel(
                base_model.resistivity_model,
                cfg,
                method_tag=base_model.method_tag,
            ).fit()
        except Exception:
            return 1e12

        misfit = 0.0
        for c in self.constraints:
            ix = _nearest_x_ix(base_model.resistivity_model, c.x)
            misfit += _constraint_residual_sq(c, result, ix)
        return float(misfit)


# ─────────────────────────────────────────────────────────────────────────────
# Residual helpers
# ─────────────────────────────────────────────────────────────────────────────

def _constraint_residual_sq(
    c: _AnyConstraint,
    result: EMHydroResult,
    ix: int,
) -> float:
    """Squared normalised residual for a single constraint."""
    if isinstance(c, WaterLevelConstraint):
        wt = float(result.water_table[ix])
        if not np.isfinite(wt):
            return 10.0   # mild penalty for undetected water table
        return ((wt - c.depth_m) / c.uncertainty_m) ** 2

    if isinstance(c, PumpingTestConstraint):
        T_model = float(result.transmissivity[ix])
        if T_model <= 0:
            return 25.0
        sigma_log = np.log10(max(c.uncertainty_factor, 1.001))
        return ((np.log10(T_model) - np.log10(c.T_m2s)) / sigma_log) ** 2

    if isinstance(c, SlugTestConstraint):
        iz      = _nearest_z_ix(result.resistivity_model, c.depth_m)
        K_model = float(result.hydraulic_K[iz, ix])
        if K_model <= 0:
            return 25.0
        sigma_log = np.log10(max(c.uncertainty_factor, 1.001))
        return ((np.log10(K_model) - np.log10(c.K_ms)) / sigma_log) ** 2

    if isinstance(c, ECConstraint):
        rho_w_obs   = float(ec_mscm_to_rho(c.ec_mscm))
        rho_w_model = result.config.rho_w
        sigma       = float(ec_mscm_to_rho(c.ec_mscm - c.uncertainty_mscm)
                            - ec_mscm_to_rho(c.ec_mscm + c.uncertainty_mscm)) / 2.0
        sigma = max(abs(sigma), 1e-4)
        return ((rho_w_model - rho_w_obs) / sigma) ** 2

    return 0.0


def _nearest_x_ix(model: ResistivityModel, x: float) -> int:
    return int(np.argmin(np.abs(model.x_centers - x)))


def _nearest_z_ix(model: ResistivityModel, z: float) -> int:
    return int(np.argmin(np.abs(model.z_centers - z)))
