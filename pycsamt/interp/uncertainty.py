# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Monte Carlo uncertainty propagation for EM hydro-geophysical models.

Petrophysical parameters (ρ_w, Archie *m*, *n*, porosity prior) are never
known exactly from the EM data alone.  This module propagates that
parametric uncertainty through the full chain

    (ρ_w, m, n, φ)  →  :class:`~pycsamt.interp.hydromodel.EMHydroModel`
    →  K, Sw, φ, water-table depth, T, S

by drawing *N* random realisations of the free parameters from user-specified
prior distributions, running :class:`~pycsamt.interp.hydromodel.EMHydroModel`
for each, and accumulating ensemble statistics (mean, std, P10, P50, P90).

EM method context
-----------------
* **TDEM** — ρ_w uncertainty dominates (water quality unknown in advance).
  Typical range: 5–100 Ω·m for freshwater targets.
* **AMT** — *m* uncertainty is significant in fractured basement
  (cementation exponent varies 1.3–2.5 in fractured vs porous media).
* **EMAP / MT** — mixed uncertainty; wide ρ_w range for basin brines.

Typical use
-----------
>>> from pycsamt.interp.uncertainty import (
...     UncertaintyBounds, MonteCarloHydro
... )
>>> bounds = UncertaintyBounds(
...     rho_w_range=(5.0, 80.0),      # fresh water range
...     m_range=(1.5, 2.2),
...     phi_prior_range=(0.15, 0.40),
... )
>>> mc = MonteCarloHydro(resistivity_model, config, bounds, n_samples=300)
>>> unc = mc.run()
>>> print(unc.p90_wt - unc.p10_wt)   # WT depth uncertainty (m)
>>> print(unc.cv_K)                    # coefficient of variation of K

References
----------
.. [1] Gómez-Treviño, E. et al. (2014). J. Appl. Geophys., 100, 1–10.
.. [2] Slater, L. (2007). *Near Surface Geophys.*, 5, 369–384.
   (Uncertainty in shallow geophysical aquifer mapping)
.. [3] Binley, A. et al. (2015). *Water Resour. Res.*, 51, 3837–3866.
"""

from __future__ import annotations

import csv
import dataclasses
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, List, Optional, Tuple, Union

import numpy as np

from ..api.property import PyCSAMTObject
from ._base import ResistivityModel
from .hydromodel import EMHydroModel, EMHydroResult, PetrophysicalConfig
from .petrophysics import ArchieModel, WaxmanSmitsModel

__all__ = [
    "UncertaintyBounds",
    "UncertaintyResult",
    "MonteCarloHydro",
]

PathLike = Union[str, Path]
_DIST_MODES = ("uniform", "normal")


# ─────────────────────────────────────────────────────────────────────────────
# Parameter bounds
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class UncertaintyBounds(PyCSAMTObject):
    """Prior distribution specification for Monte Carlo sampling.

    Each ``*_range`` parameter activates that parameter as a free variable.
    Un-specified parameters are held fixed at the central
    :class:`~pycsamt.interp.hydromodel.PetrophysicalConfig` values.

    Parameters
    ----------
    rho_w_range : (float, float), optional
        Pore-water resistivity range (Ω·m).
        *uniform*: (low, high); *normal*: (mean, std).
        Typical fresh-water range: ``(5, 100)``; brackish: ``(0.5, 20)``.
    m_range : (float, float), optional
        Archie cementation exponent.
        Typical: ``(1.3, 2.5)`` for mixed lithology.
    n_range : (float, float), optional
        Archie saturation exponent.
        Typical: ``(1.8, 2.5)`` for clean sand.
    phi_prior_range : (float, float), optional
        Porosity prior.  Typical: ``(0.10, 0.45)``.
    dist : str
        Sampling distribution: ``'uniform'`` (default, non-informative) or
        ``'normal'`` (Gaussian — *range* is interpreted as *(mean, std)*).

    Notes
    -----
    At least one ``*_range`` must be set.  ``rho_w_range`` alone covers
    the dominant source of uncertainty in most shallow EM surveys.
    """

    rho_w_range:     Optional[Tuple[float, float]] = None
    m_range:         Optional[Tuple[float, float]] = None
    n_range:         Optional[Tuple[float, float]] = None
    phi_prior_range: Optional[Tuple[float, float]] = None
    dist: str = "uniform"

    def __post_init__(self) -> None:
        if self.dist not in _DIST_MODES:
            raise ValueError(f"dist must be one of {_DIST_MODES}, got {self.dist!r}.")
        if self.n_free == 0:
            raise ValueError(
                "At least one *_range parameter must be set."
            )

    @property
    def n_free(self) -> int:
        """Number of free (uncertain) parameters."""
        return sum(r is not None for r in [
            self.rho_w_range, self.m_range,
            self.n_range, self.phi_prior_range,
        ])

    @property
    def free_names(self) -> List[str]:
        """Names of the free parameters."""
        names = []
        if self.rho_w_range is not None:     names.append("rho_w")
        if self.m_range is not None:          names.append("m")
        if self.n_range is not None:          names.append("n")
        if self.phi_prior_range is not None:  names.append("phi_prior")
        return names

    def sample(
        self,
        cfg: PetrophysicalConfig,
        n: int,
        rng: np.random.Generator,
    ) -> List[PetrophysicalConfig]:
        """Draw *n* parameter samples and return a list of configs.

        Parameters
        ----------
        cfg : PetrophysicalConfig
            Central (best-estimate) configuration — used as fallback for
            fixed parameters.
        n : int
            Number of Monte Carlo samples.
        rng : np.random.Generator

        Returns
        -------
        list of PetrophysicalConfig, length *n*
        """
        rho_w_vals   = _draw(self.rho_w_range,     cfg.rho_w,         n, rng, self.dist, (1e-3, 1e4))
        m_vals       = _draw(self.m_range,          cfg.petro.m if hasattr(cfg.petro, "m") else 1.8,
                             n, rng, self.dist, (1.0, 4.0))
        n_vals       = _draw(self.n_range,          cfg.petro.n if hasattr(cfg.petro, "n") else 2.0,
                             n, rng, self.dist, (1.0, 4.0))
        phi_vals     = _draw(self.phi_prior_range,  cfg.porosity_prior, n, rng, self.dist, (0.01, 0.80))

        configs = []
        for i in range(n):
            petro_i = _perturb_petro(cfg.petro, m=m_vals[i], n=n_vals[i])
            configs.append(dataclasses.replace(
                cfg,
                petro=petro_i,
                rho_w=float(rho_w_vals[i]),
                porosity_prior=float(phi_vals[i]),
            ))
        return configs


# ─────────────────────────────────────────────────────────────────────────────
# Uncertainty result
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class UncertaintyResult(PyCSAMTObject):
    """Ensemble statistics from a :class:`MonteCarloHydro` run.

    All 2-D arrays have shape ``(n_z, n_x)``; 1-D station arrays have
    shape ``(n_x,)``.  K values are stored on the **linear** scale (m/s).

    Attributes
    ----------
    resistivity_model : ResistivityModel
    config : PetrophysicalConfig
        Central configuration used to seed the Monte Carlo.
    bounds : UncertaintyBounds
    n_samples : int
    method_tag : str
    sampled_params : ndarray (n_samples, n_free)
        Matrix of sampled parameter values, columns in ``bounds.free_names`` order.

    Hydraulic conductivity K (m/s)
    ------------------------------
    mean_K, std_K, p10_K, p50_K, p90_K : ndarray (n_z, n_x)
    cv_K : ndarray (n_z, n_x)
        Coefficient of variation std_K / mean_K.  Dimensionless uncertainty map.

    Water saturation Sw
    -------------------
    mean_Sw, std_Sw, p10_Sw, p90_Sw : ndarray (n_z, n_x)

    Porosity φ
    ----------
    mean_phi, std_phi : ndarray (n_z, n_x)

    Water-table depth (m)
    ---------------------
    mean_wt, std_wt, p10_wt, p50_wt, p90_wt : ndarray (n_x,)
    wt_detection_rate : ndarray (n_x,)
        Fraction of samples in which the water table was detected (0–1).

    Transmissivity T (m²/s)
    -----------------------
    mean_T, std_T, p10_T, p50_T, p90_T : ndarray (n_x,)
    """

    resistivity_model:  ResistivityModel
    config:             PetrophysicalConfig
    bounds:             UncertaintyBounds
    n_samples:          int
    method_tag:         str = ""
    sampled_params:     np.ndarray = field(default_factory=lambda: np.array([]))

    # K maps
    mean_K:  np.ndarray = field(default_factory=lambda: np.array([]))
    std_K:   np.ndarray = field(default_factory=lambda: np.array([]))
    p10_K:   np.ndarray = field(default_factory=lambda: np.array([]))
    p50_K:   np.ndarray = field(default_factory=lambda: np.array([]))
    p90_K:   np.ndarray = field(default_factory=lambda: np.array([]))
    cv_K:    np.ndarray = field(default_factory=lambda: np.array([]))

    # Sw maps
    mean_Sw: np.ndarray = field(default_factory=lambda: np.array([]))
    std_Sw:  np.ndarray = field(default_factory=lambda: np.array([]))
    p10_Sw:  np.ndarray = field(default_factory=lambda: np.array([]))
    p90_Sw:  np.ndarray = field(default_factory=lambda: np.array([]))

    # porosity maps
    mean_phi: np.ndarray = field(default_factory=lambda: np.array([]))
    std_phi:  np.ndarray = field(default_factory=lambda: np.array([]))

    # water-table
    mean_wt:            np.ndarray = field(default_factory=lambda: np.array([]))
    std_wt:             np.ndarray = field(default_factory=lambda: np.array([]))
    p10_wt:             np.ndarray = field(default_factory=lambda: np.array([]))
    p50_wt:             np.ndarray = field(default_factory=lambda: np.array([]))
    p90_wt:             np.ndarray = field(default_factory=lambda: np.array([]))
    wt_detection_rate:  np.ndarray = field(default_factory=lambda: np.array([]))

    # transmissivity
    mean_T:  np.ndarray = field(default_factory=lambda: np.array([]))
    std_T:   np.ndarray = field(default_factory=lambda: np.array([]))
    p10_T:   np.ndarray = field(default_factory=lambda: np.array([]))
    p50_T:   np.ndarray = field(default_factory=lambda: np.array([]))
    p90_T:   np.ndarray = field(default_factory=lambda: np.array([]))

    # ── derived queries ────────────────────────────────────────────────────

    def prob_wt_shallower_than(
        self,
        depth_m: float,
        *,
        wt_ensemble: Optional[np.ndarray] = None,
    ) -> np.ndarray:
        """P(water-table depth < *depth_m*) per station column.

        Parameters
        ----------
        depth_m : float
            Reference depth (m, positive downward).
        wt_ensemble : ndarray (n_samples, n_x), optional
            Full water-table ensemble from :meth:`MonteCarloHydro.run_ensemble`.
            If ``None``, approximated from the stored P10/P50/P90 assuming
            a Gaussian distribution with mean=mean_wt, std=std_wt.

        Returns
        -------
        ndarray (n_x,) — probability in [0, 1]
        """
        if wt_ensemble is not None:
            return np.mean(wt_ensemble < depth_m, axis=0)
        from scipy.stats import norm as _norm
        p = np.zeros_like(self.mean_wt)
        valid = self.std_wt > 0
        p[valid] = _norm.cdf(depth_m, loc=self.mean_wt[valid], scale=self.std_wt[valid])
        mask_fixed = ~valid & np.isfinite(self.mean_wt)
        p[mask_fixed] = (self.mean_wt[mask_fixed] < depth_m).astype(float)
        return p

    def station_report(self) -> List[dict]:
        """Per-station uncertainty summary."""
        model = self.resistivity_model
        rows = []
        for ix, x in enumerate(model.x_centers):
            name = model.station_names[ix] if model.station_names and ix < len(model.station_names) else f"S{ix:03d}"
            rows.append({
                "station":          name,
                "x_m":              float(x),
                "mean_wt_m":        float(self.mean_wt[ix]),
                "std_wt_m":         float(self.std_wt[ix]),
                "p10_wt_m":         float(self.p10_wt[ix]),
                "p90_wt_m":         float(self.p90_wt[ix]),
                "wt_range_m":       float(self.p90_wt[ix] - self.p10_wt[ix]),
                "wt_detection_pct": float(self.wt_detection_rate[ix] * 100),
                "mean_T_m2s":       float(self.mean_T[ix]),
                "std_T_m2s":        float(self.std_T[ix]),
                "p10_T_m2s":        float(self.p10_T[ix]),
                "p90_T_m2s":        float(self.p90_T[ix]),
                "log10_T_range":    float(
                    np.log10(max(self.p90_T[ix], 1e-20))
                    - np.log10(max(self.p10_T[ix], 1e-20))
                ),
            })
        return rows

    def to_csv(self, path: PathLike) -> Path:
        """Write per-station uncertainty summary to CSV."""
        out = Path(path)
        out.parent.mkdir(parents=True, exist_ok=True)
        rows = self.station_report()
        if not rows:
            return out
        with out.open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)
        return out


# ─────────────────────────────────────────────────────────────────────────────
# Monte Carlo engine
# ─────────────────────────────────────────────────────────────────────────────

class MonteCarloHydro(PyCSAMTObject):
    """Monte Carlo uncertainty propagation for a quantitative EM hydro model.

    Parameters
    ----------
    resistivity_model : ResistivityModel
        Source inversion model.
    config : PetrophysicalConfig
        Central (best-estimate) petrophysical configuration.
    bounds : UncertaintyBounds
        Prior distributions for the free parameters.
    n_samples : int
        Number of Monte Carlo realisations (default 200).
        Typical: 100 for quick runs, 500–1000 for publication-quality.
    seed : int
        Random seed for reproducibility (default 42).
    method_tag : str, optional
        EM method label inherited by the output.
    verbose : bool
        Print progress every 50 iterations (default False).

    Examples
    --------
    >>> bounds = UncertaintyBounds(rho_w_range=(5.0, 100.0), m_range=(1.4, 2.4))
    >>> mc  = MonteCarloHydro(rm, cfg, bounds, n_samples=300)
    >>> unc = mc.run()
    >>> unc.cv_K.max()          # worst-case K uncertainty (fraction)
    >>> unc.p90_wt - unc.p10_wt  # WT depth spread per station (m)
    """

    def __init__(
        self,
        resistivity_model: ResistivityModel,
        config: PetrophysicalConfig,
        bounds: UncertaintyBounds,
        *,
        n_samples: int = 200,
        seed: int = 42,
        method_tag: str = "",
        verbose: bool = False,
    ) -> None:
        self.resistivity_model = resistivity_model
        self.config     = config
        self.bounds     = bounds
        self.n_samples  = int(n_samples)
        self.seed       = seed
        self.method_tag = method_tag
        self.verbose    = verbose

    # ── public ─────────────────────────────────────────────────────────────

    def run(self) -> UncertaintyResult:
        """Execute the Monte Carlo ensemble and return :class:`UncertaintyResult`.

        Returns
        -------
        UncertaintyResult
        """
        rng     = np.random.default_rng(self.seed)
        configs = self.bounds.sample(self.config, self.n_samples, rng)

        n_z = self.resistivity_model.n_z
        n_x = self.resistivity_model.n_x

        # pre-allocate ensemble arrays
        K_ens  = np.full((self.n_samples, n_z, n_x), np.nan)
        Sw_ens = np.full((self.n_samples, n_z, n_x), np.nan)
        phi_ens= np.full((self.n_samples, n_z, n_x), np.nan)
        wt_ens = np.full((self.n_samples, n_x),       np.nan)
        T_ens  = np.full((self.n_samples, n_x),       np.nan)

        # sampled parameter matrix
        param_mat = np.full((self.n_samples, self.bounds.n_free), np.nan)

        for i, cfg_i in enumerate(configs):
            if self.verbose and i % 50 == 0:
                print(f"  MC sample {i}/{self.n_samples} ...")

            try:
                result_i = EMHydroModel(
                    self.resistivity_model,
                    cfg_i,
                    method_tag=self.method_tag,
                ).fit()
            except Exception:
                continue

            K_ens[i]   = result_i.hydraulic_K
            Sw_ens[i]  = result_i.saturation
            phi_ens[i] = result_i.porosity
            wt_ens[i]  = result_i.water_table
            T_ens[i]   = result_i.transmissivity
            param_mat[i] = _extract_params(cfg_i, self.bounds)

        # ── K statistics ───────────────────────────────────────────────────
        mean_K = np.nanmean(K_ens,  axis=0)
        std_K  = np.nanstd(K_ens,   axis=0)
        p10_K  = np.nanpercentile(K_ens, 10, axis=0)
        p50_K  = np.nanpercentile(K_ens, 50, axis=0)
        p90_K  = np.nanpercentile(K_ens, 90, axis=0)
        with np.errstate(divide="ignore", invalid="ignore"):
            cv_K = np.where(mean_K > 0, std_K / mean_K, np.nan)

        # ── Sw statistics ──────────────────────────────────────────────────
        mean_Sw = np.nanmean(Sw_ens, axis=0)
        std_Sw  = np.nanstd(Sw_ens,  axis=0)
        p10_Sw  = np.nanpercentile(Sw_ens, 10, axis=0)
        p90_Sw  = np.nanpercentile(Sw_ens, 90, axis=0)

        # ── phi statistics ─────────────────────────────────────────────────
        mean_phi = np.nanmean(phi_ens, axis=0)
        std_phi  = np.nanstd(phi_ens,  axis=0)

        # ── water-table statistics ─────────────────────────────────────────
        detect_rate = np.mean(np.isfinite(wt_ens), axis=0)
        mean_wt = np.nanmean(wt_ens, axis=0)
        std_wt  = np.nanstd(wt_ens,  axis=0)
        p10_wt  = np.nanpercentile(wt_ens, 10, axis=0)
        p50_wt  = np.nanpercentile(wt_ens, 50, axis=0)
        p90_wt  = np.nanpercentile(wt_ens, 90, axis=0)

        # ── transmissivity statistics ──────────────────────────────────────
        mean_T = np.nanmean(T_ens, axis=0)
        std_T  = np.nanstd(T_ens,  axis=0)
        p10_T  = np.nanpercentile(T_ens, 10, axis=0)
        p50_T  = np.nanpercentile(T_ens, 50, axis=0)
        p90_T  = np.nanpercentile(T_ens, 90, axis=0)

        if self.verbose:
            print(f"  MC complete: {self.n_samples} samples  "
                  f"n_free={self.bounds.n_free}  "
                  f"free={self.bounds.free_names}")

        return UncertaintyResult(
            resistivity_model  = self.resistivity_model,
            config             = self.config,
            bounds             = self.bounds,
            n_samples          = self.n_samples,
            method_tag         = self.method_tag,
            sampled_params     = param_mat,
            mean_K=mean_K, std_K=std_K,
            p10_K=p10_K, p50_K=p50_K, p90_K=p90_K, cv_K=cv_K,
            mean_Sw=mean_Sw, std_Sw=std_Sw, p10_Sw=p10_Sw, p90_Sw=p90_Sw,
            mean_phi=mean_phi, std_phi=std_phi,
            mean_wt=mean_wt, std_wt=std_wt,
            p10_wt=p10_wt, p50_wt=p50_wt, p90_wt=p90_wt,
            wt_detection_rate=detect_rate,
            mean_T=mean_T, std_T=std_T,
            p10_T=p10_T, p50_T=p50_T, p90_T=p90_T,
        )

    def run_ensemble(self) -> Tuple[UncertaintyResult, np.ndarray, np.ndarray]:
        """Like :meth:`run` but also returns the raw WT and T ensembles.

        Returns
        -------
        unc : UncertaintyResult
        wt_ensemble : ndarray (n_samples, n_x)
        T_ensemble  : ndarray (n_samples, n_x)
        """
        rng     = np.random.default_rng(self.seed)
        configs = self.bounds.sample(self.config, self.n_samples, rng)
        n_z = self.resistivity_model.n_z
        n_x = self.resistivity_model.n_x

        K_ens  = np.full((self.n_samples, n_z, n_x), np.nan)
        Sw_ens = np.full((self.n_samples, n_z, n_x), np.nan)
        phi_ens= np.full((self.n_samples, n_z, n_x), np.nan)
        wt_ens = np.full((self.n_samples, n_x), np.nan)
        T_ens  = np.full((self.n_samples, n_x), np.nan)
        param_mat = np.full((self.n_samples, self.bounds.n_free), np.nan)

        for i, cfg_i in enumerate(configs):
            try:
                r = EMHydroModel(self.resistivity_model, cfg_i,
                                 method_tag=self.method_tag).fit()
            except Exception:
                continue
            K_ens[i]  = r.hydraulic_K;  Sw_ens[i]  = r.saturation
            phi_ens[i]= r.porosity;      wt_ens[i]  = r.water_table
            T_ens[i]  = r.transmissivity
            param_mat[i] = _extract_params(cfg_i, self.bounds)

        mean_K   = np.nanmean(K_ens,  axis=0)
        unc = UncertaintyResult(
            resistivity_model=self.resistivity_model, config=self.config,
            bounds=self.bounds, n_samples=self.n_samples, method_tag=self.method_tag,
            sampled_params=param_mat,
            mean_K=mean_K, std_K=np.nanstd(K_ens, axis=0),
            p10_K=np.nanpercentile(K_ens,10,axis=0),
            p50_K=np.nanpercentile(K_ens,50,axis=0),
            p90_K=np.nanpercentile(K_ens,90,axis=0),
            cv_K=np.where(mean_K>0, np.nanstd(K_ens,axis=0)/np.where(mean_K>0,mean_K,np.nan), np.nan),
            mean_Sw=np.nanmean(Sw_ens,axis=0), std_Sw=np.nanstd(Sw_ens,axis=0),
            p10_Sw=np.nanpercentile(Sw_ens,10,axis=0),
            p90_Sw=np.nanpercentile(Sw_ens,90,axis=0),
            mean_phi=np.nanmean(phi_ens,axis=0), std_phi=np.nanstd(phi_ens,axis=0),
            mean_wt=np.nanmean(wt_ens,axis=0), std_wt=np.nanstd(wt_ens,axis=0),
            p10_wt=np.nanpercentile(wt_ens,10,axis=0),
            p50_wt=np.nanpercentile(wt_ens,50,axis=0),
            p90_wt=np.nanpercentile(wt_ens,90,axis=0),
            wt_detection_rate=np.mean(np.isfinite(wt_ens),axis=0),
            mean_T=np.nanmean(T_ens,axis=0), std_T=np.nanstd(T_ens,axis=0),
            p10_T=np.nanpercentile(T_ens,10,axis=0),
            p50_T=np.nanpercentile(T_ens,50,axis=0),
            p90_T=np.nanpercentile(T_ens,90,axis=0),
        )
        return unc, wt_ens, T_ens


# ─────────────────────────────────────────────────────────────────────────────
# Internal helpers
# ─────────────────────────────────────────────────────────────────────────────

def _draw(
    rng_spec: Optional[Tuple[float, float]],
    central: float,
    n: int,
    rng: np.random.Generator,
    dist: str,
    physical_clip: Tuple[float, float],
) -> np.ndarray:
    """Draw n samples from prior; return central value repeated if not free."""
    if rng_spec is None:
        return np.full(n, central)
    lo, hi = rng_spec
    if dist == "uniform":
        samples = rng.uniform(lo, hi, size=n)
    else:
        samples = rng.normal(lo, hi, size=n)   # (mean, std) convention
    return np.clip(samples, *physical_clip)


def _perturb_petro(
    petro: Union[ArchieModel, WaxmanSmitsModel],
    m: float,
    n: float,
) -> Union[ArchieModel, WaxmanSmitsModel]:
    """Return a new petrophysical model with updated m and n."""
    return dataclasses.replace(petro, m=float(m), n=float(n))


def _extract_params(
    cfg: PetrophysicalConfig,
    bounds: UncertaintyBounds,
) -> np.ndarray:
    """Extract free-parameter values from a config as a 1-D array."""
    vals = []
    if bounds.rho_w_range is not None:     vals.append(cfg.rho_w)
    if bounds.m_range is not None:          vals.append(cfg.petro.m)
    if bounds.n_range is not None:          vals.append(cfg.petro.n)
    if bounds.phi_prior_range is not None:  vals.append(cfg.porosity_prior)
    return np.array(vals, dtype=float)
