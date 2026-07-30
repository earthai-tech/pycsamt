# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Quantitative hydro-geophysical model from EM resistivity sections.

Converts a 2-D :class:`~pycsamt.interp.ResistivityModel` (from TDEM, AMT, MT,
or EMAP inversion) into quantitative hydrogeological deliverables:

* Porosity and water-saturation maps (n_z × n_x)
* Hydraulic conductivity map via Archie → Kozeny-Carman, or cubic-law for
  fractured basement (AMT targets)
* Water-table depth per station column
* Transmissivity and storativity profiles
* Dar-Zarrouk parameters (transverse resistance TR, longitudinal conductance S)
* TDS indicator from pore-water resistivity

The main entry point is :class:`EMHydroModel` which wraps the transforms
already implemented in :mod:`pycsamt.interp.petrophysics`.  The workflow is
intentionally *deterministic* (single best-estimate); uncertainty propagation
will be added in a later phase.

Typical use
-----------
>>> from pycsamt.interp import ResistivityModel
>>> from pycsamt.interp.petrophysics import ArchieModel
>>> from pycsamt.interp.hydromodel import EMHydroModel, PetrophysicalConfig
>>>
>>> cfg = PetrophysicalConfig(
...     petro=ArchieModel(m=1.8, n=2.0),
...     rho_w=0.025,
...     porosity_prior=0.25,
... )
>>> result = EMHydroModel(resistivity_model=rm, config=cfg).fit()
>>> print(result.water_table)  # (n_x,) depth in metres
>>> print(result.transmissivity)  # (n_x,) m²/s

References
----------
.. [1] Archie, G. E. (1942). Trans. AIME, 146, 54–62.
.. [2] Kozeny, J. (1927). Sitzungsber. Akad. Wiss. Wien, 136(2a), 271–306.
.. [3] Maillet, R. (1947). Geophysics, 12, 509–519.  (Dar-Zarrouk)
.. [4] Bostick, F. X. (1977). Workshop on Electrical Methods, SEG.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Union

import numpy as np

from ..api.property import PyCSAMTObject
from ._base import ResistivityModel
from .petrophysics import (
    ArchieModel,
    WaxmanSmitsModel,
    fractured_zone_K,
    kozeny_carman_K,
    rho_w_to_tds,
    water_table_from_profile,
)

__all__ = [
    "PetrophysicalConfig",
    "EMHydroResult",
    "EMHydroModel",
]

PathLike = Union[str, Path]
_PetroModel = Union[ArchieModel, WaxmanSmitsModel]


# ─────────────────────────────────────────────────────────────────────────────
# Configuration dataclass
# ─────────────────────────────────────────────────────────────────────────────


@dataclass
class PetrophysicalConfig(PyCSAMTObject):
    """All petrophysical and hydraulic parameters needed by :class:`EMHydroModel`.

    Parameters
    ----------
    petro : ArchieModel or WaxmanSmitsModel
        Petrophysical model for ρ ↔ (φ, Sw) transforms.
        Default is ``ArchieModel(m=1.8, n=2.0)`` — appropriate for clean
        sands targeted by shallow TDEM and AMT surveys.
    rho_w : float
        Pore-water resistivity (Ω·m).  Default 20 Ω·m (EC ≈ 0.5 mS/cm,
        potable fresh water at 25 °C; range: 0.2 Ω·m seawater – 100+ Ω·m
        pristine groundwater).  Constrained by water-quality measurements
        or EC logs via :mod:`pycsamt.interp.constraints`.
    porosity_prior : float
        Reference porosity used in the vadose zone and as an upper-clip for
        unrealistic Archie-inferred porosities (default 0.25).
    Sw_water_table_threshold : float
        Saturation value that marks the water table in the Archie inverse
        (default 0.85; lower it for coarser-grained aquifers).
    d50_m : float
        Median grain size (m) for Kozeny-Carman K (default 2.5 × 10⁻⁴ m,
        fine sand).  Increase to ~1 × 10⁻³ m for coarse sand/gravel.
    kozeny_C : float
        Kozeny constant × shape factor (default 180 for spheres).
    kozeny_tortuosity : float
        Tortuosity correction in Kozeny-Carman (default 0.5).
    fracture_depth_m : float or None
        Depth (m) below which :func:`~pycsamt.interp.petrophysics.fractured_zone_K`
        replaces Kozeny-Carman.  Use for AMT/MT basement targets.  ``None``
        disables fracture K everywhere (default: None).
    fracture_rho_matrix : float
        Background resistivity of intact rock for the cubic-law fracture K
        (Ω·m; default 5 000).
    specific_storage : float
        Specific storage Ss (m⁻¹) for confined storativity (default 10⁻⁴).
    min_wt_search_depth : float
        Minimum depth (m) for water-table detection; skips near-surface noise
        (default 0.5 m).
    """

    petro: _PetroModel = field(
        default_factory=lambda: ArchieModel(m=1.8, n=2.0)
    )
    rho_w: float = 20.0
    porosity_prior: float = 0.25
    Sw_water_table_threshold: float = 0.85
    d50_m: float = 2.5e-4
    kozeny_C: float = 180.0
    kozeny_tortuosity: float = 0.5
    fracture_depth_m: float | None = None
    fracture_rho_matrix: float = 5000.0
    specific_storage: float = 1e-4
    min_wt_search_depth: float = 0.5

    def __post_init__(self) -> None:
        if not isinstance(self.petro, (ArchieModel, WaxmanSmitsModel)):
            raise TypeError("petro must be ArchieModel or WaxmanSmitsModel.")
        if self.rho_w <= 0:
            raise ValueError("rho_w must be positive.")
        if not 0 < self.porosity_prior < 1:
            raise ValueError("porosity_prior must be in (0, 1).")


# ─────────────────────────────────────────────────────────────────────────────
# Result dataclass
# ─────────────────────────────────────────────────────────────────────────────


@dataclass
class EMHydroResult(PyCSAMTObject):
    """Quantitative hydrogeological output of one EM resistivity section.

    All 2-D arrays have shape ``(n_z, n_x)`` matching the source
    :class:`~pycsamt.interp._base.ResistivityModel`; 1-D station arrays have
    shape ``(n_x,)``.

    Attributes
    ----------
    resistivity_model : ResistivityModel
    config : PetrophysicalConfig
    porosity : ndarray (n_z, n_x)
        Effective porosity per cell (dimensionless, 0–1).
        Below the water table: Archie inverse at Sw = 1.
        Above the water table: ``config.porosity_prior``.
    saturation : ndarray (n_z, n_x)
        Water saturation Sw per cell (0–1).
        Below the water table: 1.0.
        Above the water table: Archie inverse at φ = porosity_prior.
    hydraulic_K : ndarray (n_z, n_x)
        Hydraulic conductivity K (m/s).
        Porous media (z < fracture_depth_m): Kozeny-Carman from φ.
        Fractured basement (z ≥ fracture_depth_m): cubic-law from ρ contrast.
    water_table : ndarray (n_x,)
        Estimated water-table depth (m, positive downward).
        ``nan`` where detection failed (resistivity column shows no
        saturation transition above the Sw threshold).
    transmissivity : ndarray (n_x,)
        Aquifer transmissivity T = ∫K dz over the saturated zone (m²/s).
    storativity_confined : ndarray (n_x,)
        Confined storativity S = Ss × b_sat (dimensionless).
    storativity_unconfined : ndarray (n_x,)
        Unconfined storativity ≈ mean porosity in saturated zone (≈ specific yield).
    dar_zarrouk_TR : ndarray (n_x,)
        Transverse resistance TR = Σ ρᵢ hᵢ (Ω·m²) over all cells.
    dar_zarrouk_S : ndarray (n_x,)
        Longitudinal conductance S = Σ hᵢ/ρᵢ (siemens) over all cells.
    tds : float
        Estimated total dissolved solids (mg/L) from ``config.rho_w``.
    method_tag : str
        Source EM method label (e.g. ``'TDEM'``, ``'AMT'``).
    metadata : dict
        Provenance and parameter snapshot.
    """

    resistivity_model: ResistivityModel
    config: PetrophysicalConfig
    porosity: np.ndarray
    saturation: np.ndarray
    hydraulic_K: np.ndarray
    water_table: np.ndarray
    transmissivity: np.ndarray
    storativity_confined: np.ndarray
    storativity_unconfined: np.ndarray
    dar_zarrouk_TR: np.ndarray
    dar_zarrouk_S: np.ndarray
    tds: float
    method_tag: str = ""
    metadata: dict = field(default_factory=dict)

    # ── per-station summary ────────────────────────────────────────────────

    def station_report(self) -> list[dict]:
        """Return one dict per model column with key hydro indicators."""
        model = self.resistivity_model
        rows: list[dict] = []
        for ix, x in enumerate(model.x_centers):
            name = _station_name(model, ix)
            wt = float(self.water_table[ix])
            sat_mask = model.z_centers >= (wt if np.isfinite(wt) else 0.0)
            phi_sat = (
                float(np.nanmean(self.porosity[sat_mask, ix]))
                if sat_mask.any()
                else float("nan")
            )
            K_sat = (
                float(np.nanmean(self.hydraulic_K[sat_mask, ix]))
                if sat_mask.any()
                else float("nan")
            )
            rows.append(
                {
                    "station": name,
                    "x_m": float(x),
                    "water_table_m": wt,
                    "mean_porosity_sat": phi_sat,
                    "mean_K_ms": K_sat,
                    "transmissivity_m2s": float(self.transmissivity[ix]),
                    "storativity_confined": float(
                        self.storativity_confined[ix]
                    ),
                    "storativity_unconfined": float(
                        self.storativity_unconfined[ix]
                    ),
                    "dar_zarrouk_TR_ohm_m2": float(self.dar_zarrouk_TR[ix]),
                    "dar_zarrouk_S_siemens": float(self.dar_zarrouk_S[ix]),
                    "tds_mg_per_L": float(self.tds),
                }
            )
        return rows

    # ── CSV export ─────────────────────────────────────────────────────────

    def to_csv(self, path: PathLike) -> Path:
        """Write cell-level hydro parameters to a flat CSV file."""
        out = Path(path)
        out.parent.mkdir(parents=True, exist_ok=True)
        model = self.resistivity_model
        with out.open("w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(
                [
                    "station",
                    "x_m",
                    "z_m",
                    "rho_log10",
                    "rho_ohm_m",
                    "porosity",
                    "saturation",
                    "hydraulic_K_ms",
                ]
            )
            for ix, x in enumerate(model.x_centers):
                name = _station_name(model, ix)
                for iz, z in enumerate(model.z_centers):
                    rlog = float(model.rho_2d[iz, ix])
                    w.writerow(
                        [
                            name,
                            float(x),
                            float(z),
                            rlog,
                            10.0**rlog,
                            float(self.porosity[iz, ix]),
                            float(self.saturation[iz, ix]),
                            float(self.hydraulic_K[iz, ix]),
                        ]
                    )
        return out

    def station_report_csv(self, path: PathLike) -> Path:
        """Write per-station hydro summary to CSV."""
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

    def to_dataframe(self) -> Any:
        """Return station summary as a :class:`pandas.DataFrame` (requires pandas)."""
        try:
            import pandas as pd
        except ImportError as exc:
            raise ImportError(
                "pandas is required for to_dataframe()."
            ) from exc
        return pd.DataFrame(self.station_report())


# ─────────────────────────────────────────────────────────────────────────────
# Main model class
# ─────────────────────────────────────────────────────────────────────────────


class EMHydroModel(PyCSAMTObject):
    """Quantitative hydrogeological model from an EM resistivity section.

    Wires the petrophysical transforms in :mod:`pycsamt.interp.petrophysics`
    onto a full 2-D :class:`~pycsamt.interp.ResistivityModel` to produce
    spatially continuous porosity, saturation, hydraulic-conductivity, and
    water-table maps.

    The computation follows a two-pass strategy to break the Archie
    chicken-and-egg (φ needs Sw, Sw needs φ):

    1. **Water-table pass** — scan each column with the Archie-inverse Sw
       estimator to find the saturation front (Sw ≥ threshold).
    2. **Cell-property pass** — for cells **below** the water table assume
       full saturation (Sw = 1) and solve for φ; for cells **above** use
       ``porosity_prior`` and solve for Sw.

    Parameters
    ----------
    resistivity_model : ResistivityModel
        Source inversion model.
    config : PetrophysicalConfig, optional
        Full parameter bundle.  If ``None``, default parameters are used.
        Individual keyword arguments below override specific fields.
    petro : ArchieModel or WaxmanSmitsModel, optional
        Petrophysical model.
    rho_w : float, optional
        Pore-water resistivity (Ω·m).
    porosity_prior : float, optional
        Prior porosity.
    method_tag : str, optional
        EM method label (``'TDEM'``, ``'AMT'``, ``'MT'``, ``'EMAP'``).

    Examples
    --------
    >>> cfg = PetrophysicalConfig(
    ...     rho_w=0.03, porosity_prior=0.20, fracture_depth_m=300.0
    ... )
    >>> result = EMHydroModel(rm, cfg, method_tag="AMT").fit()
    >>> result.water_table  # array of depths, one per x-column
    >>> result.transmissivity  # T profile (m²/s)
    """

    def __init__(
        self,
        resistivity_model: ResistivityModel,
        config: PetrophysicalConfig | None = None,
        *,
        petro: _PetroModel | None = None,
        rho_w: float | None = None,
        porosity_prior: float | None = None,
        method_tag: str = "",
    ) -> None:
        self.resistivity_model = resistivity_model
        # Build config, applying any overrides
        base = config if config is not None else PetrophysicalConfig()
        if petro is not None:
            base = _replace_config(base, petro=petro)
        if rho_w is not None:
            base = _replace_config(base, rho_w=rho_w)
        if porosity_prior is not None:
            base = _replace_config(base, porosity_prior=porosity_prior)
        self.config = base
        self.method_tag = method_tag
        self._result: EMHydroResult | None = None

    # ── public ─────────────────────────────────────────────────────────────

    def fit(self) -> EMHydroResult:
        """Run all petrophysical transforms and return :class:`EMHydroResult`."""
        model = self.resistivity_model
        cfg = self.config
        petro = cfg.petro

        # ── pass 1: water table ────────────────────────────────────────────
        wt = self._compute_water_table()

        # ── pass 2: porosity and saturation ────────────────────────────────
        phi_map, Sw_map = self._compute_porosity_saturation(wt)

        # ── pass 3: hydraulic conductivity ─────────────────────────────────
        K_map = self._compute_hydraulic_K(phi_map)

        # ── pass 4: depth-cell thicknesses ─────────────────────────────────
        h = _cell_thicknesses(model.z_centers)  # (n_z,)

        # ── pass 5: profile quantities ─────────────────────────────────────
        T_arr, Sc_arr, Su_arr = self._compute_transmissivity_storativity(
            K_map, phi_map, wt, h
        )
        TR_arr, S_arr = self._compute_dar_zarrouk(h)
        tds_val = float(rho_w_to_tds(cfg.rho_w))

        result = EMHydroResult(
            resistivity_model=model,
            config=cfg,
            porosity=phi_map,
            saturation=Sw_map,
            hydraulic_K=K_map,
            water_table=wt,
            transmissivity=T_arr,
            storativity_confined=Sc_arr,
            storativity_unconfined=Su_arr,
            dar_zarrouk_TR=TR_arr,
            dar_zarrouk_S=S_arr,
            tds=tds_val,
            method_tag=self.method_tag or model.method,
            metadata={
                "petro": repr(petro),
                "rho_w": cfg.rho_w,
                "porosity_prior": cfg.porosity_prior,
                "fracture_depth_m": cfg.fracture_depth_m,
                "source_method": model.method,
                "rms": model.rms,
            },
        )
        self._result = result
        return result

    # ── private helpers ────────────────────────────────────────────────────

    def _compute_water_table(self) -> np.ndarray:
        """One-pass water-table depth estimate per column."""
        model = self.resistivity_model
        cfg = self.config
        petro = cfg.petro
        wt = np.full(model.n_x, np.nan)
        for ix in range(model.n_x):
            col = model.rho_2d[:, ix]
            depth = water_table_from_profile(
                col,
                model.z_centers,
                petro
                if isinstance(petro, ArchieModel)
                else _archie_from_ws(petro),
                rho_w=cfg.rho_w,
                Sw_threshold=cfg.Sw_water_table_threshold,
                min_depth=cfg.min_wt_search_depth,
            )
            if depth is not None:
                wt[ix] = depth
        return wt

    def _compute_porosity_saturation(self, wt: np.ndarray) -> tuple:
        """Two-zone porosity/saturation per cell.

        Saturated zone (z ≥ water table): assume Sw=1, invert for φ.
        Vadose zone   (z <  water table): use φ_prior, invert for Sw.
        """
        model = self.resistivity_model
        cfg = self.config
        petro = cfg.petro
        archie = (
            petro if isinstance(petro, ArchieModel) else _archie_from_ws(petro)
        )

        phi_map = np.full(model.rho_2d.shape, cfg.porosity_prior)
        Sw_map = np.ones(model.rho_2d.shape)

        for ix in range(model.n_x):
            wt_depth = wt[ix] if np.isfinite(wt[ix]) else 0.0
            for iz, z in enumerate(model.z_centers):
                rho_log = float(model.rho_2d[iz, ix])
                rho = 10.0**rho_log
                if not np.isfinite(rho) or rho <= 0:
                    continue
                if z >= wt_depth:
                    # saturated — invert for porosity
                    phi = float(archie.porosity(rho, 1.0, cfg.rho_w))
                    phi = float(
                        np.clip(phi, 1e-4, min(0.99, 3.0 * cfg.porosity_prior))
                    )
                    phi_map[iz, ix] = phi
                    Sw_map[iz, ix] = 1.0
                else:
                    # vadose — invert for saturation
                    phi_map[iz, ix] = cfg.porosity_prior
                    Sw = float(
                        archie.saturation(rho, cfg.porosity_prior, cfg.rho_w)
                    )
                    Sw_map[iz, ix] = float(np.clip(Sw, 0.0, 1.0))

        return phi_map, Sw_map

    def _compute_hydraulic_K(self, phi_map: np.ndarray) -> np.ndarray:
        """Hydraulic conductivity map.

        Below ``fracture_depth_m``: Kozeny-Carman from porosity.
        At or below ``fracture_depth_m``: cubic-law fracture K from ρ contrast.
        A logistic sigmoid over ±20 m around ``fracture_depth_m`` blends the
        two models, removing the abrupt step that caused a sharp colour band.
        """
        _BLEND_HALF_WIDTH = 20.0  # metres; sigmoid scale = half_width / 4
        model = self.resistivity_model
        cfg = self.config
        K_map = np.zeros(model.rho_2d.shape)

        for iz, z in enumerate(model.z_centers):
            rho_row = 10.0 ** model.rho_2d[iz, :]
            K_kc = kozeny_carman_K(
                phi_map[iz, :],
                d50_m=cfg.d50_m,
                C=cfg.kozeny_C,
                T=cfg.kozeny_tortuosity,
            )
            if cfg.fracture_depth_m is None:
                K_map[iz, :] = K_kc
            else:
                K_frac = fractured_zone_K(
                    rho_row,
                    rho_matrix=cfg.fracture_rho_matrix,
                )
                # sigmoid weight: 0 = pure Kozeny-Carman, 1 = pure fracture-K
                # transitions from ~0.07 to ~0.93 over the ±20 m blend window
                scale = _BLEND_HALF_WIDTH / 4.0
                w = 1.0 / (1.0 + np.exp(-(z - cfg.fracture_depth_m) / scale))
                K_map[iz, :] = (1.0 - w) * K_kc + w * K_frac
        return K_map

    def _compute_transmissivity_storativity(
        self,
        K_map: np.ndarray,
        phi_map: np.ndarray,
        wt: np.ndarray,
        h: np.ndarray,
    ) -> tuple:
        """Integrate K × dz over the saturated zone per column."""
        model = self.resistivity_model
        cfg = self.config
        n_x = model.n_x
        T_arr = np.zeros(n_x)
        Sc_arr = np.zeros(n_x)
        Su_arr = np.zeros(n_x)

        for ix in range(n_x):
            wt_depth = wt[ix] if np.isfinite(wt[ix]) else 0.0
            sat_mask = model.z_centers >= wt_depth
            if not sat_mask.any():
                continue
            K_sat = K_map[sat_mask, ix]
            phi_sat = phi_map[sat_mask, ix]
            h_sat = h[sat_mask]
            b_total = float(np.sum(h_sat))
            T_arr[ix] = float(np.sum(K_sat * h_sat))
            Sc_arr[ix] = b_total * cfg.specific_storage
            Su_arr[ix] = float(np.nanmean(phi_sat))  # specific yield ≈ φ

        return T_arr, Sc_arr, Su_arr

    def _compute_dar_zarrouk(self, h: np.ndarray) -> tuple:
        """Dar-Zarrouk parameters per column.

        TR = Σ ρᵢ hᵢ   (transverse resistance, Ω·m²)
        S  = Σ hᵢ/ρᵢ   (longitudinal conductance, siemens)
        """
        model = self.resistivity_model
        rho_2d = 10.0**model.rho_2d  # (n_z, n_x)
        h_col = h[:, np.newaxis]  # broadcast over x

        TR = np.sum(rho_2d * h_col, axis=0)
        S = np.sum(h_col / np.maximum(rho_2d, 1e-10), axis=0)
        return TR, S


# ─────────────────────────────────────────────────────────────────────────────
# Internal helpers
# ─────────────────────────────────────────────────────────────────────────────


def _cell_thicknesses(z_centers: np.ndarray) -> np.ndarray:
    """Estimate cell thicknesses from centre positions (metres)."""
    z = np.asarray(z_centers, dtype=float)
    if z.size == 1:
        return np.array([max(1.0, float(z[0]))])
    mids = 0.5 * (z[:-1] + z[1:])
    first = max(0.0, z[0] - (mids[0] - z[0]))
    last = z[-1] + (z[-1] - mids[-1])
    edges = np.r_[first, mids, last]
    return np.diff(edges)


def _station_name(model: ResistivityModel, ix: int) -> str:
    if model.station_names and ix < len(model.station_names):
        return str(model.station_names[ix])
    return f"S{ix:03d}"


def _archie_from_ws(ws: WaxmanSmitsModel) -> ArchieModel:
    """Approximate ArchieModel from WaxmanSmitsModel for water-table search."""
    return ArchieModel(m=ws.m, n=ws.n, a=ws.a)


def _replace_config(cfg: PetrophysicalConfig, **kwargs) -> PetrophysicalConfig:
    """Return a new PetrophysicalConfig with selected fields overridden."""
    import dataclasses

    return dataclasses.replace(cfg, **kwargs)
