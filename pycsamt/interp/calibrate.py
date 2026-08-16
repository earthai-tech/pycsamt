# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""ModelCalibrator — constrain a 2-D EM model with borehole ground truth.

Implements the pseudo-stratigraphic New Model (NM) construction
introduced in Kouadio et al. (2022) and generalized here to work with
any EM inversion output (Occam2D, ModEM, AI-based 1-D/2-D).

Algorithm
---------
For each station column of the Calculated Resistivity Model (CRM):

1. **Soft replace** — for every depth cell whose CRM resistivity is
   within *ptol* (default 10 %) of a borehole TRES value, replace the
   calculated value with the true value:

   .. math::

       \\xi_{\\mathrm{se}} = \\frac{|\\rho_{\\mathrm{crm}} - \\rho_k|}{\\rho_k}
       \\leq p_{\\mathrm{tol}}

2. **Autolayer** — for cells not matched in Step 1, assign the nearest
   rock from the geological database :math:`\\Gamma` in log-space.

The **misfit G (%)** between the CRM and the NM is:

.. math::

    G(\\%) = 100 \\times
    \\sqrt{\\frac{\\sum_i (\\mathrm{NM}_i - \\mathrm{CRM}_i)^2}
                 {\\sum_i \\mathrm{CRM}_i^2}}

References
----------
.. [1] Kouadio, K.L. et al. (2022). pyCSAMT: An alternative Python
   toolbox for groundwater exploration using controlled source
   audio-frequency magnetotelluric. *Journal of Applied Geophysics*,
   201, 104647. https://doi.org/10.1016/j.jappgeo.2022.104647

Example
-------
>>> from pycsamt.interp import ResistivityModel, Borehole, ModelCalibrator
>>> model = ResistivityModel.from_occam2d(result)
>>> bh = Borehole.from_csv("boreholes/Bo.csv", name="Bo", x=1050.0)
>>> cal = ModelCalibrator(ptol=0.10).fit(model, [bh])
>>> nm = cal.calibrated_model()
>>> logs = cal.stratigraphic_logs()
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np

from ..api.property import PyCSAMTObject
from ..geology.borehole import Borehole
from ..geology.lithology import RockDatabase, StratigraphicLog
from ._base import ResistivityModel

__all__ = ["ModelCalibrator"]


class ModelCalibrator(PyCSAMTObject):
    """Constrain a 2-D EM resistivity model with borehole TRES values.

    Parameters
    ----------
    ptol : float, default 0.10
        Soft-error threshold — fractional tolerance between a CRM
        resistivity value and a TRES value for accepting a replace.
        The paper uses 10 % (``ptol = 0.10``).
    max_borehole_distance : float, default 500.0
        Maximum profile-distance (metres) from a station to the nearest
        borehole.  Stations farther than this threshold are calibrated
        using the rock database only (no TRES constraint).
    db : RockDatabase, optional
        Rock physics database for autolayer assignment.
        Defaults to :meth:`~RockDatabase.default`.
    verbose : bool, default True
    """

    def __init__(
        self,
        ptol: float = 0.10,
        *,
        max_borehole_distance: float = 500.0,
        db: RockDatabase | None = None,
        verbose: bool = True,
    ) -> None:
        self.ptol = float(ptol)
        self.max_borehole_distance = float(max_borehole_distance)
        self.db = db if db is not None else RockDatabase.default()
        self.verbose = verbose

        self._model: ResistivityModel | None = None
        self._boreholes: list[Borehole] = []
        self._nm_rho_2d: np.ndarray | None = None
        self._misfit_map: np.ndarray | None = None
        self._is_fitted: bool = False

    # ------------------------------------------------------------------
    # Fit
    # ------------------------------------------------------------------

    def fit(
        self,
        model: ResistivityModel,
        boreholes: Sequence[Borehole] | None = None,
    ) -> ModelCalibrator:
        """Calibrate *model* against *boreholes*.

        Parameters
        ----------
        model : ResistivityModel
            The Calculated Resistivity Model (CRM) from any EM inversion.
        boreholes : list of Borehole, optional
            Ground-truth borehole logs.  If ``None`` or empty, the
            calibration falls back to rock-database autolayer assignment
            for every cell.

        Returns
        -------
        self
        """
        self._model = model
        self._boreholes = list(boreholes) if boreholes else []

        crm = model.rho_2d.copy()
        nm = crm.copy()

        for ix, x_sta in enumerate(model.x_centers):
            bh = self._nearest_borehole(x_sta)
            col_crm = crm[:, ix]
            if (
                bh is not None
                and abs(bh.x - x_sta) <= self.max_borehole_distance
            ):
                nm[:, ix] = self._calibrate_column(
                    col_crm, model.z_centers, bh
                )
            else:
                nm[:, ix] = self._autolayer_column(col_crm)

        self._nm_rho_2d = nm
        self._misfit_map = self._compute_misfit(crm, nm)
        self._is_fitted = True

        if self.verbose:
            g_mean = float(np.nanmean(self._misfit_map))
            print(
                f"  ModelCalibrator: fitted {model.n_x} columns, "
                f"{len(self._boreholes)} borehole(s), "
                f"mean misfit G = {g_mean:.2f} %"
            )
        return self

    # ------------------------------------------------------------------
    # Results
    # ------------------------------------------------------------------

    def calibrated_model(self) -> ResistivityModel:
        """Return a :class:`ResistivityModel` containing the NM.

        Raises
        ------
        RuntimeError
            If :meth:`fit` has not been called.
        """
        self._check_fitted()
        m = self._model
        return ResistivityModel(
            x_centers=m.x_centers.copy(),
            z_centers=m.z_centers.copy(),
            rho_2d=self._nm_rho_2d.copy(),
            station_x=m.station_x.copy(),
            station_names=list(m.station_names),
            method=m.method + "+calibrated",
            rms=m.rms,
        )

    def misfit_map(self) -> np.ndarray:
        """Return the per-column misfit G (%), shape (n_z, n_x).

        Values range 0 – 100 %.  Low values indicate good CRM–TRES
        agreement; high values flag where the inversion needed the
        most correction.
        """
        self._check_fitted()
        return self._misfit_map.copy()

    def stratigraphic_logs(
        self,
        *,
        db: RockDatabase | None = None,
        model: str = "nm",
        merge_tolerance: float = 0.2,
    ) -> list[StratigraphicLog]:
        """Build a :class:`~pycsamt.geology.lithology.StratigraphicLog`
        for every station.

        Parameters
        ----------
        db : RockDatabase, optional
            Override the classifier database.
        model : {'nm', 'crm'}
            Which model to classify: the calibrated NM or the original CRM.
        merge_tolerance : float
            Log-space merging tolerance passed to
            :meth:`~StratigraphicLog.from_column`.

        Returns
        -------
        list of StratigraphicLog
        """
        self._check_fitted()
        _db = db if db is not None else self.db
        m = self._model
        rho_2d = self._nm_rho_2d if model == "nm" else m.rho_2d

        logs: list[StratigraphicLog] = []
        for _i, (sta_x, sta_name) in enumerate(
            zip(m.station_x, m.station_names)
        ):
            ix = int(np.argmin(np.abs(m.x_centers - sta_x)))
            col = rho_2d[:, ix]
            log = StratigraphicLog.from_column(
                station_name=sta_name,
                x=float(sta_x),
                z_centers=m.z_centers,
                rho_log10=col,
                db=_db,
                merge_tolerance=merge_tolerance,
            )
            logs.append(log)

        # If no station list was provided, fall back to all x-columns
        if not logs:
            for ix, x_c in enumerate(m.x_centers):
                col = rho_2d[:, ix]
                log = StratigraphicLog.from_column(
                    station_name=f"S{ix:03d}",
                    x=float(x_c),
                    z_centers=m.z_centers,
                    rho_log10=col,
                    db=_db,
                    merge_tolerance=merge_tolerance,
                )
                logs.append(log)

        return logs

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _nearest_borehole(self, x: float) -> Borehole | None:
        if not self._boreholes:
            return None
        dists = [abs(bh.x - x) for bh in self._boreholes]
        return self._boreholes[int(np.argmin(dists))]

    def _calibrate_column(
        self,
        col_crm: np.ndarray,
        z_centers: np.ndarray,
        bh: Borehole,
    ) -> np.ndarray:
        """Apply the two-step NM construction to a single column."""
        tres_col = bh.tres_column(z_centers)  # TRES at each depth, Ω·m
        rho_crm = 10.0**col_crm  # linear Ω·m
        nm = col_crm.copy()  # initialise NM = CRM

        for i in range(len(z_centers)):
            if np.isnan(col_crm[i]):
                continue
            tres_i = tres_col[i]
            if np.isnan(tres_i):
                # No TRES at this depth — autolayer from database
                nm[i] = self._autolayer_single(rho_crm[i])
                continue

            # Step 1: soft replace
            xi_se = abs(rho_crm[i] - tres_i) / max(tres_i, 1e-6)
            if xi_se <= self.ptol:
                nm[i] = np.log10(max(tres_i, 1e-6))
            else:
                # Step 2: autolayer (nearest rock in Γ consistent with TRES)
                nm[i] = self._autolayer_single(tres_i)

        return nm

    def _autolayer_single(self, rho_ohm_m: float) -> float:
        """Return log10-rho of the database rock nearest to *rho_ohm_m*."""
        entry = self.db.classify(float(rho_ohm_m))
        return entry.log_rho_mid

    def _autolayer_column(self, col_crm: np.ndarray) -> np.ndarray:
        nm = col_crm.copy()
        for i, v in enumerate(col_crm):
            if not np.isnan(v):
                nm[i] = self._autolayer_single(float(10.0**v))
        return nm

    @staticmethod
    def _compute_misfit(crm: np.ndarray, nm: np.ndarray) -> np.ndarray:
        """Per-column G (%) misfit, broadcast to (n_z, n_x)."""
        diff_sq = (nm - crm) ** 2
        crm_sq = crm**2
        with np.errstate(divide="ignore", invalid="ignore"):
            g_col = 100.0 * np.sqrt(
                np.nansum(diff_sq, axis=0)
                / np.maximum(np.nansum(crm_sq, axis=0), 1e-12)
            )
        # Broadcast scalar per-column misfit back to a 2-D map
        g_map = np.tile(g_col, (crm.shape[0], 1))
        return g_map

    def _check_fitted(self) -> None:
        if not self._is_fitted:
            raise RuntimeError("Call fit() before accessing results.")

    # ------------------------------------------------------------------

    def __repr__(self) -> str:
        status = "fitted" if self._is_fitted else "unfitted"
        return f"ModelCalibrator(ptol={self.ptol}, {status})"
