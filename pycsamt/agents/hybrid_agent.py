# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
r"""
pycsamt.agents.hybrid_agent
============================

:class:`HybridInversionAgent` wraps
:class:`~pycsamt.ai.inversion.HybridInverter1D`,
:class:`~pycsamt.ai.inversion.HybridInverter2D`, and
:class:`~pycsamt.ai.inversion.HybridInverter3D`.

Two-stage workflow:

**Stage 1** — pre-trained supervised AI inverter
(EMInverter1D / EMInverter2D / GCNInverter3D) is
applied to produce a physically plausible initial
model in milliseconds.

**Stage 2** — physics-informed gradient descent
(same Wait 1954 loss as PINNInverter) refines the
Stage-1 model.  Convergence is faster and more
reliable than starting from a random initialisation.
"""

from __future__ import annotations

import time
from typing import Any

import numpy as np

from ._base import AgentResult, BaseAgent
from .pinn_agent import (
    _plot_loss_curves,
    _plot_pinn_section,
    _rms_from_residuals,
)

_HYBRID_SYSTEM = """\
You are an expert in hybrid AI + physics-informed
inversion for MT/CSAMT geophysics.
Given a hybrid inversion result write 4-5 sentences:
1. Compare Stage-1 (AI) and Stage-2 (physics) RMS.
2. Describe how much the physics step improved the fit.
3. State the recovered resistivity structure.
4. Flag stations where Stage-2 failed to improve on
   Stage-1 or where residuals remain high.
5. Recommend whether to retrain the AI component,
   run more physics iterations, or proceed to 2-D.
Reply in plain scientific English.
"""


class HybridInversionAgent(BaseAgent):
    r"""Two-stage AI + physics MT inversion.

    Stage 1 applies a pre-trained supervised AI
    inverter to obtain a starting model.
    Stage 2 refines it with physics-informed Adam
    gradient descent.

    Parameters
    ----------
    dim : {1, 2, 3}
        Dimensionality.  Default ``1``.
    max_iter : int
        Physics refinement iterations (Stage 2).
        Default ``200``.
    smoothness_weight : float
        Vertical regularisation weight.
        Default ``0.005``.
    lateral_weight : float
        Lateral smoothness weight (2-D only).
        Default ``0.005``.
    graph_weight : float
        Graph-Laplacian weight (3-D only).
        Default ``0.005``.
    radius : float
        Edge radius in metres for 3-D graph.
        Default ``5000.0``.
    lr : float
        Adam learning rate for Stage 2.
        Default ``5e-3``.
    solver : {"mt1d", "csamt1d"}
        Physics solver.  Default ``"mt1d"``.
    comp : str
        Impedance component (1-D only).
        Default ``"xy"``.
    n_freqs : int
        Frequency-grid size fed to the 1-D AI
        inverter.  Default ``32``.
    api_key, model, llm_provider :
        LLM configuration (optional).

    Input keys
    ----------
    ``sites`` / ``path``     observed data
    ``ai_inverter``          fitted AI inverter object
                             or path to checkpoint
    ``checkpoint``           alias for ``ai_inverter``
    ``output_dir``           optional save directory
    ``dim``, ``max_iter``,
    ``smoothness_weight``,
    ``lateral_weight``,
    ``graph_weight``         : optional overrides

    Output data keys
    ----------------
    ``inverter``         fitted HybridInverterXD
    ``section``          ndarray (n_layers, n_stations)
                         Stage-2 log10-rho section
    ``stage1_section``   ndarray — Stage-1 section
    ``models``           list of LayeredModel (1-D)
    ``stage1_models``    list of LayeredModel (1-D)
    ``n_stations``       int
    ``rms_per_station``  dict {station: float}
    ``rms_global``       float  (Stage-2)
    ``rms_stage1``       float  (Stage-1 for comparison)
    ``convergence_df``   pandas.DataFrame or None
    ``residuals_df``     pandas.DataFrame or None
    ``figures``          dict
    ``figure_paths``     dict

    Examples
    --------
    >>> from pycsamt.ai.inversion import EMInverter1D
    >>> ai = EMInverter1D.load("checkpoint.npz")
    >>> agent = HybridInversionAgent(dim=1, max_iter=100)
    >>> res = agent.execute(
    ...     {
    ...         "path": "/data/L22PLT",
    ...         "ai_inverter": ai,
    ...     }
    ... )
    >>> res["rms_global"]
    0.14
    """

    SYSTEM_PROMPT = _HYBRID_SYSTEM

    def __init__(
        self,
        *,
        dim: int = 1,
        max_iter: int = 200,
        smoothness_weight: float = 0.005,
        lateral_weight: float = 0.005,
        graph_weight: float = 0.005,
        radius: float = 5000.0,
        lr: float = 5e-3,
        solver: str = "mt1d",
        comp: str = "xy",
        n_freqs: int = 32,
        api_key: str | None = None,
        model: str | None = None,
        llm_provider: str = "claude",
    ) -> None:
        super().__init__(
            "HybridInversionAgent",
            api_key=api_key,
            model=model,
            llm_provider=llm_provider,
            section_preset="inversion",
        )
        if dim not in (1, 2, 3):
            raise ValueError(f"dim must be 1, 2, or 3; got {dim!r}.")
        self.dim = dim
        self.max_iter = max_iter
        self.smoothness_weight = smoothness_weight
        self.lateral_weight = lateral_weight
        self.graph_weight = graph_weight
        self.radius = radius
        self.lr = lr
        self.solver = solver
        self.comp = comp
        self.n_freqs = n_freqs

    # ── public entry point ────────────────────────

    def execute(self, input_data: dict[str, Any]) -> AgentResult:
        self._last_cost = 0.0
        t0 = time.time()
        warns: list[str] = []

        try:
            from ..backends import get_backend_instance

            if get_backend_instance() is None:
                raise ImportError("No DL backend available.")
        except ImportError as exc:
            return AgentResult.failed(
                f"Hybrid inversion requires PyTorch or TensorFlow: {exc}",
                hint=("pip install torch\npip install tensorflow"),
                elapsed=time.time() - t0,
            )

        dim = int(input_data.get("dim", self.dim))
        max_iter = int(input_data.get("max_iter", self.max_iter))
        sw = float(
            input_data.get(
                "smoothness_weight",
                self.smoothness_weight,
            )
        )
        lw = float(
            input_data.get(
                "lateral_weight",
                self.lateral_weight,
            )
        )
        gw = float(
            input_data.get(
                "graph_weight",
                self.graph_weight,
            )
        )
        output_dir = input_data.get("output_dir")

        ai_inv_raw = input_data.get("ai_inverter") or input_data.get(
            "checkpoint"
        )
        if ai_inv_raw is None:
            return AgentResult.failed(
                "No 'ai_inverter' or 'checkpoint' "
                "in input_data.  Provide a fitted "
                "AI inverter or a checkpoint path.",
                elapsed=time.time() - t0,
            )

        sites_raw = input_data.get("sites") or input_data.get("path")
        if sites_raw is None:
            return AgentResult.failed(
                "No 'sites' or 'path' in input_data.",
                elapsed=time.time() - t0,
            )
        try:
            from ..emtools._core import ensure_sites

            sites = ensure_sites(sites_raw, verbose=0)
        except Exception as exc:
            return AgentResult.failed(
                str(exc),
                elapsed=time.time() - t0,
            )

        try:
            (
                inv,
                mat,
                s1_mat,
                conv_df,
                res_df,
                s1_res_df,
            ) = self._run(
                dim,
                sites,
                ai_inv_raw,
                max_iter,
                sw,
                lw,
                gw,
                warns,
            )
        except Exception as exc:
            return AgentResult.failed(
                f"Hybrid fitting failed: {exc}",
                hint=(
                    "Check that the AI inverter "
                    "is trained for this dim. "
                    "Reduce max_iter if OOM."
                ),
                elapsed=time.time() - t0,
            )

        station_names = inv.stations
        n_st = inv.n_sites

        rms_per, rms_global = _rms_from_residuals(res_df)
        _, rms_s1 = _rms_from_residuals(s1_res_df)

        figures: dict[str, Any] = {}
        fig_paths: dict[str, str] = {}

        n_layers = getattr(
            inv,
            "n_layers",
            mat.shape[0] if mat is not None else 10,
        )
        depth_max = getattr(inv, "depth_max", 2000.0)

        if mat is not None and n_st > 0:
            ths = np.logspace(
                np.log10(max(depth_max / 100, 50)),
                np.log10(depth_max),
                n_layers - 1,
            )
            depths_km = np.concatenate([[0.0], np.cumsum(ths)]) / 1000.0

            try:
                fig_s2 = _plot_pinn_section(
                    mat,
                    station_names,
                    n_layers,
                    depths_km,
                    title=(f"Hybrid-{dim}D Stage-2 section"),
                )
                if fig_s2 is not None:
                    figures["section"] = fig_s2
                    p = self._save_figure(
                        fig_s2,
                        output_dir,
                        f"hybrid{dim}d_section",
                        warnings_list=warns,
                    )
                    if p:
                        fig_paths["section"] = p
            except Exception as exc:
                warns.append(f"Stage-2 section plot: {exc}")

            if s1_mat is not None:
                try:
                    fig_s1 = _plot_pinn_section(
                        s1_mat,
                        station_names,
                        n_layers,
                        depths_km,
                        title=(f"Hybrid-{dim}D Stage-1 (AI) section"),
                    )
                    if fig_s1 is not None:
                        figures["stage1_section"] = fig_s1
                        p = self._save_figure(
                            fig_s1,
                            output_dir,
                            f"hybrid{dim}d_s1_section",
                            warnings_list=warns,
                        )
                        if p:
                            fig_paths["stage1_section"] = p
                except Exception as exc:
                    warns.append(f"Stage-1 section plot: {exc}")

        if conv_df is not None and len(conv_df):
            try:
                fig_c = _plot_loss_curves(
                    conv_df,
                    title=(f"Hybrid-{dim}D Stage-2 convergence"),
                )
                if fig_c is not None:
                    figures["convergence"] = fig_c
                    p = self._save_figure(
                        fig_c,
                        output_dir,
                        f"hybrid{dim}d_convergence",
                        warnings_list=warns,
                    )
                    if p:
                        fig_paths["convergence"] = p
            except Exception as exc:
                warns.append(f"Convergence plot: {exc}")

        interp: str | None = None
        if self.api_key and n_st > 0:
            rms_s = f"{rms_global:.3f}" if not np.isnan(rms_global) else "N/A"
            rms_s1_s = f"{rms_s1:.3f}" if not np.isnan(rms_s1) else "N/A"
            n_hi = sum(1 for v in rms_per.values() if v > 0.5)
            prompt = (
                f"Hybrid-{dim}D MT inversion:\n"
                f"  Stations: {n_st}\n"
                f"  Stage-1 (AI) global RMS: "
                f"{rms_s1_s}\n"
                f"  Stage-2 (physics) global RMS: "
                f"{rms_s}\n"
                f"  High-RMS after Stage-2 "
                f"(>0.5): {n_hi}\n"
                f"  max_iter={max_iter}\n"
                "Evaluate both stages and recommend "
                "next steps."
            )
            interp = self.query_llm(prompt, max_tokens=300)

        models: list = []
        stage1_models: list = []
        if dim == 1:
            try:
                models = inv.predict()
            except Exception as exc:
                warns.append(f"predict(): {exc}")
            try:
                stage1_models = inv.stage1_models()
            except Exception as exc:
                warns.append(f"stage1_models(): {exc}")

        elapsed = time.time() - t0
        rms_str = (
            f"RMS {rms_global:.3f}" if not np.isnan(rms_global) else "RMS N/A"
        )
        return AgentResult(
            status=("success" if n_st > 0 else "needs_review"),
            summary=(
                f"Hybrid-{dim}D: {n_st} stations. "
                f"{rms_str} (Stage-2). "
                f"{len(figures)} figure(s)."
            ),
            data={
                "inverter": inv,
                "section": mat,
                "stage1_section": s1_mat,
                "models": models,
                "stage1_models": stage1_models,
                "n_stations": n_st,
                "rms_per_station": rms_per,
                "rms_global": rms_global,
                "rms_stage1": rms_s1,
                "convergence_df": conv_df,
                "residuals_df": res_df,
                "figures": figures,
                "figure_paths": fig_paths,
            },
            warnings=warns,
            llm_interpretation=interp,
            elapsed_seconds=elapsed,
            cost_estimate_usd=self._last_cost,
        )

    # ── private runners ───────────────────────────

    def _run(
        self,
        dim: int,
        sites: Any,
        ai_inv_raw: Any,
        max_iter: int,
        sw: float,
        lw: float,
        gw: float,
        warns: list[str],
    ):
        if dim == 1:
            return self._run_1d(
                sites,
                ai_inv_raw,
                max_iter,
                sw,
                warns,
            )
        if dim == 2:
            return self._run_2d(
                sites,
                ai_inv_raw,
                max_iter,
                sw,
                lw,
                warns,
            )
        return self._run_3d(
            sites,
            ai_inv_raw,
            max_iter,
            sw,
            gw,
            warns,
        )

    def _run_1d(
        self,
        sites,
        ai_inv_raw,
        max_iter,
        sw,
        warns,
    ):
        from ..ai.inversion.hybrid1d import (
            HybridInverter1D,
        )

        inv = HybridInverter1D(
            sites,
            ai_inv_raw,
            solver=self.solver,
            max_iter=max_iter,
            smoothness_weight=sw,
            lr=self.lr,
            comp=self.comp,
            n_freqs=self.n_freqs,
            verbose=0,
        )
        inv.fit(verbose=False)

        n_layers = inv._ai_inv.n_layers
        n_st = inv.n_sites

        # Stage-2 section matrix
        mat = np.full((n_layers, n_st), np.nan)
        for si, res in enumerate(inv._stage2):
            lr_arr = res["log_rho"]
            n = min(len(lr_arr), n_layers)
            mat[:n, si] = lr_arr[:n]

        # Stage-1 section matrix
        s1_mat = np.full((n_layers, n_st), np.nan)
        for si, res in enumerate(inv._stage1):
            lr_arr = res["log_rho"]
            n = min(len(lr_arr), n_layers)
            s1_mat[:n, si] = lr_arr[:n]

        try:
            conv_df = inv.convergence_curves()
        except Exception:
            conv_df = None
        try:
            res_df = inv.residuals(stage=2)
        except Exception as exc:
            warns.append(f"residuals(2): {exc}")
            res_df = None
        try:
            s1_res_df = inv.residuals(stage=1)
        except Exception:
            s1_res_df = None
        return inv, mat, s1_mat, conv_df, res_df, s1_res_df

    def _run_2d(
        self,
        sites,
        ai_inv_raw,
        max_iter,
        sw,
        lw,
        warns,
    ):
        from ..ai.inversion.hybrid2d import (
            HybridInverter2D,
        )

        inv = HybridInverter2D(
            sites,
            ai_inv_raw,
            solver=self.solver,
            max_iter=max_iter,
            smoothness_weight=sw,
            lateral_weight=lw,
            lr=self.lr,
            verbose=0,
        )
        inv.fit(verbose=False)

        try:
            mat = inv.resistivity_section(stage=2)
        except Exception:
            try:
                mat = inv.resistivity_section()
            except Exception:
                mat = None
        try:
            s1_mat = inv.resistivity_section(stage=1)
        except Exception:
            s1_mat = None
        try:
            conv_df = inv.convergence_curve()
        except Exception:
            conv_df = None
        try:
            res_df = inv.residuals(stage=2)
        except Exception as exc:
            warns.append(f"residuals(2): {exc}")
            res_df = None
        try:
            s1_res_df = inv.residuals(stage=1)
        except Exception:
            s1_res_df = None
        return inv, mat, s1_mat, conv_df, res_df, s1_res_df

    def _run_3d(
        self,
        sites,
        ai_inv_raw,
        max_iter,
        sw,
        gw,
        warns,
    ):
        from ..ai.inversion.hybrid3d import (
            HybridInverter3D,
        )

        inv = HybridInverter3D(
            sites,
            ai_inv_raw,
            max_iter=max_iter,
            smoothness_weight=sw,
            graph_weight=gw,
            radius=self.radius,
            lr=self.lr,
            verbose=0,
        )
        inv.fit(verbose=False)

        try:
            mat = inv.resistivity_volume(stage=2)
        except Exception:
            try:
                mat = inv.resistivity_volume()
            except Exception:
                mat = None
        try:
            s1_mat = inv.resistivity_volume(stage=1)
        except Exception:
            s1_mat = None
        try:
            conv_df = inv.convergence_curve()
        except Exception:
            conv_df = None
        try:
            res_df = inv.residuals(stage=2)
        except Exception as exc:
            warns.append(f"residuals(2): {exc}")
            res_df = None
        try:
            s1_res_df = inv.residuals(stage=1)
        except Exception:
            s1_res_df = None
        return inv, mat, s1_mat, conv_df, res_df, s1_res_df


__all__ = ["HybridInversionAgent"]
