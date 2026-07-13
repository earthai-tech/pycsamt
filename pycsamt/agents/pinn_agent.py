# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
r"""
pycsamt.agents.pinn_agent
=========================

:class:`PINNInversionAgent` wraps
:class:`~pycsamt.ai.inversion.PINNInverter1D`,
:class:`~pycsamt.ai.inversion.PINNInverter2D`, and
:class:`~pycsamt.ai.inversion.PINNInverter3D`.

Physics-informed optimisation requires no labelled
training data.  Model parameters are updated by
gradient descent on the Wait (1954) MT forward-
physics loss.
"""

from __future__ import annotations

import time
from typing import Any

import numpy as np

from ._base import AgentResult, BaseAgent

_PINN_SYSTEM = """\
You are an expert in physics-informed neural-network
inversion for MT/CSAMT geophysics.
Given a PINN inversion result write 4-5 sentences:
1. State dimensionality (1-D/2-D/3-D) and convergence.
2. Report final RMS (log10 rho-ohm-m) and data fit.
3. Describe the resistivity structure recovered.
4. Flag stations or regions with high residuals.
5. Recommend adjustments to epochs, regularisation,
   or whether to switch to a hybrid approach.
Reply in plain scientific English.
"""

_DEF_EPOCHS: dict[int, int] = {1: 500, 2: 300, 3: 300}


class PINNInversionAgent(BaseAgent):
    r"""PINN-based MT inversion without labelled data.

    Optimises a layered Earth by minimising a
    physics-informed loss via Adam gradient descent.
    Supports 1-D per-station, joint 2-D profile, and
    quasi-3-D graph-coupled inversion.

    Parameters
    ----------
    dim : {1, 2, 3}
        Dimensionality.  Default ``1``.
    n_layers : int
        Number of layers including the halfspace.
        Default ``10``.
    depth_max : float
        Maximum investigation depth in metres.
        Default ``2000.0``.
    smoothness_weight : float
        Vertical regularisation weight.
        Default ``0.01``.
    lateral_weight : float
        Lateral smoothness weight (2-D only).
        Default ``0.005``.
    graph_weight : float
        Graph-Laplacian spatial weight (3-D only).
        Default ``0.005``.
    radius : float
        Edge radius in metres for the 3-D graph.
        Default ``5000.0``.
    epochs : int or None
        Adam iterations.  ``None`` uses 500 for 1-D
        and 300 for 2-D / 3-D.
    lr : float
        Adam learning rate.  Default ``1e-2``.
    solver : {"mt1d", "csamt1d"}
        Physics solver.  Default ``"mt1d"``.
    comp : str
        Impedance component (1-D only).
        Default ``"xy"``.
    api_key, model, llm_provider :
        LLM configuration (optional).

    Input keys
    ----------
    ``sites`` / ``path``   observed data
    ``output_dir``         optional figure/save dir
    ``dim``, ``epochs``, ``n_layers`` : overrides

    Output data keys
    ----------------
    ``inverter``        fitted inverter object
    ``section``         ndarray (n_layers, n_stations)
                        log10-rho section matrix
    ``models``          list of LayeredModel (1-D)
    ``n_stations``      int
    ``rms_per_station`` dict {station: float}
    ``rms_global``      float
    ``loss_df``         pandas.DataFrame or None
    ``residuals_df``    pandas.DataFrame or None
    ``figures``         dict
    ``figure_paths``    dict

    Examples
    --------
    >>> agent = PINNInversionAgent(
    ...     dim=1, n_layers=10, epochs=200
    ... )
    >>> res = agent.execute(
    ...     {"path": "/data/L22PLT"}
    ... )
    >>> res["rms_global"]
    0.18
    """

    SYSTEM_PROMPT = _PINN_SYSTEM

    def __init__(
        self,
        *,
        dim: int = 1,
        n_layers: int = 10,
        depth_max: float = 2000.0,
        smoothness_weight: float = 0.01,
        lateral_weight: float = 0.005,
        graph_weight: float = 0.005,
        radius: float = 5000.0,
        epochs: int | None = None,
        lr: float = 1e-2,
        solver: str = "mt1d",
        comp: str = "xy",
        api_key: str | None = None,
        model: str | None = None,
        llm_provider: str = "claude",
    ) -> None:
        super().__init__(
            "PINNInversionAgent",
            api_key=api_key,
            model=model,
            llm_provider=llm_provider,
            section_preset="inversion",
        )
        if dim not in (1, 2, 3):
            raise ValueError(f"dim must be 1, 2, or 3; got {dim!r}.")
        self.dim = dim
        self.n_layers = n_layers
        self.depth_max = depth_max
        self.smoothness_weight = smoothness_weight
        self.lateral_weight = lateral_weight
        self.graph_weight = graph_weight
        self.radius = radius
        self.epochs = epochs
        self.lr = lr
        self.solver = solver
        self.comp = comp

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
                f"PINN requires PyTorch or TensorFlow: {exc}",
                hint=("pip install torch\npip install tensorflow"),
                elapsed=time.time() - t0,
            )

        dim = int(input_data.get("dim", self.dim))
        n_layers = int(input_data.get("n_layers", self.n_layers))
        depth_max = float(input_data.get("depth_max", self.depth_max))
        epochs = int(
            input_data.get("epochs", self.epochs) or _DEF_EPOCHS.get(dim, 300)
        )
        output_dir = input_data.get("output_dir")

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
            inv, mat, loss_df, res_df = self._run(
                dim,
                sites,
                n_layers,
                depth_max,
                epochs,
                warns,
            )
        except Exception as exc:
            return AgentResult.failed(
                f"PINN fitting failed: {exc}",
                hint=("Try fewer epochs/layers or check EDI data quality."),
                elapsed=time.time() - t0,
            )

        station_names = inv.stations
        n_st = inv.n_sites

        rms_per, rms_global = _rms_from_residuals(res_df)

        figures: dict[str, Any] = {}
        fig_paths: dict[str, str] = {}

        if mat is not None and n_st > 0:
            try:
                ths = np.logspace(
                    np.log10(max(depth_max / 100, 50)),
                    np.log10(depth_max),
                    n_layers - 1,
                )
                depths_km = np.concatenate([[0.0], np.cumsum(ths)]) / 1000.0
                fig_s = _plot_pinn_section(
                    mat,
                    station_names,
                    n_layers,
                    depths_km,
                    title=(f"PINN-{dim}D resistivity section"),
                )
                if fig_s is not None:
                    figures["section"] = fig_s
                    p = self._save_figure(
                        fig_s,
                        output_dir,
                        f"pinn{dim}d_section",
                        warnings_list=warns,
                    )
                    if p:
                        fig_paths["section"] = p
            except Exception as exc:
                warns.append(f"Section plot: {exc}")

        if loss_df is not None and len(loss_df):
            try:
                fig_c = _plot_loss_curves(
                    loss_df,
                    title=(f"PINN-{dim}D convergence"),
                )
                if fig_c is not None:
                    figures["convergence"] = fig_c
                    p = self._save_figure(
                        fig_c,
                        output_dir,
                        f"pinn{dim}d_convergence",
                        warnings_list=warns,
                    )
                    if p:
                        fig_paths["convergence"] = p
            except Exception as exc:
                warns.append(f"Convergence plot: {exc}")

        interp: str | None = None
        if self.api_key and n_st > 0:
            rms_s = f"{rms_global:.3f}" if not np.isnan(rms_global) else "N/A"
            n_hi = sum(1 for v in rms_per.values() if v > 0.5)
            prompt = (
                f"PINN-{dim}D MT inversion:\n"
                f"  n_layers={n_layers}, "
                f"depth_max={depth_max} m\n"
                f"  Stations: {n_st}\n"
                f"  Global RMS: {rms_s} "
                f"log10(Ohm·m)\n"
                f"  High-RMS stations (>0.5): "
                f"{n_hi}\n"
                f"  Epochs: {epochs}\n"
                "Evaluate and recommend next steps."
            )
            interp = self.query_llm(prompt, max_tokens=250)

        models: list = []
        if dim == 1:
            try:
                models = inv.predict()
            except Exception as exc:
                warns.append(f"predict(): {exc}")

        elapsed = time.time() - t0
        rms_str = (
            f"RMS {rms_global:.3f}" if not np.isnan(rms_global) else "RMS N/A"
        )
        return AgentResult(
            status=("success" if n_st > 0 else "needs_review"),
            summary=(
                f"PINN-{dim}D: {n_st} stations, "
                f"{n_layers} layers. {rms_str}. "
                f"{len(figures)} figure(s)."
            ),
            data={
                "inverter": inv,
                "section": mat,
                "models": models,
                "n_stations": n_st,
                "rms_per_station": rms_per,
                "rms_global": rms_global,
                "loss_df": loss_df,
                "residuals_df": res_df,
                "figures": figures,
                "figure_paths": fig_paths,
            },
            warnings=warns,
            llm_interpretation=interp,
            elapsed_seconds=elapsed,
            cost_estimate_usd=self._last_cost,
        )

    # ── private runners per dimension ─────────────

    def _run(
        self,
        dim: int,
        sites: Any,
        n_layers: int,
        depth_max: float,
        epochs: int,
        warns: list[str],
    ):
        if dim == 1:
            return self._run_1d(
                sites,
                n_layers,
                depth_max,
                epochs,
                warns,
            )
        if dim == 2:
            return self._run_2d(
                sites,
                n_layers,
                depth_max,
                epochs,
                warns,
            )
        return self._run_3d(
            sites,
            n_layers,
            depth_max,
            epochs,
            warns,
        )

    def _run_1d(
        self,
        sites,
        n_layers,
        depth_max,
        epochs,
        warns,
    ):
        from ..ai.inversion.pinn1d import (
            PINNInverter1D,
        )

        inv = PINNInverter1D(
            sites,
            solver=self.solver,
            n_layers=n_layers,
            depth_max=depth_max,
            smoothness_weight=(self.smoothness_weight),
            lr=self.lr,
            comp=self.comp,
            verbose=0,
        )
        inv.fit(epochs=epochs, verbose=False)

        n_st = inv.n_sites
        mat = np.full((n_layers, n_st), np.nan)
        for si, res in enumerate(inv._results):
            lr_arr = res["log_rho"]
            n = min(len(lr_arr), n_layers)
            mat[:n, si] = lr_arr[:n]

        try:
            loss_df = inv.loss_curves()
        except Exception:
            loss_df = None
        try:
            res_df = inv.residuals()
        except Exception as exc:
            warns.append(f"residuals(): {exc}")
            res_df = None
        return inv, mat, loss_df, res_df

    def _run_2d(
        self,
        sites,
        n_layers,
        depth_max,
        epochs,
        warns,
    ):
        from ..ai.inversion.pinn2d import (
            PINNInverter2D,
        )

        inv = PINNInverter2D(
            sites,
            n_layers=n_layers,
            depth_max=depth_max,
            smoothness_weight=(self.smoothness_weight),
            lateral_weight=self.lateral_weight,
            epochs=epochs,
            lr=self.lr,
            verbose=0,
        )
        inv.fit()
        mat = inv.resistivity_section(as_log10=True)
        try:
            loss_df = inv.convergence_curve()
        except Exception:
            loss_df = None
        try:
            res_df = inv.residuals()
        except Exception as exc:
            warns.append(f"residuals(): {exc}")
            res_df = None
        return inv, mat, loss_df, res_df

    def _run_3d(
        self,
        sites,
        n_layers,
        depth_max,
        epochs,
        warns,
    ):
        from ..ai.inversion.pinn3d import (
            PINNInverter3D,
        )

        inv = PINNInverter3D(
            sites,
            n_layers=n_layers,
            depth_max=depth_max,
            smoothness_weight=(self.smoothness_weight),
            graph_weight=self.graph_weight,
            radius=self.radius,
            epochs=epochs,
            lr=self.lr,
            verbose=0,
        )
        inv.fit()
        try:
            mat = inv.resistivity_volume()
        except Exception:
            mat = None
        try:
            loss_df = inv.convergence_curve()
        except Exception:
            loss_df = None
        try:
            res_df = inv.residuals()
        except Exception as exc:
            warns.append(f"residuals(): {exc}")
            res_df = None
        return inv, mat, loss_df, res_df


# ── module-level helpers (also used by hybrid_agent) ─


def _rms_from_residuals(df: Any):
    r"""Compute per-station and global log10-rho RMS.

    Parameters
    ----------
    df : pandas.DataFrame or None
        Must have columns ``rho_obs``, ``rho_pred``,
        ``station``.

    Returns
    -------
    rms_per : dict {station: float}
    rms_global : float
    """
    if df is None or len(df) == 0:
        return {}, np.nan
    try:
        df = df.copy()
        df["log_obs"] = np.log10(np.clip(df["rho_obs"], 1e-6, None))
        df["log_pred"] = np.log10(np.clip(df["rho_pred"], 1e-6, None))
        df["sq_err"] = (df["log_obs"] - df["log_pred"]) ** 2
        per = df.groupby("station")["sq_err"].apply(
            lambda x: float(np.sqrt(np.nanmean(x)))
        )
        return per.to_dict(), float(per.mean())
    except Exception:
        return {}, np.nan


def _plot_pinn_section(
    mat: np.ndarray,
    station_names: list[str],
    n_layers: int,
    depths_km: np.ndarray,
    *,
    title: str = "PINN resistivity section",
) -> Any:
    r"""Plot log10-rho section as colour image.

    Parameters
    ----------
    mat : ndarray (n_layers, n_stations)
        log10-resistivity values.
    station_names : list of str
    n_layers : int
    depths_km : ndarray, shape (n_layers+1,)
        Depth edges in km.
    title : str

    Returns
    -------
    matplotlib.figure.Figure or None
    """
    import matplotlib.pyplot as plt

    from ..api.section import PYCSAMT_SECTION
    from ..api.station import (
        PYCSAMT_STATION_RENDERING,
    )

    n_st = len(station_names)
    if n_st == 0:
        return None

    section = PYCSAMT_SECTION.style_for("inversion")
    fw, fh = section.figsize_for(n_stations=n_st, n_y=n_layers)
    fig, ax = plt.subplots(figsize=(fw, fh))

    valid = mat[np.isfinite(mat)]
    vmin = float(np.nanpercentile(mat, 5)) if len(valid) else 0.0
    vmax = float(np.nanpercentile(mat, 95)) if len(valid) else 4.0
    im = ax.imshow(
        mat,
        aspect="auto",
        origin="upper",
        extent=(
            -0.5,
            n_st - 0.5,
            depths_km[-1],
            depths_km[0],
        ),
        cmap="jet_r",
        vmin=vmin,
        vmax=vmax,
        interpolation="nearest",
    )
    PYCSAMT_STATION_RENDERING.apply(
        ax,
        np.arange(n_st, dtype=float),
        station_names,
        preset="inversion",
        xlim=(-0.5, n_st - 0.5),
    )
    ax.set_ylabel("Depth (km)", fontsize=9)
    ax.tick_params(axis="y", labelsize=8)
    section.add_colorbar(
        im,
        ax,
        label=r"$\log_{10}\rho$ (Ohm·m)",
    )
    ax.set_title(title, fontsize=10, fontweight="bold")
    fig.tight_layout()
    return fig


def _plot_loss_curves(
    loss_df: Any,
    *,
    title: str = "PINN convergence",
) -> Any:
    """Plot Adam loss vs epoch.

    Parameters
    ----------
    loss_df : pandas.DataFrame
        Must have columns ``epoch`` and ``loss``.
        Optional column ``station`` for 1-D
        per-station curves.
    title : str

    Returns
    -------
    matplotlib.figure.Figure or None
    """
    import matplotlib.pyplot as plt

    if loss_df is None or len(loss_df) == 0:
        return None
    fig, ax = plt.subplots(figsize=(6, 3.5))

    if "station" in loss_df.columns:
        for st, grp in loss_df.groupby("station"):
            ax.semilogy(
                grp["epoch"],
                grp["loss"],
                lw=0.8,
                label=str(st),
            )
        if loss_df["station"].nunique() <= 10:
            ax.legend(fontsize=7, ncol=2)
    else:
        ax.semilogy(
            loss_df["epoch"],
            loss_df["loss"],
            lw=1.2,
        )

    ax.set_xlabel("Epoch", fontsize=9)
    ax.set_ylabel("Loss", fontsize=9)
    ax.set_title(title, fontsize=10, fontweight="bold")
    ax.tick_params(labelsize=8)
    fig.tight_layout()
    return fig


__all__ = ["PINNInversionAgent"]
