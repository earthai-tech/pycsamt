"""
pycsamt.agents.inversion_backend
==================================

:class:`InversionBackendAgent` — Drive the ``pycsamt.inversion`` physics-based
inversion backends from within the agent system.

Supports all backends registered in :mod:`pycsamt.inversion`:

* **builtin** — pure-Python FD MT/TDEM inversion (no extra dependencies)
* **simpeg** — SimPEG MT/3D natural-source inversion (optional)
* **pygimli** — pyGIMLi 1-D MT/TDEM inversion (optional)
* **occam2d** — Occam2D data-file runner (requires Occam binary)
* **modem** — ModEM3D data-file runner (requires ModEM binary)

The agent wraps :func:`~pycsamt.inversion.run_inversion` and
:class:`~pycsamt.inversion.InversionConfig`, converts the raw
:class:`~pycsamt.inversion.InversionResult` into a standard
:class:`~pycsamt.agents._base.AgentResult`, and plots the final resistivity
section using the PYCSAMT API.
"""

from __future__ import annotations

import time
from typing import Any

import numpy as np

from ._base import AgentResult, BaseAgent

_SYSTEM_PROMPT = """\
You are an expert in MT inversion and subsurface resistivity modelling.
Given an inversion result, write 4-5 sentences that:
1. State the backend used, dimensionality, and convergence (RMS, n_iter).
2. Describe the recovered resistivity model (range, dominant structures).
3. Assess the fit quality and whether the RMS target was reached.
4. Identify stations or depth ranges with elevated misfit.
5. Recommend regularisation adjustments or mesh refinements for the next run.
Reply in plain scientific English.
"""


class InversionBackendAgent(BaseAgent):
    """Drive pycsamt.inversion physics-based backends.

    Parameters
    ----------
    api_key, model, llm_provider : str
    backend : str
        Inversion backend: ``'builtin'`` (default), ``'simpeg'``,
        ``'pygimli'``, ``'occam2d'``, ``'modem'``.
    dimension : str
        ``'1d'`` (default), ``'2d'``, or ``'3d'``.
    method : str
        ``'mt'`` (default), ``'amt'``, ``'csamt'``, or ``'tdem'``.
    n_layers : int
        Number of depth layers for 1-D inversion (default 5).
    max_iter : int
        Maximum inversion iterations (default 80).
    regularization : str
        ``'smooth'`` (default), ``'damped'``, or ``'blocky'``.
    error_floor : float
        Relative data error floor (default 0.05 = 5 %).

    Input keys
    ----------
    ``sites`` / ``path`` : Sites or str
    ``backend``, ``dimension``, ``method`` : str, optional overrides
    ``n_layers``, ``max_iter`` : int, optional overrides
    ``regularization``, ``error_floor`` : optional overrides
    ``backend_options`` : dict, optional — forwarded to InversionConfig
    ``output_dir`` : str, optional

    Output data keys
    ----------------
    ``inversion_result``  InversionResult
    ``rms``               float
    ``n_iter``            int
    ``log_rho_section``   ndarray (n_layers, n_stations)
    ``station_names``     list[str] or None
    ``backend``           str
    ``dimension``         str
    ``figures``           dict
    ``figure_paths``      dict
    """

    SYSTEM_PROMPT = _SYSTEM_PROMPT

    def __init__(
        self,
        *,
        api_key: str | None = None,
        model: str | None = None,
        llm_provider: str = "claude",
        backend: str = "builtin",
        dimension: str = "1d",
        method: str = "mt",
        n_layers: int = 5,
        max_iter: int = 80,
        regularization: str = "smooth",
        error_floor: float = 0.05,
    ) -> None:
        super().__init__(
            "InversionBackendAgent",
            api_key=api_key,
            model=model,
            llm_provider=llm_provider,
            section_preset="inversion",
        )
        self.backend = backend
        self.dimension = dimension
        self.method = method
        self.n_layers = n_layers
        self.max_iter = max_iter
        self.regularization = regularization
        self.error_floor = error_floor

    def execute(self, input_data: dict[str, Any]) -> AgentResult:
        self._last_cost = 0.0
        t0 = time.time()
        warnings: list[str] = []

        try:
            from ..inversion import (
                InversionConfig,
                available_backends,
                run_inversion,
            )
        except ImportError as exc:
            return AgentResult.failed(
                f"pycsamt.inversion not available: {exc}",
                hint="Ensure pycsamt.inversion is installed.",
                elapsed=time.time() - t0,
            )

        from ..emtools._core import ensure_sites

        sites_raw = input_data.get("sites") or input_data.get("path")
        if sites_raw is None:
            return AgentResult.failed(
                "No 'sites' or 'path'.", elapsed=time.time() - t0
            )
        try:
            sites = ensure_sites(sites_raw, verbose=0)
        except Exception as exc:
            return AgentResult.failed(str(exc), elapsed=time.time() - t0)

        backend = str(input_data.get("backend", self.backend))
        dimension = str(input_data.get("dimension", self.dimension))
        method = str(input_data.get("method", self.method))
        n_layers = int(input_data.get("n_layers", self.n_layers))
        max_iter = int(input_data.get("max_iter", self.max_iter))
        regularization = str(
            input_data.get("regularization", self.regularization)
        )
        error_floor = float(input_data.get("error_floor", self.error_floor))
        backend_options = dict(input_data.get("backend_options", {}))
        output_dir = input_data.get("output_dir", "pycsamt_inversion_output")

        # ── check backend availability ────────────────────────────────────
        avail = available_backends()
        if backend not in avail:
            warnings.append(
                f"Backend {backend!r} not in registry. "
                f"Available: {list(avail)}. Falling back to 'builtin'."
            )
            backend = "builtin"

        if avail.get(backend) is False:
            warnings.append(
                f"Backend {backend!r} is registered but its dependencies are "
                "not installed. Falling back to 'builtin'."
            )
            backend = "builtin"

        # ── build InversionConfig ─────────────────────────────────────────
        try:
            cfg = InversionConfig(
                method=method,
                dimension=dimension,
                backend=backend,
                data=sites,
                workdir=str(output_dir),
                output_dir=str(output_dir),
                n_layers=n_layers,
                error_floor=error_floor,
                regularization=regularization,
                max_iter=max_iter,
                backend_options=backend_options,
            )
        except Exception as exc:
            return AgentResult.failed(
                f"InversionConfig construction failed: {exc}",
                elapsed=time.time() - t0,
            )

        # ── run inversion ─────────────────────────────────────────────────
        self._log.info(
            "Running inversion: backend=%s, dim=%s, method=%s, n_layers=%d",
            backend,
            dimension,
            method,
            n_layers,
        )
        try:
            inv_result = run_inversion(cfg)
        except Exception as exc:
            return AgentResult.failed(
                f"Inversion failed: {exc}",
                hint=f"Check backend availability: available_backends() = {avail}",
                elapsed=time.time() - t0,
            )

        # ── extract results ───────────────────────────────────────────────
        rms = float(getattr(inv_result, "rms", float("nan")))
        n_iter = int(getattr(inv_result, "n_iter", 0))

        # log_rho_section: (n_layers, n_stations)
        log_rho_section = None
        station_names = None
        try:
            model = getattr(inv_result, "model", None)
            if model is not None:
                sec = (
                    model.get("log_rho_section")
                    if hasattr(model, "get")
                    else None
                )
                if sec is not None:
                    log_rho_section = np.asarray(sec, dtype=float)
                    if log_rho_section.ndim == 1:
                        log_rho_section = log_rho_section[:, np.newaxis]
        except Exception as exc:
            warnings.append(f"Could not extract log_rho_section: {exc}")

        try:
            station_names = list(
                getattr(inv_result, "station_names", None) or []
            )
        except Exception:
            pass

        # ── figures ───────────────────────────────────────────────────────
        figures: dict[str, Any] = {}
        fig_paths: dict[str, str] = {}

        if log_rho_section is not None:
            try:
                fig = _plot_inversion_section(
                    log_rho_section,
                    station_names,
                    rms,
                    backend,
                    dimension,
                )
                if fig is not None:
                    figures["inversion_section"] = fig
                    p = self._save_figure(
                        fig,
                        str(output_dir),
                        "inversion_backend_section",
                        warnings_list=warnings,
                    )
                    if p:
                        fig_paths["inversion_section"] = p
            except Exception as exc:
                warnings.append(f"Inversion section figure: {exc}")

        # convergence curve
        try:
            history = getattr(inv_result, "history", None)
            if history is not None:
                from ..ai.plot.convergence import (
                    plot_convergence,
                )

                fig_cv = plot_convergence(history)
                if fig_cv is not None:
                    f = (
                        fig_cv
                        if hasattr(fig_cv, "savefig")
                        else (
                            fig_cv.get_figure()
                            if hasattr(fig_cv, "get_figure")
                            else None
                        )
                    )
                    if f is not None:
                        figures["convergence"] = f
                        p = self._save_figure(
                            f,
                            str(output_dir),
                            "inversion_convergence",
                            warnings_list=warnings,
                        )
                        if p:
                            fig_paths["convergence"] = p
        except Exception as exc:
            warnings.append(f"Convergence figure: {exc}")

        # ── LLM interpretation ────────────────────────────────────────────
        interp: str | None = None
        if self.api_key:
            rms_str = f"{rms:.4f}" if not np.isnan(rms) else "N/A"
            n_sta = (
                log_rho_section.shape[1]
                if log_rho_section is not None
                else "?"
            )
            prompt = (
                f"Inversion result:\n"
                f"  Backend: {backend} | Dimension: {dimension} | Method: {method}\n"
                f"  Stations: {n_sta} | Layers: {n_layers}\n"
                f"  RMS: {rms_str} | Iterations: {n_iter}\n"
                f"  Regularization: {regularization} | Error floor: {error_floor:.0%}\n"
                f"  Warnings: {warnings[:3] if warnings else 'none'}\n\n"
                "Evaluate the inversion convergence and interpret the result."
            )
            interp = self.query_llm(prompt, max_tokens=230)

        elapsed = time.time() - t0
        rms_str = f"{rms:.4f}" if not np.isnan(rms) else "N/A"
        return AgentResult(
            status="success" if not np.isnan(rms) else "needs_review",
            summary=(
                f"Inversion ({backend}, {dimension}, {method}): "
                f"RMS={rms_str}, {n_iter} iterations. "
                f"{len(figures)} figures."
            ),
            data={
                "inversion_result": inv_result,
                "rms": rms,
                "n_iter": n_iter,
                "log_rho_section": log_rho_section,
                "station_names": station_names,
                "backend": backend,
                "dimension": dimension,
                "method": method,
                "figures": figures,
                "figure_paths": fig_paths,
            },
            warnings=warnings,
            llm_interpretation=interp,
            elapsed_seconds=elapsed,
            cost_estimate_usd=self._last_cost,
        )


# ── private helpers ───────────────────────────────────────────────────────────


def _plot_inversion_section(
    log_rho: np.ndarray,
    station_names: list[str] | None,
    rms: float,
    backend: str,
    dimension: str,
) -> Any:
    """Plot recovered log₁₀ρ section using PYCSAMT_SECTION API."""
    import matplotlib.pyplot as plt

    from ..api.section import PYCSAMT_SECTION
    from ..api.station import PYCSAMT_STATION_RENDERING

    n_layers, n_sta = log_rho.shape
    if n_sta == 0:
        return None

    names = station_names or [f"S{i + 1}" for i in range(n_sta)]
    section = PYCSAMT_SECTION.style_for("inversion")
    fig_w, fig_h = section.figsize_for(n_stations=n_sta, n_y=n_layers)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    vv = log_rho[np.isfinite(log_rho)]
    vmin = float(np.percentile(vv, 5)) if vv.size else 0.0
    vmax = float(np.percentile(vv, 95)) if vv.size else 4.0

    im = ax.imshow(
        log_rho,
        aspect="auto",
        origin="upper",
        extent=(-0.5, n_sta - 0.5, n_layers, 0),
        cmap="jet_r",
        vmin=vmin,
        vmax=vmax,
        interpolation="bilinear",
    )
    PYCSAMT_STATION_RENDERING.apply(
        ax,
        np.arange(n_sta, dtype=float),
        names,
        preset="inversion",
        xlim=(-0.5, n_sta - 0.5),
    )
    ax.set_ylabel("Layer index", fontsize=9)
    ax.tick_params(labelsize=8)
    section.add_colorbar(im, ax, label="$\\log_{10}\\rho$ (Ω·m)")
    rms_str = f"RMS = {rms:.4f}" if not np.isnan(rms) else ""
    ax.set_title(
        f"Inversion ({backend}, {dimension})  {rms_str}",
        fontsize=10,
        fontweight="bold",
    )
    fig.tight_layout()
    return fig


__all__ = ["InversionBackendAgent"]
