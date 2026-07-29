"""
pycsamt.agents.forward
======================

:class:`ForwardModelAgent` — Run 1-D, 2-D, or 3-D MT forward solvers.

Wraps :mod:`pycsamt.forward`:

1-D (``dim=1``)
    :class:`~pycsamt.forward.MT1DForward` on a
    :class:`~pycsamt.forward.LayeredModel`.

2-D (``dim=2``)
    :class:`~pycsamt.forward.MT2DForward` (finite-difference TE + TM) on a
    :class:`~pycsamt.forward.Grid2D`.  Supports halfspace, 1-D-layer, or
    embedded conductive-anomaly models.

3-D (``dim=3``)
    :class:`~pycsamt.forward.MT3DForward` (quasi-3D profile stacking) on a
    :class:`~pycsamt.forward.Grid3D`.  Supports halfspace and block-anomaly
    models.

The agent also computes data–model RMS when observed ``sites`` are provided
(1-D only), letting it act as a model-validation check before inversion.
"""

from __future__ import annotations

import time
from typing import Any

import numpy as np

from ._base import AgentResult, BaseAgent

_SYSTEM_PROMPT = """\
You are an expert in MT forward modelling and resistivity earth models.
Given a forward model result, write 3-4 sentences that:
1. Describe the model geometry (dimensionality, layers / grid, resistivity range).
2. Comment on the synthetic ρa and phase response (frequency range, lateral variation for 2D/3D).
3. If observed data are provided, interpret the data-model misfit (1-D only).
4. Suggest which model parameters to adjust to better fit the data or geology.
Reply in plain English. No bullet points or markdown.
"""

_DEFAULT_RESISTIVITIES = [100.0, 10.0, 1_000.0, 100.0]
_DEFAULT_THICKNESSES = [500.0, 1_000.0, 2_000.0]


class ForwardModelAgent(BaseAgent):
    """Run a 1-D, 2-D, or 3-D MT forward model.

    Parameters
    ----------
    api_key, model, llm_provider : str
    dim : {1, 2, 3}
        Forward solver dimensionality.
    freqs : array-like or None
        Frequencies (Hz).  Defaults to 40 log-spaced points 10⁻⁴–10³ Hz.

    Input keys
    ----------
    ``model`` : dict or LayeredModel or None
        **1-D / 2-D from 1-D layers:**
        ``{"resistivities": [...], "thicknesses": [...]}``.

        **2-D grid type override:**
        add ``"type": "halfspace" | "anomaly"`` and grid parameters such as
        ``"bg_rho"``, ``"anomaly_rho"``, ``"anomaly_bounds"``.

        **3-D grid type:**
        ``"type": "halfspace" | "block_anomaly"`` with grid parameters.

    ``dim`` : int, optional — overrides constructor dim for this call
    ``nx``, ``nz``, ``x_max``, ``z_max`` : int / float, optional (2-D grid)
    ``ny``, ``y_max``, ``nx_stations``, ``ny_stations`` : int / float (3-D)
    ``n_stations`` : int, optional — number of surface receivers (2-D)
    ``method`` : str, optional — ``"quasi3d"`` (default) for 3-D solver
    ``sites`` / ``path`` : Sites or str, optional — observed data for RMS (1-D)
    ``freqs`` : array-like, optional — overrides constructor default
    ``output_dir`` : str, optional
    ``component`` : ``"xy"`` (default) or ``"yx"`` (1-D component selection)

    Output data keys
    ----------------
    ``dim``              int
    ``layered_model``    LayeredModel (1-D / 2-D from 1-D)
    ``grid``             Grid2D or Grid3D (2-D / 3-D)
    ``response``         ForwardResponse / ForwardResponse2D / ForwardResponse3D
    ``rho_a``            ndarray — 1-D ρa
    ``phase``            ndarray — 1-D phase (°)
    ``rho_a_te``         ndarray (n_freqs, n_stations) — 2-D TE
    ``phase_te``         ndarray — 2-D TE phase
    ``rho_a_tm``         ndarray — 2-D TM
    ``phase_tm``         ndarray — 2-D TM phase
    ``rho_a_xy``         ndarray (n_freqs, n_stations) — 3-D XY
    ``phase_xy``         ndarray — 3-D XY phase
    ``rho_a_yx``         ndarray — 3-D YX
    ``phase_yx``         ndarray — 3-D YX phase
    ``freqs``            ndarray
    ``rms``              float or None
    ``figures``          dict
    ``figure_paths``     dict
    """

    SYSTEM_PROMPT = _SYSTEM_PROMPT

    def __init__(
        self,
        *,
        api_key: str | None = None,
        model: str | None = None,
        llm_provider: str = "claude",
        dim: int = 1,
        freqs: Any = None,
    ) -> None:
        super().__init__(
            "ForwardModelAgent",
            api_key=api_key,
            model=model,
            llm_provider=llm_provider,
            section_preset="pseudosection",
        )
        self.dim = int(dim)
        self._freqs_cfg = freqs

    # ── public ────────────────────────────────────────────────────────────────

    def execute(self, input_data: dict[str, Any]) -> AgentResult:
        self._last_cost = 0.0
        t0 = time.time()
        warnings: list[str] = []

        output_dir = input_data.get("output_dir")
        dim = int(input_data.get("dim", self.dim))

        freqs_raw = input_data.get("freqs") or self._freqs_cfg
        freqs = (
            np.asarray(freqs_raw, dtype=float)
            if freqs_raw is not None
            else np.logspace(-4, 3, 40)
        )

        try:
            from ..forward import (
                Grid2D,
                Grid3D,
                LayeredModel,
                MT1DForward,
                MT2DForward,
                MT3DForward,
            )
        except ImportError as exc:
            return AgentResult.failed(
                f"pycsamt.forward is not available: {exc}",
                hint="Ensure numpy and scipy are installed.",
                elapsed=time.time() - t0,
            )

        if dim == 1:
            return self._execute_1d(
                input_data,
                freqs,
                t0,
                warnings,
                MT1DForward,
                LayeredModel,
                output_dir,
            )
        if dim == 2:
            return self._execute_2d(
                input_data,
                freqs,
                t0,
                warnings,
                MT2DForward,
                LayeredModel,
                Grid2D,
                output_dir,
            )
        if dim == 3:
            return self._execute_3d(
                input_data,
                freqs,
                t0,
                warnings,
                MT3DForward,
                Grid3D,
                output_dir,
            )

        return AgentResult.failed(
            f"dim={dim} not supported. Use 1, 2, or 3.",
            elapsed=time.time() - t0,
        )

    # ── 1-D ──────────────────────────────────────────────────────────────────

    def _execute_1d(
        self, inp, freqs, t0, warnings, MT1DForward, LayeredModel, output_dir
    ):
        component = str(inp.get("component", "xy")).lower()
        ri, ci = (0, 1) if component == "xy" else (1, 0)

        layered = _build_layered_model(inp.get("model"), LayeredModel, warnings)
        if isinstance(layered, AgentResult):
            return layered

        try:
            response = MT1DForward(freqs=freqs).run(layered)
        except Exception as exc:
            return AgentResult.failed(
                f"1-D forward failed: {exc}", elapsed=time.time() - t0
            )

        rho_a = phase = None
        try:
            ra = np.asarray(response.rho_a)
            rho_a = ra[:, ri, ci] if ra.ndim == 3 else ra
            ph = np.asarray(response.phase)
            phase = ph[:, ri, ci] if ph.ndim == 3 else ph
        except Exception as exc:
            warnings.append(f"Could not extract ρa/phase: {exc}")

        rms: float | None = None
        sites_raw = inp.get("sites") or inp.get("path")
        if sites_raw is not None:
            try:
                rms = _compute_rms_1d(sites_raw, response, freqs, ri, ci)
            except Exception as exc:
                warnings.append(f"RMS computation failed: {exc}")

        figures, fig_paths = self._plot_1d(response, layered, output_dir, warnings)

        interp = self._llm_1d(layered, freqs, rms) if self.api_key else None

        n_layers = len(
            np.asarray(
                getattr(
                    layered,
                    "resistivity",
                    getattr(layered, "resistivities", []),
                )
            )
        )
        rms_str = f"RMS {rms:.3f}" if rms is not None else "no observed data"
        return AgentResult(
            status="success",
            summary=f"1-D forward: {n_layers} layers. {rms_str}. "
            f"{len(figures)} figure(s).",
            data={
                "dim": 1,
                "layered_model": layered,
                "response": response,
                "rho_a": rho_a,
                "phase": phase,
                "freqs": freqs,
                "rms": rms,
                "figures": figures,
                "figure_paths": fig_paths,
            },
            warnings=warnings,
            llm_interpretation=interp,
            elapsed_seconds=time.time() - t0,
            cost_estimate_usd=self._last_cost,
        )

    # ── 2-D ──────────────────────────────────────────────────────────────────

    def _execute_2d(
        self,
        inp,
        freqs,
        t0,
        warnings,
        MT2DForward,
        LayeredModel,
        Grid2D,
        output_dir,
    ):
        grid, layered = _build_grid_2d(inp, LayeredModel, Grid2D, warnings)
        if isinstance(grid, AgentResult):
            return grid

        try:
            response = MT2DForward(freqs, grid, verbose=False).run()
        except Exception as exc:
            return AgentResult.failed(
                f"2-D forward failed: {exc}", elapsed=time.time() - t0
            )

        figures, fig_paths = self._plot_2d(grid, response, output_dir, warnings)
        interp = self._llm_2d(grid, response, freqs) if self.api_key else None

        ns = response.rho_a_te.shape[1]
        nf = len(freqs)
        return AgentResult(
            status="success",
            summary=f"2-D forward: {grid.nx}×{grid.nz} cells, "
            f"{ns} stations, {nf} frequencies. "
            f"{len(figures)} figure(s).",
            data={
                "dim": 2,
                "layered_model": layered,
                "grid": grid,
                "response": response,
                "rho_a_te": response.rho_a_te,
                "phase_te": response.phase_te,
                "rho_a_tm": response.rho_a_tm,
                "phase_tm": response.phase_tm,
                "zxy": response.zxy,
                "zyx": response.zyx,
                "stations_x": response.stations_x,
                "freqs": freqs,
                "rms": None,
                "figures": figures,
                "figure_paths": fig_paths,
            },
            warnings=warnings,
            llm_interpretation=interp,
            elapsed_seconds=time.time() - t0,
            cost_estimate_usd=self._last_cost,
        )

    # ── 3-D ──────────────────────────────────────────────────────────────────

    def _execute_3d(self, inp, freqs, t0, warnings, MT3DForward, Grid3D, output_dir):
        method = str(inp.get("method", "quasi3d")).lower()
        grid = _build_grid_3d(inp, Grid3D, warnings)
        if isinstance(grid, AgentResult):
            return grid

        try:
            response = MT3DForward(freqs, grid, method=method, verbose=False).run()
        except Exception as exc:
            return AgentResult.failed(
                f"3-D forward failed: {exc}", elapsed=time.time() - t0
            )

        figures, fig_paths = self._plot_3d(grid, response, output_dir, warnings)
        interp = self._llm_3d(grid, response, freqs, method) if self.api_key else None

        ns = response.n_stations
        nf = len(freqs)
        return AgentResult(
            status="success",
            summary=f"3-D forward ({method}): "
            f"{grid.nx}×{grid.ny}×{grid.nz} cells, "
            f"{ns} stations, {nf} frequencies. "
            f"{len(figures)} figure(s).",
            data={
                "dim": 3,
                "grid": grid,
                "response": response,
                "rho_a_xy": response.rho_a_xy,
                "phase_xy": response.phase_xy,
                "rho_a_yx": response.rho_a_yx,
                "phase_yx": response.phase_yx,
                "zxy": response.zxy,
                "zyx": response.zyx,
                "stations_xy": response.stations_xy,
                "freqs": freqs,
                "rms": None,
                "figures": figures,
                "figure_paths": fig_paths,
            },
            warnings=warnings,
            llm_interpretation=interp,
            elapsed_seconds=time.time() - t0,
            cost_estimate_usd=self._last_cost,
        )

    # ── figure helpers ────────────────────────────────────────────────────────

    def _plot_1d(self, response, layered, output_dir, warnings):
        from ..forward import (
            plot_model_1d,
            plot_response_1d,
            plot_response_and_model_1d,
        )

        figures: dict[str, Any] = {}
        fig_paths: dict[str, str] = {}

        try:
            fig = plot_response_and_model_1d(response, layered)
            f = _unwrap_fig(fig)
            if f is not None:
                figures["response_and_model"] = f
                p = self._save_figure(
                    f,
                    output_dir,
                    "fwd1d_response_model",
                    warnings_list=warnings,
                )
                if p:
                    fig_paths["response_and_model"] = p
        except Exception as exc:
            warnings.append(f"plot_response_and_model_1d: {exc}")

        if "response_and_model" not in figures:
            for fn, key, tag in [
                (plot_response_1d, "response", "fwd1d_response"),
                (plot_model_1d, "model", "fwd1d_model"),
            ]:
                try:
                    f = _unwrap_fig(fn(response if key == "response" else layered))
                    if f is not None:
                        figures[key] = f
                        p = self._save_figure(
                            f, output_dir, tag, warnings_list=warnings
                        )
                        if p:
                            fig_paths[key] = p
                except Exception as exc:
                    warnings.append(f"{fn.__name__}: {exc}")

        return figures, fig_paths

    def _plot_2d(self, grid, response, output_dir, warnings):

        from ..forward import (
            plot_model_2d,
            plot_pseudosection_2d,
            plot_response_profiles,
        )

        figures: dict[str, Any] = {}
        fig_paths: dict[str, str] = {}

        for fn, key, tag, kwargs in [
            (
                plot_model_2d,
                "model_2d",
                "fwd2d_model",
                {"show_stations": True},
            ),
            (
                plot_pseudosection_2d,
                "pseudo_te",
                "fwd2d_pseudo_te",
                {"mode": "te", "quantity": "rho_a"},
            ),
            (
                plot_pseudosection_2d,
                "pseudo_tm",
                "fwd2d_pseudo_tm",
                {"mode": "tm", "quantity": "rho_a"},
            ),
            (plot_response_profiles, "profiles", "fwd2d_profiles", {}),
        ]:
            try:
                arg = grid if key == "model_2d" else response
                result = fn(arg, **kwargs)
                # plot_model_2d returns an Axes; others may return Figure or Axes
                f = _unwrap_fig(result)
                if f is None and hasattr(result, "get_figure"):
                    f = result.get_figure()
                if f is not None:
                    figures[key] = f
                    p = self._save_figure(f, output_dir, tag, warnings_list=warnings)
                    if p:
                        fig_paths[key] = p
            except Exception as exc:
                warnings.append(f"{fn.__name__}({key}): {exc}")

        return figures, fig_paths

    def _plot_3d(self, grid, response, output_dir, warnings):

        from ..forward import (
            plot_model_3d,
            plot_response_map_3d,
            plot_response_section_3d,
        )

        figures: dict[str, Any] = {}
        fig_paths: dict[str, str] = {}

        for fn, key, tag, kwargs in [
            (plot_model_3d, "model_3d", "fwd3d_model", {}),
            (plot_response_map_3d, "map_3d", "fwd3d_map", {"freq_idx": 0}),
            (plot_response_section_3d, "section_3d", "fwd3d_section", {}),
        ]:
            try:
                arg = grid if key == "model_3d" else response
                result = fn(arg, **kwargs)
                # plot_model_3d returns ndarray of Axes; get figure from first
                if hasattr(result, "__iter__") and not hasattr(result, "savefig"):
                    axes = list(result)
                    f = (
                        axes[0].get_figure()
                        if axes and hasattr(axes[0], "get_figure")
                        else None
                    )
                else:
                    f = _unwrap_fig(result)
                if f is not None:
                    figures[key] = f
                    p = self._save_figure(f, output_dir, tag, warnings_list=warnings)
                    if p:
                        fig_paths[key] = p
            except Exception as exc:
                warnings.append(f"{fn.__name__}({key}): {exc}")

        return figures, fig_paths

    # ── LLM helpers ───────────────────────────────────────────────────────────

    def _llm_1d(self, layered, freqs, rms):
        rhos = list(
            getattr(
                layered,
                "resistivity",
                getattr(layered, "resistivities", _DEFAULT_RESISTIVITIES),
            )
        )
        ths = list(
            getattr(
                layered,
                "thickness",
                getattr(layered, "thicknesses", _DEFAULT_THICKNESSES),
            )
        )
        rms_line = (
            f"  RMS vs observed: {rms:.3f}\n"
            if rms is not None
            else "  RMS: not computed (no observed data)\n"
        )
        prompt = (
            f"1-D MT forward model:\n"
            f"  Resistivities (Ω·m): {rhos}\n"
            f"  Thicknesses (m):     {ths}\n"
            f"  Frequency range: {freqs[0]:.2e}–{freqs[-1]:.2e} Hz\n"
            + rms_line
            + "Describe the model and interpret the misfit if available."
        )
        return self.query_llm(prompt, max_tokens=200)

    def _llm_2d(self, grid, response, freqs):
        rho_range = (
            float(grid.resistivity.min()),
            float(grid.resistivity.max()),
        )
        rho_te_mean = float(np.nanmean(response.rho_a_te))
        prompt = (
            f"2-D MT forward model (FD, TE+TM modes):\n"
            f"  Grid: {grid.nx}×{grid.nz} cells  "
            f"  ({grid.x_nodes[-1]:.0f} m × {grid.z_nodes[-1]:.0f} m)\n"
            f"  Resistivity range: {rho_range[0]:.1f}–{rho_range[1]:.1f} Ω·m\n"
            f"  Stations: {grid.n_stations}  "
            f"  Frequencies: {len(freqs)} ({freqs[0]:.2e}–{freqs[-1]:.2e} Hz)\n"
            f"  Mean ρa (TE): {rho_te_mean:.1f} Ω·m\n\n"
            "Describe the 2-D model and comment on the synthetic response."
        )
        return self.query_llm(prompt, max_tokens=200)

    def _llm_3d(self, grid, response, freqs, method):
        rho_range = (
            float(grid.resistivity.min()),
            float(grid.resistivity.max()),
        )
        rho_xy_mean = float(np.nanmean(response.rho_a_xy))
        prompt = (
            f"3-D MT forward model ({method}):\n"
            f"  Grid: {grid.nx}×{grid.ny}×{grid.nz} cells  "
            f"  ({grid.x_nodes[-1]:.0f}×{grid.y_nodes[-1]:.0f}×{grid.z_nodes[-1]:.0f} m)\n"
            f"  Resistivity range: {rho_range[0]:.1f}–{rho_range[1]:.1f} Ω·m\n"
            f"  Stations: {response.n_stations}  "
            f"  Frequencies: {len(freqs)} ({freqs[0]:.2e}–{freqs[-1]:.2e} Hz)\n"
            f"  Mean ρa (XY): {rho_xy_mean:.1f} Ω·m\n\n"
            "Describe the 3-D model and its synthetic MT response."
        )
        return self.query_llm(prompt, max_tokens=200)


# ── grid builders ─────────────────────────────────────────────────────────────


def _build_layered_model(model_raw, LayeredModel, warnings):
    """Return a LayeredModel or an AgentResult on failure."""
    if model_raw is None:
        warnings.append(
            "No model provided; using default 4-layer model. "
            "Pass model={'resistivities': [...], 'thicknesses': [...]}."
        )
        model_raw = {
            "resistivities": _DEFAULT_RESISTIVITIES,
            "thicknesses": _DEFAULT_THICKNESSES,
        }
    if not isinstance(model_raw, dict):
        return model_raw  # already a LayeredModel

    rhos = model_raw.get("resistivity", model_raw.get("resistivities", []))
    ths = model_raw.get("thickness", model_raw.get("thicknesses", []))
    try:
        return LayeredModel(
            resistivity=np.asarray(rhos, dtype=float),
            thickness=np.asarray(ths, dtype=float),
        )
    except Exception as exc:
        return AgentResult.failed(
            f"Could not build LayeredModel: {exc}",
            hint="Pass {'resistivities': [...], 'thicknesses': [...]}.",
            elapsed=0.0,
        )


def _build_grid_2d(inp, LayeredModel, Grid2D, warnings):
    """
    Build (Grid2D, layered_model_or_None) from input_data.

    Priority:
    1. ``inp["grid"]`` — pre-built Grid2D
    2. ``model["type"] == "halfspace"`` → Grid2D.halfspace
    3. ``model["type"] == "anomaly"``   → Grid2D.with_anomaly
    4. model has resistivities/thicknesses → Grid2D.from_1d_layers
    5. default halfspace ρ=100 Ω·m
    """
    if "grid" in inp and inp["grid"] is not None:
        return inp["grid"], None

    m = inp.get("model") or {}
    nx = int(inp.get("nx", m.get("nx", 40)))
    nz = int(inp.get("nz", m.get("nz", 25)))
    x_max = float(inp.get("x_max", m.get("x_max", 6_000.0)))
    z_max = float(inp.get("z_max", m.get("z_max", 3_000.0)))
    n_sta = int(inp.get("n_stations", m.get("n_stations", 10)))
    bg_rho = float(
        m.get(
            "bg_rho",
            m.get(
                "resistivity",
                m.get("resistivities", [100.0])[0]
                if isinstance(m.get("resistivities"), (list, np.ndarray))
                else 100.0,
            ),
        )
    )
    mtype = str(m.get("type", m.get("model_type", "auto"))).lower()

    layered = None

    # check if resistivities / thicknesses are supplied
    has_layers = bool(
        m.get("resistivities")
        or m.get("resistivity")
        or (isinstance(m, dict) and "thicknesses" in m)
    )

    try:
        if mtype == "halfspace":
            grid = Grid2D.halfspace(
                rho=bg_rho,
                nx=nx,
                nz=nz,
                x_max=x_max,
                z_max=z_max,
                n_stations=n_sta,
            )

        elif mtype == "anomaly":
            anomaly_rho = float(m.get("anomaly_rho", bg_rho / 50.0))
            anomaly_bounds = m.get(
                "anomaly_bounds",
                [x_max * 0.25, x_max * 0.75, z_max * 0.10, z_max * 0.30],
            )
            grid = Grid2D.with_anomaly(
                bg_rho=bg_rho,
                anomaly_rho=anomaly_rho,
                anomaly_bounds=anomaly_bounds,
                nx=nx,
                nz=nz,
                x_max=x_max,
                z_max=z_max,
                n_stations=n_sta,
            )

        elif has_layers or mtype in ("from_1d", "from_layers", "auto"):
            layered = _build_layered_model(
                m if has_layers else None, LayeredModel, warnings
            )
            if isinstance(layered, AgentResult):
                return layered, None
            grid = Grid2D.from_1d_layers(
                layered,
                nx=nx,
                x_max=x_max,
                n_stations=n_sta,
            )

        else:
            # fallback: 100 Ω·m halfspace
            warnings.append(f"Unknown model type {mtype!r}; using 100 Ω·m halfspace.")
            grid = Grid2D.halfspace(
                rho=100.0,
                nx=nx,
                nz=nz,
                x_max=x_max,
                z_max=z_max,
                n_stations=n_sta,
            )

    except Exception as exc:
        return (
            AgentResult.failed(
                f"Grid2D construction failed: {exc}",
                hint="Check nx, nz, x_max, z_max, n_stations in input_data.",
                elapsed=0.0,
            ),
            None,
        )

    return grid, layered


def _build_grid_3d(inp, Grid3D, warnings):
    """Build Grid3D from input_data."""
    if "grid" in inp and inp["grid"] is not None:
        return inp["grid"]

    m = inp.get("model") or {}
    nx = int(inp.get("nx", m.get("nx", 20)))
    ny = int(inp.get("ny", m.get("ny", 20)))
    nz = int(inp.get("nz", m.get("nz", 15)))
    x_max = float(inp.get("x_max", m.get("x_max", 6_000.0)))
    y_max = float(inp.get("y_max", m.get("y_max", 6_000.0)))
    z_max = float(inp.get("z_max", m.get("z_max", 4_000.0)))
    nx_sta = int(inp.get("nx_stations", m.get("nx_stations", 4)))
    ny_sta = int(inp.get("ny_stations", m.get("ny_stations", 4)))
    bg_rho = float(m.get("bg_rho", 100.0))
    mtype = str(m.get("type", m.get("model_type", "halfspace"))).lower()

    try:
        if mtype in ("halfspace", "auto"):
            grid = Grid3D.halfspace(
                rho=bg_rho,
                nx=nx,
                ny=ny,
                nz=nz,
                x_max=x_max,
                y_max=y_max,
                z_max=z_max,
                nx_stations=nx_sta,
                ny_stations=ny_sta,
            )

        elif mtype in ("block_anomaly", "anomaly", "block"):
            anomaly_rho = float(m.get("anomaly_rho", bg_rho / 50.0))
            grid = Grid3D.block_anomaly(
                bg_rho=bg_rho,
                anomaly_rho=anomaly_rho,
                nx=nx,
                ny=ny,
                nz=nz,
                x_max=x_max,
                y_max=y_max,
                z_max=z_max,
                nx_stations=nx_sta,
                ny_stations=ny_sta,
            )

        else:
            warnings.append(f"Unknown 3-D model type {mtype!r}; using halfspace.")
            grid = Grid3D.halfspace(
                rho=100.0,
                nx=nx,
                ny=ny,
                nz=nz,
                x_max=x_max,
                y_max=y_max,
                z_max=z_max,
                nx_stations=nx_sta,
                ny_stations=ny_sta,
            )

    except Exception as exc:
        return AgentResult.failed(
            f"Grid3D construction failed: {exc}",
            hint="Check nx, ny, nz, x_max, y_max, z_max in input_data.",
            elapsed=0.0,
        )

    return grid


# ── misc helpers ──────────────────────────────────────────────────────────────


def _unwrap_fig(obj):
    """Return a matplotlib Figure from an Axes, Figure, or ndarray of Axes."""
    if obj is None:
        return None
    if hasattr(obj, "savefig"):
        return obj
    if hasattr(obj, "get_figure"):
        return obj.get_figure()
    return None


def _compute_rms_1d(sites_raw, response, freqs_fwd, ri, ci):
    """Compute mean-log-ρa RMS between 1-D forward and observed sites."""
    from ..emtools._core import (
        _get_z_block,
        _iter_items,
        ensure_sites,
    )

    sites = ensure_sites(sites_raw, verbose=0)
    residuals: list[float] = []

    ra = np.asarray(response.rho_a)
    rho_fwd = ra[:, ri, ci] if ra.ndim == 3 else ra

    for _, ed in enumerate(_iter_items(sites)):
        _, z, fr = _get_z_block(ed)
        if z is None or fr is None:
            continue
        rho_raw = getattr(ed, "rho", None)
        rho_obs = (
            rho_raw[:, ri, ci]
            if rho_raw is not None
            else (0.2 / np.where(fr == 0, np.nan, fr)) * np.abs(z[:, ri, ci]) ** 2
        )
        per_obs = 1.0 / np.where(fr == 0, np.nan, fr)
        per_fwd = 1.0 / np.where(freqs_fwd == 0, np.nan, freqs_fwd)
        mask = np.isfinite(per_obs) & (rho_obs > 0)
        if not mask.any():
            continue
        interp = np.interp(
            np.log10(per_obs[mask]),
            np.log10(per_fwd[np.isfinite(per_fwd)]),
            np.log10(np.clip(rho_fwd[np.isfinite(per_fwd)], 1e-6, None)),
        )
        residuals.extend(
            (np.log10(np.clip(rho_obs[mask], 1e-6, None)) - interp).tolist()
        )

    if not residuals:
        raise ValueError("No valid observed data to compare.")
    return float(np.sqrt(np.mean(np.array(residuals) ** 2)))


__all__ = ["ForwardModelAgent"]
