# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Callbacks for the Inversion page (Traditional + AI Neural + PyTorch installer)."""

from __future__ import annotations

import importlib
import re
import subprocess
import sys
import threading
import time

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from dash import (
    Input,
    Output,
    State,
    clientside_callback,
    ctx,
    html,
    no_update,
)
from dash.exceptions import PreventUpdate

_INV_TRACK_IDS = ["trad", "ai", "pinn", "hybrid"]
_INV_OUT_IDS = ["result", "conv", "stats", "log", "data-fit"]

_INV_TRACK_ICON = {
    "trad": "bi-cpu",
    "ai": "bi-robot",
    "pinn": "bi-bezier2",
    "hybrid": "bi-node-plus",
}
_INV_TRACK_LABEL = {
    "trad": "Traditional",
    "ai": "AI Neural",
    "pinn": "PINN",
    "hybrid": "Hybrid",
}

from pycsamt.app.web.cache import cache_get
from pycsamt.app.web.layout import IDs
from pycsamt.app.web.utils import (
    apply_web_dark_theme,
    apply_web_light_theme,
    empty_src,
    fig_to_src,
)

_SHOW = {"display": "block"}
_HIDE = {"display": "none"}

try:
    import torch as _torch

    _TORCH_OK = True
except ImportError:
    _TORCH_OK = False


# ---

_install_lock = threading.Lock()
_install_state: dict = {
    "status": "idle",  # idle | running | done | error
    "progress": 0,
    "log_lines": [],
    "done_time": 0.0,
}


def _parse_pip_progress(lines: list[str]) -> int:
    """Estimate 0-99% progress from accumulated pip output lines."""
    progress = 5
    for line in lines:
        lo = line.lower()
        if "looking" in lo or "collecting" in lo or "resolving" in lo:
            progress = max(progress, 8)
        elif "requirement already satisfied" in lo:
            progress = max(progress, 98)
        elif "downloading" in lo:
            m = re.search(r"(\d+)%", line)
            if m:
                pct = int(m.group(1))
                # map download pct 0-100% → overall 20-78%
                progress = max(progress, 20 + int(pct * 0.58))
            else:
                progress = max(progress, 20)
        elif "installing collected" in lo:
            progress = max(progress, 82)
        elif "running setup" in lo or "building wheel" in lo:
            progress = max(progress, 88)
        elif "successfully installed" in lo:
            progress = 99  # will be set to 100 after thread confirms
    return min(progress, 99)


def _run_pip_install() -> None:
    """Run 'pip install torch' in background thread and update _install_state."""
    global _TORCH_OK

    with _install_lock:
        _install_state.update(
            status="running",
            progress=5,
            log_lines=["Starting pip install torch …"],
            done_time=0.0,
        )

    cmd = [
        sys.executable,
        "-m",
        "pip",
        "install",
        "torch",
        "--index-url",
        "https://download.pytorch.org/whl/cpu",
        "--no-cache-dir",
    ]

    try:
        proc = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
        )

        for raw_line in proc.stdout:
            line = raw_line.rstrip()
            with _install_lock:
                _install_state["log_lines"].append(line)
                _install_state["progress"] = _parse_pip_progress(
                    _install_state["log_lines"]
                )

        proc.wait()

        if proc.returncode == 0:
            # Try to import torch in this process immediately
            importlib.invalidate_caches()
            try:
                import torch  # noqa: F401

                _TORCH_OK = True
                note = "✓ torch imported in current process — ready to use."
            except ImportError:
                note = "Restart the web app to activate PyTorch."

            with _install_lock:
                _install_state["log_lines"].append(note)
                _install_state.update(
                    status="done",
                    progress=100,
                    done_time=time.time(),
                )
        else:
            with _install_lock:
                _install_state.update(
                    status="error", progress=_install_state["progress"]
                )

    except Exception as exc:
        with _install_lock:
            _install_state["log_lines"].append(f"ERROR: {exc}")
            _install_state["status"] = "error"


# ---

_LOG_RULES: list[tuple[re.Pattern, str]] = [
    (re.compile(r"^\s*(={3,}|-{3,}|\*{3,}|#{2,})"), "log-head"),
    (re.compile(r"iter|iteration|step\s+\d", re.I), "log-iter"),
    (re.compile(r"converge|misfit\s*<|rms.*<|done|finish", re.I), "log-conv"),
    (re.compile(r"warn|caution", re.I), "log-warn"),
    (re.compile(r"error|fail|traceback|exception", re.I), "log-err"),
    (re.compile(r"load|read|write|sav|export", re.I), "log-info"),
]


def _classify(line: str) -> str:
    for pat, cls in _LOG_RULES:
        if pat.search(line):
            return cls
    return ""


def _coloured_log(lines: list[str]) -> list[html.Span]:
    return [html.Span(line + "\n", className=_classify(line) or None) for line in lines]


# ---


def _build_true_model(table_data, halfspace_rho):
    from pycsamt.forward import LayeredModel

    if not table_data:
        raise ValueError("True model needs at least one finite layer.")
    halfspace_rho = float(halfspace_rho or 1000.0)
    rho = np.array([float(r["resistivity"]) for r in table_data])
    thk = np.array([float(r["thickness"]) for r in table_data])
    if np.any(thk <= 0) or np.any(rho <= 0) or halfspace_rho <= 0:
        raise ValueError("All resistivities and thicknesses must be > 0.")
    return LayeredModel(resistivity=np.append(rho, halfspace_rho), thickness=thk)


def _make_synthetic_1d(model_true, method, freqs, noise):
    from pycsamt.inversion.data import EMData

    method = (method or "mt").lower()
    noise = float(noise or 0.05)
    if method == "tdem":
        from pycsamt.forward.em1d import TEM1DForward

        resp = TEM1DForward(freqs, loop_radius=50.0).run(model_true)
        values = resp.dBz_dt * (1 + noise * np.random.randn(len(freqs)))
        return EMData(method="tdem", times=freqs, values=values)
    elif method == "csamt":
        from pycsamt.forward.em1d import CSAMT1DForward

        resp = CSAMT1DForward(freqs, source_offset=5000.0).run(model_true)
        return EMData(
            method="csamt",
            frequencies=freqs,
            rho_a=resp.rho_a * (1 + noise * np.random.randn(len(freqs))),
            phase=resp.phase + noise * 45 * np.random.randn(len(freqs)),
        )
    else:
        from pycsamt.forward.em1d import MT1DForward

        resp = MT1DForward(freqs).run(model_true)
        return EMData(
            method="mt",
            frequencies=freqs,
            rho_a=resp.rho_a * (1 + noise * np.random.randn(len(freqs))),
            phase=resp.phase + noise * 45 * np.random.randn(len(freqs)),
        )


def _make_synthetic_2d(model_true, method, freqs, noise, n_stations, profile_len):
    from pycsamt.forward.em1d import MT1DForward
    from pycsamt.inversion.data import EMData

    n_st = int(n_stations or 8)
    noise = float(noise or 0.05)
    rho_all, phase_all = [], []
    for i in range(n_st):
        np.random.seed(i)
        resp = MT1DForward(freqs).run(model_true)
        rho_all.append(resp.rho_a * (1 + noise * np.random.randn(len(freqs))))
        phase_all.append(resp.phase + noise * 45 * np.random.randn(len(freqs)))
    return EMData(
        method="mt",
        frequencies=freqs,
        rho_a=np.array(rho_all),
        phase=np.array(phase_all),
        station_names=[f"S{i:02d}" for i in range(n_st)],
        station_x=np.linspace(0.0, float(profile_len or 5000), n_st),
    )


def _plot_model_step(ax, model, dark, title="True Model"):
    rho, thk = model.resistivity, model.thickness
    tops = np.concatenate([[0.0], np.cumsum(thk)])
    d, r = [], []
    for i in range(len(rho) - 1):
        d += [tops[i], tops[i + 1]]
        r += [rho[i], rho[i]]
    bot = tops[-1] * 1.2 if tops[-1] > 0 else 500.0
    d += [tops[-1], bot]
    r += [rho[-1], rho[-1]]
    colour = "#89b4fa" if dark else "#1e66f5"
    ax.plot(r, d, color=colour, lw=2)
    ax.fill_betweenx(d, r, alpha=0.12, color=colour)
    ax.set_xscale("log")
    ax.invert_yaxis()
    ax.set_xlabel("Resistivity (Ω·m)")
    ax.set_ylabel("Depth (m)")
    ax.set_title(title)
    ax.set_xlim(left=0.1)


def _plot_1d_result(model_true, result, dark):
    from pycsamt.inversion.plot import plot_model

    has_true = model_true is not None
    ncols = 3 if has_true else 2
    fig, axes = plt.subplots(1, ncols, figsize=(5 * ncols, 5))
    idx = 0
    if has_true:
        _plot_model_step(axes[idx], model_true, dark, "True Model")
        idx += 1
    plot_model(result, ax=axes[idx])
    axes[idx].set_title(f"Recovered  (RMS = {result.rms:.3f})")
    unc = result.uncertainty
    if unc and unc.depth_confidence is not None and hasattr(result.model, "thickness"):
        tops = np.concatenate([[0.0], np.cumsum(result.model.thickness)])
        conf = np.asarray(unc.depth_confidence)
        conf = conf / (conf.max() + 1e-9)
        for d, a in zip(tops, (1 - conf[: len(tops)]) * 0.28):
            axes[idx].axhspan(d, d + 50, alpha=float(a), color="#cba6f7", lw=0)
    idx += 1
    if hasattr(result, "data") and result.data is not None:
        data = result.data
        pred = result.predicted
        c_o = "#f38ba8" if dark else "#d20f39"
        c_p = "#a6e3a1" if dark else "#40a02b"
        if hasattr(data, "frequencies") and data.frequencies is not None:
            T = 1.0 / data.frequencies
            rho_obs = np.asarray(data.rho_a)
            if rho_obs.ndim > 1:
                rho_obs = rho_obs[0]
            axes[idx].loglog(T, rho_obs, "o", color=c_o, ms=4, label="Observed")
            if pred is not None and hasattr(pred, "rho_a"):
                rho_pred = np.asarray(pred.rho_a)
                if rho_pred.ndim > 1:
                    rho_pred = rho_pred[0]
                axes[idx].loglog(T, rho_pred, "-", color=c_p, lw=2, label="Predicted")
            axes[idx].set_xlabel("Period (s)")
            axes[idx].set_ylabel("ρₐ (Ω·m)")
            axes[idx].set_title("Data Fit")
            axes[idx].legend(fontsize=9)
    fig.tight_layout()
    return fig


def _plot_2d_result(result, dark):
    from pycsamt.inversion.plot import plot_model, plot_rms

    fig, axes = plt.subplots(
        2, 1, figsize=(12, 8), gridspec_kw={"height_ratios": [3, 1]}
    )
    plot_model(result, ax=axes[0])
    axes[0].set_title(f"Recovered 2-D Section  (RMS = {result.rms:.3f})")
    plot_rms(result, ax=axes[1])
    fig.tight_layout()
    return fig


def _plot_convergence(result, dark):
    if result.history is None:
        return None
    arrs = result.history.arrays()
    iters = arrs.get("iteration", np.arange(max(len(v) for v in arrs.values())))
    c_rms = "#f38ba8" if dark else "#d20f39"
    c_m = "#89b4fa" if dark else "#1e66f5"
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    if "rms" in arrs:
        axes[0].semilogy(iters, arrs["rms"], color=c_rms, lw=2)
        axes[0].axhline(1.0, ls="--", color="#6c7086", lw=1, label="RMS=1")
        axes[0].set_xlabel("Iteration")
        axes[0].set_ylabel("RMS misfit")
        axes[0].set_title("Data Misfit")
        axes[0].legend(fontsize=9)
    if "phi_m" in arrs:
        axes[1].plot(iters, arrs["phi_m"], color=c_m, lw=2)
        axes[1].set_xlabel("Iteration")
        axes[1].set_ylabel("φ_m")
        axes[1].set_title("Model Roughness")
    elif "objective" in arrs:
        axes[1].semilogy(iters, arrs["objective"], color=c_m, lw=2)
        axes[1].set_xlabel("Iteration")
        axes[1].set_ylabel("Objective")
        axes[1].set_title("Objective Function")
    fig.tight_layout()
    return fig


def _stats_table(result, method, dim, backend, n_layers, regularize):
    status_color = "text-success" if result.converged else "text-warning"
    rows = [
        ("Status", html.Span(result.status, className=status_color)),
        ("RMS misfit", f"{result.rms:.4f}"),
        ("Iterations", str(result.n_iter)),
        ("Method", f"{method.upper()} {dim.upper()}"),
        ("Backend", backend),
        ("Regularization", regularize),
        ("N layers", str(n_layers)),
        ("Warnings", str(len(result.warnings)) if result.warnings else "0"),
    ]
    st_rms = (result.metadata or {}).get("station_rms")
    if st_rms is not None:
        a = np.asarray(st_rms)
        rows.append(("Station RMS", f"{a.mean():.3f} ± {a.std():.3f}"))
    trs = [
        html.Tr(
            [
                html.Td(
                    k,
                    style={
                        "color": "#a6adc8",
                        "paddingRight": "16px",
                        "fontSize": "11px",
                        "whiteSpace": "nowrap",
                    },
                ),
                html.Td(v, style={"color": "#cdd6f4", "fontSize": "11px"}),
            ]
        )
        for k, v in rows
    ]
    return html.Table(trs, style={"width": "100%", "borderCollapse": "collapse"})


def _fwd2d_to_tensor(resp) -> np.ndarray:
    """Convert a ForwardResponse2D → (1, 2, n_freqs, n_stations) float32 tensor.

    Uses TE mode apparent resistivity (log10) + normalised phase.
    """
    rho_a = np.asarray(resp.rho_a_te, dtype=np.float64)  # (n_freqs, n_stations)
    phase = np.asarray(resp.phase_te, dtype=np.float64)
    ch0 = np.log10(np.maximum(rho_a, 1e-6)).astype(np.float32)
    ch1 = (phase / 90.0).astype(np.float32)
    X = np.stack([ch0, ch1], axis=0)  # (2, n_freqs, n_stations)
    return X[np.newaxis]  # (1, 2, n_freqs, n_stations)


def _generate_2d_dataset(
    n_profiles: int,
    n_stations: int,
    n_layers: int,
    freqs: np.ndarray,
    noise: float,
    seed: int,
) -> tuple:
    """Generate (X, y) for EMInverter2D in-process.

    X : (n_profiles, 2, n_freqs, n_stations)   channels: log10_rho_a, norm_phase
    y : (n_profiles, n_layers, n_stations)       log10(rho) sections
    """
    from pycsamt.forward import LayeredModel
    from pycsamt.forward.em1d import MT1DForward

    n_freqs = len(freqs)
    rng = np.random.default_rng(seed)
    fwd = MT1DForward(freqs)

    X = np.zeros((n_profiles, 2, n_freqs, n_stations), dtype=np.float32)
    y = np.zeros((n_profiles, n_layers, n_stations), dtype=np.float32)

    for i in range(n_profiles):
        # Correlated log10(rho) across stations: shared background + per-station deviation
        bg_log_rho = rng.normal(2.0, 0.5, size=n_layers)
        for j in range(n_stations):
            log_rho = bg_log_rho + rng.normal(0, 0.25, size=n_layers)
            log_rho = np.clip(log_rho, 0.5, 4.0)
            thk = 10.0 ** rng.uniform(1.8, 3.3, size=n_layers - 1)
            model = LayeredModel(resistivity=10.0**log_rho, thickness=thk)
            resp = fwd.run(model)
            rho_a = np.asarray(resp.rho_a, dtype=np.float32)
            phase = np.asarray(resp.phase, dtype=np.float32)
            rho_a = np.maximum(rho_a * (1 + noise * rng.standard_normal(n_freqs)), 1e-3)
            phase = phase + noise * 45 * rng.standard_normal(n_freqs)
            X[i, 0, :, j] = np.log10(rho_a)
            X[i, 1, :, j] = phase / 90.0  # normalise to ~[-1, 1]
            y[i, :, j] = log_rho.astype(np.float32)

    return X, y


def _plot_predicted_section_2d(
    y_pred_section: np.ndarray, fwd_resp, n_layers: int, dark: bool
):
    """2-panel plot for a single AI-predicted 2-D section from a real ForwardResponse2D.

    y_pred_section : (n_layers, n_stations) — log10(rho) section
    fwd_resp       : ForwardResponse2D with .stations_x, .freqs, .rho_a_te
    """
    cmap = "viridis_r"
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    # Left: predicted resistivity section
    im0 = axes[0].imshow(y_pred_section, aspect="auto", cmap=cmap, origin="upper")
    axes[0].set_title("AI-predicted log₁₀(ρ) section")
    axes[0].set_xlabel("Station index")
    axes[0].set_ylabel("Depth layer")
    fig.colorbar(im0, ax=axes[0], label="log₁₀(ρ)")

    # Right: input pseudosection (TE rho_a) that was inverted
    rho_te = np.log10(np.maximum(np.asarray(fwd_resp.rho_a_te), 1e-6))
    im1 = axes[1].imshow(rho_te, aspect="auto", cmap=cmap, origin="upper")
    axes[1].set_title("Input: log₁₀(ρₐ TE) pseudosection")
    axes[1].set_xlabel("Station index")
    axes[1].set_ylabel("Frequency index (high→low period)")
    fig.colorbar(im1, ax=axes[1], label="log₁₀(ρₐ)")

    fig.suptitle("AI 2-D inversion of session ForwardResponse2D", fontsize=11)
    fig.tight_layout()
    return fig


def _plot_ai_section_2d(y_true: np.ndarray, y_pred: np.ndarray, dark: bool):
    """4-panel plot for EMInverter2D test results.

    y_true / y_pred : (n_profiles, n_layers, n_stations)
    """
    cmap = "viridis_r"
    fig, axes = plt.subplots(2, 2, figsize=(13, 7))

    vmin = min(float(y_true[0].min()), float(y_pred[0].min()))
    vmax = max(float(y_true[0].max()), float(y_pred[0].max()))

    im0 = axes[0, 0].imshow(
        y_true[0],
        aspect="auto",
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        origin="upper",
    )
    axes[0, 0].set_title("True section (profile 1)")
    axes[0, 0].set_ylabel("Depth layer")
    axes[0, 0].set_xlabel("Station")
    fig.colorbar(im0, ax=axes[0, 0], label="log₁₀(ρ)")

    im1 = axes[0, 1].imshow(
        y_pred[0],
        aspect="auto",
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        origin="upper",
    )
    axes[0, 1].set_title("Predicted section (profile 1)")
    axes[0, 1].set_xlabel("Station")
    fig.colorbar(im1, ax=axes[0, 1], label="log₁₀(ρ)")

    diff = y_pred[0] - y_true[0]
    lim = max(float(np.abs(diff).max()), 0.01)
    im2 = axes[1, 0].imshow(
        diff,
        aspect="auto",
        cmap="RdBu_r",
        vmin=-lim,
        vmax=lim,
        origin="upper",
    )
    axes[1, 0].set_title("Residual (pred − true)")
    axes[1, 0].set_ylabel("Depth layer")
    axes[1, 0].set_xlabel("Station")
    fig.colorbar(im2, ax=axes[1, 0], label="Δlog₁₀(ρ)")

    mae_per_profile = np.abs(y_pred - y_true).mean(axis=(1, 2))
    c = "#89b4fa" if dark else "#1e66f5"
    axes[1, 1].plot(mae_per_profile, "o-", color=c, ms=3, lw=1.5)
    axes[1, 1].set_xlabel("Test profile index")
    axes[1, 1].set_ylabel("MAE (log₁₀ ρ)")
    axes[1, 1].set_title(f"Test MAE  (mean={mae_per_profile.mean():.3f})")

    fig.tight_layout()
    return fig


def _plot_gcn_section_3d(
    y_true: np.ndarray, y_pred: np.ndarray, n_layers: int, dark: bool
):
    """3-panel plot for GCNInverter3D test results.

    y_true / y_pred : (n_samples, n_stations, 2*n_layers-1)
    Shows the first test survey's resistivity section.
    """
    cmap = "viridis_r"
    rho_t = y_true[0, :, :n_layers]  # (n_stations, n_layers)
    rho_p = y_pred[0, :, :n_layers]

    fig, axes = plt.subplots(1, 3, figsize=(14, 5))
    c_t = "#89b4fa" if dark else "#1e66f5"
    c_p = "#f38ba8" if dark else "#d20f39"

    # Mean depth profile
    axes[0].plot(rho_t.mean(axis=0), np.arange(n_layers), color=c_t, lw=2, label="True")
    axes[0].plot(
        rho_p.mean(axis=0),
        np.arange(n_layers),
        color=c_p,
        lw=2,
        ls="--",
        label="Predicted",
    )
    axes[0].invert_yaxis()
    axes[0].set_xlabel("log₁₀(ρ)")
    axes[0].set_ylabel("Layer")
    axes[0].set_title("Mean depth profile")
    axes[0].legend(fontsize=9)

    vmin = min(rho_t.min(), rho_p.min())
    vmax = max(rho_t.max(), rho_p.max())

    im1 = axes[1].imshow(
        rho_t.T,
        aspect="auto",
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        origin="upper",
    )
    axes[1].set_title("True log₁₀(ρ) section")
    axes[1].set_xlabel("Station")
    axes[1].set_ylabel("Layer")
    fig.colorbar(im1, ax=axes[1], label="log₁₀(ρ)")

    im2 = axes[2].imshow(
        rho_p.T,
        aspect="auto",
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        origin="upper",
    )
    axes[2].set_title("Predicted log₁₀(ρ) section")
    axes[2].set_xlabel("Station")
    fig.colorbar(im2, ax=axes[2], label="log₁₀(ρ)")

    fig.suptitle(
        f"GCN 3-D  ·  MAE = {np.abs(y_pred[..., :n_layers] - y_true[..., :n_layers]).mean():.3f}",
        fontsize=11,
    )
    fig.tight_layout()
    return fig


def _ai_stats_div(
    title: str,
    *,
    solver,
    n_layers,
    n_train,
    epochs,
    n_test,
    extra: list | None = None,
) -> html.Div:
    rows = [
        f"Solver: {solver}",
        f"N layers: {n_layers}",
        f"Training samples: {n_train}",
        f"Epochs: {epochs}",
        f"Test samples: {n_test}",
    ] + (extra or [])
    return html.Div(
        [
            html.Div(
                title,
                style={
                    "fontWeight": "600",
                    "color": "#a6e3a1",
                    "marginBottom": "6px",
                },
            ),
            *[html.Div(r, style={"fontSize": "11px"}) for r in rows],
        ],
        style={"padding": "10px"},
    )


def _plot_ai_convergence(history: dict, dark: bool):
    c_t = "#a6e3a1" if dark else "#40a02b"
    c_v = "#f38ba8" if dark else "#d20f39"
    t_loss = history.get("train_loss") or history.get("loss") or []
    v_loss = history.get("val_loss") or history.get("val") or []
    if not t_loss:
        return None
    epochs = np.arange(1, len(t_loss) + 1)
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.semilogy(epochs, t_loss, color=c_t, lw=2, label="Train loss")
    if v_loss:
        ax.semilogy(
            np.arange(1, len(v_loss) + 1),
            v_loss,
            color=c_v,
            lw=2,
            ls="--",
            label="Val loss",
        )
    ax.set_xlabel("Epoch")
    ax.set_ylabel("Loss (MSE)")
    ax.set_title("AI Training Convergence")
    ax.legend(fontsize=9)
    fig.tight_layout()
    return fig


def _plot_ai_comparison(y_true, y_pred, n_layers, dark):
    c_t = "#89b4fa" if dark else "#1e66f5"
    c_p = "#f38ba8" if dark else "#d20f39"
    n_show = min(4, len(y_true))
    fig, axes = plt.subplots(1, n_show, figsize=(4 * n_show, 5))
    if n_show == 1:
        axes = [axes]
    n_rho = n_layers
    n_thk = n_layers - 1
    for i, ax in enumerate(axes[:n_show]):
        yt, yp = y_true[i], y_pred[i]
        rho_t = 10 ** yt[:n_rho]
        thk_t = 10 ** yt[n_rho : n_rho + n_thk]
        rho_p = 10 ** yp[:n_rho]
        thk_p = 10 ** yp[n_rho : n_rho + n_thk]
        tops_t = np.concatenate([[0.0], np.cumsum(thk_t)])
        tops_p = np.concatenate([[0.0], np.cumsum(thk_p)])
        bot = max(tops_t[-1], tops_p[-1]) * 1.2

        def _step(rho, tops):
            d, r = [], []
            for j in range(len(rho) - 1):
                d += [tops[j], tops[j + 1]]
                r += [rho[j], rho[j]]
            d += [tops[-1], bot]
            r += [rho[-1], rho[-1]]
            return r, d

        r_t, d_t = _step(rho_t, tops_t)
        r_p, d_p = _step(rho_p, tops_p)
        ax.plot(r_t, d_t, color=c_t, lw=2, label="True")
        ax.plot(r_p, d_p, color=c_p, lw=2, ls="--", label="Predicted")
        ax.set_xscale("log")
        ax.invert_yaxis()
        ax.set_xlabel("ρ (Ω·m)")
        if i == 0:
            ax.set_ylabel("Depth (m)")
            ax.legend(fontsize=8)
        ax.set_title(f"Sample {i + 1}")
    fig.suptitle("AI 1-D Inversion: True vs Predicted", fontsize=11)
    fig.tight_layout()
    return fig


# ---


def register_inversion(app) -> None:

    # ---

    # G1. Guard Run button: disable for Hybrid
    #     when no checkpoint is loaded
    @app.callback(
        Output(IDs.BTN_INV_RUN, "disabled"),
        Output(IDs.BTN_INV_RUN, "title"),
        Input(IDs.INV_ACTIVE_TRACK, "data"),
        Input(IDs.INV_HYB_CHECKPOINT, "contents"),
    )
    def guard_run_btn(track, checkpoint):
        if (track or "trad") == "hybrid" and not checkpoint:
            return (
                True,
                "Upload a .npz AI checkpoint to enable Hybrid inversion",
            )
        return False, ""

    # G2. Clientside: show running message
    #     immediately on button click
    app.clientside_callback(
        """
        function(n, track) {
            if (!n) return "";
            var m = {
                "trad":
                    "Running Traditional inversion...",
                "ai":
                    "Running AI Neural inversion...",
                "pinn":
                    "Running PINN inversion...",
                "hybrid":
                    "Running Hybrid inversion..."
            };
            return m[track || "trad"] || "Running...";
        }
        """,
        Output(
            IDs.INV_RUN_MSG,
            "children",
            allow_duplicate=True,
        ),
        Input(IDs.BTN_INV_RUN, "n_clicks"),
        State(IDs.INV_ACTIVE_TRACK, "data"),
        prevent_initial_call=True,
    )

    # LS. Line / station selector helpers
    # ---
    # Key-name helpers: station_records may use
    # uppercase ("Line","ID") or lowercase ("line","station")

    def _rec_line(r):
        return r.get("Line") or r.get("line") or ""

    def _rec_sta(r):
        return r.get("ID") or r.get("station") or ""

    def _line_opts(records):
        seen, out = set(), []
        for r in records:
            ln = _rec_line(r)
            if ln and ln not in seen:
                seen.add(ln)
                out.append(ln)
        return [{"label": ln, "value": ln} for ln in sorted(out)]

    def _sta_opts(records, sel_lines=None):
        if sel_lines:
            records = [r for r in records if _rec_line(r) in set(sel_lines)]
        seen, out = set(), []
        for r in records:
            s = _rec_sta(r)
            if s and s not in seen:
                seen.add(s)
                out.append(s)
        return [{"label": s, "value": s} for s in sorted(out)]

    # LS-1. PINN: populate line + station from STORE_DATA
    @app.callback(
        Output(IDs.INV_PINN_LINE_SEL, "options"),
        Output(
            IDs.INV_PINN_STATION_SEL,
            "options",
            allow_duplicate=True,
        ),
        Input(IDs.STORE_DATA, "data"),
        prevent_initial_call=True,
    )
    def populate_pinn_selectors(store_data):
        recs = (store_data or {}).get("station_records", [])
        return _line_opts(recs), _sta_opts(recs)

    # LS-2. PINN: filter stations by selected lines
    @app.callback(
        Output(
            IDs.INV_PINN_STATION_SEL,
            "options",
            allow_duplicate=True,
        ),
        Input(IDs.INV_PINN_LINE_SEL, "value"),
        State(IDs.STORE_DATA, "data"),
        prevent_initial_call=True,
    )
    def filter_pinn_stations(sel_lines, store_data):
        recs = (store_data or {}).get("station_records", [])
        return _sta_opts(recs, sel_lines)

    # LS-3. Hybrid: populate line + station from STORE_DATA
    @app.callback(
        Output(IDs.INV_HYB_LINE_SEL, "options"),
        Output(
            IDs.INV_HYB_STATION_SEL,
            "options",
            allow_duplicate=True,
        ),
        Input(IDs.STORE_DATA, "data"),
        prevent_initial_call=True,
    )
    def populate_hyb_selectors(store_data):
        recs = (store_data or {}).get("station_records", [])
        return _line_opts(recs), _sta_opts(recs)

    # LS-4. Hybrid: filter stations by selected lines
    @app.callback(
        Output(
            IDs.INV_HYB_STATION_SEL,
            "options",
            allow_duplicate=True,
        ),
        Input(IDs.INV_HYB_LINE_SEL, "value"),
        State(IDs.STORE_DATA, "data"),
        prevent_initial_call=True,
    )
    def filter_hyb_stations(sel_lines, store_data):
        recs = (store_data or {}).get("station_records", [])
        return _sta_opts(recs, sel_lines)

    # LS-5. Data-fit: populate line options from
    #        stations that were actually inverted
    @app.callback(
        Output(IDs.INV_DATAFIT_LINE, "options"),
        Output("inv-datafit-selectors", "style"),
        Output("inv-datafit-hint", "style"),
        Input(IDs.INV_DATAFIT_STA, "options"),
        State(IDs.STORE_DATA, "data"),
        prevent_initial_call=True,
    )
    def populate_datafit_lines(sta_opts, store_data):
        if not sta_opts:
            return (
                [],
                {"display": "none"},
                {},
            )
        sta_names = {o["value"] for o in sta_opts}
        recs = (store_data or {}).get("station_records", [])
        active_recs = [r for r in recs if _rec_sta(r) in sta_names]
        lines = _line_opts(active_recs)
        sel_style = {} if lines else {"display": "none"}
        hint_style = {"display": "none"} if sta_opts else {}
        return lines, sel_style, hint_style

    # LS-6. Data-fit: filter stations by selected line
    @app.callback(
        Output(
            IDs.INV_DATAFIT_STA,
            "options",
            allow_duplicate=True,
        ),
        Input(IDs.INV_DATAFIT_LINE, "value"),
        State(IDs.INV_DATAFIT_STA, "options"),
        State(IDs.STORE_DATA, "data"),
        prevent_initial_call=True,
    )
    def filter_datafit_stations(sel_line, current_opts, store_data):
        if not sel_line or not current_opts:
            return no_update
        recs = (store_data or {}).get("station_records", [])
        line_stas = {_rec_sta(r) for r in recs if _rec_line(r) == sel_line}
        filtered = [o for o in current_opts if o["value"] in line_stas]
        return filtered or no_update

    # T1. Track pill buttons → INV_ACTIVE_TRACK store
    @app.callback(
        Output(IDs.INV_ACTIVE_TRACK, "data"),
        [Input(f"inv-track-{t}", "n_clicks") for t in _INV_TRACK_IDS],
        prevent_initial_call=True,
    )
    def switch_inv_track(*_):
        tid = ctx.triggered_id
        if not tid:
            return no_update
        return tid.replace("inv-track-", "")

    # T2. Output-tab buttons → INV_ACTIVE_OUT store
    @app.callback(
        Output(IDs.INV_ACTIVE_OUT, "data"),
        [Input(f"inv-out-tab-{t}", "n_clicks") for t in _INV_OUT_IDS],
        prevent_initial_call=True,
    )
    def switch_inv_out(*_):
        tid = ctx.triggered_id
        if not tid:
            return no_update
        return tid.replace("inv-out-tab-", "")

    # T3. Store → sidebar ctrl panels + track btn classes
    clientside_callback(
        """
        function(active_track) {
            if (!active_track)
                return Array(8).fill(
                    window.dash_clientside.no_update
                );
            var tracks = ["trad","ai","pinn","hybrid"];
            var blk  = {display:"block"};
            var none = {display:"none"};
            var styles = tracks.map(function(t){
                return active_track === t ? blk : none;
            });
            var classes = tracks.map(function(t){
                return "fwd-dim-btn"
                    + (active_track === t
                        ? " active" : "");
            });
            return styles.concat(classes);
        }
        """,
        Output("inv-ctrl-trad", "style"),
        Output("inv-ctrl-ai", "style"),
        Output("inv-ctrl-pinn", "style"),
        Output("inv-ctrl-hybrid", "style"),
        Output("inv-track-trad", "className"),
        Output("inv-track-ai", "className"),
        Output("inv-track-pinn", "className"),
        Output("inv-track-hybrid", "className"),
        Input(IDs.INV_ACTIVE_TRACK, "data"),
        prevent_initial_call=False,
    )

    # T4. Output-store → view panels + output tab classes
    clientside_callback(
        """
        function(active_out) {
            if (!active_out)
                return Array(10).fill(
                    window.dash_clientside.no_update
                );
            var tabs = [
                "result","conv","stats",
                "log","data-fit"
            ];
            var flex = {
                display:"flex",
                flexDirection:"column",
                height:"100%",
                minHeight:"0",
                flex:"1"
            };
            var blk  = {display:"block"};
            var none = {display:"none"};
            var noFlex = ["stats","log","data-fit"];
            var styles = tabs.map(function(t){
                if (active_out !== t) return none;
                return noFlex.indexOf(t) >= 0
                    ? blk : flex;
            });
            var classes = tabs.map(function(t){
                return "prof-tab-btn"
                    + (active_out === t
                        ? " active" : "");
            });
            return styles.concat(classes);
        }
        """,
        Output("inv-panel-result", "style"),
        Output("inv-panel-conv", "style"),
        Output("inv-panel-stats", "style"),
        Output("inv-panel-log", "style"),
        Output("inv-panel-data-fit", "style"),
        Output("inv-out-tab-result", "className"),
        Output("inv-out-tab-conv", "className"),
        Output("inv-out-tab-stats", "className"),
        Output("inv-out-tab-log", "className"),
        Output("inv-out-tab-data-fit", "className"),
        Input(IDs.INV_ACTIVE_OUT, "data"),
        prevent_initial_call=False,
    )

    # T5. Context info bar
    @app.callback(
        Output(IDs.INV_CTX_BAR, "children"),
        Input(IDs.INV_ACTIVE_TRACK, "data"),
        Input(IDs.INV_METHOD, "value"),
        Input(IDs.INV_AI_DIM, "value"),
        Input(IDs.INV_PINN_DIM, "value"),
        Input(IDs.INV_HYB_DIM, "value"),
    )
    def update_inv_ctx_bar(
        active_track,
        method,
        ai_dim,
        pinn_dim,
        hyb_dim,
    ):
        active_track = active_track or "trad"
        icon = _INV_TRACK_ICON.get(active_track, "bi-cpu")
        label = _INV_TRACK_LABEL.get(active_track, "Traditional")
        parts = [
            html.Span(
                [
                    html.I(className=f"bi {icon} me-1"),
                    label,
                ],
                className="adv-info-group",
            ),
        ]
        if active_track == "trad" and method:
            parts += [
                html.Span(".", className="adv-info-sep"),
                html.Span(
                    method.upper(),
                    className="adv-info-plot",
                ),
            ]
        elif active_track == "ai" and ai_dim:
            dim_label = {
                "1d": "1-D ResNet",
                "2d": "2-D UNet",
                "3d": "3-D GCN",
            }.get(ai_dim, ai_dim)
            parts += [
                html.Span(".", className="adv-info-sep"),
                html.Span(
                    dim_label,
                    className="adv-info-plot",
                ),
            ]
        elif active_track == "pinn" and pinn_dim:
            parts += [
                html.Span(".", className="adv-info-sep"),
                html.Span(
                    f"PINN {pinn_dim.upper()}",
                    className="adv-info-plot",
                ),
            ]
        elif active_track == "hybrid" and hyb_dim:
            parts += [
                html.Span(".", className="adv-info-sep"),
                html.Span(
                    f"Hybrid {hyb_dim.upper()}",
                    className="adv-info-plot",
                ),
            ]
        parts += [
            html.Span(".", className="adv-info-sep"),
            html.Span(
                "Configure parameters then click Run Inversion",
                style={
                    "color": "var(--sub0)",
                    "fontSize": "10.5px",
                },
            ),
        ]
        return parts

    # ---

    # 1. Show/hide data-source panels
    @app.callback(
        Output(IDs.INV_SYN_PANEL, "style"),
        Output(IDs.INV_SESS_PANEL, "style"),
        Input(IDs.INV_DATA_SRC, "value"),
    )
    def sync_data_src(src):
        src = src or "synthetic"
        return (
            _SHOW if src == "synthetic" else _HIDE,
            _SHOW if src == "session" else _HIDE,
        )

    # 2. Show/hide 2D station/profile params
    @app.callback(
        Output(IDs.INV_2D_PARAMS, "style"),
        Input(IDs.INV_T_DIM, "value"),
    )
    def sync_t_dim(dim):
        return _SHOW if (dim or "1d") == "2d" else _HIDE

    # 3. Update AI architecture options + show/hide dim-specific panels
    @app.callback(
        Output(IDs.INV_AI_ARCH, "options"),
        Output(IDs.INV_AI_ARCH, "value"),
        Output("inv-ai-arch-note", "children"),
        Output(IDs.INV_AI_2D_PANEL, "style"),
        Output(IDs.INV_AI_3D_PANEL, "style"),
        Input(IDs.INV_AI_DIM, "value"),
    )
    def sync_ai_dim(ai_dim):
        ai_dim = ai_dim or "1d"
        show_2d = _SHOW if ai_dim == "2d" else _HIDE
        show_3d = _SHOW if ai_dim == "3d" else _HIDE
        if ai_dim == "1d":
            opts = [
                {"label": "ResNet1D  [recommended]", "value": "resnet"},
                {"label": "CNN1D", "value": "cnn1d"},
                {"label": "FCN1D", "value": "fcn"},
            ]
            return (
                opts,
                "resnet",
                "1-D sounding inversion from MT/CSAMT/TEM response.",
                show_2d,
                show_3d,
            )
        elif ai_dim == "2d":
            opts = [{"label": "UNet2D  (encoder–decoder)", "value": "unet2d"}]
            return (
                opts,
                "unet2d",
                "2-D section from multi-station MT response (U-Net encoder–decoder).",
                show_2d,
                show_3d,
            )
        else:
            opts = [{"label": "GCN3D  (Graph Conv. Net)", "value": "gcn3d"}]
            return (
                opts,
                "gcn3d",
                "3-D quasi-inversion via graph-based station connectivity (GCN).",
                show_2d,
                show_3d,
            )

    # 4. Add layer to true model table
    @app.callback(
        Output(IDs.INV_TRUE_TABLE, "data"),
        Input(IDs.BTN_INV_ADD_LAYER, "n_clicks"),
        State(IDs.INV_TRUE_TABLE, "data"),
        prevent_initial_call=True,
    )
    def add_true_layer(n_clicks, data):
        if not n_clicks:
            raise PreventUpdate
        data = list(data or [])
        last = float(data[-1]["resistivity"]) if data else 100.0
        data.append({"layer": len(data) + 1, "resistivity": last, "thickness": 500})
        return data

    # 5. ── PyTorch install: open confirm modal when AI run without torch ────────
    # (handled via run_inversion returning INV_INSTALL_MODAL is_open=True)

    # 6. Start PyTorch installation in background thread
    @app.callback(
        Output(IDs.INV_INSTALL_INTERVAL, "disabled"),
        Output(IDs.INV_INSTALL_PROG_AREA, "style"),
        Output(IDs.BTN_INV_INSTALL_YES, "disabled"),
        Output(IDs.BTN_INV_INSTALL_NO, "disabled"),
        Output(IDs.INV_INSTALL_LABEL, "children"),
        Input(IDs.BTN_INV_INSTALL_YES, "n_clicks"),
        prevent_initial_call=True,
    )
    def start_install(n_clicks):
        if not n_clicks:
            raise PreventUpdate
        # Launch background installer thread
        t = threading.Thread(target=_run_pip_install, daemon=True)
        t.start()
        return (
            False,  # enable interval
            _SHOW,  # show progress area
            True,  # disable Install button
            True,  # disable Cancel button (install running)
            "Resolving dependencies…",
        )

    # 7. Cancel / close install modal
    @app.callback(
        Output(IDs.INV_INSTALL_MODAL, "is_open", allow_duplicate=True),
        Input(IDs.BTN_INV_INSTALL_NO, "n_clicks"),
        prevent_initial_call=True,
    )
    def cancel_install(n_clicks):
        if not n_clicks:
            raise PreventUpdate
        return False

    # 8. Poll install progress (fires every 800ms while interval is enabled)
    @app.callback(
        Output(IDs.INV_INSTALL_PROGRESS, "value"),
        Output(IDs.INV_INSTALL_PROGRESS, "animated"),
        Output(IDs.INV_INSTALL_PROGRESS, "striped"),
        Output(IDs.INV_INSTALL_PROGRESS, "color"),
        Output(IDs.INV_INSTALL_LOG, "children"),
        Output(IDs.INV_INSTALL_LABEL, "children", allow_duplicate=True),
        Output(IDs.BTN_INV_INSTALL_YES, "disabled", allow_duplicate=True),
        Output(IDs.BTN_INV_INSTALL_NO, "disabled", allow_duplicate=True),
        Output(IDs.INV_INSTALL_MODAL, "is_open", allow_duplicate=True),
        Output(IDs.INV_INSTALL_INTERVAL, "disabled", allow_duplicate=True),
        Output(IDs.INV_FEEDBACK, "children", allow_duplicate=True),
        Output(IDs.TOAST_SUCCESS, "is_open", allow_duplicate=True),
        Output(IDs.TOAST_SUCCESS_BODY, "children", allow_duplicate=True),
        Input(IDs.INV_INSTALL_INTERVAL, "n_intervals"),
        prevent_initial_call=True,
    )
    def poll_install(n_intervals):
        with _install_lock:
            state = dict(_install_state)

        status = state.get("status", "idle")
        if status == "idle":
            raise PreventUpdate

        progress = state.get("progress", 0)
        log_lines = state.get("log_lines", [])
        log_text = "\n".join(log_lines[-40:])  # last 40 lines

        if status == "running":
            label = f"Installing…  {progress}%"
            return (
                progress,
                True,
                True,
                "info",
                log_text,
                label,
                True,
                True,  # keep buttons disabled
                no_update,
                no_update,  # keep modal open, keep interval running
                no_update,
                False,
                "",
            )

        elif status == "done":
            # Show success state briefly, then auto-close
            elapsed = time.time() - state.get("done_time", 0)
            label = "✓ PyTorch installed successfully!"
            if elapsed < 2.5:
                # Still in "show success" window — progress bar full + green
                return (
                    100,
                    False,
                    False,
                    "success",
                    log_text,
                    label,
                    True,
                    True,  # buttons still disabled
                    no_update,
                    no_update,
                    no_update,
                    False,
                    "",
                )
            else:
                # Close modal, disable interval, show success toast + feedback
                feedback = "✓ PyTorch ready — AI inversion enabled. Run again to start."
                return (
                    100,
                    False,
                    False,
                    "success",
                    log_text,
                    label,
                    False,
                    False,  # re-enable buttons
                    False,
                    True,  # close modal, stop interval
                    feedback,
                    True,
                    "PyTorch installed successfully! Restart the app if inference fails.",
                )

        else:  # error
            label = "✗ Installation failed. See log below."
            return (
                progress,
                False,
                False,
                "danger",
                log_text,
                label,
                False,
                False,  # re-enable buttons for retry
                no_update,
                True,  # keep modal open, stop interval
                no_update,
                False,
                "",
            )

    # 9. Main inversion runner (trad + AI tracks)
    @app.callback(
        Output(IDs.IMG_INV, "src", allow_duplicate=True),
        Output(IDs.IMG_INV_CONV, "src", allow_duplicate=True),
        Output(IDs.INV_LOG, "children", allow_duplicate=True),
        Output(IDs.INV_STATS, "children", allow_duplicate=True),
        Output(IDs.INV_SPINNER, "children", allow_duplicate=True),
        Output(IDs.INV_FEEDBACK, "children", allow_duplicate=True),
        Output(IDs.INV_ACTIVE_OUT, "data", allow_duplicate=True),
        Output(IDs.TOAST_ERROR, "is_open", allow_duplicate=True),
        Output(IDs.TOAST_BODY, "children", allow_duplicate=True),
        Output(IDs.INV_INSTALL_MODAL, "is_open", allow_duplicate=True),
        Output(IDs.INV_RUN_MSG, "children", allow_duplicate=True),
        Input(IDs.BTN_INV_RUN, "n_clicks"),
        State(IDs.INV_ACTIVE_TRACK, "data"),
        State(IDs.INV_METHOD, "value"),
        State(IDs.INV_T_DIM, "value"),
        State(IDs.INV_BACKEND, "value"),
        State(IDs.INV_DATA_SRC, "value"),
        State(IDs.INV_TRUE_TABLE, "data"),
        State(IDs.INV_HALFSPACE_RHO, "value"),
        State(IDs.INV_NOISE_LEVEL, "value"),
        State(IDs.INV_N_STATIONS, "value"),
        State(IDs.INV_PROFILE_LEN, "value"),
        State(IDs.INV_N_LAYERS, "value"),
        State(IDs.INV_REGULARIZE, "value"),
        State(IDs.INV_MAX_ITER, "value"),
        State(IDs.INV_ERROR_FLOOR, "value"),
        State(IDs.INV_INCLUDE_PHASE, "value"),
        State(IDs.INV_FREQ_MIN, "value"),
        State(IDs.INV_FREQ_MAX, "value"),
        State(IDs.INV_N_FREQ, "value"),
        State(IDs.INV_AI_DIM, "value"),
        State(IDs.INV_AI_ARCH, "value"),
        State(IDs.INV_AI_SOLVER, "value"),
        State(IDs.INV_AI_N_LAYERS, "value"),
        State(IDs.INV_AI_N_SAMPLES, "value"),
        State(IDs.INV_AI_EPOCHS, "value"),
        State(IDs.INV_AI_BATCH, "value"),
        State(IDs.INV_AI_LR, "value"),
        State(IDs.INV_AI_NOISE, "value"),
        # 2D / 3D specific
        State(IDs.INV_AI_N_COMPONENTS, "value"),
        State(IDs.INV_AI_N_DEPTH, "value"),
        State(IDs.INV_AI_N_STATIONS, "value"),
        State(IDs.INV_AI_N_STATIONS_3D, "value"),
        State(IDs.INV_AI_UNET_DEPTH, "value"),
        State(IDs.SESSION_ID, "data"),
        State(IDs.STORE_THEME, "data"),
        prevent_initial_call=True,
    )
    def run_inversion(
        n_clicks,
        track,
        method,
        t_dim,
        backend,
        data_src,
        true_table,
        halfspace_rho,
        noise_level,
        n_stations,
        profile_len,
        n_layers,
        regularize,
        max_iter,
        error_floor,
        include_phase,
        freq_min,
        freq_max,
        n_freq,
        ai_dim,
        ai_arch,
        ai_solver,
        ai_n_layers,
        ai_n_samples,
        ai_epochs,
        ai_batch,
        ai_lr,
        ai_noise,
        ai_n_components,
        ai_n_depth,
        ai_n_stations,
        ai_n_stations_3d,
        ai_unet_depth,
        session_id,
        theme,
    ):
        if not n_clicks:
            raise PreventUpdate

        track = track or "trad"
        if track in ("pinn", "hybrid"):
            raise PreventUpdate

        dark = (theme or "dark") == "dark"
        if dark:
            apply_web_dark_theme()
        else:
            apply_web_light_theme()

        _empty = empty_src(dark=dark)

        # Convenience: 10-tuple returns with modal=False (no install prompt)
        def _ok(
            src,
            src_conv,
            spans,
            stats,
            msg,
            active_tab="result",
        ):
            return (
                src,
                src_conv,
                spans,
                stats,
                "",
                msg,
                active_tab,
                False,
                "",
                False,
                "",
            )

        def _err(spans, msg):
            return (
                _empty,
                _empty,
                spans,
                no_update,
                "",
                f"Error: {msg}",
                "log",
                True,
                msg,
                False,
                "",
            )

        def _open_install(spans, msg):
            return (
                _empty,
                _empty,
                spans,
                no_update,
                "",
                msg,
                "log",
                False,
                "",
                True,
                "",
            )

        log_lines: list[str] = []

        def _log(m):
            log_lines.append(m)

        f_min = float(freq_min if freq_min is not None else -3)
        f_max = float(freq_max if freq_max is not None else 4)
        n_pts = int(n_freq or 20)
        freqs = np.logspace(f_min, f_max, n_pts)
        np.random.seed(0)

        # ── AI track ──────────────────────────────────────────────────────
        if track == "ai":
            _log("=" * 58)
            _log("  pycsamt AI Neural Inversion")
            _log("=" * 58)

            if not _TORCH_OK:
                _log("INFO  PyTorch not found in current environment.")
                _log("      Click 'Install PyTorch' to install automatically.")
                spans = _coloured_log(log_lines)
                return _open_install(
                    spans,
                    "PyTorch not installed — click Install PyTorch to continue.",
                )

            try:
                ai_dim = (ai_dim or "1d").lower()
                ai_arch = ai_arch or "resnet"
                ai_solver = (ai_solver or "mt1d").lower()
                ai_n_lay = int(ai_n_layers or 4)
                n_samp = int(ai_n_samples or 1000)
                epochs = int(ai_epochs or 40)
                batch_sz = int(ai_batch or 128)
                lr = float(ai_lr or 1e-3)
                noise = float(ai_noise or 0.05)
                n_sta_ai = int(n_stations or 16)

                # ── 1-D ──────────────────────────────────────────────────
                if ai_dim == "1d":
                    from pycsamt.ai.inversion import (
                        EMInverter1D,
                    )
                    from pycsamt.forward.batch import (
                        generate_dataset,
                    )

                    _log(f"Generating {n_samp} 1-D training samples …")
                    ds_train = generate_dataset(
                        solver=ai_solver,
                        n_samples=n_samp,
                        freqs=freqs,
                        n_layers=ai_n_lay,
                        noise_level=noise,
                        seed=42,
                        verbose=False,
                    )
                    ds_test = generate_dataset(
                        solver=ai_solver,
                        n_samples=30,
                        freqs=freqs,
                        n_layers=ai_n_lay,
                        noise_level=noise,
                        seed=99,
                        verbose=False,
                    )
                    _log(f"X_train {ds_train.X.shape}  X_test {ds_test.X.shape}")
                    _log(
                        f"Training EMInverter1D(arch={ai_arch}) "
                        f"— {epochs} epochs, batch {batch_sz}, lr {lr} …"
                    )
                    inv = EMInverter1D(
                        arch=ai_arch, n_layers=ai_n_lay, solver=ai_solver
                    )
                    inv.fit(
                        ds_train.X,
                        ds_train.y,
                        epochs=epochs,
                        batch_size=batch_sz,
                        lr=lr,
                        verbose=False,
                        seed=42,
                    )
                    _log("Training complete.")

                    y_pred = inv.predict(ds_test.X, as_log_rho=True)
                    _log(f"Predicted {len(y_pred)} test samples.")

                    fig_res = _plot_ai_comparison(ds_test.y, y_pred, ai_n_lay, dark)
                    src_res = fig_to_src(fig_res)
                    plt.close(fig_res)

                    src_conv = _empty
                    hist = getattr(inv, "history_", None) or getattr(
                        inv, "_history", None
                    )
                    if hist:
                        fig_c = _plot_ai_convergence(hist, dark)
                        if fig_c:
                            src_conv = fig_to_src(fig_c)
                            plt.close(fig_c)

                    _log("Done.")
                    stats = _ai_stats_div(
                        f"AI 1-D  ({ai_arch.upper()})",
                        solver=ai_solver.upper(),
                        n_layers=ai_n_lay,
                        n_train=n_samp,
                        epochs=epochs,
                        n_test=len(y_pred),
                    )
                    return _ok(
                        src_res,
                        src_conv,
                        _coloured_log(log_lines),
                        stats,
                        "AI 1-D inversion complete",
                    )

                # ── 2-D ──────────────────────────────────────────────────
                elif ai_dim == "2d":
                    from pycsamt.ai.inversion import (
                        EMInverter2D,
                    )
                    from pycsamt.app.web.cache import (
                        cache_get_fwd,
                    )

                    # Resolve 2D-specific user params
                    n_comp = int(ai_n_components or 2)
                    n_dep = int(ai_n_depth or ai_n_lay)
                    unet_d = (
                        None
                        if (ai_unet_depth or "auto") == "auto"
                        else int(ai_unet_depth)
                    )

                    # Check for a cached ForwardResponse2D from the forward page
                    fwd_resp = cache_get_fwd(session_id, "2d")
                    if fwd_resp is not None:
                        n_freqs = fwd_resp.n_freqs
                        n_sta_ai = fwd_resp.n_stations
                        _log(
                            f"Using cached ForwardResponse2D: "
                            f"{n_freqs} freqs × {n_sta_ai} stations."
                        )
                        X_te_fwd = _fwd2d_to_tensor(fwd_resp)  # (1, 2, nf, ns)
                        has_fwd_test = True
                    else:
                        n_freqs = len(freqs)
                        n_sta_ai = int(ai_n_stations or 8)
                        has_fwd_test = False
                        _log("No forward model in session — using synthetic data.")

                    _log(
                        f"Generating {n_samp} synthetic 2-D profiles "
                        f"({n_sta_ai} stations, {n_dep} depth cells, "
                        f"{n_comp} components) …"
                    )
                    X_tr, y_tr = _generate_2d_dataset(
                        n_samp, n_sta_ai, n_dep, freqs, noise, seed=42
                    )
                    X_te, y_te = _generate_2d_dataset(
                        20, n_sta_ai, n_dep, freqs, noise, seed=99
                    )
                    _log(
                        f"X_train {X_tr.shape}  y_train {y_tr.shape}  "
                        f"X_test {X_te.shape}"
                    )

                    inv = EMInverter2D(
                        n_components=n_comp,
                        n_depth=n_dep,
                        n_stations=n_sta_ai,
                        n_freqs=n_freqs,
                        unet_depth=unet_d,
                        dropout=0.1,
                    )
                    unet_info = (
                        f"auto ({len(inv._channels) - 1} stages)"
                        if unet_d is None
                        else f"{unet_d} stages (manual)"
                    )
                    _log(
                        f"UNet channels: {inv._channels}  [{unet_info}]  "
                        f"input {n_freqs}f × {n_sta_ai}s"
                    )
                    inv.fit(
                        X_tr,
                        y_tr,
                        epochs=epochs,
                        batch_size=batch_sz,
                        lr=lr,
                        patience=max(8, epochs // 5),
                        verbose=False,
                        seed=42,
                    )
                    _log("Training complete.")

                    y_pred = inv.predict(X_te, as_log_rho=True)
                    mae_synth = float(np.abs(y_pred - y_te).mean())
                    _log(
                        f"Predicted {len(y_pred)} synthetic test profiles. "
                        f"MAE = {mae_synth:.3f}"
                    )

                    # If we have a real forward response, also predict on it
                    extra_lines = [
                        f"Input channels (n_components): {n_comp}",
                        f"Stations per profile: {n_sta_ai}",
                        f"Depth cells (n_depth): {n_dep}",
                        f"Frequencies: {n_freqs}",
                        f"UNet: {inv._channels}  [{unet_info}]",
                        f"Synthetic test MAE: {mae_synth:.3f}",
                    ]
                    if has_fwd_test:
                        y_pred_fwd = inv.predict(X_te_fwd, as_log_rho=True)
                        _log(
                            f"Predicted on session ForwardResponse2D "
                            f"— shape {y_pred_fwd.shape}"
                        )
                        extra_lines.append("Session forward response: predicted ✓")
                        fig_res = _plot_predicted_section_2d(
                            y_pred_fwd[0], fwd_resp, n_dep, dark
                        )
                    else:
                        fig_res = _plot_ai_section_2d(y_te, y_pred, dark)
                    src_res = fig_to_src(fig_res)
                    plt.close(fig_res)

                    src_conv = _empty
                    hist = getattr(inv, "_history", None)
                    if hist:
                        fig_c = _plot_ai_convergence(hist, dark)
                        if fig_c:
                            src_conv = fig_to_src(fig_c)
                            plt.close(fig_c)

                    _log("Done.")
                    stats = _ai_stats_div(
                        f"AI 2-D  (UNet2D  {len(inv._channels) - 1}-stage)",
                        solver="MT 2-D synthetic + "
                        + ("session" if has_fwd_test else "synthetic")
                        + " test",
                        n_layers=n_dep,
                        n_train=n_samp,
                        epochs=epochs,
                        n_test=len(y_pred),
                        extra=extra_lines,
                    )
                    return _ok(
                        src_res,
                        src_conv,
                        _coloured_log(log_lines),
                        stats,
                        "AI 2-D inversion complete",
                    )

                # ── 3-D ──────────────────────────────────────────────────
                else:  # ai_dim == "3d"
                    from pycsamt.ai.inversion import (
                        GCNInverter3D,
                    )
                    from pycsamt.ai.nets.gcn import (
                        build_adjacency,
                    )
                    from pycsamt.forward.batch import (
                        generate_dataset_3d,
                    )

                    n_sta_3d = int(ai_n_stations_3d or 20)
                    _log(
                        f"Generating {n_samp} pseudo-3D surveys "
                        f"({n_sta_3d} stations, {ai_n_lay} layers) …"
                    )
                    ds = generate_dataset_3d(
                        solver="mt1d",
                        n_surveys=n_samp + 20,
                        n_stations=n_sta_3d,
                        n_layers=ai_n_lay,
                        freqs=freqs,
                        noise_level=noise,
                        seed=42,
                        verbose=False,
                    )
                    X_tr, y_tr = ds.X[:n_samp], ds.y[:n_samp]
                    X_te, y_te = ds.X[n_samp:], ds.y[n_samp:]
                    coords = ds.coords
                    _log(f"X_train {X_tr.shape}  X_test {X_te.shape}")

                    A = build_adjacency(
                        coords,
                        radius=float(
                            np.percentile(
                                np.linalg.norm(
                                    coords[:, None] - coords[None, :], axis=-1
                                )[np.triu_indices(n_sta_3d, k=1)],
                                30,
                            )
                        ),
                    )
                    n_features = X_tr.shape[-1]
                    _log(
                        f"Training GCNInverter3D  n_features={n_features} "
                        f"n_layers={ai_n_lay}  {epochs} epochs …"
                    )
                    inv = GCNInverter3D(
                        n_features=n_features,
                        n_layers=ai_n_lay,
                        hidden=(128, 64, 32),
                    )
                    inv.fit(
                        X_tr,
                        y_tr,
                        adjacency=A,
                        epochs=epochs,
                        batch_size=batch_sz,
                        lr=lr,
                        patience=max(8, epochs // 5),
                        verbose=False,
                        seed=42,
                    )
                    _log("Training complete.")

                    y_pred = inv.predict(X_te, adjacency=A, as_log_rho=True)
                    mae = float(
                        np.abs(y_pred[..., :ai_n_lay] - y_te[..., :ai_n_lay]).mean()
                    )
                    _log(
                        f"Predicted {len(y_pred)} test surveys. "
                        f"MAE(log10 ρ) = {mae:.3f}"
                    )

                    fig_res = _plot_gcn_section_3d(y_te, y_pred, ai_n_lay, dark)
                    src_res = fig_to_src(fig_res)
                    plt.close(fig_res)

                    src_conv = _empty
                    hist = getattr(inv, "_history", None)
                    if hist:
                        fig_c = _plot_ai_convergence(hist, dark)
                        if fig_c:
                            src_conv = fig_to_src(fig_c)
                            plt.close(fig_c)

                    _log("Done.")
                    stats = _ai_stats_div(
                        "AI 3-D  (GCN3D)",
                        solver="MT GCN pseudo-3D",
                        n_layers=ai_n_lay,
                        n_train=n_samp,
                        epochs=epochs,
                        n_test=len(y_pred),
                        extra=[
                            f"Stations: {n_sta_3d}",
                            f"Features/station: {n_features}",
                            f"Test MAE: {mae:.3f}",
                        ],
                    )
                    return _ok(
                        src_res,
                        src_conv,
                        _coloured_log(log_lines),
                        stats,
                        "AI 3-D inversion complete",
                    )

            except Exception as exc:
                import traceback as _tb

                _log(f"ERROR  {exc}")
                _log(_tb.format_exc()[-1200:])
                return _err(_coloured_log(log_lines), str(exc))

        # ── Traditional track ──────────────────────────────────────────────
        method = (method or "mt").lower()
        t_dim = (t_dim or "1d").lower()
        backend = (backend or "builtin").lower()
        data_src = (data_src or "synthetic").lower()
        n_layers = int(n_layers or 4)
        regularize = regularize or "smooth"
        max_iter = int(max_iter or 60)
        error_flr = float(error_floor or 0.05)
        inc_phase = bool(include_phase)

        _log("=" * 58)
        _log(
            f"  Traditional Inversion  "
            f"method={method.upper()} dim={t_dim} backend={backend}"
        )
        _log("=" * 58)

        if backend != "builtin":
            msg = (
                f"'{backend}' requires an external binary. "
                "Only the built-in (SciPy) backend is available here."
            )
            _log(f"WARNING  {msg}")
            return _ok(
                _empty,
                _empty,
                _coloured_log(log_lines),
                html.Div(msg, style={"color": "#f9e2af", "padding": "10px"}),
                msg,
                "log",
            )

        try:
            from pycsamt.inversion.config import (
                InversionConfig,
            )
            from pycsamt.inversion.workflow import (
                InversionWorkflow,
            )

            model_true = None

            if data_src == "session":
                sites = cache_get(session_id)
                if not sites:
                    msg = "No session data. Load EDI files first."
                    _log(f"WARNING  {msg}")
                    return _ok(
                        _empty,
                        _empty,
                        _coloured_log(log_lines),
                        html.Div(msg, style={"color": "#f9e2af", "padding": "10px"}),
                        msg,
                        "log",
                    )
                from pycsamt.inversion.data import EMData

                all_rho, all_phase, st_x, st_names = [], [], [], []
                for i, site in enumerate(sites[:20]):
                    try:
                        freqs_s = getattr(site, "freqs", freqs)
                        rho_s = getattr(site, "rho_a", None)
                        phi_s = getattr(site, "phase", None)
                        if rho_s is None:
                            continue
                        all_rho.append(
                            np.interp(np.log10(freqs), np.log10(freqs_s), rho_s)
                        )
                        all_phase.append(
                            np.interp(np.log10(freqs), np.log10(freqs_s), phi_s)
                        )
                        st_x.append(getattr(site, "x", i * 200.0))
                        st_names.append(getattr(site, "station", f"S{i:02d}"))
                    except Exception:
                        continue
                if not all_rho:
                    return _err(
                        _coloured_log(log_lines),
                        "Could not extract MT response from session data.",
                    )
                rho_arr = np.array(all_rho)
                phase_arr = np.array(all_phase)
                em_data = EMData(
                    method="mt",
                    frequencies=freqs,
                    rho_a=rho_arr if rho_arr.ndim > 1 else rho_arr[0],
                    phase=phase_arr if phase_arr.ndim > 1 else phase_arr[0],
                    station_names=st_names,
                    station_x=np.array(st_x) if len(st_x) > 1 else None,
                )
                _log(f"Session data: {em_data.n_stations} station(s).")

            else:
                model_true = _build_true_model(true_table, halfspace_rho)
                _log(
                    f"True model: {model_true.n_layers} layers "
                    f"ρ={list(np.round(model_true.resistivity, 1))}"
                )
                if t_dim == "2d":
                    em_data = _make_synthetic_2d(
                        model_true,
                        method,
                        freqs,
                        noise_level,
                        n_stations,
                        profile_len,
                    )
                    _log(f"Synthetic 2-D: {em_data.n_stations} stations.")
                else:
                    em_data = _make_synthetic_1d(model_true, method, freqs, noise_level)
                    _log("Synthetic 1-D data generated.")

            _log(
                f"InversionConfig: n_layers={n_layers} reg={regularize} "
                f"max_iter={max_iter} …"
            )
            cfg = InversionConfig(
                method=method,
                dimension=t_dim,
                backend=backend,
                data=em_data,
                n_layers=n_layers,
                regularization=regularize,
                max_iter=max_iter,
                error_floor=error_flr,
                include_phase=inc_phase,
            )
            result = InversionWorkflow(cfg).run()
            try:
                from pycsamt.app.web.cache import (
                    cache_set_inversion_result,
                )

                cache_set_inversion_result(session_id, result)
            except Exception:
                pass

            _log(
                f"Status: {result.status}   RMS: {result.rms:.4f}   "
                f"N_iter: {result.n_iter}"
            )
            for w in result.warnings or []:
                _log(f"WARN  {w}")
            _log("=" * 58 + "\nDone.")

            fig_res = (
                _plot_2d_result(result, dark)
                if t_dim == "2d"
                else _plot_1d_result(model_true, result, dark)
            )
            src_res = fig_to_src(fig_res)
            plt.close(fig_res)

            fig_cv = _plot_convergence(result, dark)
            src_cv = (fig_to_src(fig_cv), plt.close(fig_cv))[0] if fig_cv else _empty

            stats = _stats_table(result, method, t_dim, backend, n_layers, regularize)
            msg = (
                f"{method.upper()} {t_dim.upper()} converged — "
                f"RMS={result.rms:.3f}, {result.n_iter} iterations"
            )
            return _ok(src_res, src_cv, _coloured_log(log_lines), stats, msg)

        except Exception as exc:
            import traceback as _tb

            _log(f"ERROR  {exc}")
            _log(_tb.format_exc()[-800:])
            return _err(_coloured_log(log_lines), str(exc))

    # ---

    @app.callback(
        Output("inv-pinn-solver-card", "style"),
        Output(IDs.INV_PINN_2D_PANEL, "style"),
        Output(IDs.INV_PINN_3D_PANEL, "style"),
        Output("inv-pinn-lamx-row", "style"),
        Output("inv-pinn-lamg-row", "style"),
        Input(IDs.INV_PINN_DIM, "value"),
    )
    def sync_pinn_dim(dim):
        dim = dim or "1d"
        solver = _SHOW if dim == "1d" else _HIDE
        p2d = _SHOW if dim == "2d" else _HIDE
        p3d = _SHOW if dim == "3d" else _HIDE
        lamx = _SHOW if dim != "1d" else _HIDE
        lamg = _SHOW if dim == "3d" else _HIDE
        return solver, p2d, p3d, lamx, lamg

    @app.callback(
        Output(IDs.INV_HYB_2D_PANEL, "style"),
        Output(IDs.INV_HYB_3D_PANEL, "style"),
        Output("inv-hyb-lamx-row", "style"),
        Output("inv-hyb-lamg-row", "style"),
        Input(IDs.INV_HYB_DIM, "value"),
    )
    def sync_hyb_dim(dim):
        dim = dim or "1d"
        p2d = _SHOW if dim == "2d" else _HIDE
        p3d = _SHOW if dim == "3d" else _HIDE
        lamx = _SHOW if dim != "1d" else _HIDE
        lamg = _SHOW if dim == "3d" else _HIDE
        return p2d, p3d, lamx, lamg

    # ---

    @app.callback(
        Output(IDs.INV_HYB_MODEL_INFO, "children"),
        Input(IDs.INV_HYB_CHECKPOINT, "contents"),
        State(IDs.INV_HYB_CHECKPOINT, "filename"),
        prevent_initial_call=True,
    )
    def hyb_checkpoint_info(contents, filename):
        if not contents:
            raise PreventUpdate
        try:
            import numpy as _np

            from pycsamt.app.web.utils_pinn import (
                decode_npz_checkpoint,
            )

            path = decode_npz_checkpoint(contents)
            with _np.load(path) as npz:
                keys = list(npz.keys())
            import os

            os.unlink(path)
            return html.Span(
                f"{filename}  —  keys: " + ", ".join(keys[:6]),
                style={
                    "fontSize": "11px",
                    "color": "#a6e3a1",
                },
            )
        except Exception as exc:
            return html.Span(
                f"Could not read: {exc}",
                style={
                    "fontSize": "11px",
                    "color": "#f38ba8",
                },
            )

    # ---

    @app.callback(
        Output(IDs.IMG_INV, "src", allow_duplicate=True),
        Output(IDs.IMG_INV_CONV, "src", allow_duplicate=True),
        Output(IDs.INV_LOG, "children", allow_duplicate=True),
        Output(IDs.INV_STATS, "children", allow_duplicate=True),
        Output(IDs.INV_SPINNER, "children", allow_duplicate=True),
        Output(IDs.INV_FEEDBACK, "children", allow_duplicate=True),
        Output(IDs.INV_ACTIVE_OUT, "data", allow_duplicate=True),
        Output(IDs.TOAST_ERROR, "is_open", allow_duplicate=True),
        Output(IDs.TOAST_BODY, "children", allow_duplicate=True),
        Output(IDs.INV_INSTALL_MODAL, "is_open", allow_duplicate=True),
        Output(IDs.INV_PINN_RESULT_STORE, "data"),
        Output(IDs.INV_HYB_RESULT_STORE, "data"),
        Output(IDs.INV_DATAFIT_STA, "options"),
        Output(IDs.INV_RUN_MSG, "children", allow_duplicate=True),
        Input(IDs.BTN_INV_RUN, "n_clicks"),
        State(IDs.INV_ACTIVE_TRACK, "data"),
        # PINN params
        State(IDs.INV_PINN_DIM, "value"),
        State(IDs.INV_PINN_DATA_SRC, "value"),
        State(IDs.INV_PINN_STATION_SEL, "value"),
        State(IDs.INV_PINN_SOLVER, "value"),
        State(IDs.INV_PINN_MODE, "value"),
        State(IDs.INV_PINN_N_LAYERS, "value"),
        State(IDs.INV_PINN_DEPTH_MAX, "value"),
        State(IDs.INV_PINN_N_FREQS, "value"),
        State(IDs.INV_PINN_EPOCHS, "value"),
        State(IDs.INV_PINN_LR, "value"),
        State(IDs.INV_PINN_LOG_EVERY, "value"),
        State(IDs.INV_PINN_DEVICE, "value"),
        State(IDs.INV_PINN_LAM_Z, "value"),
        State(IDs.INV_PINN_LAM_X, "value"),
        State(IDs.INV_PINN_LAM_G, "value"),
        State(IDs.INV_PINN_RADIUS, "value"),
        State(IDs.INV_PINN_SPC, "value"),
        State(IDs.INV_PINN_COMP_TE, "value"),
        State(IDs.INV_PINN_COMP_TM, "value"),
        State(IDs.INV_BACKEND_SELECT, "value"),
        # Hybrid params
        State(IDs.INV_HYB_DIM, "value"),
        State(IDs.INV_HYB_DATA_SRC, "value"),
        State(IDs.INV_HYB_STATION_SEL, "value"),
        State(IDs.INV_HYB_CHECKPOINT, "contents"),
        State(IDs.INV_HYB_MODE, "value"),
        State(IDs.INV_HYB_N_LAYERS, "value"),
        State(IDs.INV_HYB_N_FREQS, "value"),
        State(IDs.INV_HYB_EPOCHS, "value"),
        State(IDs.INV_HYB_LR, "value"),
        State(IDs.INV_HYB_DEVICE, "value"),
        State(IDs.INV_HYB_LAM_Z, "value"),
        State(IDs.INV_HYB_COMP_TE, "value"),
        State(IDs.INV_HYB_COMP_TM, "value"),
        State(IDs.SESSION_ID, "data"),
        State(IDs.STORE_THEME, "data"),
        prevent_initial_call=True,
    )
    def run_pinn_hybrid_inversion(
        n_clicks,
        track,
        pinn_dim,
        pinn_src,
        pinn_stations,
        pinn_solver,
        pinn_mode,
        pinn_n_layers,
        pinn_depth_max,
        pinn_n_freqs,
        pinn_epochs,
        pinn_lr,
        pinn_log_every,
        pinn_device,
        pinn_lam_z,
        pinn_lam_x,
        pinn_lam_g,
        pinn_radius,
        pinn_spc,
        pinn_comp_te,
        pinn_comp_tm,
        pinn_backend,
        hyb_dim,
        hyb_src,
        hyb_stations,
        hyb_ckpt,
        hyb_mode,
        hyb_n_layers,
        hyb_n_freqs,
        hyb_epochs,
        hyb_lr,
        hyb_device,
        hyb_lam_z,
        hyb_comp_te,
        hyb_comp_tm,
        session_id,
        theme,
    ):
        if not n_clicks:
            raise PreventUpdate
        track = track or "trad"
        if track not in ("pinn", "hybrid"):
            raise PreventUpdate

        dark = (theme or "dark") == "dark"
        if dark:
            apply_web_dark_theme()
        else:
            apply_web_light_theme()

        _empty = empty_src(dark=dark)

        def _ok13(
            src,
            src_conv,
            spans,
            stats,
            msg,
            tab="result",
            pinn_store=None,
            hyb_store=None,
            sta_opts=None,
        ):
            return (
                src,
                src_conv,
                spans,
                stats,
                "",
                msg,
                tab,
                False,
                "",
                False,
                pinn_store,
                hyb_store,
                sta_opts or [],
                "",
            )

        def _err13(spans, msg):
            return (
                _empty,
                _empty,
                spans,
                no_update,
                "",
                f"Error: {msg}",
                "log",
                True,
                msg,
                False,
                no_update,
                no_update,
                no_update,
                "",
            )

        log_lines: list = []

        def _log(m):
            log_lines.append(m)

        from pycsamt.app.web.utils_pinn import (
            build_hybrid_inv,
            build_pinn_inv,
            decode_npz_checkpoint,
            pinn_stats_div,
            plot_pinn_convergence,
            plot_pinn_section,
            session_to_obs_1d,
            session_to_obs_2d,
        )

        # ── PINN track ────────────────────────────────

        if track == "pinn":
            _log("=" * 50)
            _log("  pycsamt PINN Inversion")
            _log("=" * 50)

            dim = (pinn_dim or "1d").lower()
            data_src = (pinn_src or "session").lower()
            solver = (pinn_solver or "mt1d").lower()
            mode = (pinn_mode or "te").lower()
            n_layers = int(pinn_n_layers or 5)
            depth_max = float(pinn_depth_max or 2000.0)
            n_freqs = int(pinn_n_freqs or 32)
            epochs = int(pinn_epochs or 300)
            lr = float(pinn_lr or 0.01)
            log_every = int(pinn_log_every or 50)
            device = pinn_device or None
            lam_z = float(pinn_lam_z or 0.01)
            lam_x = float(pinn_lam_x or 0.005)
            lam_g = float(pinn_lam_g or 0.005)
            radius = float(pinn_radius or 5000.0)
            spc = float(pinn_spc or 500.0)
            comp_te = pinn_comp_te or "xy"
            comp_tm = pinn_comp_tm or "yx"

            try:
                if not _TORCH_OK:
                    _log("PyTorch not found.")
                    return (
                        _empty,
                        _empty,
                        _coloured_log(log_lines),
                        no_update,
                        "",
                        "PyTorch required for PINN.",
                        "log",
                        False,
                        "",
                        True,
                        no_update,
                        no_update,
                        no_update,
                    )

                _log(f"dim={dim} src={data_src} epochs={epochs} lr={lr}")

                if data_src == "session":
                    if not session_id:
                        return _err13(
                            _coloured_log(log_lines),
                            "No session. Load EDI files first.",
                        )
                    _log("Loading session obs …")
                    if dim == "1d":
                        obs = session_to_obs_1d(
                            session_id,
                            pinn_stations or None,
                            comp=comp_te,
                        )
                    else:
                        obs = session_to_obs_2d(
                            session_id,
                            pinn_stations or None,
                            comp_te=comp_te,
                            comp_tm=comp_tm,
                        )
                    _log(f"Loaded {len(obs)} stations.")
                else:
                    if dim != "1d":
                        return _err13(
                            _coloured_log(log_lines),
                            "Synthetic data is only"
                            " supported for 1-D PINN."
                            " Use session data for"
                            " 2-D / 3-D.",
                        )
                    from pycsamt.ai.inversion._sites_bridge import (  # noqa: E501
                        SiteObs1D as _SiteObs1D,
                    )
                    from pycsamt.forward.em1d import (
                        MT1DForward,
                    )
                    from pycsamt.forward.synthetic import (
                        LayeredModel,
                    )

                    freqs = np.logspace(-2, 4, n_freqs)
                    rng = np.random.default_rng(42)
                    obs = []
                    for i in range(8):
                        log_r = rng.uniform(0.5, 3.5, n_layers)
                        thk = 10.0 ** rng.uniform(1.5, 3.2, n_layers - 1)
                        m = LayeredModel(
                            resistivity=10.0**log_r,
                            thickness=thk,
                        )
                        resp = MT1DForward(freqs).run(m)
                        obs.append(
                            _SiteObs1D(
                                name=f"SYN{i:02d}",
                                freq=freqs,
                                rho_obs=resp.rho_a,
                                phase_obs=resp.phase,
                            )
                        )
                    _log(f"Generated {len(obs)} synthetic 1-D obs.")

                _log("Building PINN inverter …")
                t0 = time.time()
                inv = build_pinn_inv(
                    obs,
                    dim,
                    solver=solver,
                    mode=mode,
                    n_layers=n_layers,
                    depth_max=depth_max,
                    n_freqs=n_freqs,
                    epochs=epochs,
                    lr=lr,
                    comp=comp_te,
                    comp_te=comp_te,
                    comp_tm=comp_tm,
                    smoothness_weight=lam_z,
                    lateral_weight=lam_x,
                    graph_weight=lam_g,
                    radius=radius,
                    station_spacing=spc,
                    device=device,
                )
                _log("Fitting …")
                if dim == "1d":
                    inv.fit(
                        epochs=epochs,
                        verbose=False,
                        log_every=log_every,
                    )
                else:
                    inv.fit(
                        verbose=False,
                        log_every=log_every,
                    )
                elapsed = time.time() - t0
                _log(f"Done in {elapsed:.1f}s.")

                src_res = plot_pinn_section(inv, dim, dark=dark)
                src_conv = plot_pinn_convergence(inv, dark=dark, label="PINN")
                stats = pinn_stats_div(
                    inv,
                    dim,
                    elapsed_s=elapsed,
                    label="PINN",
                )

                try:
                    resid_json = inv.residuals().to_dict("records")
                    sta_names = list(inv.stations)
                except Exception:
                    resid_json = []
                    sta_names = []

                pinn_store = {
                    "residuals": resid_json,
                    "stations": sta_names,
                    "dim": dim,
                }
                sta_opts = [{"label": s, "value": s} for s in sta_names]

                msg = (
                    f"PINN {dim.upper()} complete"
                    f" — {len(obs)} stations"
                    f" in {elapsed:.1f}s"
                )
                return _ok13(
                    src_res,
                    src_conv,
                    _coloured_log(log_lines),
                    stats,
                    msg,
                    pinn_store=pinn_store,
                    sta_opts=sta_opts,
                )

            except Exception as exc:
                import traceback as _tb

                _log(f"ERROR  {exc}")
                _log(_tb.format_exc()[-1200:])
                return _err13(_coloured_log(log_lines), str(exc))

        # ── Hybrid track ──────────────────────────────

        _log("=" * 50)
        _log("  pycsamt Hybrid Inversion")
        _log("=" * 50)

        dim = (hyb_dim or "1d").lower()
        data_src = (hyb_src or "session").lower()
        mode = (hyb_mode or "te").lower()
        n_layers = int(hyb_n_layers or 5)
        n_freqs = int(hyb_n_freqs or 32)
        epochs = int(hyb_epochs or 150)
        lr = float(hyb_lr or 0.005)
        device = hyb_device or None
        lam_z = float(hyb_lam_z or 0.005)
        comp_te = hyb_comp_te or "xy"
        comp_tm = hyb_comp_tm or "yx"

        try:
            if not _TORCH_OK:
                _log("PyTorch not found.")
                return (
                    _empty,
                    _empty,
                    _coloured_log(log_lines),
                    no_update,
                    "",
                    "PyTorch required for Hybrid.",
                    "log",
                    False,
                    "",
                    True,
                    no_update,
                    no_update,
                    no_update,
                )

            if not hyb_ckpt:
                return _err13(
                    _coloured_log(log_lines),
                    "Upload a .npz AI checkpoint before running Hybrid.",
                )

            _log("Decoding AI checkpoint …")
            ckpt_path = decode_npz_checkpoint(hyb_ckpt)

            _log("Loading session obs …")
            if dim == "1d":
                obs = session_to_obs_1d(
                    session_id,
                    hyb_stations or None,
                    comp=comp_te,
                )
            else:
                obs = session_to_obs_2d(
                    session_id,
                    hyb_stations or None,
                    comp_te=comp_te,
                    comp_tm=comp_tm,
                )
            _log(f"Loaded {len(obs)} stations.")

            _log("Building Hybrid inverter …")
            t0 = time.time()
            inv = build_hybrid_inv(
                obs,
                dim,
                ckpt_path,
                n_layers=n_layers if n_layers else None,
                n_freqs=n_freqs,
                mode=mode,
                smoothness_weight=lam_z,
                max_iter=epochs,
                lr=lr,
                comp=comp_te,
                comp_te=comp_te,
                comp_tm=comp_tm,
                device=device,
            )

            import os as _os

            _os.unlink(ckpt_path)

            _log("Fitting (Stage 1 AI + Stage 2 PINN)…")
            inv.fit(verbose=False, log_every=20)
            elapsed = time.time() - t0
            _log(f"Done in {elapsed:.1f}s.")

            src_res = plot_pinn_section(inv, dim, dark=dark)
            src_conv = plot_pinn_convergence(inv, dark=dark, label="Hybrid")
            stats = pinn_stats_div(
                inv,
                dim,
                elapsed_s=elapsed,
                label="Hybrid",
            )

            try:
                resid_json = inv.residuals().to_dict("records")
                sta_names = list(inv.stations)
            except Exception:
                resid_json = []
                sta_names = []

            hyb_store = {
                "residuals": resid_json,
                "stations": sta_names,
                "dim": dim,
            }
            sta_opts = [{"label": s, "value": s} for s in sta_names]
            msg = (
                f"Hybrid {dim.upper()} complete"
                f" — {len(obs)} stations"
                f" in {elapsed:.1f}s"
            )
            return _ok13(
                src_res,
                src_conv,
                _coloured_log(log_lines),
                stats,
                msg,
                hyb_store=hyb_store,
                sta_opts=sta_opts,
            )

        except Exception as exc:
            import traceback as _tb

            _log(f"ERROR  {exc}")
            _log(_tb.format_exc()[-1200:])
            return _err13(_coloured_log(log_lines), str(exc))

    # ---

    @app.callback(
        Output(IDs.IMG_INV_DATA_FIT, "src"),
        Input(IDs.INV_DATAFIT_STA, "value"),
        State(IDs.INV_PINN_RESULT_STORE, "data"),
        State(IDs.INV_HYB_RESULT_STORE, "data"),
        State(IDs.STORE_THEME, "data"),
        prevent_initial_call=True,
    )
    def update_data_fit(station, pinn_store, hyb_store, theme):
        if not station:
            raise PreventUpdate

        dark = (theme or "dark") == "dark"

        store = pinn_store or hyb_store
        if not store:
            raise PreventUpdate

        resid_records = store.get("residuals", [])
        if not resid_records:
            raise PreventUpdate

        import pandas as pd

        from pycsamt.app.web.utils_pinn import (
            plot_pinn_data_fit,
        )

        df = pd.DataFrame(resid_records)
        return plot_pinn_data_fit(
            None,
            station,
            dark=dark,
            df=df,
        )
