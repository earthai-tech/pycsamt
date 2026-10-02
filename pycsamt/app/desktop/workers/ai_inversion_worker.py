# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
AIInversionWorker — QThread for AI-based 1-D / 2-D / 3-D MT inversion.

Workflow
--------
1. Generate synthetic training dataset   (generate_dataset / generate_dataset_3d)
2. Instantiate and fit the inverter      (EMInverter1D / EMInverter2D / GCNInverter3D)
3. Predict on observed data              (inverter.predict)
4. Emit result as a plain dict ready for plotting

Signals
-------
progress(int)      0-100 estimated progress
log_line(str)      one line of text for the live log
finished(dict)     result dict with keys 'dim', 'model', 'inverter', ...
error(str)         human-readable error on failure
"""

from __future__ import annotations

import logging
from typing import Any

import numpy as np
from PySide6.QtCore import QThread, Signal

logger = logging.getLogger(__name__)


class AIInversionWorker(QThread):
    """Background thread for training + predicting with an AI inverter."""

    progress = Signal(int)
    log_line = Signal(str)
    finished = Signal(dict)
    error = Signal(str)

    def __init__(self, params: dict[str, Any], parent=None) -> None:
        super().__init__(parent)
        self._params = params

    # ── QThread entry point ───────────────────────────────────────────────────

    def run(self) -> None:
        dim = self._params.get("dim", "1D")
        try:
            if dim == "1D":
                result = self._run_1d()
            elif dim == "2D":
                result = self._run_2d()
            else:
                result = self._run_3d()
            self.finished.emit(result)
        except Exception as exc:
            logger.exception("AI inversion worker error")
            self.error.emit(str(exc))

    # ── helpers ───────────────────────────────────────────────────────────────

    def _log(self, msg: str) -> None:
        logger.info(msg)
        self.log_line.emit(msg)

    # ── 1-D  (EMInverter1D) ───────────────────────────────────────────────────

    def _run_1d(self) -> dict:
        import numpy as np

        from pycsamt.ai.inversion.inv1d import EMInverter1D
        from pycsamt.forward.batch import generate_dataset

        p = self._params
        arch = p.get("arch", "resnet")
        n_layers = int(p.get("n_layers", 5))
        solver = p.get("solver", "mt1d")
        epochs = int(p.get("epochs", 60))
        batch_size = int(p.get("batch_size", 256))
        lr = float(p.get("lr", 1e-3))
        n_samples = int(p.get("n_samples", 2000))
        n_freq = int(p.get("n_freq", 30))
        f_min = max(float(p.get("f_min", 1e-3)), 1e-6)
        f_max = max(float(p.get("f_max", 1e3)), f_min * 10)
        noise_lvl = float(p.get("noise_level", 0.05))
        geology = p.get("geology", None) or None

        freqs = np.logspace(np.log10(f_min), np.log10(f_max), n_freq)

        self._log(f"Generating {n_samples} 1-D training samples…")
        self.progress.emit(5)
        dataset = generate_dataset(
            solver=solver,
            n_samples=n_samples,
            freqs=freqs,
            n_layers=(max(2, n_layers - 1), n_layers + 1),
            noise_level=noise_lvl,
            geology=geology,
            verbose=False,
        )
        X_train = dataset.X  # (n_samples, 2*n_freq)
        y_train = dataset.y  # (n_samples, 2*n_layers-1)

        self._log(f"Training EMInverter1D  arch={arch}  epochs={epochs}…")
        self.progress.emit(20)
        inv = EMInverter1D(arch=arch, n_layers=n_layers, solver=solver)
        inv.fit(
            X_train,
            y_train,
            epochs=epochs,
            batch_size=batch_size,
            lr=lr,
            verbose=False,
        )
        self.progress.emit(80)

        # Predict on observed data if supplied, otherwise use last training sample
        X_obs = p.get("X_obs")
        if X_obs is not None:
            X_obs = np.asarray(X_obs, dtype=float)
        else:
            X_obs = X_train[:5]  # demo: first 5 synthetic samples
        self._log(f"Predicting on {len(X_obs)} station(s)…")
        y_pred = inv.predict(X_obs, as_log_rho=True)  # (n, 2*n_layers-1)

        self.progress.emit(100)
        self._log("AI 1-D inversion complete.")

        return {
            "dim": "1D",
            "y_pred": y_pred,
            "X_obs": X_obs,
            "n_layers": n_layers,
            "freqs": freqs,
            "inverter": inv,
            # Real station names when X_obs came from the loaded survey;
            # empty for the synthetic demo samples.
            "stations": list(p.get("stations") or []) if p.get("X_obs")
            is not None else [],
        }

    # ── 2-D  (EMInverter2D via Inv2DAgent) ──────────────────────────────────────
    #
    # Delegates to pycsamt.agents.inv2d_agent.Inv2DAgent rather than building
    # its own dataset/training pipeline. That agent already implements a
    # validated "physics" dispatch (mt1d tiled-1-D vs. mt2d real 2-D
    # finite-difference training data via generate_2d_maxwell_dataset), trains
    # and predicts directly against the real Sites survey, and produces its
    # own section figure -- the hand-rolled version this replaced duplicated
    # that logic with a real, confirmed bug (constructing EMInverter2D with
    # n_components=4 while only ever building a 2-channel array, guaranteeing
    # a shape-mismatch crash on first real use).

    def _run_2d(self) -> dict:
        from pycsamt.agents.inv2d_agent import Inv2DAgent

        p = self._params
        sites = p.get("sites")
        if sites is None:
            raise RuntimeError("No sites loaded for 2-D AI inversion.")
        physics = p.get("physics", "mt1d")
        n_freqs = int(p.get("n_freq", 32))
        f_min = max(float(p.get("f_min", 1e-3)), 1e-6)
        f_max = max(float(p.get("f_max", 1e2)), f_min * 10)
        freqs = np.logspace(np.log10(f_min), np.log10(f_max), n_freqs)

        self._log(f"Running Inv2DAgent (physics={physics})…")
        self.progress.emit(10)
        agent = Inv2DAgent(
            n_depth=int(p.get("n_depth", 40)),
            n_freqs=n_freqs,
            n_components=int(p.get("n_components", 2)),
            n_train_profiles=int(p.get("n_samples", 500)),
            n_stations_per_profile=int(p.get("n_stations", 20)),
            epochs=int(p.get("epochs", 40)),
            physics=physics,
            verbose=False,
        )
        result = agent.execute({"sites": sites, "freqs": freqs})
        self.progress.emit(90)

        for w in result.warnings:
            self._log(f"WARNING: {w}")
        if result.status != "success":
            raise RuntimeError(result.error or "Inv2DAgent failed.")

        self.progress.emit(100)
        self._log(result.summary)
        return {"dim": "2D", "agent_result": result}

    # ── 3-D  (GCNInverter3D via Inv3DAgent) ─────────────────────────────────────
    #
    # Same rationale as _run_2d: delegates to the already-validated
    # pycsamt.agents.inv3d_agent.Inv3DAgent, which trains/predicts on the
    # real Sites survey (including its actual station coordinates for the
    # GCN adjacency graph) instead of a synthetic duplicated-coordinate
    # station grid, and offers the same real mt1d/mt3d physics dispatch.

    def _run_3d(self) -> dict:
        from pycsamt.agents.inv3d_agent import Inv3DAgent

        p = self._params
        sites = p.get("sites")
        if sites is None:
            raise RuntimeError("No sites loaded for 3-D AI inversion.")
        physics = p.get("physics", "mt1d")
        n_freq = int(p.get("n_freq", 20))
        f_min = max(float(p.get("f_min", 1e-3)), 1e-6)
        f_max = max(float(p.get("f_max", 1e1)), f_min * 10)
        freqs = np.logspace(np.log10(f_min), np.log10(f_max), n_freq)

        self._log(f"Running Inv3DAgent (physics={physics})…")
        self.progress.emit(10)
        agent = Inv3DAgent(
            n_layers=int(p.get("n_layers", 5)),
            n_freqs=n_freq,
            hidden=tuple(p.get("hidden", [256, 128, 64])),
            dropout=float(p.get("dropout", 0.1)),
            n_train_profiles=int(p.get("n_samples", 300)),
            epochs=int(p.get("epochs", 40)),
            radius=float(p.get("radius", 5000.0)),
            physics=physics,
            n_mc=0,  # MC-dropout uncertainty not surfaced by this panel yet
            verbose=False,
        )
        result = agent.execute({"sites": sites, "freqs": freqs})
        self.progress.emit(90)

        for w in result.warnings:
            self._log(f"WARNING: {w}")
        if result.status != "success":
            raise RuntimeError(result.error or "Inv3DAgent failed.")

        self.progress.emit(100)
        self._log(result.summary)
        return {"dim": "3D", "agent_result": result}
