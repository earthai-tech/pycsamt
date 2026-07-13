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
        }

    # ── 2-D  (EMInverter2D) ───────────────────────────────────────────────────

    def _run_2d(self) -> dict:
        from pycsamt.ai.inversion.inv2d import EMInverter2D
        from pycsamt.forward.batch import generate_dataset

        p = self._params
        n_components = int(p.get("n_components", 4))
        n_depth = int(p.get("n_depth", 40))
        n_stations = int(p.get("n_stations", 20))
        n_freqs = int(p.get("n_freq", 32))
        epochs = int(p.get("epochs", 40))
        batch_size = int(p.get("batch_size", 16))
        lr = float(p.get("lr", 1e-3))
        n_samples = int(p.get("n_samples", 500))

        f_min = max(float(p.get("f_min", 1e-3)), 1e-6)
        f_max = max(float(p.get("f_max", 1e2)), f_min * 10)
        freqs = np.logspace(np.log10(f_min), np.log10(f_max), n_freqs)

        self._log(f"Generating {n_samples} 2-D training samples…")
        self.progress.emit(5)

        # For 2D, generate per-station 1D datasets and stack into profiles
        dataset = generate_dataset(
            solver="mt1d",
            n_samples=n_samples * n_stations,
            freqs=freqs,
            n_layers=(3, 6),
            noise_level=0.05,
            verbose=False,
        )
        # Reshape to (n_samples, n_components=2, n_freqs, n_stations)
        # Using only rho_a and phase (2 components) for simplicity
        n_total = n_samples * n_stations
        X_flat = dataset.X[:n_total]  # (n*ns, 2*nf)
        y_flat = dataset.y[:n_total]  # (n*ns, 2*nl-1)
        n_used = (len(X_flat) // n_stations) * n_stations
        X_2d = X_flat[:n_used].reshape(-1, n_stations, 2, n_freqs)
        # → (n_samples, n_stations, 2, n_freqs) → (n_samples, 2, n_freqs, n_stations)
        X_2d = X_2d.transpose(0, 2, 3, 1)
        y_2d = y_flat[:n_used].reshape(-1, n_stations, y_flat.shape[-1])

        self._log(f"Training EMInverter2D  epochs={epochs}…")
        self.progress.emit(20)
        inv = EMInverter2D(
            n_components=n_components,
            n_depth=n_depth,
            n_stations=n_stations,
            n_freqs=n_freqs,
        )
        inv.fit(
            X_2d,
            y_2d,
            epochs=epochs,
            batch_size=batch_size,
            lr=lr,
            verbose=False,
        )
        self.progress.emit(80)

        X_obs = p.get("X_obs")
        if X_obs is not None:
            X_obs = np.asarray(X_obs, dtype=float)
        else:
            X_obs = X_2d[:2]
        self._log(f"Predicting on {len(X_obs)} profile(s)…")
        y_pred = inv.predict(X_obs, as_log_rho=True)

        self.progress.emit(100)
        self._log("AI 2-D inversion complete.")
        return {
            "dim": "2D",
            "y_pred": y_pred,
            "X_obs": X_obs,
            "n_stations": n_stations,
            "n_depth": n_depth,
            "freqs": freqs,
            "inverter": inv,
        }

    # ── 3-D  (GCNInverter3D) ──────────────────────────────────────────────────

    def _run_3d(self) -> dict:
        from pycsamt.ai.inversion.inv3d import GCNInverter3D
        from pycsamt.forward.batch import generate_dataset

        p = self._params
        n_features = int(p.get("n_features", 40))
        n_layers = int(p.get("n_layers", 5))
        hidden = p.get("hidden", [256, 128, 64])
        dropout = float(p.get("dropout", 0.1))
        epochs = int(p.get("epochs", 40))
        batch_size = int(p.get("batch_size", 16))
        lr = float(p.get("lr", 1e-3))
        n_samples = int(p.get("n_samples", 300))
        n_sta = int(p.get("n_sta", 16))
        radius = float(p.get("radius", 5000.0))

        f_min = max(float(p.get("f_min", 1e-3)), 1e-6)
        f_max = max(float(p.get("f_max", 1e1)), f_min * 10)
        n_freq = int(p.get("n_freq", 20))
        freqs = np.logspace(np.log10(f_min), np.log10(f_max), n_freq)

        self._log(
            f"Generating {n_samples} 3-D training samples ({n_sta} stations each)…"
        )
        self.progress.emit(5)

        dataset = generate_dataset(
            solver="mt1d",
            n_samples=n_samples * n_sta,
            freqs=freqs,
            n_layers=(3, 6),
            noise_level=0.05,
            verbose=False,
        )
        n_total = n_samples * n_sta
        X_flat = dataset.X[:n_total]
        y_flat = dataset.y[:n_total]
        n_used = (len(X_flat) // n_sta) * n_sta
        X_3d = X_flat[:n_used].reshape(
            -1, n_sta, X_flat.shape[-1]
        )  # (N, ns, nf*2)
        X_3d = X_3d.reshape(
            -1, X_flat.shape[-1]
        )  # flatten back to (N*ns, nf*2)
        y_3d = y_flat[:n_used]

        # Synthetic station coordinates for training
        grid_side = int(np.ceil(np.sqrt(n_sta)))
        xs = np.linspace(0, (grid_side - 1) * 1000, grid_side)
        ys = np.linspace(0, (grid_side - 1) * 1000, grid_side)
        gx, gy = np.meshgrid(xs, ys)
        coords_base = np.column_stack(
            [gx.ravel()[:n_sta], gy.ravel()[:n_sta]]
        )
        coords = np.tile(coords_base, (n_samples, 1))[:n_used]

        self._log(f"Training GCNInverter3D  epochs={epochs}…")
        self.progress.emit(20)
        inv = GCNInverter3D(
            n_features=n_features,
            n_layers=n_layers,
            hidden=hidden,
            dropout=dropout,
        )
        inv.fit(
            X_3d,
            y_3d,
            coords=coords,
            radius=radius,
            epochs=epochs,
            batch_size=batch_size,
            lr=lr,
            verbose=False,
        )
        self.progress.emit(80)

        X_obs = p.get("X_obs")
        coords_obs = p.get("coords_obs")
        if X_obs is not None:
            X_obs = np.asarray(X_obs, dtype=float)
        else:
            X_obs = X_3d[:n_sta]
            coords_obs = coords[:n_sta]
        self._log(f"Predicting on {len(X_obs)} station(s)…")
        y_pred = inv.predict(
            X_obs, coords=coords_obs, radius=radius, as_log_rho=True
        )

        self.progress.emit(100)
        self._log("AI 3-D inversion complete.")
        return {
            "dim": "3D",
            "y_pred": y_pred,
            "X_obs": X_obs,
            "coords": coords_obs,
            "freqs": freqs,
            "inverter": inv,
        }
