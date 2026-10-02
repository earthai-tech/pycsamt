# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
AI inversion page of the Inversion Studio.

Three learned inverters, trained on synthetic physics then applied to the
loaded survey:

* **EMInverter1D** — ResNet / CNN / FCN, one layered model per station;
* **EMInverter2D** — U-Net section through ``Inv2DAgent``;
* **GCNInverter3D** — graph network through ``Inv3DAgent``.

Results: model, data fit, training loss, and — for 1-D — a comparison
with a classical Occam1D run loaded in *View Results* (same stations).
"""

from __future__ import annotations

import numpy as np
from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import (
    QButtonGroup,
    QComboBox,
    QDoubleSpinBox,
    QFormLayout,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QPushButton,
    QRadioButton,
    QScrollArea,
    QSizePolicy,
    QSpinBox,
    QSplitter,
    QStackedWidget,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

from pycsamt.app.desktop.controllers.inversion_engines import (
    PlotUnavailable,
    site_list,
    sites_to_1d_features,
)
from pycsamt.app.desktop.widgets.canvas_stack import CanvasResultView
from pycsamt.app.desktop.windows._base import make_group
from pycsamt.app.desktop.windows.inversion.monitor import RunMonitor

AI_MODELS = [
    ("inv1d", "1D", "EMInverter1D", "ResNet / CNN / FCN — one layered model "
     "per station, in seconds after training."),
    ("inv2d", "2D", "EMInverter2D", "U-Net resistivity section along the "
     "profile (Inv2DAgent)."),
    ("inv3d", "3D", "GCNInverter3D", "Graph neural network over the station "
     "network (Inv3DAgent)."),
]


def _spin(value, lo, hi, decimals=0, step=1.0) -> QDoubleSpinBox:
    sb = QDoubleSpinBox()
    sb.setRange(lo, hi)
    sb.setDecimals(decimals)
    sb.setSingleStep(step)
    sb.setValue(value)
    sb.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
    return sb


def _ispin(value, lo, hi) -> QSpinBox:
    sb = QSpinBox()
    sb.setRange(lo, hi)
    sb.setValue(value)
    sb.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
    return sb


def draw_into(view: CanvasResultView, fn) -> bool:
    """Run ``fn(figure)`` on *view*'s canvas; show a reason card on
    :class:`PlotUnavailable` or failure.  ``fn`` may return a new Figure."""
    canvas = view.canvas
    canvas.figure.clear()
    try:
        out = fn(canvas.figure)
    except PlotUnavailable as exc:
        view.show_unavailable(exc.title, exc.reason, exc.guidance)
        return False
    except Exception as exc:  # a library plot failed on this data
        view.show_unavailable("Plot failed", f"{type(exc).__name__}: {exc}")
        return False
    if out is not None and out is not canvas.figure:
        canvas.show_figure(out)
    else:
        canvas.draw()
    view.show_canvas()
    return True


class AIInversionPage(QWidget):
    """Train & predict with an AI inverter; show its results."""

    log = Signal(str)
    run_state = Signal(str, int, bool)  # text, percent (-1 busy), running
    result_ready = Signal(dict)

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._sites = None
        self._worker = None
        self._result: dict | None = None
        self._classical = None  # LoadedRun (Occam1D) for "Compare"

        split = QSplitter(Qt.Orientation.Horizontal, self)
        split.setChildrenCollapsible(False)
        lay = QVBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.addWidget(split)

        # ── left: model + settings ────────────────────────────────────
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setHorizontalScrollBarPolicy(
            Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        scroll.setMinimumWidth(270)
        scroll.setMaximumWidth(380)
        left = QWidget()
        lv = QVBoxLayout(left)
        lv.setContentsMargins(6, 6, 6, 6)
        lv.setSpacing(6)
        grp, gl = make_group("AI model")
        self._model_group = QButtonGroup(self)
        self._model_rbs: dict[str, QRadioButton] = {}
        for key, dim, label, desc in AI_MODELS:
            rb = QRadioButton(f"{label}  ·  {dim}")
            rb.setProperty("engine", key)
            rb.setToolTip(desc)
            self._model_group.addButton(rb)
            self._model_rbs[key] = rb
            gl.addWidget(rb)
            d = QLabel(desc)
            d.setObjectName("InfoLabel")
            d.setWordWrap(True)
            d.setContentsMargins(22, 0, 0, 4)
            gl.addWidget(d)
        self._rb_inv1d = self._model_rbs["inv1d"]
        self._rb_inv2d = self._model_rbs["inv2d"]
        self._rb_inv3d = self._model_rbs["inv3d"]
        self._rb_inv1d.setChecked(True)
        lv.addWidget(grp)

        grp, gl = make_group("Training && network")
        self._ai_stack = QStackedWidget()
        self._ai_stack.addWidget(self._build_cfg_inv1d())
        self._ai_stack.addWidget(self._build_cfg_inv2d())
        self._ai_stack.addWidget(self._build_cfg_inv3d())
        gl.addWidget(self._ai_stack)
        lv.addWidget(grp)
        self._model_group.buttonToggled.connect(self._on_model_changed)

        grp, gl = make_group("Data")
        f = QFormLayout()
        f.setSpacing(4)
        self._f_min = _spin(1e-3, 1e-5, 1e3, 5, 1e-4)
        self._f_max = _spin(1e3, 1e-3, 1e6, 1, 10.0)
        self._f_min.setSuffix(" Hz")
        self._f_max.setSuffix(" Hz")
        f.addRow("f min:", self._f_min)
        f.addRow("f max:", self._f_max)
        gl.addLayout(f)
        self._data_lbl = QLabel("No data loaded")
        self._data_lbl.setObjectName("InfoLabel")
        self._data_lbl.setWordWrap(True)
        gl.addWidget(self._data_lbl)
        lv.addWidget(grp)

        row = QHBoxLayout()
        self._btn_run = QPushButton("▶  Train && predict")
        self._btn_run.setObjectName("RunButton")
        self._btn_run.clicked.connect(self.run)
        row.addWidget(self._btn_run)
        lv.addLayout(row)
        lv.addStretch(1)
        scroll.setWidget(left)
        split.addWidget(scroll)

        # ── right: monitor + results ──────────────────────────────────
        right = QWidget()
        rv = QVBoxLayout(right)
        rv.setContentsMargins(0, 0, 0, 0)
        rv.setSpacing(4)
        self.monitor = RunMonitor()
        self.monitor.stop_requested.connect(self.stop)
        rv.addWidget(self.monitor)
        self._tabs = QTabWidget()
        self._tab_model_view = CanvasResultView(
            empty_title="No prediction yet",
            empty_reason="Train an AI inverter to see its resistivity model.")
        self._tab_fit_view = CanvasResultView(
            empty_title="No data yet",
            empty_reason="The observed soundings appear after a run.")
        self._tab_convergence_view = CanvasResultView(
            empty_title="No training yet",
            empty_reason="The training loss appears after a run.")
        self._tab_compare_view = CanvasResultView(
            empty_title="Nothing to compare yet",
            empty_reason="Run EMInverter1D on the loaded stations and open "
                         "an Occam1D run in Classical ▸ View Results.")
        self._tab_model = self._tab_model_view.canvas
        self._tab_fit = self._tab_fit_view.canvas
        self._tab_convergence = self._tab_convergence_view.canvas
        self._tabs.addTab(self._tab_model_view, "Model")
        self._tabs.addTab(self._tab_fit_view, "Data")
        self._tabs.addTab(self._tab_convergence_view, "Training loss")
        self._tabs.addTab(self._tab_compare_view, "Compare with classical")
        rv.addWidget(self._tabs, 1)
        split.addWidget(right)
        split.setStretchFactor(1, 1)
        split.setSizes([300, 900])

    # ── configs (attribute names kept from the pre-2.6 window) ───────────
    def _build_cfg_inv1d(self) -> QWidget:
        w = QWidget()
        f = QFormLayout(w)
        f.setContentsMargins(0, 0, 0, 0)
        f.setSpacing(4)
        self._ai1_arch = QComboBox()
        self._ai1_arch.addItems(["resnet", "cnn1d", "fcn"])
        self._ai1_solver = QComboBox()
        self._ai1_solver.addItems(["mt1d", "csamt1d", "tem1d"])
        self._ai1_n_layers = _ispin(5, 2, 20)
        self._ai1_n_samples = _ispin(2000, 100, 50000)
        self._ai1_epochs = _ispin(60, 5, 500)
        self._ai1_batch = _ispin(256, 8, 2048)
        self._ai1_lr = _spin(1e-3, 1e-5, 1e-1, 5, 1e-4)
        self._ai1_noise = _spin(0.05, 0.0, 0.5, 3, 0.01)
        self._ai1_geology = QComboBox()
        self._ai1_geology.addItems([
            "(none)", "sedimentary", "crystalline", "geothermal", "marine",
            "permafrost", "basement", "coastal", "evaporite", "hydrothermal",
            "laterite", "mineralized", "porphyry", "volcanic"])
        f.addRow("Architecture:", self._ai1_arch)
        f.addRow("Physics:", self._ai1_solver)
        f.addRow("Layers:", self._ai1_n_layers)
        f.addRow("Training samples:", self._ai1_n_samples)
        f.addRow("Epochs:", self._ai1_epochs)
        f.addRow("Batch size:", self._ai1_batch)
        f.addRow("Learning rate:", self._ai1_lr)
        f.addRow("Noise level:", self._ai1_noise)
        f.addRow("Geology prior:", self._ai1_geology)
        return w

    def _build_cfg_inv2d(self) -> QWidget:
        # Inv2DAgent fixes batch size and learning rate internally, so no
        # widgets are offered for them.
        w = QWidget()
        f = QFormLayout(w)
        f.setContentsMargins(0, 0, 0, 0)
        f.setSpacing(4)
        self._ai2_physics = QComboBox()
        self._ai2_physics.addItems(["mt1d (fast, tiled 1-D)",
                                    "mt2d (2-D finite-difference)"])
        self._ai2_n_comp = _ispin(2, 1, 8)
        self._ai2_n_depth = _ispin(40, 10, 200)
        self._ai2_n_sta = _ispin(20, 2, 100)
        self._ai2_n_freq = _ispin(32, 8, 128)
        self._ai2_n_samples = _ispin(500, 50, 10000)
        self._ai2_epochs = _ispin(40, 5, 300)
        f.addRow("Physics:", self._ai2_physics)
        f.addRow("Components:", self._ai2_n_comp)
        f.addRow("Depth cells:", self._ai2_n_depth)
        f.addRow("Stations/profile:", self._ai2_n_sta)
        f.addRow("Frequencies:", self._ai2_n_freq)
        f.addRow("Training profiles:", self._ai2_n_samples)
        f.addRow("Epochs:", self._ai2_epochs)
        return w

    def _build_cfg_inv3d(self) -> QWidget:
        # Inv3DAgent derives n_features from n_freqs and takes the station
        # graph from the real survey coordinates.
        w = QWidget()
        f = QFormLayout(w)
        f.setContentsMargins(0, 0, 0, 0)
        f.setSpacing(4)
        self._ai3_physics = QComboBox()
        self._ai3_physics.addItems(["mt1d (fast, tiled 1-D)",
                                    "mt3d (3-D finite-difference, small "
                                    "grid)"])
        self._ai3_n_layers = _ispin(5, 2, 20)
        self._ai3_hidden = QLineEdit("256,128,64")
        self._ai3_hidden.setToolTip("Hidden layer sizes, comma-separated")
        self._ai3_dropout = _spin(0.1, 0.0, 0.9, 2, 0.05)
        self._ai3_n_samples = _ispin(300, 50, 5000)
        self._ai3_epochs = _ispin(40, 5, 300)
        self._ai3_radius = _spin(5000.0, 100.0, 1e6, 0, 500.0)
        self._ai3_radius.setSuffix(" m")
        f.addRow("Physics:", self._ai3_physics)
        f.addRow("Layers:", self._ai3_n_layers)
        f.addRow("Hidden layers:", self._ai3_hidden)
        f.addRow("Dropout:", self._ai3_dropout)
        f.addRow("Training profiles:", self._ai3_n_samples)
        f.addRow("Epochs:", self._ai3_epochs)
        f.addRow("Graph radius:", self._ai3_radius)
        return w

    # ── state ─────────────────────────────────────────────────────────
    def current_engine(self) -> str:
        btn = self._model_group.checkedButton()
        return btn.property("engine") if btn else "inv1d"

    def _on_model_changed(self, *_):
        idx = {"inv1d": 0, "inv2d": 1, "inv3d": 2}[self.current_engine()]
        self._ai_stack.setCurrentIndex(idx)

    def set_sites(self, sites, selected: list[str] | None = None) -> None:
        self._sites = sites
        self._selected = list(selected or [])
        n = len(site_list(sites))
        self._data_lbl.setText(
            f"{n} station(s) loaded — predictions use the stations ticked in "
            "Classical ▸ Data." if n else "No data loaded — the 1-D model "
            "then predicts five synthetic samples as a demo.")

    def _selected_sites(self):
        if self._sites is None:
            return None
        names = getattr(self, "_selected", [])
        if not names:
            return self._sites
        try:
            return self._sites.select(names=names)
        except Exception:
            return self._sites

    def set_classical_run(self, run) -> None:
        """An Occam1D LoadedRun to compare 1-D predictions with."""
        self._classical = run
        self._draw_compare()

    @property
    def running(self) -> bool:
        return self._worker is not None and self._worker.isRunning()

    # ── params ────────────────────────────────────────────────────────
    def _build_ai_params(self, engine: str, dim: str | None = None) -> dict:
        dim = dim or {"inv1d": "1D", "inv2d": "2D", "inv3d": "3D"}[engine]
        p: dict = {"dim": dim, "f_min": self._f_min.value(),
                   "f_max": self._f_max.value()}
        if engine == "inv1d":
            geo = self._ai1_geology.currentText()
            p.update({
                "arch": self._ai1_arch.currentText(),
                "solver": self._ai1_solver.currentText(),
                "n_layers": self._ai1_n_layers.value(),
                "n_samples": self._ai1_n_samples.value(),
                "epochs": self._ai1_epochs.value(),
                "batch_size": self._ai1_batch.value(),
                "lr": self._ai1_lr.value(),
                "noise_level": self._ai1_noise.value(),
                "geology": None if geo == "(none)" else geo,
                "n_freq": 30,
            })
            X, names = self._build_X_obs(engine, dim, p)
            p["X_obs"] = X
            p["stations"] = names
            return p
        if engine == "inv2d":
            p.update({
                "physics": "mt2d" if self._ai2_physics.currentIndex() == 1
                else "mt1d",
                "n_components": self._ai2_n_comp.value(),
                "n_depth": self._ai2_n_depth.value(),
                "n_stations": self._ai2_n_sta.value(),
                "n_freq": self._ai2_n_freq.value(),
                "n_samples": self._ai2_n_samples.value(),
                "epochs": self._ai2_epochs.value(),
                "sites": self._selected_sites(),
            })
            return p
        try:
            hidden = [int(x) for x in self._ai3_hidden.text().split(",")
                      if x.strip()]
        except ValueError:
            hidden = [256, 128, 64]
        p.update({
            "physics": "mt3d" if self._ai3_physics.currentIndex() == 1
            else "mt1d",
            "n_layers": self._ai3_n_layers.value(),
            "hidden": hidden or [256, 128, 64],
            "dropout": self._ai3_dropout.value(),
            "n_samples": self._ai3_n_samples.value(),
            "epochs": self._ai3_epochs.value(),
            "radius": self._ai3_radius.value(),
            "sites": self._selected_sites(),
        })
        return p

    def _build_X_obs(self, engine: str, dim: str, params: dict):
        """Real soundings on the training frequency grid, or (None, [])."""
        sites = self._selected_sites()
        if sites is None:
            return None, []
        f_min = max(float(params.get("f_min", 1e-3)), 1e-6)
        f_max = max(float(params.get("f_max", 1e3)), f_min * 10)
        freqs = np.logspace(np.log10(f_min), np.log10(f_max),
                            int(params.get("n_freq", 30)))
        return sites_to_1d_features(sites, freqs)

    # ── run ───────────────────────────────────────────────────────────
    def run(self) -> None:
        if self.running:
            return
        from pycsamt.app.desktop.workers.ai_inversion_worker import (
            AIInversionWorker,
        )

        engine = self.current_engine()
        if engine != "inv1d" and self._sites is None:
            self.monitor.finish("error", "Load EDI data first — the 2-D and "
                                "3-D inverters predict on the real survey.")
            return
        params = self._build_ai_params(engine)
        label = {k: lab for k, _d, lab, _t in AI_MODELS}[engine]
        self._worker = AIInversionWorker(params, parent=self)
        self._worker.log_line.connect(self.log)
        self._worker.progress.connect(self._on_progress)
        self._worker.finished.connect(self._on_finished)
        self._worker.error.connect(self._on_error)
        self.monitor.start(f"Training {label}", max_iter=0, target=None)
        self._btn_run.setEnabled(False)
        self.log.emit(f"── AI inversion: {label} ──")
        if engine == "inv1d":
            n = len(params.get("stations") or [])
            self.log.emit(f"Predicting on {n} real station(s)." if n else
                          "No stations loaded: predicting 5 synthetic "
                          "samples (demo).")
        self.run_state.emit(f"AI: {label}", 0, True)
        self._worker.start()

    def stop(self) -> None:
        if self._worker is not None and self._worker.isRunning():
            # Training runs inside library code with no cancel hook; the
            # thread is abandoned and its result ignored.
            self._worker.finished.disconnect()
            self._worker.error.disconnect()
            self._worker = None
        self._btn_run.setEnabled(True)
        self.monitor.finish("stopped", "Stopped — the training thread is "
                            "left to finish in the background.")
        self.run_state.emit("", 0, False)

    def _on_progress(self, pct: int) -> None:
        self.monitor.bar.setRange(0, 1000)
        self.monitor.bar.setValue(int(pct) * 10)
        self.run_state.emit("AI inversion", int(pct), True)

    def _on_finished(self, result: dict) -> None:
        self._btn_run.setEnabled(True)
        self._result = result
        self.monitor.finish("done", "AI inversion finished")
        self.run_state.emit("", 100, False)
        self._plot_ai_result(result)
        self.result_ready.emit({"engine": self.current_engine(),
                                "result": result})

    def _on_error(self, msg: str) -> None:
        self._btn_run.setEnabled(True)
        self.monitor.finish("error", msg.splitlines()[0][:160])
        self.log.emit(f"ERROR: {msg}")
        self.run_state.emit("", 0, False)

    # ── plots ─────────────────────────────────────────────────────────
    def _plot_ai_result(self, result: dict) -> None:
        if result.get("dim", "1D") in ("2D", "3D"):
            self._plot_agent_result(result)
        else:
            draw_into(self._tab_model_view,
                      lambda fig: self._draw_model_1d(fig, result))
            draw_into(self._tab_fit_view,
                      lambda fig: self._draw_data_1d(fig, result))
            self._plot_loss(result.get("inverter"))
            self._draw_compare()
        self._tabs.setCurrentIndex(0)

    def _plot_agent_result(self, result: dict) -> None:
        agent = result.get("agent_result")
        figures = {} if agent is None else dict(agent.get("figures") or {})

        def model(fig):
            if not figures:
                raise PlotUnavailable("No figure returned",
                                      "The inversion agent returned no "
                                      "section figure.")
            return next(iter(figures.values()))

        draw_into(self._tab_model_view, model)

        def summary(fig):
            ax = fig.add_subplot(111)
            ax.set_axis_off()
            lines = []
            if agent is not None:
                if agent.get("physics"):
                    lines.append(f"Physics: {agent.get('physics')}")
                if agent.get("rms_global") is not None:
                    lines.append(f"Data-space RMS: {agent.get('rms_global'):.3f}")
                if agent.warnings:
                    lines.append(f"{len(agent.warnings)} warning(s) — see "
                                 "the console")
            ax.text(0.5, 0.5, "\n".join(lines) or "No summary available",
                    transform=ax.transAxes, ha="center", va="center",
                    fontsize=11)
            return fig

        draw_into(self._tab_fit_view, summary)
        self._plot_loss(None if agent is None else agent.get("inverter"))
        self._tab_compare_view.show_unavailable(
            "Comparison is 1-D only",
            "Compare with classical overlays EMInverter1D profiles on "
            "Occam1D models of the same stations.")

    def _plot_loss(self, inv) -> None:
        def draw(fig):
            hist = getattr(inv, "loss_history_", None)
            if inv is None or hist is None or not len(hist):
                raise PlotUnavailable(
                    "Loss history not available",
                    "This inverter did not record a training-loss history.")
            ax = fig.add_subplot(111)
            ax.plot(np.arange(1, len(hist) + 1), hist, "-", color="#1864ab",
                    lw=1.6)
            ax.set_yscale("log")
            ax.set_xlabel("Epoch")
            ax.set_ylabel("Loss")
            ax.set_title("Training loss", fontsize=10)
            ax.grid(True, alpha=0.3)
            return fig

        draw_into(self._tab_convergence_view, draw)

    @staticmethod
    def _profiles(result: dict):
        y = result.get("y_pred")
        n_layers = int(result.get("n_layers", 5))
        if y is None:
            return []
        out = []
        for row in np.asarray(y, float):
            rho = 10.0 ** row[:n_layers]
            thick = np.abs(row[n_layers:])
            tops = np.concatenate([[0.0], np.cumsum(thick)])
            out.append((tops, rho))
        return out

    def _draw_model_1d(self, fig, result: dict):
        profiles = self._profiles(result)
        if not profiles:
            raise PlotUnavailable("No model prediction returned",
                                  "The inversion did not return a "
                                  "resistivity prediction.")
        names = result.get("stations") or [f"sample {i + 1}"
                                            for i in range(len(profiles))]
        ax = fig.add_subplot(111)
        for (tops, rho), name in list(zip(profiles, names))[:12]:
            depth = np.append(tops, tops[-1] * 1.5 if tops[-1] else 1.0)
            ax.step(np.append(rho, rho[-1]), np.maximum(depth, 1.0),
                    where="post", lw=1.4, alpha=0.85, label=name)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.invert_yaxis()
        ax.set_xlabel("Resistivity (Ω·m)")
        ax.set_ylabel("Depth (m)")
        ax.set_title("EMInverter1D predicted models"
                     + ("" if result.get("stations") else " (synthetic demo)"),
                     fontsize=10)
        ax.grid(True, which="both", alpha=0.25)
        ax.legend(fontsize=7, ncol=2)
        return fig

    def _draw_data_1d(self, fig, result: dict):
        X, freqs = result.get("X_obs"), result.get("freqs")
        if X is None or freqs is None:
            raise PlotUnavailable("No observed data",
                                  "The run returned no soundings.")
        n = len(freqs)
        period = 1.0 / np.asarray(freqs, float)
        names = result.get("stations") or [f"sample {i + 1}"
                                            for i in range(len(X))]
        ax1 = fig.add_subplot(211)
        ax2 = fig.add_subplot(212, sharex=ax1)
        for row, name in list(zip(np.asarray(X, float), names))[:8]:
            ax1.loglog(period, 10.0 ** row[:n], "o-", ms=2.5, lw=1,
                       label=name)
            ax2.semilogx(period, row[n:2 * n], "o-", ms=2.5, lw=1)
        ax1.set_ylabel("ρa (Ω·m)")
        ax2.set_ylabel("φ (°)")
        ax2.set_xlabel("Period (s)")
        ax1.set_title("Soundings given to the network", fontsize=10)
        ax1.legend(fontsize=7, ncol=2)
        for ax in (ax1, ax2):
            ax.grid(True, which="both", alpha=0.25)
        return fig

    def _draw_compare(self) -> None:
        def draw(fig):
            result, run = self._result, self._classical
            if not result or result.get("dim", "1D") != "1D":
                raise PlotUnavailable(
                    "No AI 1-D result",
                    "Run EMInverter1D on the loaded stations first.")
            names = result.get("stations") or []
            if not names:
                raise PlotUnavailable(
                    "Synthetic prediction",
                    "The last AI run used demo samples, not real stations.")
            if run is None or run.engine != "occam1d":
                raise PlotUnavailable(
                    "No classical run loaded",
                    "Open an Occam1D run folder in Classical ▸ View Results.")
            from pycsamt.app.desktop.controllers.inversion_engines import (
                _o1d_load_station,
            )

            dirs = run.extra.get("dirs", {})
            common = [n for n in names if n in dirs][:4]
            if not common:
                raise PlotUnavailable(
                    "No common stations",
                    "The AI prediction and the Occam1D run share no station "
                    "names.")
            profiles = dict(zip(names, self._profiles(result)))
            for i, name in enumerate(common, 1):
                ax = fig.add_subplot(1, len(common), i)
                tops, rho = profiles[name]
                ax.step(np.append(rho, rho[-1]),
                        np.maximum(np.append(tops, tops[-1] * 1.5), 1.0),
                        where="post", color="#c92a2a", lw=1.6, label="AI")
                st = _o1d_load_station(dirs[name])
                if st["result"] is not None:
                    d = np.maximum(st["depth"], 1.0)
                    ax.step(st["result"].final.resistivity, d, where="post",
                            color="#1864ab", lw=1.4, label="Occam1D")
                ax.set_xscale("log")
                ax.set_yscale("log")
                ax.invert_yaxis()
                ax.set_title(name, fontsize=9)
                ax.set_xlabel("ρ (Ω·m)")
                if i == 1:
                    ax.set_ylabel("Depth (m)")
                    ax.legend(fontsize=7)
                ax.grid(True, which="both", alpha=0.25)
            return fig

        draw_into(self._tab_compare_view, draw)


__all__ = ["AIInversionPage", "draw_into"]
