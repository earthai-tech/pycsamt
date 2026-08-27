# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Pcsf3DWindow — the desktop app's first 3-D/volume panel.

Phase 7 of the PCSF format plan: renders any ``.pcsf`` inversion-result
file directly (Occam2D/DUHI ``grid2d``, ModEM ``grid3d``, MARE2DEM
``mesh_unstructured``, or a Phase-5 ``multiline`` stack) with no
backend-specific logic in this window — every geometry kind is
flattened the same way by
:func:`pycsamt.format.pointcloud.pcsf_to_point_cloud`, so this window
only ever has to know how to scatter-plot points.

Deliberately independent of the main window's loaded EDI session (no
``set_sites``): the whole point of a persisted ``.pcsf`` file is that a
view can be rebuilt from it alone.

Left params panel
──────────────────
  File          [Load .pcsf file…] + loaded-file summary
  Display       Max points, colormap, [Render]

Right content
─────────────
  MplCanvas — a Matplotlib 3-D (mplot3d) scatter plot, reusing the same
  canvas widget every other desktop panel already uses (no new
  rendering stack — pyqtgraph/OpenGL has no existing precedent
  anywhere in this app and is not introduced here).
"""

from __future__ import annotations

from pathlib import Path

from PySide6.QtWidgets import (
    QComboBox,
    QFileDialog,
    QLabel,
    QSizePolicy,
    QSpinBox,
    QVBoxLayout,
    QWidget,
)

from pycsamt.app.desktop.widgets.mpl_canvas import MplCanvas
from pycsamt.app.desktop.windows._base import (
    PanelWindow,
    icon_button,
    make_group,
)

_CMAPS = ["viridis", "plasma", "turbo", "jet", "RdBu_r", "terrain"]


class Pcsf3DWindow(PanelWindow):
    """Floating window that loads and renders a ``.pcsf`` file in 3-D."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(
            title="PCSF 3D Viewer",
            session_key="pcsf3d",
            params_width=260,
            icon_name="3d",
            parent=parent,
        )
        self.resize(1100, 800)
        self._model = None
        self._cloud = None
        self._path: Path | None = None

    # ── Params panel ──────────────────────────────────────────────────

    def _build_params(self, layout: QVBoxLayout) -> None:
        grp_file, lay_file = make_group("File")
        self._btn_load = icon_button(
            "📂  Load .pcsf file…", "3d", "Open a .pcsf inversion-result file"
        )
        self._btn_load.clicked.connect(self._on_load)
        lay_file.addWidget(self._btn_load)

        self._file_lbl = QLabel("No file loaded.")
        self._file_lbl.setWordWrap(True)
        self._file_lbl.setObjectName("InfoLabel")
        lay_file.addWidget(self._file_lbl)
        layout.addWidget(grp_file)

        grp_disp, lay_disp = make_group("Display")
        lay_disp.addWidget(QLabel("Max points"))
        self._spin_max_points = QSpinBox()
        self._spin_max_points.setRange(1_000, 500_000)
        self._spin_max_points.setSingleStep(5_000)
        self._spin_max_points.setValue(50_000)
        self._spin_max_points.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed
        )
        lay_disp.addWidget(self._spin_max_points)

        lay_disp.addWidget(QLabel("Colormap"))
        self._combo_cmap = QComboBox()
        self._combo_cmap.addItems(_CMAPS)
        self._combo_cmap.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed
        )
        lay_disp.addWidget(self._combo_cmap)

        self._btn_render = icon_button(
            "↻  Render", "3d", "Render the loaded model"
        )
        self._btn_render.clicked.connect(self._on_render)
        self._btn_render.setEnabled(False)
        lay_disp.addWidget(self._btn_render)
        layout.addWidget(grp_disp)

        self._status_lbl = QLabel("")
        self._status_lbl.setObjectName("InfoLabel")
        self._status_lbl.setWordWrap(True)
        layout.addWidget(self._status_lbl)

    # ── Content panel ─────────────────────────────────────────────────

    def _build_content(self, layout: QVBoxLayout) -> None:
        self._canvas = MplCanvas(self, toolbar=True)
        layout.addWidget(self._canvas)

    # ── Public API ────────────────────────────────────────────────────

    def load_pcsf(self, path: str) -> None:
        """Load and render *path* — the programmatic entry point behind
        the "Load .pcsf file…" button, also usable directly (e.g. by
        tests or a future "open with" integration)."""
        from pycsamt.format import read_pcsf

        self._model = read_pcsf(path)
        self._path = Path(path)
        n_lines = (
            len(self._model.geometry.lines)
            if self._model.kind == "multiline"
            else 1
        )
        self._file_lbl.setText(
            f"{self._path.name}\n"
            f"kind={self._model.kind}  backend={self._model.source_backend}"
            + (f"  lines={n_lines}" if self._model.kind == "multiline" else "")
        )
        self._btn_render.setEnabled(True)
        self._on_render()

    # ── Slots ─────────────────────────────────────────────────────────

    def _on_load(self) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self, "Load PCSF file", "", "PCSF files (*.pcsf)"
        )
        if not path:
            return
        try:
            self.load_pcsf(path)
        except Exception as exc:
            self._status_lbl.setText(f"Load error: {exc}")
            self._btn_render.setEnabled(False)

    def _on_render(self) -> None:
        if self._model is None:
            return
        try:
            from pycsamt.format.pointcloud import pcsf_to_point_cloud

            cloud = pcsf_to_point_cloud(
                self._model, max_points=int(self._spin_max_points.value())
            )
            self._cloud = cloud
            self._draw_cloud(cloud)
            self._status_lbl.setText(cloud.label)
        except Exception as exc:
            self._status_lbl.setText(f"Render error: {exc}")

    # ── Rendering ─────────────────────────────────────────────────────

    def _draw_cloud(self, cloud) -> None:
        # Registers the "3d" projection with add_subplot below.
        from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

        fig = self._canvas.figure
        fig.clear()
        ax = fig.add_subplot(111, projection="3d")
        scatter = ax.scatter(
            cloud.x,
            cloud.y,
            cloud.z,
            c=cloud.value,
            cmap=self._combo_cmap.currentText(),
            s=4,
            depthshade=True,
        )
        ax.set_xlabel("x (m)")
        ax.set_ylabel("y / line offset (m)")
        ax.set_zlabel("z (m)")
        fig.colorbar(scatter, ax=ax, label="log10 ρ (Ω·m)", shrink=0.6)
        self._canvas.draw()
