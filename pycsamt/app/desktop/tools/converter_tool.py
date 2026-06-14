# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
FormatConverterDialog — re-export a survey folder to EDI or CSV.

Supported output formats
------------------------
EDI   — write one .edi file per station using the loaded survey's EDI writer
CSV   — station list with name, lat, lon, n_freq, T_min, T_max, has_z_err
JSON  — same metadata as CSV but structured as a JSON array

Workflow
--------
1. Source: the already-loaded survey  OR  pick a folder of EDI/AVG/J files.
2. Choose output format and target folder.
3. Click Convert → background worker iterates stations, progress bar updates.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from PySide6.QtCore import Qt, QThread, Signal
from PySide6.QtWidgets import (
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QFileDialog,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QProgressBar,
    QPushButton,
    QTextEdit,
    QVBoxLayout,
    QWidget,
)

_FORMATS = ["CSV (station metadata)", "JSON (station metadata)", "EDI (re-export)"]


# ── Background worker ─────────────────────────────────────────────────────────

class _ConvertWorker(QThread):
    progress = Signal(int, int, str)   # current, total, message
    done     = Signal(str)             # summary
    error    = Signal(str)

    def __init__(self, sites, out_dir: Path, fmt: str):
        super().__init__()
        self._sites   = sites
        self._out_dir = out_dir
        self._fmt     = fmt

    def run(self):
        try:
            from pycsamt.emtools._core import _iter_items, _unwrap
            items = list(_iter_items(self._sites))
        except Exception:
            try:
                items = list(self._sites) if hasattr(self._sites, "__iter__") else [self._sites]
            except Exception as exc:
                self.error.emit(f"Cannot iterate survey: {exc}")
                return

        self._out_dir.mkdir(parents=True, exist_ok=True)
        total = len(items)
        rows  = []

        for idx, ed in enumerate(items):
            try:
                ed = _unwrap(ed)
            except Exception:
                pass

            # Collect metadata
            name = (
                getattr(ed, "station", None)
                or getattr(ed, "id", f"S{idx:03d}")
            )
            lat  = getattr(ed, "lat", None) or getattr(ed, "latitude",  None)
            lon  = getattr(ed, "lon", None) or getattr(ed, "longitude", None)
            try:
                lat = float(lat) if lat is not None else float("nan")
                lon = float(lon) if lon is not None else float("nan")
            except Exception:
                lat = lon = float("nan")

            Z_obj = getattr(ed, "Z", None) or getattr(
                getattr(ed, "edi", None), "Z", None
            )
            n_freq = 0
            t_min = t_max = float("nan")
            has_err = False
            if Z_obj is not None:
                freqs = getattr(Z_obj, "freq", None)
                if freqs is not None and len(freqs) > 0:
                    fa = np.asarray(freqs, dtype=float)
                    fa = fa[fa > 0]
                    n_freq = len(fa)
                    if n_freq:
                        periods = 1.0 / fa
                        t_min, t_max = float(periods.min()), float(periods.max())
                z_err = getattr(Z_obj, "z_err", None)
                if z_err is not None:
                    has_err = bool(np.any(np.isfinite(np.asarray(z_err))))

            row = dict(
                station=str(name),
                lat=lat, lon=lon,
                n_freq=n_freq,
                t_min=round(t_min, 6) if np.isfinite(t_min) else None,
                t_max=round(t_max, 6) if np.isfinite(t_max) else None,
                has_z_err=has_err,
            )
            rows.append(row)

            # EDI re-export (best-effort via mtpy-style write if available)
            if "EDI" in self._fmt:
                try:
                    edi_obj = getattr(ed, "edi", ed)
                    write_fn = getattr(edi_obj, "write_edi_file", None)
                    if write_fn:
                        edi_path = self._out_dir / f"{name}.edi"
                        write_fn(str(edi_path))
                except Exception:
                    pass

            self.progress.emit(idx + 1, total, str(name))

        # Write tabular output
        try:
            if "CSV" in self._fmt:
                import csv
                out_path = self._out_dir / "survey_stations.csv"
                with open(out_path, "w", newline="", encoding="utf-8") as fh:
                    writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
                    writer.writeheader()
                    writer.writerows(rows)
            elif "JSON" in self._fmt:
                import json
                out_path = self._out_dir / "survey_stations.json"
                out_path.write_text(json.dumps(rows, indent=2))
        except Exception as exc:
            self.error.emit(f"Write error: {exc}")
            return

        self.done.emit(
            f"Converted {total} stations → {self._out_dir}"
        )


# ── Dialog ────────────────────────────────────────────────────────────────────

class FormatConverterDialog(QDialog):
    """
    Survey format converter dialog.

    Parameters
    ----------
    sites : any
        Already-loaded survey (from MainWindow._controller.sites).
        If None the dialog still opens and the user can pick a source folder
        (conversion will use the loaded survey once data is loaded).
    parent : QWidget, optional
    """

    def __init__(self, sites, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setWindowTitle("Format Converter")
        self.setMinimumSize(520, 420)
        self._sites  = sites
        self._worker = None
        self._build_ui()

    # ── Build ─────────────────────────────────────────────────────────────────

    def _build_ui(self) -> None:
        root = QVBoxLayout(self)
        root.setSpacing(10)

        # ── Source ────────────────────────────────────────────────────────────
        grp_src = QGroupBox("Source")
        form_src = QFormLayout(grp_src)
        n = getattr(self._sites, "__len__", None)
        src_lbl = "Loaded survey"
        if self._sites is not None:
            try:
                from pycsamt.emtools._core import _iter_items
                cnt = sum(1 for _ in _iter_items(self._sites))
                src_lbl = f"Loaded survey  ({cnt} stations)"
            except Exception:
                pass
        self._src_lbl = QLabel(src_lbl if self._sites else "No survey loaded")
        form_src.addRow("Input:", self._src_lbl)
        root.addWidget(grp_src)

        # ── Output ────────────────────────────────────────────────────────────
        grp_out = QGroupBox("Output")
        form_out = QFormLayout(grp_out)
        form_out.setSpacing(8)

        self._fmt_combo = QComboBox()
        self._fmt_combo.addItems(_FORMATS)
        form_out.addRow("Format:", self._fmt_combo)

        dir_row = QHBoxLayout()
        self._dir_edit = QLineEdit(str(Path.home() / "pycsamt_export"))
        btn_browse = QPushButton("Browse…")
        btn_browse.setFixedWidth(70)
        btn_browse.clicked.connect(self._pick_dir)
        dir_row.addWidget(self._dir_edit)
        dir_row.addWidget(btn_browse)
        form_out.addRow("Output folder:", dir_row)

        root.addWidget(grp_out)

        # ── Progress ──────────────────────────────────────────────────────────
        self._progress = QProgressBar()
        self._progress.setValue(0)
        root.addWidget(self._progress)

        self._log = QTextEdit()
        self._log.setReadOnly(True)
        self._log.setFixedHeight(100)
        root.addWidget(self._log)

        # ── Buttons ───────────────────────────────────────────────────────────
        btn_row = QHBoxLayout()
        self._run_btn = QPushButton("Convert")
        self._run_btn.clicked.connect(self._on_run)
        btn_row.addWidget(self._run_btn)
        btn_row.addStretch()
        box = QDialogButtonBox(QDialogButtonBox.StandardButton.Close)
        box.rejected.connect(self.reject)
        btn_row.addWidget(box)
        root.addLayout(btn_row)

    # ── Slots ─────────────────────────────────────────────────────────────────

    def _pick_dir(self) -> None:
        d = QFileDialog.getExistingDirectory(
            self, "Select output folder", self._dir_edit.text()
        )
        if d:
            self._dir_edit.setText(d)

    def _on_run(self) -> None:
        if self._sites is None:
            self._log.append("Load survey data first.")
            return
        out_dir = Path(self._dir_edit.text().strip())
        fmt     = self._fmt_combo.currentText()
        self._run_btn.setEnabled(False)
        self._progress.setValue(0)
        self._log.clear()
        self._log.append(f"Converting → {out_dir}  [{fmt}]")

        self._worker = _ConvertWorker(self._sites, out_dir, fmt)
        self._worker.progress.connect(self._on_progress)
        self._worker.done.connect(self._on_done)
        self._worker.error.connect(self._on_error)
        self._worker.start()

    def _on_progress(self, cur: int, total: int, name: str) -> None:
        self._progress.setMaximum(total)
        self._progress.setValue(cur)
        self._log.append(f"  [{cur}/{total}]  {name}")

    def _on_done(self, msg: str) -> None:
        self._run_btn.setEnabled(True)
        self._log.append(f"\nDone — {msg}")

    def _on_error(self, msg: str) -> None:
        self._run_btn.setEnabled(True)
        self._log.append(f"\nError: {msg}")
