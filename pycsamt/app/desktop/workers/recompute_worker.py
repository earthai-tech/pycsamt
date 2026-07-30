# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
RecomputeWorker — QThread worker that drives EDIRecomputer station-by-station,
emitting Qt signals for per-station progress updates.
"""

from __future__ import annotations

from typing import Any

from PySide6.QtCore import QThread, Signal


class RecomputeWorker(QThread):
    """Background EDI-recompute thread.

    Drives ``EDIRecomputer`` one station at a time so the UI receives
    per-station progress signals rather than a single blocking call.

    Signals
    -------
    station_done(int, int, str, str)
        (done_count, total, station_name, status)  after each station.
    finished(object)
        :class:`pycsamt.site.recompute.EDIRecomputeResult` on success.
    error(str)
        Human-readable error message if setup or an unhandled exception fails.
    """

    station_done = Signal(
        int, int, str, str, str
    )  # done, total, name, status, message
    finished = Signal(object)  # EDIRecomputeResult
    error = Signal(str)

    def __init__(
        self,
        source: Any,
        *,
        output_root: str | None = None,
        output_name: str = "recomputed_edis",
        preserve_line_dirs: bool = True,
        template: str = "{source_stem}.edi",
        overwrite: bool = False,
        write: bool = True,
        manifest_csv: bool = True,
        rotate_angle: float | None = None,
        rotate_components: tuple = ("Z", "Tip"),
        fmin: float | None = None,
        fmax: float | None = None,
        fill_missing_values: str | None = None,
        recompute_resphase: bool = True,
        datatype: str | None = None,
        synthesize_spectra: bool = False,
        recursive: bool = True,
        strict: bool = False,
        on_dup: str = "replace",
        verbose: int = 0,
        parent=None,
    ) -> None:
        super().__init__(parent)
        self._source = source
        self._kwargs = dict(
            output_root=output_root,
            output_name=output_name,
            preserve_line_dirs=preserve_line_dirs,
            template=template,
            overwrite=overwrite,
            write=write,
            manifest_csv=manifest_csv,
            rotate_angle=rotate_angle,
            rotate_components=rotate_components,
            fmin=fmin,
            fmax=fmax,
            fill_missing_values=fill_missing_values,
            recompute_resphase=recompute_resphase,
            datatype=datatype,
            synthesize_spectra=synthesize_spectra,
            recursive=recursive,
            strict=strict,
            on_dup=on_dup,
            verbose=verbose,
        )
        self.result = None

    def run(self) -> None:
        try:
            from pycsamt.site.recompute import EDIRecomputer

            kw = self._kwargs.copy()

            def _progress(
                done: int,
                total: int,
                station: str,
                status: str,
                message: str = "",
            ) -> None:
                if self.isInterruptionRequested():
                    raise InterruptedError("Recompute cancelled.")
                self.station_done.emit(done, total, station, status, message)

            recomputer = EDIRecomputer(
                progress=False,
                copy=True,
                progress_callback=_progress,
                **kw,
            )
            result = recomputer.run(self._source)

            self.result = result
            self.finished.emit(result)

        except InterruptedError as exc:
            self.error.emit(str(exc))
        except Exception as exc:
            self.error.emit(str(exc))
