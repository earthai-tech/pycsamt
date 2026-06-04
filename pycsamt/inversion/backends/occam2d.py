# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Occam2D adapter backend."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from ..base import BaseInversionBackend
from ..results import InversionResult

__all__ = ["Occam2DBackend"]


class Occam2DBackend(BaseInversionBackend):
    """Prepare, optionally run, or load an Occam2D inversion directory."""

    name = "occam2d"
    supports = (
        ("mt", "2d"),
        ("amt", "2d"),
        ("csamt", "2d"),
        ("emap", "2d"),
    )

    def run(self, data: Any | None = None) -> InversionResult:
        self.check_supported()
        cfg = self.config
        source = cfg.data if data is None else data
        workdir = Path(cfg.workdir)
        workdir.mkdir(parents=True, exist_ok=True)
        warnings: list[str] = []
        files: dict[str, str] = {}
        native = None
        status = "prepared"

        try:
            from ...models.occam2d import InputBuilder, InversionResult as OccamResult
        except ImportError as exc:
            raise ImportError("Occam2D backend requires pycsamt.models.occam2d.") from exc

        if source is not None:
            builder = InputBuilder(source, workdir=workdir)
            opts = dict(cfg.backend_options)
            builder.build(
                error_floor_rho=opts.pop("error_floor_rho", cfg.error_floor),
                error_floor_phase=opts.pop("error_floor_phase", cfg.phase_error),
                **opts,
            )
            native = builder
            if getattr(builder, "is_ready", False):
                for key, obj in {
                    "data": builder.data,
                    "mesh": builder.mesh,
                    "model": builder.model,
                    "startup": builder.startup,
                }.items():
                    path = getattr(obj, "path", None)
                    if path is not None:
                        files[key] = str(path)

        if cfg.run_external:
            try:
                from ...models.occam2d import OccamRunner

                runner = OccamRunner(workdir)
                runner.run(**cfg.backend_options.get("runner", {}))
                native = runner
            except Exception as exc:
                warnings.append(f"Occam2D runner failed: {exc}")
                status = "needs_review"

        try:
            loaded = OccamResult(workdir)
            if getattr(loaded, "rho_2d", None) is not None:
                native = loaded
                status = "loaded"
                rms = float(getattr(loaded, "final_rms", float("nan")))
            else:
                rms = float("nan")
        except Exception as exc:
            warnings.append(f"Occam2D result loading skipped: {exc}")
            rms = float("nan")

        return InversionResult(
            method=cfg.method,
            dimension=cfg.dimension,
            backend=self.name,
            status=status,
            data=source,
            rms=rms,
            workdir=str(workdir),
            files=files,
            native=native,
            warnings=warnings,
            metadata=cfg.metadata,
        )
