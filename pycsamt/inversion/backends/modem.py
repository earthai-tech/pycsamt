# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""ModEM adapter backend."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from ..base import BaseInversionBackend
from ..results import InversionResult

__all__ = ["ModEMBackend"]


class ModEMBackend(BaseInversionBackend):
    """Prepare or load ModEM 2-D/3-D workflow files."""

    name = "modem"
    supports = (
        ("mt", "2d"),
        ("mt", "3d"),
        ("amt", "2d"),
        ("amt", "3d"),
    )

    def run(self, data: Any | None = None) -> InversionResult:
        self.check_supported()
        cfg = self.config
        source = cfg.data if data is None else data
        workdir = Path(cfg.workdir)
        workdir.mkdir(parents=True, exist_ok=True)
        warnings: list[str] = []
        native = None

        try:
            from ...models import modem
        except ImportError as exc:
            raise ImportError("ModEM backend requires pycsamt.models.modem.") from exc

        builder_cls = getattr(modem, "InputBuilder", None)
        if builder_cls is not None and source is not None:
            try:
                modem_cfg = cfg.backend_options.get("config")
                builder = builder_cls(config=modem_cfg)
                build_options = {
                    key: value for key, value in cfg.backend_options.items()
                    if key != "config"
                }
                files = builder.build(source, workdir=workdir, **build_options)
                native = builder
                file_map = {key: str(value) for key, value in files.items()}
                status = "prepared"
            except Exception as exc:
                warnings.append(f"ModEM preparation failed: {exc}")
                file_map = {}
                status = "needs_review"
        else:
            status = "prepared"
            file_map = {}
            warnings.append("ModEM InputBuilder is not available; created workdir only.")

        result_cls = getattr(modem, "InversionResult", None)
        rms = float("nan")
        if result_cls is not None:
            try:
                loaded = result_cls(workdir)
                native = loaded
                status = "loaded"
                rms = float(getattr(loaded, "final_rms", getattr(loaded, "rms", float("nan"))))
            except Exception as exc:
                warnings.append(f"ModEM result loading skipped: {exc}")

        return InversionResult(
            method=cfg.method,
            dimension=cfg.dimension,
            backend=self.name,
            status=status,
            data=source,
            rms=rms,
            workdir=str(workdir),
            files=file_map,
            native=native,
            warnings=warnings,
            metadata=cfg.metadata,
        )
