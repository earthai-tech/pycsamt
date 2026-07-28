"""
pycsamt.agents.inversion_prep
==============================

:class:`InversionPrepAgent` — Prepare Occam2D / ModEM data files.

Phase 2 stub: validates the Sites object and writes inversion-ready
data files in the requested format.  Full Occam2D write support
uses pycsamt's existing model writers; ModEM is planned for Phase 3.
"""

from __future__ import annotations

import os
import time
from typing import Any

from ._base import AgentResult, BaseAgent

_SYSTEM_PROMPT = """\
You are an expert MT inversion specialist.
Given a dataset summary and chosen inversion code, recommend:
1. Appropriate period band and error floor.
2. Mesh geometry (cell sizes, padding, depth extent).
3. Regularisation starting values.
4. Any pre-processing steps still needed before inversion.
Reply in 4–5 sentences, plain English.
"""


class InversionPrepAgent(BaseAgent):
    """Prepare MT data files for 2-D / 3-D inversion codes.

    Currently supported output formats
    -----------------------------------
    * ``"occam2d"`` — Occam2D DataFile format
    * ``"modem"``   — ModEM3D data file (Phase 3)

    Input keys
    ----------
    ``sites`` / ``path`` : Sites or str
    ``code`` : str — ``"occam2d"`` (default) or ``"modem"``
    ``period_range`` : [T_min, T_max], optional
    ``component`` : str — ``"xy"``, ``"yx"``, or ``"both"``
    ``error_floor`` : float — minimum relative error floor (default 0.05)
    ``output_dir`` : str

    Output data keys
    ----------------
    ``data_file_path``   str — path to the written data file
    ``n_periods``        int
    ``n_stations``       int
    ``code``             str
    """

    SYSTEM_PROMPT = _SYSTEM_PROMPT

    def __init__(
        self,
        *,
        api_key: str | None = None,
        model: str | None = None,
        llm_provider: str = "claude",
        code: str = "occam2d",
        error_floor: float = 0.05,
    ) -> None:
        super().__init__(
            "InversionPrepAgent",
            api_key=api_key,
            model=model,
            llm_provider=llm_provider,
        )
        self.code = code
        self.error_floor = error_floor

    def execute(self, input_data: dict[str, Any]) -> AgentResult:
        self._last_cost = 0.0
        t0 = time.time()
        warnings: list[str] = []

        from ..emtools._core import (
            _get_z_block,
            _iter_items,
            ensure_sites,
        )

        sites_raw = input_data.get("sites") or input_data.get("path")
        if sites_raw is None:
            return AgentResult.failed("No 'sites' or 'path'.", elapsed=time.time() - t0)
        try:
            sites = ensure_sites(sites_raw, verbose=0)
        except Exception as exc:
            return AgentResult.failed(str(exc), elapsed=time.time() - t0)

        code = str(input_data.get("code", self.code)).lower()
        output_dir = input_data.get("output_dir", "pycsamt_inversion_prep")
        per_range = input_data.get("period_range")
        component = str(input_data.get("component", "both")).lower()
        err_floor = float(input_data.get("error_floor", self.error_floor))
        os.makedirs(output_dir, exist_ok=True)

        # count stations + periods
        n_st = 0
        per_all = []
        import numpy as np

        for _i, ed in enumerate(_iter_items(sites)):
            n_st += 1
            Z_obj, z, fr = _get_z_block(ed)
            if fr is not None:
                per = 1.0 / np.where(fr == 0, np.nan, fr)
                if per_range:
                    mask = (per >= per_range[0]) & (per <= per_range[1])
                    per = per[mask]
                per_all.extend(per[np.isfinite(per)].tolist())

        n_per = len(set(f"{p:.4e}" for p in per_all))

        data_file_path: str | None = None
        if code == "occam2d":
            try:
                from ..models.occam2d import (
                    write_occam2d_data,
                )

                out_path = os.path.join(output_dir, "OccamDataFile.dat")
                write_occam2d_data(
                    sites,
                    out_path,
                    period_range=per_range,
                    component=component,
                    error_floor=err_floor,
                )
                data_file_path = out_path
            except ImportError:
                warnings.append(
                    "pycsamt.models.occam2d not available. "
                    "Data file was not written."
                )
            except Exception as exc:
                warnings.append(f"write_occam2d_data failed: {exc}")
        elif code == "modem":
            warnings.append(
                "ModEM data file writer is planned for Phase 3. " "No file was written."
            )
        else:
            warnings.append(f"Unknown inversion code {code!r}.")

        interp: str | None = None
        if self.api_key:
            prompt = (
                f"Inversion prep summary:\n"
                f"  Code: {code}, stations: {n_st}, periods: {n_per}\n"
                f"  Period range: {per_range}, error floor: {err_floor}\n"
                "Give inversion recommendations."
            )
            interp = self.query_llm(prompt, max_tokens=200)

        elapsed = time.time() - t0
        return AgentResult(
            status="success" if data_file_path else "needs_review",
            summary=(
                f"Inversion prep ({code}): {n_st} stations, ~{n_per} periods. "
                + (
                    f"Data file: {data_file_path}"
                    if data_file_path
                    else "No data file written."
                )
            ),
            data={
                "data_file_path": data_file_path,
                "n_periods": n_per,
                "n_stations": n_st,
                "code": code,
                "sites": sites,
            },
            warnings=warnings,
            llm_interpretation=interp,
            elapsed_seconds=elapsed,
            cost_estimate_usd=self._last_cost,
        )


__all__ = ["InversionPrepAgent"]
