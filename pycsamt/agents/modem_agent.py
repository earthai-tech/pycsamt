"""
pycsamt.agents.modem_agent
===========================

:class:`ModEmAgent` — Write a ModEM3D inversion data file.

Wraps :class:`~pycsamt.models.modem.ModEmData`:

* :meth:`~pycsamt.models.modem.ModEmData.from_edi` — build data object
  from any EDI source.
* :meth:`~pycsamt.models.modem.ModEmData.write` — write the ModEM-format
  data file.

The ModEM data file format stores MT impedances in the standardised
ModEM3D block format that the ModEM Fortran binary reads directly.
"""

from __future__ import annotations

import os
import time
from pathlib import Path
from typing import Any

from ._base import AgentResult, BaseAgent

_SYSTEM_PROMPT = """\
You are an expert in 3-D MT inversion setup using ModEM3D.
Given a ModEM data file summary, write 3–4 sentences that:
1. Confirm the data contains the expected stations, periods, and components.
2. Recommend an initial model (background resistivity, layer structure).
3. Suggest suitable covariance parameters (smoothing length, roughness).
4. Flag any data issues that could cause convergence problems.
Reply in plain English.
"""


class ModEmAgent(BaseAgent):
    """Write a ModEM3D MT data file from EDI sources.

    Parameters
    ----------
    api_key, model, llm_provider : str
    component_types : list of str or None
        Impedance components to include.
        E.g. ``["Full_Impedance", "Full_Vertical_Components"]``.
        Default: ``None`` → ModEmData auto-selects.
    error_floor : float
        Minimum relative error floor (default 0.05 = 5 %).

    Input keys
    ----------
    ``sites`` / ``path`` : Sites or str
    ``output_dir`` : str
    ``period_range`` : [T_min, T_max], optional
    ``component_types`` : list, optional
    ``error_floor`` : float, optional

    Output data keys
    ----------------
    ``data_path``    Path — ModEM data file
    ``n_stations``   int
    ``n_periods``    int
    ``output_dir``   str

    Examples
    --------
    >>> agent  = ModEmAgent()
    >>> result = agent.execute({
    ...     "path":       "/data/WILLY_DATA",
    ...     "output_dir": "/out/modem",
    ... })
    >>> print(result["data_path"])
    /out/modem/ModEM_Data.dat
    """

    SYSTEM_PROMPT = _SYSTEM_PROMPT

    def __init__(
        self,
        *,
        api_key: str | None = None,
        model: str | None = None,
        llm_provider: str = "claude",
        component_types: list[str] | None = None,
        error_floor: float = 0.05,
    ) -> None:
        super().__init__(
            "ModEmAgent",
            api_key=api_key,
            model=model,
            llm_provider=llm_provider,
            section_preset="inversion",
        )
        self.component_types = component_types
        self.error_floor     = error_floor

    def execute(self, input_data: dict[str, Any]) -> AgentResult:
        self._last_cost = 0.0
        t0 = time.time()
        warnings: list[str] = []

        # ── resolve sites ─────────────────────────────────────────────────────
        from ..emtools._core import ensure_sites
        sites_raw = input_data.get("sites") or input_data.get("path")
        if sites_raw is None:
            return AgentResult.failed("No 'sites' or 'path'.", elapsed=time.time()-t0)
        try:
            sites = ensure_sites(sites_raw, verbose=0)
        except Exception as exc:
            return AgentResult.failed(str(exc), elapsed=time.time()-t0)

        output_dir      = input_data.get("output_dir", "pycsamt_modem")
        per_range       = input_data.get("period_range")
        comp_types      = input_data.get("component_types", self.component_types)
        error_floor     = float(input_data.get("error_floor", self.error_floor))
        os.makedirs(output_dir, exist_ok=True)

        # ── import ModEM layer ────────────────────────────────────────────────
        try:
            from ..models.modem import ModEmData
            from ..models.modem.config import ModEmConfig
        except ImportError as exc:
            return AgentResult.failed(
                f"pycsamt.models.modem not available: {exc}",
                elapsed=time.time() - t0,
            )

        # ── build ModEmData ───────────────────────────────────────────────────
        data_path: Path | None = None
        n_stations = n_periods = 0

        try:
            cfg_kwargs: dict[str, Any] = {
                "error_floor_z": error_floor,
            }
            if comp_types:
                # accept list or single string; ModEmConfig takes a single string
                ct = comp_types[0] if isinstance(comp_types, list) else str(comp_types)
                cfg_kwargs["component_type"] = ct

            cfg      = ModEmConfig(**cfg_kwargs)
            modem_d  = ModEmData.from_edi(sites, config=cfg)

            n_stations = int(modem_d.n_sites)
            n_periods  = int(modem_d.n_periods)

            # apply period range filter if requested
            if per_range is not None:
                from ..emtools._core import _iter_items
                self._log.info(
                    "Period range filter (%s – %s s) is applied via config.",
                    per_range[0], per_range[1],
                )
                warnings.append(
                    "Period-range filtering should be applied to the Sites "
                    "before calling ModEmAgent for precise control."
                )

            data_path = Path(output_dir) / "ModEM_Data.dat"
            modem_d.write(data_path)
            self._log.info(
                "ModEM_Data.dat written: %d stations × %d periods",
                n_stations, n_periods,
            )
        except Exception as exc:
            return AgentResult.failed(
                f"ModEmData.from_edi failed: {exc}",
                hint=(
                    "Check that all EDIs have valid coordinates (lat/lon) and "
                    "impedance data.  Pass period_range to limit the band."
                ),
                elapsed=time.time() - t0,
            )

        # ── validate ──────────────────────────────────────────────────────────
        if data_path.exists():
            size_kb = data_path.stat().st_size // 1024
            self._log.info("ModEM data file: %s (%d KB)", data_path.name, size_kb)
        else:
            warnings.append(f"ModEM data file not found after write: {data_path}")

        # ── LLM interpretation ────────────────────────────────────────────────
        interp: str | None = None
        if self.api_key and n_stations:
            prompt = (
                f"ModEM3D data file summary:\n"
                f"  Stations: {n_stations}, periods: {n_periods}\n"
                f"  Components: {comp_types or 'auto'}\n"
                f"  Error floor: {error_floor:.0%}\n"
                f"  Warnings: {warnings[:3] if warnings else 'none'}\n\n"
                "Advise on the 3-D inversion setup."
            )
            interp = self.query_llm(prompt, max_tokens=200)

        elapsed = time.time() - t0
        return AgentResult(
            status="success" if data_path and data_path.exists() else "needs_review",
            summary=(
                f"ModEM3D prep: {n_stations} stations × {n_periods} periods. "
                f"Data file: {data_path}."
            ),
            data={
                "data_path":  data_path,
                "n_stations": n_stations,
                "n_periods":  n_periods,
                "output_dir": output_dir,
            },
            warnings=warnings,
            llm_interpretation=interp,
            elapsed_seconds=elapsed,
            cost_estimate_usd=self._last_cost,
        )


__all__ = ["ModEmAgent"]
