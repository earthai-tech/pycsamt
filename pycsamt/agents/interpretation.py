"""
pycsamt.agents.interpretation
==============================

:class:`InterpretationAgent` — Translate a resistivity model into geology.

This agent is primarily LLM-driven: it takes the resistivity model output
(from inversion or forward modelling) plus optional geological context and
uses the LLM to write a formation-level interpretation.  It also maps
standard resistivity ranges to lithology classes (after Keller 1988 /
Palacky 1987) as a rule-based fallback when no LLM key is available.
"""

from __future__ import annotations

import time
from typing import Any

from ._base import AgentResult, BaseAgent

_SYSTEM_PROMPT = """\
You are an expert hydrogeologist and applied geophysicist specialising in \
electrical resistivity interpretation.
Given a 1-D or 2-D resistivity model from MT/AMT inversion, write a \
geological interpretation that:
1. Identifies the likely lithological units (e.g. weathered zone, basement,
   aquifer, clay layer) based on resistivity ranges.
2. Estimates formation depths and thicknesses.
3. Discusses implications for groundwater, mineral exploration, or hazard
   assessment depending on the survey context.
4. States uncertainties and what additional data would reduce them.
Reply in plain scientific English, 5–8 sentences. No bullet points.
"""

# ── Keller (1988) / Palacky (1987) resistivity-lithology table ───────────────
_LITHO_TABLE = [
    (0, 2, "sea water / saline brine"),
    (2, 10, "marine clay / saturated shale"),
    (10, 50, "weathered basement / clay-rich soil"),
    (50, 200, "fractured rock / sandy alluvium"),
    (200, 1000, "limestone / sandstone / crystalline basement"),
    (1000, 5000, "compact limestone / granite / quartzite"),
    (5000, 1e9, "massive resistive basement / dry crystalline rock"),
]


def resistivity_to_lithology(rho: float) -> str:
    """Return a likely lithological description for a resistivity value (Ω·m)."""
    for lo, hi, label in _LITHO_TABLE:
        if lo <= rho < hi:
            return label
    return "unknown lithology"


class InterpretationAgent(BaseAgent):
    """Interpret a resistivity model in terms of geological formations.

    Parameters
    ----------
    api_key, model, llm_provider : str
    context : str
        Optional geological context passed to the LLM
        (e.g. ``"Semi-arid Precambrian terrain, looking for aquifers"``).

    Input keys
    ----------
    ``model`` : dict or LayeredModel
        ``{"resistivities": [...], "thicknesses": [...]}``
    ``rms`` : float, optional
    ``context`` : str, optional — overrides constructor default
    ``sites`` / ``path`` : optional — for period range context
    ``output_dir`` : str, optional

    Output data keys
    ----------------
    ``layer_interpretations``   list of dict
    ``summary_text``            str
    ``dominant_lithology``      str
    ``formation_depths_m``      list[float]
    """

    SYSTEM_PROMPT = _SYSTEM_PROMPT

    def __init__(
        self,
        *,
        api_key: str | None = None,
        model: str | None = None,
        llm_provider: str = "claude",
        context: str = "",
    ) -> None:
        super().__init__(
            "InterpretationAgent",
            api_key=api_key,
            model=model,
            llm_provider=llm_provider,
        )
        self.context = context

    def execute(self, input_data: dict[str, Any]) -> AgentResult:
        self._last_cost = 0.0
        t0 = time.time()
        warnings: list[str] = []

        model_raw = input_data.get("model")
        if model_raw is None:
            model_raw = input_data.get("layered_model")
        geo_context = input_data.get("context", self.context)
        rms = input_data.get("rms")

        if model_raw is None:
            return AgentResult.failed(
                "No model available for "
                "interpretation. Run inversion "
                "first and ensure EDI data is "
                "loaded.",
                elapsed=time.time() - t0,
            )

        # ── extract resistivities and thicknesses ─────────────────────────────
        if isinstance(model_raw, dict):
            rhos = list(
                model_raw.get(
                    "resistivity", model_raw.get("resistivities", [])
                )
            )
            ths = list(
                model_raw.get("thickness", model_raw.get("thicknesses", []))
            )
        else:
            rhos = list(
                getattr(
                    model_raw,
                    "resistivity",
                    getattr(model_raw, "resistivities", []),
                )
            )
            ths = list(
                getattr(
                    model_raw,
                    "thickness",
                    getattr(model_raw, "thicknesses", []),
                )
            )

        if not rhos:
            return AgentResult.failed(
                "Model has no resistivities.",
                elapsed=time.time() - t0,
            )

        # ── rule-based layer interpretation ───────────────────────────────────
        depths: list[float] = [0.0]
        for h in ths:
            depths.append(depths[-1] + float(h))
        # Pad depths if thicknesses shorter than layers
        while len(depths) <= len(rhos):
            depths.append(depths[-1] + 50.0)

        layer_interps: list[dict] = []
        for k, rho in enumerate(rhos):
            top = depths[k]
            bot = depths[k + 1] if k + 1 < len(depths) else None
            litho = resistivity_to_lithology(float(rho))
            layer_interps.append(
                {
                    "layer": k + 1,
                    "resistivity": float(rho),
                    "depth_top_m": top,
                    "depth_bot_m": bot,
                    "thickness_m": float(ths[k]) if k < len(ths) else None,
                    "lithology": litho,
                }
            )

        dominant = max(
            layer_interps, key=lambda d: d.get("thickness_m") or 1e9
        )["lithology"]

        # ── LLM interpretation ────────────────────────────────────────────────
        interp_text = ""
        if self.api_key:

            def _bot_str(d):
                b = d.get("depth_bot_m")
                return "basement" if b is None else f"{b:.0f}"

            layer_str = "\n".join(
                f"  Layer {d['layer']}: "
                f"{d['resistivity']:.0f} Ohm.m, "
                f"depth {d['depth_top_m']:.0f}"
                f"-{_bot_str(d)} m"
                f" -- {d['lithology']}"
                for d in layer_interps
            )
            rms_line = f"RMS misfit: {rms:.3f}\n" if rms else ""
            prompt = (
                f"Geological context: "
                f"{geo_context or 'not specified'}\n"
                f"{rms_line}"
                f"Resistivity model:\n"
                f"{layer_str}\n\n"
                "Write a geological interpretation."
            )
            interp_text = self.query_llm(prompt, max_tokens=350) or ""
        else:
            # fallback: build a simple rule-based summary
            lines = [
                f"Layer {d['layer']}: {d['resistivity']:.0f} Ω·m "
                f"({d['depth_top_m']:.0f} m) — {d['lithology']}"
                for d in layer_interps
            ]
            interp_text = ". ".join(lines) + "."

        elapsed = time.time() - t0
        return AgentResult(
            status="success",
            summary=(
                f"Interpreted {len(layer_interps)}-layer model. "
                f"Dominant lithology: {dominant}."
            ),
            data={
                "layer_interpretations": layer_interps,
                "summary_text": interp_text,
                "dominant_lithology": dominant,
                "formation_depths_m": depths,
            },
            warnings=warnings,
            llm_interpretation=interp_text,
            elapsed_seconds=elapsed,
            cost_estimate_usd=self._last_cost,
        )


__all__ = ["InterpretationAgent"]
