# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt.agents.distortion_agent
=================================

:class:`DistortionClassificationAgent` — Triage galvanic-distortion
regime per station.

Wraps :class:`~pycsamt.ai.processing.distortion.DistortionTypeClassifier`:
an MLP self-trained on rule-based labels (clean / static-shift-only /
distorted) that routes each station toward the correction tool that
actually solves its physics --
:func:`~pycsamt.emtools.ss.correct_ss_ama` for static-shift-only
stations, :func:`~pycsamt.emtools.gb.apply_groom_bailey` for fully
distorted ones -- rather than re-solving the distortion physics itself.

**Known, documented instability (not a bug fixed by this agent):**
:meth:`DistortionTypeClassifier.fit`'s own docstring records that on a
small survey (tens of stations, typical for a real MT line), repeated
fits at low epoch counts can disagree on which regime is the *majority*
class, not just on a minority label appearing or not. This agent
defaults to ``epochs=200`` (the class's own low default of 80 is a
"quick look" setting, not a converged one) and always attaches an
explicit warning to the result -- read one run's regime table as one
plausible smoothing of the rule boundary, not a converged ground truth.
"""

from __future__ import annotations

import time
from typing import Any

from ._base import AgentResult, BaseAgent

_SELF_TRAINING_WARNING = (
    "DistortionTypeClassifier is self-trained on rule-based labels. On a "
    "small survey (tens of stations), repeated fits can disagree on which "
    "regime is the majority class, not just on a minority label appearing "
    "or not -- this is a documented instability of the estimator, not a "
    "one-off bug. Treat this run's regime_label column as one plausible "
    "smoothing of the rule boundary, not a converged ground truth; a "
    "different seed or epoch count may relabel some stations."
)

_SYSTEM_PROMPT = """\
You are an expert in MT/AMT galvanic-distortion diagnosis (static shift vs.
full Groom-Bailey distortion).
Given a distortion-triage result summary, write 3-4 sentences that:
1. State the regime counts (clean / static-shift-only / distorted) and
   which correction tool each non-clean group was routed to.
2. Explicitly flag that this classifier's labels come from self-training on
   rule-based thresholds and can shift the majority-class assignment between
   runs on a small survey -- treat the regime table as a plausible triage,
   not a final answer, and recommend a human sanity-check on borderline
   stations before committing to a correction.
3. Note the held-out/training convergence (val_loss) if available.
4. Recommend a concrete next step (e.g. re-run with more epochs, inspect
   feature-space plot for borderline stations).
Reply in plain English.
"""


class DistortionClassificationAgent(BaseAgent):
    """Triage each station's galvanic-distortion regime and route it.

    Parameters
    ----------
    api_key, model, llm_provider : str
    epochs : int
        Training epochs (default 200 -- higher than
        ``DistortionTypeClassifier``'s own low "quick look" default of
        80; see this module's docstring on the self-training
        instability this guards against).
    route_corrections : bool
        When ``True`` (default), also run
        ``correct_ss_ama``/``apply_groom_bailey`` on the routed
        stations and return the corrected sites.

    Input keys
    ----------
    ``sites`` / ``path`` : Sites or str
    ``output_dir`` : str, optional

    Output data keys
    ----------------
    ``table``               pandas DataFrame — station, regime_label, confidence
    ``regime_counts``       dict {regime_label: count}
    ``ss_corrected_sites``  Sites or None — static_shift_only stations, corrected
    ``distorted_corrected_sites`` Sites or None — distorted stations, corrected
    ``figures``             dict
    ``figure_paths``        dict
    """

    SYSTEM_PROMPT = _SYSTEM_PROMPT

    def __init__(
        self,
        *,
        api_key: str | None = None,
        model: str | None = None,
        llm_provider: str = "claude",
        epochs: int = 200,
        route_corrections: bool = True,
    ) -> None:
        super().__init__(
            "DistortionClassificationAgent",
            api_key=api_key,
            model=model,
            llm_provider=llm_provider,
            section_preset="pseudosection",
        )
        self.epochs = epochs
        self.route_corrections = route_corrections

    def execute(self, input_data: dict[str, Any]) -> AgentResult:
        self._last_cost = 0.0
        t0 = time.time()
        warnings: list[str] = [_SELF_TRAINING_WARNING]

        try:
            from ..ai.processing import (
                DistortionTypeClassifier,
                build_distortion_features_table,
            )
            from ..backends import get_backend_instance

            if get_backend_instance() is None:
                raise ImportError("No DL backend.")
        except ImportError as exc:
            return AgentResult.failed(
                f"DistortionClassificationAgent requires PyTorch or TensorFlow: {exc}",
                hint="pip install torch  or  pip install tensorflow",
                elapsed=time.time() - t0,
            )

        from ..emtools._core import ensure_sites

        sites_raw = input_data.get("sites") or input_data.get("path")
        if sites_raw is None:
            return AgentResult.failed(
                "No 'sites' or 'path'.", elapsed=time.time() - t0
            )
        try:
            sites = ensure_sites(sites_raw, verbose=0)
        except Exception as exc:
            return AgentResult.failed(str(exc), elapsed=time.time() - t0)

        output_dir = input_data.get("output_dir")

        try:
            feats = build_distortion_features_table(sites)
        except Exception as exc:
            return AgentResult.failed(
                f"build_distortion_features_table failed: {exc}",
                elapsed=time.time() - t0,
            )
        if len(feats) < 5:
            return AgentResult.failed(
                "Need >= 5 stations to classify distortion regime.",
                elapsed=time.time() - t0,
            )

        try:
            clf = DistortionTypeClassifier.from_features_table(
                feats, epochs=self.epochs, seed=0, verbose=False
            )
            table = clf.predict_table(sites)
        except Exception as exc:
            return AgentResult.failed(
                f"DistortionTypeClassifier fit/predict failed: {exc}",
                elapsed=time.time() - t0,
            )

        regime_counts = table["regime_label"].value_counts().to_dict()

        # ── route to the real correction tools ──────────────────────────────
        ss_corrected_sites = None
        distorted_corrected_sites = None
        if self.route_corrections:
            try:
                from ..emtools.gb import apply_groom_bailey
                from ..emtools.ss import correct_ss_ama

                ss_only = table.loc[
                    table["regime_label"] == "static_shift_only", "station"
                ].tolist()
                distorted = table.loc[
                    table["regime_label"] == "distorted", "station"
                ].tolist()
                if ss_only:
                    ss_corrected_sites = correct_ss_ama(
                        sites.select(names=ss_only), inplace=False
                    )
                if distorted:
                    distorted_corrected_sites = apply_groom_bailey(
                        sites.select(names=distorted), inplace=False
                    )
            except Exception as exc:
                warnings.append(f"Routed correction failed: {exc}")

        # ── figures ───────────────────────────────────────────────────────
        figures: dict[str, Any] = {}
        fig_paths: dict[str, str] = {}
        try:
            from ..ai.processing import plot_distortion_summary

            fig = plot_distortion_summary(
                table,
                history=clf,
                suptitle="DistortionTypeClassifier -- regime triage",
            )
            figures["distortion_summary"] = fig
            p = self._save_figure(
                fig, output_dir, "distortion_summary", warnings_list=warnings
            )
            if p:
                fig_paths["distortion_summary"] = p
        except Exception as exc:
            warnings.append(f"Distortion summary figure: {exc}")

        # ── LLM interpretation ──────────────────────────────────────────────
        interp: str | None = None
        if self.llm_available:
            prompt = (
                f"Distortion triage summary:\n"
                f"  Stations: {len(table)}\n"
                f"  Regime counts: {regime_counts}\n"
                f"  Routed to correct_ss_ama: "
                f"{len(ss_corrected_sites) if ss_corrected_sites is not None else 0}\n"
                f"  Routed to apply_groom_bailey: "
                f"{len(distorted_corrected_sites) if distorted_corrected_sites is not None else 0}\n\n"
                f"IMPORTANT CAVEAT: {_SELF_TRAINING_WARNING}\n\n"
                "Diagnose the survey's distortion landscape and recommend follow-up."
            )
            interp = self.query_llm(prompt, max_tokens=240)

        elapsed = time.time() - t0
        return AgentResult(
            status="success",
            summary=(
                f"Triaged {len(table)} station(s): {regime_counts}. "
                f"{len(figures)} figure(s). Self-training instability "
                "disclaimer attached (see warnings)."
            ),
            data={
                "table": table,
                "regime_counts": regime_counts,
                "ss_corrected_sites": ss_corrected_sites,
                "distorted_corrected_sites": distorted_corrected_sites,
                "figures": figures,
                "figure_paths": fig_paths,
            },
            warnings=warnings,
            llm_interpretation=interp,
            elapsed_seconds=elapsed,
            cost_estimate_usd=self._last_cost,
        )


__all__ = ["DistortionClassificationAgent"]
