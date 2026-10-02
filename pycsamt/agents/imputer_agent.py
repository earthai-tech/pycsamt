# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt.agents.imputer_agent
=============================

:class:`ImputerAgent` — Fill genuinely missing MT impedance cells.

Wraps :class:`~pycsamt.ai.processing.imputer.EMImputer`: a 1-D
convolutional autoencoder trained by synthetically hiding observed
cells, then applied at inference time to reconstruct cells that are
*actually* missing (``NaN``) on the loaded survey -- distinct from
:class:`~pycsamt.agents.denoising.DenoisingAgent`'s ``EMDenoiser``,
which suppresses noise on cells that are present but corrupted.

Requires PyTorch **or** TensorFlow; degrades to a frequency-axis
linear-interpolation fallback (no cross-station information) when
neither is available -- see :class:`EMImputer`'s own docstring.
"""

from __future__ import annotations

import time
from typing import Any

import numpy as np

from ._base import AgentResult, BaseAgent

_SYSTEM_PROMPT = """\
You are an expert in MT/AMT data completeness and gap-filling.
Given a gap-filling result summary, write 3-4 sentences that:
1. State how many stations/frequency rows were genuinely missing and now filled.
2. Note which impedance components were reconstructed vs. left untouched
   (n_components=4 only models the off-diagonal Zxy/Zyx; the diagonal
   Zxx/Zyy at a missing row stays NaN unless n_components=8 was used).
3. Flag any station with an unusually large gap fraction as needing review.
4. Recommend whether the filled data is safe to use directly in downstream
   processing/inversion, or should be treated as lower-confidence.
Reply in plain English.
"""


class ImputerAgent(BaseAgent):
    """Fill genuinely missing impedance cells in MT data.

    Parameters
    ----------
    api_key, model, llm_provider : str
    n_components : {4, 8}
        4 (default) models the off-diagonal Zxy/Zyx only; 8 also models
        the diagonal Zxx/Zyy.
    mask_frac : float
        Fraction of *observed* cells synthetically hidden on every
        training step to supply a self-supervised reconstruction
        target (default 0.15).
    epochs : int
        Training epochs (default 80).
    channels : tuple of int
        Encoder/decoder channel widths (default ``(32, 64, 32)``).

    Input keys
    ----------
    ``sites`` / ``path`` : Sites or str
    ``output_dir`` : str, optional

    Output data keys
    ----------------
    ``filled_sites``     Sites with genuinely-missing cells reconstructed
    ``gap_table``        pandas DataFrame {station, freq, missing}
    ``n_missing_rows``   int — total genuinely-missing frequency rows
    ``n_gapped_stations`` int — stations with at least one missing row
    ``gapped_stations``  list[str]
    ``figures``          dict
    ``figure_paths``     dict
    """

    SYSTEM_PROMPT = _SYSTEM_PROMPT

    def __init__(
        self,
        *,
        api_key: str | None = None,
        model: str | None = None,
        llm_provider: str = "claude",
        n_components: int = 4,
        mask_frac: float = 0.15,
        epochs: int = 80,
        channels: tuple[int, ...] = (32, 64, 32),
    ) -> None:
        super().__init__(
            "ImputerAgent",
            api_key=api_key,
            model=model,
            llm_provider=llm_provider,
            section_preset="pseudosection",
        )
        # Cast defensively: the registry UI's "n_components" combo (4/8)
        # passes its value through as the widget's raw text (a str), same
        # as every QComboBox-backed registry param -- do not rely on
        # EMImputer's own int(n_components) cast to paper over that.
        self.n_components = int(n_components)
        self.mask_frac = mask_frac
        self.epochs = epochs
        self.channels = channels

    def execute(self, input_data: dict[str, Any]) -> AgentResult:
        self._last_cost = 0.0
        t0 = time.time()
        warnings: list[str] = []

        try:
            from ..ai.processing import EMImputer, prepare_z_features
            from ..backends import get_backend_instance

            if get_backend_instance() is None:
                raise ImportError("No DL backend.")
        except ImportError as exc:
            return AgentResult.failed(
                f"ImputerAgent requires PyTorch or TensorFlow: {exc}",
                hint="pip install torch  or  pip install tensorflow",
                elapsed=time.time() - t0,
            )

        from ..emtools._core import _get_z_block, _iter_items, _name, ensure_sites

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
        labels = [_name(ed, i) for i, ed in enumerate(_iter_items(sites))]

        try:
            X = prepare_z_features(sites, n_components=self.n_components)
        except Exception as exc:
            return AgentResult.failed(
                f"prepare_z_features failed: {exc}", elapsed=time.time() - t0
            )
        n_sites, _n_comp, n_freq = X.shape

        freq = None
        for ed in _iter_items(sites):
            result = _get_z_block(ed, with_errors=False)
            freq = result[2] if len(result) >= 3 else None
            if freq is not None:
                break
        if freq is None:
            return AgentResult.failed(
                "Could not resolve a frequency grid from the loaded sites.",
                elapsed=time.time() - t0,
            )

        missing_mask = ~np.isfinite(X)
        row_missing = missing_mask.any(axis=1)  # (n_sites, n_freq)
        per_site = row_missing.sum(axis=1)
        gapped_stations = [labels[i] for i in range(n_sites) if per_site[i] > 0]
        n_missing_rows = int(row_missing.sum())

        if n_missing_rows == 0:
            return AgentResult(
                status="success",
                summary="No genuinely missing cells found -- nothing to fill.",
                data={
                    "filled_sites": sites,
                    "gap_table": None,
                    "n_missing_rows": 0,
                    "n_gapped_stations": 0,
                    "gapped_stations": [],
                    "figures": {},
                    "figure_paths": {},
                },
                elapsed_seconds=time.time() - t0,
            )

        import pandas as pd

        gap_table = pd.DataFrame(
            {
                "station": np.repeat(labels, n_freq),
                "freq": np.tile(freq, n_sites),
                "missing": row_missing.astype(int).ravel(),
            }
        )

        # ── fit + fill ────────────────────────────────────────────────────
        try:
            imputer = EMImputer(
                n_components=self.n_components, channels=self.channels
            )
            imputer.fit(
                X, mask_frac=self.mask_frac, epochs=self.epochs, verbose=False
            )
        except Exception as exc:
            return AgentResult.failed(
                f"EMImputer.fit failed: {exc}", elapsed=time.time() - t0
            )

        try:
            X_filled = imputer.transform(X)
            filled_sites = imputer.apply(sites, inplace=False)
        except Exception as exc:
            return AgentResult.failed(
                f"EMImputer.transform/apply failed: {exc}",
                elapsed=time.time() - t0,
            )

        # ── figures ───────────────────────────────────────────────────────
        figures: dict[str, Any] = {}
        fig_paths: dict[str, str] = {}
        try:
            from ..ai.processing import plot_imputer_summary

            fig = plot_imputer_summary(
                gap_table,
                freq,
                X,
                X_filled,
                history=imputer,
                station_labels=labels,
                n_show=2,
                suptitle="EMImputer -- gap-filling summary",
            )
            figures["imputer_summary"] = fig
            p = self._save_figure(
                fig, output_dir, "imputer_summary", warnings_list=warnings
            )
            if p:
                fig_paths["imputer_summary"] = p
        except Exception as exc:
            warnings.append(f"Imputer summary figure: {exc}")

        # ── LLM interpretation ──────────────────────────────────────────────
        interp: str | None = None
        if self.llm_available:
            top5 = sorted(
                zip(labels, per_site.tolist()), key=lambda x: -x[1]
            )[:5]
            prompt = (
                f"Gap-filling summary:\n"
                f"  Stations: {n_sites}, gapped: {len(gapped_stations)}\n"
                f"  Missing frequency rows filled: {n_missing_rows}\n"
                f"  n_components modelled: {self.n_components} "
                f"({'off-diagonal only' if self.n_components == 4 else 'all 4 components'})\n"
                f"  Largest gaps: {top5}\n\n"
                "Assess data completeness and recommend follow-up."
            )
            interp = self.query_llm(prompt, max_tokens=220)

        elapsed = time.time() - t0
        return AgentResult(
            status="success",
            summary=(
                f"Filled {n_missing_rows} missing frequency row(s) across "
                f"{len(gapped_stations)}/{n_sites} station(s). "
                f"{len(figures)} figure(s)."
            ),
            data={
                "filled_sites": filled_sites,
                "gap_table": gap_table,
                "n_missing_rows": n_missing_rows,
                "n_gapped_stations": len(gapped_stations),
                "gapped_stations": gapped_stations,
                "figures": figures,
                "figure_paths": fig_paths,
            },
            warnings=warnings,
            llm_interpretation=interp,
            elapsed_seconds=elapsed,
            cost_estimate_usd=self._last_cost,
        )


__all__ = ["ImputerAgent"]
