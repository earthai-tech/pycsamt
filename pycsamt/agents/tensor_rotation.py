"""
pycsamt.agents.tensor_rotation
================================

:class:`TensorRotationAgent` — Apply a strike rotation to impedance tensors
and write corrected EDI files.

Closes the loop between :class:`~pycsamt.agents.PhaseAnalysisAgent` (which
identifies the optimal strike angle) and the inversion-preparation step:
rotated EDIs are written to disk and are ready to be consumed by Occam2D,
ModEM, or any external inversion code.

Rotation convention
-------------------
The standard two-sided rotation ``Z' = R(θ) Z R(θ)ᵀ`` is applied per
frequency using :func:`~pycsamt.seg.ops.rotate_impedance`.  Tipper vectors
are rotated with :func:`~pycsamt.seg.ops.rotate_tipper`.  A positive angle
θ rotates the measurement frame counter-clockwise (geological azimuth
convention: N → E positive).
"""

from __future__ import annotations

import os
import time
from typing import Any

import numpy as np

from ._base import AgentResult, BaseAgent

_SYSTEM_PROMPT = """\
You are an expert in MT tensor rotation and coordinate-frame correction.
Given a rotation result, write 3-4 sentences that:
1. State the rotation angle applied and the original coordinate frame.
2. Describe how the off-diagonal impedances (Zxy, Zyx) changed after rotation.
3. Assess whether the rotation removed apparent 2-D coupling from the diagonal terms.
4. Recommend whether further refinement (per-frequency rotation, decomposition) is needed.
Reply in plain English.
"""


class TensorRotationAgent(BaseAgent):
    """Rotate impedance tensors and tipper vectors by a fixed strike angle.

    Parameters
    ----------
    api_key, model, llm_provider : str
    strike_deg : float
        Default rotation angle (degrees). Overridden by ``input_data["strike_deg"]``.

    Input keys
    ----------
    ``sites`` / ``path`` : Sites or str
    ``strike_deg`` : float — counter-clockwise rotation angle (degrees)
    ``output_dir`` : str — directory for rotated EDI files
    ``overwrite`` : bool — allow overwriting existing files (default False)
    ``file_suffix`` : str — appended to station name, e.g. ``"_rot"``

    Output data keys
    ----------------
    ``strike_deg``        float — angle applied
    ``written_paths``     list[str] — successfully written EDI paths
    ``failed_stations``   list[str]
    ``n_written``         int
    ``z_diag_reduction``  float — mean |Zxx/Zxy| before − after (proxy for rotation quality)
    ``figures``           dict
    ``figure_paths``      dict

    Examples
    --------
    >>> agent = TensorRotationAgent(strike_deg=42.0)
    >>> r = agent.execute({
    ...     "path":       "/data/WILLY_EDIs",
    ...     "output_dir": "/data/WILLY_rotated",
    ... })
    >>> print(r["n_written"], "EDIs written")
    """

    SYSTEM_PROMPT = _SYSTEM_PROMPT

    def __init__(
        self,
        *,
        api_key: str | None = None,
        model: str | None = None,
        llm_provider: str = "claude",
        strike_deg: float = 0.0,
    ) -> None:
        super().__init__(
            "TensorRotationAgent",
            api_key=api_key,
            model=model,
            llm_provider=llm_provider,
            section_preset="pseudosection",
        )
        self.strike_deg = float(strike_deg)

    def execute(self, input_data: dict[str, Any]) -> AgentResult:
        self._last_cost = 0.0
        t0 = time.time()
        warnings: list[str] = []

        try:
            from ..seg.ops import (
                rotate_impedance,
                rotate_tipper,
            )
        except ImportError as exc:
            return AgentResult.failed(
                f"pycsamt.seg.ops not available: {exc}",
                elapsed=time.time() - t0,
            )

        from ..emtools._core import (
            _get_z_block,
            _iter_items,
            _name,
            ensure_sites,
        )

        sites_raw = input_data.get("sites") or input_data.get("path")
        if sites_raw is None:
            return AgentResult.failed(
                "No 'sites' or 'path'.", elapsed=time.time() - t0
            )
        try:
            sites = ensure_sites(sites_raw, verbose=0)
        except Exception as exc:
            return AgentResult.failed(str(exc), elapsed=time.time() - t0)

        theta = float(input_data.get("strike_deg", self.strike_deg))
        output_dir = input_data.get("output_dir")
        overwrite = bool(input_data.get("overwrite", False))
        suffix = str(input_data.get("file_suffix", "_rot"))

        if output_dir is None:
            return AgentResult.failed(
                "No 'output_dir' specified.",
                hint="Pass output_dir='/path/to/rotated_edis'.",
                elapsed=time.time() - t0,
            )
        os.makedirs(output_dir, exist_ok=True)

        written_paths: list[str] = []
        failed_stations: list[str] = []
        diag_ratios_before: list[float] = []
        diag_ratios_after: list[float] = []

        for i, ed in enumerate(_iter_items(sites)):
            nm = _name(ed, i)
            _, z, fr = _get_z_block(ed)
            if z is None:
                warnings.append(f"{nm}: no Z data — skipped.")
                failed_stations.append(nm)
                continue

            # ── diagonal-suppression metric before rotation ─────────────────
            mask = np.isfinite(z[:, 0, 0]) & (np.abs(z[:, 0, 1]) > 1e-30)
            if mask.any():
                ratio_b = float(
                    np.nanmean(np.abs(z[mask, 0, 0]) / np.abs(z[mask, 0, 1]))
                )
                diag_ratios_before.append(ratio_b)

            # ── rotate Z ───────────────────────────────────────────────────
            try:
                z_rot = rotate_impedance(z, theta)  # (n_freq, 2, 2)
                if z_rot.ndim == 2:
                    z_rot = z_rot[np.newaxis]
            except Exception as exc:
                warnings.append(f"{nm}: rotate_impedance failed: {exc}")
                failed_stations.append(nm)
                continue

            if mask.any():
                ratio_a = float(
                    np.nanmean(
                        np.abs(z_rot[mask, 0, 0]) / np.abs(z_rot[mask, 0, 1])
                    )
                )
                diag_ratios_after.append(ratio_a)

            # ── rotate tipper (if present) ────────────────────────────────
            tip_rot = None
            try:
                tip_raw = getattr(getattr(ed, "Tip", None), "tipper", None)
                if tip_raw is not None:
                    t_arr = np.asarray(tip_raw)
                    if t_arr.ndim == 3 and t_arr.shape[1] == 1:
                        t_arr = t_arr[:, 0, :]
                    if t_arr.shape == (len(fr), 2):
                        tip_rot = rotate_tipper(t_arr, theta)
            except Exception as exc:
                warnings.append(f"{nm}: tipper rotation failed: {exc}")

            # ── write rotated EDI ─────────────────────────────────────────
            out_path = os.path.join(output_dir, f"{nm}{suffix}.edi")
            if os.path.exists(out_path) and not overwrite:
                warnings.append(
                    f"{nm}: {out_path} exists — skipped (set overwrite=True)."
                )
                failed_stations.append(nm)
                continue

            try:
                written = _write_rotated_edi(
                    ed, z_rot, tip_rot, out_path, nm, theta, warnings
                )
                if written:
                    written_paths.append(written)
                else:
                    failed_stations.append(nm)
            except Exception as exc:
                warnings.append(f"{nm}: EDI write failed: {exc}")
                failed_stations.append(nm)

        n_written = len(written_paths)

        # ── diagonal-suppression summary ──────────────────────────────────
        diag_reduction = float("nan")
        if diag_ratios_before and diag_ratios_after:
            diag_reduction = float(
                np.mean(diag_ratios_before) - np.mean(diag_ratios_after)
            )

        # ── figure: before/after Zxy pseudosection ─────────────────────────
        figures: dict[str, Any] = {}
        fig_paths: dict[str, str] = {}
        try:
            fig = _plot_rotation_summary(sites, theta, warnings)
            if fig is not None:
                figures["rotation_summary"] = fig
                p = self._save_figure(
                    fig,
                    output_dir,
                    "tensor_rotation_summary",
                    warnings_list=warnings,
                )
                if p:
                    fig_paths["rotation_summary"] = p
        except Exception as exc:
            warnings.append(f"Rotation summary figure: {exc}")

        # ── LLM interpretation ────────────────────────────────────────────
        interp: str | None = None
        if self.api_key and n_written:
            dr_str = (
                f"{diag_reduction:+.3f}"
                if not np.isnan(diag_reduction)
                else "N/A"
            )
            prompt = (
                f"Tensor rotation summary:\n"
                f"  Rotation angle: {theta:.1f}°\n"
                f"  EDIs written: {n_written}\n"
                f"  Failed stations: {len(failed_stations)}\n"
                f"  Mean |Zxx/Zxy| reduction: {dr_str} (positive = improvement)\n"
                f"  Warnings: {warnings[:3] if warnings else 'none'}\n\n"
                "Evaluate the rotation result and advise on coordinate-frame correction."
            )
            interp = self.query_llm(prompt, max_tokens=200)

        elapsed = time.time() - t0
        dr_disp = (
            f"{diag_reduction:+.3f}"
            if not np.isnan(diag_reduction)
            else "N/A"
        )
        return AgentResult(
            status="success" if n_written > 0 else "needs_review",
            summary=(
                f"Rotation {theta:.1f}°: {n_written} EDIs written, "
                f"{len(failed_stations)} failed. "
                f"Diagonal suppression: {dr_disp}."
            ),
            data={
                "strike_deg": theta,
                "written_paths": written_paths,
                "failed_stations": failed_stations,
                "n_written": n_written,
                "z_diag_reduction": diag_reduction,
                "figures": figures,
                "figure_paths": fig_paths,
            },
            warnings=warnings,
            llm_interpretation=interp,
            elapsed_seconds=elapsed,
            cost_estimate_usd=self._last_cost,
        )


# ── private helpers ───────────────────────────────────────────────────────────


def _write_rotated_edi(
    ed: Any,
    z_rot: np.ndarray,
    tip_rot: np.ndarray | None,
    out_path: str,
    station_name: str,
    theta: float,
    warnings: list[str],
) -> str | None:
    """Write rotated impedance to a new EDI file. Returns path or None."""
    try:
        # Try write_new_edi if the EDI object supports it
        if (
            hasattr(ed, "write_new_edi")
            and hasattr(ed, "Z")
            and ed.Z is not None
        ):
            import copy

            Z_copy = copy.deepcopy(ed.Z)
            Z_copy.z = z_rot
            Tip_copy = None
            if (
                tip_rot is not None
                and hasattr(ed, "Tip")
                and ed.Tip is not None
            ):
                Tip_copy = copy.deepcopy(ed.Tip)
                if hasattr(Tip_copy, "tipper"):
                    Tip_copy.tipper = tip_rot[:, np.newaxis, :]
            return ed.write_new_edi(
                edi_fn=out_path,
                Z=Z_copy,
                Tipper=Tip_copy,
            )

        # Fallback: patch Z on a deepcopy and call write()
        if hasattr(ed, "write"):
            import copy

            ed_copy = copy.deepcopy(ed)
            if hasattr(ed_copy, "Z") and ed_copy.Z is not None:
                ed_copy.Z.z = z_rot
            if (
                tip_rot is not None
                and hasattr(ed_copy, "Tip")
                and ed_copy.Tip is not None
            ):
                ed_copy.Tip.tipper = tip_rot[:, np.newaxis, :]
            return ed_copy.write(new_edifn=out_path)

        warnings.append(f"{station_name}: EDI object has no write method.")
        return None

    except Exception as exc:
        warnings.append(f"{station_name}: _write_rotated_edi: {exc}")
        return None


def _plot_rotation_summary(
    sites: Any, theta: float, warnings: list[str]
) -> Any:
    """Simple bar chart: mean |Zxx/Zxy| per station before rotation."""
    import matplotlib.pyplot as plt

    from ..emtools._core import (
        _get_z_block,
        _iter_items,
        _name,
    )
    from ..seg.ops import rotate_impedance

    station_names, ratios_before, ratios_after = [], [], []

    for i, ed in enumerate(_iter_items(sites)):
        nm = _name(ed, i)
        _, z, fr = _get_z_block(ed)
        if z is None:
            continue
        try:
            mask = np.isfinite(z[:, 0, 0]) & (np.abs(z[:, 0, 1]) > 1e-30)
            if not mask.any():
                continue
            rb = float(
                np.nanmean(np.abs(z[mask, 0, 0]) / np.abs(z[mask, 0, 1]))
            )
            z_r = rotate_impedance(z, theta)
            if z_r.ndim == 2:
                z_r = z_r[np.newaxis]
            ra = float(
                np.nanmean(np.abs(z_r[mask, 0, 0]) / np.abs(z_r[mask, 0, 1]))
            )
            station_names.append(nm)
            ratios_before.append(rb)
            ratios_after.append(ra)
        except Exception:
            continue

    if not station_names:
        return None

    n = len(station_names)
    x = np.arange(n)
    fig, ax = plt.subplots(figsize=(max(6, n * 0.6), 4))
    ax.bar(
        x - 0.18,
        ratios_before,
        width=0.35,
        label="Before",
        color="#3498db",
        alpha=0.8,
    )
    ax.bar(
        x + 0.18,
        ratios_after,
        width=0.35,
        label=f"After θ={theta:.1f}°",
        color="#e74c3c",
        alpha=0.8,
    )
    ax.set_xticks(x)
    ax.set_xticklabels(station_names, rotation=90, fontsize=7)
    ax.set_ylabel("|Zxx| / |Zxy|  (lower = better suppression)", fontsize=8)
    ax.set_title(
        f"Tensor rotation — diagonal suppression  (θ = {theta:.1f}°)",
        fontsize=9,
        fontweight="bold",
    )
    ax.legend(fontsize=8)
    ax.tick_params(labelsize=7)
    fig.tight_layout()
    return fig


__all__ = ["TensorRotationAgent"]
