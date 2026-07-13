"""
pycsamt.agents.edi_export
==========================

:class:`EDIExportAgent` — Write processed Sites back to EDI files on disk.

Allows any corrected or processed dataset (static-shift corrected, denoised,
frequency-decimated) to be exported as standard SEG EDI files, making the
results reusable by external tools such as Occam2D GUI, WinGLink, WALDIM,
or any other MT processing software.

Wraps :meth:`~pycsamt.seg.collection.EDICollection.export` when Sites can be
cast to an EDICollection, and falls back to per-item EDI writing otherwise.
"""

from __future__ import annotations

import os
import time
from typing import Any

from ._base import AgentResult, BaseAgent

_SYSTEM_PROMPT = """\
You are an expert in MT data management and EDI file format conventions.
Given an EDI export result, write 2-3 sentences that:
1. Confirm how many files were written and their location.
2. Note any stations that failed and the likely cause.
3. Recommend next steps for the exported data (e.g. inversion, external QC).
Reply in plain English.
"""


class EDIExportAgent(BaseAgent):
    """Write processed Sites to EDI files on disk.

    Parameters
    ----------
    api_key, model, llm_provider : str
    file_pattern : str
        Filename format using ``{station}`` placeholder.
        Default ``"{station}.edi"``.
    overwrite : bool
        Overwrite existing files (default False).

    Input keys
    ----------
    ``sites`` / ``path`` : Sites or str
    ``output_dir`` : str — target directory (created if absent)
    ``file_pattern`` : str, optional — override default
    ``overwrite`` : bool, optional — override default

    Output data keys
    ----------------
    ``written_paths``   list[str]
    ``failed``          list[tuple[str, str]]  — (station, error message)
    ``n_written``       int
    ``n_failed``        int
    ``output_dir``      str

    Examples
    --------
    >>> agent = EDIExportAgent()
    >>> r = agent.execute({
    ...     "path":       "/data/WILLY_EDIs",
    ...     "output_dir": "/out/willy_corrected",
    ... })
    >>> print(r["n_written"], "EDIs exported")
    """

    SYSTEM_PROMPT = _SYSTEM_PROMPT

    def __init__(
        self,
        *,
        api_key: str | None = None,
        model: str | None = None,
        llm_provider: str = "claude",
        file_pattern: str = "{station}.edi",
        overwrite: bool = False,
    ) -> None:
        super().__init__(
            "EDIExportAgent",
            api_key=api_key,
            model=model,
            llm_provider=llm_provider,
        )
        self.file_pattern = file_pattern
        self.overwrite = overwrite

    def execute(self, input_data: dict[str, Any]) -> AgentResult:
        self._last_cost = 0.0
        t0 = time.time()
        warnings: list[str] = []

        from ..emtools._core import (
            ensure_sites,
        )

        sites_raw = input_data.get("sites") or input_data.get("path")
        if sites_raw is None:
            return AgentResult.failed(
                "No 'sites' or 'path'.", elapsed=time.time() - t0
            )

        output_dir = input_data.get("output_dir")
        if output_dir is None:
            return AgentResult.failed(
                "No 'output_dir' specified.",
                hint="Pass output_dir='/path/to/export_dir'.",
                elapsed=time.time() - t0,
            )

        os.makedirs(output_dir, exist_ok=True)

        try:
            sites = ensure_sites(sites_raw, verbose=0)
        except Exception as exc:
            return AgentResult.failed(str(exc), elapsed=time.time() - t0)

        pattern = str(input_data.get("file_pattern", self.file_pattern))
        overwrite = bool(input_data.get("overwrite", self.overwrite))

        written_paths: list[str] = []
        failed: list[tuple[str, str]] = []

        # ── try EDICollection.export() first ──────────────────────────────
        if hasattr(sites, "export"):
            try:
                result = sites.export(
                    output_dir,
                    file_pattern=pattern,
                )
                written_paths = [str(p) for p in result.get("successful", [])]
                for sid, err in result.get("failed", []):
                    failed.append((str(sid), str(err)))
                    warnings.append(f"{sid}: export failed — {err}")
            except Exception as exc:
                warnings.append(
                    f"EDICollection.export() failed ({exc}); "
                    "falling back to per-item write."
                )
                written_paths, failed = _per_item_export(
                    sites, output_dir, pattern, overwrite, warnings
                )
        else:
            written_paths, failed = _per_item_export(
                sites, output_dir, pattern, overwrite, warnings
            )

        n_written = len(written_paths)
        n_failed = len(failed)

        # ── LLM interpretation ────────────────────────────────────────────
        interp: str | None = None
        if self.api_key:
            prompt = (
                f"EDI export summary:\n"
                f"  Output directory: {output_dir}\n"
                f"  Files written: {n_written}\n"
                f"  Files failed: {n_failed}\n"
                f"  Pattern: {pattern!r}\n"
                f"  Warnings: {warnings[:3] if warnings else 'none'}\n\n"
                "Confirm the export and recommend next steps."
            )
            interp = self.query_llm(prompt, max_tokens=150)

        elapsed = time.time() - t0
        status = "success" if n_written > 0 else "failed"
        return AgentResult(
            status=status,
            summary=(
                f"EDI export: {n_written} files written to '{output_dir}'. "
                f"{n_failed} failed."
            ),
            data={
                "written_paths": written_paths,
                "failed": failed,
                "n_written": n_written,
                "n_failed": n_failed,
                "output_dir": output_dir,
                "figures": {},
                "figure_paths": {},
            },
            warnings=warnings,
            llm_interpretation=interp,
            elapsed_seconds=elapsed,
            cost_estimate_usd=self._last_cost,
        )


# ── private helpers ───────────────────────────────────────────────────────────


def _per_item_export(
    sites: Any,
    output_dir: str,
    pattern: str,
    overwrite: bool,
    warnings: list[str],
) -> tuple[list[str], list[tuple[str, str]]]:
    """Write each item individually when bulk export is unavailable."""
    from ..emtools._core import _iter_items, _name

    written: list[str] = []
    failed: list[tuple[str, str]] = []

    for i, ed in enumerate(_iter_items(sites)):
        nm = _name(ed, i)
        try:
            filename = pattern.format(station=nm)
        except (KeyError, ValueError):
            filename = f"{nm}.edi"

        out_path = os.path.join(output_dir, filename)

        if os.path.exists(out_path) and not overwrite:
            warnings.append(
                f"{nm}: {out_path} exists — skipped (set overwrite=True)."
            )
            failed.append((nm, "file exists"))
            continue

        try:
            if hasattr(ed, "write"):
                path = ed.write(new_edifn=out_path)
                written.append(str(path or out_path))
            elif hasattr(ed, "write_new_edi"):
                path = ed.write_new_edi(edi_fn=out_path)
                written.append(str(path or out_path))
            else:
                failed.append((nm, "no write method"))
                warnings.append(f"{nm}: no write method available.")
        except Exception as exc:
            failed.append((nm, str(exc)))
            warnings.append(f"{nm}: write failed — {exc}")

    return written, failed


__all__ = ["EDIExportAgent"]
