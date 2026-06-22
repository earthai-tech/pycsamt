# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
r"""
pycsamt.assistant.rag.config
============================

What to index, what to exclude, and how to tag/prioritise it.

Keep this declarative — the indexers and ingest walker consult these
constants and helpers so the policy lives in one place.
"""

from __future__ import annotations

from pathlib import PurePosixPath

__all__ = [
    "EXCLUDE_PATTERNS",
    "INDEX_ROOTS",
    "ROOT_DOCS",
    "INDEX_SUFFIXES",
    "HIGH_PRIORITY_PATHS",
    "WORKFLOW_KEYWORDS",
    "is_excluded",
    "should_index",
    "priority_for",
    "infer_workflow",
]

# ── exclusions (noise that pollutes retrieval) ──────────────────────────────────
EXCLUDE_PATTERNS = (
    "__pycache__",
    ".pyc",
    "docs/source/_static",
    "docs/source/images",
    "pycsamt/app/",                 # Dash/desktop UI glue, not the science API
    "pycsamt/app/desktop/resources",
    "pycsamt/app/web/assets",
    ".pytest_cache",
    ".ruff_cache",
    ".git/",
    "build/",
    "dist/",
    ".egg-info",
    "/paper/",
    "cag-paper",
)

# Only files under these roots (relative to the repo root) are indexed,
# plus the root-level docs in ROOT_DOCS.
INDEX_ROOTS = (
    "pycsamt",
    "docs/source",
    "examples",
    "assistant_recipes",
)
ROOT_DOCS = ("README.md", "AGENTS.md")

INDEX_SUFFIXES = (".py", ".rst", ".md")

# Implementation areas worth boosting at retrieval time.
HIGH_PRIORITY_PATHS = (
    "pycsamt/agents",
    "pycsamt/emtools",
    "pycsamt/api",
    "pycsamt/site",
    "pycsamt/seg",
    "pycsamt/z",
    "pycsamt/pipeline",
    "pycsamt/forward",
    "pycsamt/ai",
    "assistant_recipes",
)

# ── workflow tagging (keep aligned with pycsamt.agents._workflows) ──────────────
WORKFLOW_KEYWORDS: dict[str, tuple[str, ...]] = {
    "static_shift": (
        "static shift", "static-shift", "galvanic",
        "estimate_ss_ama", "correct_ss_ama", "staticshift",
    ),
    "qc": (
        "quality control", "dataqc", "qc_flags",
        "build_qc_table", "dead band", "snr",
    ),
    "phase_analysis": (
        "phase tensor", "phaseanalysis", "strike",
        "skew", "dimensionality", "mohr", "argand",
    ),
    "denoise": ("denoise", "denoising", "rpca", "hampel"),
    "tipper": ("tipper", "induction arrow", "wiese", "parkinson"),
    "sensitivity": ("sensitivity", "bostick", "depth of investigation"),
    "rotation": ("tensor rotation", "rotate", "strike rotation"),
    "ai_inversion": (
        "ai inversion", "eminverter", "aiinversion",
        "neural inversion", "cnn inversion",
    ),
    "inv2d": ("unet", "inv2d", "2d inversion"),
    "inv3d": ("gcn", "inv3d", "3d inversion"),
    "pre_inversion": ("occam2d", "occam", "mesh", "startup"),
    "modem": ("modem",),
    "forward": ("forward model", "mt1dforward", "layeredmodel"),
    "report": ("survey report", "reportagent"),
    "code_gen": ("code generation", "codegenerationagent"),
}


# ── helpers ─────────────────────────────────────────────────────────────────────

def _posix(rel_path: str) -> str:
    return PurePosixPath(rel_path).as_posix()


def is_excluded(rel_path: str) -> bool:
    """True when *rel_path* (repo-relative) matches any exclusion."""
    p = _posix(rel_path)
    return any(pat in p for pat in EXCLUDE_PATTERNS)


def should_index(rel_path: str) -> bool:
    """True when *rel_path* should be ingested.

    Must have an indexable suffix, not be excluded, and live under an
    :data:`INDEX_ROOTS` root (or be a root-level doc in :data:`ROOT_DOCS`).
    """
    p = _posix(rel_path)
    if not p.endswith(INDEX_SUFFIXES):
        return False
    if is_excluded(p):
        return False
    if p in ROOT_DOCS:
        return True
    return any(
        p == root or p.startswith(root + "/")
        for root in INDEX_ROOTS
    )


def priority_for(rel_path: str) -> int:
    """Static priority: 3 for high-priority areas, 2 for docs, else 1."""
    p = _posix(rel_path)
    if any(p.startswith(h) for h in HIGH_PRIORITY_PATHS):
        return 3
    if p.startswith("docs/source") or p in ROOT_DOCS:
        return 2
    return 1


def infer_workflow(*texts: str) -> str | None:
    """Best-effort workflow tag from path / symbol / docstring text.

    Returns the first workflow whose keyword appears in the lowered,
    concatenated *texts*, else ``None``.
    """
    blob = " ".join(t for t in texts if t).lower()
    if not blob:
        return None
    for wf, kws in WORKFLOW_KEYWORDS.items():
        if any(kw in blob for kw in kws):
            return wf
    return None
