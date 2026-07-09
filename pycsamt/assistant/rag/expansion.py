# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
r"""
pycsamt.assistant.rag.expansion
===============================

Deterministic **query expansion** — the zero-dependency "semantic" layer.

A lexical retriever (BM25) can only match words the user actually typed.
But a geophysicist rarely types the API vocabulary: they ask *"clean up
spikes and outliers"* (the code says ``rpca`` / ``hampel`` / ``denoise``)
or *"how far down can I trust my model"* (the code says ``sensitivity`` /
``bostick`` / depth of investigation). This module bridges that gap
**without embeddings**: it maps natural-language triggers to the domain
terms that actually appear in the corpus, so the retriever can find the
right chunk even when the query and the source share no surface words.

It is intentionally high-precision — a trigger only fires on a specific
phrase, so an off-topic query expands to nothing and retrieval stays
empty rather than drifting. The expansion terms are scored at a lower
weight than the user's own words (see :mod:`~pycsamt.assistant.rag.
retriever`), so they *nudge* ranking rather than dominate it.
"""

from __future__ import annotations

import re

__all__ = ["QUERY_EXPANSIONS", "expand_query"]


# Each entry maps a trigger (matched as a lowercased substring / word) to
# extra domain terms — all of which are real tokens found in the pyCSAMT
# corpus (workflow ids, method names, recipe/doc vocabulary). Ordered for
# readability; every matching trigger contributes its terms.
QUERY_EXPANSIONS: tuple[tuple[str, tuple[str, ...]], ...] = (
    # ── denoising / outliers ────────────────────────────────────────────
    ("spike", ("denoise", "denoising", "rpca", "hampel", "outlier")),
    ("outlier", ("denoise", "denoising", "rpca", "hampel")),
    ("despike", ("denoise", "rpca", "hampel")),
    ("glitch", ("denoise", "outlier", "hampel")),
    ("clean up", ("denoise", "denoising", "outlier")),
    ("remove noise", ("denoise", "denoising", "rpca")),
    ("noise removal", ("denoise", "denoising", "rpca")),
    # ── quality control ─────────────────────────────────────────────────
    ("noisy", ("qc", "snr", "denoise")),
    ("bad frequenc", ("qc", "frequency", "dead", "band", "snr")),
    ("dead band", ("qc", "snr", "frequency")),
    ("drop frequenc", ("qc", "frequency", "decimation")),
    ("data quality", ("qc", "snr", "flags")),
    ("flag", ("qc", "snr")),
    # ── static shift / galvanic ─────────────────────────────────────────
    ("shifted", ("static", "shift", "galvanic", "ama")),
    ("vertical offset", ("static", "shift", "galvanic")),
    ("offset vertical", ("static", "shift", "galvanic")),
    ("parallel curve", ("static", "shift", "galvanic")),
    ("level shift", ("static", "shift", "galvanic")),
    ("distortion", ("galvanic", "static", "shift", "correction")),
    # ── phase tensor / strike / dimensionality ──────────────────────────
    ("out of quadrant", ("phase", "tensor", "dimensionality")),
    ("quadrant", ("phase", "tensor")),
    ("phase wrap", ("phase", "tensor")),
    ("which direction", ("strike", "dimensionality", "azimuth")),
    ("geology run", ("strike", "azimuth")),
    ("orientation", ("strike", "azimuth", "rotation")),
    ("geoelectric strike", ("strike", "dimensionality")),
    ("1d 2d 3d", ("dimensionality", "skew", "phase", "tensor")),
    ("1-d 2-d 3-d", ("dimensionality", "skew", "phase", "tensor")),
    ("dimensional", ("dimensionality", "skew", "strike")),
    # ── sensitivity / depth of investigation ────────────────────────────
    ("how far down", ("sensitivity", "bostick", "depth", "investigation")),
    ("how deep", ("sensitivity", "bostick", "depth", "investigation")),
    ("penetrat", ("sensitivity", "bostick", "depth", "investigation")),
    ("trust", ("sensitivity", "bostick", "investigation")),
    ("reliable depth", ("sensitivity", "bostick", "investigation")),
    ("investigation depth", ("sensitivity", "bostick")),
    # ── inversion / modelling ───────────────────────────────────────────
    ("without a mesh", ("aiinversion", "eminverter", "neural", "inversion")),
    ("no mesh", ("aiinversion", "eminverter", "inversion")),
    ("resistivity model", ("inversion", "occam", "model")),
    ("depth model", ("inversion", "sensitivity", "bostick")),
    ("third party", ("occam2d", "modem", "mare2dem", "export")),
    ("2d inversion", ("occam2d", "modem", "inv2d", "mesh")),
    ("3d inversion", ("modem", "inv3d")),
    # ── tipper / induction ──────────────────────────────────────────────
    ("induction arrow", ("tipper", "wiese", "parkinson")),
    ("tipper", ("tipper", "induction", "wiese")),
    # ── profile / geometry ──────────────────────────────────────────────
    ("combine", ("profile", "section", "collection")),
    ("merge line", ("profile", "section", "collection")),
    ("stitch", ("profile", "section")),
    # ── loading ─────────────────────────────────────────────────────────
    ("load edi", ("ensure_sites", "sites", "edi", "collection")),
    ("read edi", ("ensure_sites", "sites", "edi")),
)

_WORD = re.compile(r"[a-z0-9]+")


def expand_query(query: str) -> list[str]:
    """Return extra domain terms implied by *query* (possibly empty).

    High-precision: a term is only added when its trigger phrase appears
    in the lowercased query, and duplicates / terms already present in the
    query are dropped so expansion never double-counts a word the user
    already typed.
    """
    ql = (query or "").lower()
    present = set(_WORD.findall(ql))
    out: list[str] = []
    seen: set[str] = set()
    for trigger, terms in QUERY_EXPANSIONS:
        if trigger in ql:
            for t in terms:
                if t not in seen and t not in present:
                    seen.add(t)
                    out.append(t)
    return out
