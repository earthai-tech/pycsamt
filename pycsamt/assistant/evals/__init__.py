# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
r"""
pycsamt.assistant.evals
=======================

Evaluation harness for the RAG assistant: score intent/workflow/line
classification, retrieval & symbol recall, and hallucination guards over
JSONL suites in ``evals/suites/``. See
``PYCSAMT-V2-RAG-IMPLEMENTATION.md`` §17.
"""

from __future__ import annotations

from .harness import EvalReport, evaluate, load_suite

__all__ = ["EvalReport", "evaluate", "load_suite"]
