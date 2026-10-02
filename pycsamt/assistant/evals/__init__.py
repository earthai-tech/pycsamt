# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
r"""
pycsamt.assistant.evals
=======================

Evaluation harnesses for the assistant.

- :mod:`.harness` scores intent/workflow/line classification, retrieval
  and symbol recall, and hallucination guards over JSONL suites in
  ``evals/suites/``. See ``PYCSAMT-V2-RAG-IMPLEMENTATION.md`` §17.
- :mod:`.answers` scores Agent Master final answers and generated
  scripts against the benchmark in ``evals/answer_suites/``;
  :mod:`.runner` runs that benchmark through the real chat path for the
  offline, local and cloud providers.
"""

from __future__ import annotations

from .answers import AnswerReport, load_answer_suite, score_case, score_run
from .harness import EvalReport, evaluate, load_suite

__all__ = [
    "AnswerReport",
    "EvalReport",
    "evaluate",
    "load_answer_suite",
    "load_suite",
    "score_case",
    "score_run",
]
