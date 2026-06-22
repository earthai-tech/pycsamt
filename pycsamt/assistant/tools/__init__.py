# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
r"""
pycsamt.assistant.tools
=======================

Executable tools the assistant can call (as opposed to RAG, which only
retrieves context). Step 2 adds the project registry; later steps add
data / static-shift / QC / plot / inversion / validation tools.
"""

from __future__ import annotations

from .project_registry import ProjectRegistry, resolve_line
from .validation_tools import validate_generated_code
from .workflow_tools import (
    resolve_target,
    run_phase_analysis,
    run_qc,
    run_static_shift,
    run_workflow,
)

__all__ = [
    "ProjectRegistry",
    "resolve_line",
    "validate_generated_code",
    "resolve_target",
    "run_workflow",
    "run_static_shift",
    "run_qc",
    "run_phase_analysis",
]
