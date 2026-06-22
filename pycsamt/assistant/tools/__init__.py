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

__all__ = ["ProjectRegistry", "resolve_line"]
