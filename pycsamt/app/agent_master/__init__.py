# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
r"""
pycsamt.app.agent_master
========================

pyCSAMT Agent Master GUI — the single
entry-point to all pyCSAMT sub-agents.

Launch
------
From Python::

    from pycsamt.app.agent_master import launch
    launch()   # opens http://localhost:8765

From the CLI::

    python -m pycsamt.app.agent_master

Features
--------
- Natural-language workflow dispatch to the
  full pyCSAMT agent chain.
- Top-bar EDI loader (files or folder) with
  3-mode survey-line detection.
- Persistent session save / reload.
- Settings panel: Claude / OpenAI / Gemini /
  DeepSeek / MiniMax API keys, model selection,
  default export format.
- Inline figure cards with full-screen viewer
  and export to PNG / SVG / EPS / PDF.
- Light and dark themes (CSS custom properties,
  persisted in browser localStorage).
"""

from __future__ import annotations

from .app import create_app, launch

__all__ = ["create_app", "launch"]
