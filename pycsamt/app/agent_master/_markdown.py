# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
r"""
pycsamt.app.agent_master._markdown
==================================

Pure, dependency-free parsing of the lightweight markdown subset that agents
emit (headings, bullets, fenced code blocks, ``**bold**``).

Kept dash-free on purpose: :mod:`pycsamt.app.agent_master.callbacks.chat`
turns these tokens into Dash ``html`` components, but the *parsing* logic lives
here so it can be unit-tested without importing Dash (or the GUI package).
"""

from __future__ import annotations

import re

__all__ = ["split_inline_bold", "parse_markdown"]


def split_inline_bold(text: str) -> list[tuple[bool, str]]:
    """Split *text* into ``(is_bold, chunk)`` segments on ``**...**``.

    Examples
    --------
    >>> split_inline_bold("a **b** c")
    [(False, 'a '), (True, 'b'), (False, ' c')]
    >>> split_inline_bold("plain")
    [(False, 'plain')]
    >>> split_inline_bold("")
    [(False, '')]
    """
    out: list[tuple[bool, str]] = []
    parts = re.split(r"\*\*(.+?)\*\*", text)
    for i, chunk in enumerate(parts):
        if chunk == "":
            continue
        out.append((i % 2 == 1, chunk))
    return out or [(False, "")]


def parse_markdown(text: str) -> list[tuple]:
    """Tokenise lightweight markdown into block tokens.

    Returns a list of tuples whose first element is the block kind:
    ``("code", lang, body)``, ``("heading", text)``, ``("bullet", text)``,
    ``("para", text)``, or ``("blank",)``.

    Only the small subset agents actually emit is supported: fenced code
    blocks, ``#``-headings, ``-``/``*`` bullets, and paragraphs.

    Examples
    --------
    >>> parse_markdown("# Title\\n- a\\n```py\\nx=1\\n```")
    [('heading', 'Title'), ('bullet', 'a'), ('code', 'py', 'x=1')]
    >>> parse_markdown("hello **world**")
    [('para', 'hello **world**')]
    """
    tokens: list[tuple] = []
    lines = (text or "").split("\n")
    i = 0
    n = len(lines)
    while i < n:
        line = lines[i]
        stripped = line.strip()
        if stripped.startswith("```"):
            lang = stripped[3:].strip() or "python"
            body: list[str] = []
            i += 1
            while i < n and not lines[i].strip().startswith("```"):
                body.append(lines[i])
                i += 1
            i += 1  # skip closing fence
            tokens.append(("code", lang, "\n".join(body)))
            continue
        if not stripped:
            tokens.append(("blank",))
        elif stripped.startswith("#"):
            tokens.append(("heading", stripped.lstrip("#").strip()))
        elif stripped[:2] in ("- ", "* "):
            tokens.append(("bullet", stripped[2:].strip()))
        else:
            tokens.append(("para", stripped))
        i += 1
    return tokens
