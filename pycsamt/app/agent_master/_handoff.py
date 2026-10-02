# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Desktop → Agent Master survey hand-off (no Qt, no Dash).

When Agent Master is launched from the desktop with a survey loaded, the
desktop writes the survey *as it is in memory* -- corrections, edits and
pipeline output included -- to a hand-off folder, one sub-folder per line,
with a ``session.json`` in Agent Master's own ``STORE_EDI`` shape.  The
browser opens ``/?handoff=<token>`` and Agent Master starts on that data
instead of the empty welcome screen, so a long desktop session is never
re-read from the original files (or lost).

Folders live in ``~/.pycsamt/agent_master/handoff/<token>/`` (override:
``PYCSAMT_AGENT_HANDOFF_DIR``); only the newest few are kept.
"""

from __future__ import annotations

import datetime as _dt
import json
import os
import re
import secrets
import shutil
from pathlib import Path
from typing import Any

__all__ = [
    "HANDOFF_PARAM",
    "handoff_root",
    "read_handoff",
    "token_from_search",
    "write_handoff",
]

HANDOFF_PARAM = "handoff"
_TOKEN_RE = re.compile(r"^[0-9a-f]{12,40}$")
_SESSION = "session.json"


def handoff_root() -> Path:
    env = os.environ.get("PYCSAMT_AGENT_HANDOFF_DIR")
    return Path(env) if env else (Path.home() / ".pycsamt" / "agent_master"
                                  / "handoff")


def _folder_name(line: str) -> str:
    name = re.sub(r"[^\w.-]+", "_", str(line).strip()) or "line"
    return name[:60]


def _site_name(site) -> str:
    return str(getattr(site, "name", "") or getattr(site, "station", "")
               or "site")


def write_handoff(sites, lines: dict | None = None, *, edited: bool = False,
                  label: str = "", keep: int = 5) -> str:
    """Write *sites* for Agent Master; returns the hand-off token.

    Parameters
    ----------
    sites : Sites or iterable of Site
        The survey as the desktop holds it now.
    lines : dict, optional
        ``{station: line}``; stations without a line go to ``"Survey"``.
    edited : bool
        The survey was changed since it was loaded (shown in Agent
        Master's note).
    label : str
        Where the data came from (e.g. the loaded folder).
    keep : int
        Number of recent hand-offs to keep.
    """
    from pycsamt.site.base import Sites

    items = list(sites)
    if not items:
        raise ValueError("no stations to hand over")
    root = handoff_root()
    # timestamp (to the microsecond) first: name order = creation order
    token = _dt.datetime.now().strftime("%y%m%d%H%M%S%f") + secrets.token_hex(4)
    folder = root / token
    groups: dict[str, list] = {}
    for s in items:
        line = str((lines or {}).get(_site_name(s)) or "Survey")
        if line in ("nan", "None"):
            line = "Survey"
        groups.setdefault(line, []).append(s)
    written: dict[str, list[str]] = {}
    for line, members in groups.items():
        paths = Sites(members).write(folder / _folder_name(line),
                                     exist_ok=True)
        written[line] = [str(p) for p in paths]
    n = sum(len(v) for v in written.values())
    session = {
        "path": str(folder),
        "groups": written,
        "n_edi": n,
        "mode": "desktop",
        "source": "desktop",
        "edited": bool(edited),
        "label": str(label or ""),
        "created": _dt.datetime.now().isoformat(timespec="seconds"),
    }
    (folder / _SESSION).write_text(json.dumps(session, indent=2),
                                   encoding="utf-8")
    _prune(root, keep)
    return token


def _prune(root: Path, keep: int) -> None:
    try:
        olds = sorted((p for p in root.iterdir()
                       if p.is_dir() and _TOKEN_RE.match(p.name)),
                      key=lambda p: p.name)
    except OSError:
        return
    for p in olds[:-max(1, keep)]:
        shutil.rmtree(p, ignore_errors=True)


def token_from_search(search: str | None) -> str:
    """The hand-off token of a URL query string (``?handoff=…``)."""
    if not search:
        return ""
    from urllib.parse import parse_qs

    vals = parse_qs(str(search).lstrip("?")).get(HANDOFF_PARAM) or [""]
    tok = vals[0].strip()
    return tok if _TOKEN_RE.match(tok) else ""


def read_handoff(token: str) -> dict[str, Any] | None:
    """The ``STORE_EDI``-shaped session of *token*, or ``None``.

    Only tokens of the expected form, inside the hand-off root, whose
    files still exist are accepted.
    """
    if not token or not _TOKEN_RE.match(token):
        return None
    root = handoff_root().resolve()
    folder = (root / token).resolve()
    if root not in folder.parents:
        return None
    try:
        data = json.loads((folder / _SESSION).read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None
    groups = {str(k): [str(f) for f in v if Path(f).exists()]
              for k, v in (data.get("groups") or {}).items()}
    groups = {k: v for k, v in groups.items() if v}
    if not groups:
        return None
    data["groups"] = groups
    data["n_edi"] = sum(len(v) for v in groups.values())
    data["path"] = str(folder)
    return data
