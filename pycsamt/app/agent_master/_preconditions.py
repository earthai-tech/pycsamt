# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Deterministic request preconditions checked before any dispatch.

Each check returns ``(kind, message)`` when the request cannot be served as
asked — a named pycsamt API that does not exist, an input path that is
absent, or literal correction factors that are not finite and positive —
and ``None`` otherwise. They never run a workflow or a model.
"""

from __future__ import annotations

import math
import re
from pathlib import Path

__all__ = ["check_request", "no_data_message"]

_API_REF = re.compile(r"(?<![\w/.])pycsamt(?:\.[A-Za-z_]\w*)+")
_PATH = re.compile(
    r"(?<![\w/\\])(?:[A-Za-z]:[\\/]|/)[^\s'\"`<>|?*]+"
    r"\.(?:edi|xml|j|dat|avg|csv|txt|h5|pcsf)\b",
    re.IGNORECASE,
)
_FACTORS = re.compile(
    r"\b(?:apply|use|correct\w*)\b[^.]*?\bfactors?\b\s*(?:[:=]|of|are)?\s*(.+)",
    re.IGNORECASE,
)
_WORD_VALUES = {"nan": math.nan, "inf": math.inf, "infinity": math.inf,
                "zero": 0.0, "one": 1.0, "two": 2.0, "three": 3.0}
_LINE_TOKEN = re.compile(r"\bL\d+[A-Z]*\b")
_INVERSIONS = {"ai_inversion", "inv1d", "inv2d", "inv3d", "pinn_inversion",
               "hybrid_inversion", "ensemble_inversion"}
_SOLVERS = {"modem": "ModEM", "mare2dem": "MARE2DEM", "pre_inversion": "Occam2D"}


def _unresolved_api(text):
    refs = sorted({m.rstrip(".") for m in _API_REF.findall(text)})
    if not refs:
        return None
    from pycsamt.assistant.tools._static_api import StaticAPI

    api = StaticAPI()
    missing = [r for r in refs if api.resolve(r)[0] == "failed"]
    if not missing:
        return None
    names = ", ".join(f"`{m}`" for m in missing)
    return ("answer",
            f"{names} does not exist in this pyCSAMT checkout, so I will not "
            "import or call it. To process a survey, load the data and ask for "
            "a registered workflow (for example quality control, static-shift "
            "correction or phase tensor analysis), or ask me to write a script "
            "using verified pyCSAMT functions.")


def _missing_path(text):
    for match in _PATH.finditer(text):
        path = match.group(0).rstrip(".,;:)")
        if path.startswith("//") or text[:match.start()].endswith(":"):
            continue  # URL, not a local path
        if not Path(path).exists():
            return ("clarify",
                    f"The input `{path}` does not exist, so nothing was read "
                    "and no station statistics can be reported. Which existing "
                    "EDI file, folder or registered survey line should I use?")
    return None


def _factor_values(fragment):
    values = []
    for token in re.findall(r"minus\s+\w+|-?\d+(?:\.\d+)?(?:e-?\d+)?|[A-Za-z]+",
                            fragment, re.IGNORECASE):
        low = token.lower()
        sign = 1.0
        if low.startswith("minus"):
            sign, low = -1.0, low.split()[-1]
        if low in _WORD_VALUES:
            values.append(sign * _WORD_VALUES[low])
        else:
            try:
                values.append(sign * float(low))
            except ValueError:
                continue
    return values


def _invalid_factors(text):
    match = _FACTORS.search(text)
    if not match:
        return None
    values = _factor_values(match.group(1))
    bad = [v for v in values if not (math.isfinite(v) and v > 0)]
    if not bad:
        return None
    return ("answer",
            "These static-shift factors cannot be applied. Correction factors "
            "multiply resistivity, so each must be finite and strictly "
            "positive: NaN or infinite values are undefined, zero removes the "
            "signal, and negative values have no physical meaning. No "
            "correction was applied and the source data are unchanged. "
            "Re-estimate the factors (for example with AMA) or supply finite "
            "positive values.")


def check_request(text):
    """Return ``(kind, message)`` for an unservable request, else ``None``."""
    for check in (_unresolved_api, _missing_path, _invalid_factors):
        found = check(text or "")
        if found:
            return found
    return None


_ANGLE = re.compile(r"(?<![\w.])([-+]?\d+(?:\.\d+)?)\s*(?:°|deg(?:ree)?s?\b)", re.I)
_CLOCKWISE = re.compile(r"\bclockwise\b|\beast\s+of\s+north\b|\bazimuth\b|\bN\s*→?\s*E\b", re.I)
_ANTICLOCKWISE = re.compile(
    r"\b(?:counter[- ]?clockwise|anti[- ]?clockwise)\b|\bwest\s+of\s+north\b", re.I)
_TO_STRIKE = re.compile(r"\b(?:to|onto|along|into)\s+(?:the\s+)?(?:estimated\s+|geoelectric\s+)?"
                        r"(?:strike|principal\s+axes)\b", re.I)
ROTATION_CONVENTION = (
    "In pyCSAMT a positive angle rotates the measurement axes clockwise in map view, "
    "from north toward east (x toward y), relative to the data's current frame.")


def rotation_request(text):
    """Resolve the angle of a rotation request, or the question to ask.

    Returns ``{"angle": float | None, "question": str | None,
    "suggestions": list}``; suggestions are one-click replies to the question. ``angle`` is
    in pyCSAMT's convention (positive = clockwise from north); ``None`` with
    no question means "rotate to the estimated strike". An explicit angle
    without a stated sense is ambiguous, so it is asked about rather than
    guessed.
    """
    text = text or ""
    match = _ANGLE.search(text)
    if _TO_STRIKE.search(text) or not match:
        return {"angle": None, "question": None, "suggestions": []}
    value = float(match.group(1))
    anticlockwise = _ANTICLOCKWISE.search(text)
    if anticlockwise:
        return {"angle": -abs(value), "question": None, "suggestions": []}
    if _CLOCKWISE.search(text):
        return {"angle": abs(value), "question": None, "suggestions": []}
    return {"angle": None, "question": (
        f"In which sense should I rotate by {abs(value):g}°? {ROTATION_CONVENTION} "
        f"Reply, for example, \"rotate by {abs(value):g} degrees clockwise\", "
        f"\"rotate by {abs(value):g} degrees counter-clockwise\", or \"rotate to the "
        "strike\" to use the estimated geoelectric strike. The original EDI files are "
        "kept; rotated copies are written."),
        "suggestions": [
            {"label": f"{abs(value):g}° clockwise",
             "reply": f"Rotate the data by {abs(value):g} degrees clockwise."},
            {"label": f"{abs(value):g}° counter-clockwise",
             "reply": f"Rotate the data by {abs(value):g} degrees counter-clockwise."},
            {"label": "Rotate to the strike", "reply": "Rotate the data to the strike."},
        ]}


def line_suggestions(text, registry=None):
    """Replies that retry *text* on registered lines when it names an unknown one."""
    known = list(registry.lines()) if registry is not None else []
    unknown = [t for t in _LINE_TOKEN.findall(text or "")
               if t.upper() not in {k.upper() for k in known}]
    if not unknown or not known:
        return []
    return [{"label": line, "reply": re.sub(re.escape(unknown[0]), line, text, count=1)}
            for line in known[:4]]


def no_data_message(workflow, text, label, registry=None):
    """Explain what a data-dependent *workflow* needs, specific to the request."""
    known = list(registry.lines()) if registry is not None else []
    named = [t for t in _LINE_TOKEN.findall(text or "")
             if t.upper() not in {k.upper() for k in known}]
    if named:
        listing = f" Known lines: {', '.join(known[:12])}." if known else ""
        return (f"Line {named[0]} is not in the project registry, so no data "
                f"were loaded and nothing was run.{listing} Which line or EDI "
                "folder should I use?")
    if workflow in _SOLVERS:
        solver = _SOLVERS[workflow]
        return (f"Running {solver} needs loaded station data and an installed "
                f"{solver} executable (see `pycsamt build`); nothing was run. "
                f"Once data are loaded I can prepare the {solver} input files, "
                "and run the solver only if it is installed. Which dataset "
                "should I use?")
    if workflow in _INVERSIONS:
        return ("Which data should I invert, and at what scope? Options include "
                "a 1-D AI inversion per station, 2-D inversion preparation "
                "(Occam2D, MARE2DEM) or 3-D ModEM preparation. Load an EDI or "
                "XML-TF dataset or name a registered line; I will not start an "
                "inversion with invented data or settings.")
    return (f"{label[:1].upper() + label[1:]} needs station data, and none is "
            "loaded. Load an EDI or XML-TF dataset with Load Data, or name a "
            "registered survey line, then ask again.")
