# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
r"""
pycsamt.assistant.tools.project_registry
========================================

Resolve a survey-line name (e.g. ``"L22PLT"``) to a real EDI directory
plus default processing parameters.

Without this, the assistant can only guess what ``line L22PLT`` means.
With it, a user phrase like *"process line L22PLT"* resolves to::

    {
      "line": "L22PLT",
      "edi_dir": "<abs path>",
      "exists": True,
      "n_edi_files": 25,
      "sort_by": "lon",
      "output_root": "results/willy/L22PLT",
      "default_workflows": [...],
      "static_shift": {...defaults...},
      "plot": {...defaults...},
    }

Registries are YAML; see ``projects/willy_project_registry.yml`` and
``projects/example_project_registry.yml``.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any

import yaml

__all__ = ["ProjectRegistry", "resolve_line", "find_default_registry"]

# Recognise a line token in free text: "L22PLT", "line L22PLT", "l22plt".
_LINE_TOKEN = re.compile(r"\b([A-Za-z]{1,3}\d{1,3}[A-Za-z]{0,4})\b")


def _repo_root() -> Path:
    here = Path(__file__).resolve()
    for parent in here.parents:
        if (parent / "pycsamt").is_dir() and (
            parent / "pyproject.toml"
        ).exists():
            return parent
    return here.parents[3]


def find_default_registry(root: Path | None = None) -> Path | None:
    """Locate a registry under ``projects/`` (prefers ``willy_*``).

    Returns the first match or ``None`` if no real registry is present
    (the ``example_*`` template is ignored).
    """
    root = root or _repo_root()
    proj = root / "projects"
    if not proj.is_dir():
        return None
    willy = proj / "willy_project_registry.yml"
    if willy.exists():
        return willy
    for p in sorted(proj.glob("*_project_registry.yml")):
        if p.name != "example_project_registry.yml":
            return p
    return None


def _norm(name: str) -> str:
    """Normalise a line name for matching: strip 'line', upper, alnum."""
    s = (name or "").strip().lower()
    s = re.sub(r"^(survey\s+)?line\s+", "", s)
    s = re.sub(r"[^a-z0-9]", "", s)
    return s.upper()


class ProjectRegistry:
    """Load a project registry and resolve survey lines.

    Parameters
    ----------
    path : path-like
        Path to a registry YAML file.
    root : path-like, optional
        Base directory for resolving relative ``edi_dir`` paths.
        Defaults to the repository root.
    """

    def __init__(
        self, path: str | Path, *, root: str | Path | None = None
    ) -> None:
        self.path = Path(path)
        self.root = Path(root) if root is not None else _repo_root()
        self.data: dict[str, Any] = (
            yaml.safe_load(self.path.read_text(encoding="utf-8")) or {}
        )
        self._lines: dict[str, Any] = self.data.get("lines", {}) or {}
        # normalised key → canonical line name (+ aliases)
        self._index: dict[str, str] = {_norm(k): k for k in self._lines}
        for alias, target in (self.data.get("aliases", {}) or {}).items():
            if target in self._lines:
                self._index[_norm(alias)] = target

    # ── lookup ──────────────────────────────────────────────────────────────
    @classmethod
    def from_default(cls, root: Path | None = None) -> ProjectRegistry | None:
        """Build from the auto-discovered registry, or ``None``."""
        reg = find_default_registry(root)
        return cls(reg, root=root) if reg else None

    def lines(self) -> list[str]:
        """Sorted canonical line names."""
        return sorted(self._lines)

    def __contains__(self, name: str) -> bool:
        return _norm(name) in self._index

    def canonical(self, name: str) -> str | None:
        """Canonical line name for *name* (handles case / 'line ' / alias)."""
        return self._index.get(_norm(name))

    def find_line_in_text(self, text: str) -> str | None:
        """Return the first registry line mentioned in *text*, if any."""
        if not text:
            return None
        # direct: any known line token appears
        for tok in _LINE_TOKEN.findall(text):
            canon = self.canonical(tok)
            if canon:
                return canon
        # fallback: whole-text normalised contains a key
        norm_blob = _norm(text)
        for key_norm, canon in self._index.items():
            if key_norm and key_norm in norm_blob:
                return canon
        return None

    def resolve_line(self, line_name: str) -> dict[str, Any]:
        """Resolve *line_name* to a path + parameters dict.

        Raises
        ------
        KeyError
            If the line is unknown (message lists the known lines).
        """
        canon = self.canonical(line_name)
        if canon is None:
            known = ", ".join(self.lines()) or "(none)"
            raise KeyError(
                f"Unknown survey line {line_name!r}. Known lines: {known}"
            )

        line = dict(self._lines[canon])
        project = self.data.get("project", {}) or {}
        out_root = project.get("default_output_root", "results")

        edi_dir_raw = line.get("edi_dir", "")
        edi_path = Path(edi_dir_raw)
        if not edi_path.is_absolute():
            edi_path = (self.root / edi_dir_raw).resolve()
        exists = edi_path.is_dir()
        n_edi = len(list(edi_path.glob("*.edi"))) if exists else 0

        output_root = line.get("output_root", f"{out_root}/{canon}")

        return {
            "project": project.get("name"),
            "line": canon,
            "edi_dir": str(edi_path),
            "edi_dir_raw": edi_dir_raw,
            "exists": exists,
            "n_edi_files": n_edi,
            "station_pattern": line.get("station_pattern"),
            "sort_by": line.get("sort_by", "lon"),
            "output_root": output_root,
            "default_workflows": line.get("default_workflows", []),
            "static_shift": dict(
                self.data.get("static_shift_defaults", {}) or {}
            ),
            "plot": dict(self.data.get("plot_defaults", {}) or {}),
        }


def resolve_line(
    line_name: str,
    *,
    registry_path: str | Path | None = None,
    root: Path | None = None,
) -> dict[str, Any]:
    """Convenience: resolve *line_name* via a registry.

    Uses *registry_path* if given, else the auto-discovered default.
    Raises :class:`FileNotFoundError` when no registry is available and
    :class:`KeyError` when the line is unknown.
    """
    reg = (
        ProjectRegistry(registry_path, root=root)
        if registry_path is not None
        else ProjectRegistry.from_default(root)
    )
    if reg is None:
        raise FileNotFoundError(
            "No project registry found. Create "
            "projects/willy_project_registry.yml or pass registry_path."
        )
    return reg.resolve_line(line_name)
