"""Config-file loaders for the pyCSAMT pipeline.

Supports three config formats:

YAML
    ``Pipeline.from_yaml("workflow.yaml")``

JSON
    ``Pipeline.from_json("workflow.json")``

Python module
    ``Pipeline.from_py("workflow.py")``

    The Python file must expose a module-level variable named
    ``pipeline_config`` that is a :class:`dict` with the same schema as the
    YAML/JSON format.

Config schema
-------------
::

    name: "My Workflow"        # optional pipeline label
    output_dir: "results/"     # default outdir when Pipeline.run() is called
    preset: "basic_qc"         # optional: seed steps from a preset before
                               # applying the steps list below
    steps:
      - name: notch            # user label for the step (optional)
        code: NR001            # step registry code OR name
        params:                # keyword arguments (optional)
          mains_hz: 50
          n_harm: 30
      - name: select_band
        code: FREQ001
        params:
          band_hz: [0.001, 10000]
      - code: FREQ004          # name is optional — auto-set from registry
      - code: SS001

The ``steps`` list is processed in order; if ``preset`` is specified, its
steps come *first* and the explicit ``steps`` entries are appended.
"""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path
from typing import Any

# ---------------------------------------------------------------------------
# Raw-dict → step list helper
# ---------------------------------------------------------------------------


def _steps_from_list(
    raw_steps: list[dict[str, Any]],
) -> list[tuple[str, Any]]:
    """Convert a list of step dicts to ``[(label, Step), ...]``."""
    from ._steps import Step

    result = []
    for entry in raw_steps:
        code = entry.get("code") or entry.get("name")
        if not code:
            raise ValueError(
                f"Each step entry must have a 'code' field.  Got: {entry!r}"
            )
        params = entry.get("params") or {}
        step = Step(code, **params)
        label = entry.get("name") or step.spec.name
        result.append((label, step))
    return result


def _pipeline_from_dict(cfg: dict[str, Any]) -> Any:
    """Build a :class:`~pycsamt.emtools.pipe.Pipeline` from a raw config dict.

    This function lives here (not in _pipeline.py) to avoid circular imports.
    It is called by the ``from_yaml`` / ``from_json`` / ``from_py``
    classmethods.
    """
    from ._presets import get_preset

    name = cfg.get("name", "unnamed")
    output_dir = cfg.get("output_dir", None)
    preset_name = cfg.get("preset", None)
    raw_steps = cfg.get("steps") or []

    steps: list[tuple[str, Any]] = []

    if preset_name:
        preset = get_preset(preset_name)
        steps.extend(preset.steps)

    steps.extend(_steps_from_list(raw_steps))

    return steps, name, output_dir


# ---------------------------------------------------------------------------
# Format-specific loaders
# ---------------------------------------------------------------------------


def load_yaml(path: str | Path) -> dict[str, Any]:
    """Parse a YAML pipeline config file and return the raw dict."""
    try:
        import yaml
    except ImportError as exc:
        raise ImportError(
            "PyYAML is required to load YAML pipeline configs.  "
            "Install it with:  pip install pyyaml"
        ) from exc

    with open(path, encoding="utf-8") as fh:
        data = yaml.safe_load(fh)
    if not isinstance(data, dict):
        raise ValueError(
            f"YAML config must be a mapping at the top level.  Got: {type(data)}"
        )
    return data


def load_json(path: str | Path) -> dict[str, Any]:
    """Parse a JSON pipeline config file and return the raw dict."""
    with open(path, encoding="utf-8") as fh:
        data = json.load(fh)
    if not isinstance(data, dict):
        raise ValueError(
            f"JSON config must be an object at the top level.  Got: {type(data)}"
        )
    return data


def load_py(path: str | Path) -> dict[str, Any]:
    """Import a Python config file and return its ``pipeline_config`` dict."""
    path = Path(path)
    spec = importlib.util.spec_from_file_location("_pipe_cfg", path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot import pipeline config from {path}")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)  # type: ignore[attr-defined]

    cfg = getattr(mod, "pipeline_config", None)
    if cfg is None:
        raise AttributeError(
            f"Python pipeline config file {path} must define a "
            "'pipeline_config' variable."
        )
    if not isinstance(cfg, dict):
        raise TypeError(
            f"'pipeline_config' in {path} must be a dict.  Got: {type(cfg)}"
        )
    return cfg


# ---------------------------------------------------------------------------
# Pipeline serialiser (YAML string output)
# ---------------------------------------------------------------------------


def _coerce_for_yaml(obj: Any) -> Any:
    """Recursively convert tuples → lists so yaml.safe_load can reload."""
    if isinstance(obj, dict):
        return {k: _coerce_for_yaml(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_coerce_for_yaml(v) for v in obj]
    return obj


def pipeline_to_yaml(
    steps: list[tuple[str, Any]],
    name: str = "pipeline",
    output_dir: str | None = None,
) -> str:
    """Serialise a list of ``(label, Step)`` tuples to a YAML string.

    Produces the canonical format that :func:`load_yaml` can reload.
    Falls back to a plain-text representation when PyYAML is unavailable.
    """
    step_list = []
    for label, step in steps:
        entry: dict[str, Any] = {"name": label, "code": step.spec.code}
        if step.params:
            entry["params"] = step.params
        step_list.append(entry)

    data: dict[str, Any] = {"name": name, "steps": step_list}
    if output_dir is not None:
        data["output_dir"] = output_dir

    try:
        import yaml

        return yaml.dump(
            _coerce_for_yaml(data),
            default_flow_style=False,
            sort_keys=False,
        )
    except ImportError:
        # Minimal fallback: hand-craft YAML without pyyaml
        lines = [f"name: {name!r}"]
        if output_dir:
            lines.append(f"output_dir: {output_dir!r}")
        lines.append("steps:")
        for entry in step_list:
            lines.append(f"  - name: {entry['name']!r}")
            lines.append(f"    code: {entry['code']!r}")
            if entry.get("params"):
                lines.append("    params:")
                for k, v in entry["params"].items():
                    lines.append(f"      {k}: {v!r}")
        return "\n".join(lines) + "\n"


__all__ = [
    "load_yaml",
    "load_json",
    "load_py",
    "pipeline_to_yaml",
    "_pipeline_from_dict",
]
