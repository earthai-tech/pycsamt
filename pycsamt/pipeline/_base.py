"""Base class for pyCSAMT pipeline objects.

:class:`PipelineBase` is the common ancestor for :class:`~.Pipeline` (and any
future pipeline variant).  It provides three groups of functionality that are
domain-specific to the pipeline subsystem and would not fit naturally into
:class:`~pycsamt.api.property.PyCSAMTObject` (which targets lightweight
config/style objects).

Why not inherit PyCSAMTObject?
-------------------------------
``PyCSAMTObject`` gives generic ``__repr__``, ``to_dict``, and ``summary``
built for dataclass-style config holders.  ``Pipeline`` needs a specialized
scikit-learn-style repr and its own serialization contract (``to_yaml``,
``to_json``, ``to_py``).  Inheriting from two incompatible repr contracts
creates confusion; a clean separation is better.

Groups provided
---------------
1. **Registry access** (static / class methods)
   Explore what steps are available without constructing a pipeline.

2. **Export** (instance methods)
   ``to_py(path)`` — new format that produces an editable Python script
   (YAML and JSON are on :class:`~.Pipeline` itself).

3. **Scaffold** (class method)
   ``scaffold(path, fmt, preset, name)`` — generate a ready-to-edit
   starter config in YAML, JSON, or Python.
"""

from __future__ import annotations

import json
import textwrap
from pathlib import Path
from typing import Any

from ._registry import (
    STEP_REGISTRY,
    categories,
    list_steps,
    lookup_step,
)

# ── formatting helpers ────────────────────────────────────────────────────────

_SEP  = "─" * 68
_SEP2 = "═" * 68


def _defaults_str(defaults: dict) -> str:
    if not defaults:
        return ""
    return "  " + "  ".join(f"{k}={v!r}" for k, v in defaults.items())


def _params_repr(params: dict) -> str:
    """Python repr for a params dict, suitable for embedding in source."""
    if not params:
        return ""
    inner = ", ".join(f"{k}={v!r}" for k, v in params.items())
    return f"dict({inner})"


# ── base class ────────────────────────────────────────────────────────────────

class PipelineBase:
    """Base class for pyCSAMT pipeline objects.

    Inherit from this to gain registry introspection, ``to_py`` export, and
    ``scaffold`` template generation.  :class:`~.Pipeline` inherits this.
    """

    # =========================================================================
    # 1. Registry access — static / class methods
    # =========================================================================

    @staticmethod
    def available_steps(category: str | None = None) -> list:
        """Return all registered :class:`~._registry.StepSpec` objects.

        Parameters
        ----------
        category:
            When supplied, restrict to steps in that category
            (e.g. ``"noise_removal"``, ``"frequency"``).

        Examples
        --------
        >>> Pipeline.available_steps()
        >>> Pipeline.available_steps("static_shift")
        """
        return list_steps(category)

    @staticmethod
    def available_categories() -> list[str]:
        """Return a sorted list of step category names.

        Examples
        --------
        >>> Pipeline.available_categories()
        ['dimensionality', 'frequency', 'noise_removal', ...]
        """
        return categories()

    @staticmethod
    def step_info(code_or_name: str) -> str:
        """Return a formatted info block for a single step.

        Parameters
        ----------
        code_or_name:
            Registry code (``"NR001"``) or snake-case name
            (``"notch_powerline"``).

        Examples
        --------
        >>> print(Pipeline.step_info("NR001"))
        >>> print(Pipeline.step_info("correct_ss_ama"))
        """
        spec = lookup_step(code_or_name)
        qc_names = [fn_name for _, fn_name in spec.qc_defs] if spec.qc_defs else []
        defs_str = (
            "  ".join(f"{k}={v!r}" for k, v in spec.defaults.items())
            if spec.defaults
            else "—"
        )
        qc_str = ", ".join(qc_names) if qc_names else "—"
        fn_path = (
            f"{spec.mod}.{spec.fn_name}" if spec.mod and spec.fn_name
            else "(built-in wrapper)"
        )
        lines = [
            _SEP,
            f"  {spec.code}  {spec.label}  [{spec.category}]",
            _SEP,
            f"  name     : {spec.name}",
            f"  function : {fn_path}",
            f"  defaults : {defs_str}",
            f"  qc plots : {qc_str}",
            f"  returns  : {'Sites (transform)' if spec.returns_sites else 'Sites unchanged (diagnostic)'}",
            "",
            "  Usage:",
            f'    Step("{spec.code}")                    # defaults',
        ]
        if spec.defaults:
            first_k, first_v = next(iter(spec.defaults.items()))
            lines.append(f'    Step("{spec.code}", {first_k}={first_v!r})   # override')
        lines.append(f'    Step("{spec.name}")              # by name')
        lines.append(_SEP)
        return "\n".join(lines)

    @classmethod
    def catalogue(cls, category: str | None = None) -> str:
        """Return a full formatted catalogue of all available steps.

        Parameters
        ----------
        category:
            Restrict output to one category, or ``None`` for all.

        Examples
        --------
        >>> print(Pipeline.catalogue())
        >>> print(Pipeline.catalogue("tensor"))
        """
        cats = [category] if category else categories()
        lines = []

        total = len(STEP_REGISTRY)
        n_cats = len(categories())
        header = (
            f"  {total} steps across {n_cats} categories"
            if category is None
            else f"  Category: {category}"
        )
        lines += [_SEP2, "  Available Pipeline Steps", header, _SEP2, ""]

        col_code  = 8
        col_name  = 22
        col_label = 34

        for cat in cats:
            specs = list_steps(cat)
            if not specs:
                continue
            lines.append(f"  {cat.upper()}  ─ {len(specs)} step{'s' if len(specs) != 1 else ''}")
            lines.append("  " + "─" * 64)
            for sp in specs:
                dstr = _defaults_str(sp.defaults)
                diag = "  ★ diagnostic" if not sp.returns_sites else ""
                lines.append(
                    f"  {sp.code:<{col_code}}  "
                    f"{sp.name:<{col_name}}  "
                    f"{sp.label:<{col_label}}"
                    f"{dstr}{diag}"
                )
            lines.append("")

        return "\n".join(lines)

    # =========================================================================
    # 2. Export — to_py (instance method; to_yaml / to_json stay on Pipeline)
    # =========================================================================

    def to_py(self, path: str | Path | None = None) -> str:
        """Serialize this pipeline to a Python config script.

        Produces a ``pipeline_config`` dict in a clean, readable ``.py``
        file that can be edited and reloaded with
        :meth:`~.Pipeline.from_py`.  Steps are grouped by category and
        commented with their human-readable labels.

        Parameters
        ----------
        path:
            If given, write the output to this file path.

        Returns
        -------
        str
            The generated Python source as a string (always returned,
            whether or not *path* is provided).

        Examples
        --------
        >>> src = pipe.to_py()
        >>> print(src)
        >>> pipe.to_py("config/my_workflow.py")
        """
        name = getattr(self, "name", "pipeline")
        outdir = getattr(self, "_output_dir", None) or "pipe_results"
        steps = list(getattr(self, "_steps", []))

        # Group steps by category for comments
        lines = [
            f"# pyCSAMT Pipeline: {name}",
            "# Generated by pycsamt v2 — Pipeline.to_py()",
            "# Reload : pipe = Pipeline.from_py(__file__)",
            "# Run    : result = pipe.run(sites)",
            "",
            "pipeline_config = dict(",
            f'    name       = "{name}",',
            f'    output_dir = "{outdir}",',
            "    steps      = [",
        ]

        current_cat: str | None = None
        for label, step in steps:
            cat = step.spec.category
            if cat != current_cat:
                lines.append(f"        # ── {cat} {'─' * max(0, 56 - len(cat))}")
                current_cat = cat
            params_part = (
                f", params={_params_repr(step.params)}"
                if step.params
                else ""
            )
            lines.append(
                f'        dict(name={label!r:<22}, code="{step.spec.code}"'
                f'{params_part}),'
            )

        lines += [
            "    ],",
            ")",
        ]

        src = "\n".join(lines) + "\n"

        if path is not None:
            p = Path(path)
            p.parent.mkdir(parents=True, exist_ok=True)
            p.write_text(src, encoding="utf-8")

        return src

    # =========================================================================
    # 3. Scaffold — generate an editable starter config
    # =========================================================================

    @classmethod
    def scaffold(
        cls,
        path: str | Path | None = None,
        *,
        fmt: str = "yaml",
        preset: str | None = None,
        name: str = "my_workflow",
        outdir: str = "pipe_results",
    ) -> str:
        """Generate a ready-to-edit starter pipeline config.

        The scaffold contains:
        - Active steps from *preset* (or a sensible default set)
        - All other available steps commented out with descriptions

        Parameters
        ----------
        path:
            Write the output to this file.  Extension is inferred from *fmt*
            if the path has none.
        fmt:
            Output format: ``"yaml"`` (default), ``"json"``, or ``"py"``.
        preset:
            Name of a built-in preset to use as the active step set
            (e.g. ``"basic_qc"``, ``"full_processing"``).  When ``None``
            a minimal ``basic_qc``-style set is used.
        name:
            Pipeline name written into the config.
        outdir:
            Default output directory written into the config.

        Returns
        -------
        str
            The generated config content as a string.

        Examples
        --------
        >>> print(Pipeline.scaffold())
        >>> Pipeline.scaffold("config/starter.yaml", preset="full_processing")
        >>> Pipeline.scaffold("config/starter.py",   fmt="py")
        """
        fmt = fmt.lower()
        if fmt not in {"yaml", "json", "py"}:
            raise ValueError(f"fmt must be 'yaml', 'json', or 'py'.  Got {fmt!r}")

        # Resolve active steps from preset or default
        active_steps = _resolve_scaffold_steps(preset)

        if fmt == "yaml":
            content = _scaffold_yaml(name, outdir, active_steps)
        elif fmt == "json":
            content = _scaffold_json(name, outdir, active_steps)
        else:
            content = _scaffold_py(name, outdir, active_steps, preset)

        if path is not None:
            p = Path(path)
            if not p.suffix:
                p = p.with_suffix(f".{fmt}")
            p.parent.mkdir(parents=True, exist_ok=True)
            p.write_text(content, encoding="utf-8")

        return content


# ── internal scaffold helpers ─────────────────────────────────────────────────

def _resolve_scaffold_steps(preset_name: str | None) -> list:
    """Return (label, StepSpec, params) triples for the scaffold active steps."""
    if preset_name is not None:
        from ._presets import get_preset
        preset = get_preset(preset_name)
        return [(lbl, step.spec, step.params) for lbl, step in preset.steps]
    # minimal default: FREQ002 → FREQ001 → FREQ004 → NR001 → SS001 → TZ001 → QC001
    defaults = [
        ("drop_dup",    "FREQ002", {}),
        ("select_band", "FREQ001", {"band_hz": (1e-3, 1e4)}),
        ("align",       "FREQ004", {}),
        ("notch",       "NR001",   {"mains_hz": 50}),
        ("correct_ss",  "SS001",   {}),
        ("rotate",      "TZ001",   {"method": "swift"}),
        ("qc",          "QC001",   {}),
    ]
    return [
        (lbl, lookup_step(code), params)
        for lbl, code, params in defaults
    ]


def _step_to_dict(label: str, spec: Any, params: dict) -> dict:
    d: dict = {"name": label, "code": spec.code}
    if params:
        # convert tuples → lists for JSON/YAML safety
        safe_params = {
            k: list(v) if isinstance(v, tuple) else v
            for k, v in params.items()
        }
        d["params"] = safe_params
    return d


def _scaffold_yaml(name: str, outdir: str, active: list) -> str:
    """Generate the YAML scaffold."""
    try:
        import yaml as _yaml
    except ImportError:
        _yaml = None  # type: ignore[assignment]

    active_codes = {spec.code for _, spec, _ in active}

    # Build header
    header = textwrap.dedent(f"""\
        # pyCSAMT Pipeline Configuration
        # Generated by: Pipeline.scaffold("{name}.yaml")
        # ─────────────────────────────────────────────────────────────
        # Usage :
        #   pipe   = Pipeline.from_yaml("this_file.yaml")
        #   result = pipe.run(sites)
        # Help  : Pipeline.catalogue()   |   Pipeline.step_info("NR001")
        # ─────────────────────────────────────────────────────────────

        name: {name}
        output_dir: {outdir}

        steps:
        """)

    step_lines: list[str] = []
    for cat in categories():
        specs = list_steps(cat)
        step_lines.append(f"  # ── {cat} ({'─' * max(0, 56 - len(cat))})")
        for sp in specs:
            _step_to_dict("", sp, sp.defaults)
            entry = f"  - {{name: {sp.name}, code: {sp.code}"
            if sp.defaults:
                defs_yaml = ", ".join(
                    f"{k}: {list(v) if isinstance(v, tuple) else v!r}"
                    for k, v in sp.defaults.items()
                )
                entry += f", params: {{{defs_yaml}}}"
            entry += f"}}  # {sp.label}"

            if sp.code in active_codes:
                step_lines.append(entry)
            else:
                step_lines.append(f"  # {entry.lstrip()}")

        step_lines.append("")

    return header + "\n".join(step_lines)


def _scaffold_json(name: str, outdir: str, active: list) -> str:
    """Generate the JSON scaffold."""
    data = {
        "_comment": (
            f"pyCSAMT Pipeline — {name}. "
            "Load: Pipeline.from_json(path)"
        ),
        "name": name,
        "output_dir": outdir,
        "steps": [
            _step_to_dict(lbl, spec, params)
            for lbl, spec, params in active
        ],
        "_available_steps_hint": (
            "Run Pipeline.catalogue() to see all 33 available steps."
        ),
    }
    return json.dumps(data, indent=2)


def _scaffold_py(
    name: str,
    outdir: str,
    active: list,
    preset_name: str | None,
) -> str:
    """Generate the Python scaffold."""
    subtitle = f"Preset: {preset_name}" if preset_name else f"Pipeline: {name}"
    header = textwrap.dedent(f"""\
        # pyCSAMT Pipeline Configuration
        # {subtitle}
        # ─────────────────────────────────────────────────────────────
        # Reload : pipe = Pipeline.from_py(__file__)
        # Run    : result = pipe.run(sites)
        # Help   : Pipeline.catalogue()  |  Pipeline.step_info("NR001")
        # ─────────────────────────────────────────────────────────────

        """)

    active_codes = {spec.code for _, spec, _ in active}

    body_lines: list[str] = [
        "pipeline_config = dict(",
        f'    name       = "{name}",',
        f'    output_dir = "{outdir}",',
        "    steps      = [",
    ]

    current_cat: str | None = None
    for cat in categories():
        specs = list_steps(cat)
        for sp in specs:
            if cat != current_cat:
                body_lines.append(
                    f"        # ── {cat} {'─' * max(0, 52 - len(cat))}"
                )
                current_cat = cat

            use_params = sp.defaults
            params_part = (
                f", params={_params_repr(use_params)}"
                if use_params
                else ""
            )
            entry = (
                f'        dict(name={sp.name!r:<24}, code="{sp.code}"{params_part}),'
                f"  # {sp.label}"
            )
            if sp.code in active_codes:
                body_lines.append(entry)
            else:
                body_lines.append(f"        # {entry.strip()}")

    body_lines += [
        "    ],",
        ")",
    ]

    return header + "\n".join(body_lines) + "\n"
