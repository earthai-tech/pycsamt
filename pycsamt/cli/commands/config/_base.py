# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt.cli.commands.config._base
==================================

Shared helpers for ``pycsamt config``:

* TOML read / write  (``~/.pycsamt.toml``)
* Dot-notation key parsing   ``control.rho.view`` → (section, ``rho__view``)
* Per-section appliers        map TOML dicts to the relevant ``configure_*()``
* Per-section singleton getters  for reading current effective values
* ``load_all_config()``       called at CLI startup to hydrate all singletons

TOML layout
-----------
::

    [control]
    rho__view        = "log10"
    x__view          = "log10_period"
    phase__wrap      = true

    [style]
    preset           = "publication"
    multiline__mode  = "gradient"

    [plot]
    fmt              = "pdf"
    dpi              = 300
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import click

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

TOML_PATH: Path = Path.home() / ".pycsamt.toml"

_KNOWN_SECTIONS: frozenset[str] = frozenset(
    {
        "cli",
        "control",
        "style",
        "section_view",
        "station",
        "interp",
        "plot",
        "agent",
        "pipe",
        "view",
    }
)

# Human-readable descriptions (shown in config list / show)
_SECTION_LABELS: dict[str, str] = {
    "cli": "CLI behaviour (verbosity, color, output format)",
    "control": "Plot axis controls (rho scale, phase range, x-axis)",
    "style": "Visual style (MT components, multiline, PT ellipses)",
    "section_view": "Section-plot layout (figure size, colorbar, axis direction)",
    "station": "Station-tick rendering (density, rotation, marker)",
    "interp": "Hydro-interp colormaps and figure sizes",
    "plot": "Figure export (format, DPI, save directory)",
    "agent": "LLM provider configuration (model, budget)",
    "pipe": "Pipeline runtime settings (error mode, progress, output)",
    "view": "DataFrame backend (pycsamt / pandas)",
}

# Environment variables recognised by pycsamt (shown in config env)
_ENV_VARS: list[tuple[str, str, str]] = [
    # var,                          section,  description
    ("PYCSAMT_VERBOSE", "cli", "Verbosity level (0=quiet, 1=info, 2=debug)"),
    ("PYCSAMT_NO_COLOR", "cli", "Disable ANSI color (any non-empty value)"),
    ("PYCSAMT_OUTPUT", "cli", "Output format: text | json | csv"),
    ("PYCSAMT_OUTPUT_DIR", "cli", "Default output directory"),
    ("PYCSAMT_JOBS", "cli", "Parallel worker count"),
    ("PYCSAMT_FMT", "plot", "Export format: png | pdf | svg | jpg"),
    ("PYCSAMT_BASE_FMT", "plot", "Base format for additive tokens"),
    ("PYCSAMT_DPI", "plot", "Raster DPI"),
    ("PYCSAMT_SAVEDIR", "plot", "Default figure save directory"),
    (
        "PYCSAMT_TRANSPARENT",
        "plot",
        "Transparent background (any non-empty value)",
    ),
    ("PYCSAMT_BBOX", "plot", "bbox_inches value (e.g. 'tight')"),
    ("ANTHROPIC_API_KEY", "agent", "Claude / Anthropic API key"),
    ("PYCSAMT_CLAUDE_API_KEY", "agent", "Claude API key (pycsamt alias)"),
    ("OPENAI_API_KEY", "agent", "OpenAI API key"),
    ("PYCSAMT_OPENAI_API_KEY", "agent", "OpenAI API key (pycsamt alias)"),
    ("GOOGLE_API_KEY", "agent", "Google Gemini API key"),
    (
        "GOOGLE_GENERATIVEAI_API_KEY",
        "agent",
        "Google Generative AI API key (alias)",
    ),
    ("PYCSAMT_GEMINI_API_KEY", "agent", "Gemini API key (pycsamt alias)"),
    ("PYCSAMT_API_VIEW", "view", "DataFrame backend: pycsamt | pandas"),
]

# ---------------------------------------------------------------------------
# TOML I/O
# ---------------------------------------------------------------------------


def _read_toml() -> dict[str, Any]:
    """Load ``~/.pycsamt.toml``; return empty dict when absent or unreadable."""
    if not TOML_PATH.exists():
        return {}
    try:
        try:
            import tomllib  # Python 3.11+
        except ImportError:
            try:
                import tomli as tomllib  # type: ignore[no-redef]
            except ImportError:
                return {}
        with open(TOML_PATH, "rb") as fh:
            return tomllib.load(fh)
    except Exception:  # noqa: BLE001
        return {}


def _write_toml(data: dict[str, Any]) -> None:
    """Write *data* to ``~/.pycsamt.toml`` using tomli_w or manual formatting."""
    TOML_PATH.parent.mkdir(parents=True, exist_ok=True)
    try:
        import tomli_w

        with open(TOML_PATH, "wb") as fh:
            tomli_w.dump(data, fh)
    except ImportError:
        # Fallback: hand-write the TOML (handles simple types only)
        lines: list[str] = []
        for section, kv in data.items():
            if not isinstance(kv, dict):
                continue
            if not kv:
                continue
            lines.append(f"\n[{section}]")
            for k, v in kv.items():
                if isinstance(v, bool):
                    lines.append(f"{k} = {str(v).lower()}")
                elif isinstance(v, str):
                    lines.append(f'{k} = "{v}"')
                elif isinstance(v, (int, float)):
                    lines.append(f"{k} = {v}")
                else:
                    lines.append(f'{k} = "{v}"')
        with open(TOML_PATH, "w", encoding="utf-8") as fh:
            fh.write("\n".join(lines).lstrip() + "\n")


# ---------------------------------------------------------------------------
# Key parsing
# ---------------------------------------------------------------------------


def parse_key(key: str) -> tuple[str, str]:
    """Parse a dot-notation key into ``(section, toml_key)``.

    Examples
    --------
    ``"control.rho.view"``   → ``("control", "rho__view")``
    ``"style.preset"``       → ``("style",   "preset")``
    ``"style.mt.xy.color"``  → ``("style",   "mt__xy__color")``
    ``"pipe.output_root"``   → ``("pipe",    "output_root")``
    """
    parts = key.split(".")
    if len(parts) < 2:
        raise click.BadParameter(
            f"{key!r} must be at least <section>.<attribute>  "
            f"(e.g. plot.dpi, control.rho.view)",
            param_hint="KEY",
        )
    section = parts[0]
    if section not in _KNOWN_SECTIONS:
        raise click.BadParameter(
            f"Unknown section {section!r}.  "
            f"Valid sections: {', '.join(sorted(_KNOWN_SECTIONS))}",
            param_hint="KEY",
        )
    toml_key = "__".join(parts[1:])
    return section, toml_key


def toml_key_to_dotted(section: str, toml_key: str) -> str:
    """Reverse of ``parse_key``: ``("control", "rho__view")`` → ``"control.rho.view"``."""
    return section + "." + toml_key.replace("__", ".")


# ---------------------------------------------------------------------------
# Value coercion
# ---------------------------------------------------------------------------


def coerce(value: str) -> bool | int | float | str:
    """Coerce a CLI string value to an appropriate Python scalar."""
    if value.lower() in ("true", "yes", "on"):
        return True
    if value.lower() in ("false", "no", "off"):
        return False
    try:
        return int(value)
    except ValueError:
        pass
    try:
        return float(value)
    except ValueError:
        pass
    return value


# ---------------------------------------------------------------------------
# Per-section appliers   (TOML dict  →  configure_*() call)
# ---------------------------------------------------------------------------


def apply_section(section: str, kwargs: dict[str, Any]) -> None:
    """Apply a TOML section's key-value pairs to the matching API singleton.

    Parameters
    ----------
    section:
        One of the known section names.
    kwargs:
        Flat dict of ``{toml_key: value}`` pairs, exactly as read from TOML.
        Keys use ``__`` to separate sub-attributes (e.g. ``rho__view``).
    """
    if not kwargs:
        return

    kw = dict(kwargs)  # local copy so we can pop without affecting caller

    try:
        if section == "cli":
            from pycsamt.api.cli.config import (
                configure_cli,  # noqa: PLC0415
            )

            _apply_flat_or_nested(kw, configure_cli)

        elif section == "control":
            from pycsamt.api.control import (
                configure_control,  # noqa: PLC0415
            )

            configure_control(**kw)

        elif section == "style":
            preset = kw.pop("preset", None)
            if preset:
                from pycsamt.api.style import (
                    use_style,  # noqa: PLC0415
                )

                use_style(preset)
            if kw:
                from pycsamt.api.style import (
                    configure_style,  # noqa: PLC0415
                )

                configure_style(**kw)

        elif section == "section_view":
            from pycsamt.api.section import (
                configure_section,  # noqa: PLC0415
            )

            configure_section(**kw)

        elif section == "station":
            from pycsamt.api.station import (
                configure_station_rendering,  # noqa: PLC0415
            )

            configure_station_rendering(**kw)

        elif section == "interp":
            preset = kw.pop("preset", None)
            if preset:
                from pycsamt.api.interp import (
                    use_interp,  # noqa: PLC0415
                )

                use_interp(preset)
            if kw:
                from pycsamt.api.interp import (
                    configure_interp,  # noqa: PLC0415
                )

                configure_interp(**kw)

        elif section == "plot":
            from pycsamt.api.plot import (  # noqa: PLC0415
                PLOT_CONFIG,
                set_dpi,
                set_fmt,
                set_savedir,
            )

            fmt = kw.pop("fmt", None)
            dpi = kw.pop("dpi", None)
            savedir = kw.pop("savedir", None)
            base_fmt = kw.pop("base_fmt", None)
            if fmt is not None:
                set_fmt(str(fmt))
            if dpi is not None:
                set_dpi(int(dpi))
            if savedir is not None:
                set_savedir(str(savedir))
            if base_fmt is not None:
                PLOT_CONFIG.base_fmt = str(base_fmt)
            for k, v in kw.items():
                if hasattr(PLOT_CONFIG, k):
                    setattr(PLOT_CONFIG, k, v)

        elif section == "agent":
            from pycsamt.api.agents import (
                AGENT_CONFIG,  # noqa: PLC0415
            )

            provider = kw.pop("provider", None)
            model = kw.pop("model", None)
            budget_usd = kw.pop("budget_usd", None)
            if provider:
                try:
                    AGENT_CONFIG.switch(provider, model=model)
                except Exception:  # noqa: BLE001
                    pass
            elif model:
                AGENT_CONFIG.model = model
            if budget_usd is not None:
                AGENT_CONFIG.set_budget(usd=float(budget_usd))

        elif section == "pipe":
            from pycsamt.api.pipe.config import (
                configure_pipe,  # noqa: PLC0415
            )

            configure_pipe(**kw)

        elif section == "view":
            from pycsamt.api.view.config import (
                configure_api_view,  # noqa: PLC0415
            )

            configure_api_view(**kw)

    except Exception as exc:  # noqa: BLE001
        # Non-fatal at load time; report and continue
        click.echo(
            f"  Warning: could not apply [{section}] config: {exc}",
            err=True,
        )


def _apply_flat_or_nested(
    kw: dict[str, Any],
    configure_fn: Any,
) -> None:
    """Handle both flat (``log__level=0``) and nested (``{"log": {"level": 0}}``) formats."""
    for k, v in kw.items():
        if isinstance(v, dict):
            # Legacy nested TOML: [cli.log] level = 0
            for sub_k, sub_v in v.items():
                configure_fn(**{f"{k}__{sub_k}": sub_v})
        else:
            configure_fn(**{k: v})


# ---------------------------------------------------------------------------
# Apply all TOML sections  (called at CLI startup)
# ---------------------------------------------------------------------------


def load_all_config(toml_data: dict[str, Any] | None = None) -> None:
    """Apply every section in ``~/.pycsamt.toml`` to its singleton.

    Called once at CLI startup after loading the TOML file.
    ``apply_section`` is wrapped in a broad except so a bad config
    file never prevents the CLI from starting.
    """
    data = toml_data if toml_data is not None else _read_toml()
    for section in _KNOWN_SECTIONS:
        section_data = data.get(section)
        if not section_data or not isinstance(section_data, dict):
            continue
        apply_section(section, dict(section_data))


# ---------------------------------------------------------------------------
# Per-section singleton getters   (for config get / config show)
# ---------------------------------------------------------------------------

_SINGLETON_MAP: dict[str, tuple[str, str]] = {
    "cli": ("pycsamt.api.cli.config", "PYCSAMT_CLI"),
    "control": ("pycsamt.api.control", "PYCSAMT_CONTROL"),
    "style": ("pycsamt.api.style", "PYCSAMT_STYLE"),
    "section_view": ("pycsamt.api.section", "PYCSAMT_SECTION"),
    "station": ("pycsamt.api.station", "PYCSAMT_STATION_RENDERING"),
    "interp": ("pycsamt.api.interp", "PYCSAMT_INTERP"),
    "plot": ("pycsamt.api.plot", "PLOT_CONFIG"),
    "agent": ("pycsamt.api.agents", "AGENT_CONFIG"),
    "pipe": ("pycsamt.api.pipe.config", "PYCSAMT_PIPE"),
    "view": ("pycsamt.api.view.config", "PYCSAMT_API_VIEW"),
}


def get_singleton(section: str) -> Any:
    """Return the live singleton for *section*."""
    import importlib  # noqa: PLC0415

    module_path, attr = _SINGLETON_MAP[section]
    mod = importlib.import_module(module_path)
    return getattr(mod, attr)


def get_value(section: str, toml_key: str) -> Any:
    """Return the current effective value of *section*.*toml_key* from the singleton.

    The key ``rho__view`` is traversed as ``obj.rho.view`` on the singleton.
    """
    singleton = get_singleton(section)
    parts = toml_key.replace("__", ".").split(".")
    obj = singleton
    for part in parts:
        obj = getattr(obj, part)
    return obj


# ---------------------------------------------------------------------------
# Section summary helpers  (used by config show)
# ---------------------------------------------------------------------------


def section_summary(section: str) -> str:
    """Return a human-readable summary string for the given section."""
    try:
        singleton = get_singleton(section)
    except Exception:  # noqa: BLE001
        return "(unavailable)"

    # Most singletons expose summary()
    if hasattr(singleton, "summary"):
        try:
            return singleton.summary()
        except Exception:  # noqa: BLE001
            pass

    # Fallback: repr
    try:
        return repr(singleton)
    except Exception:  # noqa: BLE001
        return str(singleton)


__all__ = [
    "TOML_PATH",
    "_KNOWN_SECTIONS",
    "_SECTION_LABELS",
    "_ENV_VARS",
    "_SINGLETON_MAP",
    "parse_key",
    "toml_key_to_dotted",
    "coerce",
    "apply_section",
    "load_all_config",
    "get_singleton",
    "get_value",
    "section_summary",
    "_read_toml",
    "_write_toml",
]
