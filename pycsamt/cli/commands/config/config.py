# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt config — root Click group plus list / get / set / unset / reset / show / env.

Architecture
------------
* ``list [SECTION]``   — show all persisted keys (or one section) from TOML.
* ``get  KEY``         — print the current *effective* value of a key.
* ``set  KEY VALUE``   — write to TOML and apply immediately.
* ``unset KEY``        — remove a key from TOML.
* ``reset [SECTION]``  — remove a section (or clear the entire TOML).
* ``show [SECTION]``   — rich panel of the current effective configuration.
* ``env``              — table of environment variables recognised by pycsamt.
* ``style PRESET``     — shortcut to activate a named visual-style preset.
* ``interp PRESET``    — shortcut to activate a named interp preset.
* ``agent``            — agent sub-commands (status, set-key guidance).
"""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path
from typing import Any

import click

from . import _base as _b
from ._base import (
    _ENV_VARS,
    _KNOWN_SECTIONS,
    _SECTION_LABELS,
    apply_section,
    coerce,
    get_value,
    parse_key,
    section_summary,
    toml_key_to_dotted,
)

# Convenience accessors that always resolve through _b so monkeypatch works:
def _TOML_PATH():       return _b.TOML_PATH
def _read_toml():       return _b._read_toml()
def _write_toml(data):  return _b._write_toml(data)

# ---------------------------------------------------------------------------
# Optional rich formatting
# ---------------------------------------------------------------------------

try:
    from rich.console import Console
    from rich.panel import Panel
    from rich.table import Table
    _RICH = True
except ImportError:
    _RICH = False


# ---------------------------------------------------------------------------
# Root group
# ---------------------------------------------------------------------------

@click.group("config")
@click.pass_context
def config(ctx: click.Context) -> None:
    """Read and write pyCSAMT configuration.

    Settings are persisted in ``~/.pycsamt.toml`` and applied automatically
    at every CLI invocation.  Use ``pycsamt config show`` to view the full
    current configuration, or ``pycsamt config set KEY VALUE`` to change a
    setting permanently.

    \b
    Key naming
    ----------
    Keys use dot-notation: ``<section>.<attribute>``  (or deeper).
    Sections: cli, control, style, section_view, station, interp, plot,
              agent, pipe, view.

    \b
    Quick examples
    --------------
      pycsamt config set  plot.dpi 300
      pycsamt config set  plot.fmt pdf
      pycsamt config set  control.rho.view log10
      pycsamt config set  style.preset publication
      pycsamt config set  agent.provider claude
      pycsamt config set  agent.budget_usd 5.0
      pycsamt config get  plot.dpi
      pycsamt config list plot
      pycsamt config show
      pycsamt config env
    """
    ctx.ensure_object(dict)


# ---------------------------------------------------------------------------
# config list
# ---------------------------------------------------------------------------

@config.command("list")
@click.argument(
    "section",
    required=False,
    default=None,
    metavar="[SECTION]",
)
@click.option(
    "--format", "output_format",
    type=click.Choice(["text", "json"], case_sensitive=False),
    default="text",
    show_default=True,
)
@click.pass_context
def list_cmd(ctx: click.Context, section: str | None, output_format: str) -> None:
    """Show all persisted keys in ``~/.pycsamt.toml``.

    With no SECTION argument all sections are shown.  Pass a section name
    (e.g. ``plot``, ``style``) to show only that group.

    \b
    Examples:
      pycsamt config list
      pycsamt config list plot
      pycsamt config list agent --format json
    """
    data = _read_toml()

    if section is not None:
        section = section.lower()
        if section not in _KNOWN_SECTIONS:
            click.echo(
                f"Unknown section {section!r}.  "
                f"Valid: {', '.join(sorted(_KNOWN_SECTIONS))}",
                err=True,
            )
            sys.exit(1)
        sections_to_show = [section]
    else:
        sections_to_show = sorted(_KNOWN_SECTIONS)

    if output_format == "json":
        subset = {s: data.get(s, {}) for s in sections_to_show}
        click.echo(json.dumps(subset, indent=2))
        return

    if not data:
        click.echo(
            f"No configuration file found at {_TOML_PATH()}\n"
            "  Run 'pycsamt config set <key> <value>' to create one."
        )
        return

    if _RICH:
        console = Console()
        for s in sections_to_show:
            kv = data.get(s, {})
            if not kv:
                continue
            tbl = Table.grid(padding=(0, 2))
            tbl.add_column(style="cyan", no_wrap=True)
            tbl.add_column(style="green")
            for k, v in sorted(kv.items()):
                dotted = toml_key_to_dotted(s, k)
                tbl.add_row(dotted, str(v))
            console.print(Panel(
                tbl,
                title=f"[bold]{s}[/bold]  [dim]{_SECTION_LABELS.get(s, '')}[/dim]",
                border_style="blue",
                expand=False,
            ))
    else:
        for s in sections_to_show:
            kv = data.get(s, {})
            if not kv:
                continue
            click.echo(f"\n[{s}]")
            for k, v in sorted(kv.items()):
                click.echo(f"  {toml_key_to_dotted(s, k)} = {v!r}")


# ---------------------------------------------------------------------------
# config get
# ---------------------------------------------------------------------------

@config.command("get")
@click.argument("key")
@click.option(
    "--format", "output_format",
    type=click.Choice(["text", "json"], case_sensitive=False),
    default="text",
)
def get_cmd(key: str, output_format: str) -> None:
    """Print the current effective value of KEY.

    The returned value reflects env-var overrides and any runtime
    changes — not just what is persisted in the TOML file.

    \b
    Examples:
      pycsamt config get plot.dpi
      pycsamt config get control.rho.view
      pycsamt config get agent.provider
    """
    try:
        section, toml_key = parse_key(key)
    except click.BadParameter as exc:
        click.echo(str(exc), err=True)
        sys.exit(1)

    try:
        value = get_value(section, toml_key)
    except AttributeError:
        click.echo(
            f"Key {key!r} not found in the {section!r} configuration.",
            err=True,
        )
        sys.exit(1)
    except Exception as exc:  # noqa: BLE001
        click.echo(f"Could not read {key!r}: {exc}", err=True)
        sys.exit(1)

    if output_format == "json":
        click.echo(json.dumps({key: value}, default=str))
    else:
        click.echo(f"{key} = {value!r}")


# ---------------------------------------------------------------------------
# config set
# ---------------------------------------------------------------------------

@config.command("set")
@click.argument("key")
@click.argument("value")
@click.option(
    "--dry-run",
    is_flag=True,
    default=False,
    help="Validate the value without writing to disk.",
)
def set_cmd(key: str, value: str, dry_run: bool) -> None:
    """Set KEY to VALUE and persist it in ``~/.pycsamt.toml``.

    The value is applied immediately (to validate it) and written to TOML
    so that all future invocations inherit it.

    Value types are inferred automatically:
    ``true``/``false`` → bool, integer strings → int,
    decimal strings → float, everything else → str.

    \b
    Examples:
      pycsamt config set plot.dpi 300
      pycsamt config set plot.fmt pdf
      pycsamt config set control.rho.view log10
      pycsamt config set style.preset publication
      pycsamt config set agent.provider claude
      pycsamt config set agent.budget_usd 5.0
      pycsamt config set view.backend pandas
      pycsamt config set pipe.on_step_error warn
      pycsamt config set control.phase.wrap true
    """
    try:
        section, toml_key = parse_key(key)
    except click.BadParameter as exc:
        click.echo(str(exc), err=True)
        sys.exit(1)

    coerced = coerce(value)

    # Validate by applying immediately (raises on bad values)
    try:
        apply_section(section, {toml_key: coerced})
    except Exception as exc:  # noqa: BLE001
        click.echo(f"Invalid value for {key!r}: {exc}", err=True)
        sys.exit(1)

    if dry_run:
        click.echo(f"[dry-run] {key} = {coerced!r}  (not written)")
        return

    # Persist
    data = _read_toml()
    if section not in data:
        data[section] = {}
    data[section][toml_key] = coerced
    _write_toml(data)

    click.echo(f"Set  {key} = {coerced!r}")
    click.echo(f"     (written to {_TOML_PATH()})")


# ---------------------------------------------------------------------------
# config unset
# ---------------------------------------------------------------------------

@config.command("unset")
@click.argument("key")
def unset_cmd(key: str) -> None:
    """Remove KEY from ``~/.pycsamt.toml``.

    The runtime value is not changed (the singleton keeps its current
    value until the next CLI invocation, when it will revert to the
    package default or env-var override).

    \b
    Examples:
      pycsamt config unset plot.dpi
      pycsamt config unset style.preset
    """
    try:
        section, toml_key = parse_key(key)
    except click.BadParameter as exc:
        click.echo(str(exc), err=True)
        sys.exit(1)

    data = _read_toml()
    section_data = data.get(section, {})

    if toml_key not in section_data:
        click.echo(
            f"{key!r} is not set in {_TOML_PATH()}  (nothing to unset).",
            err=True,
        )
        return

    del data[section][toml_key]
    if not data[section]:
        del data[section]
    _write_toml(data)
    click.echo(f"Unset  {key}  (removed from {_TOML_PATH()})")


# ---------------------------------------------------------------------------
# config reset
# ---------------------------------------------------------------------------

@config.command("reset")
@click.argument("section", required=False, default=None, metavar="[SECTION]")
@click.option(
    "--yes", "-y",
    is_flag=True,
    default=False,
    help="Skip the confirmation prompt.",
)
def reset_cmd(section: str | None, yes: bool) -> None:
    """Remove a SECTION from TOML, or clear the entire config file.

    Without SECTION, resets *all* pycsamt configuration back to defaults.
    This does not reset env-variable overrides — only the TOML file is
    affected.

    \b
    Examples:
      pycsamt config reset plot
      pycsamt config reset style --yes
      pycsamt config reset          # clears ~/.pycsamt.toml entirely
    """
    if section is not None:
        section = section.lower()
        if section not in _KNOWN_SECTIONS:
            click.echo(
                f"Unknown section {section!r}.  "
                f"Valid: {', '.join(sorted(_KNOWN_SECTIONS))}",
                err=True,
            )
            sys.exit(1)
        target = f"[{section}] section of {_TOML_PATH()}"
    else:
        target = str(_TOML_PATH())

    if not yes:
        click.confirm(
            f"Reset {target}?",
            default=False,
            abort=True,
        )

    data = _read_toml()

    if section is not None:
        if section in data:
            del data[section]
            _write_toml(data)
            click.echo(f"Reset [{section}] in {_TOML_PATH()}")
        else:
            click.echo(f"[{section}] was not set — nothing to reset.")
    else:
        if _TOML_PATH().exists():
            _TOML_PATH().unlink()
        click.echo(f"Removed {_TOML_PATH()}  (all settings reset to defaults)")


# ---------------------------------------------------------------------------
# config show
# ---------------------------------------------------------------------------

@config.command("show")
@click.argument("section", required=False, default=None, metavar="[SECTION]")
@click.option(
    "--format", "output_format",
    type=click.Choice(["text", "json"], case_sensitive=False),
    default="text",
    show_default=True,
)
def show_cmd(section: str | None, output_format: str) -> None:
    """Show the current *effective* configuration (including env-var overrides).

    Unlike ``config list``, this shows the live values held by the API
    singletons — not just what is in the TOML file.

    \b
    Examples:
      pycsamt config show
      pycsamt config show plot
      pycsamt config show agent --format json
    """
    if section is not None:
        section = section.lower()
        if section not in _KNOWN_SECTIONS:
            click.echo(
                f"Unknown section {section!r}.  "
                f"Valid: {', '.join(sorted(_KNOWN_SECTIONS))}",
                err=True,
            )
            sys.exit(1)
        sections_to_show = [section]
    else:
        sections_to_show = sorted(_KNOWN_SECTIONS)

    if output_format == "json":
        out: dict[str, Any] = {}
        for s in sections_to_show:
            try:
                from ._base import get_singleton  # noqa: PLC0415
                singleton = get_singleton(s)
                if hasattr(singleton, "to_dict"):
                    out[s] = singleton.to_dict()
                elif hasattr(singleton, "__dict__"):
                    out[s] = {k: str(v) for k, v in vars(singleton).items()
                              if not k.startswith("_")}
                else:
                    out[s] = str(singleton)
            except Exception:  # noqa: BLE001
                out[s] = "(unavailable)"
        click.echo(json.dumps(out, indent=2, default=str))
        return

    if _RICH:
        console = Console()
        console.print()
        for s in sections_to_show:
            summary = _rich_section_summary(s)
            console.print(Panel(
                summary or "[dim](no detail available)[/dim]",
                title=(
                    f"[bold cyan]{s}[/bold cyan]  "
                    f"[dim]{_SECTION_LABELS.get(s, '')}[/dim]"
                ),
                border_style="bright_blue",
                expand=False,
            ))
        console.print()
    else:
        for s in sections_to_show:
            click.echo(f"\n=== {s} — {_SECTION_LABELS.get(s, '')} ===")
            click.echo(section_summary(s))


def _rich_section_summary(section: str) -> str:
    """Build a rich-markup string for one config section."""
    raw = section_summary(section)
    # Summary text is plain — just return it (rich renders it as text)
    return raw


# ---------------------------------------------------------------------------
# config env
# ---------------------------------------------------------------------------

@config.command("env")
@click.option("--section", "-s", default=None, help="Filter by section name.")
def env_cmd(section: str | None) -> None:
    """Show environment variables recognised by pycsamt.

    These override TOML settings and take effect for every invocation.
    Set them in your shell profile (``~/.bashrc``, ``~/.zshrc``, etc.) to
    make them permanent.

    \b
    Examples:
      pycsamt config env
      pycsamt config env --section plot
      pycsamt config env --section agent
    """
    rows = [
        (var, sec, desc)
        for var, sec, desc in _ENV_VARS
        if section is None or sec == section.lower()
    ]

    if _RICH:
        console = Console()
        tbl = Table(
            "Variable",
            "Section",
            "Description",
            "Current value",
            title="pyCSAMT environment variables",
            show_header=True,
            header_style="bold cyan",
            border_style="bright_blue",
        )
        for var, sec, desc in rows:
            value = os.environ.get(var, "")
            masked = _mask(var, value)
            tbl.add_row(
                f"[yellow]{var}[/yellow]",
                f"[dim]{sec}[/dim]",
                desc,
                f"[green]{masked}[/green]" if value else "[dim](not set)[/dim]",
            )
        console.print()
        console.print(tbl)
        console.print()
        console.print(
            "[dim]Tip: set in ~/.bashrc or ~/.zshrc — "
            "env vars always override ~/.pycsamt.toml[/dim]"
        )
        console.print()
    else:
        click.echo(f"\n{'Variable':<40}  {'Section':<12}  {'Set?':<6}  Description")
        click.echo("-" * 100)
        for var, sec, desc in rows:
            value = os.environ.get(var, "")
            masked = _mask(var, value) if value else ""
            is_set = f"[{masked[:8]}…]" if masked else "no"
            click.echo(f"  {var:<38}  {sec:<12}  {is_set:<6}  {desc}")


def _mask(var: str, value: str) -> str:
    """Mask API keys so only the first 8 chars are shown."""
    if "KEY" in var.upper() or "SECRET" in var.upper():
        return value[:8] + "…" if len(value) > 8 else "***"
    return value


# ---------------------------------------------------------------------------
# config style  (preset shortcut)
# ---------------------------------------------------------------------------

@config.command("style")
@click.argument(
    "preset",
    type=click.Choice(
        ["pycsamt", "publication", "dark", "modem"], case_sensitive=False
    ),
)
@click.option(
    "--no-persist",
    is_flag=True,
    default=False,
    help="Apply in this session only; do not write to ~/.pycsamt.toml.",
)
def style_cmd(preset: str, no_persist: bool) -> None:
    """Activate a named visual-style preset.

    Shortcut for ``pycsamt config set style.preset PRESET``.

    Available presets
    -----------------
    pycsamt      Default colours and markers.
    publication  Thin lines, no fill markers, B&W-friendly.
    dark         Bright colours on dark backgrounds.
    modem        Optimised for ModEM-like section overlays.

    \b
    Examples:
      pycsamt config style publication
      pycsamt config style dark --no-persist
    """
    from pycsamt.api.style import use_style  # noqa: PLC0415
    use_style(preset)

    if not no_persist:
        data = _read_toml()
        data.setdefault("style", {})["preset"] = preset
        _write_toml(data)
        click.echo(f"Style preset set to {preset!r}  (saved to {_TOML_PATH()})")
    else:
        click.echo(f"Style preset set to {preset!r}  (session only)")


# ---------------------------------------------------------------------------
# config interp  (preset shortcut)
# ---------------------------------------------------------------------------

@config.command("interp")
@click.argument(
    "preset",
    type=click.Choice(
        ["default", "publication", "dark", "accessible"], case_sensitive=False
    ),
)
@click.option(
    "--no-persist",
    is_flag=True,
    default=False,
    help="Apply in this session only; do not write to ~/.pycsamt.toml.",
)
def interp_cmd(preset: str, no_persist: bool) -> None:
    """Activate a named hydro-interp style preset.

    Shortcut for ``pycsamt config set interp.preset PRESET``.

    Available presets
    -----------------
    default       Standard colormaps and figure sizes.
    publication   Muted palette suitable for journal figures.
    dark          High-contrast on dark backgrounds.
    accessible    Colorblind-safe palette (Okabe-Ito / CARTO).

    \b
    Examples:
      pycsamt config interp accessible
      pycsamt config interp publication --no-persist
    """
    from pycsamt.api.interp import use_interp  # noqa: PLC0415
    use_interp(preset)

    if not no_persist:
        data = _read_toml()
        data.setdefault("interp", {})["preset"] = preset
        _write_toml(data)
        click.echo(f"Interp preset set to {preset!r}  (saved to {_TOML_PATH()})")
    else:
        click.echo(f"Interp preset set to {preset!r}  (session only)")


# ---------------------------------------------------------------------------
# config agent  (sub-group)
# ---------------------------------------------------------------------------

@config.group("agent")
def agent_group() -> None:
    """LLM provider configuration.

    API keys are never stored in ``~/.pycsamt.toml``.
    Use environment variables — run  ``pycsamt config env --section agent``
    for the full list.
    """


@agent_group.command("status")
@click.option(
    "--format", "output_format",
    type=click.Choice(["text", "json"], case_sensitive=False),
    default="text",
)
def agent_status(output_format: str) -> None:
    """Show the current LLM agent configuration and API-key status.

    \b
    Examples:
      pycsamt config agent status
      pycsamt config agent status --format json
    """
    from pycsamt.api.agents import AGENT_CONFIG  # noqa: PLC0415

    rows = {
        "provider":      AGENT_CONFIG.provider,
        "model":         AGENT_CONFIG.model,
        "is_configured": AGENT_CONFIG.is_configured,
        "spent_usd":     AGENT_CONFIG.spent_usd,
        "remaining_usd": AGENT_CONFIG.remaining_usd,
    }

    if output_format == "json":
        click.echo(json.dumps(rows, default=str))
        return

    if _RICH:
        console = Console()
        tbl = Table.grid(padding=(0, 2))
        tbl.add_column(style="cyan", no_wrap=True)
        tbl.add_column(style="green")
        for k, v in rows.items():
            tbl.add_row(k, str(v) if v is not None else "[dim](not set)[/dim]")
        console.print(Panel(
            tbl,
            title="[bold]LLM agent status[/bold]",
            border_style="bright_blue",
            expand=False,
        ))
    else:
        for k, v in rows.items():
            click.echo(f"  {k:<20} {v}")


@agent_group.command("set-key")
@click.argument(
    "provider",
    type=click.Choice(["claude", "openai", "gemini"], case_sensitive=False),
)
def agent_set_key(provider: str) -> None:
    """Show instructions for setting the API key for PROVIDER.

    pycsamt never stores API keys in ``~/.pycsamt.toml`` — always use the
    environment variable shown below.

    \b
    Examples:
      pycsamt config agent set-key claude
      pycsamt config agent set-key openai
      pycsamt config agent set-key gemini
    """
    instructions: dict[str, tuple[str, str]] = {
        "claude": (
            "ANTHROPIC_API_KEY",
            "https://console.anthropic.com/  →  API Keys",
        ),
        "openai": (
            "OPENAI_API_KEY",
            "https://platform.openai.com/  →  API keys",
        ),
        "gemini": (
            "GOOGLE_API_KEY",
            "https://aistudio.google.com/  →  Get API key",
        ),
    }
    env_var, url = instructions[provider.lower()]

    if _RICH:
        console = Console()
        console.print()
        console.print(Panel(
            f"[bold]Step 1[/bold]  Obtain your key from:\n"
            f"  [link]{url}[/link]\n\n"
            f"[bold]Step 2[/bold]  Export the key in your shell profile:\n"
            f"  [yellow]export {env_var}=<your-key>[/yellow]\n\n"
            f"[bold]Step 3[/bold]  Activate the provider:\n"
            f"  [cyan]pycsamt config set agent.provider {provider}[/cyan]",
            title=f"[bold]Setting up {provider.title()} API key[/bold]",
            border_style="green",
            expand=False,
        ))
        console.print()
    else:
        click.echo(f"\nSetting up {provider.title()} API key")
        click.echo(f"  1. Obtain your key: {url}")
        click.echo(f"  2. Add to your shell profile:")
        click.echo(f"       export {env_var}=<your-key>")
        click.echo(f"  3. Activate the provider:")
        click.echo(f"       pycsamt config set agent.provider {provider}")
        click.echo()
