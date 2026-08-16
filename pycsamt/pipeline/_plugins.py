"""Entry-point discovery for third-party pipeline steps.

A plugin package registers pipeline steps by declaring a zero-argument
callable under the ``pycsamt.pipeline.steps`` entry-point group, e.g. in its
``pyproject.toml``::

    [project.entry-points."pycsamt.pipeline.steps"]
    my_plugin = "my_package.pipeline_steps:register"

where ``my_package.pipeline_steps.register()`` calls
:func:`pycsamt.pipeline.register_step` for whatever it contributes.

Discovery is **never** run automatically on ``import pycsamt.pipeline`` — no
third-party code should execute just because a package happens to be
installed.  Call :func:`discover_plugins` explicitly, or use the
``pycsamt pipe`` CLI, which calls it for you before every command.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass
from importlib import metadata as _importlib_metadata

ENTRY_POINT_GROUP = "pycsamt.pipeline.steps"


@dataclass
class PluginLoadResult:
    """Outcome of loading one ``pycsamt.pipeline.steps`` entry point.

    Attributes
    ----------
    name:
        The entry point's registered name (left-hand side in
        ``[project.entry-points."pycsamt.pipeline.steps"]``).
    ok:
        Whether the entry point loaded and ran without raising.
    error:
        ``str(exception)`` when ``ok`` is ``False``, else ``None``.
    """

    name: str
    ok: bool
    error: str | None = None


def _iter_entry_points(group: str):
    try:
        return _importlib_metadata.entry_points(group=group)
    except TypeError:
        # Python 3.9: entry_points() takes no arguments, returns a
        # dict[str, list[EntryPoint]] keyed by group name.
        return _importlib_metadata.entry_points().get(group, [])


def discover_plugins(
    *, group: str = ENTRY_POINT_GROUP, on_error: str = "warn"
) -> list[PluginLoadResult]:
    """Load every ``pycsamt.pipeline.steps`` entry point and run it.

    Each entry point must resolve to a zero-argument callable; it is
    expected to call :func:`pycsamt.pipeline.register_step` itself for
    whatever steps it contributes.  Never called automatically — see the
    module docstring.

    Parameters
    ----------
    group:
        Entry-point group to scan. Defaults to ``ENTRY_POINT_GROUP``.
    on_error:
        ``"warn"`` (default): a plugin that fails to load or raises is
        reported as a failed :class:`PluginLoadResult` and a
        :class:`UserWarning`; discovery continues with the remaining
        plugins. ``"raise"``: the first failure propagates immediately.

    Returns
    -------
    list[PluginLoadResult]
        One entry per discovered entry point, in discovery order.
    """
    if on_error not in ("warn", "raise"):
        raise ValueError(f"on_error must be 'warn' or 'raise', got {on_error!r}")

    results: list[PluginLoadResult] = []
    for ep in _iter_entry_points(group):
        try:
            register_fn = ep.load()
            register_fn()
        except Exception as exc:  # noqa: BLE001 - one bad plugin must not sink the rest
            if on_error == "raise":
                raise
            warnings.warn(
                f"pyCSAMT pipeline plugin {ep.name!r} failed to load: {exc}",
                stacklevel=2,
            )
            results.append(PluginLoadResult(name=ep.name, ok=False, error=str(exc)))
        else:
            results.append(PluginLoadResult(name=ep.name, ok=True))
    return results


__all__ = ["ENTRY_POINT_GROUP", "PluginLoadResult", "discover_plugins"]
