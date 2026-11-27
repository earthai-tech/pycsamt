"""Plugin interfaces for the pyCSAMT agent loop.

The agent loop (router, orchestrator, provenance) depends on three
domain-specific pieces. This module names them as
:class:`typing.Protocol` classes so another package can supply its own
domain without forking :mod:`pycsamt.agents`:

* :class:`WorkflowRegistry` - which workflows exist, how free text is
  classified into one, and which agent steps run for it;
* :class:`PlanValidator` - which plans are blocked, which need human
  review, and the deterministic route cross-check;
* :class:`RecordWriter` - the run record (events, snapshot,
  persistence) and the retrieval-corpus manifest hook.

A :class:`DomainPlugin` bundles one of each. Plugins are registered
with :func:`register_plugin` or advertised through the entry-point
group ``pycsamt.agent_plugins`` and found by :func:`discover_plugins`.
:func:`get_plugin` returns one by name; with no name it returns the
built-in ``"pycsamt"`` plugin, whose adapters delegate to the existing
tables and functions, so the default agent behaviour is unchanged.

The module itself uses the standard library only; the built-in
adapters import :mod:`pycsamt.agents` internals lazily, at call time,
so that patching those module attributes (as the test-suite does)
keeps working.

Writing a plugin
----------------
Implement the three protocols with plain classes (no inheritance is
needed; the protocols are structural), bundle them, and register:

>>> from pycsamt.agents.plugins import (
...     DomainPlugin, get_plugin, register_plugin)
>>> base = get_plugin()
>>> demo = DomainPlugin(
...     name="demo", registry=base.registry,
...     validator=base.validator,
...     record_factory=base.record_factory)
>>> register_plugin(demo, replace=True) is demo
True
>>> get_plugin("demo").name
'demo'

To ship it from another distribution, declare an entry point that
resolves to a :class:`DomainPlugin` instance or to a zero-argument
callable returning one::

    [project.entry-points."pycsamt.agent_plugins"]
    mtpy = "my_package.agent_plugin:make_plugin"

See ``examples/agent_plugins/mtpy_example`` in the repository for a
complete plugin.
"""
from __future__ import annotations

import logging
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass
from importlib import metadata
from pathlib import Path
from typing import Any, Protocol, runtime_checkable

__all__ = [
    "ENTRY_POINT_GROUP",
    "DEFAULT_PLUGIN",
    "WorkflowRegistry",
    "PlanValidator",
    "RecordWriter",
    "DomainPlugin",
    "register_plugin",
    "get_plugin",
    "available_plugins",
    "discover_plugins",
]

_log = logging.getLogger(__name__)

#: Entry-point group scanned by :func:`discover_plugins`.
ENTRY_POINT_GROUP = "pycsamt.agent_plugins"
#: Name of the built-in plugin returned by ``get_plugin()``.
DEFAULT_PLUGIN = "pycsamt"

#: One dispatch step: ``(step_name, agent_class_name, input_fn,
#: description[, required])``, the tuple shape the orchestrator runs.
Step = tuple


# ── protocols ─────────────────────────────────────────────────────


@runtime_checkable
class WorkflowRegistry(Protocol):
    """Catalogue of workflows a domain offers to the agent loop.

    Attributes
    ----------
    name : str
        Short identifier of the registry (normally the plugin name).
    descriptions : Mapping[str, str]
        Every known workflow id mapped to a one-line description.
        Ids absent from this mapping are unknown and get blocked.
    valid_workflows : frozenset of str
        Ids that go through the plan gate (validation and review).
        Known ids outside this set are dispatched on other paths
        (plots, tools) without a plan.
    dispatchable : frozenset of str
        Ids that have a step chain, i.e. for which
        :meth:`steps` returns a sequence.

    Notes
    -----
    Implementations must be deterministic: the same text always
    classifies to the same id, because the cross-check and the
    paper's offline router rows rely on it.
    """

    name: str

    @property
    def descriptions(self) -> Mapping[str, str]: ...

    @property
    def valid_workflows(self) -> frozenset[str]: ...

    @property
    def dispatchable(self) -> frozenset[str]: ...

    def classify(
        self, text: str, *, default: str | None = None
    ) -> str | None:
        """Return the workflow id *text* asks for, or *default*."""
        ...

    def family(self, workflow: str | None) -> str | None:
        """Return the analysis family of *workflow* (``None`` in)."""
        ...

    def normalise(self, workflow: str) -> str:
        """Return the canonical id for an alias or spelling."""
        ...

    def steps(self, workflow: str) -> Sequence[Step] | None:
        """Return the dispatch steps of *workflow*, or ``None``."""
        ...


@runtime_checkable
class PlanValidator(Protocol):
    """Gate that decides whether a plan may run.

    The orchestrator calls :meth:`errors` first (any error blocks
    execution), then :meth:`needs_review` (``True`` requires the
    caller to approve the exact plan hash).
    """

    def errors(self, plan: Any) -> list[str]:
        """Return human-readable validation errors (empty = valid)."""
        ...

    def needs_review(self, plan: Any) -> bool:
        """Return ``True`` when the plan needs human approval."""
        ...

    def cross_check(
        self, text: str, intent: str, workflow: str | None
    ) -> str | None:
        """Router-level check of an LLM route against a keyword route.

        Returns a reason string when the two routes disagree, else
        ``None``.
        """
        ...

    def route_flags(self, text: str, workflow: str) -> list[str]:
        """Orchestrator-level cross-check, as plan risk flags.

        Any returned flag is added to the plan's ``risk_flags`` and
        therefore sends the plan to the review gate.
        """
        ...


@runtime_checkable
class RecordWriter(Protocol):
    """Append-only run record written next to the outputs."""

    def emit(self, event: str, **fields: Any) -> None:
        """Append one event with JSON-serialisable *fields*."""
        ...

    def snapshot(self) -> dict[str, Any]:
        """Return the record as a redacted, JSON-ready dictionary."""
        ...

    def persist(self, path: str | Path) -> None:
        """Write the snapshot to *path* (and keep it updated)."""
        ...

    def corpus_manifest(self) -> Any:
        """Return the manifest of the retrieval index, or ``None``."""
        ...


@dataclass(frozen=True)
class DomainPlugin:
    """Bundle of the three pieces one domain supplies.

    Parameters
    ----------
    name : str
        Registry key, e.g. ``"pycsamt"`` or ``"mtpy_example"``.
    registry : WorkflowRegistry
        Workflow catalogue and dispatch steps.
    validator : PlanValidator
        Plan gate and cross-check.
    record_factory : callable
        ``record_factory(payload) -> RecordWriter``; *payload* is the
        request dictionary the run started from.
    description : str, default ""
        One-line summary shown by :func:`available_plugins`.

    Raises
    ------
    TypeError
        If a member does not satisfy its protocol.
    """

    name: str
    registry: WorkflowRegistry
    validator: PlanValidator
    record_factory: Callable[[dict], RecordWriter]
    description: str = ""

    def __post_init__(self) -> None:
        if not isinstance(self.name, str) or not self.name:
            raise TypeError("DomainPlugin.name must be a non-empty "
                            "str.")
        if not isinstance(self.registry, WorkflowRegistry):
            raise TypeError("registry does not satisfy "
                            "WorkflowRegistry.")
        if not isinstance(self.validator, PlanValidator):
            raise TypeError("validator does not satisfy "
                            "PlanValidator.")
        if not callable(self.record_factory):
            raise TypeError("record_factory must be callable.")


# ── built-in pyCSAMT adapters (delegate, never re-implement) ──────


class _PycsamtRegistry:
    """Adapter over :mod:`._workflows` and orchestrator steps."""

    name = DEFAULT_PLUGIN

    @property
    def descriptions(self) -> Mapping[str, str]:
        from . import _workflows

        return _workflows.WORKFLOW_DESCRIPTIONS

    @property
    def valid_workflows(self) -> frozenset[str]:
        from . import _workflow_plan

        return _workflow_plan.VALID_WORKFLOWS

    @property
    def dispatchable(self) -> frozenset[str]:
        from . import orchestrator

        return frozenset(orchestrator._WORKFLOW_STEPS)

    def classify(
        self, text: str, *, default: str | None = None
    ) -> str | None:
        from . import _workflows

        return _workflows.classify_workflow(text, default=default)

    def family(self, workflow: str | None) -> str | None:
        from . import _workflows

        return _workflows.workflow_family(workflow)

    def normalise(self, workflow: str) -> str:
        from . import _workflows

        return _workflows.normalise_workflow(workflow)

    def steps(self, workflow: str) -> Sequence[Step] | None:
        from . import orchestrator

        return orchestrator._WORKFLOW_STEPS.get(workflow)


class _PycsamtValidator:
    """Adapter over :mod:`._workflow_plan` and the route checks."""

    def errors(self, plan: Any) -> list[str]:
        from . import _workflow_plan

        return _workflow_plan.validate_workflow_plan(plan)

    def needs_review(self, plan: Any) -> bool:
        return plan.requires_human_review

    def cross_check(
        self, text: str, intent: str, workflow: str | None
    ) -> str | None:
        from . import router

        return router.route_disagreement(text, intent, workflow)

    def route_flags(self, text: str, workflow: str) -> list[str]:
        from . import orchestrator

        return orchestrator._route_cross_check(text, workflow)


def _pycsamt_record(payload: dict) -> RecordWriter:
    """Return a :class:`~._provenance.RunRecord` with the hook."""
    from ._provenance import RunRecord

    class PycsamtRunRecord(RunRecord):
        def corpus_manifest(self) -> Any:
            from ..assistant.rag.index_store import read_manifest

            return read_manifest()

    return PycsamtRunRecord(payload)


# ── plugin registry ───────────────────────────────────────────────

_PLUGINS: dict[str, DomainPlugin] = {}
_DISCOVERED = False


def _ensure_builtin() -> None:
    if DEFAULT_PLUGIN not in _PLUGINS:
        _PLUGINS[DEFAULT_PLUGIN] = DomainPlugin(
            name=DEFAULT_PLUGIN,
            registry=_PycsamtRegistry(),
            validator=_PycsamtValidator(),
            record_factory=_pycsamt_record,
            description="Built-in pyCSAMT CSAMT/MT workflows.",
        )


def register_plugin(
    plugin: DomainPlugin, *, replace: bool = False
) -> DomainPlugin:
    """Register *plugin* under ``plugin.name``.

    Parameters
    ----------
    plugin : DomainPlugin
        The bundle to register.
    replace : bool, default False
        Overwrite an existing plugin of the same name.

    Returns
    -------
    DomainPlugin
        *plugin*, so the call can be used as an expression.

    Raises
    ------
    TypeError
        If *plugin* is not a :class:`DomainPlugin`.
    ValueError
        If the name is taken and *replace* is false, or if the name
        is the built-in ``"pycsamt"``.

    Examples
    --------
    >>> from pycsamt.agents.plugins import get_plugin, register_plugin
    >>> register_plugin(get_plugin())
    Traceback (most recent call last):
    ...
    ValueError: The built-in 'pycsamt' plugin cannot be replaced.
    """
    if not isinstance(plugin, DomainPlugin):
        raise TypeError("register_plugin expects a DomainPlugin.")
    _ensure_builtin()
    if plugin.name == DEFAULT_PLUGIN:
        raise ValueError(
            f"The built-in {DEFAULT_PLUGIN!r} plugin cannot be "
            "replaced."
        )
    if plugin.name in _PLUGINS and not replace:
        raise ValueError(
            f"A plugin named {plugin.name!r} is already registered; "
            "pass replace=True to overwrite it."
        )
    _PLUGINS[plugin.name] = plugin
    return plugin


def discover_plugins() -> dict[str, str]:
    """Load plugins advertised in :data:`ENTRY_POINT_GROUP`.

    Each entry point must resolve to a :class:`DomainPlugin` or to a
    zero-argument callable returning one. Already-registered names
    are left alone. Discovery runs once per process; later calls
    only report.

    Returns
    -------
    dict of str to str
        Entry-point name mapped to ``"loaded"``, ``"skipped"`` (name
        already registered) or the error message of a failed load.
        Failures are also logged at WARNING level; they never stop
        the agent loop.
    """
    global _DISCOVERED
    _ensure_builtin()
    status: dict[str, str] = {}
    if _DISCOVERED:
        return status
    _DISCOVERED = True
    for ep in metadata.entry_points(group=ENTRY_POINT_GROUP):
        try:
            obj = ep.load()
            plugin = obj if isinstance(obj, DomainPlugin) else obj()
            if plugin.name in _PLUGINS:
                status[ep.name] = "skipped"
                continue
            register_plugin(plugin)
            status[ep.name] = "loaded"
        except Exception as exc:  # noqa: BLE001 - third-party code
            _log.warning("Agent plugin %r failed to load: %s",
                         ep.name, exc)
            status[ep.name] = f"{type(exc).__name__}: {exc}"
    return status


def get_plugin(name: str | None = None) -> DomainPlugin:
    """Return the plugin called *name* (built-in when ``None``).

    Parameters
    ----------
    name : str or None, default None
        Plugin name. ``None`` returns the built-in ``"pycsamt"``
        plugin. An unknown name triggers :func:`discover_plugins`
        once before failing.

    Returns
    -------
    DomainPlugin

    Raises
    ------
    KeyError
        If no plugin of that name is registered or discoverable.

    Examples
    --------
    >>> from pycsamt.agents.plugins import get_plugin
    >>> p = get_plugin()
    >>> p.name, p.registry.normalise("PINN")
    ('pycsamt', 'pinn_inversion')
    """
    _ensure_builtin()
    key = DEFAULT_PLUGIN if name is None else name
    if key not in _PLUGINS:
        discover_plugins()
    try:
        return _PLUGINS[key]
    except KeyError:
        raise KeyError(
            f"No agent plugin named {key!r}; available: "
            f"{sorted(_PLUGINS)}"
        ) from None


def available_plugins() -> dict[str, str]:
    """Return registered plugin names mapped to their descriptions."""
    _ensure_builtin()
    return {k: p.description for k, p in sorted(_PLUGINS.items())}
