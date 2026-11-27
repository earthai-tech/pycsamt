.. _agents-plugins:

Agent Plugins
=============

The agent loop in :mod:`pycsamt.agents` is general (route a request,
build a plan, gate it, run a chain of agents, record what happened),
but three of its pieces are specific to CSAMT/MT processing: the
catalogue of workflows, the rules that decide whether a plan may run,
and the run record. :mod:`pycsamt.agents.plugins` names those pieces as
:class:`typing.Protocol` classes, so another package can supply its own
domain without forking pyCSAMT. A bundle of the three is an
:term:`agent plugin`.

.. seealso::

   :doc:`review_gate`
       How the validator of a plugin is used before execution.

   :doc:`orchestrator`
       The dispatcher that looks up the active plugin.

   :mod:`pycsamt.agents.plugins`
       API reference.


The Three Protocols
-------------------

The protocols are structural: a plain class that has the listed
attributes and methods satisfies them, no inheritance is needed.

.. list-table::
   :header-rows: 1
   :widths: 22 78

   * - Protocol
     - Members and contract
   * - ``WorkflowRegistry``
     - ``name``; ``descriptions`` (every known workflow id mapped to one
       line; ids absent here are unknown and get blocked);
       ``valid_workflows`` (ids that go through the plan gate);
       ``dispatchable`` (ids that have a step chain);
       ``classify(text, default=None)``, ``normalise(workflow)``,
       ``family(workflow)`` and ``steps(workflow)``. ``steps`` returns
       ``(step_name, agent_class_name, input_fn, description[, required])``
       tuples, the shape the orchestrator runs. ``classify`` must be
       deterministic, because the route cross-check relies on it.
   * - ``PlanValidator``
     - ``errors(plan)`` (any error blocks the run),
       ``needs_review(plan)`` (``True`` requires approval of the plan
       hash), ``cross_check(text, intent, workflow)`` (router-level
       comparison, returns a reason or ``None``) and
       ``route_flags(text, workflow)`` (orchestrator-level comparison,
       returned as plan risk flags).
   * - ``RecordWriter``
     - ``emit(event, **fields)``, ``snapshot()``, ``persist(path)`` and
       ``corpus_manifest()``. Subclassing
       ``pycsamt.agents._provenance.RunRecord`` provides the first three.

:class:`~pycsamt.agents.plugins.DomainPlugin` bundles a registry, a
validator and a ``record_factory(payload) -> RecordWriter`` under a
name. It checks each member against its protocol when it is built and
raises :class:`TypeError` when one does not fit.


The Built-In Plugin
-------------------

``get_plugin()`` with no name returns the built-in ``"pycsamt"`` plugin.
Its adapters delegate to the existing workflow tables, plan validation
and route checks, so the default behaviour of the agents is the same as
before the interface existed.

.. code-block:: pycon
   :linenos:

   >>> from pycsamt.agents.plugins import get_plugin
   >>> base = get_plugin()
   >>> base.name, base.registry.normalise("PINN")
   ('pycsamt', 'pinn_inversion')
   >>> base.registry.classify("compute the phase tensor ellipses")
   'pt_psection'
   >>> base.registry.family("pt_psection")
   'phase_analysis'

The built-in name is reserved: ``register_plugin`` refuses a plugin
called ``"pycsamt"``.

In this release, :class:`~pycsamt.agents.WorkflowOrchestratorAgent`
consults the active plugin for the list of dispatchable workflows, the
step chain, the orchestrator-level route flags, plan validation and the
review decision. The intent router and the run record still use the
built-in pieces, and the agent class names in a step chain are resolved
against the pyCSAMT agent registry. A plugin whose steps name agents of
its own can therefore be registered, classified and validated, but its
chain cannot yet be executed by the orchestrator.


Registering A Plugin
--------------------

``examples/agent_plugins/mtpy_example`` in the repository is a complete
plugin for three MTpy-style workflows: reading EDI files and plotting
apparent resistivity and phase, computing the phase tensor, and
estimating strike from it. MTpy is needed only by its step functions;
importing, registering, classifying and validating work without it.
Run from the repository root:

.. code-block:: pycon
   :linenos:

   >>> import sys
   >>> sys.path.insert(0, "examples/agent_plugins")
   >>> from mtpy_example import make_plugin
   >>> from pycsamt.agents.plugins import available_plugins, register_plugin

   >>> plugin = register_plugin(make_plugin(), replace=True)
   >>> sorted(available_plugins())
   ['mtpy_example', 'pycsamt']
   >>> plugin.registry.classify("estimate strike for site 12")
   'strike'
   >>> plugin.registry.family("strike"), plugin.registry.normalise("PT")
   ('phase_tensor', 'phase_tensor')
   >>> plugin.validator.cross_check(
   ...     "plot the phase tensor ellipses", "workflow", "read_edi_plot_rhophi")
   "LLM workflow 'read_edi_plot_rhophi' disagrees with keyword route 'phase_tensor'."
   >>> plugin.validator.route_flags("plot the phase tensor ellipses", "strike")
   []

``strike`` belongs to the ``phase_tensor`` family of this plugin, so a
language model that chose ``strike`` for a phase-tensor request is not
flagged, while one that chose the resistivity plot is. The example
validator also checks its own parameters: ``freq_band`` must lie in
:math:`10^{-5}` to :math:`10^{5}` Hz and increase,
``strike_tolerance_deg`` must be in :math:`(0, 90]` degrees and
``n_periods`` must be an integer between 1 and 500.

``replace=True`` overwrites a plugin of the same name; without it a
second registration raises :class:`ValueError`.

To make the orchestrator use a plugin, set its ``plugin_name`` class
attribute, for example in a subclass:

.. code-block:: python

   from pycsamt.agents import WorkflowOrchestratorAgent

   class MTpyOrchestrator(WorkflowOrchestratorAgent):
       plugin_name = "mtpy_example"


Shipping A Plugin
-----------------

A third-party distribution advertises its plugin in the
``pycsamt.agent_plugins`` entry-point group, distinct from the
``pycsamt.pipeline.steps`` group used by pipeline step plugins. The entry
point may resolve to a ``DomainPlugin`` instance or to a zero-argument
callable that returns one:

.. code-block:: toml

   [project.entry-points."pycsamt.agent_plugins"]
   mtpy_example = "mtpy_example:make_plugin"

After the package is installed, ``get_plugin("mtpy_example")`` finds it:
an unknown name triggers :func:`~pycsamt.agents.plugins.discover_plugins`
once before :class:`KeyError` is raised. Discovery runs once per
process, never replaces a plugin that is already registered, and logs
a failing entry point at WARNING level without stopping the agent loop.
``discover_plugins()`` returns the status of each entry point
(``"loaded"``, ``"skipped"`` or the error message), which is the first
thing to check when a plugin does not appear.
