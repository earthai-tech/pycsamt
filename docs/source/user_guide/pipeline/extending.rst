.. _pipeline-extending:

Extending The Pipeline
======================

:doc:`steps` describes the *built-in* extension path: a new operation is
implemented in ``pycsamt.emtools``, given a :class:`~pycsamt.pipeline.StepSpec`
entry in ``pycsamt/pipeline/_registry.py``, and reviewed with tests and
documentation before it ships in a pyCSAMT release. That path is intentional
for anything that should become part of the shared scientific catalogue.

It is the wrong path for a step that only makes sense for one project, one
instrument vendor, or one organisation's internal correction. For that case,
pyCSAMT accepts a second kind of step: a :term:`pipeline plugin`, registered
at runtime with ``register_step`` and never touching pyCSAMT's own source
tree. Once registered, a plugin step behaves exactly like a built-in one --
it can be used from ``Step``, ``Pipeline``, presets, configuration files, and
the CLI.

Two Extension Paths
-------------------

.. list-table::
   :header-rows: 1
   :widths: 22 39 39

   * -
     - Built-in step
     - Plugin step
   * - Where it lives
     - ``pycsamt/pipeline/_registry.py``
     - Any installed package, or a project script
   * - How it ships
     - Reviewed pull request, released with pyCSAMT
     - ``pip install`` of a third-party package, or a runtime call
   * - Registered by
     - The literal ``StepSpec`` list in the registry module
     - ``register_step`` from Python, or a
       :term:`plugin entry-point group` callable
   * - ``StepSpec.origin``
     - ``"builtin"``
     - ``"plugin"``
   * - Good fit for
     - Operations useful across surveys and organisations
     - Site-specific corrections, vendor formats, internal QC rules

Registering A Step
------------------

``register_step`` takes a fully-built :class:`~pycsamt.pipeline.StepSpec` and
inserts it into the same :term:`step registry` that the 47 built-in steps
live in:

.. code-block:: pycon
   :linenos:

   >>> from pycsamt.pipeline import StepSpec, register_step, lookup_step
   >>> def scale_amplitude(sites, factor: float = 2.0):
   ...     """Multiply every impedance value in *sites* by *factor* (toy transform)."""
   ...     return sites
   ...
   >>> spec = register_step(
   ...     StepSpec(
   ...         code="DEMO001",
   ...         name="scale_amplitude",
   ...         label="Demo Amplitude Scale",
   ...         category="demo",
   ...         override_fn=scale_amplitude,
   ...         defaults={"factor": 2.0},
   ...     )
   ... )
   >>> spec.origin
   'plugin'
   >>> lookup_step("DEMO001") is spec
   True

``origin`` is always stamped ``"plugin"`` by ``register_step`` itself,
regardless of what the caller passed -- a plugin author cannot accidentally
label a step ``"builtin"``. Once registered, the step is usable exactly like
any other:

.. code-block:: pycon
   :linenos:

   >>> from pycsamt.pipeline import Step
   >>> Step("DEMO001", factor=3.0).params
   {'factor': 3.0}

Registration is a mutation of one process-wide registry, the same way
``pycsamt.forward.maxwell``'s :ref:`backend registry <forward_maxwell_backends>`
is process-wide. A second registration under the same code or name fails
rather than silently swapping the implementation, because a silent swap could
change numerical behaviour somewhere else in a long-running session:

.. code-block:: pycon
   :linenos:

   >>> register_step(
   ...     StepSpec(
   ...         code="DEMO001",
   ...         name="scale_amplitude",
   ...         label="duplicate",
   ...         category="demo",
   ...         override_fn=scale_amplitude,
   ...     )
   ... )
   Traceback (most recent call last):
   ...
   ValueError: Pipeline step code='DEMO001' or name='scale_amplitude' is already registered.  Pass replace_existing=True to overwrite it.

Passing ``replace_existing=True`` allows a deliberate overwrite, including of
a built-in step -- useful for a site that wants to patch one operation's
defaults without forking pyCSAMT. Formally, if :math:`R` is the current
registry mapping codes to specs and a caller registers spec :math:`s` under
code :math:`c`,

.. math::

   R' =
   \begin{cases}
   R \cup \{c: s\}, & c \notin R,\\
   R \setminus \{c: R[c]\} \cup \{c: s\}, & c \in R \text{ and } \texttt{replace\_existing},\\
   \text{raise } \texttt{ValueError}, & c \in R \text{ and not } \texttt{replace\_existing}.
   \end{cases}

``register_step`` also validates the spec before inserting it: it calls
``spec.get_fn()`` once, so a typo'd module path fails immediately at
registration time rather than three steps into a pipeline run:

.. code-block:: pycon
   :linenos:

   >>> register_step(
   ...     StepSpec(
   ...         code="DEMO_BAD",
   ...         name="demo_bad",
   ...         label="Bad",
   ...         category="demo",
   ...         mod="my_package.pipeline_steps",
   ...         fn_name="not_a_real_function",
   ...     )
   ... )
   Traceback (most recent call last):
   ...
   ModuleNotFoundError: No module named 'my_package'

The registry is left untouched when validation fails; ``DEMO_BAD`` never
appears in ``lookup_step`` or ``list_steps``.

Removing A Step
---------------

``unregister_step`` reverses a registration, by code or by name:

.. code-block:: pycon
   :linenos:

   >>> from pycsamt.pipeline import unregister_step
   >>> unregister_step("DEMO001")
   >>> lookup_step("DEMO001")
   Traceback (most recent call last):
   ...
   KeyError: "No pipeline step found for 'DEMO001'.  Call list_steps() to see all available steps."

``unregister_step("SOME_CODE", missing_ok=True)`` is a safe no-op when the
step was never registered -- the form to reach for in a test fixture's
teardown, so a test does not fail just because an earlier test already
cleaned up the same code.

Packaging A Plugin
------------------

A plugin package announces itself through the ``pycsamt.pipeline.steps``
:term:`plugin entry-point group`. The package's own ``pyproject.toml`` points
at a zero-argument callable that performs the registration:

.. code-block:: toml
   :linenos:

   [project.entry-points."pycsamt.pipeline.steps"]
   demo = "demo_pipe_plugin:register"

``demo`` is the plugin's display name -- it is what shows up in
``pycsamt pipe plugins`` output. ``demo_pipe_plugin:register`` is
``module:callable``, resolved with :mod:`importlib.metadata` the same way
``console_scripts`` entry points resolve for the ``pycsamt`` command itself.
The referenced module calls ``register_step`` for whatever it contributes:

.. code-block:: python
   :linenos:

   """Toy pyCSAMT pipeline plugin used to smoke-test the plugin discovery API."""

   from pycsamt.pipeline import StepSpec, register_step


   def scale_amplitude(sites, factor: float = 2.0):
       """Multiply every impedance value in *sites* by *factor* (toy transform)."""
       return sites


   def register() -> None:
       register_step(
           StepSpec(
               code="DEMO001",
               name="scale_amplitude",
               label="Demo Amplitude Scale",
               category="demo",
               mod="demo_pipe_plugin",
               fn_name="scale_amplitude",
               defaults={"factor": 2.0},
           )
       )

Installing the package with ``pip install demo-pipe-plugin`` (or
``pip install -e .`` during development) is enough for
:mod:`importlib.metadata` to see the entry point. It is not, by itself,
enough for the step to appear in ``STEP_REGISTRY`` -- see the next section.

Explicit Discovery
------------------

Installing a plugin package never runs its ``register`` function by itself.
:term:`Plugin discovery` -- the act of scanning the entry-point group and
calling every callable found there -- only happens when
``pycsamt.pipeline.discover_plugins`` is called:

.. code-block:: pycon
   :linenos:

   >>> from pycsamt.pipeline import discover_plugins
   >>> discover_plugins()
   [PluginLoadResult(name='demo', ok=True, error=None)]

This is deliberate, for two reasons. First, running arbitrary third-party
code merely because a package happens to be installed is a real security
surface pyCSAMT does not want to open by default. Second, the scan itself is
not free: :func:`importlib.metadata.entry_points` has to enumerate every
installed distribution to find the ones that declare
``pycsamt.pipeline.steps``, and on the Anaconda environment used to build
this page -- several hundred installed packages -- a single scan measured
several seconds. Paying that cost on every ``import pycsamt.pipeline``, or on
every single ``pycsamt pipe`` invocation, would make the common case (no
plugins at all) noticeably slower for everyone in order to serve a case most
users never hit.

A plugin that fails to load is reported, not fatal to the others -- and this
demo plugin's own ``register()`` is itself a case in point. It calls
``register_step`` unconditionally with no ``replace_existing``, so calling
``discover_plugins`` a *second* time in the same process reports it as
failed, not because anything broke, but because ``DEMO001`` is already
registered from the first call:

.. code-block:: pycon
   :linenos:

   >>> results = discover_plugins()
   >>> [(r.name, r.ok) for r in results]
   [('demo', False)]

This is why pyCSAMT's own CLI calls ``discover_plugins`` at most once per
process -- either inside ``pycsamt pipe plugins`` or from the ``--with-plugins``
group flag, never both (see :ref:`CLI Plugin Discovery <pipeline-extending-cli>`
below). A plugin author who expects ``register()`` to be called more than
once per process should make it tolerate that, either by checking
``spec.code not in STEP_REGISTRY`` first or by passing
``replace_existing=True``. Pass ``on_error="raise"`` to stop at the first
failure instead of collecting a :class:`PluginLoadResult` for it -- useful in
a CI job that should fail loudly on a broken plugin rather than merely warn.

.. _pipeline-extending-cli:

CLI Plugin Discovery
--------------------

``pycsamt pipe plugins`` always discovers, because discovery is the entire
point of that command:

.. code-block:: console
   :linenos:

   pycsamt pipe plugins

.. code-block:: text
   :linenos:

   Discovered 1 pipeline plugin(s):
     demo                     ok

   Registered 1 plugin step(s):
     DEMO001    scale_amplitude              [demo]  Demo Amplitude Scale

Every other ``pycsamt pipe`` subcommand leaves discovery off by default, for
the same latency reason explained above. Trying to use a plugin step's code
without opting in fails with a hint toward the fix:

.. code-block:: console
   :linenos:

   pycsamt pipe run --steps DEMO001 --survey ./edis/ --dry-run

.. code-block:: text
   :linenos:

   Error: Invalid value for '--steps': Unknown step 'DEMO001'.  Run  pycsamt pipe steps  to see all available steps.  If this is a plugin step, pass  pipe --with-plugins  or run  pycsamt pipe plugins  first.

``--with-plugins`` is a flag on the ``pipe`` group itself, not on individual
subcommands such as ``run`` -- it has to run before Click parses
``--steps``, because ``--steps`` is validated against the registry at parse
time. Setting the ``PYCSAMT_PIPELINE_LOAD_PLUGINS`` environment variable has
the same effect without retyping the flag on every invocation. With
discovery opted into, the same plugin step resolves in one shot:

.. code-block:: console
   :linenos:

   pycsamt pipe --with-plugins run --steps DEMO001 --survey ./edis/ --dry-run

.. code-block:: text
   :linenos:

   Pipeline  'cli_pipeline'  ───────────────────────────────────────────────  1 step
     ( 1) scale_amplitude  [DEMO001]  Demo Amplitude Scale  factor=2.0
   ────────────────────────────────────────────────────────────────────────────────

   Sites   : 25
   Steps   : 1
   Out dir : pipe_results (default)

   Dry run — no processing performed.

Captured against a real 25-site WILLY survey directory. ``pycsamt pipe steps
--info`` (also under ``--with-plugins``) formats a plugin step the same way
it formats a built-in one:

.. code-block:: console
   :linenos:

   pycsamt pipe --with-plugins steps --info DEMO001

.. code-block:: text
   :linenos:

   DEMO001  Demo Amplitude Scale  [demo]
     name     : scale_amplitude
     function : demo_pipe_plugin.scale_amplitude
     defaults : factor=2.0
     qc plots : —
     returns  : Sites (transform)

Plugin Origin
-------------

Every :class:`~pycsamt.pipeline.StepSpec` carries an ``origin`` field so
built-in and plugin steps can be told apart programmatically:

.. code-block:: pycon
   :linenos:

   >>> from pycsamt.pipeline import STEP_REGISTRY
   >>> sorted(s.code for s in STEP_REGISTRY.values() if s.origin == "plugin")
   ['DEMO001']
   >>> lookup_step("NR001").origin
   'builtin'

A script that wants to report only what a survey's specific pipeline
configuration relies on beyond the shared catalogue can filter on this field
rather than maintaining a separate list of "steps I added."

First-Party AI Steps
--------------------

pyCSAMT ships one opt-in step of its own, built on exactly the mechanism
this page teaches for third-party plugins: ``AI001`` /
``audit_survey``, wrapping :func:`pycsamt.ai.domain_gap.audit.audit_survey`
as a diagnostic pipeline step. It is not a built-in the way the other 50
steps are, and it is not discovered through the entry-point mechanism above
either -- it is registered directly by
:func:`pycsamt.pipeline.register_ai_steps`:

.. code-block:: pycon
   :linenos:

   >>> from pycsamt.pipeline import register_ai_steps, lookup_step
   >>> registered = register_ai_steps()
   >>> [s.code for s in registered]
   ['AI001']
   >>> lookup_step("AI001").origin
   'plugin'

Why this needs its own opt-in, distinct from ``register_step`` collisions or
entry-point scanning: resolving ``AI001`` imports ``pycsamt.ai``, and
``pycsamt.ai``'s own package ``__init__`` eagerly imports
``pycsamt.ai.nets.drcnn``, which imports ``torch`` at module level (a
deliberate choice there so its classes stay picklable for checkpointing --
see that module's own comment). That is a real, multi-second cost the CLI
must not force on every pipeline user merely for ``pycsamt pipe run``. The
``pipe`` group's ``--with-ai-steps`` flag (or the
``PYCSAMT_PIPELINE_LOAD_AI_STEPS`` environment variable) opts in the same
way ``--with-plugins`` does, and for the same structural reason: it is a
*group*-level flag so it runs before ``--steps`` is parsed and validated:

.. code-block:: console
   :linenos:

   pycsamt pipe --with-ai-steps steps --info AI001

.. code-block:: text
   :linenos:

   AI001  AI Domain-Gap Survey Audit  [ai]
     name     : audit_survey
     function : pycsamt.pipeline.ai_steps.qc_audit_survey
     defaults : —
     qc plots : —
     returns  : Sites unchanged (diagnostic)

Captured running the step for real against the 25-station WILLY L22 line:

.. code-block:: console
   :linenos:

   pycsamt pipe --with-ai-steps run --steps AI001 \
       --survey data/AMT/WILLY_DATA/L22PLT --out results/line22_audit

.. code-block:: text
   :linenos:

   Survey audit (generated 2026-08-14T08:05:48Z)
     Stations: 25 input, 25 included, 0 excluded
     Frequency grid: matched
     Impedance coverage: 100.0%
     Frequency range: 1.008-10400 Hz
     Declared error / |Z|: p05=0.0214, p50=0.0839, p95=0.3492
     Station spacing (m): min=59.9, median=100.0, max=120.5
     CRS declared: False
     Elevation coverage: 100.0%
     Dimensionality: n=1325, 1D=9.5%, 2D=18.0%, 3D=72.5%
     Strike (consensus): -29.0 deg (IQR 118.9 deg)
     Static shift log10 sigma: 0.1219
     Distortion sigma: gain(log10)=0.0000, twist_deg=14.63, shear=0.3680, anisotropy=0.1234

``qc_audit_survey`` (:mod:`pycsamt.pipeline.ai_steps`) also accepts a
``report_path`` parameter that writes the full report as JSON via
:meth:`~pycsamt.ai.domain_gap.audit.SurveyAuditReport.write_json`, for
runs where the structured report should be kept alongside the rest of the
pipeline's output.

Testing Plugin Steps
--------------------

Plugin registration is process-wide mutable state, the same way the Maxwell
:ref:`backend registry <forward_maxwell_backends>` is. A test that registers
a plugin step must clean it up, or it leaks into unrelated tests that run
later in the same process:

.. code-block:: pycon
   :linenos:

   >>> import pytest
   >>> from pycsamt.pipeline import StepSpec, register_step, unregister_step
   >>>
   >>> @pytest.fixture()
   ... def demo_plugin_step():
   ...     spec = register_step(
   ...         StepSpec(
   ...             code="DEMO001",
   ...             name="scale_amplitude",
   ...             label="Demo Amplitude Scale",
   ...             category="demo",
   ...             override_fn=scale_amplitude,
   ...         )
   ...     )
   ...     yield spec
   ...     unregister_step("DEMO001", missing_ok=True)
   ...
   >>> def test_plugin_step_usable_in_a_pipeline(demo_plugin_step, sites):
   ...     from pycsamt.pipeline import Pipeline, Step
   ...     pipe = Pipeline([("scale", Step("DEMO001"))])
   ...     result = pipe.run(sites, outdir=None, save_plots=False)
   ...     assert result.ok

``missing_ok=True`` in the teardown matters: if the test body itself already
called ``unregister_step``, or failed before registration completed, a
strict teardown would raise a second, unrelated error that masks the real
failure.

Troubleshooting
---------------

``ValueError: ... is already registered``
    Another step already claims this code or name. Pass
    ``replace_existing=True`` if overwriting it is intentional, or choose a
    less generic code -- a vendor or project prefix such as
    ``ACME_DEMO001`` avoids colliding with both the built-in catalogue and
    other plugins.

A plugin step is "unknown" from the CLI
    Discovery did not run. Pass ``pycsamt pipe --with-plugins`` before the
    subcommand, run ``pycsamt pipe plugins`` first, or set
    ``PYCSAMT_PIPELINE_LOAD_PLUGINS=1``.

The first ``pycsamt pipe plugins`` call feels slow
    That is :func:`importlib.metadata.entry_points` enumerating every
    installed package, not the plugin's own code. It is a one-time cost per
    process, not per step.

A plugin's ``register`` function raised
    ``discover_plugins`` reports it as a failed :class:`PluginLoadResult`
    and warns, rather than crashing the whole discovery pass. Check
    ``pycsamt pipe plugins`` output, or ``PluginLoadResult.error``, for the
    underlying exception message.

Related Pages
-------------

* :doc:`steps` explains the built-in step registry and the reviewed path for
  adding a step to pyCSAMT itself.
* :doc:`concepts` explains how registered steps compose into a pipeline run.
* :doc:`cli_pipe` explains the ``pycsamt pipe`` command group in full,
  including ``--with-plugins``.
* :doc:`caching` explains ``cache=``/``--cache``, including how it applies
  to the opt-in ``AI001`` step introduced on this page.
