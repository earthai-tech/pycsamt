.. _development-api-policy:

API Policy
==========

This page defines how the pyCSAMT v2 API is structured, how contributors
should expose new functionality, and what compatibility promises users can
expect.  It is the development contract for the project: before adding a new
module, function, class, CLI command, agent, or pipeline step, check this page.

pyCSAMT v2 is larger than a single numerical package.  It contains scientific
data containers, EM processing tools, inversion interfaces, AI models,
AI-assisted agents, declarative pipelines, plotting helpers, a command-line
interface, and application-facing view objects.  The API policy keeps those
parts discoverable without making every internal implementation detail
permanent.


API design goals
----------------

The v2 API follows six principles.

.. list-table::
   :header-rows: 1
   :widths: 24 76

   * - Principle
     - Meaning for contributors
   * - Stable imports
     - Users should import from documented public namespaces, not from private
       implementation files.
   * - Scientific clarity
     - Public functions must make units, coordinate assumptions, array shapes,
       and impedance/resistivity conventions explicit.
   * - Composable workflows
     - Processing, inversion, agents, and CLI commands should share structured
       inputs and outputs so they can be chained.
   * - Optional heavy dependencies
     - AI, LLM, plotting, GIS, and backend-specific features must fail lazily
       and clearly when optional dependencies are missing.
   * - Reproducibility
     - User-facing results should carry enough metadata to recreate the
       operation: parameters, versions, backend, warnings, output paths, and
       diagnostics.
   * - Compatibility discipline
     - Public names should not be renamed, moved, or removed casually.  When a
       change is required, use a documented deprecation path.


Public namespace map
--------------------

Use this table when deciding where a new API belongs.

.. list-table::
   :header-rows: 1
   :widths: 25 35 40

   * - Namespace
     - Role
     - Stability expectation
   * - ``pycsamt``
     - Top-level convenience imports such as version, backend helpers, and the
       most common pipeline symbols.
     - Stable.  Keep small and intentional.
   * - ``pycsamt.api``
     - Public configuration, result wrappers, style controls, CLI option
       helpers, plotting configuration, and lightweight view helpers.
     - Stable public facade.
   * - ``pycsamt.agents``
     - AI-assisted workflow agents, lazy-loaded to avoid hard LLM dependency at
       import time.
     - Public, but individual agent capabilities may be marked experimental.
   * - ``pycsamt.pipeline``
     - Declarative processing engine, steps, presets, and pipeline runtime
       configuration.
     - Stable public workflow API.
   * - ``pycsamt.cli``
     - Click-based command-line application and command groups.
     - Stable at command and option level once documented.
   * - ``pycsamt.inversion``
     - Physics-based inversion interfaces, models, results, backends, and
       workflow helpers.
     - Stable for documented classes and functions.
   * - ``pycsamt.models``
     - File builders, runners, result readers, and configuration objects for
       engines such as Occam2D and ModEM.
     - Public for documented engine-facing objects.
   * - ``pycsamt.ai``
     - Neural inversion models, training utilities, AI processing tools, and AI
       plotting helpers.
     - Experimental-to-stable depending on page and docstring status.
   * - ``pycsamt.emtools``
     - Electromagnetic processing, diagnostics, source effects, tensor tools,
       plots, and legacy-compatible helpers.
     - Stable for documented functions; legacy aliases follow deprecation
       rules.
   * - ``pycsamt.seg``, ``pycsamt.jones``, ``pycsamt.zonge``, ``pycsamt.tdem``
     - Format-specific readers, writers, parsers, transforms, and survey
       objects.
     - Stable for documented readers and data objects.
   * - ``pycsamt.site``
     - Site, survey, profile, location, selection, export, and reporting
       helpers.
     - Stable for documented site and survey interfaces.
   * - ``pycsamt.interp``
     - Geological, lithological, hydrological, petrophysical, uncertainty, and
       export interpretation APIs.
     - Stable for documented interpretation workflows.
   * - ``pycsamt.backends``
     - Backend detection and backend selection for AI/neural runtimes.
     - Stable for selection helpers; backend implementations may evolve.
   * - ``pycsamt.utils`` and ``pycsamt.compat``
     - Shared utilities, compatibility shims, and migration helpers.
     - Public only when explicitly documented or exported through
       ``__all__``.


Public, private, experimental, deprecated
-----------------------------------------

Every API must fall into one of these categories.

.. list-table::
   :header-rows: 1
   :widths: 18 30 28 24

   * - Category
     - How to recognize it
     - Compatibility promise
     - Developer action
   * - Public stable
     - Documented in ``docs/source/api`` or a user/developer guide, exported in
       ``__all__``, and importable from a public namespace.
     - Keep compatible across v2 minor releases.
     - Add tests, docstring, docs entry, and release note.
   * - Public experimental
     - Documented as experimental or provisional.  Common for young AI,
       backend, and agent features.
     - May change, but changes must be explained.
     - Add warnings in docs and prefer keyword-only extension points.
   * - Private
     - Name starts with ``_`` or lives in an underscored module such as
       ``pycsamt.pipeline._registry``.
     - No compatibility promise.
     - Do not document as an import target.  Re-export stable wrappers instead.
   * - Deprecated
     - Emits ``FutureWarning`` or is listed in migration/release notes.
     - Must remain available until the announced removal version.
     - Provide replacement guidance and tests for the warning.
   * - Removed
     - Not importable, or guarded by an explanatory ``__getattr__`` error.
     - No runtime support.
     - Keep migration guidance where useful.


The import rule
---------------

Users should not need to know the internal file layout.  Public imports should
come from stable package-level facades.

Preferred public imports:

.. code-block:: python
   :linenos:

   from pycsamt.api import read_edis, configure_style
   from pycsamt.pipeline import Pipeline, Step, list_steps
   from pycsamt.agents import DataQCAgent, WorkflowOrchestratorAgent
   from pycsamt.backends import set_backend, get_backend

Avoid exposing internal implementation paths as the primary documented import:

.. code-block:: python
   :linenos:

   # Avoid in user documentation unless the module itself is the public API.
   from pycsamt.pipeline._registry import STEP_REGISTRY
   from pycsamt.agents._base import BaseAgent

There are valid exceptions.  Some submodules, such as
``pycsamt.models.occam2d.builder`` or ``pycsamt.inversion.workflow``, are
domain-specific public entry points.  When that is intended, document the
module and include its public names in ``__all__``.


Top-level ``pycsamt`` policy
----------------------------

The top-level package is for the most common cross-package shortcuts only.  It
currently exposes:

* ``__version__``
* core subpackages such as ``interp``, ``inversion``, and ``tdem``
* pipeline shortcuts such as ``Pipeline``, ``Step``, ``configure_pipe``,
  ``reset_pipe``, and ``PYCSAMT_PIPE``
* backend helpers such as ``auto_detect``, ``get_backend``,
  ``get_backend_instance``, ``list_backends``, and ``set_backend``

Do not add every new class to ``pycsamt.__all__``.  Add a top-level shortcut
only when all of these are true:

* the object is broadly useful across workflows;
* importing it does not force optional heavy dependencies;
* the name is stable enough to support through v2;
* the object already has a natural package-level home.

Good example:

.. code-block:: python
   :linenos:

   from pycsamt import Pipeline, Step

More specific example:

.. code-block:: python
   :linenos:

   from pycsamt.models.occam2d import InputBuilder, OccamRunner


The ``pycsamt.api`` facade
--------------------------

``pycsamt.api`` is the public facade for cross-cutting configuration and
application-facing helpers.  Use it for stable objects that are not tied to one
scientific algorithm but affect how pyCSAMT behaves or returns data.

Examples include:

* global style configuration;
* plotting/export configuration;
* CLI option and parameter helpers;
* result wrappers such as ``APIResult``;
* data-view helpers such as ``read_edi`` and ``read_edis``;
* agent LLM configuration such as ``configure_agents``;
* pipeline runtime configuration such as ``configure_pipe``.

When adding to ``pycsamt.api``:

* keep objects lightweight at import time;
* avoid importing AI, GIS, plotting, or inversion engines unless needed;
* define an explicit ``__all__`` in the source module;
* re-export from ``pycsamt/api/__init__.py`` only for stable public names;
* include examples in the relevant guide page.

Pattern:

.. code-block:: python
   :linenos:

   # pycsamt/api/example.py
   from dataclasses import dataclass

   @dataclass
   class ExampleConfig:
       enabled: bool = True

   PYCSAMT_EXAMPLE = ExampleConfig()

   def configure_example(**kwargs):
       for key, value in kwargs.items():
           setattr(PYCSAMT_EXAMPLE, key, value)
       return PYCSAMT_EXAMPLE

   def reset_example():
       global PYCSAMT_EXAMPLE
       PYCSAMT_EXAMPLE = ExampleConfig()
       return PYCSAMT_EXAMPLE

   __all__ = ["ExampleConfig", "PYCSAMT_EXAMPLE",
              "configure_example", "reset_example"]


Agents API policy
-----------------

Agents are public workflow components.  They must be usable from Python,
coordinators, orchestrators, notebooks, and eventually web/CLI workflows.

Public agents live under ``pycsamt.agents`` and are lazy-loaded through
``pycsamt/agents/__init__.py``.  New agents should follow the existing
``BaseAgent`` and ``AgentResult`` contract.

Required behavior for every public agent:

* inherit from ``BaseAgent`` unless there is a strong architectural reason not
  to;
* implement ``execute(...)`` with a dictionary-like input contract;
* return an ``AgentResult`` or a compatible mapping with ``success``, ``data``,
  ``warnings``, ``figures``, and ``metadata`` semantics;
* support ``dry_run`` when the operation can be planned without execution;
* keep LLM calls optional;
* report missing optional dependencies through clear exceptions or failed
  ``AgentResult`` objects;
* document input keys, output data keys, examples, and typical chained usage.

Recommended minimal agent skeleton:

.. code-block:: python
   :linenos:

   from pycsamt.agents import AgentResult, BaseAgent

   class MyAgent(BaseAgent):
       """Short scientific purpose of the agent."""

       def execute(self, inputs, *, dry_run=False, **kwargs):
           config = dict(inputs or {})

           if dry_run:
               return AgentResult(
                   success=True,
                   data={"planned": True, "steps": ["load", "qc"]},
                   metadata={"agent": self.__class__.__name__},
               )

           # Run the concrete operation here.
           return AgentResult(
               success=True,
               data={"result": "..."},
               warnings=[],
               metadata={"agent": self.__class__.__name__},
           )

   __all__ = ["MyAgent"]

After creating an agent:

* add it to the lazy map in ``pycsamt.agents``;
* add it to the catalogue docs;
* add focused tests under ``pycsamt/agents/tests``;
* update coordinator or orchestrator routing only when the agent is intended
  for automatic workflow selection.


Pipeline API policy
-------------------

The pipeline package is the declarative processing layer.  Its public API is
``Pipeline``, ``Step``, registry discovery helpers, presets, and runtime
configuration.

Pipeline steps must be stable, inspectable, and serializable.  A user should be
able to move between Python, YAML/JSON configuration, CLI execution, and agent
execution without changing the meaning of a step.

Recommended user-facing pattern:

.. code-block:: python
   :linenos:

   from pycsamt.pipeline import Pipeline, Step

   pipe = Pipeline([
       ("notch", Step("NR001", mains_hz=50)),
       ("band", Step("FREQ001", fmin=1e-3, fmax=1.0)),
       ("static_shift", Step("SS001", method="spatial_median")),
   ])

   result = pipe.run(sites, outdir="results/willy")

When adding pipeline functionality:

* register steps with stable codes and descriptive names;
* keep step parameters JSON/YAML serializable where possible;
* include category, description, required inputs, outputs, and failure modes;
* add the step to presets only after it has focused tests;
* preserve backward compatibility for existing step codes.

Do not rename a step code once documented.  If a step must be replaced, keep
the old code as an alias, emit a deprecation warning, and document the new
code.


Inversion and model API policy
------------------------------

Inversion APIs sit at the boundary between pyCSAMT objects and external
engines.  They must be explicit about what is computed inside pyCSAMT and what
is delegated to an external executable or backend.

Public inversion objects should document:

* supported data type: AMT, CSAMT, MT, EMAP, TDEM, or mixed;
* dimensionality: 1-D, 2-D, or 3-D;
* required input format and units;
* mesh/model assumptions;
* backend or executable requirements;
* output files and result object fields;
* reproducibility metadata.

Typical import style:

.. code-block:: python
   :linenos:

   from pycsamt.inversion import InversionConfig
   from pycsamt.models.occam2d import InputBuilder, OccamRunner
   from pycsamt.models.modem import ModEmConfig, ModEmRunner

Backend wrappers must not hide external failures.  If Occam2D, ModEM, PyGIMLi,
PyTorch, or TensorFlow is missing, raise a clear error or return a failed
structured result that names the missing backend.


AI and backend policy
---------------------

AI APIs are allowed to depend on optional packages, but importing pyCSAMT must
not require those packages.  Heavy imports should happen inside methods,
factory functions, or backend adapters.

Rules:

* use ``pycsamt.backends`` for neural backend selection where possible;
* keep model configuration serializable;
* expose training history and inference metadata;
* document tensor shapes and scaling conventions;
* never silently switch backend in a way that changes numerical results without
  recording it in metadata;
* keep pretrained model downloads explicit.

Example:

.. code-block:: python
   :linenos:

   from pycsamt.backends import set_backend
   from pycsamt.agents import Inv2DAgent

   set_backend("torch")
   result = Inv2DAgent().execute({
       "edi_dir": "data/willy",
       "output_dir": "results/willy_inv2d",
   })


CLI API policy
--------------

The command-line interface is a public API.  Users may depend on command names,
options, exit codes, output files, and machine-readable output.

When adding or changing a CLI command:

* keep command groups aligned with package domains, such as ``avg``, ``edi``,
  ``forward``, ``invert``, ``pipe``, ``site``, ``tdem``, and ``transform``;
* reuse shared option decorators and parameter types from ``pycsamt.api.cli``;
* support ``--help`` with concrete examples where possible;
* return non-zero exit codes for failures;
* do not change existing option names without deprecation;
* keep Python API and CLI behavior consistent.

Preferred structure:

.. code-block:: python
   :linenos:

   # pycsamt/cli/commands/example/run.py
   import click

   from pycsamt.api.cli import output_dir_option, verbose_option

   @click.command()
   @output_dir_option
   @verbose_option
   def run(output_dir, verbose):
       """Run the example workflow."""
       ...


Data and result contracts
-------------------------

Public workflows should return structured objects rather than unlabelled
tuples.  pyCSAMT uses several result styles, including ``AgentResult``,
``PipelineResult``, ``APIResult``, inversion result objects, and dictionaries
with explicit keys.

For new public APIs, prefer one of these patterns:

* dataclass result for domain objects with a stable schema;
* ``AgentResult`` for agents;
* ``PipelineResult`` or ``StepResult`` for pipeline execution;
* ``APIResult`` for application-facing wrappers;
* xarray, pandas, NumPy, or pathlib objects only when the data type is obvious
  from the function name and docstring.

Avoid returning long positional tuples from public functions.  If a tuple is
unavoidable for compatibility, document the order and provide a named
alternative.

Result metadata should include:

* pyCSAMT version when available;
* input paths or source identifiers;
* major parameter values;
* backend and optional dependency versions when relevant;
* warnings and quality flags;
* output paths for generated files;
* random seed for stochastic workflows.


Docstring and documentation requirements
----------------------------------------

A public API is not considered stable until it has documentation.

Minimum documentation for a public function:

* one-sentence summary;
* parameter section with types, units, and defaults;
* return section with shape and type;
* raises/warns section when important;
* at least one example for workflow-level APIs;
* API page entry through autosummary or explicit documentation.

Minimum documentation for a public class:

* summary of purpose;
* constructor parameters;
* important attributes;
* method examples;
* notes about optional dependencies;
* stability label if experimental.

Use the conventions in :ref:`development-docstring-style`.


Deprecation policy
------------------

Deprecation is allowed, but silent breakage is not.  Use deprecation when a
public name, argument, command, return field, or behavior must be replaced.

Required deprecation information:

* old name or behavior;
* replacement;
* version where deprecation starts;
* planned removal version;
* migration note in changelog or release notes;
* test that the warning is emitted.

Runtime warnings should use ``FutureWarning`` for user-facing deprecations.

Example:

.. code-block:: python
   :linenos:

   import warnings

   def old_reader(*args, **kwargs):
       warnings.warn(
           "old_reader is deprecated since v2.0.0 and will be removed in "
           "v2.2.0. Use read_edis instead.",
           FutureWarning,
           stacklevel=2,
       )
       return read_edis(*args, **kwargs)

Deprecation windows:

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - API type
     - Minimum support window
   * - Python function, class, or argument
     - At least one minor release after warning appears.
   * - CLI command or option
     - At least one minor release, preferably two for common commands.
   * - Pipeline step code
     - Keep as alias for at least two minor releases.
   * - File format reader behavior
     - Keep compatibility unless the old behavior is scientifically wrong or
       unsafe.  Document migration carefully.
   * - Private API
     - No deprecation window required.


Optional dependency policy
--------------------------

Optional dependencies must be imported lazily and fail with actionable
messages.  This applies especially to:

* LLM clients: Anthropic, OpenAI, Google Gemini;
* neural frameworks: PyTorch, TensorFlow;
* GIS/export libraries: GDAL, rasterio, geopandas, shapely;
* plotting stacks beyond base matplotlib;
* external inversion executables.

Good pattern:

.. code-block:: python
   :linenos:

   def run_torch_model(config):
       try:
           import torch
       except ImportError as exc:
           raise ImportError(
               "This operation requires PyTorch. Install pycsamt with the "
               "AI extras or install torch manually."
           ) from exc

       ...

Do not import optional heavy dependencies in ``pycsamt.__init__`` or in a
facade module unless the dependency is required for base installation.


Testing requirements
--------------------

New public APIs need tests at the correct level.

.. list-table::
   :header-rows: 1
   :widths: 32 68

   * - Change type
     - Required tests
   * - New pure function or data class
     - Unit tests for normal input, edge cases, and invalid input.
   * - New reader/writer
     - Round-trip or fixture-based tests with small sample files.
   * - New agent
     - Dry-run test, success-path test with small or mocked data, missing
       optional dependency test, and LLM-disabled test.
   * - New pipeline step
     - Registry lookup test, parameter validation test, execution test, and
       serialization test if applicable.
   * - New CLI command
     - ``CliRunner`` help test, success test, and failure test.
   * - New inversion/backend feature
     - Backend availability test plus a small deterministic run or a mocked
       backend test.
   * - Deprecation
     - Warning test and replacement behavior test.


Contributor checklist
---------------------

Before opening a pull request for a public API change, verify:

.. code-block:: text
   :linenos:

   [ ] The public import path is intentional.
   [ ] The name is included in __all__ only where it should be public.
   [ ] Optional dependencies are imported lazily.
   [ ] Inputs, units, shapes, and defaults are documented.
   [ ] Return type is structured and documented.
   [ ] Errors and warnings are actionable.
   [ ] Tests cover normal, edge, and failure paths.
   [ ] CLI behavior is consistent with Python API behavior, if applicable.
   [ ] Agent or pipeline metadata is reproducible, if applicable.
   [ ] Deprecation notes exist for renamed or removed public behavior.
   [ ] Autosummary/API docs build without new warnings.
   [ ] Release notes or changelog mention user-visible changes.


Decision guide
--------------

Use this quick guide when unsure where to put a new feature.

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Feature
     - Put it here
   * - Global plotting, style, CLI, agent, or pipeline configuration
     - ``pycsamt.api``
   * - A reusable processing operation that can be sequenced
     - ``pycsamt.emtools`` plus a ``pycsamt.pipeline`` step if needed
   * - A natural-language or automated workflow component
     - ``pycsamt.agents``
   * - A declarative workflow composition
     - ``pycsamt.pipeline``
   * - A command-line entry point
     - ``pycsamt.cli.commands``
   * - A specific external inversion engine interface
     - ``pycsamt.models.<engine>`` or ``pycsamt.inversion.backends``
   * - A scientific inversion abstraction independent of one engine
     - ``pycsamt.inversion``
   * - A format-specific reader or writer
     - ``pycsamt.seg``, ``pycsamt.jones``, ``pycsamt.zonge``, ``pycsamt.tdem``,
       or ``pycsamt.io``
   * - Site, station, profile, coordinate, or survey organization
     - ``pycsamt.site`` or ``pycsamt.metadata``
   * - A compatibility alias from older code
     - ``pycsamt.compat`` or a documented legacy module with warnings


In short
--------

The public pyCSAMT v2 API should be small at the top level, clear at package
boundaries, rich in documentation, conservative about compatibility, and
honest about optional dependencies.  Internal modules can evolve quickly, but
documented public imports are a promise to users building scientific workflows
on top of pyCSAMT.
