.. _development-docstring-style:

Docstring Style
===============

pyCSAMT uses NumPy-style docstrings for public Python APIs.  The goal is not
only pretty generated reference pages.  Good docstrings are part of the v2 API
contract: they explain scientific assumptions, units, array shapes, optional
dependencies, agent input/output schemas, and workflow reproducibility.

This page defines the expected style for public functions, classes, agents,
pipeline steps, CLI helpers, inversion objects, and compatibility aliases.


Why docstrings matter
---------------------

pyCSAMT users often work across notebooks, scripts, CLI workflows, agents, and
generated API pages.  A docstring may be read in any of these places:

* ``help(pycsamt.api.read_edis)`` in a Python shell;
* an IDE hover card;
* an autosummary-generated API page;
* a tutorial that links to the API reference;
* a developer review of a new public function;
* an agent or pipeline page explaining structured outputs.

For this reason, public docstrings must be complete enough to stand alone, but
short enough that the API reference remains scannable.


Relationship to the API policy
------------------------------

The :ref:`development-api-policy` defines which names are public.  This page
defines how those public names should be documented.

In short:

* public stable APIs need complete NumPy-style docstrings;
* public experimental APIs need complete docstrings plus a stability note;
* private helpers may have short implementation docstrings;
* deprecated aliases must document the replacement and removal plan.


General rules
-------------

Use these rules everywhere unless a more specific section below says otherwise.

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Rule
     - Requirement
   * - Start with one sentence
     - The first line should describe the object in plain language and fit on
       one line when reasonable.
   * - Use NumPy sections
     - Prefer recognized sections such as ``Parameters``, ``Returns``,
       ``Attributes``, ``Raises``, ``Warns``, ``Notes``, ``Examples``, and
       ``See Also``.
   * - Be explicit about science
     - Document units, coordinate systems, tensor conventions, frequency versus
       period, resistivity units, and array shapes.
   * - Document optional dependencies
     - Name optional packages required at execution time, not only at install
       time.
   * - Prefer keyword clarity
     - Public functions should document defaults and accepted strings.
   * - Keep examples runnable
     - Examples should use public imports and realistic small workflows.
   * - Avoid private import paths
     - Examples should not teach users to import from underscored modules.
   * - No custom top-level sections
     - Do not add headings such as ``Input Keys`` or ``Output Data Keys`` at the
       same level as ``Parameters``.  Put those details inside ``Notes``.


Recognized section order
------------------------

Use this order for most public docstrings.

.. code-block:: text
   :linenos:

   One-line summary.

   Optional extended summary.  Explain the scientific purpose, not the
   implementation line by line.

   Parameters
   ----------
   name : type
       Description.

   Returns
   -------
   type
       Description.

   Raises
   ------
   ErrorType
       When and why this is raised.

   Warns
   -----
   WarningType
       When and why this warning is emitted.

   See Also
   --------
   related_function : Short description.

   Notes
   -----
   Scientific assumptions, units, shapes, algorithms, or structured schemas.

   Examples
   --------
   >>> from pycsamt.api import read_edis
   >>> survey = read_edis("data/edis")

Not every docstring needs every section.  For a simple helper, a summary,
``Parameters``, ``Returns``, and a short example may be enough.


Sections to avoid
-----------------

numpydoc warns when it sees unknown sections.  Avoid custom top-level section
headings like these:

.. code-block:: text
   :linenos:

   Input Keys
   ----------
   ...

   Output Data Keys
   ----------------
   ...

   Resolution Rules
   ----------------
   ...

   Recognised Variables
   --------------------
   ...

Instead, put the same content under ``Notes`` using paragraphs, bullet lists,
or simple tables.

Preferred:

.. code-block:: text
   :linenos:

   Notes
   -----
   Input mapping keys:

   ``sites`` : Sites
       Survey object to process.
   ``path`` : str or path-like
       Directory or file used when ``sites`` is not supplied.

   Output data keys:

   ``qc_table`` : pandas.DataFrame
       Per-station quality metrics.
   ``n_flagged`` : int
       Number of stations that failed the configured thresholds.


Function docstrings
-------------------

A public function docstring must explain what the function does, what it
accepts, what it returns, and what can fail.

Template:

.. code-block:: python
   :linenos:

   def read_edis(
       sources,
       *,
       recursive=True,
       strict=False,
       on_dup="replace",
       progress="auto",
       verbose=0,
   ):
       """Read EDI files and return a public survey view.

       Parameters
       ----------
       sources : path-like, sequence of path-like, or file-like
           EDI file, directory, glob-compatible collection, or explicit
           sequence of EDI paths.
       recursive : bool, default=True
           If True, search directories recursively.
       strict : bool, default=False
           If True, raise when a file cannot be parsed.  If False, keep
           successful files and record parse errors on the parser metadata.
       on_dup : {"replace", "keep"}, default="replace"
           Policy used when several files resolve to the same station name.
       progress : bool or {"auto"}, default="auto"
           Whether to display progress while reading many files.
       verbose : int, default=0
           Verbosity level passed to the low-level parser.

       Returns
       -------
       APISurvey
           Public survey wrapper containing an EDI collection and parser
           metadata.

       Raises
       ------
       OSError
           If an input path cannot be read.
       ValueError
           If ``on_dup`` is not a recognized duplicate policy.

       Notes
       -----
       The returned object is the public view layer used by notebooks, CLI
       commands, and higher-level workflows.  Low-level parser details remain
       available through metadata when needed.

       Examples
       --------
       >>> from pycsamt.api import read_edis
       >>> survey = read_edis("data/willy/edis")
       >>> survey.name
       'edi_survey'
       """


Parameter style
---------------

Follow these conventions in ``Parameters``.

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Situation
     - Style
   * - Default values
     - Put the default in the type line: ``verbose : int, default=0``.
   * - Accepted strings
     - Use braces: ``method : {"composite", "snr"}, default="composite"``.
   * - Optional values
     - Use ``type or None`` or ``optional`` consistently.
   * - Path-like values
     - Prefer ``path-like`` unless the function truly requires ``pathlib.Path``.
   * - Arrays
     - Include shape and units: ``rho : array-like of shape (n_freq, n_site)``.
   * - Dictionaries
     - Document important keys under the parameter description or in ``Notes``.
   * - Forwarded kwargs
     - Explain where they go: ``**kwargs : dict`` followed by the callee.

Good examples:

.. code-block:: text
   :linenos:

   period_range : tuple of float, optional
       Two-element ``(min_period, max_period)`` window in seconds.
   components : {"xy", "yx", "det", "ssq"}, default="det"
       Impedance component or invariant used to compute apparent resistivity.
   station_spacing : float, optional
       Station spacing in metres.  Required when station coordinates are not
       available.


Return style
------------

Prefer structured return documentation.  Avoid vague lines such as
``Returns result``.

Examples:

.. code-block:: text
   :linenos:

   Returns
   -------
   APISurvey
       Survey wrapper containing parsed EDI files, station metadata, and parser
       diagnostics.

.. code-block:: text
   :linenos:

   Returns
   -------
   rho : ndarray of shape (n_depth, n_station)
       Estimated resistivity section in ohm-m.
   depth : ndarray of shape (n_depth,)
       Cell-centre depths in metres.

For multiple named return values, either return a dataclass/result object or
document each output name.  Do not return long positional tuples from new
public APIs unless there is a compatibility reason.


Units, shapes, and conventions
------------------------------

Scientific docstrings must record the assumptions that affect interpretation.

Document:

* frequency in hertz or period in seconds;
* apparent resistivity in ohm-m;
* phase in degrees or radians;
* distance, elevation, and depth in metres unless otherwise stated;
* longitude/latitude in decimal degrees;
* projected coordinates and EPSG/UTM assumptions when applicable;
* impedance tensor component order;
* shape names such as ``(n_freq, n_station)`` or ``(n_depth, n_x)``;
* whether arrays are linear scale, log10 scale, or normalized.

Example:

.. code-block:: text
   :linenos:

   Notes
   -----
   ``rho`` is expected in ohm-m on linear scale with shape
   ``(n_frequency, n_station)``.  ``phase`` is expected in degrees with the
   same shape.  Frequencies are sorted from high to low by the reader, but
   plotting functions may display period increasing downward.


Class docstrings
----------------

Class docstrings should document construction and user-facing state.  Do not
list every private attribute.

Template:

.. code-block:: python
   :linenos:

   class ForwardConfig:
       """Configuration for forward electromagnetic modelling.

       Parameters
       ----------
       dimensionality : {1, 2, 3}, default=2
           Model dimensionality.
       frequencies : array-like of float, optional
           Frequencies in hertz.
       resistivity_background : float, default=100.0
           Background resistivity in ohm-m.

       Attributes
       ----------
       dimensionality : int
           Validated model dimensionality.
       frequencies : ndarray or None
           Frequencies used by the forward solver.

       Notes
       -----
       This object stores configuration only.  It does not run a solver or
       import optional backend packages.
       """

Use ``Attributes`` for values users are expected to inspect or modify.  Use
``Notes`` for implementation details, algorithms, and scientific assumptions.


Dataclass docstrings
--------------------

For dataclasses, document constructor fields in ``Parameters`` when users
instantiate the class directly.  Document computed or post-init fields in
``Attributes``.

Example:

.. code-block:: python
   :linenos:

   @dataclass
   class StepResult:
       """Record produced after one pipeline step runs.

       Parameters
       ----------
       step_idx : int
           One-based position of the step in the pipeline.
       step_name : str
           User-supplied step label.
       step_code : str
           Stable pipeline registry code, such as ``"NR001"``.

       Attributes
       ----------
       error : Exception or None
           Exception captured when the pipeline continued after a failed step.
       """


Agent docstrings
----------------

Agents need richer docstrings because they are both Python objects and workflow
building blocks.  They should document constructor parameters, execution input
mapping, output mapping, LLM behavior, and examples.

Use recognized numpydoc sections.  Put agent schemas inside ``Notes`` rather
than custom top-level ``Input Keys`` sections.

Template:

.. code-block:: python
   :linenos:

   class DataQCAgent(BaseAgent):
       """Run data quality control on an MT/AMT/CSAMT survey.

       Parameters
       ----------
       api_key : str, optional
           Provider API key used for optional LLM interpretation.
       model : str, optional
           Provider model name.  If omitted, the global agent configuration is
           used.
       llm_provider : {"claude", "openai", "gemini"}, default="claude"
           LLM provider used when interpretation is enabled.
       method : {"composite", "presence", "snr", "spatial"}, default="composite"
           Confidence scoring method.
       min_frac_ok : float, default=0.6
           Minimum fraction of acceptable frequencies required for a station to
           pass QC.

       Returns
       -------
       AgentResult
           Structured result.  The ``data`` mapping contains QC tables,
           flagged stations, figures, and saved figure paths.

       Notes
       -----
       Execution input mapping:

       ``sites`` : Sites, optional
           Survey object to evaluate.
       ``path`` : path-like, optional
           EDI directory or file used when ``sites`` is not supplied.
       ``output_dir`` : path-like, optional
           Directory where figures are written.
       ``period_range`` : tuple of float, optional
           Period window in seconds.

       Output ``data`` mapping:

       ``qc_table`` : pandas.DataFrame
           Per-station metrics.
       ``flags`` : pandas.DataFrame
           Pass/fail flags per station.
       ``n_flagged`` : int
           Number of stations that failed configured thresholds.
       ``figures`` : dict
           Matplotlib figures keyed by plot name.

       LLM interpretation is optional.  When no provider key is configured,
       the agent still runs deterministic QC and leaves
       ``llm_interpretation`` empty.

       Examples
       --------
       >>> from pycsamt.agents import DataQCAgent
       >>> agent = DataQCAgent(method="composite")
       >>> result = agent.execute({
       ...     "path": "data/willy/edis",
       ...     "output_dir": "results/willy_qc",
       ... })
       >>> result.success
       True
       """

Agent examples should use public imports:

.. code-block:: python
   :linenos:

   from pycsamt.agents import DataQCAgent

   result = DataQCAgent().execute({
       "path": "data/willy/edis",
       "output_dir": "results/qc",
   })

   if result.success:
       print(result.data["n_flagged"])


Coordinator and orchestrator docstrings
---------------------------------------

Coordinator and orchestrator docstrings should describe workflow behavior, not
only constructor arguments.

Document:

* how steps are selected or registered;
* what ``dry_run`` returns;
* how failed steps are represented;
* how LLM routing behaves when no key is configured;
* which keys are required in the input mapping;
* whether output is deterministic.

Example ``Notes`` content:

.. code-block:: text
   :linenos:

   Notes
   -----
   ``dry_run=True`` returns a plan without executing workflow steps.  The plan
   includes the selected workflow name, step order, and input keys expected by
   each step.  No LLM call is required for dry-run routing when
   ``workflow_type`` is supplied explicitly.


Pipeline docstrings
-------------------

Pipeline docstrings should make serialization and reproducibility clear.

For ``Step``:

* document accepted registry identifiers;
* mention that parameter overrides are merged over registry defaults;
* document serialization with ``to_dict`` or pipeline config files.

For ``Pipeline``:

* document accepted step declarations;
* document input survey type;
* document output directory behavior;
* document error policy;
* document generated plots/reports;
* include a short example using ``Pipeline`` and ``Step`` from
  ``pycsamt.pipeline``.

Example:

.. code-block:: python
   :linenos:

   from pycsamt.pipeline import Pipeline, Step

   pipe = Pipeline([
       ("notch", Step("NR001", mains_hz=50)),
       ("band", Step("FREQ001", fmin=1e-3, fmax=1.0)),
   ])
   result = pipe.run(sites, outdir="results/willy")


Pipeline registry docstrings
----------------------------

Registry objects such as step specifications should document user-facing
fields.  Step codes are public once documented.

Preferred style:

.. code-block:: text
   :linenos:

   Parameters
   ----------
   code : str
       Stable registry code, for example ``"NR001"``.
   name : str
       Snake-case registry name.
   category : str
       Processing category used by discovery helpers.
   defaults : dict
       JSON/YAML-serializable default parameters.

   Notes
   -----
   Step codes are part of the public pipeline API.  Do not rename a documented
   code without adding a deprecated alias.


Inversion and model docstrings
------------------------------

Inversion APIs must document scientific and runtime assumptions carefully.

Include:

* dimensionality;
* supported survey type;
* required input object or file format;
* units and shapes;
* mesh/model assumptions;
* external executable or backend requirements;
* output files;
* reproducibility metadata.

Example:

.. code-block:: text
   :linenos:

   Notes
   -----
   The model grid uses metres in projected coordinates.  Resistivity values are
   stored in ohm-m on linear scale unless ``log10=True`` is passed.  The
   Occam2D executable is not invoked by this builder; it only writes validated
   input files.


AI model docstrings
-------------------

AI model and training docstrings should be precise about tensors and backend
requirements.

Document:

* backend: PyTorch, TensorFlow, or backend-agnostic;
* input tensor shape and channel order;
* target tensor shape;
* normalization/scaling;
* stochastic behavior and seeds;
* checkpoint format;
* GPU/CPU assumptions when relevant.

Example:

.. code-block:: text
   :linenos:

   Parameters
   ----------
   X : array-like of shape (n_sample, n_frequency, n_station, n_channel)
       Normalized input features.  Channel order is ``rho_xy``, ``phase_xy``,
       ``rho_yx``, ``phase_yx``.
   y : array-like of shape (n_sample, n_depth, n_station)
       Target resistivity sections in log10 ohm-m.

   Notes
   -----
   This class imports PyTorch lazily when the model is constructed.  Importing
   :mod:`pycsamt.ai` does not require PyTorch.


CLI helper docstrings
---------------------

CLI command functions often appear in generated docs and tests.  Their
docstrings should describe user behavior, not Click internals.

Good:

.. code-block:: python
   :linenos:

   def run(output_dir, verbose):
       """Run the configured inversion workflow from the command line."""

For reusable Click option decorators and parameter types, document:

* accepted CLI spelling;
* Python value returned by Click;
* validation behavior;
* examples if parsing is non-trivial.


Optional dependency docstrings
------------------------------

When a public API uses optional packages, state when the dependency is needed.

Example:

.. code-block:: text
   :linenos:

   Notes
   -----
   This function imports ``geopandas`` only when GeoPackage export is
   requested.  CSV and JSON export do not require GIS extras.

Do not claim that an optional dependency is required to import pyCSAMT unless
that is actually true.


Warnings, errors, and deprecations
----------------------------------

Use ``Raises`` for exceptions and ``Warns`` for runtime warnings.

Example:

.. code-block:: text
   :linenos:

   Raises
   ------
   ImportError
       If the selected backend requires an optional package that is not
       installed.
   ValueError
       If ``method`` is not one of the supported correction methods.

   Warns
   -----
   FutureWarning
       If a deprecated argument alias is used.

Deprecated APIs must include the replacement and planned removal in the
docstring:

.. code-block:: python
   :linenos:

   def old_reader(path):
       """Read EDI data using the legacy reader.

       .. deprecated:: 2.0.0
          Use :func:`pycsamt.api.read_edis` instead.  This alias is planned
          for removal in pyCSAMT 2.2.0.
       """


Examples style
--------------

Examples should be short, public, and realistic.

Use doctest prompts for small expressions:

.. code-block:: text
   :linenos:

   Examples
   --------
   >>> from pycsamt.pipeline import Step
   >>> Step("NR001", mains_hz=50).to_dict()["code"]
   'NR001'

Use narrative examples in user guides when setup is too heavy for doctest.
Inside docstrings, avoid examples that require large data downloads, API keys,
or long-running inversion jobs.

For LLM examples, show both configured and no-LLM behavior when relevant:

.. code-block:: text
   :linenos:

   Examples
   --------
   >>> from pycsamt.agents import ContextInputAgent
   >>> agent = ContextInputAgent()  # no API key: deterministic fallback
   >>> result = agent.execute({"request": "Load EDIs and run QC"})
   >>> "workflow" in result.data["config"]
   True


Cross-references
----------------

Use Sphinx roles for public names when the target is stable:

.. code-block:: text
   :linenos:

   See Also
   --------
   pycsamt.api.read_edis : Read EDI files into an API survey view.
   pycsamt.pipeline.Pipeline : Declarative processing workflow.
   pycsamt.agents.WorkflowOrchestratorAgent : Route natural-language requests
       to agent workflows.

In prose, prefer:

.. code-block:: text
   :linenos:

   The returned object is compatible with
   :class:`pycsamt.pipeline.Pipeline`.

Avoid fragile references to private helpers unless the docstring is for that
private helper.


Module docstrings
-----------------

Public modules should start with a short module docstring explaining the module
purpose and main entry points.

Good module docstring:

.. code-block:: python
   :linenos:

   """Public readers for EDI files and site collections.

   This module provides lightweight facade functions used by notebooks, CLI
   commands, and application views.  Low-level parser classes remain in the
   format-specific packages.
   """

Private implementation modules may use shorter module docstrings, but should
still explain non-obvious architecture.


Private helper docstrings
-------------------------

Private helpers do not need full NumPy-style documentation.  A short docstring
is enough when the function is local and simple.

Example:

.. code-block:: python
   :linenos:

   def _to_figure(obj):
       """Return a matplotlib Figure from a Figure-like or Axes-like object."""

Use fuller documentation for private helpers only when the algorithm is
complex or the helper is widely used internally.


Style details
-------------

Follow these formatting details to keep generated docs consistent.

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Item
     - Convention
   * - Line length
     - Keep source lines readable, usually below 88-96 characters.
   * - Markup
     - Use double backticks for literals and Sphinx roles for references.
   * - Bullets
     - Use bullets in ``Notes`` for readable lists.
   * - Tables
     - Avoid complex reST tables inside docstrings.  They are fragile in
       autosummary output.
   * - Unicode
     - Prefer ASCII in new code docstrings unless the symbol is scientifically
       important or already used nearby.
   * - Ellipses
     - Use ``...`` in examples, not the single Unicode ellipsis.
   * - External links
     - Put long links in user guides or references pages rather than API
       docstrings.


Common fixes for current warnings
---------------------------------

When Sphinx emits numpydoc warnings, the fix is usually one of these.

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Warning pattern
     - Fix
   * - ``Unknown section Input Keys``
     - Move the content under ``Notes`` as ``Execution input mapping``.
   * - ``Unknown section Output Data Keys``
     - Move the content under ``Notes`` as ``Output data mapping``.
   * - ``Unknown section Resolution Rules``
     - Move the content under ``Notes`` or ``Examples``.
   * - ``Unknown section Recognised Variables``
     - Move the variable list under ``Notes``.
   * - ``potentially wrong underline length``
     - Make the underline exactly as long as the heading text, or remove the
       custom heading.
   * - ``Citation is not referenced``
     - Move citations to the project references page or cite them explicitly in
       prose.


Review checklist
----------------

Use this checklist before exposing a new public API.

.. code-block:: text
   :linenos:

   [ ] First line is a clear one-sentence summary.
   [ ] Parameters include types, defaults, accepted strings, units, and shapes.
   [ ] Returns section names the structured result and important fields.
   [ ] Raises and Warns sections describe user-visible failures.
   [ ] Notes document scientific assumptions and optional dependencies.
   [ ] Agent input/output mappings are under Notes, not custom sections.
   [ ] Examples use public imports and small realistic workflows.
   [ ] Deprecated APIs include replacement and planned removal.
   [ ] Cross-references point to stable public names.
   [ ] Sphinx builds without new numpydoc warnings from this docstring.


Minimal examples by API type
----------------------------

Function:

.. code-block:: python
   :linenos:

   def estimate_skin_depth(rho, frequency):
       """Estimate electromagnetic skin depth.

       Parameters
       ----------
       rho : float or array-like
           Resistivity in ohm-m.
       frequency : float or array-like
           Frequency in hertz.

       Returns
       -------
       float or ndarray
           Skin depth in metres.
       """

Class:

.. code-block:: python
   :linenos:

   class SurveyMeta:
       """Metadata describing a survey acquisition.

       Parameters
       ----------
       name : str
           Survey name.
       crs : str, optional
           Coordinate reference system, such as ``"EPSG:32629"``.
       """

Agent:

.. code-block:: python
   :linenos:

   class ReportAgent(BaseAgent):
       """Assemble workflow outputs into Markdown, HTML, or PDF reports.

       Parameters
       ----------
       format : {"markdown", "html", "pdf"}, default="html"
           Report format to generate.

       Returns
       -------
       AgentResult
           Structured result containing report paths and collected warnings.

       Notes
       -----
       Execution input mapping:

       ``results`` : mapping
           Workflow results to summarize.
       ``output_dir`` : path-like
           Directory where report files are written.
       """

Pipeline step:

.. code-block:: python
   :linenos:

   class Step:
       """Configured pipeline step.

       Parameters
       ----------
       code_or_name : str
           Stable registry code or snake-case registry name.
       **params : dict
           Parameter overrides merged over registry defaults.

       Returns
       -------
       Step
           Configured step ready to place in a pipeline.
       """


In short
--------

A pyCSAMT docstring should answer four questions quickly:

* What does this object do?
* What inputs does it require, with which units and shapes?
* What structured result does it return?
* What assumptions, optional dependencies, warnings, or reproducibility details
  must the user know?

When those answers are clear, the generated API reference becomes useful
documentation rather than a mechanical listing of names.
