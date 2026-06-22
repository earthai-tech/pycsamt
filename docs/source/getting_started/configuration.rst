.. _getting-started-configuration:

Configuration
=============

pyCSAMT can be used with almost no configuration for a first survey, but a
real project usually needs a few choices to be explicit:

* where outputs are written;
* how progress and CLI output should look;
* whether dataframe-like results are wrapped in pyCSAMT view objects;
* which plot style should be used;
* which AI backend should be selected;
* whether AI-assisted agents may call an LLM provider;
* how much LLM cost is allowed in one session.

Most v2 configuration follows the same pattern:

.. code-block:: python
   :linenos:

   configure_something(...)
   reset_something()

   with PYCSAMT_SOMETHING.context(...):
       ...

Use global configuration for a notebook or application session.  Use context
managers when one workflow needs temporary settings.


Quick project setup
-------------------

This is a reasonable starting configuration for a local project.

.. code-block:: python
   :linenos:

   from pycsamt.api import (
       configure_api_view,
       configure_cli,
       configure_pipe,
       configure_style,
   )
   from pycsamt.backends import set_backend

   configure_api_view(backend="pycsamt")

   configure_cli(
       log__level=1,
       output__format="text",
       output__dir="results",
       build__n_jobs=4,
   )

   configure_pipe(
       output_root="results/pipeline",
       plot_dpi=200,
       plot_fmt="png",
       show_progress=True,
       on_step_error="warn",
   )

   configure_style(
       multiline__mode="gradient",
       multiline__base_color="blue",
       mt__xy__color="#003f88",
       mt__yx__color="#d62828",
   )

   # Optional: only needed for AI/model-zoo workflows.
   set_backend("auto")

This does not configure LLM agents.  Agents can still run deterministic steps
where supported, but LLM interpretation remains disabled until a provider key
is configured.


Configuration layers
--------------------

pyCSAMT has several configuration layers because different parts of the
package solve different problems.

.. list-table::
   :header-rows: 1
   :widths: 24 34 42

   * - Layer
     - Main object
     - What it controls
   * - API view
     - ``PYCSAMT_API_VIEW``
     - Whether dataframe-like results are returned as pyCSAMT view objects or
       raw pandas-style objects.
   * - CLI
     - ``PYCSAMT_CLI``
     - Verbosity, color, output format, output directory, cache, and job count
       for command-line workflows.
   * - Pipeline
     - ``PYCSAMT_PIPE``
     - Pipeline output directories, progress style, plot format, report format,
       and step error policy.
   * - Style
     - ``PYCSAMT_STYLE``
     - Plot colors, line styles, rose diagrams, component colors, and
       publication-style visual defaults.
   * - AI backend
     - ``pycsamt.backends``
     - PyTorch/TensorFlow backend selection for AI inversion and model-zoo
       workflows.
   * - Agents
     - ``AGENT_CONFIG``
     - LLM provider, model, API keys, pricing, and session budget.


Inspect current settings
------------------------

The global configuration objects are normal Python objects.  You can print
them or inspect their attributes.

.. code-block:: python
   :linenos:

   from pycsamt.api import PYCSAMT_CLI, PYCSAMT_PIPE
   from pycsamt.api.view import PYCSAMT_API_VIEW
   from pycsamt.backends import get_backend, list_backends

   print(PYCSAMT_CLI)
   print(PYCSAMT_PIPE)
   print(PYCSAMT_API_VIEW)
   print(get_backend())
   print(list_backends())

For agents:

.. code-block:: python
   :linenos:

   from pycsamt.agents import AGENT_CONFIG

   print(AGENT_CONFIG.info())


Temporary configuration
-----------------------

Context managers are useful when a single block of code needs different
settings without changing the rest of the session.

.. code-block:: python
   :linenos:

   from pycsamt.api import PYCSAMT_PIPE
   from pycsamt.pipeline import Pipeline, Step

   pipe = Pipeline([
       ("notch", Step("NR001", mains_hz=50)),
       ("band", Step("FREQ001")),
   ])

   with PYCSAMT_PIPE.context(
       output_root="results/publication",
       plot_dpi=300,
       plot_fmt="pdf",
       show_progress=False,
   ):
       result = pipe.run(sites)

   # The previous pipeline settings are restored here.


API view configuration
----------------------

Many public helpers return dataframe-like data.  The API view layer controls
whether those results are wrapped in pyCSAMT's lightweight view objects.

Default behavior uses pyCSAMT wrappers:

.. code-block:: python
   :linenos:

   from pycsamt.api import configure_api_view, read_edis

   configure_api_view(backend="pycsamt")
   survey = read_edis("data/willy/edis")

To receive raw pandas-style objects where applicable:

.. code-block:: python
   :linenos:

   from pycsamt.api import configure_api_view

   configure_api_view(backend="pandas")

Accepted values:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Value
     - Meaning
   * - ``"pycsamt"``, ``"api"``, ``"view"``, ``True``
     - Return pyCSAMT API view wrappers.
   * - ``"pandas"``, ``"raw"``, ``"none"``, ``False``
     - Return the raw dataframe-like object.
   * - callable
     - Call ``wrapper(data, **metadata)`` and return the custom wrapper result.

Environment variable:

.. code-block:: bash
   :linenos:

   export PYCSAMT_API_VIEW=pandas


CLI configuration
-----------------

The CLI configuration uses section-style keys with double underscores:

.. code-block:: python
   :linenos:

   from pycsamt.api import PYCSAMT_CLI, configure_cli, reset_cli

   configure_cli(
       log__level=1,
       log__color=True,
       output__format="json",
       output__dir="results/cli",
       output__overwrite=False,
       build__n_jobs=4,
       build__cache=True,
   )

   print(PYCSAMT_CLI.summary())

   reset_cli()

Important keys:

.. list-table::
   :header-rows: 1
   :widths: 30 28 42

   * - Key
     - Values
     - Meaning
   * - ``log__level``
     - ``0``, ``1``, ``2``
     - Quiet, info, or debug verbosity.
   * - ``log__color``
     - ``True`` or ``False``
     - Enable colored terminal output.
   * - ``output__format``
     - ``"text"``, ``"json"``, ``"csv"``
     - Preferred command output format.
   * - ``output__dir``
     - path-like
     - Default output directory.
   * - ``output__overwrite``
     - ``True`` or ``False``
     - Whether existing output may be overwritten.
   * - ``build__n_jobs``
     - integer >= 1
     - Number of jobs for CLI workflows that support parallel execution.
   * - ``build__cache``
     - ``True`` or ``False``
     - Enable command cache where supported.

Environment variables loaded at import time:

.. list-table::
   :header-rows: 1
   :widths: 34 66

   * - Variable
     - Meaning
   * - ``PYCSAMT_VERBOSE``
     - Integer verbosity level, usually ``0``, ``1``, or ``2``.
   * - ``PYCSAMT_NO_COLOR``
     - Any non-empty value disables color.
   * - ``PYCSAMT_OUTPUT``
     - Output format: ``text``, ``json``, or ``csv``.
   * - ``PYCSAMT_OUTPUT_DIR``
     - Default output directory.
   * - ``PYCSAMT_JOBS``
     - Integer job count.

Example shell setup:

.. code-block:: bash
   :linenos:

   export PYCSAMT_VERBOSE=1
   export PYCSAMT_OUTPUT=json
   export PYCSAMT_OUTPUT_DIR=results/cli
   export PYCSAMT_JOBS=4


Pipeline configuration
----------------------

Pipeline configuration controls the default behavior of
``pycsamt.pipeline.Pipeline.run``.

.. code-block:: python
   :linenos:

   from pycsamt.api import PYCSAMT_PIPE, configure_pipe, reset_pipe

   configure_pipe(
       output_root="results/pipeline",
       processed_subdir="processed",
       plots_subdir="plots",
       on_step_error="warn",
       save_intermediate=False,
       show_progress=True,
       progress_style="bar",
       plot_dpi=150,
       plot_fmt="png",
       report_formats=("html", "txt"),
   )

   print(PYCSAMT_PIPE)

   reset_pipe()

Key settings:

.. list-table::
   :header-rows: 1
   :widths: 28 32 40

   * - Setting
     - Default
     - Meaning
   * - ``output_root``
     - ``"pipe_results"``
     - Root directory used when ``Pipeline.run`` is called without ``outdir``.
   * - ``processed_subdir``
     - ``"processed"``
     - Subdirectory for processed EDI files.
   * - ``plots_subdir``
     - ``"plots"``
     - Subdirectory for QC and diagnostic figures.
   * - ``on_step_error``
     - ``"warn"``
     - Step failure policy: ``"raise"``, ``"warn"``, or ``"skip"``.
   * - ``save_intermediate``
     - ``False``
     - Save EDI snapshots after each step.
   * - ``show_progress``
     - ``True``
     - Show progress during pipeline execution.
   * - ``progress_style``
     - ``"bar"``
     - Progress style: ``"bar"``, ``"log"``, or ``"silent"``.
   * - ``plot_dpi``
     - ``150``
     - DPI for saved figures.
   * - ``plot_fmt``
     - ``"png"``
     - Figure format, for example ``"png"``, ``"pdf"``, or ``"svg"``.
   * - ``report_formats``
     - ``("html", "txt")``
     - Pipeline report formats.

Use stricter error behavior when validating a new workflow:

.. code-block:: python
   :linenos:

   from pycsamt.api import PYCSAMT_PIPE

   with PYCSAMT_PIPE.context(on_step_error="raise"):
       result = pipe.run(sites, outdir="results/debug")


Plot style configuration
------------------------

Plot style configuration keeps figures consistent across notebooks,
tutorials, reports, and agent-generated outputs.

Use named style presets where available:

.. code-block:: python
   :linenos:

   from pycsamt.api import use_style

   use_style("publication")

Override specific style attributes with double-underscore paths:

.. code-block:: python
   :linenos:

   from pycsamt.api import configure_style

   configure_style(
       multiline__mode="gradient",
       multiline__base_color="teal",
       multiline__lw=1.8,
       mt__xy__color="#003f88",
       mt__yx__color="#d62828",
       correction__before__color="#808080",
       correction__after__color="#005f73",
   )

Temporary style:

.. code-block:: python
   :linenos:

   from pycsamt.api import PYCSAMT_STYLE

   with PYCSAMT_STYLE.context("dark"):
       fig = plot_function(...)

The exact style sections are documented in the API reference, but the most
common sections are:

* ``multiline`` for multi-station and multi-profile lines;
* ``mt`` for MT component colors;
* ``rose`` for rose diagrams;
* ``correction`` for before/after correction figures;
* ``raw`` for diagnostic raw-data traces.


AI backend configuration
------------------------

AI inversion and model-zoo workflows can use PyTorch or TensorFlow.  Neither
backend is required for the base package.  Backend packages are loaded lazily
when AI functionality is used.

Inspect available backends:

.. code-block:: python
   :linenos:

   from pycsamt.backends import get_backend, list_backends

   print(list_backends())
   print(get_backend())

Select a backend for this session:

.. code-block:: python
   :linenos:

   from pycsamt.backends import set_backend

   set_backend("torch")
   # or
   set_backend("tensorflow")
   # or
   set_backend("auto")

Persist a backend choice to ``~/.pycsamt/config.json``:

.. code-block:: python
   :linenos:

   from pycsamt.backends import set_backend

   set_backend("torch", persist=True)

Environment variable:

.. code-block:: bash
   :linenos:

   export PYCSAMT_AI_BACKEND=torch

Backend resolution order:

1. Explicit ``set_backend(...)`` call in the current Python session.
2. ``PYCSAMT_AI_BACKEND`` environment variable.
3. ``~/.pycsamt/config.json`` key ``"ai_backend"``.
4. Auto-detection from installed frameworks.


Agent and LLM configuration
---------------------------

pyCSAMT agents can run deterministic scientific steps without an LLM when the
agent supports fallback behavior.  LLM configuration is only needed for
natural-language interpretation, report wording, and assistant-style workflow
guidance.

Basic setup:

.. code-block:: python
   :linenos:

   from pycsamt.agents import configure_agents

   configure_agents(
       provider="claude",
       api_key="sk-ant-...",
       model="claude-sonnet-4-6",
   )

Configure with environment variables instead of hard-coding keys:

.. code-block:: bash
   :linenos:

   export ANTHROPIC_API_KEY="sk-ant-..."
   export OPENAI_API_KEY="sk-..."
   export GOOGLE_API_KEY="AIza..."

Then select the provider in Python:

.. code-block:: python
   :linenos:

   from pycsamt.agents import configure_agents

   configure_agents(provider="openai")

Provider environment variables:

.. list-table::
   :header-rows: 1
   :widths: 28 72

   * - Provider
     - Variables checked
   * - ``"claude"``
     - ``ANTHROPIC_API_KEY``, ``PYCSAMT_CLAUDE_API_KEY``
   * - ``"openai"``
     - ``OPENAI_API_KEY``, ``PYCSAMT_OPENAI_API_KEY``
   * - ``"gemini"``
     - ``GOOGLE_API_KEY``, ``GOOGLE_GENERATIVEAI_API_KEY``,
       ``PYCSAMT_GEMINI_API_KEY``

Use a session budget before experimenting with LLM-assisted agents:

.. code-block:: python
   :linenos:

   from pycsamt.agents import AGENT_CONFIG

   AGENT_CONFIG.set_budget(usd=2.0)
   print(AGENT_CONFIG.remaining_usd)

Switch providers:

.. code-block:: python
   :linenos:

   from pycsamt.agents import AGENT_CONFIG

   AGENT_CONFIG.set_key("claude", "sk-ant-...")
   AGENT_CONFIG.set_key("openai", "sk-...")

   AGENT_CONFIG.switch("claude")
   # run Claude-assisted workflow

   AGENT_CONFIG.switch("openai")
   # run OpenAI-assisted workflow

Temporary LLM override:

.. code-block:: python
   :linenos:

   from pycsamt.agents import AGENT_CONFIG, DataQCAgent

   with AGENT_CONFIG.using(provider="gemini", api_key="AIza..."):
       result = DataQCAgent().execute({"path": "data/willy/edis"})

   # Original agent settings are restored here.

For the full agent configuration guide, see
:ref:`agents-llm-configuration`.


Configuration files
-------------------

Most day-to-day configuration is done in Python or environment variables.
The AI backend layer can also persist the selected backend to:

.. code-block:: text
   :linenos:

   ~/.pycsamt/config.json

Example content:

.. code-block:: json
   :linenos:

   {
     "ai_backend": "torch"
   }

Pipeline workflows may also be stored as YAML or JSON pipeline configuration
files.  Those files describe processing steps, not global runtime settings.
See the pipeline configuration guide for details:
:ref:`pipeline-configuration-files`.


Logging behavior
----------------

CLI verbosity is controlled by ``PYCSAMT_CLI`` and the ``PYCSAMT_VERBOSE``
environment variable.  Python logging is initialized by the package at import
time so command-line and library workflows have a consistent baseline.

For most users:

.. code-block:: bash
   :linenos:

   export PYCSAMT_VERBOSE=1

For Python scripts:

.. code-block:: python
   :linenos:

   from pycsamt.api import configure_cli

   configure_cli(log__level=2)

Use ``log__level=2`` for debugging configuration and parsing problems.  Use
``log__level=0`` for quiet batch runs.


Recommended setups
------------------

Notebook exploration:

.. code-block:: python
   :linenos:

   from pycsamt.api import configure_api_view, configure_pipe, use_style

   configure_api_view(backend="pycsamt")
   configure_pipe(output_root="notebook_results", show_progress=True)
   use_style("publication")

Batch processing:

.. code-block:: python
   :linenos:

   from pycsamt.api import configure_cli, configure_pipe

   configure_cli(log__level=1, output__format="json", build__n_jobs=8)
   configure_pipe(
       output_root="batch_results",
       progress_style="log",
       on_step_error="warn",
       save_intermediate=True,
   )

AI inversion:

.. code-block:: python
   :linenos:

   from pycsamt.backends import set_backend
   from pycsamt.agents import AGENT_CONFIG

   set_backend("auto")
   AGENT_CONFIG.set_budget(usd=5.0)

   # Add provider only when LLM explanation is needed.
   # AGENT_CONFIG.configure(provider="claude")

Report/publication figures:

.. code-block:: python
   :linenos:

   from pycsamt.api import configure_pipe, use_style

   use_style("publication")
   configure_pipe(plot_dpi=300, plot_fmt="pdf")


Reset configuration
-------------------

Reset helpers restore package defaults for the current Python session.

.. code-block:: python
   :linenos:

   from pycsamt.api import reset_api_view, reset_cli, reset_pipe, reset_style
   from pycsamt.agents import reset_agents

   reset_api_view()
   reset_cli()
   reset_pipe()
   reset_style()
   reset_agents()

Backend persistence is separate.  If you used ``set_backend(..., persist=True)``,
edit or remove ``~/.pycsamt/config.json`` to clear the saved backend choice.


Common mistakes
---------------

.. list-table::
   :header-rows: 1
   :widths: 36 64

   * - Problem
     - Fix
   * - Agent does not produce LLM interpretation
     - Configure a provider and key, or accept deterministic no-LLM fallback.
   * - AI inversion says no backend is installed
     - Install PyTorch or TensorFlow, then run ``set_backend("auto")``.
   * - Pipeline writes to an unexpected directory
     - Pass ``outdir=...`` to ``Pipeline.run`` or configure
       ``output_root``.
   * - CLI ignores environment changes
     - Set environment variables before importing pyCSAMT or call
       ``PYCSAMT_CLI.load_env()``.
   * - Plots are inconsistent across notebooks
     - Configure style once at the start of the session.
   * - Raw pandas objects are needed
     - Use ``configure_api_view(backend="pandas")``.


In short
--------

For a first project, configure output directories, progress behavior, and
plot style.  Add AI backend selection only when using AI inversion.  Add LLM
agent configuration only when natural-language interpretation or assistant
guidance is required.
