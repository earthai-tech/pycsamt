.. _api-configuration:

Configuration
=============

Every runtime behaviour of pyCSAMT — where outputs land, how figures look,
what tables are returned, how the CLI logs, how much an AI agent may spend —
is controlled from :mod:`pycsamt.api` through one repeated pattern.  This
page walks through that pattern and gives a worked example for each family.

For a "set up my first session" walkthrough, see
:doc:`../getting_started/configuration`.  This page is the systematic tour.

The Dotted-Path Convention
--------------------------

``configure_*`` functions accept keyword arguments whose double underscores
descend into nested settings, so ``configure_style(mt__xy__color="#003f88")``
sets ``style.mt.xy.color``.  Every family also exposes its live singleton —
print it to inspect the current state:

.. code-block:: python

   from pycsamt.api import PYCSAMT_STYLE, PYCSAMT_PIPE

   print(PYCSAMT_PIPE)     # current pipeline settings
   print(PYCSAMT_STYLE)    # current style settings

View Layer
----------

Decide what dataframe-returning functions give you when ``api=True``:

.. code-block:: python

   from pycsamt.api import configure_api_view, reset_api_view

   configure_api_view(backend="pycsamt")   # APIFrame / APIResult (default)
   configure_api_view(backend="pandas")    # plain DataFrames everywhere
   reset_api_view()

The full story — ``APIFrame``, multi-table results, custom wrappers — is on
the :doc:`views` page.

Pipeline Outputs
----------------

Control where pipeline runs write results and how they report progress:

.. code-block:: python

   from pycsamt.api import configure_pipe

   configure_pipe(
       output_root="results/pipeline",
       plot_dpi=200,
       plot_fmt="png",
       show_progress=True,
       on_step_error="warn",
   )

Batch runs often prefer log-style progress and intermediate saves:

.. code-block:: python

   configure_pipe(
       output_root="batch_results",
       progress_style="log",
       on_step_error="warn",
       save_intermediate=True,
   )

CLI Defaults
------------

The same settings the ``pycsamt`` command reads from the terminal can be
pre-configured in Python:

.. code-block:: python

   from pycsamt.api import configure_cli

   configure_cli(
       log__level=1,
       output__format="text",
       output__dir="results",
       build__n_jobs=4,
   )

Use ``log__level=0`` for quiet batch runs, ``output__format="json"`` for
machine-readable output.

Plot Styles
-----------

Named presets cover the common cases; dotted paths tune individual elements:

.. code-block:: python

   from pycsamt.api import use_style, configure_style

   use_style("publication")        # or "pycsamt" (default), "dark"

   configure_style(
       multiline__mode="gradient",
       multiline__base_color="blue",
       mt__xy__color="#003f88",
       mt__yx__color="#d62828",
   )

Figure Output
-------------

Global saving defaults apply to every figure pyCSAMT writes:

.. code-block:: python

   from pycsamt.api import set_dpi, set_fmt, set_savedir, save_fig

   set_dpi(300)                    # 150 screen, 300 print
   set_fmt("png", "pdf")           # save every figure in both formats
   set_savedir("figures/")

   paths = save_fig(fig, "response_S17")   # honours the settings above

View Controls
-------------

Axis conventions for apparent-resistivity and phase views, including phase
wrapping and frequency-axis direction:

.. code-block:: python

   from pycsamt.api import configure_control, reset_control

   configure_control(...)          # dotted-path arguments, see reference
   reset_control()

Sections, Stations, Interpretation, Topography
----------------------------------------------

The remaining plot families follow the identical pattern:

.. code-block:: python

   from pycsamt.api import (
       configure_section,             # resistivity-section figures
       configure_station_rendering,   # station map markers and axes
       configure_interp, use_interp,  # hydrogeological styles (with presets)
       configure_topo,                # topography and y-axis conventions
   )

See :doc:`../api/api` for every accepted dotted path.

Agents
------

Cap what AI-assisted workflows may spend, and pick the LLM provider only
when needed:

.. code-block:: python

   from pycsamt.agents import AGENT_CONFIG

   AGENT_CONFIG.set_budget(usd=5.0)
   # AGENT_CONFIG.configure(provider="claude")   # enable LLM explanations

Environment Variables
---------------------

Settings can be fixed before Python starts, which is useful for CI and
batch environments:

.. code-block:: bash

   PYCSAMT_API_VIEW=pandas python workflow.py
   PYCSAMT_API_VIEW=pycsamt python workflow.py

Recommended Setups
------------------

**Notebook exploration**

.. code-block:: python

   from pycsamt.api import configure_api_view, configure_pipe, use_style

   configure_api_view(backend="pycsamt")
   configure_pipe(output_root="notebook_results", show_progress=True)
   use_style("publication")

**Publication figures**

.. code-block:: python

   from pycsamt.api import configure_pipe, use_style

   use_style("publication")
   configure_pipe(plot_dpi=300, plot_fmt="pdf")

Reset Everything
----------------

Each family has a ``reset_*`` helper; together they restore a clean session:

.. code-block:: python

   from pycsamt.api import reset_api_view, reset_cli, reset_pipe, reset_style
   from pycsamt.agents import reset_agents

   reset_api_view()
   reset_cli()
   reset_pipe()
   reset_style()
   reset_agents()
