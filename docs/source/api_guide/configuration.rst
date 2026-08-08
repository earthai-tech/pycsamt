.. _api-configuration:

Configuration
=============

Every runtime behaviour of pyCSAMT -- where outputs land, how figures look,
what tables are returned, how the CLI logs, how much an AI agent may spend --
is controlled from :mod:`pycsamt.api` through one repeated pattern (see
:doc:`overview` if that pattern is new). This page is the systematic tour:
one real, runnable example per family, then a pointer to that family's own
page for the full depth.

For a "set up my first session" walkthrough instead, see
:doc:`../getting_started/configuration`.

The Dotted-Path Convention
--------------------------

``configure_*`` functions accept keyword arguments whose double underscores
descend into nested settings, so ``configure_style(mt__xy__color="#003f88")``
sets ``style.mt.xy.color``. Every family also exposes its live singleton --
print it to inspect the current state:

.. code-block:: pycon

   >>> from pycsamt.api import PYCSAMT_STYLE, PYCSAMT_PIPE

   >>> print(PYCSAMT_PIPE)
   PipelineAPIConfig
     output_root: 'pipe_results'
     processed_subdir: 'processed'
     plots_subdir: 'plots'
     on_step_error: 'warn'
     save_intermediate: False
     show_progress: True
     progress_style: 'bar'
     repr_width: 80
     plot_dpi: 150
     plot_fmt: 'png'
     report_formats: ('html', 'txt')

View Layer
----------

Decide what dataframe-returning functions give you when ``api=True``:

.. code-block:: pycon

   >>> from pycsamt.api import configure_api_view, reset_api_view, PYCSAMT_API_VIEW

   >>> print(PYCSAMT_API_VIEW)
   APIViewConfig(backend='pycsamt')
   >>> configure_api_view(backend="pandas")   # plain DataFrames everywhere
   >>> print(PYCSAMT_API_VIEW)
   APIViewConfig(backend='pandas')
   >>> reset_api_view()

The full story -- ``APIFrame``, multi-table results, custom wrappers -- is
on the :doc:`views` page.

Pipeline Outputs
----------------

Control where pipeline runs write results and how they report progress:

.. code-block:: pycon

   >>> from pycsamt.api import configure_pipe, reset_pipe

   >>> configure_pipe(
   ...     output_root="results/pipeline",
   ...     plot_dpi=200,
   ...     plot_fmt="png",
   ...     show_progress=True,
   ...     on_step_error="warn",
   ... )
   >>> print(PYCSAMT_PIPE)
   PipelineAPIConfig
     output_root: 'results/pipeline'
     processed_subdir: 'processed'
     plots_subdir: 'plots'
     on_step_error: 'warn'
     save_intermediate: False
     show_progress: True
     progress_style: 'bar'
     repr_width: 80
     plot_dpi: 200
     plot_fmt: 'png'
     report_formats: ('html', 'txt')

   >>> reset_pipe()

Batch runs often prefer log-style progress and intermediate saves:
``configure_pipe(progress_style="log", save_intermediate=True)``.

Site Ordering
-------------

Set the station-ordering policy once for the current Python process. Calls
that normalize their inputs through :func:`pycsamt.emtools.ensure_sites`, and
direct calls to :meth:`pycsamt.site.Sites.ordered` without a ``by`` argument,
then use the same policy:

.. code-block:: pycon

   >>> from pathlib import Path
   >>> from pycsamt.api import configure_ordering, PYCSAMT_ORDERING, reset_ordering
   >>> from pycsamt.emtools import ensure_sites

   >>> _ = configure_ordering(mode="auto")
   >>> sites = ensure_sites(Path("data/AMT/WILLY_DATA/L18PLT"))
   >>> sites.ordering["applied"], sites.ordering["n_sites"]
   ('chainage', 28)

   >>> reset_ordering()

``auto`` is the recommended default for survey lines. It converts latitude
and longitude to local metre coordinates, finds the principal profile axis,
and sorts stations by projected chainage -- applying that spatial order only
when the coordinates describe a credible approximately straight line (28 of
28 WILLY stations qualified above). If coordinates are missing or fail the
geometry checks, input order is preserved instead of guessing from station
names or from only one coordinate component.

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Mode
     - Behaviour
   * - ``auto``
     - Use validated coordinate-derived chainage; otherwise preserve input
       order.
   * - ``chainage``
     - Force projection along the coordinate-derived profile axis. Sites
       without usable coordinates remain at the end in their input order.
   * - ``input``
     - Preserve the order received from the loader or caller.
   * - ``station``
     - Natural numeric station-name order, for example ``S2`` before ``S10``.
   * - ``latitude``
     - Sort by latitude only.
   * - ``longitude``
     - Sort by longitude only.

The conservative acceptance thresholds for ``auto`` can be adjusted for a
known survey geometry
(``configure_ordering(mode="auto", min_linearity=0.95, max_cross_track_ratio=0.15)``),
an explicit per-call strategy remains authoritative and does not change the
global setting (``ensure_sites(path, order_by="input")``), and
``PYCSAMT_ORDERING.context(mode="station")`` scopes an override to one block.
The configuration is process-local: set it near the start of each script,
notebook kernel, worker process, or application startup, and restore package
defaults with :func:`~pycsamt.api.reset_ordering`.

CLI Defaults
------------

The same settings the ``pycsamt`` command reads from the terminal can be
pre-configured in Python:

.. code-block:: pycon

   >>> from pycsamt.api import configure_cli, PYCSAMT_CLI, reset_cli

   >>> configure_cli(
   ...     log__level=1,
   ...     output__format="text",
   ...     output__dir="results",
   ...     build__n_jobs=4,
   ... )
   >>> print(PYCSAMT_CLI)
   PyCSAMTCLI
     log.level     = 1  (info)
     log.color     = True
     log.file      = None
     output.format = 'text'
     output.dir    = 'results'
     output.overwrite = False
     build.n_jobs  = 4
     build.cache   = True
     build.cache_dir = None

   >>> reset_cli()

Use ``log__level=0`` for quiet batch runs, ``output__format="json"`` for
machine-readable output.

Plot Styles
-----------

Named presets cover the common cases; dotted paths tune individual elements:

.. code-block:: pycon

   >>> from pycsamt.api import use_style, configure_style, reset_style

   >>> use_style("publication")        # or "pycsamt" (default), "dark", "modem"
   >>> configure_style(mt__xy__color="#003f88", mt__yx__color="#d62828")
   >>> PYCSAMT_STYLE.mt.xy.color
   '#003f88'
   >>> reset_style()

MT component colours, multiline gradients, correction pairs, raw-data
style, phase-tensor ellipses, and rose diagrams -- plus the presets'
compounding behaviour and the rose functions' ``style="pycsamt"`` literal
gotcha -- are covered in depth on :doc:`style`.

Figure Output
-------------

Global saving defaults apply to every figure pyCSAMT writes:

.. code-block:: pycon

   >>> from pycsamt.api import set_dpi, set_fmt, set_savedir, PLOT_CONFIG, reset_plot_config

   >>> set_dpi(300)                    # 150 screen, 300 print
   >>> set_fmt("png", "pdf")           # save every figure in both formats
   >>> set_savedir("figures/")
   >>> print(PLOT_CONFIG)
   PlotConfig
     fmt              = ['png', 'pdf']
     resolved formats = ['png', 'pdf']
     base_fmt         = 'png'
     dpi              = 300
     bbox_inches      = 'tight'
     transparent      = False
     facecolor        = 'white'
     savedir          = 'figures/'
     close_after_save = False
     verbose          = True

   >>> reset_plot_config()

``save_fig(fig, "response_S17")`` honours the settings above -- and returns
the paths actually written, one per configured format.

View Controls
-------------

Apparent-resistivity scale, phase wrapping, and the frequency/period axis
convention -- the only family with no named presets:

.. code-block:: pycon

   >>> from pycsamt.api import configure_control, PYCSAMT_CONTROL, reset_control

   >>> configure_control(rho__view="linear", phase__range=(0.0, 360.0))
   >>> PYCSAMT_CONTROL.rho.view, PYCSAMT_CONTROL.phase.range
   ('linear', (0.0, 360.0))
   >>> reset_control()

The log10-vs-linear tradeoff, correct error propagation for either, and why
``use_log_scale()`` must be checked rather than assumed from the view name,
are covered on :doc:`view_controls`.

Section Layout
--------------

Figure sizing (fixed or data-aware ``"dynamic"``), axis direction, and
topography-awareness for every station-by-depth or station-by-period plot:

.. code-block:: pycon

   >>> from pycsamt.api import configure_section, PYCSAMT_SECTION, reset_section

   >>> configure_section(publication__colorbar__max_ticks=3)
   >>> PYCSAMT_SECTION.publication.colorbar.max_ticks
   3
   >>> reset_section()

Six presets (``"pseudosection"``, ``"inversion"``, ``"publication"``,
``"compact"``, ``"dashboard"``, ``"dynamic"``), the dynamic-sizing clamp
math, and the topography gate are covered on :doc:`section`.

Station Rendering
-----------------

Tick marks, adaptive label thinning, and marker glyphs along a profile's
station axis:

.. code-block:: pycon

   >>> from pycsamt.api import configure_station_rendering, PYCSAMT_STATION_RENDERING, reset_station_rendering

   >>> configure_station_rendering(inversion__marker__facecolor="crimson")
   >>> PYCSAMT_STATION_RENDERING.inversion.marker.facecolor
   'crimson'
   >>> reset_station_rendering()

The three presets, the "nice step" adaptive label-thinning algorithm, and
the terrain-following ``topo_elev=`` marker mode are covered on
:doc:`station_rendering`.

Interpretation
--------------

Colours and figure geometry for every hydrogeophysical plot -- sections,
water-table profiles, uncertainty panels:

.. code-block:: pycon

   >>> from pycsamt.api import use_interp, configure_interp, PYCSAMT_INTERP, reset_interp

   >>> use_interp("accessible")        # or "default", "publication", "dark"
   >>> configure_interp(section__cmap_K="plasma", section__wt_color="white")
   >>> PYCSAMT_INTERP.default.section.cmap_K
   'plasma'
   >>> reset_interp()

Unlike the style families above, every :mod:`pycsamt.interp.plot` class
defaults ``style=None``, which always resolves to this live singleton --
no ``style=PYCSAMT_INTERP.default`` hand-off needed. See :doc:`interpretation`
for the full preset comparison and that contrast explained.

Topography
----------

Whether and how terrain elevation is drawn on depth-like sections --
period/frequency sections are silently skipped regardless of this setting:

.. code-block:: pycon

   >>> from pycsamt.api import configure_topo, PYCSAMT_TOPO, reset_topo

   >>> configure_topo(enabled=True, exaggeration=2.0)
   >>> PYCSAMT_TOPO.enabled, PYCSAMT_TOPO.exaggeration
   (True, 2.0)
   >>> reset_topo()

:doc:`../user_guide/topo/concepts` covers this singleton in depth, including
:meth:`~pycsamt.topo.config.TopoConfig.is_active_for` -- the same
depth-vs-period gate that :doc:`section` and :doc:`station_rendering` build
on.

Mesh Display
------------

Rectilinear and triangular meshes share one preset system --
``"filled"`` (colour only), ``"review"`` (colour + cell edges), or
``"diagram"`` (edges only):

.. code-block:: pycon

   >>> from pycsamt.api import configure_mesh, PYCSAMT_MESH, reset_mesh

   >>> configure_mesh(review__edge__alpha=0.4, review__edge__linewidth=0.25)
   >>> PYCSAMT_MESH.review.edge.alpha
   0.4
   >>> reset_mesh()

The full story -- both mesh families, worked ``draw_mesh``/``draw_tri_mesh``
examples, and the context-manager form -- is on the :doc:`mesh` page.

Agents
------

Cap what AI-assisted workflows may spend, and pick the LLM provider only
when needed:

.. code-block:: pycon

   >>> from pycsamt.agents import AGENT_CONFIG

   >>> _ = AGENT_CONFIG.set_budget(usd=5.0)
   >>> AGENT_CONFIG.spent_usd, AGENT_CONFIG.remaining_usd
   (0.0, 5.0)
   >>> _ = AGENT_CONFIG.reset_budget(cap=True)

``AGENT_CONFIG`` breaks from the dotted-path convention on purpose -- its
settings are plain keyword arguments, not a nested style tree. The full
guide, including credential resolution order and ``.env.local``, lives at
:doc:`../user_guide/agents/llm_configuration`; :doc:`agent_config` is the
short version in this family's usual place.

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

   from pycsamt.api import (
       reset_api_view,
       reset_cli,
       reset_control,
       reset_interp,
       reset_mesh,
       reset_ordering,
       reset_pipe,
       reset_plot_config,
       reset_section,
       reset_station_rendering,
       reset_style,
       reset_topo,
   )
   from pycsamt.agents import reset_agents

   reset_api_view()
   reset_cli()
   reset_control()
   reset_interp()
   reset_mesh()
   reset_ordering()
   reset_pipe()
   reset_plot_config()
   reset_section()
   reset_station_rendering()
   reset_style()
   reset_topo()
   reset_agents()

Next Steps
----------

* :doc:`overview` for the shared ``PYCSAMT_<FAMILY>`` / ``configure_*`` /
  ``reset_*`` / ``use_*`` contract itself.
* :doc:`style`, :doc:`interpretation`, :doc:`agent_config`, :doc:`section`,
  :doc:`station_rendering`, :doc:`view_controls`, :doc:`views`, and
  :doc:`mesh` for each family's own worked depth.
* :doc:`../getting_started/configuration` for the settings most users touch
  on day one, without the full family-by-family tour.
