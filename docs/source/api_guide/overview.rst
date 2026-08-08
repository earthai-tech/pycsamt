.. _api-guide-overview:

API configuration
=================

:mod:`pycsamt.api` is the front door of the package.  It is a single import
point that gathers three things every workflow needs:

1. **Readers** that turn field data into survey objects —
   :func:`~pycsamt.api.read_edis`, :func:`~pycsamt.api.read_edi`,
   :func:`~pycsamt.api.read_sites`.
2. **Views** that present tabular results for interactive work — the
   opt-in :class:`~pycsamt.api.APIFrame` / :class:`~pycsamt.api.APIResult`
   layer behind ``api=True``.
3. **Configuration** for every runtime behaviour of the package — output
   directories, plot styles, axis conventions, CLI defaults, agent budgets —
   through one consistent ``configure_*`` / ``reset_*`` pattern.

Everything documented here is importable directly from ``pycsamt.api``:

.. code-block:: pycon

   >>> from pycsamt.api import (
   ...     read_edis,            # data in
   ...     configure_api_view,   # what api=True returns
   ...     configure_pipe,       # where pipeline outputs go
   ...     configure_ordering,   # how survey stations are ordered
   ...     use_style,            # how figures look
   ... )

   >>> configure_api_view(backend="pycsamt")
   >>> _ = configure_pipe(output_root="results/run01", plot_dpi=200)
   >>> _ = configure_ordering(mode="auto")
   >>> use_style("publication")

   >>> survey = read_edis("data/AMT/WILLY_DATA/L18PLT", progress="auto")
   >>> print(survey.summary())        # APIFrame: compact, metadata-rich
   APIFrame: edi_survey_summary
   kind: edi.summary
   shape: 28 rows x 6 columns
   columns: station, path, n_freq, tipper, spectra, ts
   numeric: 1 columns
   missing: 0.0%
   source: data/AMT/WILLY_DATA/L18PLT

One Pattern Everywhere
----------------------

Each behaviour family follows the same contract, so learning one family
teaches you all of them:

.. list-table::
   :header-rows: 1
   :widths: 26 74

   * - Entry point
     - Role
   * - ``PYCSAMT_<FAMILY>``
     - Module-level singleton holding the current settings.  Print it to see
       the active configuration.
   * - ``configure_<family>(**kw)``
     - Update settings with dotted-path keywords, e.g.
       ``configure_style(mt__xy__color="#003f88")``.
   * - ``reset_<family>()``
     - Restore package defaults for the current session.
   * - ``use_<family>(preset)``
     - Where presets exist (styles, interpretation), apply a named preset in
       one call, e.g. ``use_style("publication")``.

Configuration Families
----------------------

.. list-table::
   :header-rows: 1
   :widths: 22 38 40

   * - Family
     - Entry points
     - Controls
   * - View layer
     - :func:`~pycsamt.api.configure_api_view`, ``PYCSAMT_API_VIEW``
     - What ``api=True`` returns: pyCSAMT views, pandas, or a custom wrapper.
   * - Pipeline
     - :func:`~pycsamt.api.configure_pipe`, ``PYCSAMT_PIPE``
     - Output roots, plot DPI/format, progress display, step-error policy.
   * - Site ordering
     - :func:`~pycsamt.api.configure_ordering`, ``PYCSAMT_ORDERING``
     - Package-wide station order: validated coordinate chainage, input,
       natural station name, latitude, or longitude.
   * - CLI
     - :func:`~pycsamt.api.configure_cli`, ``PYCSAMT_CLI``
     - Logging level, output format and directory, parallel build jobs.
   * - Plot styles
     - :func:`~pycsamt.api.configure_style`, :func:`~pycsamt.api.use_style`,
       ``PYCSAMT_STYLE`` -- :doc:`style`
     - MT component colours, multiline gradients, correction, rose, and
       phase-tensor ellipse styles.  Presets: ``pycsamt``, ``publication``,
       ``dark``, ``modem``.
   * - Figure output
     - :func:`~pycsamt.api.save_fig`, :func:`~pycsamt.api.set_dpi`,
       :func:`~pycsamt.api.set_fmt`, :func:`~pycsamt.api.set_savedir`,
       ``PLOT_CONFIG``
     - Global figure saving defaults: DPI, formats, output directory.
   * - View controls
     - :func:`~pycsamt.api.configure_control`,
       :func:`~pycsamt.api.wrap_phase`, ``PYCSAMT_CONTROL`` --
       :doc:`view_controls`
     - Apparent-resistivity and phase axis behaviour, frequency-axis
       direction, phase wrapping.
   * - Sections
     - :func:`~pycsamt.api.configure_section`, ``PYCSAMT_SECTION`` --
       :doc:`section`
     - Resistivity-section figure, axis, and colourbar styling.
   * - Station rendering
     - :func:`~pycsamt.api.configure_station_rendering`,
       ``PYCSAMT_STATION_RENDERING`` -- :doc:`station_rendering`
     - Station map markers and axis styling.
   * - Interpretation
     - :func:`~pycsamt.api.configure_interp`,
       :func:`~pycsamt.api.use_interp`, ``PYCSAMT_INTERP`` --
       :doc:`interpretation`
     - Hydrogeological profile and section styles.
   * - Topography
     - :func:`~pycsamt.api.configure_topo`, ``PYCSAMT_TOPO`` --
       :doc:`../user_guide/topo/concepts`
     - Topography handling and depth/frequency y-axis conventions.
   * - Mesh display
     - :func:`~pycsamt.api.configure_mesh`, :func:`~pycsamt.api.draw_mesh`,
       :func:`~pycsamt.api.draw_tri_mesh`, ``PYCSAMT_MESH`` -- :doc:`mesh`
     - Rectilinear and triangular mesh rendering: filled, reviewed (fill +
       edges), or diagram (edges only) presets, shared by both mesh
       families.
   * - Agents
     - :func:`~pycsamt.api.configure_agents`, ``AGENT_CONFIG`` --
       :doc:`agent_config`
     - LLM provider selection and spending budget for AI agents.

Besides the families above, ``pycsamt.api`` also exposes shared axis-label
constants (``STATION_LABEL``, ``FREQUENCY_LABEL``, ``PERIOD_LABEL``, ...) and
the base objects (:class:`~pycsamt.api.PyCSAMTObject`,
:class:`~pycsamt.api.MetadataMixin`) that carry metadata through workflows.
To build your own classes on these bases — and inherit the MT math and
EDI interop that come with them — see :doc:`/development/extending`.

In This Section
---------------

:doc:`views`
    The opt-in view layer: ``api=True``, :class:`~pycsamt.api.APIFrame`,
    multi-table :class:`~pycsamt.api.APIResult` containers, global backend
    policy, custom wrappers, and the public readers.

:doc:`mesh`
    Rendering computational and inversion meshes: the shared ``"filled"``/
    ``"review"``/``"diagram"`` preset system, rectilinear vs. triangular
    mesh drawing, and dotted-path style configuration.

:doc:`style`
    Plot styles in depth: MT component colours, multiline gradients,
    correction pairs, raw-data style, phase-tensor ellipses, rose
    diagrams, and the four named presets.

:doc:`interpretation`
    Hydrogeophysical plot styles: section and profile colours, the four
    presets, and why these plot classes track live configuration changes
    without an explicit style hand-off.

:doc:`agent_config`
    ``AGENT_CONFIG`` in the family lineup, with a fast quick-reference;
    points to :doc:`../user_guide/agents/llm_configuration` for the full
    guide.

:doc:`section`
    Section-plot figure sizing, axis direction, and the topography-
    awareness gate shared by every pseudosection and inversion section.

:doc:`station_rendering`
    Station tick marks, the adaptive label-thinning algorithm, and
    terrain-following marker placement.

:doc:`view_controls`
    Apparent-resistivity scale, phase wrapping, and the frequency/period
    axis convention read by most :mod:`pycsamt.emtools` plotting
    functions.

:doc:`configuration`
    The configuration system in depth: the dotted-path convention, worked
    examples for each family, recommended setups, environment variables,
    and how to reset everything.

.. seealso::

   :doc:`../getting_started/configuration`
       First-session setup: the settings most users touch on day one.

   :doc:`../api/api`
       Auto-generated reference for every ``pycsamt.api`` function and class.
