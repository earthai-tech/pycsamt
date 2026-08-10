.. _getting-started-configuration:

Configure a first session
=========================

pyCSAMT works with its defaults, so configuration is not a prerequisite for
loading a first survey. For a reproducible project, however, make three choices
explicit near the beginning of the script or notebook:

* how stations are ordered;
* where generated figures are written;
* which visual style is used.

These settings are process-local. They affect subsequent calls in the current
Python process but do not rewrite the field data. The complete configuration
system—including views, pipelines, meshes, sections, interpretation, CLI, and
agents—is documented in :doc:`../api_guide/configuration`.

Start with a small explicit setup
---------------------------------

The following setup preserves the loader's station order and gives figures a
predictable destination and appearance:

.. code-block:: pycon

   >>> from pycsamt.api import (
   ...     PLOT_CONFIG,
   ...     PYCSAMT_ORDERING,
   ...     configure_ordering,
   ...     set_dpi,
   ...     set_savedir,
   ...     use_style,
   ... )
   >>> configure_ordering(mode="input")
   SiteOrderingConfig(mode='input', min_linearity=0.95, max_cross_track_ratio=0.15, min_coordinate_fraction=0.6)
   >>> use_style("pycsamt")
   >>> set_savedir("results/figures")
   >>> set_dpi(200)
   >>> PYCSAMT_ORDERING.mode
   'input'
   >>> PLOT_CONFIG.savedir, PLOT_CONFIG.dpi
   ('results/figures', 200)

``mode="input"`` is a conservative first-survey choice because the returned
row order follows file discovery. It is reproducible only when the input file
list itself is reproducible. When filenames do not reflect profile position,
use the ordering modes described below after validating the coordinates.

``set_savedir`` establishes a default for plotting functions that honor the
global figure configuration. It does not create a scientific project layout or
force every third-party Matplotlib call into that directory. Pass an explicit
output path when a function documents one.

Choose station ordering deliberately
------------------------------------

Station order controls how survey rows and profile plots correspond to
physical locations. The main choices are:

.. list-table::
   :header-rows: 1
   :widths: 22 78

   * - Mode
     - Use it when
   * - ``"input"``
     - The supplied file or object sequence already has the required order.
   * - ``"station"``
     - Station identifiers encode the desired numeric order, such as ``S2``
       before ``S10``.
   * - ``"latitude"`` or ``"longitude"``
     - A cardinal coordinate direction represents the profile adequately.
   * - ``"chainage"``
     - Valid coordinates define a survey line and order should follow distance
       projected along it.
   * - ``"auto"``
     - pyCSAMT may use coordinate-derived chainage when the geometry passes its
       single-line checks, otherwise preserving input order.

Automatic ordering is convenient, but it cannot detect every coordinate or
survey-layout error. Compare the resulting station sequence with the field
manifest before using row position as profile distance. See
:func:`pycsamt.emtools.ensure_sites` and :doc:`../user_guide/data_loading` for
the processing boundary where ordering is applied.

Choose a plotting preset
------------------------

Named styles provide consistent colors and line conventions without requiring
new users to configure individual components:

.. code-block:: pycon

   >>> from pycsamt.api import use_style
   >>> use_style("publication")

Use ``"pycsamt"`` for the normal project appearance, ``"publication"`` for
print-oriented figures, and ``"dark"`` only when the surrounding medium also
uses a dark background. Plot style changes presentation, not the underlying
resistivity, phase, uncertainty, or inversion result.

For multi-format output, set the formats and resolution together:

.. code-block:: pycon

   >>> from pycsamt.api import set_dpi, set_fmt, set_savedir
   >>> set_fmt("png", "pdf")
   >>> set_dpi(300)
   >>> set_savedir("results/publication")

PNG is convenient for web pages and notebooks; PDF preserves vector content
for publication when the plotting elements support it. Higher DPI improves
raster resolution but does not add information absent from the data.

Keep API result views at their default initially
------------------------------------------------

Public table-producing functions normally return pyCSAMT view objects such as
``APIFrame``. These retain metadata and provide readable summaries while still
allowing conversion to pandas:

.. code-block:: pycon

   >>> from pycsamt.api import read_edis
   >>> survey = read_edis(
   ...     "data/AMT/WILLY_DATA/L18PLT",
   ...     recursive=False,
   ...     progress=False,
   ... )
   >>> summary = survey.summary()
   >>> type(summary).__name__
   'APIFrame'
   >>> summary.to_pandas(copy=True).shape
   (28, 6)

There is usually no reason to change this setting during Getting Started. If
an integration requires raw pandas outputs globally, use
:func:`pycsamt.api.configure_api_view` as described in
:doc:`../api_guide/views`.

Use temporary settings for isolated work
----------------------------------------

Configuration singletons provide context managers when one block needs a
temporary override. The previous values are restored even if the block raises
an exception:

.. code-block:: pycon

   >>> from pycsamt.api import PYCSAMT_ORDERING
   >>> before = PYCSAMT_ORDERING.mode
   >>> with PYCSAMT_ORDERING.context(mode="station"):
   ...     print(PYCSAMT_ORDERING.mode)
   station
   >>> PYCSAMT_ORDERING.mode == before
   True

This is safer in notebooks and reusable libraries than changing a global
setting for one operation and relying on a later cell to restore it.

Reset the settings you changed
------------------------------

Reset helpers restore package defaults for the current process:

.. code-block:: pycon

   >>> from pycsamt.api import (
   ...     reset_ordering,
   ...     reset_plot_config,
   ...     reset_style,
   ... )
   >>> reset_ordering()
   >>> reset_plot_config()
   >>> reset_style()

Resetting does not remove output files and does not undo transformations
already applied to data objects.

Configuration that can wait
---------------------------

Do not configure every subsystem before the first survey. Add settings when a
workflow actually needs them:

* pipeline output and failure policies: :doc:`../user_guide/pipeline/index`;
* mesh, section, station, and interpretation rendering:
  :doc:`../api_guide/configuration`;
* AI framework selection: :doc:`../user_guide/ai_inversion/index`;
* agent providers, credentials, and spending limits:
  :doc:`../user_guide/agents/llm_configuration`;
* command-line defaults: :doc:`../cli/index`.

.. important::

   Never place provider API keys in documentation, notebooks committed to
   version control, pipeline configuration files, or command history. Use the
   provider's supported environment variable or a local secret store.

With input format and session defaults established, continue to
:doc:`first_survey` to load and inspect the first survey line.
