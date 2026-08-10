.. _api-contour-configuration:

Contour overlays
================

Contour lines can make boundaries in a continuous color field easier to
follow, but excessive levels or heavy lines can conceal the underlying data.
pyCSAMT therefore separates the color renderer from a shared contour style.
Plots that support this API resolve their contour defaults through
:data:`pycsamt.api.PYCSAMT_CONTOUR`, while an explicit function argument
remains authoritative.

The live default is the ``review`` style used by the Getting Started survey
fingerprint:

.. code-block:: pycon

   >>> from pycsamt.api import PYCSAMT_CONTOUR
   >>> print(PYCSAMT_CONTOUR)
   PyCSAMTContour
     enabled     = True
     levels      = 7
     colors      = '#202020'
     linewidths  = 0.8
     linestyles  = 'solid'
     alpha       = 0.8
     labels      = False

These are rendering defaults, not scientific thresholds. Seven automatically
spaced levels summarize the displayed numeric range; they do not identify
geological units or confidence classes.

Choose a preset
---------------

Four named alternatives cover common presentation contexts:

.. list-table::
   :header-rows: 1
   :widths: 22 78

   * - Preset
     - Intended use
   * - ``"review"``
     - Visible dark lines for interactive inspection; this is the package
       default.
   * - ``"subtle"``
     - Thin, translucent lines when the color field should dominate.
   * - ``"publication"``
     - Moderately weighted, high-opacity lines with numeric labels enabled.
   * - ``"off"``
     - Disable contour overlays while retaining the configured style values.

Select a preset for subsequent compatible plots:

.. code-block:: pycon

   >>> from pycsamt.api import PYCSAMT_CONTOUR, use_contour
   >>> use_contour("subtle")
   >>> PYCSAMT_CONTOUR.default.linewidths
   0.35
   >>> PYCSAMT_CONTOUR.default.alpha
   0.45

Use :func:`pycsamt.api.reset_contour` to restore ``review`` as the live
default.

Configure the live style
------------------------

:func:`pycsamt.api.configure_contour` accepts direct attributes for the live
default and dotted paths for named presets:

.. code-block:: pycon

   >>> from pycsamt.api import configure_contour, reset_contour
   >>> configure_contour(
   ...     levels=9,
   ...     colors="white",
   ...     linewidths=1.0,
   ...     linestyles="solid",
   ...     alpha=0.9,
   ... )
   >>> PYCSAMT_CONTOUR.default.levels
   9
   >>> configure_contour(publication__label_fmt="%.1f")
   >>> reset_contour()

``levels`` may be an integer or an explicit tuple of numeric boundaries. Use
an integer for exploratory display. Use explicit values only when the quantity
and thresholds have a defensible physical meaning, and state the units in the
caption or colorbar.

Temporary overrides
-------------------

Use a context when one figure needs a different style:

.. code-block:: pycon

   >>> with PYCSAMT_CONTOUR.context(
   ...     "publication",
   ...     levels=(0.2, 0.4, 0.6, 0.8),
   ...     colors="#111111",
   ... ):
   ...     fig = plot_function(...)

The previous configuration is restored when the block exits, including when
plotting raises an exception.

Per-call control
----------------

:func:`pycsamt.emtools.plot_survey_fingerprint` demonstrates the precedence
used by compatible plots:

1. ``PYCSAMT_CONTOUR.default`` supplies the base style;
2. ``contours=True`` or ``False`` overrides only the enabled state;
3. ``contour_kws`` overrides Matplotlib contour keywords for that call.

.. code-block:: pycon

   >>> fig = plot_survey_fingerprint(
   ...     sites,
   ...     render="imshow",
   ...     contours=True,
   ...     contour_kws={
   ...         "levels": (-20, -10, 0, 10, 20),
   ...         "colors": "white",
   ...         "linewidths": 1.0,
   ...         "linestyles": "dashed",
   ...         "alpha": 0.9,
   ...     },
   ... )

Passing ``contours=False`` suppresses the overlay without changing the global
configuration. Passing ``None``—the function default—uses the configured
enabled state.

Align stacked panels by station
-------------------------------

Contours follow values and may cross several station columns. When multiple
fingerprint panels share the same station axis, ``station_grid=True`` draws a
vertical guide through every station centre on every panel:

.. code-block:: pycon

   >>> fig = plot_survey_fingerprint(
   ...     sites,
   ...     render="imshow",
   ...     station_grid=True,
   ... )

The guides are disabled by default because dense surveys may become visually
busy. Their built-in style is a light dotted line with moderate transparency.
Override Matplotlib line properties for one figure with
``station_grid_kws``:

.. code-block:: pycon

   >>> fig = plot_survey_fingerprint(
   ...     sites,
   ...     station_grid=True,
   ...     station_grid_kws={
   ...         "color": "white",
   ...         "linewidth": 0.9,
   ...         "linestyle": "--",
   ...         "alpha": 0.8,
   ...         "zorder": 3,
   ...     },
   ... )

The controls are those accepted by
:meth:`matplotlib.axes.Axes.axvline`. Guides are placed at station centres
(``0.5, 1.5, ...``), identically for ``imshow`` and ``pcolormesh``. They show
sampling alignment only; they do not imply station spacing, interpolation
support, or a geological boundary.

Labels and scientific interpretation
------------------------------------

The ``publication`` preset enables labels through ``labels=True``. Label
format, font size, and inline placement are controlled by ``label_fmt``,
``label_fontsize``, and ``label_inline``. Labels should remain sparse enough
to preserve the color field and station annotations.

Contours connect equal values in the plotting grid. With an ``imshow`` view,
bilinear interpolation may smooth the color display, while contour vertices
still derive from the gridded values supplied to Matplotlib. Neither operation
adds measurements between stations or periods. Use ``pcolormesh`` and explicit
levels when exact cell boundaries and thresholds matter more than visual
continuity.

Reset configuration
-------------------

.. code-block:: pycon

   >>> from pycsamt.api import reset_contour
   >>> reset_contour()

Resetting affects future plots only. It does not modify figures already
created or files already saved.
