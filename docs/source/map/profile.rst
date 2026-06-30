Profile Maps And Pseudosections
===============================

Profile maps show apparent resistivity or phase as a function of
station/profile position and period.  They are quick-look data
sections, not inversion sections.

Function API
------------

.. code-block:: python

   from pycsamt.map import ProfileMapOptions, plot_pseudosection

   fig = plot_pseudosection(
       "data/AMT/WILLY_DATA/L18PLT",
       options=ProfileMapOptions(
           quantity="rho",
           components=("xy", "yx"),
           period_range=(0.001, 10.0),
           x_axis="distance",
       ),
   )

Builder API
-----------

.. code-block:: python

   from pycsamt.map import ProfileMap

   fig = (
       ProfileMap("data/AMT/WILLY_DATA/L18PLT")
       .with_quantity("phase")
       .with_component("xy")
       .pseudosection()
   )

Quantities And Components
-------------------------

``quantity`` accepts ``"rho"`` / ``"resistivity"`` and
``"phase"`` / ``"phi"``.  Components use the shared component parser:
``xx``, ``xy``, ``yx``, ``yy``, ``det``, and ``avg``.

Axes
----

``x_axis="station"``
   Use station IDs along the horizontal axis.

``x_axis="distance"``
   Use approximate station distance in kilometres when coordinates are
   available.

Backends
--------

The Plotly backend creates interactive heatmaps.  The Matplotlib
backend creates static report figures.
