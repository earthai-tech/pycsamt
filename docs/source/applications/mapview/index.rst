.. _applications-mapview:

MapView
=======

MapView is the dedicated map workbench of pyCSAMT: a Dash/Plotly browser
application focused entirely on **seeing a survey in space** — basemap station
maps, profile overlays, Surfer-style contours, and interactive 3-D resistivity
scenes with fence sections, blocks, depth slices, iso-surfaces, and draped
topography. It is the point-and-click layer over the :doc:`Python mapping API
</user_guide/map/index>`: everything the app draws comes from the same
:class:`~pycsamt.map.MapView` façade you can drive from code, so a scene you
build here can be reproduced exactly in a script.

.. figure:: ../../_static/applications/mapview/mapview-walkthrough.gif
   :alt: Animated tour of the pyCSAMT MapView workbench
   :class: pycsamt-screenshot

   An animated tour of MapView: from the welcome screen to a survey on a
   basemap, Surfer-style contours, interactive 3-D fence and block scenes,
   isolating conductive structure, and an imported ModEM inversion model.
   Each control is documented on the pages below.

A good first pass, in order, is **Overview → Installation → Loading &
Sessions → Views & Controls → Exports → Troubleshooting**.

.. toctree::
   :numbered: 4
   :maxdepth: 1
   :class: pycsamt-guide-toc

   overview
   installation
   loading_and_sessions
   views
   exports
   troubleshooting
