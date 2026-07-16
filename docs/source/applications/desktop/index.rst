.. _applications-desktop:

Desktop GUI
===========

The pyCSAMT desktop GUI is the local application for interactive CSAMT, AMT,
and MT survey work. It is designed for the day-to-day workflow of loading EDI
data, checking station geometry, inspecting profiles and maps, running QC,
testing corrections, preparing forward and inversion models, exporting figures
and processed EDIs, and documenting the processing chain.

The desktop app uses the same scientific package as the Python API. The GUI
is not a separate science layer; it is an interactive surface over pyCSAMT
loaders, ``Sites`` objects, processing functions, plotting tools, inversion
builders, and export helpers.

.. figure:: ../../_static/applications/desktop/desktop-walkthrough.gif
   :alt: Animated tour of the pyCSAMT desktop GUI
   :class: pycsamt-screenshot

   An animated tour of the desktop app: launch the main window, load EDI survey
   lines, inspect the station map and per-station responses, run quality
   control, correct static shift, and build a forward model. Each step is
   documented on the pages below.

A good first pass, in order, is **Overview → Installation → Loading &
Sessions → Workspace → Maps & Profiles → Processing Workflows → Exports →
Troubleshooting**.

.. toctree::
   :maxdepth: 3
   :class: pycsamt-guide-toc

   overview
   installation
   loading_and_sessions
   workspace
   maps_and_profiles
   processing_workflows
   exports
   troubleshooting
