.. _applications:

Applications
============

pyCSAMT v2 includes application interfaces for users who prefer an
interactive environment over direct Python or command-line workflows.
These applications sit on top of the same scientific packages documented
elsewhere in the project: data loading, site management, quality control,
forward modelling, inversion, interpretation, pipelines, and AI-assisted
agents.

.. note::

   This section is a documentation stub. The desktop and web applications
   are already present in :mod:`pycsamt.app`, but their user manuals will be
   expanded after the core workflow documentation is stabilized.

Application Surfaces
--------------------

.. list-table::
   :header-rows: 1
   :widths: 22 32 46

   * - Interface
     - Entry point
     - Purpose
   * - Desktop application
     - ``pycsamt-gui`` or ``python -m pycsamt.app.desktop``
     - A PySide6 desktop workspace for survey loading, QC, maps, profiles,
       corrections, forward modelling, inversion, interpretation, pipelines,
       and agent-assisted tasks.
   * - Web application
     - ``python -m pycsamt.app.web``
     - A Dash/Plotly browser interface that mirrors the main desktop
       workflow areas for local or network-accessible interactive use.

Install Application Dependencies
--------------------------------

The base pyCSAMT installation does not force GUI and web dependencies on
all users. Install the application extra when you want to run either
interface:

.. code-block:: bash
   :linenos:

   pip install "pycsamt[app]"

For local development from the source tree:

.. code-block:: bash
   :linenos:

   pip install -e ".[app,dev]"

Relationship To Other Documentation
-----------------------------------

The applications should be read as interface layers, not separate
scientific implementations. When documenting an application feature, link
back to the underlying workflow documentation:

.. list-table::
   :header-rows: 1
   :widths: 26 74

   * - Application area
     - Canonical documentation to reference
   * - Data loading and survey context
     - :doc:`../getting_started/data_formats`,
       :doc:`../getting_started/first_survey`, :doc:`../cli/survey`
   * - Quality control and processing
     - :doc:`../user_guide/processing`, :doc:`../pipeline/steps`
   * - Forward modelling
     - :doc:`../cli/forward`, :doc:`../api/forward`
   * - Inversion
     - :doc:`../user_guide/inversion`, :doc:`../cli/invert`,
       :doc:`../api/inversion`
   * - Pipelines
     - :doc:`../pipeline/index`, :doc:`../pipeline/cli_pipe`
   * - AI-assisted workflows
     - :doc:`../agents/index`, :doc:`../agents/agent_catalogue`

Planned Pages
-------------

The application documentation will eventually include:

* installation and platform notes for GUI dependencies;
* opening surveys and managing sessions;
* map, profile, station, and QC panels;
* correction, forward, inversion, interpretation, and pipeline windows;
* agent panel usage and background-worker behavior;
* export workflows and reproducibility notes;
* deployment notes for the web application.

.. toctree::
   :maxdepth: 1

   desktop
   web
