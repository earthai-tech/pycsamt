.. _applications-desktop-installation:

Installation And Launch
=======================

The desktop application is the local Qt interface for pyCSAMT.  It is meant
for interactive work: loading EDI surveys, inspecting station metadata,
opening map and profile views, running quality-control and correction tools,
and preparing forward or inversion workflows from one workspace.

Install The Desktop Extra
-------------------------

From a released package, install the desktop extra:

.. code-block:: bash

   pip install "pycsamt[desktop]"

If you also want the Dash web app in the same environment, install the
application extra instead:

.. code-block:: bash

   pip install "pycsamt[app]"

For development from a source checkout, install the package in editable mode
with the desktop and documentation dependencies:

.. code-block:: bash

   git clone https://github.com/earthai-tech/pycsamt.git
   cd pycsamt
   pip install -e ".[desktop,docs]"

The desktop extra installs the Qt runtime and plotting dependencies used by
the GUI:

.. list-table::
   :header-rows: 1
   :widths: 28 72

   * - Dependency
     - Used for
   * - ``PySide6``
     - the desktop windowing toolkit, dialogs, menus, and toolbar controls
   * - ``pyqtgraph``
     - responsive station and signal views
   * - ``contextily``
     - optional map tile support for geospatial views

Launch The Application
----------------------

After installation, start the desktop GUI from the active environment:

.. code-block:: bash

   pycsamt-desktop

The older command remains available as an alias:

.. code-block:: bash

   pycsamt-gui

From a source checkout, the module entry point is equivalent:

.. code-block:: bash

   python -m pycsamt.app.desktop

On first launch the application opens with an empty survey workspace.  The
left panel is reserved for stations, the center panel summarizes the survey,
and the toolbar exposes the main workflow windows.

.. figure:: ../../_static/applications/desktop/home.png
   :alt: Empty pyCSAMT desktop workspace
   :class: pycsamt-screenshot

   Empty desktop workspace before survey data are loaded.

Load A First EDI Survey
-----------------------

Use **Open / Load EDI...** to load one or more EDI files.  The loader accepts
individual files or a folder containing a survey line.  For the sample WILLY
data, a typical folder is one of the ``L18PLT``, ``L22PLT``, ``L26PLT``,
``L30PLT``, or ``L34PLT`` line folders under ``data/AMT/WILLY_DATA``.

.. figure:: ../../_static/applications/desktop/load-edis.png
   :alt: Selecting EDI files in the desktop loader
   :class: pycsamt-screenshot

   Loading a line by selecting multiple ``.edi`` files from the survey folder.

After loading, the station table, overview panel, profile view, map view, QC
tools, correction tools, forward modelling, inversion, interpretation, and
agent windows operate on the active session.

Recommended Environment
-----------------------

Use a dedicated environment for the desktop app so Qt and geospatial packages
do not conflict with other Python projects:

.. code-block:: bash

   conda create -n pycsamt-v2 python=3.12
   conda activate pycsamt-v2
   pip install -e ".[desktop,docs]"

For day-to-day work in this repository, activate the environment before
launching or building documentation:

.. code-block:: bash

   conda activate pycsamt-v2
   pycsamt-desktop

Platform Notes
--------------

Windows users should run the command from Anaconda Prompt, PowerShell, or Git
Bash after activating the environment.  If the command is not found, verify
that the environment's ``Scripts`` directory is on ``PATH``.

Linux users may need system Qt/OpenGL libraries in minimal containers or
headless servers.  The desktop app is intended for an interactive display; use
the Python API, CLI, or web app for non-interactive runs.

macOS users should launch from a terminal attached to the active environment.
If a gatekeeper prompt appears for local Python or Qt components, approve the
environment that owns the installed dependencies.

Quick Verification
------------------

Check that the command is available:

.. code-block:: bash

   pycsamt-desktop --help
   pycsamt-gui --help

Then launch the GUI:

.. code-block:: bash

   pycsamt-desktop

If the application reports that ``PySide6`` is missing, reinstall with the
desktop extra:

.. code-block:: bash

   pip install -e ".[desktop]"
