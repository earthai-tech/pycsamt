.. _installation:

Installation
============

pyCSAMT ships as a small core with optional feature groups, so you install
only what your workflow needs.  This page is the complete installation
reference: requirements, every optional extra, console commands, compiled
solvers, and verification.  For a guided, step-by-step environment setup,
see :doc:`getting_started/installation` instead.

Requirements
------------

* **Python 3.9 or later** (3.9 – 3.13 are tested in CI).
* Any platform — Linux, macOS, and Windows are supported.
* The core installation pulls in a deliberately small scientific stack:

.. list-table::
   :header-rows: 1
   :widths: 30 20 50

   * - Package
     - Minimum
     - Used for
   * - NumPy
     - 1.22
     - Array mathematics throughout the package
   * - SciPy
     - 1.8
     - Signal processing, interpolation, optimisation
   * - Matplotlib
     - 3.5
     - All plotting
   * - pandas
     - 1.4
     - Tabular results and the API view layer
   * - PyYAML
     - 5.4
     - Pipeline and configuration files
   * - tqdm
     - 4.60
     - Progress bars
   * - click / rich
     - 8.1 / 13.0
     - The ``pycsamt`` command-line interface

Standard Install
----------------

.. code-block:: bash

   pip install pycsamt              # core: I/O, processing, plotting, CLI
   pip install "pycsamt[full]"      # everything below in one command

``full`` bundles ``torch``, ``geo``, ``dev``, ``docs``, ``app``, and
``agents``.  It intentionally prefers the PyTorch backend — add
``tensorflow`` explicitly if you need Keras models.

Optional Feature Groups
-----------------------

Every group can be combined freely, e.g.
``pip install "pycsamt[torch,geo,agents]"``.

.. list-table::
   :header-rows: 1
   :widths: 18 40 42

   * - Extra
     - Installs
     - Enables
   * - ``torch``
     - PyTorch ≥ 1.13
     - PINN and hybrid deep-learning inverters (recommended backend)
   * - ``tensorflow``
     - TensorFlow ≥ 2.10
     - The Keras/TensorFlow model backend
   * - ``geo``
     - pyproj, xarray, contextily
     - Coordinate reprojection, gridded data, web basemaps for station maps
   * - ``agents``
     - anthropic, openai, google-generativeai
     - LLM-driven agents: Claude, OpenAI (and DeepSeek via the OpenAI SDK),
       and Gemini providers
   * - ``desktop``
     - PySide6, pyqtgraph, contextily
     - The native desktop application
   * - ``web``
     - Dash, dash-bootstrap-components, Plotly, diskcache, Pillow
     - The Dash web dashboard
   * - ``app``
     - ``desktop`` + ``web``
     - Both interactive applications
   * - ``agent-master``
     - ``agents`` + Dash 4 stack
     - The Agent Master web application (chat-driven workflows)
   * - ``dev``
     - pytest, pytest-cov, ruff, pre-commit
     - Running the test suite and contributing
   * - ``docs``
     - Sphinx, PyData theme, numpydoc, MyST, sphinx-design, ...
     - Building this documentation locally

Console Commands
----------------

Installing pyCSAMT registers these entry points (application commands
require the matching extra):

.. list-table::
   :header-rows: 1
   :widths: 28 30 42

   * - Command
     - Requires
     - Launches
   * - ``pycsamt``
     - core
     - The command-line interface (``pycsamt --help``)
   * - ``pycsamt-desktop`` / ``pycsamt-gui``
     - ``desktop``
     - The native desktop application (both names are equivalent)
   * - ``pycsamt-web``
     - ``web``
     - The Dash web dashboard
   * - ``pycsamt-agent``
     - ``agent-master``
     - The Agent Master web application
   * - ``pycsamt-mapview``
     - ``web``
     - The map-view workbench

Conda Environments
------------------

pyCSAMT itself installs from PyPI, but conda users can create the
environment first and then pip-install into it.  The repository ships an
``environment.yml`` for a reproducible setup:

.. code-block:: bash

   conda env create -f environment.yml
   conda activate pycsamt
   pip install -e ".[full]"

A minimal manual equivalent:

.. code-block:: bash

   conda create -n pycsamt python=3.11
   conda activate pycsamt
   pip install "pycsamt[full]"

Install From Source
-------------------

Track the active development branch (``v2``):

.. code-block:: bash

   git clone https://github.com/earthai-tech/pycsamt.git
   cd pycsamt
   git checkout v2
   pip install -e ".[full]"

An editable (``-e``) install picks up local code changes without
reinstalling — pair it with the ``dev`` extra when contributing, and read
:doc:`contributing` before opening a pull request.

Compiled Inversion Solvers
--------------------------

The classical inversion solvers (Occam2D, ModEM) are distributed as Fortran
source under ``pycsamt/models/*/_source/`` and are compiled on first use via
``f2py``.  A Fortran compiler must be available:

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Platform
     - Compiler
   * - Linux
     - ``sudo apt install gfortran`` (Debian/Ubuntu) or the distribution
       equivalent
   * - macOS
     - ``brew install gcc`` (provides ``gfortran``)
   * - Windows
     - MinGW-w64 (e.g. via ``conda install m2w64-toolchain``) or WSL

Pure-Python workflows — processing, QC, plotting, PINN inversion, agents,
and apps — do not need a Fortran compiler.

Verify The Installation
-----------------------

.. code-block:: bash

   python -c "import pycsamt; print(pycsamt.__version__)"
   pycsamt --help

Check optional pieces only if you installed them:

.. code-block:: python

   import torch                    # [torch]
   import pyproj                   # [geo]
   import anthropic                # [agents]

A fuller verification checklist — including backends, apps, and agents — is
in :doc:`getting_started/installation`.

Upgrade Or Remove
-----------------

.. code-block:: bash

   pip install --upgrade pycsamt        # latest release
   pip uninstall pycsamt                # remove (leaves your data untouched)

Troubleshooting
---------------

**pip resolves an old version.**
Upgrade the installer first: ``python -m pip install --upgrade pip``.

**PyTorch or TensorFlow wheels fail to install.**
Install the backend on its own first, following the selector on
`pytorch.org <https://pytorch.org/get-started/locally/>`_ or
`tensorflow.org <https://www.tensorflow.org/install>`_, then install
pyCSAMT without that extra.

**Qt platform errors when launching** ``pycsamt-desktop`` **on Linux.**
Install the system libraries PySide6 needs, e.g.
``sudo apt install libxcb-cursor0``.

**Fortran solver compilation fails.**
Confirm ``gfortran --version`` works in the same shell; on Windows prefer
the conda toolchain or WSL.

More scenarios are covered in the
:doc:`getting started troubleshooting section <getting_started/installation>`.

Next Steps
----------

* :doc:`getting_started/data_formats` — identify your field data format.
* :doc:`getting_started/configuration` — configure outputs and styles.
* :doc:`getting_started/first_survey` — load and QC a first survey.

.. toctree::
   :maxdepth: 1
   :hidden:

   getting_started/index
