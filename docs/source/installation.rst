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

* **Python 3.9 or later** (3.9–3.12 are tested in CI).
* Linux, macOS, and Windows are supported for the Python package. External
  solver toolchains have additional platform requirements described below.
* The core installation pulls in a deliberately small scientific stack:

.. list-table::
   :header-rows: 1
   :widths: 27 20 53

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
   * - Triangle
     - 20220202
     - Quality-graded triangular meshes for forward modelling
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
   * - scikit-learn
     - 1.1
     - Shared compatibility layer and machine-learning utilities

Python 3.9 installations constrain NumPy to the 1.x line for compatibility
with the supported Matplotlib and desktop stack. pip applies this marker
automatically from the package metadata.

Standard Install
----------------

.. code-block:: bash

   python -m pip install pycsamt          # core: I/O, processing, plotting, CLI
   python -m pip install "pycsamt[full]"  # broad development installation

``full`` bundles ``torch``, ``geo``, ``dev``, ``docs``, ``app``, and
``agents``. It intentionally prefers the PyTorch backend. It does **not**
include ``tensorflow`` or the separate ``agent-master`` application stack;
add either one explicitly when required.

Optional Feature Groups
-----------------------

Every group can be combined freely, for example
``python -m pip install "pycsamt[torch,geo,agents]"``.

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
     - pyproj, xarray, contextily, h5py
     - Reprojection, gridded data, web basemaps, and HDF5 elevation data
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
   * - ``full``
     - ``torch`` + ``geo`` + ``dev`` + ``docs`` + ``app`` + ``agents``
     - Broad development setup; excludes ``tensorflow`` and ``agent-master``

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

pyCSAMT itself installs with pip, but conda users can create and manage the
environment first. The repository's ``environment.yml`` creates an editable
development environment with the ``dev``, ``docs``, ``geo``, ``web``, and
``agents`` extras:

.. code-block:: bash

   conda env create -f environment.yml
   conda activate pycsamt

Run those commands from a source checkout. The package is already installed
by the environment file; do not install ``.[full]`` again unless you also
want its additional desktop and PyTorch dependencies.

A minimal manual equivalent:

.. code-block:: bash

   conda create -n pycsamt python=3.11
   conda activate pycsamt
   python -m pip install "pycsamt[full]"

Install From Source
-------------------

Clone the repository and install it in editable mode:

.. code-block:: bash

   git clone https://github.com/earthai-tech/pycsamt.git
   cd pycsamt
   git checkout v2
   python -m pip install -e ".[dev,docs]"

The ``v2`` branch is the active development branch. An editable (``-e``)
install picks up local code changes without reinstalling. Add the workflow
extras you need—for example ``geo``, ``torch``, or ``app``—and read
:doc:`contributing` before opening a pull request.

Compiled Inversion Solvers
--------------------------

Occam2D and ModEM source code is vendored under
``pycsamt/models/*/_source/``. These are external executables: installing or
importing pyCSAMT does not compile them, and they are not built with ``f2py``.
Build them explicitly from a source checkout with the dispatcher:

.. code-block:: bash

   bash pycsamt/models/_solver_build/build.sh occam2d
   bash pycsamt/models/_solver_build/build.sh modem2d
   bash pycsamt/models/_solver_build/build.sh modem3d

Use ``--help`` on the dispatcher or an individual script before enabling
options such as ``--auto-install``. Occam2D and the serial ModEM builds need
``gfortran`` and ``make``; ModEM also needs linkable LAPACK/BLAS libraries.

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Platform
     - Typical toolchain
   * - Linux
     - gfortran, make, and LAPACK/BLAS development packages from the
       distribution package manager
   * - macOS
     - Homebrew ``gcc`` (which provides gfortran), ``make``, and OpenBLAS
   * - Windows
     - The build scripts can create an isolated conda MinGW-w64 environment;
       WSL is also supported

MARE2DEM is different: its source is downloaded separately and its build
requires an Intel MPI compiler toolchain and MKL on Linux, macOS, or WSL. See
``pycsamt/models/_solver_build/README.md`` in the source checkout for the
solver matrix, build options, and platform-specific details.

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

For a backend-aware check, use the public backend registry:

.. code-block:: bash

   python -c "from pycsamt.backends import list_backends; print(list_backends())"

Application commands can also be checked with ``--help`` after installing
their matching extras.

Upgrade Or Remove
-----------------

.. code-block:: bash

   python -m pip install --upgrade pycsamt  # latest release
   python -m pip uninstall pycsamt          # leaves your data untouched

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

**A solver build fails.**
Run the relevant build script with ``--help``, then confirm that the required
compiler, ``make``, and numerical libraries are visible in the same shell.
Consult ``pycsamt/models/_solver_build/README.md`` before changing compilers
or attempting a manual build.

Next Steps
----------

* :doc:`getting_started/data_formats` — identify your field data format.
* :doc:`getting_started/configuration` — configure outputs and styles.
* :doc:`getting_started/first_survey` — load and QC a first survey.
