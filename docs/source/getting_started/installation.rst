.. _getting-started-installation:

Installation
============

This page gives the recommended installation paths for pyCSAMT v2.  The short
version is:

.. code-block:: console
   :linenos:

   python -m pip install pycsamt
   pycsamt --help

For active development from the repository, use:

.. code-block:: console
   :linenos:

   git clone https://github.com/earthai-tech/pycsamt.git
   cd pycsamt
   python -m pip install -e ".[dev,docs]"
   pycsamt --help

The rest of this page explains which environment to choose, which optional
extras to install, and how to verify that the command line, Python API,
plotting stack, AI backends, desktop app, and documentation tools are ready.

Supported Python Versions
-------------------------

The package metadata requires Python ``>=3.9``.  For most scientific
installations, Python ``3.10``, ``3.11``, or ``3.12`` is the most comfortable
choice because the geospatial, GUI, and machine-learning dependencies usually
publish wheels for those versions first.

Use Python ``3.13`` only if every optional dependency you need already
supports it in your environment.  In particular, TensorFlow, PySide, and
geospatial packages can lag behind the newest Python release.

Check your interpreter before installing:

.. code-block:: console
   :linenos:

   python --version
   python -m pip --version

Recommended Environment
-----------------------

Use an isolated environment for every pyCSAMT project.  pyCSAMT sits on a
large scientific stack, and isolation prevents one project from accidentally
changing NumPy, SciPy, Matplotlib, PySide, TensorFlow, or PyTorch versions for
another project.

With ``venv``:

.. code-block:: console
   :linenos:

   python -m venv .venv
   source .venv/bin/activate
   python -m pip install --upgrade pip setuptools wheel

On Windows PowerShell:

.. code-block:: powershell
   :linenos:

   py -3.11 -m venv .venv
   .\.venv\Scripts\Activate.ps1
   python -m pip install --upgrade pip setuptools wheel

With conda or mamba:

.. code-block:: console
   :linenos:

   mamba create -n pycsamt python=3.11
   mamba activate pycsamt
   python -m pip install --upgrade pip setuptools wheel

Conda is often convenient when you need geospatial or GUI packages because it
can provide compiled system libraries.  ``venv`` is usually enough for core
data loading, tables, plotting, command-line use, and documentation work.

Install The Stable Package
--------------------------

Install the core package from PyPI:

.. code-block:: console
   :linenos:

   python -m pip install pycsamt

The core install includes the packages needed for common command-line and
Python workflows:

* ``numpy`` for arrays;
* ``scipy`` for numerical algorithms;
* ``matplotlib`` for figures;
* ``pandas`` for tables;
* ``pyyaml`` for configuration files;
* ``tqdm`` for progress bars;
* ``click`` and ``rich`` for the CLI.

After installation, verify that the package and command-line entry point are
available:

.. code-block:: console
   :linenos:

   python -c "import pycsamt; print(pycsamt.__version__)"
   pycsamt --help

Install From Source
-------------------

Use a source install when you are developing pyCSAMT, testing unreleased
features, editing documentation, or working from a branch.

.. code-block:: console
   :linenos:

   git clone https://github.com/earthai-tech/pycsamt.git
   cd pycsamt
   python -m pip install -e .

An editable install means Python imports the package directly from the cloned
repository.  Code changes are visible immediately without reinstalling the
package.

For normal development, install the development and documentation extras:

.. code-block:: console
   :linenos:

   python -m pip install -e ".[dev,docs]"

This adds ``pytest``, ``pytest-cov``, ``ruff``, ``pre-commit``, Sphinx, the
PyData Sphinx theme, and the documentation extensions used by this docs tree.

Install Optional Feature Groups
-------------------------------

pyCSAMT v2 keeps optional features out of the core install.  Choose only the
extras you need.

.. list-table::
   :header-rows: 1
   :widths: 18 34 48

   * - Extra
     - Install command
     - Use it when
   * - Core
     - ``python -m pip install pycsamt``
     - You need EDI reading, tables, plotting, CLI commands, processing, and
       pipeline basics.
   * - ``torch``
     - ``python -m pip install "pycsamt[torch]"``
     - You want PyTorch-backed AI inversion or neural processing.
   * - ``tensorflow``
     - ``python -m pip install "pycsamt[tensorflow]"``
     - You want TensorFlow/Keras-backed AI inversion or neural processing.
   * - ``geo``
     - ``python -m pip install "pycsamt[geo]"``
     - You need projection, map, xarray, or contextily-based geospatial
       workflows.
   * - ``app``
     - ``python -m pip install "pycsamt[app]"``
     - You want the desktop GUI, pyqtgraph views, or Dash web app support.
   * - ``docs``
     - ``python -m pip install "pycsamt[docs]"``
     - You want to build the documentation locally.
   * - ``dev``
     - ``python -m pip install "pycsamt[dev]"``
     - You want test, lint, and contribution tooling.
   * - ``full``
     - ``python -m pip install "pycsamt[full]"``
     - You want every declared optional group in one environment.

When installing from a cloned repository, keep the editable flag:

.. code-block:: console
   :linenos:

   python -m pip install -e ".[geo,app,dev,docs]"

.. note::

   Quote extras on shells such as ``zsh``.  For example, use
   ``"pycsamt[geo]"`` rather than ``pycsamt[geo]``.

AI And Deep-Learning Backends
-----------------------------

pyCSAMT's AI modules can use PyTorch or TensorFlow.  Install one backend
first, then verify that pyCSAMT can detect it.

PyTorch:

.. code-block:: console
   :linenos:

   python -m pip install "pycsamt[torch]"
   python -c "from pycsamt.backends import list_backends; print(list_backends())"

TensorFlow:

.. code-block:: console
   :linenos:

   python -m pip install "pycsamt[tensorflow]"
   python -c "from pycsamt.backends import list_backends; print(list_backends())"

Select the backend explicitly in a script:

.. code-block:: python
   :linenos:

   from pycsamt.backends import get_backend, set_backend

   set_backend("torch")          # or "tensorflow"
   print(get_backend())

Or set it for one command:

.. code-block:: console
   :linenos:

   PYCSAMT_AI_BACKEND=torch python my_script.py

For GPU-enabled PyTorch or TensorFlow, follow the official installation
instructions for your CUDA, ROCm, or Apple Silicon setup.  GPU packages are
platform-specific, so pyCSAMT does not try to guess the correct accelerator
wheel for you.

AI Agents And LLM Providers
---------------------------

The agent system can run many deterministic workflows with the normal
pyCSAMT installation.  LLM-backed agents require provider SDKs and API keys,
but those provider SDKs are intentionally lazy imports: install only the
provider you plan to use.

OpenAI:

.. code-block:: console
   :linenos:

   python -m pip install openai
   export OPENAI_API_KEY="..."

Anthropic:

.. code-block:: console
   :linenos:

   python -m pip install anthropic
   export ANTHROPIC_API_KEY="..."

Google Gemini:

.. code-block:: console
   :linenos:

   python -m pip install google-generativeai
   export GOOGLE_API_KEY="..."

The exact environment variables used by an agent workflow are documented in
:doc:`../agents/llm_configuration`.

Geospatial Support
------------------

Install geospatial extras when you need coordinate transforms, map tiles,
xarray-backed grids, or map-style plots:

.. code-block:: console
   :linenos:

   python -m pip install "pycsamt[geo]"

This installs ``pyproj``, ``xarray``, and ``contextily``.  Some systems may
also need compiled geospatial libraries such as PROJ or GDAL from the system
package manager or conda-forge.

With conda:

.. code-block:: console
   :linenos:

   mamba install -c conda-forge proj gdal pyproj
   python -m pip install "pycsamt[geo]"

If you see a warning about ``GDAL_DATA`` during import, non-geospatial
features can still work.  Fix the geospatial environment before relying on
map projections, raster layers, or GIS exports.

Desktop And Web App Support
---------------------------

Install the application extras if you want to launch the desktop GUI or use
the web application components:

.. code-block:: console
   :linenos:

   python -m pip install "pycsamt[app]"

Verify the desktop entry point:

.. code-block:: console
   :linenos:

   pycsamt-desktop --help

Verify the web launcher:

.. code-block:: console
   :linenos:

   pycsamt-web --no-browser --port 0

The ``app`` extra installs ``PySide6``, ``pyqtgraph``, ``dash``,
``dash-bootstrap-components``, ``diskcache``, ``multiprocess``, ``plotly``,
``Pillow``, and ``contextily``.  On Linux, Qt may also require system
libraries for X11, Wayland, OpenGL, and font rendering.  If the GUI starts
but the window is blank, first test a minimal PySide6 example in the same
environment.

Documentation Build Dependencies
--------------------------------

To build the documentation locally, install the documentation extra from the
repository:

.. code-block:: console
   :linenos:

   python -m pip install -e ".[docs]"

Then build HTML from the ``docs`` directory:

.. code-block:: console
   :linenos:

   cd docs
   sphinx-build -b html source build/html

For live documentation editing:

.. code-block:: console
   :linenos:

   cd docs
   sphinx-autobuild source build/html

For the detailed documentation workflow, warning policy, and selected-page
build commands, see :doc:`../development/documentation_build`.

Developer Installation
----------------------

For contributors, install pyCSAMT in editable mode with development tools:

.. code-block:: console
   :linenos:

   git clone https://github.com/earthai-tech/pycsamt.git
   cd pycsamt
   python -m pip install -e ".[dev,docs]"
   pre-commit install

Run the test suite:

.. code-block:: console
   :linenos:

   pytest

Run only a small package area while working:

.. code-block:: console
   :linenos:

   pytest pycsamt/api
   pytest pycsamt/agents
   pytest pycsamt/pipeline

Run linting:

.. code-block:: console
   :linenos:

   ruff check pycsamt

The public API policy and docstring conventions are documented in
:doc:`../development/api_policy` and
:doc:`../development/docstring_style`.

Verify A Complete Installation
------------------------------

After installation, run a few checks.  Start with the package import and CLI:

.. code-block:: console
   :linenos:

   python -c "import pycsamt; print(pycsamt.__version__)"
   pycsamt --help
   pycsamt config show

Check the public API facade:

.. code-block:: console
   :linenos:

   python - <<'PY'
   from pycsamt.api import configure_api_view, read_edis
   from pycsamt.backends import list_backends

   configure_api_view(backend="pycsamt")
   print("Backends:", list_backends())
   print("Reader:", read_edis)
   PY

Check the pipeline package:

.. code-block:: console
   :linenos:

   python - <<'PY'
   from pycsamt.pipeline import Pipeline, list_presets

   print("Presets:", list_presets())
   print(Pipeline.from_preset("basic_qc"))
   PY

If you installed the docs extra, check Sphinx:

.. code-block:: console
   :linenos:

   sphinx-build --version

If you installed the app extra, check the GUI entry point:

.. code-block:: console
   :linenos:

   pycsamt-desktop --help

Install For Common Roles
------------------------

Field-data inspection
    Install the core package.  Add ``geo`` if you need map projections or
    basemap context.

.. code-block:: console
   :linenos:

   python -m pip install "pycsamt[geo]"

Processing and inversion preparation
    Install the core package.  Add ``app`` if you want GUI inspection tools.

.. code-block:: console
   :linenos:

   python -m pip install "pycsamt[app]"

AI inversion experiments
    Install one deep-learning backend and set it explicitly.

.. code-block:: console
   :linenos:

   python -m pip install "pycsamt[torch]"
   PYCSAMT_AI_BACKEND=torch python my_training_script.py

Documentation contributor
    Install from source with development and documentation tools.

.. code-block:: console
   :linenos:

   python -m pip install -e ".[dev,docs]"

Full local workstation
    Install the full declared dependency set.  Use this for integration work,
    not for a minimal production environment.

.. code-block:: console
   :linenos:

   python -m pip install -e ".[full]"

Troubleshooting
---------------

``pycsamt`` command is not found
    The environment that installed pyCSAMT is not active, or the Python
    scripts directory is not on ``PATH``.  Try ``python -m pip show pycsamt``
    and reactivate the environment.

``ModuleNotFoundError: No module named 'pycsamt'``
    Install into the same interpreter that runs your script.  Prefer
    ``python -m pip install ...`` over a bare ``pip install ...`` so the pip
    executable is tied to the active interpreter.

An optional import fails
    Install the matching extra.  For example, use ``pycsamt[geo]`` for
    geospatial workflows, ``pycsamt[app]`` for the GUI, and ``pycsamt[torch]``
    or ``pycsamt[tensorflow]`` for neural AI workflows.

TensorFlow or PyTorch does not install
    Check your Python version and platform.  Machine-learning wheels are
    platform-specific and often support fewer Python versions than pyCSAMT
    itself.

Qt or PySide6 fails on Linux
    Install the required desktop system libraries, then retry the ``app``
    extra.  Conda-forge environments can be easier for GUI stacks.

Geospatial warnings mention ``GDAL_DATA`` or PROJ
    Core pyCSAMT workflows can still run, but map/projection workflows may be
    incomplete.  Install geospatial libraries from conda-forge or your system
    package manager.

The documentation build cannot fetch intersphinx inventories
    This usually means the machine is offline or DNS is restricted.  The
    local docs can still build, but links to external Python, NumPy, SciPy,
    pandas, Matplotlib, or scikit-learn objects may be incomplete.

Next Steps
----------

After installation:

* open :doc:`data_formats` to understand supported input files;
* open :doc:`configuration` to set API, plotting, pipeline, and agent
  defaults;
* open :doc:`first_survey` to load your first EDI directory and run basic QC;
* open :doc:`../agents/llm_configuration` if you plan to use LLM-assisted
  agents.
