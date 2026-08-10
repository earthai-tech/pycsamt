.. _getting-started-installation:

Install pyCSAMT
===============

This page takes you from an available Python interpreter to a verified core
pyCSAMT installation. The core package is enough for the first survey workflow
in this section. You can add optional features later, when a workflow actually
needs them.

For supported versions, every optional dependency group, source and conda
installs, application commands, compiled solvers, and troubleshooting, use the
complete :doc:`installation reference <../installation>`.

Check Python
------------

pyCSAMT requires Python 3.9 or later. Python 3.11 or 3.12 is a practical
choice for a new environment because both versions are covered by the test
suite and are widely supported by scientific Python packages.

Check the interpreter that will create the environment:

.. code-block:: console

   python --version

If that command reports an older version, install a supported Python release
before continuing.

Create An Isolated Environment
------------------------------

An isolated environment keeps pyCSAMT and its scientific dependencies from
changing packages used by other projects. Create one inside your project or
working directory.

On Linux or macOS:

.. code-block:: console

   python -m venv .venv
   source .venv/bin/activate
   python -m pip install --upgrade pip

On Windows PowerShell:

.. code-block:: powershell

   py -3.11 -m venv .venv
   .\.venv\Scripts\Activate.ps1
   python -m pip install --upgrade pip

Your shell prompt will usually show ``(.venv)`` after activation. Confirm that
``python`` and pip now belong to the same environment:

.. code-block:: console

   python -m pip --version

Install The Core Package
------------------------

Install the stable package from PyPI:

.. code-block:: console

   python -m pip install pycsamt

Using ``python -m pip`` ties the installer to the active interpreter and
avoids a common problem in which ``pip`` installs into a different Python
environment.

Verify The Installation
-----------------------

First verify the Python package:

.. code-block:: console

   python -c "import pycsamt; print(pycsamt.__version__)"

The command should print a version without a traceback. Then verify the
command-line interface:

.. code-block:: console

   pycsamt --help

The help output should list the available command groups. If ``pycsamt`` is
not found, reactivate ``.venv`` and confirm the install location with:

.. code-block:: console

   python -m pip show pycsamt

Do not continue until both the import and command-line checks succeed. This
prevents environment problems from being mistaken for data-loading or survey
errors later in the workflow.

Add Optional Features Only When Needed
--------------------------------------

The core installation is the recommended starting point. Add an extra only
when the next workflow requires it:

.. list-table::
   :header-rows: 1
   :widths: 24 38 38

   * - Need
     - Command
     - Extra
   * - Coordinate transforms and map layers
     - ``python -m pip install "pycsamt[geo]"``
     - ``geo``
   * - PyTorch-backed AI workflows
     - ``python -m pip install "pycsamt[torch]"``
     - ``torch``
   * - Desktop and web applications
     - ``python -m pip install "pycsamt[app]"``
     - ``app``

These are examples, not the complete dependency matrix. See the optional
feature groups in :doc:`the installation reference <../installation>` before
setting up TensorFlow, LLM providers, documentation tools, a full development
environment, or an external inversion solver.

Continue With Your Data
-----------------------

Your environment is ready when:

* the virtual environment is active;
* importing ``pycsamt`` prints a version;
* ``pycsamt --help`` displays the command groups.

Continue to :doc:`data_formats` to identify the files from your instrument or
processing system. After that, :doc:`configuration` establishes output and
style defaults, and :doc:`first_survey` loads and checks a first survey.

If you are contributing to pyCSAMT rather than installing it for use, follow
the source-install instructions in :doc:`../installation` and then read
:doc:`../contributing`.
