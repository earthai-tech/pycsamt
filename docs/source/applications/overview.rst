.. _applications_overview:

Overview
========

pyCSAMT includes four user-facing applications built on the same
scientific core — the same readers, processing tools, plotting
conventions, and configuration system you use from Python. They differ
only in working style: a native desktop GUI for local interactive
review, a browser app for shared servers and team demonstrations, a
conversational surface that delegates workflows to the pyCSAMT agents,
and a dedicated map workbench for seeing a survey in space. Results are
interchangeable — a survey processed in one surface can be picked up in
any other, or in plain Python.

Install The App Extra
------------------------

The application surfaces use optional GUI and web dependencies:

.. code-block:: bash

   pip install "pycsamt[app]"

For development from a source checkout:

.. code-block:: bash

   pip install -e ".[app,dev]"
