.. _applications_overview:

Overview
========

pyCSAMT includes five user-facing applications built on the same
scientific core — the same readers, processing tools, plotting
conventions, and configuration system you use from Python. Four of them
differ only in working style: a native desktop GUI for local interactive
review, a browser app for shared servers and team demonstrations, a
conversational surface that delegates workflows to the pyCSAMT agents,
and a dedicated map workbench for seeing a survey in space. Results are
interchangeable — a survey processed in one surface can be picked up in
any other, or in plain Python.

The fifth, the standalone :doc:`Format Converter <converter>`, is
different in kind rather than in working style: it has no survey-loading
or processing features at all, only file-format conversion (inversion
results, EDI/EMTF-XML, boreholes, geology legends, structural evidence,
points of interest), which keeps it small enough to freeze into a
double-click binary for someone who does not want a Python environment.

Install The App Extra
------------------------

The application surfaces use optional GUI and web dependencies:

.. code-block:: bash

   pip install "pycsamt[app]"

For development from a source checkout:

.. code-block:: bash

   pip install -e ".[app,dev]"
