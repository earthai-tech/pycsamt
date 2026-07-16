.. _cli:

CLI
===

The pyCSAMT command-line interface is the shell entry point for the same
survey workflows exposed by the Python API: data inspection, format
conversion, site inventory, processing pipelines, forward modelling,
inversion preparation, interpretation, mapping, and TDEM utilities.

Use the CLI when you need repeatable commands for a field dataset,
batch-friendly outputs for automation, or a quick way to validate that a
survey directory can be loaded before writing Python code.

.. code-block:: console

   pycsamt --help
   pycsamt edi validate data/edis/
   pycsamt site info data/edis/
   pycsamt map stations data/edis/

.. toctree::
   :numbered: 4
   :maxdepth: 3
   :class: pycsamt-guide-toc

   overview
   config
   survey
   edi
   avg
   jones
   site
   transform
   forward
   invert
   map
   interp
   pipe
   tdem
