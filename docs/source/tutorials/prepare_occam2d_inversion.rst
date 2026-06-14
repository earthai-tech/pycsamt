Prepare an Occam2D Inversion
============================

This tutorial outlines the pyCSAMT v2 path for preparing a 2-D inversion using
the Occam2D backend. Occam2D preparation is normally done after loading,
quality control, frequency editing, and static-shift review.

.. note::

   Occam2D is an external inversion workflow. pyCSAMT can prepare and manage
   configuration and files, but actual execution depends on the external
   Occam2D program being available in the working environment.

Load and prepare survey data
----------------------------

Start with a loaded EDI survey:

.. code-block:: python

   from pycsamt.api import read_edis

   survey = read_edis("data/edis")
   sites = survey.collection

Run first-pass QC and optional static-shift correction:

.. code-block:: python

   from pycsamt.emtools.qc import station_confidence_table
   from pycsamt.emtools.ss import correct_ss_ama

   confidence = station_confidence_table(sites, api=True)
   print(confidence)

   corrected = correct_ss_ama(sites, inplace=False)

Create an inversion workflow
----------------------------

Use the backend-neutral inversion workflow API to select the Occam2D backend:

.. code-block:: python

   from pycsamt.inversion.workflow import InversionWorkflow

   workflow = InversionWorkflow(
       method="mt",
       dimension="2d",
       backend="occam2d",
       data=corrected,
       workdir="occam2d_work",
       run_external=False,
   )

   result = workflow.run()

``run_external=False`` is the conservative preparation mode for external
backends. It allows configuration and file generation without launching an
external program. Backend-specific options can be passed through
``backend_options``:

.. code-block:: python

   workflow = InversionWorkflow(
       method="mt",
       dimension="2d",
       backend="occam2d",
       data=corrected,
       workdir="occam2d_work",
       run_external=False,
       backend_options={
           "profile_name": "line01",
           "target_rms": 1.0,
       },
   )

Expected outputs
----------------

The exact file names depend on the backend configuration, but an Occam2D
preparation workflow typically creates:

- data file
- mesh or model file
- startup/configuration file
- run metadata
- warnings or validation messages

Inspect the backend-neutral result object:

.. code-block:: python

   print(result)
   print(result.files)
   print(result.warnings)

CLI equivalent
--------------

The inversion CLI provides build/run/status/result commands:

.. code-block:: bash

   pycsamt invert build --backend occam2d --input data/edis --output occam2d_work
   pycsamt invert status occam2d_work
   pycsamt invert results occam2d_work

See Also
--------
:doc:`read_edi_survey`
    Load input EDI files.
:doc:`inspect_and_qc_survey`
    Review data quality before inversion.
:doc:`correct_static_shift`
    Review static-shift correction.
:doc:`../api/inversion`
    Inversion API reference.
:doc:`../cli/invert`
    Inversion CLI reference.
