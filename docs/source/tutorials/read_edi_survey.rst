Read an EDI Survey
==================

This tutorial shows the shortest path from a directory of EDI files to a
survey object that can be inspected, summarized, and passed into processing
workflows.

Use case
--------

You have one directory containing MT, AMT, or CSAMT EDI files and want to check
what pyCSAMT can read before doing quality control or inversion preparation.

Read the files
--------------

Use :func:`pycsamt.api.read_edis` for the public v2 API:

.. code-block:: python

   from pycsamt.api import read_edis

   survey = read_edis(
       "data/edis/**/*.edi",
       recursive=True,
       strict=False,
       on_dup="replace",
       progress="auto",
   )

   print(survey)

The returned object is an ``APISurvey``. It wraps the underlying EDI
collection, so both interactive use and lower-level processing remain
available:

.. code-block:: python

   collection = survey.collection
   raw_data = survey.data

Summarize the survey
--------------------

Generate a compact station summary:

.. code-block:: python

   summary = survey.summary()
   print(summary)

If a plain dataframe is needed:

.. code-block:: python

   df = summary.df if hasattr(summary, "df") else summary
   print(df.head())

Turn progress off
-----------------

For scripts and tests, disable progress output:

.. code-block:: python

   survey = read_edis("data/edis", progress=False)

CLI equivalent
--------------

The CLI section documents command-line workflows for EDI inspection:

.. code-block:: bash

   pycsamt edi info data/edis
   pycsamt edi stations data/edis
   pycsamt edi validate data/edis

Next steps
----------

- Run quality-control tables with :doc:`inspect_and_qc_survey`.
- Correct static shift with :doc:`correct_static_shift`.
- Prepare inversion files with :doc:`prepare_occam2d_inversion`.

See Also
--------
:doc:`../user_guide/data_loading`
    Data loading concepts.
:doc:`../api/api`
    Public API reference.
:doc:`../cli/edi`
    EDI command reference.
