.. _user_guide_data:

Data loading
============

Data loading is the first step in most pyCSAMT workflows. In v2, the preferred
public entry points return API-facing survey objects while preserving access to
the lower-level file collections used by the processing modules.

EDI surveys
-----------

EDI is the primary interchange format for MT, AMT, and CSAMT transfer-function
data. Use :func:`pycsamt.api.read_edis` to read a directory, glob pattern, or
list of files:

.. code-block:: python

   from pycsamt.api import read_edis

   survey = read_edis("data/edis/**/*.edi", recursive=True)

   print(survey)
   sites = survey.collection

The returned survey is an API view around the underlying EDI collection. This
means it is convenient for notebooks and reporting, but it still gives access
to the collection expected by lower-level processing functions.

Inspecting a survey
-------------------

Most first checks should answer three questions:

- Which stations were loaded?
- How many frequencies are available?
- Are transfer-function and tipper values present?

Use the survey summary when you want a compact table:

.. code-block:: python

   summary = survey.summary()
   print(summary)

When the API view layer is enabled, table helpers can return rich
``APIFrame`` objects. Use the ``.df`` attribute when a plain pandas dataframe
is needed:

.. code-block:: python

   table = survey.summary()
   df = table.df if hasattr(table, "df") else table

Duplicate stations and strict parsing
-------------------------------------

Field exports sometimes contain repeated station names, incomplete files, or
mixed directory contents. The public reader exposes conservative controls:

.. code-block:: python

   survey = read_edis(
       "data/edis",
       recursive=True,
       strict=False,
       on_dup="replace",
       progress="auto",
   )

``strict=False`` keeps a workflow moving when non-critical files are malformed.
Use ``strict=True`` when preparing a reproducible production workflow and you
want malformed input to fail early.

Other supported formats
-----------------------

pyCSAMT v2 also includes packages for AVG, J, SEG-style files, Zonge exports,
and TDEM workflows. These are documented in the format-specific sections:

- :mod:`pycsamt.jones`
- :mod:`pycsamt.zonge`
- :mod:`pycsamt.seg`
- :mod:`pycsamt.tdem`
- :mod:`pycsamt.transformers`

See Also
--------
:doc:`../tutorials/read_edi_survey`
    First EDI loading tutorial.
:doc:`../api/api`
    Public API helpers, including ``read_edis``.
:doc:`../api/seg`
    EDI and SEG-style parsing internals.
