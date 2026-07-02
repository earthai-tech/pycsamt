.. _api-view:

API Views
=========

pyCSAMT keeps its scientific functions pandas-friendly by default.  Functions
that historically returned :class:`pandas.DataFrame` still do so unless you
explicitly ask for the public API view layer.

The view layer is an opt-in presentation and workflow layer around tabular
results.  It gives you compact printing, metadata, table profiling, dataframe
access, NumPy conversion, and multi-table result containers without removing
pandas behavior.

Default Behavior Stays Pandas
-----------------------------

Internal and workflow functions keep their old return type:

.. code-block:: python

   from pycsamt.emtools.qc import station_confidence_table

   df = station_confidence_table(sites)
   type(df).__name__
   # 'DataFrame'

Opt In With ``api=True``
------------------------

Use ``api=True`` when you want a richer object for display or downstream
interactive work:

.. code-block:: python

   view = station_confidence_table(sites, api=True)

   print(view)
   # APIFrame: station_confidence_table
   # kind: emtools.qc.station_confidence
   # shape: ...

   view.df          # the underlying pandas DataFrame
   view.data        # NumPy array
   view.columns     # pandas columns
   view.stats       # compact pyCSAMT profile
   view.confidence  # column access when the column name is safe

``APIFrame`` delegates ordinary dataframe operations to pandas where possible:

.. code-block:: python

   high = view.query("confidence >= 0.9")
   grouped = view.groupby("method")["confidence"].median()

Common Table Functions
----------------------

The API view layer is available on selected dataframe-returning functions:

.. code-block:: python

   from pycsamt.emtools.inspect import sites_summary
   from pycsamt.emtools.qc import (
       build_qc_table,
       station_confidence_table,
       frequency_confidence_table,
   )
   from pycsamt.emtools.frequency import (
       frequency_edit_report,
       frequency_edit_decision_table,
   )
   from pycsamt.metadata.quality import quality_dataframe
   from pycsamt.metadata.geology import CATALOG

   sites_summary(sites, api=True)
   build_qc_table(sites, api=True)
   station_confidence_table(sites, api=True)
   frequency_confidence_table(sites, api=True)
   quality_dataframe(sites, api=True)
   CATALOG.to_dataframe(api=True)

Multi-Table Results
-------------------

Workflows that naturally produce several related outputs can return
``APIResult`` when ``api=True``.  For example, frequency editing returns the
edited sites plus report and decision tables:

.. code-block:: python

   from pycsamt.emtools.frequency import edit_frequencies_by_confidence

   result = edit_frequencies_by_confidence(
       sites,
       mode="recover",
       method="presence",
       api=True,
   )

   print(result)
   # APIResult: frequency_edit
   # items: sites, report, decisions, ...

   result.sites       # edited site collection
   result.report      # APIFrame
   result.decisions   # APIFrame
   result.n_dropped
   result.n_recovered

The default call still returns the original workflow result:

.. code-block:: python

   result = edit_frequencies_by_confidence(sites)
   type(result).__name__
   # 'FrequencyEditResult'

Global Configuration
--------------------

The global API view policy controls what ``api=True`` means.

The default backend is pyCSAMT's built-in ``APIFrame`` / ``APIResult`` view:

.. code-block:: python

   from pycsamt.api import PYCSAMT_API_VIEW

   PYCSAMT_API_VIEW
   # APIViewConfig(backend='pycsamt')

Disable wrapping globally when you want old behavior even in code paths that
pass ``api=True``:

.. code-block:: python

   from pycsamt.api import configure_api_view, reset_api_view

   configure_api_view(backend=False)

   out = station_confidence_table(sites, api=True)
   type(out).__name__
   # 'DataFrame'

   reset_api_view()

You can also use an environment variable before importing pyCSAMT:

.. code-block:: bash

   PYCSAMT_API_VIEW=pandas python workflow.py
   PYCSAMT_API_VIEW=pycsamt python workflow.py

Custom Wrappers
---------------

Advanced users can provide their own dataframe wrapper.  The callable receives
the dataframe-like object and the metadata passed by pyCSAMT:

.. code-block:: python

   from pycsamt.api import configure_api_view

   def my_table(df, **meta):
       return {
           "dataframe": df,
           "name": meta.get("name"),
           "kind": meta.get("kind"),
           "source": meta.get("source"),
       }

   configure_api_view(wrapper=my_table)

   table = station_confidence_table(sites, api=True)
   table["kind"]
   # 'emtools.qc.station_confidence'

Reading EDI Collections
-----------------------

The public reader helpers return API-facing survey objects while preserving the
underlying EDI collection:

.. code-block:: python

   from pycsamt.api import read_edis

   survey = read_edis("data/**/*.edi", progress="auto")

   print(survey)
   survey.collection    # underlying EDICollection
   survey.data          # same collection
   survey.summary()     # APIFrame by default

Progress display is controlled per call:

.. code-block:: python

   read_edis("data/*.edi", progress=True, leave=True)
   read_edis("data/*.edi", progress=False)

Practical Rule
--------------

Use plain calls when you are writing library code or scripts that expect pandas:

.. code-block:: python

   df = station_confidence_table(sites)

Use ``api=True`` when you are presenting, inspecting, exporting, or composing
workflow results interactively:

.. code-block:: python

   result = edit_frequencies_by_confidence(sites, api=True)
   print(result.report)
