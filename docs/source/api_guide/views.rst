.. _api-view:

API Views
=========

Most dataframe-returning functions in pyCSAMT accept an ``api`` keyword
that controls what they return: the raw :class:`pandas.DataFrame`, or
:class:`~pycsamt.api.APIFrame` -- a thin, opt-in wrapper adding compact
printing, metadata, table profiling, and multi-table result containers on
top of the same underlying data. The parameter always defaults to
``api=None``, and what ``None`` resolves to is worth internalising before
anything else on this page:

.. code-block:: pycon

   >>> from pycsamt.api import read_edis
   >>> from pycsamt.emtools.qc import station_confidence_table

   >>> survey = read_edis(
   ...     "data/AMT/WILLY_DATA/L18PLT",
   ...     recursive=False,
   ...     strict=False,
   ...     progress=False,
   ... )
   >>> sites = survey.collection

   >>> default_out = station_confidence_table(sites)
   >>> type(default_out).__name__
   'APIFrame'

A bare call, with no ``api=`` argument at all, already returns an
``APIFrame`` -- not a plain ``DataFrame``. This is the package-wide
default (:data:`~pycsamt.api.PYCSAMT_API_VIEW`'s own default backend is
``"pycsamt"``, not ``"pandas"``), and it is the opposite of what "pandas
by default" would suggest, so it is worth stating plainly rather than
leaving a reader to discover it by checking a type. ``api=False`` is the
one value that reliably returns the raw dataframe, regardless of the
global setting:

.. code-block:: pycon

   >>> plain = station_confidence_table(sites, api=False)
   >>> type(plain).__name__
   'DataFrame'
   >>> plain[["station", "confidence", "coverage"]].head(3).to_string(index=False)
   'station  confidence  coverage\n18-001A    0.709038       1.0\n18-002U    0.774634       1.0\n18-003A    0.713303       1.0'

Opt In Explicitly With ``api=True``
-----------------------------------

Passing ``api=True`` also produces an ``APIFrame``, and does so
regardless of the global setting -- the difference from the bare call
above only matters once the global default has been changed (see
`Global Configuration`_ below):

.. code-block:: pycon

   >>> view = station_confidence_table(sites, api=True)
   >>> view.name, view.kind, view.df.shape
   ('station_confidence_table', 'emtools.qc.station_confidence', (28, 13))
   >>> view.description
   'Station-level composite confidence scores.'

The same object exposes the underlying data through several angles at
once, without giving up pandas access:

.. code-block:: pycon

   >>> type(view.df).__name__, type(view.data).__name__
   ('DataFrame', 'ndarray')
   >>> list(view.columns)[:4]
   ['station', 'distance_m', 'confidence', 'confidence_err']
   >>> view.stats.rows, view.stats.columns, view.stats.missing_fraction
   (28, 13, 0.0)
   >>> view.confidence.to_numpy()[:3].round(6)
   array([0.709038, 0.774634, 0.713303])

``view.df`` is the real pandas object, ``view.data`` its NumPy values,
``view.stats`` a compact :class:`~pycsamt.api.FrameProfile`, and a safe
column name (``view.confidence``, not ``view.n/freq`` -- that one still
needs ``view["n_freq"]``) is available as a direct attribute.

Ordinary pandas methods pass straight through to the underlying
dataframe, but the result they return is the *plain* pandas object, not
another ``APIFrame`` -- ``query``/``groupby`` are pandas operations that
know nothing about the wrapper around their input:

.. code-block:: pycon

   >>> high = view.query("confidence >= 0.7")
   >>> type(high).__name__, len(high)
   ('DataFrame', 12)
   >>> view.groupby("method")["confidence"].median()
   method
   composite    0.671735
   Name: confidence, dtype: float64

Wrap the result again explicitly (``station_confidence_table(sites, api=True).df.query(...)``
piped back through :func:`~pycsamt.api.wrap_frame` if one is needed) when
a chained operation's output should stay an ``APIFrame`` too; it does not
happen automatically.

Common Table Functions
----------------------

The same ``api=True`` pattern works identically across every table
function in :mod:`pycsamt.emtools`, :mod:`pycsamt.metadata`, and beyond
-- learning it once on ``station_confidence_table`` covers all of them:

.. code-block:: pycon

   >>> from pycsamt.emtools.inspect import sites_summary
   >>> from pycsamt.emtools.qc import build_qc_table, frequency_confidence_table
   >>> from pycsamt.metadata.quality import quality_dataframe
   >>> from pycsamt.metadata.geology import CATALOG

   >>> sites_summary(sites, api=True).df.shape
   (28, 7)
   >>> build_qc_table(sites, api=True).df.shape
   (28, 11)
   >>> frequency_confidence_table(sites, api=True).df.shape
   (1484, 18)
   >>> quality_dataframe(sites, api=True).df.shape
   (28, 6)
   >>> CATALOG.to_dataframe(api=True).df.shape
   (13, 10)

``frequency_confidence_table`` has far more rows than the others
(1484 = 28 stations x 53 frequencies) because it reports one row per
station-frequency pair rather than one row per station.

Multi-Table Results
-------------------

Workflows that naturally produce several related outputs return
:class:`~pycsamt.api.APIResult` instead, when ``api=True``. Frequency
editing is a good example -- it returns the edited sites alongside report
and decision tables in one object:

.. code-block:: pycon

   >>> from pycsamt.emtools.frequency import edit_frequencies_by_confidence

   >>> result = edit_frequencies_by_confidence(
   ...     sites, mode="recover", method="presence", api=True,
   ... )
   >>> result.kind
   'emtools.frequency.edit'
   >>> type(result.sites).__name__
   'Sites'
   >>> type(result.report).__name__, type(result.decisions).__name__
   ('APIFrame', 'APIFrame')
   >>> result.report.df.shape, result.decisions.df.shape
   ((28, 18), (1484, 10))
   >>> result.n_dropped, result.n_recovered
   (0, 0)

Both counts are genuinely zero for this survey -- ``L18PLT`` already has
full frequency coverage at every station, so there is nothing for
``mode="recover"`` to drop or recover here; the real value of running it
is confirming that, not producing a dramatic before/after. ``result.sites``
carries the (unmodified, in this case) edited collection forward as an
ordinary :class:`~pycsamt.site.Sites` object, not wrapped, since that is
what downstream processing steps expect regardless of the API view
setting. A bare call, with no ``api=`` argument, behaves exactly like the
first ``station_confidence_table`` example on this page rather than
reverting to the function's original, pre-API-view return type:

.. code-block:: pycon

   >>> result2 = edit_frequencies_by_confidence(sites, mode="recover", method="presence")
   >>> type(result2).__name__
   'APIResult'

That last line is not a typo: this particular function's own default
already resolves through the same enabled-by-default global config as
every dataframe function above, so the *bare* call is an ``APIResult``
here too, for the same reason the very first example on this page was
already an ``APIFrame``. Pass ``api=False`` if the plain,
non-``APIResult`` return type is what a script actually needs.

Global Configuration
--------------------

:data:`~pycsamt.api.PYCSAMT_API_VIEW` is the single switch every function
above reads when ``api=None`` (the default):

.. code-block:: pycon

   >>> from pycsamt.api import PYCSAMT_API_VIEW, configure_api_view, reset_api_view

   >>> print(PYCSAMT_API_VIEW)
   APIViewConfig(backend='pycsamt')

Disable wrapping globally to make bare calls return plain pandas again --
useful in a script written before the API view layer existed, or one that
explicitly wants the old behaviour everywhere without touching every call
site:

.. code-block:: pycon

   >>> configure_api_view(backend=False)
   >>> out = station_confidence_table(sites)
   >>> type(out).__name__
   'DataFrame'

An explicit ``api=True`` still overrides a disabled global setting --
"force wrapping" means exactly that, not "wrap only if the global switch
allows it":

.. code-block:: pycon

   >>> still_forced = station_confidence_table(sites, api=True)
   >>> type(still_forced).__name__
   'APIFrame'

   >>> reset_api_view()
   >>> print(PYCSAMT_API_VIEW)
   APIViewConfig(backend='pycsamt')

The same switch is available as an environment variable, read once at
import time, for CI and batch environments that should not need a
``configure_api_view()`` call in every entry-point script:

.. code-block:: bash

   PYCSAMT_API_VIEW=pandas python workflow.py
   PYCSAMT_API_VIEW=pycsamt python workflow.py

Custom Wrappers
---------------

Advanced users can replace the ``APIFrame`` wrapper entirely with their
own callable. It receives the dataframe and the same metadata pyCSAMT
would otherwise attach to an ``APIFrame`` (``name``, ``kind``, ``source``,
and others) as keyword arguments:

.. code-block:: pycon

   >>> def my_table(df, **meta):
   ...     return {
   ...         "dataframe": df,
   ...         "name": meta.get("name"),
   ...         "kind": meta.get("kind"),
   ...         "source": meta.get("source"),
   ...     }
   ...
   >>> configure_api_view(wrapper=my_table)
   >>> table = station_confidence_table(sites, api=True)
   >>> type(table).__name__, table["kind"]
   ('dict', 'emtools.qc.station_confidence')

   >>> reset_api_view()

A configured custom wrapper applies to both an ``api=None`` bare call and
an explicit ``api=True`` call alike -- only ``api=False`` bypasses it, by
returning the raw dataframe directly with no wrapper involved at all.

Reading EDI Collections
-----------------------

:func:`~pycsamt.api.read_edis` is the public entry point for turning field
data into a survey object, and follows the exact same ``api``-aware
pattern as every function above -- ``survey`` here is already the
default-wrapped result:

.. code-block:: pycon

   >>> print(survey)
   APISurvey: edi_survey
   sites: 28
   stations: 23-18-001A, 23-18-002U, 23-18-003A, 23-18-004A, 23-18-005U, 23-18-006A, 23-18-007U, 23-18-008U, ...
   source: data/AMT/WILLY_DATA/L18PLT

   >>> survey.collection is survey.data
   True
   >>> type(survey.collection).__name__
   'EDICollection'

``survey.collection`` and ``survey.data`` are the same
:class:`~pycsamt.seg.collection.EDICollection` object under two names --
``collection`` reads better in EDI-specific code, ``data`` matches the
generic accessor every other API view container also exposes. Progress
display is controlled per call, independently of the ``api`` setting:

.. code-block:: pycon

   >>> _ = read_edis("data/AMT/WILLY_DATA/L18PLT", progress=False)

Practical Rule
--------------

Use ``api=False`` explicitly in library code and scripts that must work
regardless of whatever a caller's environment or notebook session has
configured globally:

.. code-block:: pycon

   >>> df = station_confidence_table(sites, api=False)
   >>> type(df).__name__
   'DataFrame'

Use ``api=True`` explicitly when presenting, inspecting, exporting, or
composing workflow results interactively, and rely on the bare
``api=None`` default only in throwaway scripts where the current global
setting -- whatever it happens to be -- is an acceptable outcome either
way.

Next Steps
----------

* :doc:`overview` for how the view layer fits alongside every other
  :mod:`pycsamt.api` configuration family.
* :doc:`configuration` for the systematic tour of every other family,
  including ``PYCSAMT_API_VIEW``'s own row in the family table.
* :doc:`../api/api` for the complete generated reference of
  :class:`~pycsamt.api.APIFrame`, :class:`~pycsamt.api.APIResult`, and
  every function this page used.
