Inspect and QC a Survey
=======================

This tutorial shows how to produce first-pass quality-control tables from a
loaded EDI survey.

Load a survey
-------------

.. code-block:: python

   from pycsamt.api import read_edis

   survey = read_edis("data/edis", progress="auto")
   sites = survey.collection

Build a QC table
----------------

Use :func:`pycsamt.emtools.qc.build_qc_table` for a compact station-level
diagnostic table:

.. code-block:: python

   from pycsamt.emtools.qc import build_qc_table

   qc = build_qc_table(sites, api=True)
   print(qc)

   df = qc.df
   print(df[["station", "n_freq", "frac_ok", "snr_med"]].head())

The exact columns depend on available transfer-function and tipper data, but
the table is intended to answer:

- how many frequencies each station contains
- how many transfer-function rows are finite
- whether tipper rows are available
- whether skew diagnostics can be estimated

Station confidence
------------------

For a normalized confidence view, use
:func:`pycsamt.emtools.qc.station_confidence_table`:

.. code-block:: python

   from pycsamt.emtools.qc import station_confidence_table

   confidence = station_confidence_table(sites, api=True)
   print(confidence)

Use higher-confidence stations for first inversion trials, and review low
confidence stations before applying automated processing.

CLI equivalent
--------------

.. code-block:: bash

   pycsamt edi info data/edis
   pycsamt edi validate data/edis
   pycsamt pipe run qc --input data/edis --output qc_results

The exact pipeline preset names are documented in :doc:`../pipeline/presets`.

See Also
--------
:doc:`read_edi_survey`
    Load an EDI survey.
:doc:`correct_static_shift`
    Apply a common first correction.
:doc:`../api/emtools`
    EMTools API reference.
