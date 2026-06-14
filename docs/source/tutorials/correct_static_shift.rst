Correct Static Shift
====================

This tutorial shows a conservative static-shift correction workflow using the
AMA-style helpers in :mod:`pycsamt.emtools.ss`.

Load and inspect data first
---------------------------

Static-shift correction should follow basic data inspection. Start with a
loaded survey and review station quality:

.. code-block:: python

   from pycsamt.api import read_edis
   from pycsamt.emtools.qc import station_confidence_table

   survey = read_edis("data/edis")
   sites = survey.collection

   confidence = station_confidence_table(sites, api=True)
   print(confidence)

Estimate correction factors
---------------------------

Estimate station-level correction factors:

.. code-block:: python

   from pycsamt.emtools.ss import estimate_ss_ama

   factors = estimate_ss_ama(
       sites,
       sort_by="lon",
       half_window=3,
       weights="tri",
       pband=None,
       max_skew=6.0,
       api=True,
   )

   print(factors)

The important output columns are:

- ``station``: station identifier
- ``delta_log10_rho``: estimated log-resistivity shift
- ``fac_rho``: apparent-resistivity scale factor
- ``fac_z``: impedance scale factor
- ``n_used``: number of samples used in the estimate

Apply the correction
--------------------

Apply the estimated correction to a copy of the survey:

.. code-block:: python

   from pycsamt.emtools.ss import correct_ss_ama

   corrected = correct_ss_ama(
       sites,
       sort_by="lon",
       half_window=3,
       weights="tri",
       inplace=False,
   )

Use ``inplace=False`` during exploration. Switch to ``inplace=True`` only when
you intentionally want to mutate the loaded site objects.

Review before inversion
-----------------------

Static-shift correction changes apparent-resistivity level. Always compare
before and after correction before preparing inversion files.

See Also
--------
:doc:`inspect_and_qc_survey`
    Inspect station quality before correction.
:doc:`prepare_occam2d_inversion`
    Use corrected data for inversion preparation.
:doc:`../theory/static_shift`
    Scientific background.
:doc:`../api/emtools`
    EMTools API reference.
