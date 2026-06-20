.. _tutorial_inspect_and_qc_survey:

Inspect and QC a Survey
=======================

This tutorial shows how to perform a first-pass quality-control review of an
EDI survey before correction, recomputation, modelling, or inversion. The goal
is not to decide the final geophysical interpretation. The goal is to answer a
more practical question:

*Which stations and frequency bands are reliable enough to use first, and which
ones need manual review?*

The workflow follows the style used throughout pyCSAMT v2:

- read the survey once
- inspect the station inventory
- build compact station-level QC tables
- compute confidence scores
- identify stations and frequency bands that need attention
- save review tables and diagnostic plots
- use CLI commands for quick checks where they are available

What You Will Learn
-------------------

After this tutorial you should be able to:

- load one EDI file, one survey directory, or a recursive survey tree
- inspect station names, frequency coverage, and missing transfer functions
- create a station-level quality-control table with
  :func:`pycsamt.emtools.qc.build_qc_table`
- generate simple flags with :func:`pycsamt.emtools.qc.qc_flags`
- compute station and frequency confidence tables
- export QC results for a field notebook, spreadsheet, or processing report
- choose which stations to keep, review, recompute, or remove from the first
  inversion trial

Input Assumptions
-----------------

The examples below assume that the EDI files are stored in ``data/edis``. The
path can be:

- a single ``.edi`` file
- a directory containing EDI files
- a survey directory with line subdirectories

For example:

.. code-block:: text

   data/
     edis/
       L18PLT/
         S001.edi
         S002.edi
       L22PLT/
         S101.edi
         S102.edi

When a directory contains several lines, use ``recursive=True`` so pyCSAMT
discovers EDI files below the top-level folder.

Load the Survey
---------------

Start by reading the EDI survey through the public API. ``strict=False`` is
useful during the first inspection because it allows the reader to continue
past recoverable metadata issues.

.. code-block:: python
   :linenos:

   from pycsamt.api import read_edis

   survey = read_edis(
       "data/edis",
       recursive=True,
       strict=False,
       progress="auto",
   )

   sites = survey.collection
   print(survey.summary())

``survey`` is the public survey object returned by the API. ``sites`` is the
site collection used by most station editing, computation, and QC helpers.

If you only need a pandas-friendly station inventory, inspect the survey
summary table:

.. code-block:: python
   :linenos:

   inventory = survey.df.to_pandas(copy=True)
   print(inventory.head())
   print(inventory.columns)

At this stage, check that the number of stations is what you expected and that
station names are not duplicated unintentionally.

Create an Output Folder
-----------------------

QC usually produces several small tables and plots. Put them in one folder so
the review is reproducible.

.. code-block:: python
   :linenos:

   from pathlib import Path

   qc_dir = Path("qc_review")
   qc_dir.mkdir(exist_ok=True)

Build the Station QC Table
--------------------------

The first compact table is produced by
:func:`pycsamt.emtools.qc.build_qc_table`. It summarizes the usable impedance
rows per station and, when possible, adds simple phase-tensor skew diagnostics.

.. code-block:: python
   :linenos:

   from pycsamt.emtools.qc import build_qc_table

   qc = build_qc_table(
       sites,
       include_skew=True,
       recursive=False,
       api=True,
   )

   qc_df = qc.to_pandas(copy=True)
   print(qc)
   print(qc_df.head())

The most useful columns are:

``station``
    Station name inferred from the EDI object.

``n_freq``
    Number of frequency rows available at the station.

``n_ok``
    Number of rows where the impedance tensor is finite enough for QC.

``frac_ok``
    Fraction of usable impedance rows, from ``0`` to ``1``.

``n_tip`` and ``n_tip_ok``
    Number of tipper rows and usable tipper rows when tipper data are present.

``snr_med``
    Median signal-to-noise proxy derived from impedance values and impedance
    errors when error tensors are available.

``pmin`` and ``pmax``
    Shortest and longest periods represented by the station.

``skew_med`` and ``skew_iqr``
    Median absolute phase-tensor skew and its interquartile spread. These
    columns are only added when ``include_skew=True``.

Export the table:

.. code-block:: python
   :linenos:

   qc_df.to_csv(qc_dir / "station_qc.csv", index=False)

Sort Stations by Review Priority
--------------------------------

A practical first review is to sort stations by coverage and signal quality.

.. code-block:: python
   :linenos:

   columns = [
       name for name in (
           "station",
           "n_freq",
           "frac_ok",
           "snr_med",
           "skew_med",
           "pmin",
           "pmax",
       )
       if name in qc_df.columns
   ]

   review_order = qc_df.sort_values(
       by=["frac_ok", "snr_med"],
       ascending=[True, True],
   )

   print(review_order[columns].head(15))

Stations near the top of this sorted table should be checked before automated
static-shift correction, dimensionality analysis, or inversion preparation.

Create Simple QC Flags
----------------------

Use :func:`pycsamt.emtools.qc.qc_flags` to attach simple rule-based labels.
The thresholds below are intentionally conservative for a first review. Adjust
them for local data quality, acquisition style, and survey objectives.

.. code-block:: python
   :linenos:

   from pycsamt.emtools.qc import qc_flags

   flagged = qc_flags(
       sites,
       min_frac_ok=0.75,
       min_snr_med=3.0,
       max_skew_med=6.0,
       recursive=False,
   )

   flagged.to_csv(qc_dir / "station_qc_flags.csv", index=False)
   print(flagged[["station", "frac_ok", "snr_med", "flags"]].head(20))

Typical flags are:

``low_coverage``
    Too few usable impedance rows. This can indicate incomplete spectra,
    parsing problems, bad frequency windows, or severe masking.

``low_snr``
    The median signal-to-noise proxy is below the selected threshold.

``high_skew``
    Phase-tensor skew is high enough to deserve dimensionality review.

A flag is not a deletion instruction. It is a review instruction. In field
datasets, low confidence can be caused by a real local conductor, poor electrode
contact, cultural noise, incorrect coordinate metadata, a rotation mismatch, or
format conversion from another software package.

Compute Station Confidence
--------------------------

The QC table is compact. The confidence table is more diagnostic. It combines
several indicators into a normalized confidence score between ``0`` and ``1``.

.. code-block:: python
   :linenos:

   from pycsamt.emtools.qc import station_confidence_table

   station_ci = station_confidence_table(
       sites,
       method="composite",
       relerr_threshold=0.20,
       offdiag_tolerance_log10=0.35,
       diagonal_leakage_max=0.35,
       phase_jump_tolerance_deg=90.0,
       spatial_tolerance_log10=0.60,
       spacing_m=200.0,
       recursive=False,
       api=True,
   )

   station_ci_df = station_ci.to_pandas(copy=True)
   station_ci_df.to_csv(qc_dir / "station_confidence.csv", index=False)
   print(station_ci_df.head())

Important confidence columns include:

``confidence``
    Composite score. Values close to ``1`` are usually safer for first-pass
    modelling. Values near ``0`` need review.

``confidence_err``
    Uncertainty proxy for the confidence score.

``coverage``
    Fraction of finite impedance rows.

``uncertainty``
    Score derived from impedance error tensors when available.

``offdiag``
    Consistency score between the two off-diagonal impedance components.

``diagonal``
    Score based on diagonal leakage relative to the off-diagonal components.

``phase``
    Score based on abrupt phase jumps.

``spatial``
    Neighbor-coherence score along the station profile.

Review Low-Confidence Stations
------------------------------

For first-pass work, a common pattern is to review stations below ``0.6`` and
start inversion tests with stations above ``0.8``. These are practical defaults,
not universal geophysical laws.

.. code-block:: python
   :linenos:

   low_ci = station_ci_df[station_ci_df["confidence"] < 0.60]
   high_ci = station_ci_df[station_ci_df["confidence"] >= 0.80]

   print("stations to review")
   print(low_ci[["station", "confidence", "coverage", "phase", "spatial"]])

   print("stations suitable for first trials")
   print(high_ci[["station", "confidence", "coverage"]].head())

   low_ci.to_csv(qc_dir / "stations_to_review.csv", index=False)

If many neighboring stations have low confidence in the same frequency band,
the problem may be survey-wide cultural noise or a real geologic response. If
one isolated station is poor across nearly all frequencies, inspect acquisition
metadata, electrode layout, orientation, and file conversion history.

Inspect Frequency-Level Confidence
----------------------------------

Station averages can hide narrow-band problems. Use
:func:`pycsamt.emtools.qc.frequency_confidence_table` to inspect every
station-frequency sample.

.. code-block:: python
   :linenos:

   from pycsamt.emtools.qc import frequency_confidence_table

   freq_ci = frequency_confidence_table(
       sites,
       method="composite",
       ci_hi=0.95,
       ci_lo=0.50,
       recursive=False,
       api=True,
   )

   freq_ci_df = freq_ci.to_pandas(copy=True)
   freq_ci_df.to_csv(qc_dir / "frequency_confidence.csv", index=False)

   weak_freq = freq_ci_df[freq_ci_df["confidence"] < 0.50]
   print(weak_freq[[
       "station",
       "frequency_hz",
       "period_s",
       "confidence",
       "flags",
   ]].head(20))

The frequency table is useful when you want to:

- mask a narrow noisy band rather than remove a whole station
- compare confidence between short-period and long-period data
- find stations with repeated phase jumps
- identify frequencies that may be affected by source-field instability,
  dead-band behavior, or cultural noise

Create a Confidence Profile Plot
--------------------------------

A line plot helps decide whether poor confidence is isolated, clustered, or
profile-wide.

.. code-block:: python
   :linenos:

   import matplotlib.pyplot as plt
   from pycsamt.emtools.qc import plot_confidence_profile

   ax = plot_confidence_profile(
       sites,
       method="composite",
       ci_hi=0.95,
       ci_lo=0.50,
       station_labels=True,
       spacing_m=200.0,
       recursive=False,
   )

   ax.figure.savefig(
       qc_dir / "station_confidence_profile.png",
       dpi=200,
       bbox_inches="tight",
   )
   plt.close(ax.figure)

Use this plot to distinguish three common cases:

- one or two isolated low-confidence stations, often worth manual repair
- a continuous low-confidence interval, often worth comparing with geology,
  topography, and acquisition notes
- profile-wide low confidence, often caused by import, rotation, unit, or
  configuration problems

Select Stations for the Next Step
---------------------------------

After QC, create explicit station lists. This makes later correction and
inversion runs easier to reproduce.

.. code-block:: python
   :linenos:

   keep_stations = (
       station_ci_df.loc[station_ci_df["confidence"] >= 0.60, "station"]
       .astype(str)
       .tolist()
   )

   review_stations = (
       station_ci_df.loc[station_ci_df["confidence"] < 0.60, "station"]
       .astype(str)
       .tolist()
   )

   (qc_dir / "keep_stations.txt").write_text(
       "\n".join(keep_stations) + "\n",
       encoding="utf-8",
   )
   (qc_dir / "review_stations.txt").write_text(
       "\n".join(review_stations) + "\n",
       encoding="utf-8",
   )

These lists can be used by later site selection, recomputation, static-shift
correction, or inversion preparation steps.

CLI Quick Checks
----------------

The Python API gives the richest QC tables. The CLI is useful for quick
terminal checks before opening a notebook.

.. code-block:: bash
   :linenos:

   pycsamt edi info data/edis
   pycsamt edi validate data/edis --deep
   pycsamt site info data/edis --format csv > qc_review/site_inventory.csv
   pycsamt site compute strike data/edis --format csv > qc_review/strike.csv
   pycsamt site compute tipper data/edis --format csv > qc_review/tipper.csv

Use the CLI outputs as supporting diagnostics. Use the Python QC tables when
you need confidence thresholds, frequency-level masking, or custom reporting.

What to Do With Poor Stations
-----------------------------

The QC result should lead to a processing decision. Common decisions are:

``keep``
    The station has high confidence and acceptable frequency coverage.

``review``
    The station has local problems but may be recoverable after frequency
    masking, component rotation, or static-shift analysis.

``recompute``
    The station was exported by another program or has inconsistent component
    orientation, naming, or metadata. Use pyCSAMT site recomputation before
    modelling.

``exclude from first trial``
    The station is too incomplete or unstable for the first inversion run.
    Keep it documented so it can be revisited after the first model explains
    the main survey response.

Troubleshooting
---------------

No stations are loaded
    Check the input path, file extension, and ``recursive=True`` when EDI files
    are inside line subdirectories.

All stations have low coverage
    Inspect the EDI format and component names. The files may need
    recomputation or conversion before QC.

Confidence is low only at long periods
    This may indicate noise, weak source field, dead-band behavior, or a real
    long-period response. Compare neighboring stations before masking.

Confidence is low only for one station
    Check station coordinates, orientation, contact resistance notes, and
    whether the station was exported differently from the rest of the line.

Skew is high across a profile segment
    Do not automatically delete those stations. Compare with structural
    geology, dimensionality diagnostics, and inversion residuals.

See Also
--------

:doc:`read_edi_survey`
    Load EDI files and inspect the survey object.

:doc:`correct_static_shift`
    Apply a common first correction after QC.

:doc:`../site/recompute`
    Recompute and rewrite EDI files imported from external software.

:doc:`../site/computed_diagnostics`
    Compute strike, resistivity, phase-slope, and tipper diagnostics.

:doc:`../api/emtools`
    EMTools API reference.
