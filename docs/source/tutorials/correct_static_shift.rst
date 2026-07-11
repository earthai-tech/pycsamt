.. _tutorial_correct_static_shift:

Correct Static Shift
====================

This tutorial shows a conservative static-shift correction workflow for EDI
surveys. It is written for the common profile-processing case: you have already
loaded and inspected a survey, you suspect that some apparent-resistivity
curves are shifted vertically, and you want to estimate station-level
correction factors before preparing inversion files.

Static-shift correction should be handled carefully. It changes impedance
amplitudes and apparent-resistivity levels. It does not fix bad phases,
frequency-dependent source effects, wrong station coordinates, incorrect
component orientation, or a poorly converted EDI file. Always inspect before
and after correction.

What You Will Learn
-------------------

After this tutorial you should be able to:

- diagnose when static-shift correction is worth trying
- estimate station-level correction factors with AMA
- inspect the meaning of ``delta_log10_rho``, ``fac_rho``, and ``fac_z``
- apply the correction without mutating the original survey
- compare before and after correction using static-shift QC plots
- export corrected EDI files and a correction manifest
- choose between AMA, LOESS, reference-median, and bilateral correction styles
- run static-shift correction from a pipeline configuration when needed

Prerequisites
-------------

Run first-pass survey inspection before correcting static shift:

.. code-block:: python
   :linenos:

   from pycsamt.api import read_edis
   from pycsamt.emtools.qc import build_qc_table, station_confidence_table

   survey = read_edis(
       "data/edis",
       recursive=True,
       strict=False,
       progress="auto",
   )
   sites = survey.collection

   qc = build_qc_table(sites, api=True)
   confidence = station_confidence_table(sites, method="composite", api=True)

   print(qc)
   print(confidence)

Static-shift correction is usually a second step. If the survey did not pass
basic loading, station-order, frequency-coverage, and quality checks, solve
those problems first. See :doc:`inspect_and_qc_survey`.

When Static Shift Is a Good Candidate
-------------------------------------

Static shift is plausible when:

- one station is shifted up or down relative to neighboring stations;
- the apparent-resistivity curve shape is similar to neighboring curves;
- phase does not show an equivalent anomaly;
- the offset is approximately frequency independent over the selected band;
- the station belongs to a profile where neighboring stations can define a
  local spatial trend.

Static-shift correction is risky when:

- the anomaly is frequency dependent;
- phase changes strongly at the same station;
- the station order or coordinates are wrong;
- several neighboring stations are shifted together because the geology is
  genuinely different;
- the data contain strong CSAMT source effects or near-field behavior;
- the EDI file came from another software package and still needs
  recomputation or component normalization.

For the scientific background, read :doc:`../theory/static_shift`.

Create a Review Folder
----------------------

Keep factors, plots, and corrected files together.

.. code-block:: python
   :linenos:

   from pathlib import Path

   review_dir = Path("static_shift_review")
   plot_dir = review_dir / "plots"
   corrected_dir = review_dir / "corrected_edis"

   plot_dir.mkdir(parents=True, exist_ok=True)
   corrected_dir.mkdir(parents=True, exist_ok=True)

Estimate AMA Correction Factors
-------------------------------

The standard interactive correction path uses
:func:`pycsamt.emtools.ss.estimate_ss_ama`. AMA means adaptive moving average:
each station is compared with its neighbors in log apparent-resistivity space.

.. code-block:: python
   :linenos:

   from pycsamt.emtools.ss import estimate_ss_ama

   factors = estimate_ss_ama(
       sites,
       sort_by="lon",
       half_window=3,
       weights="tri",
       pband=None,
       max_skew=6.0,
       robust_freq="median",
       robust_overall="median",
       recursive=False,
       api=True,
   )

   factor_df = factors.to_pandas(copy=True)
   factor_df.to_csv(review_dir / "static_shift_factors_ama.csv", index=False)
   print(factors)
   print(factor_df)

The factor table contains:

``station``
    Station identifier used to match the factor back to the EDI object.

``delta_log10_rho``
    Estimated log10 apparent-resistivity offset before correction. Positive
    values mean the station is high relative to its local trend. Negative
    values mean it is low.

``fac_rho``
    Apparent-resistivity scale factor. This is the multiplicative correction
    for apparent resistivity.

``fac_z``
    Impedance scale factor. pyCSAMT applies this to the impedance tensor
    because apparent resistivity is proportional to impedance amplitude
    squared.

``n_used``
    Number of frequency samples used to estimate the station factor.

Choose the Spatial Order
------------------------

AMA needs the profile order to make physical sense. The ``sort_by`` parameter
controls that order:

``sort_by="lon"``
    Use longitude. This is often reasonable for east-west lines.

``sort_by="lat"``
    Use latitude. This is often reasonable for north-south lines.

``sort_by="name"``
    Use station name order. This is useful when station names already encode
    line order and coordinates are missing or unreliable.

If the line is oblique, inspect the coordinates before deciding. A wrong order
can make the algorithm compare stations that are not neighbors in the field.

Use a Period Band
-----------------

Sometimes only part of the spectrum is appropriate for estimating static
shift. Use ``pband=(min_period, max_period)`` in seconds to restrict the
factor estimate.

.. code-block:: python
   :linenos:

   factors_mid_band = estimate_ss_ama(
       sites,
       sort_by="name",
       half_window=3,
       weights="tri",
       pband=(0.01, 10.0),
       max_skew=6.0,
       recursive=False,
       api=True,
   )

   factors_mid_band.to_pandas(copy=True).to_csv(
       review_dir / "static_shift_factors_ama_0p01s_10s.csv",
       index=False,
   )

Choose a band where the curves appear shifted but not strongly distorted. Avoid
using a band dominated by dead-band noise, transmitter source effects, or
obvious phase jumps.

Inspect Large Corrections
-------------------------

Large static-shift factors should be reviewed before they are applied.

.. code-block:: python
   :linenos:

   large = factor_df[
       (factor_df["fac_rho"] < 0.5)
       | (factor_df["fac_rho"] > 2.0)
       | (factor_df["n_used"] < 5)
   ]

   print("stations requiring manual review")
   print(large[["station", "delta_log10_rho", "fac_rho", "fac_z", "n_used"]])

A factor outside ``0.5`` to ``2.0`` is not automatically wrong, but it should
make you ask why the station is so different from its neighbors.

Apply Factors Manually
----------------------

The safest workflow is to estimate factors first, inspect them, then apply
them to a copy of the site collection with
:func:`pycsamt.emtools.ss.apply_ss_factors`.

.. code-block:: python
   :linenos:

   from pycsamt.emtools.ss import apply_ss_factors

   corrected = apply_ss_factors(
       sites,
       factors,
       key="fac_z",
       inplace=False,
       recursive=False,
       verbose=1,
   )

Use ``inplace=False`` during exploration. This returns a corrected copy and
keeps ``sites`` unchanged for comparison. Use ``inplace=True`` only inside a
controlled pipeline or after you have saved the original state.

Apply AMA in One Step
---------------------

For routine processing, :func:`pycsamt.emtools.ss.correct_ss_ama` combines
factor estimation and application:

.. code-block:: python
   :linenos:

   from pycsamt.emtools.ss import correct_ss_ama

   corrected = correct_ss_ama(
       sites,
       sort_by="lon",
       half_window=3,
       weights="tri",
       pband=None,
       max_skew=6.0,
       robust_freq="median",
       robust_overall="median",
       inplace=False,
       recursive=False,
       verbose=1,
   )

The one-step function is convenient, but the two-step factor workflow is better
when you are still learning the dataset.

Compare Before and After
------------------------

After correction, compare the original and corrected collections. A good
correction should reduce station-level vertical jumps without creating
frequency-dependent artifacts.

.. code-block:: python
   :linenos:

   import matplotlib.pyplot as plt
   from pycsamt.emtools.ss import (
       plot_ss_delta_profile,
       plot_ss_delta_psection,
       plot_ss_station_curves,
   )

   ax = plot_ss_delta_profile(
       sites,
       corrected,
       pband=None,
       robust="median",
   )
   ax.figure.savefig(plot_dir / "delta_profile.png", dpi=200, bbox_inches="tight")
   plt.close(ax.figure)

   ax = plot_ss_delta_psection(
       sites,
       corrected,
       axis_y="logperiod",
       pband=None,
   )
   ax.figure.savefig(plot_dir / "delta_psection.png", dpi=200, bbox_inches="tight")
   plt.close(ax.figure)

   ax = plot_ss_station_curves(
       sites,
       corrected,
       station=str(factor_df["station"].iloc[0]),
       pband=None,
   )
   ax.figure.savefig(plot_dir / "station_curve_before_after.png", dpi=200)
   plt.close(ax.figure)

The profile plot shows one correction value per station. The pseudo-section
shows whether the correction is approximately uniform across period. The
station-curve plot is useful for manual review of individual stations.

Use One-Call QC Wrappers
------------------------

The ``ss_qc_*`` helpers estimate, apply, and plot in one call. They are useful
for fast experiments. Use ``return_sites=True`` when you want the corrected
collection returned with the plot.

.. code-block:: python
   :linenos:

   from pycsamt.emtools.ss import ss_qc_psection, ss_qc_profile

   ax, corrected_preview = ss_qc_profile(
       sites,
       method="ama",
       return_sites=True,
       sort_by="lon",
       half_window=3,
       max_skew=6.0,
   )
   ax.figure.savefig(plot_dir / "ama_preview_profile.png", dpi=200)
   plt.close(ax.figure)

   ax = ss_qc_psection(
       sites,
       method="ama",
       sort_by="lon",
       half_window=3,
       max_skew=6.0,
       axis_y="logperiod",
   )
   ax.figure.savefig(plot_dir / "ama_preview_psection.png", dpi=200)
   plt.close(ax.figure)

These wrappers are excellent for exploration. For a formal processing record,
also save the factor table produced by ``estimate_ss_ama``.

Try Alternative Estimators
--------------------------

AMA is the recommended starting point, but pyCSAMT also exposes other factor
estimators:

.. code-block:: python
   :linenos:

   from pycsamt.emtools.ss import (
       estimate_ss_bilateral,
       estimate_ss_loess,
       estimate_ss_refmedian,
   )

   loess = estimate_ss_loess(
       sites,
       half_window=3,
       poly=1,
       it=2,
       pband=(0.01, 10.0),
       max_skew=6.0,
       api=True,
   )

   bilateral = estimate_ss_bilateral(
       sites,
       half_window=4,
       pband=(0.01, 10.0),
       max_skew=6.0,
       summary="median",
       api=True,
   )

   refmedian = estimate_ss_refmedian(
       sites,
       pband=(0.01, 10.0),
       max_skew=6.0,
       api=True,
   )

   loess.to_pandas(copy=True).to_csv(review_dir / "factors_loess.csv", index=False)
   bilateral.to_pandas(copy=True).to_csv(
       review_dir / "factors_bilateral.csv",
       index=False,
   )
   refmedian.to_pandas(copy=True).to_csv(
       review_dir / "factors_refmedian.csv",
       index=False,
   )

Use the alternatives as checks:

``estimate_ss_loess``
    Useful when you want a local smooth trend rather than a simple moving
    neighbor median.

``estimate_ss_refmedian``
    Useful when the whole profile can be compared to a survey-wide reference.

``estimate_ss_bilateral``
    Useful when you want spatial weighting that also respects value similarity.

If all methods point to the same shifted stations, confidence in the correction
increases. If the methods disagree strongly, inspect the data before applying a
large correction.

Export Corrected EDI Files
--------------------------

Once the correction is accepted, export the corrected collection as new EDI
files. Do not overwrite the raw EDI directory during tutorial or exploratory
work.

.. code-block:: python
   :linenos:

   from pycsamt.site.export import write_sites

   written = write_sites(
       corrected,
       corrected_dir,
       template="{station}_static_shift_corrected.edi",
       exist_ok=True,
       manifest_csv=review_dir / "corrected_edi_manifest.csv",
   )

   print(f"wrote {len(written)} corrected EDI files")

The manifest records which station was written to which file. Keep it with the
factor table and QC plots.

Run Through a Pipeline
----------------------

For repeatable production work, put static-shift correction in a pipeline
configuration. ``SS001`` is the AMA correction step.

.. code-block:: yaml
   :linenos:

   name: static_shift_review
   input: data/edis
   output_dir: results/static_shift_review

   steps:
     - name: qc_before
       code: QC001
     - name: static_shift
       code: SS001
       params:
         sort_by: lon
         half_window: 3
         weights: tri
         max_skew: 6.0
     - name: qc_after
       code: QC001

Inspect the registered static-shift steps from the command line:

.. code-block:: bash
   :linenos:

   pycsamt pipe steps --category static_shift

The pipeline registry currently includes:

``SS001``
    AMA correction with :func:`pycsamt.emtools.ss.correct_ss_ama`.

``SS002``
    LOESS factor estimation followed by factor application.

``SS003``
    Reference-median factor estimation followed by factor application.

``SS004``
    Bilateral factor estimation followed by factor application.

Use the pipeline when the same correction recipe must be shared, reviewed, or
re-run on several lines.

Checklist Before Inversion
--------------------------

Before using corrected data for inversion preparation, confirm that:

- the original EDI files are still available;
- the factor table has been saved;
- large correction factors have been manually reviewed;
- phase behavior is still geologically reasonable;
- correction plots show station-level offsets, not frequency-dependent
  distortion;
- corrected EDI files were exported into a separate directory;
- the inversion input preparation step points to the corrected directory, not
  the raw directory by accident.

Troubleshooting
---------------

All factors are exactly one
    The selected band may not contain usable samples, all stations may already
    be aligned, or the estimator could not build a neighbor trend.

One station receives an extreme factor
    Check whether it has enough frequencies, whether its coordinates are
    correct, and whether the station is an actual geologic outlier.

The pseudo-section correction is not period independent
    The problem may not be static shift. Inspect source effects, dead-band
    behavior, phase jumps, and frequency-dependent noise.

Corrected curves look worse
    Revisit station order, period band, ``half_window``, and skew threshold.
    Try ``sort_by="name"`` if coordinates are unreliable.

Correction changes interpretation too strongly
    Treat the correction as provisional. Run one inversion with raw data and
    one with corrected data, then compare residuals and geologic plausibility.

See Also
--------

:doc:`inspect_and_qc_survey`
    Inspect station quality before correction.

:doc:`prepare_occam2d_inversion`
    Prepare corrected data for inversion.

:doc:`../theory/static_shift`
    Scientific background and interpretation cautions.

:doc:`../user_guide/pipeline/steps`
    Pipeline static-shift step catalogue.

:doc:`../user_guide/site/recompute`
    Recompute and rewrite EDI files before correction when imports are
    inconsistent.

:doc:`../api/emtools`
    EMTools API reference.
