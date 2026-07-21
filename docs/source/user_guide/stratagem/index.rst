.. _user_guide_stratagem:

Stratagem Surveys
=================

The :mod:`pycsamt.stratagem` package covers the post-acquisition workflow
for AMT surveys collected with Geometrics/EMI Stratagem hardware, from raw
instrument output through to processed, georeferenced EDI files ready for
:doc:`../emtools/index` or inversion. Stratagem hardware does not export
EDI directly: each station is recorded as a set of whitespace-separated
19-column ASCII component files (frequency, stack count, cross-spectral
values), and a separate desktop tool (WinGLink) converts those into EDI
files with placeholder coordinates (``LAT=0:00:00.00``). Everything from
that point on — reading the raw files for hardware-level QC, injecting
real station coordinates from a GPS table, statistical QC, static-shift
correction, frequency editing, noise removal, and final export/renaming —
is what this package handles.

.. code-block:: text

   raw hardware files (X*.NNN, Y*.NNN, Z*.NNN)
     ↓  WinGLink (external tool, required once)
       ↓  EDIBatch / Sites        load & natural-sort the EDI export
         ↓  CoordinateInjector    inject GPS coords from CSV / XLS / XLSX
           ↓  QualityController   per-station QC report (optional)
             ↓  StaticShiftCorrector  AMA static-shift removal
               ↓  FrequencyFilter     band select + incoherent-freq masking
                 ↓  NoiseRemover      powerline notch + Hampel + smoothing
                   ↓  export() / rename()   write corrected, renamed EDIs

:class:`~pycsamt.stratagem.survey.StratagemSurvey` composes every step
above into a single fluent pipeline, but each stage is also usable on its
own — useful when a survey needs a non-standard station-count reconciliation,
a custom coordinate table layout, or only a subset of the processing steps.

Real Stratagem deliveries are rarely clean: GPS tables often don't have one
row per station (calibration shots, repeat "check point" measurements,
stations the instrument failed to record), coordinate columns can be
mislabeled or shifted, and a bad GPS fix can sit undetected among hundreds
of otherwise-good rows. This guide works through those cases directly
rather than only the happy path.

Use this section when ingesting a new Stratagem field delivery, when GPS
coordinates need reconciling against a field report before
:class:`~pycsamt.stratagem.gis_correct.CoordinateInjector` can run, or when
composing the full raw-to-corrected-EDI pipeline with
:class:`~pycsamt.stratagem.survey.StratagemSurvey`. :doc:`tutorial` runs
that same workflow as one visual, reproducible case study, using
:mod:`pycsamt.map` and :mod:`pycsamt.emtools` to check the result at
each step rather than only reading printed numbers.

.. toctree::
   :maxdepth: 3
   :class: pycsamt-guide-toc

   concepts
   loading
   coordinates
   quality_control
   processing
   export_rename
   pipeline
   tutorial
