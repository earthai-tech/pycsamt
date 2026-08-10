.. _getting-started-first-survey:

Inspect a first survey
======================

A first survey pass should answer a limited set of questions before any
correction or inversion begins: Did the intended files load? Are station
identities and order plausible? Do neighboring stations have comparable
frequency coverage? Which optional transfer-function sections are absent? Do
the response diagnostics change systematically along the line?

This page uses the bundled WILLY L18PLT AMT line. Every console result and
figure below was generated from the 28 EDI files in
``data/AMT/WILLY_DATA/L18PLT``. If your input is AVG, J-format, spectral,
time-series, or TEM/TDEM data, first identify and convert it as described in
:doc:`data_formats`.

Load the survey
---------------

Use :func:`pycsamt.api.read_edis` for the public survey view. The example keeps
the search inside one known line, requires at least one readable input, and
suppresses progress output so that the captured result is stable:

.. code-block:: pycon

   >>> from pycsamt.api import read_edis
   >>> survey = read_edis(
   ...     "data/AMT/WILLY_DATA/L18PLT",
   ...     recursive=False,
   ...     strict=True,
   ...     on_dup="replace",
   ...     progress=False,
   ... )
   >>> print(survey)
   APISurvey: edi_survey
   sites: 28
   stations: 23-18-001A, 23-18-002U, 23-18-003A, 23-18-004A, 23-18-005U, 23-18-006A, 23-18-007U, 23-18-008U, ...
   source: data/AMT/WILLY_DATA/L18PLT

The count must be compared with the field manifest; ``28`` is correct only for
this bundled line. ``on_dup="replace"`` retains the later occurrence when the
EDI parser encounters a repeated station identity. It does not determine which
field file is scientifically authoritative, so investigate duplicates rather
than relying on discovery order in a controlled workflow.

Review source files and parser errors
-------------------------------------

Keep station identities tied to their source files while diagnosing a load:

.. code-block:: pycon

   >>> from pathlib import Path
   >>> survey.n_sites
   28
   >>> survey.stations[:5]
   ['23-18-001A', '23-18-002U', '23-18-003A', '23-18-004A', '23-18-005U']
   >>> [Path(path).name for path in survey.paths[:5]]
   ['18-001A.edi', '18-002U.edi', '18-003A.edi', '18-004A.edi', '18-005U.edi']
   >>> len(survey.errors())
   0

An empty error list means every discovered file used in this load was parsed;
it does not validate coordinates or transfer-function physics. With
``strict=False``, readable files can still be returned while parser failures
remain available through ``survey.errors()``. Always inspect that list before
treating a partial load as the complete survey.

Build the station inventory
---------------------------

``survey.summary()`` returns an ``APIFrame`` with one row per parsed EDI file.
Convert a copy to pandas for ordinary tabular inspection:

.. code-block:: pycon

   >>> summary = survey.summary()
   >>> print(summary)
   APIFrame: edi_survey_summary
   kind: edi.summary
   shape: 28 rows x 6 columns
   columns: station, path, n_freq, tipper, spectra, ts
   numeric: 1 columns
   missing: 0.0%
   source: data/AMT/WILLY_DATA/L18PLT
   >>> inventory = summary.to_pandas(copy=True)
   >>> inventory[["station", "n_freq", "tipper"]].head(3)
         station  n_freq  tipper
   0  23-18-001A      53   False
   1  23-18-002U      53   False
   2  23-18-003A      53   False

All 28 stations contain 53 impedance-frequency rows, while tipper, spectra,
and time-series sections are absent. Equal row counts are encouraging, but a
count alone cannot reveal whether stations sampled the same bands. Use the
inventory overview to inspect both facts together:

.. code-block:: pycon

   >>> from pycsamt.emtools import plot_survey_inventory_overview
   >>> labels = [name.replace("23-", "") for name in survey.stations]
   >>> fig = plot_survey_inventory_overview(
   ...     survey.collection,
   ...     station_labels=labels,
   ...     recursive=False,
   ...     title="L18PLT acquisition inventory and period coverage",
   ... )

.. figure:: /images/tutorials/read_edi_survey/survey_inventory_overview.png
   :align: center
   :width: 95%
   :alt: Frequency-row counts and sampled period coverage across the L18PLT stations.

   Inventory of the loaded L18PLT line. Station triangles and labels identify
   the shared columns. The upper markers report 53 rows per station, while the
   lower map confirms that those rows occupy the same period bands.

The uniform upper profile rules out isolated loss of complete frequency rows;
the uninterrupted lower map additionally rules out hidden station-specific
band gaps. The inventory table still shows that optional tipper, spectra, and
time-series sections are absent throughout this delivery. Methods requiring
those observations are therefore inappropriate unless they are obtained from
another source.

Inspect one station
-------------------

Before computing line-wide diagnostics, inspect one parsed EDI object:

.. code-block:: pycon

   >>> first = survey[0]
   >>> first.station, Path(first.path).name
   ('23-18-001A', '18-001A.edi')
   >>> len(first.Z.freq)
   53
   >>> survey.get_site("23-18-001A").station
   '23-18-001A'
   >>> survey.get_site("missing-station") is None
   True

The EDI header identity and source filename are related but not identical.
Keep both columns when auditing a survey, and do not reconstruct station
identity by stripping or adding a filename prefix unless the field metadata
documents that convention.

Plot a survey fingerprint
-------------------------

A compact first diagnostic compares phase-sensitive quantities across station
and period. The following public function produces the captured figure:

.. code-block:: pycon

   >>> from pycsamt.emtools import plot_survey_fingerprint
   >>> fig = plot_survey_fingerprint(
   ...     survey.collection,
   ...     quantities=["skew", "ellipt", "s1"],
   ...     render="imshow",
   ...     plot_kws={"interpolation": "bilinear"},
   ...     station_grid=True,
   ...     period_range=(1e-4, 1.0),
   ...     recursive=False,
   ...     title="L18PLT quick survey fingerprint",
   ...     figsize=(11.2, 7.6),
   ... )
   >>> len(fig.axes) >= 3
   True

The public API defaults to three panels—``skew``, ``ellipt``, and ``s1``
(:math:`\phi_{\max}`)—rendered with ``pcolormesh``. This Getting Started view
requests ``imshow`` with bilinear interpolation for a visibly continuous
raster presentation. Pass a
single name such as ``quantities="skew"`` for one panel, or any ordered subset
such as ``["ellipticity", "phi_max"]`` for two. ``cmaps`` accepts either one
colormap or a quantity-to-colormap mapping, while ``plot_kws`` and
``quantity_kws`` forward global and per-panel options to the selected
Matplotlib renderer. Contours resolve through the package-wide
:data:`pycsamt.api.PYCSAMT_CONTOUR` review style by default. Pass
``contours=False`` to hide them, or use ``contour_kws`` for a one-call
override. The complete configuration pattern is documented in
:doc:`../api_guide/contour`.

The dotted vertical guides pass through the same station centre in all three
panels. Set ``station_grid=False`` to hide them, or pass
``station_grid_kws`` to control their colour, width, line style, opacity, and
drawing order.

Bilinear interpolation smooths the displayed pixels only; it does not add
measured stations or periods. Use the default ``pcolormesh`` view when exact
cell boundaries are important for diagnosis.

.. figure:: /images/tutorials/read_edi_survey/survey_fingerprint.png
   :align: center
   :width: 100%
   :alt: Skew, ellipticity, and maximum phase-tensor value across L18PLT.

   Phase-sensitive fingerprint of the L18PLT survey over periods from
   :math:`10^{-4}` to :math:`1` s. Columns are stations and rows are period
   samples; color represents the quantity shown on each panel.

The left part of the line is comparatively subdued at short periods, whereas
stations in the latter half show strong, coherent changes in skew sign,
ellipticity, and the maximum phase-tensor value. Because several quantities
change in the same station-period region, the pattern deserves targeted
phase-tensor, strike, and source-effect review. It should not be labeled a bad
station or a geological boundary from this quick-look figure alone: distortion,
three-dimensional structure, noise, and acquisition effects can produce
overlapping signatures.

The deterministic generator for the inventory and fingerprint figures is
available below. It uses public pyCSAMT interfaces and writes the images into
the documentation image tree; user workflows should call those public APIs
directly rather than importing this documentation script.

.. code-dropdown:: /../../scripts/generate_tutorial_read_edi.py
   :language: python
   :linenos:
   :title: View the figure-generation source

Decide whether to continue
--------------------------

Before processing this or another survey, verify:

* the loaded count matches the expected field stations;
* every parser error is understood;
* station identities map unambiguously to source files;
* coordinates and elevation have the expected datum and units;
* frequency values—not only their counts—are compatible across stations;
* required impedance or tipper components are present;
* unusual station-period regions have been reviewed with more specific
  diagnostics.

.. warning::

   A complete inventory and an attractive diagnostic plot do not establish
   inversion readiness. Static shift, dead bands, source effects, coordinate
   errors, uncertainty estimates, and dimensionality assumptions require
   separate checks.

Continue to :doc:`../tutorials/read_edi_survey` for deeper EDI inspection or
:doc:`../tutorials/inspect_and_qc_survey` for quality-control diagnostics. The
reusable loading and normalization rules are summarized in
:doc:`../user_guide/data_loading`; processing methods are organized under
:doc:`../user_guide/emtools/index`.
