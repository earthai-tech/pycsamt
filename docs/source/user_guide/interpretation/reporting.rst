.. _interpretation_reporting:

Export and reporting
====================

Interpretation reporting turns analysis into a reviewable and reproducible
project record. A useful report does more than display a polished resistivity
section: it identifies the source data, processing and inversion decisions,
interpretation assumptions, calibration evidence, uncertainty, limitations,
and the exact products being delivered.

pyCSAMT provides focused exporters and plotting classes rather than a single
``make_report()`` function for interpretation. This separation is deliberate.
CSV, XYZ, LAS, VTK, figures, configuration files, and narrative conclusions
serve different audiences and must be assembled according to the project's
review and governance requirements.

.. admonition:: Reporting principle
   :class: important

   A deliverable is not reproducible merely because it can be reopened. A
   reviewer must be able to identify its inputs, units, coordinate convention,
   method, assumptions, software version, uncertainty, and approval state.

Reporting objectives
--------------------

A complete interpretation package should allow another qualified reviewer to:

* understand the survey and decision being supported;
* trace every final product back to a reviewed inversion model;
* distinguish measured, inverted, calibrated, and interpreted quantities;
* reproduce the principal calculations and figures;
* inspect calibration residuals and withheld validation results;
* identify unsupported areas and alternative explanations;
* know which files are preliminary, superseded, approved, or authoritative;
* reuse machine-readable outputs without guessing units or coordinates.

Recommended reporting workflow
------------------------------

#. define the audience, decision, and approval level;
#. freeze the reviewed source model and interpretation configuration;
#. assign stable identifiers and output status;
#. generate machine-readable tables and grids;
#. generate figures using consistent scales and labels;
#. write the methods, results, uncertainty, and limitations narrative;
#. validate exported content independently of the in-memory objects;
#. assemble a manifest and checksums;
#. complete scientific and technical review;
#. publish an immutable approved package while preserving working files.

1. Define audience and reporting level
--------------------------------------

Different audiences need different products:

Scientific reviewer
   Requires data quality, inversion diagnostics, parameter rationale,
   calibration residuals, uncertainty, and alternative interpretations.

Project decision-maker
   Requires the decision question, principal findings, confidence, risks,
   limitations, and recommended next actions in plain language.

GIS or modeling specialist
   Requires coordinate reference, geometry, units, null handling, field
   definitions, and machine-readable files.

Field team
   Requires station names, profile direction, target coordinates, depth
   reference, uncertainty, access constraints, and unambiguous maps.

Regulator or client
   May require named standards, approval signatures, data lineage, controlled
   revisions, and explicit statements of professional responsibility.

Define the reporting level before export:

``working``
   Internal exploratory output. It may change and must not be used for a final
   decision.

``review``
   Frozen candidate package submitted for scientific or technical review.

``approved``
   Versioned package that has passed the project's acceptance procedure.

``superseded``
   Previously issued package retained for audit but no longer authoritative.

Put the status in the report, manifest, and directory name—not only in an
email or surrounding conversation.

2. Separate evidence classes
----------------------------

Every table and figure should make clear which kind of quantity it contains:

Measured
   Field or laboratory observations such as impedance, water level, EC,
   lithology, pumping-test transmissivity, or slug-test conductivity.

Processed
   Corrected or derived observations such as filtered impedance, static-shift
   corrected responses, or QC metrics.

Inverted
   Model properties estimated by fitting the geophysical observations, such as
   the calculated resistivity model (CRM).

Calibrated
   A model modified or parameterized using borehole or field constraints, such
   as the calibrated new model (NM).

Interpreted
   Geological or hydrogeological labels and derived properties based on
   assumptions and evidence.

Predicted
   Values calculated at validation locations or under scenarios.

Do not combine these classes in one column named ``value`` without an origin
field. A geological contact inferred from resistivity is not a measured
borehole contact, even when they agree.

3. Freeze provenance before export
----------------------------------

Assign stable identifiers to the source survey, processing run, inversion run,
interpretation run, and reporting package. A simple naming pattern is:

.. code-block:: text

   <project>_<line>_<stage>_<YYYYMMDD>_<revision>

For example:

.. code-block:: text

   willy_L18_interpretation_20260712_r01

Before generating final files, record:

* project and survey-line identifiers;
* input file inventory or upstream manifest;
* processing and inversion run identifiers;
* source-model method, RMS, mesh, and reliable depth range;
* coordinate reference system, profile origin, azimuth, and vertical datum;
* pyCSAMT and Python versions;
* interpretation configuration and rock-database version;
* boreholes and constraints used for calibration;
* observations reserved for validation;
* uncertainty bounds, sample count, seed, and failure rate;
* author, reviewer, organization, timestamp, and status.

The source model and configuration should remain unchanged while the review
package is generated. If either changes, issue a new interpretation run rather
than silently replacing files inside the existing package.

4. Use a controlled directory structure
---------------------------------------

A practical interpretation package can use:

.. code-block:: text

   willy_L18_interpretation_20260712_r01/
   ├── README.md
   ├── manifest.yml
   ├── CHANGELOG.md
   ├── source/
   │   ├── inversion_manifest.yml
   │   ├── model_snapshot.npz
   │   └── residual_summary.csv
   ├── configuration/
   │   ├── interpretation.yml
   │   ├── petrophysics.yml
   │   └── rock_database.csv
   ├── evidence/
   │   ├── borehole_inventory.csv
   │   ├── constraints.csv
   │   └── validation_observations.csv
   ├── tables/
   │   ├── stratigraphic_logs.csv
   │   ├── hydro_cells.csv
   │   ├── hydro_by_station.csv
   │   ├── uncertainty_by_station.csv
   │   └── validation_residuals.csv
   ├── grids/
   │   └── calibrated_resistivity.vtk
   ├── logs/
   │   └── S017.las
   ├── gis/
   │   └── profile.xyz
   ├── figures/
   │   ├── crm_nm_misfit.png
   │   ├── hydraulic_K.png
   │   └── uncertainty_profile.png
   ├── report/
   │   └── technical_report.pdf
   └── checksums.sha256

The exact structure may follow organizational standards. The important rule is
to separate immutable sources, configuration, evidence, machine-readable
outputs, visual products, and narrative reports.

5. Export stratigraphic logs to CSV
----------------------------------

:func:`pycsamt.interp.export.to_csv` writes all
:class:`pycsamt.interp.StratigraphicLog` objects to a flat table:

.. code-block:: python

   from pathlib import Path
   from pycsamt.interp import export

   root = Path("willy_L18_interpretation_20260712_r01")
   table_path = export.to_csv(
       logs,
       root / "tables" / "stratigraphic_logs.csv",
   )
   print(table_path)

The output fields are:

``station``
   Station identifier.

``x_m``
   Along-profile position in metres.

``depth_m``
   Depth below the model surface in metres.

``rho_log10``
   :math:`\log_{10}` resistivity where linear resistivity is in ohm metres.

``rho_ohm_m``
   Linear resistivity in ohm metres.

``lithology``
   Interpreted lithology label assigned to the depth cell.

The current exporter writes both resistivity columns. Its ``log_rho`` argument
is retained in the signature, but does not remove either column. Consumers
should select the explicitly named field rather than infer units from values.

Validate the CSV after writing:

.. code-block:: python

   import csv

   with table_path.open(newline="", encoding="utf-8") as stream:
       reader = csv.DictReader(stream)
       required = {
           "station", "x_m", "depth_m",
           "rho_log10", "rho_ohm_m", "lithology",
       }
       missing = required.difference(reader.fieldnames or [])
       if missing:
           raise ValueError(f"Missing CSV fields: {sorted(missing)}")
       rows = list(reader)

   if not rows:
       raise ValueError("The stratigraphic export is empty.")

Never use a spreadsheet's automatic formatting as the authoritative copy.
Station names may be converted to dates or numbers, and scientific notation or
decimal separators can change across locales.

6. Export Oasis Montaj XYZ
--------------------------

:func:`pycsamt.interp.export.to_oasis_montaj_xyz` writes each station log as a
``/ Line`` block:

.. code-block:: python

   import numpy as np

   elevation_m = np.array([238.4, 239.1, 240.0, 239.6])

   xyz_path = export.to_oasis_montaj_xyz(
       logs,
       root / "gis" / "profile.xyz",
       y=0.0,
       elevation=elevation_m,
       log_rho=False,
   )

Coordinate behavior must be documented:

* X is ``log.station_x``, normally an along-profile distance, not necessarily
  an easting;
* the ``y`` argument is one scalar assigned to every point;
* without ``elevation``, Z is negative depth;
* with ``elevation``, Z is surface elevation minus depth;
* the elevation array must correspond to the log order;
* ``log_rho=True`` writes log10 resistivity; ``False`` writes linear ohm-m
  resistivity;
* spaces in lithology labels are converted to underscores.

This exporter does not attach a coordinate reference system or vertical datum.
Include those in the manifest and, where possible, in a companion metadata
file. Do not label profile distance as easting unless it has actually been
converted into the project CRS.

Custom channel names can be supplied:

.. code-block:: python

   export.to_oasis_montaj_xyz(
       logs,
       root / "gis" / "profile_linear_rho.xyz",
       log_rho=False,
       channels=["PROFILE_X_M", "PROFILE_Y_M", "ELEV_M",
                 "RHO_OHM_M", "LITHOLOGY"],
   )

The custom header changes labels, not coordinate transformation or data
semantics.

7. Export individual LAS logs
-----------------------------

:func:`pycsamt.interp.export.to_las` writes one stratigraphic log as LAS 2.0:

.. code-block:: python

   las_path = export.to_las(
       logs[0],
       root / "logs" / "S017.las",
       well_name="S017",
       company="Example Hydrogeophysics Project",
       null_value=-9999.25,
       log_rho=False,
   )

The depth curve is in metres. ``log_rho=True`` writes log10 resistivity;
``False`` writes linear resistivity. Review the LAS header and curves in the
receiving application before delivery.

.. warning::

   The current LAS exporter encodes lithology as ``hash(lithology) % 1000``.
   Python hash randomization means these integer codes are not guaranteed to
   remain stable across processes. Do not treat them as a durable corporate
   lithology dictionary. Deliver an explicit station/depth/lithology CSV and a
   controlled code table when stable codes are required.

LAS output is a station-column interpretation, not a drilled well log unless
the station actually represents a borehole and the interpretation has been
validated accordingly. State ``EM-derived`` prominently in the report and
curve description.

8. Export calibrated models to VTK
----------------------------------

:func:`pycsamt.interp.export.to_vtk` writes a
:class:`pycsamt.interp.ResistivityModel` as an ASCII rectilinear grid:

.. code-block:: python

   vtk_path = export.to_vtk(
       calibrated_model,
       root / "grids" / "calibrated_resistivity.vtk",
       log_rho=False,
       field_name="rho_ohm_m",
   )

Important format details:

* X coordinates are model ``x_centers``;
* model depths are written as the VTK Y coordinates;
* the VTK Z dimension contains one coordinate at zero;
* values are written as point data;
* missing resistivity is written as ``-9999.0``;
* no CRS, vertical datum, topography, or interpretation confidence is embedded;
* the function exports resistivity only, not porosity, saturation, K, or
  lithological labels.

Use ``field_name`` to encode units explicitly. A file named
``calibrated_model.vtk`` is still ambiguous unless the manifest states whether
the field is CRM or NM, log10 or linear, and which depth and coordinate
conventions apply.

9. Export deterministic hydro results
-------------------------------------

An :class:`pycsamt.interp.EMHydroResult` provides two CSV levels:

.. code-block:: python

   hydro_result.to_csv(root / "tables" / "hydro_cells.csv")
   hydro_result.station_report_csv(
       root / "tables" / "hydro_by_station.csv"
   )

The cell-level file includes station, profile position, depth, log10 and linear
resistivity, porosity, saturation, and hydraulic conductivity. The station
summary includes water-table depth, saturated-zone porosity and K summaries,
transmissivity, storativity, Dar–Zarrouk parameters, and the TDS indicator.

Document these interpretation qualifications:

* water table is threshold-derived and may be ``nan``;
* K is petrophysically derived, not measured;
* transmissivity integrates the represented saturated model interval;
* unconfined storativity is approximated from porosity;
* TDS is based on configured scalar pore-water resistivity;
* columns with failed water-table detection require special review, as
  explained in :doc:`hydrogeophysics`.

If pandas is installed, ``hydro_result.to_dataframe()`` supports further
review, but any derived table should retain the original field names and units.

10. Export qualitative hydro interpretation
--------------------------------------------

A :class:`pycsamt.interp.HydroGeophysicalModel` can write cell classifications
and interpreted zones:

.. code-block:: python

   qualitative_model.to_csv(
       root / "tables" / "hydro_units.csv"
   )
   qualitative_model.zones_to_csv(
       root / "tables" / "aquifer_zones.csv"
   )

The cell table contains hydro-unit labels and confidence values. The zone table
contains station, position, top, bottom, thickness, mean resistivity,
confidence, and zone type.

Report the rule set, thresholds, context, rock database, and evidence used to
produce these categories. A numerical confidence emitted by a rule-based
classifier is not automatically a calibrated probability.

11. Export uncertainty summaries
--------------------------------

:class:`pycsamt.interp.UncertaintyResult` writes a per-station summary:

.. code-block:: python

   uncertainty.to_csv(
       root / "tables" / "uncertainty_by_station.csv"
   )

The table includes water-table mean, standard deviation, P10, P90, P90–P10
range, detection percentage, and transmissivity summaries. Preserve the
corresponding bounds, distribution type, free parameter order, sample count,
seed, and failure diagnostics in configuration or manifest files.

The raw water-table and transmissivity ensembles returned by
``MonteCarloHydro.run_ensemble()`` are not written by this CSV method. Archive
them separately in an appropriate binary array format when empirical
probabilities or distribution plots must be reproduced.

State that the intervals are conditional on the sampled parameters and fixed
source resistivity model. See :doc:`uncertainty` for uncertainty sources that
remain outside this ensemble.

12. Report calibration residuals
--------------------------------

When quantitative field constraints are used, retain per-constraint residuals:

.. code-block:: python

   residuals = calibrator.constraint_residuals(calibrated_result)

The method returns dictionaries rather than writing a file. Save them using a
transparent table writer:

.. code-block:: python

   import csv

   residual_path = root / "tables" / "calibration_residuals.csv"
   residual_path.parent.mkdir(parents=True, exist_ok=True)

   if residuals:
       fields = sorted({key for row in residuals for key in row})
       with residual_path.open("w", newline="", encoding="utf-8") as stream:
           writer = csv.DictWriter(stream, fieldnames=fields)
           writer.writeheader()
           writer.writerows(residuals)

Report individual residuals, not only the optimizer's total objective. Clearly
separate calibration observations from withheld validation observations.

13. Generate figures consistently
---------------------------------

Interpretation figures are review evidence and communication products. Useful
classes include:

* :class:`pycsamt.interp.plot.PlotCalibratedModel`;
* :class:`pycsamt.interp.plot.PlotStratigraphicLog`;
* :class:`pycsamt.interp.plot.PlotFenceDiagram`;
* :class:`pycsamt.interp.plot.PlotHydroSection`;
* :class:`pycsamt.interp.plot.PlotWaterTableProfile`;
* :class:`pycsamt.interp.plot.PlotPetrophysicalCrossPlot`;
* :class:`pycsamt.interp.plot.PlotAquiferCharacterization`;
* :class:`pycsamt.interp.plot.PlotUncertaintySection`;
* :class:`pycsamt.interp.plot.PlotUncertaintyProfile`;
* :class:`pycsamt.interp.plot.PlotUncertaintyHistogram`.

For example:

.. code-block:: python

   from pycsamt.interp import plot as iplot

   figures = root / "figures"
   figures.mkdir(parents=True, exist_ok=True)

   fig = iplot.PlotCalibratedModel(
       crm,
       calibrated_model,
       calibrator.misfit_map(),
       vmin_rho=1.0,
       vmax_rho=4.5,
   ).plot()
   fig.savefig(figures / "crm_nm_misfit.png", dpi=300,
               bbox_inches="tight")

   fig = iplot.PlotHydroSection(
       hydro_result,
       quantity="K",
       vmin=-10.0,
       vmax=-3.0,
       depth_max=200.0,
   ).plot()
   fig.savefig(figures / "hydraulic_K.png", dpi=300,
               bbox_inches="tight")

Figure rules
~~~~~~~~~~~~

Every figure should identify:

* project, line, method, and model status;
* profile direction and horizontal coordinate;
* depth or elevation reference and unit;
* property and unit, including log transformation;
* station locations where relevant;
* consistent color scale across compared scenarios;
* missing or masked cells;
* water-table detection gaps;
* uncertainty representation;
* run or figure identifier.

Do not use different automatic color limits to compare scenarios. The same
structure can appear stronger or weaker solely because the color normalization
changed. Avoid rainbow palettes where they obscure ordering or accessibility,
and check grayscale and color-vision readability when required by the project.

14. Write the technical narrative
---------------------------------

A concise but complete interpretation report normally contains:

Executive summary
   Decision, principal findings, confidence, limitations, and recommended next
   action. Avoid unexplained software or inversion terminology.

Objectives and scope
   Survey area, question, target depth, methods, exclusions, and reporting
   status.

Data and processing
   Acquisition inventory, data quality, exclusions, corrections, coordinate
   handling, and unresolved artifacts.

Inversion
   Backend, dimensionality, mesh, errors, regularization, convergence,
   residuals, sensitivity, model scenarios, and reliable interpretation depth.

Interpretation method
   CRM normalization, boreholes, rock database, calibration tolerance,
   classification logic, hydrogeophysical equations, and configuration.

Results
   Observed patterns and derived quantities stated separately from geological
   hypotheses.

Calibration and validation
   Evidence roles, residuals, withheld results, scale compatibility, matches,
   and mismatches.

Uncertainty
   Data, inversion, petrophysical, calibration, and interpretive uncertainty;
   intervals and detection rates; assumptions not propagated.

Conclusions and recommendations
   Answers to the stated questions, confidence-qualified targets, rejected
   alternatives, and specific follow-up measurements.

Limitations
   Scientific, spatial, computational, and operational constraints that affect
   use of the deliverables.

Appendices
   Configuration tables, file manifest, symbols and units, residuals,
   additional scenarios, and reviewer record.

Keep observation and inference linguistically separate. For example:

* observation: "A conductive zone occurs between profile distances 600 and
  900 m below approximately 40 m depth.";
* interpretation: "The zone is consistent with saturated weathered material,
  but clay-rich material remains a plausible alternative.";
* validation: "BH03 intersects weathered granite in this interval; BH03 was
  withheld from calibration.";
* uncertainty: "The boundary varies from 35 to 58 m across accepted scenarios."

15. Build a machine-readable manifest
-------------------------------------

The manifest is the package's index. YAML or JSON is suitable. A minimal YAML
structure might be:

.. code-block:: yaml

   schema_version: 1
   package_id: willy_L18_interpretation_20260712_r01
   status: review
   project: willy
   survey_line: L18
   created_utc: 2026-07-12T12:00:00Z
   software:
     package: pycsamt
     version: 2.0.0
   coordinates:
     horizontal_reference: profile_distance
     horizontal_unit: m
     vertical_reference: depth_below_surface
     vertical_positive: down
     vertical_unit: m
   resistivity:
     linear_unit: ohm_m
     model_storage: log10_ohm_m
   source_runs:
     processing: processing_run_id
     inversion: inversion_run_id
   interpretation:
     run_id: interpretation_run_id
     calibration_boreholes: [BH01, BH02]
     validation_boreholes: [BH03]
   uncertainty:
     n_samples: 500
     seed: 42
     interval: P10_P90
   files:
     - path: tables/stratigraphic_logs.csv
       role: interpreted_station_cells
     - path: figures/crm_nm_misfit.png
       role: calibration_review

Extend this with checksums, sizes, media types, CRS identifiers, configuration
hashes, and approval metadata according to project requirements. Do not put
secrets, personal data, or machine-specific absolute paths in a deliverable
manifest.

16. Add checksums and validate files
-----------------------------------

Checksums detect accidental change after approval. Generate them with an
organizationally approved tool and store paths relative to the package root.
Checksums prove file integrity, not scientific correctness.

Validation should include:

Structural checks
   Required files exist, are non-empty, and match manifest entries.

Schema checks
   CSV headers, units, data types, null conventions, and unique identifiers are
   correct.

Numerical checks
   Exported ranges and row counts agree with in-memory results; log and linear
   resistivity correspond; P10 ≤ P50 ≤ P90 where finite.

Coordinate checks
   Profile distance, station order, elevation/depth convention, CRS, and datum
   agree across figures and files.

Visual checks
   Figures open, labels are legible, color scales match comparisons, masked
   regions are visible, and no plotting layer was clipped.

Round-trip checks
   Open each format in at least one target consumer when practical. Verify LAS
   curves, XYZ columns, VTK orientation, and CSV encoding.

Scientific checks
   Conclusions match the approved tables and figures, calibration/validation
   roles are correct, and limitations are not omitted.

17. Review and approval
-----------------------

Use separate review roles when project scale permits:

Scientific review
   Tests geophysical, geological, hydrogeological, and uncertainty reasoning.

Technical review
   Tests code paths, units, file schemas, reproducibility, and internal
   consistency.

Editorial review
   Tests clarity, terminology, captions, accessibility, and audience fit.

Approval
   Confirms that the package meets the project's governance and release
   requirements.

Track comments and dispositions. If a review changes inputs, parameters, or
conclusions, increment the revision and regenerate dependent outputs. Do not
edit an approved binary or CSV in place.

18. Handle revisions and superseded products
--------------------------------------------

Maintain a changelog with:

* revision identifier and date;
* author and approver;
* files added, removed, or replaced;
* scientific reason for change;
* impact on conclusions and downstream users;
* identifier of the superseded package.

Never reuse the same approved package identifier for different content.
Preserve superseded packages in read-only archival storage with an obvious
status marker. Notify downstream users when a revision changes target
locations, depths, confidence, or safety-relevant conclusions.

19. Protect sensitive information
---------------------------------

Interpretation packages can contain private well locations, infrastructure,
landowner details, water-quality information, or commercially sensitive
targets. Before release:

* classify each file according to project policy;
* remove unnecessary personal or machine-specific information;
* limit coordinate precision when authorized and appropriate;
* avoid embedding credentials, API keys, or local absolute paths;
* verify image metadata and document properties;
* separate public summaries from controlled technical appendices.

Do not reduce coordinate precision in a scientific archive unless the precise
authoritative coordinates are preserved in an appropriately controlled source.

Complete export example
-----------------------

The following example assembles core machine-readable products and figures.
It assumes ``logs``, ``crm``, ``calibrated_model``, ``calibrator``,
``hydro_result``, and ``uncertainty`` were created and reviewed as described in
the preceding guides:

.. code-block:: python

   from pathlib import Path
   import csv

   from pycsamt.interp import export, plot as iplot

   root = Path("willy_L18_interpretation_20260712_r01")
   tables = root / "tables"
   figures = root / "figures"
   grids = root / "grids"
   logs_dir = root / "logs"
   gis = root / "gis"

   for directory in (tables, figures, grids, logs_dir, gis):
       directory.mkdir(parents=True, exist_ok=True)

   # Geological interpretation products.
   export.to_csv(logs, tables / "stratigraphic_logs.csv")
   export.to_oasis_montaj_xyz(
       logs, gis / "profile.xyz", log_rho=False
   )
   export.to_las(
       logs[0], logs_dir / f"{logs[0].station_name}.las",
       log_rho=False,
   )
   export.to_vtk(
       calibrated_model,
       grids / "calibrated_resistivity.vtk",
       log_rho=False,
       field_name="rho_ohm_m",
   )

   # Hydrogeophysical and uncertainty products.
   hydro_result.to_csv(tables / "hydro_cells.csv")
   hydro_result.station_report_csv(tables / "hydro_by_station.csv")
   uncertainty.to_csv(tables / "uncertainty_by_station.csv")

   # Calibration residuals.
   residuals = calibrator.constraint_residuals(hydro_result)
   if residuals:
       fields = sorted({key for row in residuals for key in row})
       with (tables / "calibration_residuals.csv").open(
           "w", newline="", encoding="utf-8"
       ) as stream:
           writer = csv.DictWriter(stream, fieldnames=fields)
           writer.writeheader()
           writer.writerows(residuals)

   # Review figures with explicit comparison scales.
   fig = iplot.PlotCalibratedModel(
       crm,
       calibrated_model,
       calibrator.misfit_map(),
       vmin_rho=1.0,
       vmax_rho=4.5,
   ).plot()
   fig.savefig(figures / "crm_nm_misfit.png", dpi=300,
               bbox_inches="tight")

   fig = iplot.PlotUncertaintyProfile(uncertainty).plot()
   fig.savefig(figures / "uncertainty_profile.png", dpi=300,
               bbox_inches="tight")

This script creates files; it does not by itself create the narrative report,
manifest, checksums, review record, or approval. Those are required parts of a
controlled reporting workflow.

Delivery checklist
------------------

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Check
     - Acceptance evidence
   * - Scope and status are explicit
     - Audience, decision, project, line, revision, and working/review/approved
       status.
   * - Inputs are traceable
     - Survey, processing, inversion, interpretation, and configuration IDs.
   * - Quantity classes are separated
     - Measured, inverted, calibrated, interpreted, and predicted fields are
       labeled.
   * - Coordinates are unambiguous
     - CRS or profile reference, origin, direction, units, datum, and vertical
       sign convention.
   * - Units are explicit
     - Linear/log resistivity, depth/elevation, K, T, storativity, EC, and TDS.
   * - Exports are validated
     - Schema, row count, numerical range, null handling, target-application
       round trip, and manifest match.
   * - Figures are comparable
     - Stable color limits, visible missing data, readable labels, and figure
       identifiers.
   * - Calibration is reviewable
     - Constraints, fitted parameters, bounds, restarts, and residuals.
   * - Validation is independent
     - Withheld evidence, prediction intervals, mismatches, and dispositions.
   * - Uncertainty is conditional and complete
     - Bounds, samples, seed, failures, detection rates, scenarios, and omitted
       uncertainty sources.
   * - Limitations and alternatives are stated
     - Reliable depth, resolution, conceptual ambiguity, and usage limits.
   * - Package integrity is controlled
     - Manifest, checksums, changelog, reviewer, approver, and immutable copy.

Common reporting mistakes
-------------------------

Avoid these errors:

* delivering only a color image without source values or provenance;
* labeling along-profile distance as easting;
* mixing depth below surface with elevation above datum;
* omitting whether resistivity is log10 or linear;
* presenting calibrated lithology as directly observed geology;
* treating LAS hash-based lithology codes as stable identifiers;
* calling derived hydraulic conductivity a measured field value;
* hiding water-table non-detections or Monte Carlo failures;
* reporting P10–P90 as the complete range of possible outcomes;
* changing color scales between scenarios;
* using calibration wells again as independent validation;
* manually editing exported tables without recording the change;
* replacing an approved package without incrementing its revision;
* assuming checksums establish scientific correctness;
* omitting known alternative explanations because they complicate the summary.

Next steps
----------

Use this page with:

* :doc:`workflow` for geological interpretation and calibration;
* :doc:`hydrogeophysics` for deterministic hydrogeophysical products;
* :doc:`uncertainty` for conditional intervals and validation;
* :doc:`../map/index` for spatial context and mapping exports;
* :doc:`../inversion/index` for source-model provenance and diagnostics.
