.. _getting-started-data-formats:

Data Formats
============

pyCSAMT v2 works with several electromagnetic field-data and workflow formats.
The safest starting point is simple:

* use EDI when you already have frequency-domain MT/AMT/CSAMT transfer
  functions;
* use Zonge AVG/AMTAVG readers when you are working directly from Zonge field
  exports;
* use J-format readers when your survey is stored in Jones-style transfer
  function files;
* use TDEM readers and converters when you need to transform time-domain EM
  soundings into frequency-domain objects;
* use transformers when the goal is to convert AVG, J, or spectra-like data
  into EDI-compatible products;
* use site tables when the task is station selection, survey metadata,
  coordinates, or profile editing.

This page gives the practical map.  It does not replace the API reference, but
it tells you which package to reach for first.


Format overview
---------------

.. list-table::
   :header-rows: 1
   :widths: 18 25 28 29

   * - Format
     - Typical file/content
     - Main package
     - First user-facing entry point
   * - EDI
     - ``.edi`` transfer-function files for MT, AMT, CSAMT, or EMAP-style
       workflows.
     - ``pycsamt.api`` and ``pycsamt.seg``
     - :func:`pycsamt.api.read_edis`
   * - SEG-style EDI sections
     - EDI headers, spectra sections, time-series sections, and survey
       sections.
     - ``pycsamt.seg``
     - :class:`pycsamt.seg.EDIFile`
   * - Zonge AVG / AMTAVG
     - Zonge apparent resistivity, phase, frequency, station, and QC exports.
     - ``pycsamt.zonge``
     - :func:`pycsamt.zonge.load_avg`
   * - J-format
     - Jones-style transfer-function files.
     - ``pycsamt.jones``
     - :class:`pycsamt.jones.JFile`
   * - TDEM / TEM
     - Time gates, decay curves, TEM AVG exports, logs, coordinates, and
       waveform metadata.
     - ``pycsamt.tdem``
     - :class:`pycsamt.tdem.TEMReader` or
       :func:`pycsamt.tdem.read_temavg_survey`
   * - Spectra
     - Spectral estimates that can be transformed into EDI-compatible objects.
     - ``pycsamt.seg`` and ``pycsamt.transformers``
     - :class:`pycsamt.transformers.SpectraToEDI`
   * - Site and station tables
     - Coordinates, station selection files, profile definitions, reports, and
       export tables.
     - ``pycsamt.site`` and ``pycsamt.metadata``
     - ``pycsamt.site`` helpers and :func:`pycsamt.api.read_sites`
   * - Inversion files
     - Occam2D, ModEM, mesh, model, data, startup, response, and result files.
     - ``pycsamt.models`` and ``pycsamt.inversion``
     - Engine-specific builders such as ``pycsamt.models.occam2d``


Recommended first path
----------------------

If you have a directory of EDI files, start with the public API:

.. code-block:: python
   :linenos:

   from pycsamt.api import read_edis

   survey = read_edis(
       "data/willy/edis",
       recursive=True,
       strict=False,
       on_dup="replace",
       progress="auto",
   )

   print(survey)
   print(survey.summary())

This returns an ``APISurvey`` view around the underlying EDI collection.  It is
the easiest object to pass into tutorials, pipeline workflows, agents, and
interactive inspection.

Use ``strict=True`` when validating data before a production run:

.. code-block:: python
   :linenos:

   survey = read_edis("data/willy/edis", strict=True)


EDI files
---------

EDI is the main interchange format for frequency-domain MT, AMT, and CSAMT
workflows.  It usually contains station metadata, frequencies or periods,
impedance tensor components, tipper components, errors, and survey headers.

Use the high-level API for batches:

.. code-block:: python
   :linenos:

   from pycsamt.api import read_edis

   survey = read_edis("data/edis/**/*.edi", recursive=True)

Read one file with the lower-level SEG parser:

.. code-block:: python
   :linenos:

   from pycsamt.seg import EDIFile

   edi = EDIFile("data/edis/S001.edi")
   print(edi.station)

Validate EDI-like files:

.. code-block:: python
   :linenos:

   from pycsamt.seg import IsEdi

   validator = IsEdi("data/edis/S001.edi")

CLI equivalents:

.. code-block:: bash
   :linenos:

   pycsamt edi info data/edis
   pycsamt edi stations data/edis
   pycsamt edi validate data/edis

Use EDI when:

* your data are already transfer functions;
* you want to run processing, QC, static-shift correction, tensor analysis, or
  inversion preparation;
* you need a portable format accepted by many geophysical tools.


SEG-style parser package
------------------------

``pycsamt.seg`` is the low-level package for EDI-like sections and related
survey structures.  It exposes objects such as:

* ``EDIFile`` for a parsed EDI file;
* ``Spectra`` and ``SpectraSECT`` for spectra sections;
* ``TimeSeries`` and ``TSect`` for time-series sections;
* ``EDIProfile``, ``Stations``, and ``Topography`` for survey geometry;
* ``build_dataset`` for xarray-style data assembly where supported.

Most users should start with :func:`pycsamt.api.read_edis`.  Use
``pycsamt.seg`` directly when you are building parsers, inspecting individual
sections, debugging file content, or assembling lower-level datasets.


Zonge AVG and AMTAVG
--------------------

Zonge AVG/AMTAVG files are common field-processing outputs containing station,
frequency, apparent resistivity, phase, component, QC, and survey metadata.

Read an AVG file:

.. code-block:: python
   :linenos:

   from pycsamt.zonge import load_avg

   avg = load_avg("data/zonge/line01.avg")
   print(avg)

Object-oriented entry points:

.. code-block:: python
   :linenos:

   from pycsamt.zonge import AVG, AMTAVG

   avg = AVG("data/zonge/line01.avg")
   amt = AMTAVG("data/zonge/line01.amtavg")

Write AVG-compatible output where supported:

.. code-block:: python
   :linenos:

   from pycsamt.zonge import write_avg

   write_avg(avg, "results/line01_clean.avg")

Zonge processing helpers include adaptive moving averages, tensor utilities,
phase/resistivity objects, station metadata, and QC alias maps.

Use Zonge readers when:

* you are close to the field-export stage;
* you need to preserve Zonge column aliases and QC columns;
* you want to convert AVG-like data into EDI-compatible products.


J-format files
--------------

The ``pycsamt.jones`` package supports Jones-style transfer-function files.
Use it when your data source is J-format rather than EDI.

Read one J file:

.. code-block:: python
   :linenos:

   from pycsamt.jones import JFile

   jfile = JFile("data/jones/S001.j")
   print(jfile)

Validate a J file:

.. code-block:: python
   :linenos:

   from pycsamt.jones import is_j_file

   if is_j_file("data/jones/S001.j"):
       print("valid J-format file")

Read a collection:

.. code-block:: python
   :linenos:

   from pycsamt.jones import JCollection

   collection = JCollection("data/jones")

J-format is useful when:

* legacy MT/AMT processing software produced J files;
* a collaborator delivered Jones-format transfer functions;
* you need to convert J-format data toward EDI workflows.


Transformers: AVG, J, and spectra to EDI
----------------------------------------

Transformers are conversion objects.  They are useful when the input data are
valid, but the next pyCSAMT workflow expects EDI-like transfer functions.

Available public transformers include:

* ``AVGtoEDI`` for AVG-like inputs;
* ``JtoEDI`` for J-format inputs;
* ``SpectraToEDI`` for spectra-like inputs;
* ``TransformResult`` for structured conversion results.

Example:

.. code-block:: python
   :linenos:

   from pycsamt.transformers import AVGtoEDI

   transformer = AVGtoEDI()
   result = transformer.transform(
       "data/zonge/line01.avg",
       output_dir="results/edi_from_avg",
   )

   print(result)

J-format conversion:

.. code-block:: python
   :linenos:

   from pycsamt.transformers import JtoEDI

   result = JtoEDI().transform(
       "data/jones",
       output_dir="results/edi_from_j",
   )

Spectra conversion:

.. code-block:: python
   :linenos:

   from pycsamt.transformers import SpectraToEDI

   result = SpectraToEDI().transform(
       "data/spectra",
       output_dir="results/edi_from_spectra",
   )

After conversion, use :func:`pycsamt.api.read_edis` on the output directory to
enter the standard pyCSAMT survey workflow.


TDEM and TEM data
-----------------

``pycsamt.tdem`` handles time-domain electromagnetic data.  A typical TDEM
workflow starts from time gates or decay curves, attaches geometry and waveform
metadata, and then transforms the sounding into frequency-domain objects that
can join the rest of the pyCSAMT workflow.

Create a synthetic sounding:

.. code-block:: python
   :linenos:

   import numpy as np
   from pycsamt.tdem import TEMSounding, TEMtoEDI

   t = np.logspace(-5, -2, 30)
   dBdt = 5e-5 * t ** (-2.5)

   sounding = TEMSounding(
       t,
       dBdt,
       current=8.0,
       tx_area=100.0 ** 2,
       data_type="dBdt",
       station_name="S01",
       x=1000.0,
       y=500.0,
   )

   collection = TEMtoEDI(method="late_time").transform(sounding)

Read TEM AVG survey-style data:

.. code-block:: python
   :linenos:

   from pycsamt.tdem import read_temavg_survey

   survey = read_temavg_survey("data/tdem/line01.temavg")

Transform a TEM AVG survey:

.. code-block:: python
   :linenos:

   from pycsamt.tdem import transform_temavg_survey

   conversion = transform_temavg_survey(
       "data/tdem/line01.temavg",
       output_dir="results/tdem_to_edi",
   )

TDEM support also includes coordinate tables, logs, waveforms, survey maps,
decay plots, gate profiles, and transformed apparent-resistivity plots.


Site and station data
---------------------

Site and station data are not always a transfer-function format, but they are
essential for real surveys.  The ``pycsamt.site`` package handles survey
organization tasks such as:

* station objects and station collections;
* coordinate editing and location utilities;
* profile orientation and line selection;
* site export tables;
* site reports.

Use site utilities when:

* station coordinates need correction;
* you need to select a subset of stations;
* topography or profile order must be attached before inversion;
* outputs need station tables for reporting or GIS handoff.

The API-level alias :func:`pycsamt.api.read_sites` currently follows the EDI
survey reading path, so use it when your site collection is represented by EDI
sources:

.. code-block:: python
   :linenos:

   from pycsamt.api import read_sites

   sites = read_sites("data/edis")


Inversion and model files
-------------------------

Inversion formats are workflow outputs rather than raw field inputs.  pyCSAMT
groups them mainly under ``pycsamt.models`` and ``pycsamt.inversion``.

Common examples:

.. list-table::
   :header-rows: 1
   :widths: 26 34 40

   * - Engine/workflow
     - Package
     - Typical files
   * - Occam2D
     - ``pycsamt.models.occam2d``
     - data, mesh, model, startup, iteration, response, and log files.
   * - ModEM
     - ``pycsamt.models.modem``
     - data, covariance, control, model, mesh, response, and log files.
   * - Generic inversion
     - ``pycsamt.inversion``
     - configuration, mesh, model, result, uncertainty, export, and backend
       objects.

For a first inversion path, start from EDI survey data, then use the inversion
preparation tutorial:

.. code-block:: python
   :linenos:

   from pycsamt.api import read_edis

   survey = read_edis("data/edis")

See :doc:`../tutorials/prepare_occam2d_inversion` for the next step.


Choosing the right entry point
------------------------------

.. list-table::
   :header-rows: 1
   :widths: 42 58

   * - You have...
     - Start with...
   * - A folder of ``.edi`` files
     - :func:`pycsamt.api.read_edis`
   * - One EDI file and need to inspect sections
     - :class:`pycsamt.seg.EDIFile`
   * - Zonge AVG or AMTAVG exports
     - :func:`pycsamt.zonge.load_avg` or :class:`pycsamt.zonge.AVG`
   * - Jones J-format files
     - :class:`pycsamt.jones.JFile` or :class:`pycsamt.jones.JCollection`
   * - AVG/J/spectra data that should become EDI
     - ``pycsamt.transformers``
   * - TEM/TDEM gates or decay curves
     - ``pycsamt.tdem``
   * - Station tables, profile edits, or coordinate cleanup
     - ``pycsamt.site``
   * - Occam2D or ModEM preparation
     - ``pycsamt.models.occam2d`` or ``pycsamt.models.modem``


Strict versus permissive reading
--------------------------------

When loading field data, decide whether the reader should fail immediately or
keep partial results.

For EDI batches:

.. code-block:: python
   :linenos:

   from pycsamt.api import read_edis

   # Exploration: keep good files and record problems.
   survey = read_edis("data/edis", strict=False)

   # Validation: fail on unreadable or malformed files.
   survey = read_edis("data/edis", strict=True)

Use permissive mode during early exploration and strict mode before production
processing, inversion preparation, or archiving.


Duplicate stations
------------------

Field directories often contain repeated station names: reprocessed files,
alternate components, or old exports.  For EDI batches, ``on_dup`` controls the
policy:

.. code-block:: python
   :linenos:

   from pycsamt.api import read_edis

   latest = read_edis("data/edis", on_dup="replace")
   first = read_edis("data/edis", on_dup="keep")

Use ``"replace"`` when the last discovered file should win.  Use ``"keep"``
when the first discovered file should be preserved.


Units and conventions
---------------------

Always check units before mixing formats.  In pyCSAMT documentation and public
APIs, the expected defaults are:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Quantity
     - Default convention
   * - Frequency
     - Hertz.
   * - Period
     - Seconds.
   * - Apparent resistivity
     - Ohm-m.
   * - Phase
     - Degrees unless a function explicitly says radians.
   * - Distance, depth, elevation
     - Metres.
   * - Longitude/latitude
     - Decimal degrees.
   * - Projected coordinates
     - CRS/EPSG or UTM zone should be recorded in metadata.

When converting formats, keep the original files and write converted outputs
to a new directory.  This makes provenance clearer and keeps inversion
preparation reproducible.


CLI format commands
-------------------

The CLI mirrors the main format families.

.. code-block:: bash
   :linenos:

   pycsamt edi info data/edis
   pycsamt edi validate data/edis
   pycsamt avg info data/zonge/line01.avg
   pycsamt jones info data/jones/S001.j
   pycsamt tdem info data/tdem/line01.temavg
   pycsamt transform avg data/zonge/line01.avg --output-dir results/edi

See the CLI pages for command-specific options:

* :doc:`../cli/edi`
* :doc:`../cli/avg`
* :doc:`../cli/jones`
* :doc:`../cli/tdem`
* :doc:`../cli/transform`


Next steps
----------

After choosing a data format:

* read an EDI survey with :doc:`../tutorials/read_edi_survey`;
* inspect and QC a survey with :doc:`../tutorials/inspect_and_qc_survey`;
* configure output/style behavior with :doc:`configuration`;
* learn processing concepts in :doc:`../user_guide/data_loading`;
* prepare inversion files with :doc:`../tutorials/prepare_occam2d_inversion`.


In short
--------

EDI is the common path through pyCSAMT.  AVG, J, spectra, and TDEM support let
you enter that path from field and legacy formats.  Site utilities keep station
metadata coherent, and model packages write the engine-specific files needed
for inversion.
