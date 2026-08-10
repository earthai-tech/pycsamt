.. _getting-started-data-formats:

Identify your data format
=========================

Before loading a survey, identify what the files contain rather than relying
only on their extension. pyCSAMT can read frequency-domain EDI transfer
functions directly. Instrument exports, cross-power spectra, field time
series, and time-domain electromagnetic decays need a format-specific reader
or conversion step before they enter the same processing workflow.

Keep the original field delivery unchanged. Write converted or corrected data
to a separate directory and record the source files, coordinate reference
system, units, processing settings, and pyCSAMT version. Format conversion
changes representation; it does not prove that the observations are complete
or scientifically valid.

Choose the nearest input family
-------------------------------

.. list-table::
   :header-rows: 1
   :widths: 19 31 25 25

   * - Input family
     - How to recognize it
     - Start with
     - Continue with
   * - EDI transfer functions
     - Text files, usually ``.edi``, with EDI sections and impedance or
       resistivity/phase values indexed by frequency.
     - :func:`pycsamt.api.read_edis`
     - :doc:`first_survey`
   * - Zonge AVG or AMTAVG
     - Instrument-averaged station/component rows containing frequency,
       apparent resistivity, phase, and quality fields.
     - :class:`pycsamt.zonge.avg.AVG` or
       :class:`pycsamt.zonge.avg.AMTAVG`
     - :doc:`../user_guide/transformers`
   * - Jones J-format
     - A station-oriented transfer-function file with Jones-format headers and
       component blocks; the suffix is not always ``.j``.
     - :class:`pycsamt.jones.j.JFile`
     - :doc:`../user_guide/transformers`
   * - Spectral EDI
     - An EDI container with cross-power spectral blocks but no final
       impedance tensor.
     - :class:`pycsamt.transformers.SpectraToEDI`
     - :doc:`../user_guide/transformers`
   * - Electromagnetic time series
     - Synchronized electric- and magnetic-field samples rather than
       frequency-domain transfer functions.
     - :class:`pycsamt.transformers.TStoEDI`
     - :doc:`../user_guide/transformers`
   * - TEM/TDEM decay data
     - Time gates and voltage or :math:`\mathrm{d}B/\mathrm{d}t` decay values,
       often accompanied by waveform and transmitter geometry.
     - :mod:`pycsamt.tdem`
     - :doc:`../user_guide/transformers`
   * - Stratagem field delivery
     - Geometrics/EMI raw component files, hardware exports, coordinate tables,
       or a WinGLink hand-off.
     - :class:`pycsamt.stratagem.StratagemSurvey`
     - :doc:`../user_guide/stratagem/index`
   * - Station or coordinate table
     - CSV or tabular station identifiers, longitude/latitude, projected
       coordinates, elevation, or profile position without transfer functions.
     - :mod:`pycsamt.site`
     - :doc:`../user_guide/site/index`
   * - Inversion input or result
     - Mesh, model, data, startup, response, iteration, covariance, or solver
       log files.
     - :mod:`pycsamt.models`
     - :doc:`../user_guide/inversion/index`

Station tables are normally auxiliary metadata, not electromagnetic
observations. They become scientifically useful only after station identities
are matched unambiguously to transfer functions or time-domain soundings.
Likewise, inversion files are downstream products; do not pass them to a field
data reader simply because they contain station names or resistivity values.

Load EDI data directly
----------------------

EDI is the common interchange format for MT, AMT, and CSAMT
frequency-domain responses. A useful EDI file normally identifies the station
and frequencies and contains an impedance tensor or an equivalent
apparent-resistivity/phase representation. Tipper and spectral sections are
optional and should not be assumed present.

The bundled WILLY L18PLT line provides a reproducible first check:

.. code-block:: pycon

   >>> from pycsamt.api import read_edis
   >>> survey = read_edis(
   ...     "data/AMT/WILLY_DATA/L18PLT",
   ...     recursive=False,
   ...     strict=True,
   ...     on_dup="replace",
   ...     progress=False,
   ... )
   >>> survey.n_sites
   28
   >>> survey.stations[:3]
   ['23-18-001A', '23-18-002U', '23-18-003A']
   >>> survey.summary().shape
   (28, 6)

``recursive=False`` keeps discovery inside the selected survey-line directory.
Set it to ``True`` only when nested directories belong to the same intended
load. ``strict=True`` makes an unresolved input fail instead of returning an
empty survey, while ``on_dup="replace"`` retains the later file when the EDI
parser encounters the same station identity more than once.

.. warning::

   A successful EDI parse establishes structural readability, not data
   quality. Before processing, compare the station count with the field
   manifest and inspect coordinates, component availability, frequency
   coverage, errors, and duplicate identities. The complete loading checks are
   described in :doc:`../user_guide/data_loading`.

Convert transfer-function formats when necessary
------------------------------------------------

AVG and J-format files already describe frequency-domain responses, but their
metadata and component layout differ from EDI. Convert them through the public
transformers rather than renaming the files or manually rearranging columns:

.. code-block:: pycon

   >>> from pycsamt.transformers import AVGtoEDI, JtoEDI
   >>> from pycsamt.zonge.avg import AVG
   >>> avg = AVG.from_file("data/avg/K1.AVG")
   >>> edi_collection = AVGtoEDI().transform(avg)
   >>> len(edi_collection)
   47

For J-format input, :class:`~pycsamt.transformers.JtoEDI` follows the same
conversion contract. The suffix alone is unreliable: the bundled example
``data/j/kb0-s001.txt`` is a J-format file despite its ``.txt`` extension.

Conversion requires more than reproducing numeric columns. Station identity,
component convention, frequency ordering, coordinates, missing-value rules,
and error estimates must remain traceable. The complete, verified AVG, J,
spectra, time-series, and TEM examples—including their limitations—are in
:doc:`../user_guide/transformers`.

Distinguish spectra from impedance EDI
--------------------------------------

An EDI file can contain raw or processed cross-power spectra instead of a
finished impedance tensor. Spectral conversion estimates the transfer
function from electric/magnetic cross-spectra and magnetic auto-spectra; it is
not a file-copy operation. Use :class:`~pycsamt.transformers.SpectraToEDI` only
when a genuine spectral block is present.

Conversely, do not send an ordinary impedance EDI through the spectral
transformer. Load it directly with :func:`pycsamt.api.read_edis`. The bundled
``data/MT/SPECTRA/spectra01.edi`` file is available for the worked spectral
example in :doc:`../user_guide/transformers`.

Treat time-domain data as a separate measurement family
--------------------------------------------------------

Time series and TEM/TDEM decays are not interchangeable with
frequency-domain EDI values. A field time series records electric and magnetic
channels as functions of sample time; estimating an impedance requires
windowing, spectral estimation, robust processing, and uncertainty choices.
A TEM/TDEM sounding instead records a transient response over time gates and
requires its transmitter geometry, current, waveform, units, and component
convention.

Do not infer missing acquisition metadata from a filename. An incorrect coil
area, current, time unit, or response convention can scale the transformed
result while leaving it numerically smooth. Use :class:`~pycsamt.transformers.TStoEDI`
for supported field time series and :mod:`pycsamt.tdem` for transient data,
then validate the result against the acquisition record.

Check units and spatial references
----------------------------------

Before combining files, establish the conventions used by each source:

.. list-table::
   :header-rows: 1
   :widths: 27 30 43

   * - Quantity
     - Common documentation convention
     - Check before conversion
   * - Frequency / period
     - Hz / s
     - Whether the file stores frequency, period, or both, and their ordering.
   * - Apparent resistivity
     - :math:`\Omega\,\mathrm{m}`
     - Whether values are linear or logarithmic and which component they
       represent.
   * - Phase
     - degrees
     - Sign, quadrant, component convention, and whether radians are used.
   * - Distance and elevation
     - m
     - Unit scale, vertical datum, and whether elevation is measured or
       interpolated.
   * - Geographic coordinates
     - decimal degrees
     - Longitude/latitude order, hemisphere, datum, and valid ranges.
   * - Projected coordinates
     - m
     - CRS/EPSG identifier, UTM zone, hemisphere, and axis order.
   * - Time-domain response
     - API-specific SI units
     - Voltage versus field response, :math:`B` versus
       :math:`\mathrm{d}B/\mathrm{d}t`, time unit, and normalization.

The table gives common pyCSAMT documentation conventions, not permission to
assume them for an unknown field file. Preserve the instrument header and
processing report, and state every explicit conversion in the project record.

Continue from the identified format
-----------------------------------

* If the data are EDI, continue to :doc:`first_survey` and then
  :doc:`../user_guide/data_loading`.
* If the data are AVG, J-format, spectral EDI, time series, or TEM/TDEM, use
  :doc:`../user_guide/transformers` before loading the converted collection.
* If the delivery comes from Stratagem hardware, follow
  :doc:`../user_guide/stratagem/index`.
* If the files are inversion inputs or results, begin with
  :doc:`../user_guide/inversion/index` rather than a survey reader.

Once the data family is clear, :doc:`configuration` establishes output and
display defaults, and :doc:`first_survey` performs the first survey-level
inspection.
