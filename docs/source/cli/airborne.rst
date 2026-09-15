Airborne Commands
==================

``pycsamt airborne`` inspects and diagnoses airborne EM datasets — ZTEM,
MobileMT, and AFMAG/AirMt. Unlike the ground-based ``edi``/``site`` groups,
it has no raw-vendor-file reader: pyCSAMT never decodes a proprietary
airborne binary directly, only the EMTF-XML a processing vendor or a prior
``pycsamt`` step already produced (see
:doc:`/applications/converter` for one way to get EDI/AI results into a
shape other tools can read, and :doc:`/user_guide/emtf/xml` for the
EMTF-XML schema itself). Point ``airborne`` at one ``.xml`` file or a
directory of them and it does the rest: coercing them into an
:class:`~pycsamt.airborne.site.AirborneSites` collection, auto-detecting
which of the three technologies produced them, and building that
technology's canonical diagnostic table.

Command Map
-----------

.. list-table::
   :header-rows: 1
   :widths: 26 34 40

   * - Command
     - Purpose
     - Main output
   * - ``pycsamt airborne info``
     - Dataset or single-site summary report.
     - A per-site table (or one site's full record), text/JSON/CSV.
   * - ``pycsamt airborne diagnose``
     - Technology-aware diagnostic table: ZTEM total-divergence, MobileMT
       admittance, or AFMAG/AirMt tilt.
     - One ``pandas.DataFrame``, printed and optionally saved to CSV.

Both commands share one loader: :func:`~pycsamt.airborne.site.ensure_asites`
coerces ``SOURCE`` (a single EMTF-XML file, or a directory searched
recursively by default) into an
:class:`~pycsamt.airborne.site.AirborneSites` collection — the airborne
counterpart of :class:`~pycsamt.site.base.Sites`. Every site records its own
:attr:`~pycsamt.airborne.site.AirborneSite.technology`
(``ztem``/``mobilemt``/``afmag_original``/``afmag_airmt``), read straight out
of the EMTF-XML rather than guessed from a filename.

Info
----

Usage:

.. code-block:: console

   pycsamt airborne info SOURCE [--site NAME] [--no-recursive] [-f {text,json,csv}]

With no ``--site``, ``info`` summarizes every site in ``SOURCE`` — count,
the technologies and flight-line ids present, and one row per site:

.. code-block:: console

   pycsamt airborne info data/ZTEM/gold_springs_nv

which prints:

.. code-block:: text

   Sites:        105
   Technologies: ztem
   Lines:        -

        name line_id sample_id technology  nfreq       lat         lon   elev  tipper  admittance
   GO_L1_001    None GO_L1_001       ztem      6 37.984042 -114.150000 1750.0    True       False
   GO_L1_002    None GO_L1_002       ztem      6 37.984042 -114.149430 1750.0    True       False
   GO_L1_003    None GO_L1_003       ztem      6 37.984042 -114.148860 1750.0    True       False
   …

``line_id`` is ``None`` here because generic EMTF-XML ingestion never
populates it — only sites built explicitly through
:meth:`~pycsamt.airborne.site.AirborneSites.from_line` carry a real flight-line
id (see the ``--line`` note under :ref:`cli-airborne-diagnose` below for what
this means for ZTEM).

Pass ``--site`` to see one station's full record instead of the collection
table:

.. code-block:: console

   pycsamt airborne info data/ZTEM/gold_springs_nv --site GO_L1_001

which prints:

.. code-block:: text

   AirborneSite(name='GO_L1_001', technology='ztem', nfreq=6, coords=(37.98404,-114.15000,1750.0))
     name         GO_L1_001
     line_id      None
     sample_id    GO_L1_001
     technology   ztem
     nfreq        6
     lat          37.9840424002875
     lon          -114.15
     elev         1750.0
     tipper       True
     admittance   False

An unknown ``--site`` name exits non-zero and lists the sites that *are*
available, rather than a bare traceback. ``--format json`` gives the same
information machine-readably — the collection form additionally nests every
site's summary under ``sites``:

.. code-block:: console

   pycsamt airborne info data/AFMAG/abitibi_on --format json

which prints:

.. code-block:: text

   {
     "n_sites": 13,
     "technologies": [
       "afmag_original"
     ],
     "line_ids": [],
     "sites": [
       {
         "name": "AB_001",
         "line_id": null,
         "sample_id": "AB_001",
         "technology": "afmag_original",
         "nfreq": 2,
         "lat": 48.25,
         "lon": -79.3,
         "elev": 305.0,
         "tipper": false,
         "admittance": false
       },
       …

.. _cli-airborne-diagnose:

Diagnose
--------

Usage:

.. code-block:: console

   pycsamt airborne diagnose SOURCE [--technology {ztem,mobilemt,afmag_original,afmag_airmt}] \
       [--line LINE_ID] [--component {tzx,tzy}] [--spacing-m METRES] \
       [--no-recursive] [-o FILE] [-f {text,json,csv}]

``diagnose`` auto-detects the technology from the loaded sites unless
``--technology`` forces one, then dispatches to that technology's canonical
``emtools`` table function:

.. list-table::
   :header-rows: 1
   :widths: 18 30 52

   * - ``--technology``
     - Table function
     - What it measures
   * - ``ztem``
     - :func:`~pycsamt.emtools.ztem.total_divergence_table`
     - Along-line first-difference of the selected tipper component between
       adjacent stations, at every frequency — the "Total Divergence" (Lo
       and Zang, 2008) / VLF-style "Peaker" (Pedersen et al., 1994) quantity
       used to locate lateral conductivity contrasts under a flight line.
   * - ``mobilemt``
     - :func:`~pycsamt.emtools.mobilemt.admittance_table`
     - The complex admittance tensor (:math:`Y_{xx}, Y_{xy}, Y_{yx},
       Y_{yy}`, plus :math:`Y_{hzx}, Y_{hzy}`) at every frequency, with a
       native apparent-conductivity column.
   * - ``afmag_original``
     - :func:`~pycsamt.emtools.afmag.original_afmag_tilt_table`
     - Scalar tilt angle (degrees) per station/frequency — the classical
       AFMAG comparator-tilt observable (Ward, 1959).
   * - ``afmag_airmt``
     - :func:`~pycsamt.emtools.afmag.airmt_tilt_angles`
     - Real/imaginary tilt magnitude and azimuth, plus a resultant tilt —
       the richer tensor-AFMAG ("AirMt") observable, a genuinely different
       response shape from the scalar original.

A dataset with more than one technology present (a mixed survey folder) has
to be told which one to build, or ``diagnose`` refuses with ``Multiple
technologies present (…). Pass --technology to select one.`` — auto-detection
only fires when exactly one technology is found across every loaded site.

MobileMT and both AFMAG variants need no further options; running against
the bundled sample surveys:

.. code-block:: console

   pycsamt airborne diagnose data/mobileMT/flammefjeld_greenland

which prints:

.. code-block:: text

   line_id sample_id      x_m   freq_hz  period_s   Yxx_real   Yxx_imag  Yxy_real  Yxy_imag  …  apparent_conductivity_native_Sm
      L001    FL_001  0.000000  25.000000  0.040000  0.000103   0.000040 -0.002350  0.002277  …                          0.001330
      L001    FL_001  0.000000  49.642735  0.020144 -0.000060   0.000002 -0.001664  0.001617  …                          0.001235
      L001    FL_002  100.248599  25.000000  0.040000 -0.000105 -0.000123 -0.002469  0.002471  …                          0.001434
      …

.. code-block:: console

   pycsamt airborne diagnose data/AFMAG/abitibi_on

which prints a scalar tilt per station/frequency:

.. code-block:: text

   station  freq   period  tilt_deg
    AB_001 150.0 0.006667 -0.140700
    AB_001 510.0 0.001961  0.367173
    AB_002 150.0 0.006667  0.264722
    …

The tensor AFMAG/AirMt variant instead reports real/imaginary tilt
magnitude and azimuth plus a resultant:

.. code-block:: console

   pycsamt airborne diagnose data/AFMAG/yulong_belt_cn

which prints:

.. code-block:: text

   station       freq   period  tilt_real_deg  tilt_real_azimuth_deg  tilt_imag_deg  tilt_imag_azimuth_deg  tilt_resultant_deg
    YU_001  25.000000 0.040000       0.786468             178.656090       0.725171             113.584294            1.069707
    YU_001  47.204376 0.021184       0.337883             155.976973       0.706125             -46.444678            0.782786
    …

ZTEM is the one technology that needs an extra decision. Its table is a
*profile* quantity — an along-line finite difference between adjacent
stations — so it is only meaningful for one flight line at a time.
Generic EMTF-XML ingestion never populates a real ``line_id`` (see
**Info** above), so ``--line`` accepts either a real line id when one
exists, or a 0-based index into flight lines detected geometrically from
station coordinates. Running against a multi-line survey without
``--line`` refuses outright rather than differentiating across a line
boundary silently:

.. code-block:: console

   pycsamt airborne diagnose data/ZTEM/gold_springs_nv

which refuses outright:

.. code-block:: text

   Error: Detected 7 flight lines geometrically from station coordinates --
   differentiating the total-divergence table across a line boundary is not
   physically meaningful.  Pass --line with a 0-based group index (0..6) to
   select one, e.g. --line 0.

Passing the suggested index selects that one flight line and produces the
table:

.. code-block:: console

   pycsamt airborne diagnose data/ZTEM/gold_springs_nv --line 0

which prints:

.. code-block:: text

   station_a station_b        x_m      dx_m  freq_hz  period_s  divergence_real  divergence_imag  divergence_abs
   GO_L7_001 GO_L7_002  25.042321 50.084643     30.0  0.033333        -0.000499        -0.000110         0.000511
   GO_L7_002 GO_L7_003  75.126972 50.084658     30.0  0.033333         0.000367         0.000184         0.000410
   GO_L7_003 GO_L7_004 125.211638 50.084674     30.0  0.033333        -0.000393        -0.000142         0.000418
   …

``--component`` (``tzx`` default, or ``tzy``) picks which tipper component
the divergence is taken of; ``--spacing-m`` (default 200 m) only matters as
a fall-back when station coordinates are unavailable to compute real
inter-station distance directly.

Saving And Machine-Readable Output
-----------------------------------

``-o FILE`` writes the table to CSV regardless of ``--format``, printing a
confirmation line; ``--format`` (``text``/``json``/``csv``) additionally
controls what — if anything beyond that confirmation — prints to the
console:

.. code-block:: console

   pycsamt airborne diagnose data/AFMAG/abitibi_on -o afmag_table.csv --format csv

which prints:

.. code-block:: text

   Table saved → afmag_table.csv
   station,freq,period,tilt_deg
   AB_001,149.99999999999994,0.006666666666666669,-0.14069988133109063
   AB_001,510.0,0.00196078431372549,0.3671731486153681
   …

Common Failures
----------------

``Path '...' does not exist``
    ``SOURCE`` is a Click path argument checked before anything else runs —
    a typo or a relative path resolved from the wrong working directory is
    the usual cause.

``No airborne sites found under ...``
    ``SOURCE`` exists but contains no EMTF-XML files
    :func:`~pycsamt.airborne.site.ensure_asites` recognises (or the search
    was not recursive — check ``--no-recursive`` was not passed by
    accident for a nested directory).

``Site '...' not found.  Available: ...``
    ``info --site`` was given a name that does not match any loaded site;
    the message lists what *is* available so there is no need to re-run
    without ``--site`` just to check.

``Multiple technologies present (...).  Pass --technology to select one.``
    ``diagnose`` auto-detects only when the whole loaded collection agrees
    on one technology. A folder mixing surveys needs ``--technology``, or
    should be split into one ``SOURCE`` per technology first.

``Could not detect a technology for any site under ...``
    None of the loaded sites reported a recognised
    :attr:`~pycsamt.airborne.site.AirborneSite.technology` — pass
    ``--technology`` explicitly if you know which one produced the data.

``Detected N flight lines geometrically … Pass --line with a 0-based group index``
    ZTEM-specific — see :ref:`cli-airborne-diagnose` above. Re-run with
    ``--line`` set to one of the suggested indices (or a real ``line_id``
    if the dataset has one).

``No sites found for --line '...' (checked both real line_id metadata and the N geometrically detected flight line(s) …)``
    The ``--line`` value matched neither a real flight-line id nor a valid
    geometric group index for this dataset.

``No diagnostic rows produced from this dataset.``
    The technology's table function ran without raising but returned an
    empty table — usually a dataset with the right technology tag but no
    usable frequencies/tipper data for that table (a warning, not a hard
    failure; the command still exits ``0``).

Python Equivalents
-------------------

The CLI is a thin layer over :mod:`pycsamt.airborne.site` and the
per-technology ``emtools`` table functions:

.. code-block:: python

   from pycsamt.airborne.site import ensure_asites
   from pycsamt.emtools.afmag import original_afmag_tilt_table

   sites = ensure_asites("data/AFMAG/abitibi_on")
   print(sites.technologies, len(sites))          # ('afmag_original',) 13
   table = original_afmag_tilt_table(sites)

For ZTEM, pre-select one flight line the same way ``--line`` does before
calling the table function directly — passing a multi-line collection
straight to :func:`~pycsamt.emtools.ztem.total_divergence_table` produces a
physically meaningless value at every line-to-line join rather than
raising:

.. code-block:: python

   from pycsamt.airborne.site import ensure_asites
   from pycsamt.emtools.ztem import total_divergence_table

   sites = ensure_asites("data/ZTEM/gold_springs_nv")
   line = sites.select(predicate=lambda s: s.line_id == "L7")  # or slice by index
   table = total_divergence_table(line, spacing_m=200.0, component="tzx")

.. seealso::

   :doc:`/applications/converter`
       ``pycsamt format edi-to-xml``/``xml-to-edi`` and the standalone GUI —
       the other half of getting data into or out of EMTF-XML.

   :doc:`/user_guide/emtf/xml`
       The EMTF-XML schema these commands read.
