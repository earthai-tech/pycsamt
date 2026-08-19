.. _user_guide_metadata_frequency_bands:

Frequency Bands
===============

.. currentmodule:: pycsamt.metadata.frequency

Every EM method operates over a characteristic frequency range, and
knowing that range matters for survey design, quality review, and
depth-of-investigation expectations. :class:`FrequencyBand` captures
one method's operational frequency range together with its period
bounds and a rule-of-thumb depth-of-investigation estimate from the
skin-depth relationship. :data:`MT_BANDS` is a ready-made registry
covering the principal EM methods used in applied geophysics, and
:func:`register_band` lets a project add its own named bands without
modifying the package. Unlike most classes covered earlier in this
section, ``FrequencyBand`` is a plain dataclass, not a
:class:`~pycsamt.api.property.PyCSAMTObject` subclass, and it is
general-purpose rather than tied to the EMTF document model -- nothing
here maps to an EMTF XML element.

Built-In Bands
--------------

:data:`MT_BANDS` and the mutable :data:`REGISTRY` (seeded from it)
currently hold eight named bands:

.. code-block:: pycon

   >>> from pycsamt.metadata.frequency import MT_BANDS, REGISTRY

   >>> print(sorted(REGISTRY.keys()))
   ['AMT', 'BBMT', 'CSAMT', 'CSEM', 'LAMT', 'LMT', 'MT', 'TEM']

   >>> amt = MT_BANDS["AMT"]
   >>> print(amt.f_min, amt.f_max)
   10.0 100000.0
   >>> print(amt.period_range)
   (1e-05, 0.1)
   >>> print(round(amt.n_decades, 2))
   4.0

:meth:`~FrequencyBand.doi_range_m` turns the band into a rough
depth-of-investigation range via the skin-depth formula
:math:`\delta \approx 503.3\sqrt{\rho/f}`, evaluated at ``f_max``
(shallow) and ``f_min`` (deep). The reference resistivity defaults to
100 :math:`\Omega \cdot \text{m}` but can be overridden -- a more
conductive earth attenuates the same frequency faster, so a lower
resistivity always narrows both DOI bounds:

.. code-block:: pycon

   >>> print(amt.doi_range_m())
   (15.9, 1591.6)
   >>> print(amt.doi_range_m(rho=10.0))
   (5.0, 503.3)
   >>> print(amt.doi_range_m(rho=1000.0))
   (50.3, 5033.0)
   >>> print(repr(amt))
   FrequencyBand('AMT'  10–1e+05 Hz  DOI≈16–1592 m  [AMT])

Which Band Covers a Frequency
-----------------------------

:func:`band_for_frequency` returns every band containing a query
frequency, narrowest first. Applied to the real frequency endpoints of
the two stations used throughout this section, it shows how a single
sounding can straddle several named bands at once:

.. code-block:: pycon

   >>> from pathlib import Path
   >>> from pycsamt.seg.edi import EDIFile
   >>> from pycsamt.site.base import Site
   >>> from pycsamt.metadata.frequency import band_for_frequency

   >>> site_willy = Site(EDIFile(Path("data/AMT/WILLY_DATA/L18PLT/18-001A.edi")))
   >>> freq_willy = site_willy.freq
   >>> print(freq_willy.min(), freq_willy.max(), freq_willy.size)
   1.008 10400.0 53

   >>> print([b.name for b in band_for_frequency(freq_willy.min())])
   ['LAMT', 'CSAMT', 'TEM', 'MT', 'BBMT']
   >>> print([b.name for b in band_for_frequency(freq_willy.max())])
   ['AMT']

WILLY's lowest frequency, 1.008 Hz, falls inside five different named
bands at once -- it is squarely in the overlap where LAMT, CSAMT, TEM,
MT, and BBMT all claim coverage. Its highest frequency, 10.4 kHz, is
narrow enough to belong to ``AMT`` alone. The Gabbs Valley station
used in :doc:`processing_and_quality` spans an even wider range and
shows the opposite extreme at its low end:

.. code-block:: pycon

   >>> site_gv = Site(EDIFile(Path("data/gv_data/gv_final_edi/gv100.edi")))
   >>> freq_gv = site_gv.freq
   >>> print(freq_gv.min(), freq_gv.max(), freq_gv.size)
   0.0004882812 767.9902 48

   >>> print([b.name for b in band_for_frequency(freq_gv.min())])
   ['LMT', 'TEM', 'MT', 'BBMT']
   >>> print([b.name for b in band_for_frequency(freq_gv.max())])
   ['AMT', 'LAMT', 'CSAMT', 'MT', 'BBMT']

.. note::

   ``gv100.edi`` above is public-domain USGS data; see
   :doc:`provenance_and_bibliography` for the required citation.

Neither result names the survey's own acquisition method -- both files
are conventional broadband MT -- ``band_for_frequency`` answers "which
*named bands* contain this frequency," not "which method recorded it."

Overlap and Clipping
--------------------

:meth:`~FrequencyBand.overlaps` and :meth:`~FrequencyBand.intersection`
compare two bands directly, and :meth:`~FrequencyBand.clip_frequencies`
filters a real frequency array down to one band's range:

.. code-block:: pycon

   >>> csamt = MT_BANDS["CSAMT"]
   >>> mt = MT_BANDS["MT"]
   >>> print(amt.overlaps(csamt), amt.intersection(csamt))
   True (10.0, 10000.0)
   >>> print(amt.overlaps(mt), amt.intersection(mt))
   True (10.0, 1000.0)

   >>> clipped = amt.clip_frequencies(freq_willy)
   >>> print(clipped.size, freq_willy.size)
   40 53

Thirteen of WILLY's 53 frequencies -- everything below 10 Hz -- fall
outside the ``AMT`` band and are dropped by
:meth:`~FrequencyBand.clip_frequencies`. That is expected: WILLY is a
broadband AMT-class survey extending below the nominal AMT floor, not
a sign that a quarter of the sounding is somehow invalid. Use
clipping to restrict a plot or an inversion sounding to one named
band's nominal range, not as a data-quality filter --
:doc:`processing_and_quality`'s :class:`~pycsamt.metadata.quality.DataQuality`
is the right tool for that.

Extending the Registry
----------------------

:func:`register_band` adds a custom :class:`FrequencyBand` to the
global :data:`REGISTRY` without touching the package. The committed
:mod:`pycsamt.airborne` package (see the top-level user guide) has no
frequency-band registry of its own, so a project working with a
MobileMT-style system might register its two frequency stacks
directly:

.. code-block:: pycon

   >>> from pycsamt.metadata.frequency import FrequencyBand, register_band, frequency_range

   >>> low = FrequencyBand(
   ...     name="MOBILEMT_LOW", label="MobileMT low stack",
   ...     f_min=25.0, f_max=6000.0, method="AIRBORNE",
   ... )
   >>> high = FrequencyBand(
   ...     name="MOBILEMT_HIGH", label="MobileMT high stack",
   ...     f_min=6000.0, f_max=115000.0, method="AIRBORNE",
   ... )
   >>> register_band(low)
   >>> register_band(high)
   >>> print(frequency_range("AIRBORNE"))
   (25.0, 115000.0)

That result is the method-level union :func:`frequency_range`'s
docstring promises: two bands sharing ``method="AIRBORNE"``, neither
named exactly ``"AIRBORNE"``, so the query falls through to the union
branch. That branch is easy to miss in the *built-in* registry,
because it never actually fires there -- every ``MT_BANDS`` band's
``method`` string is also the name of a band in its own right (``AMT``
is both a band and a method, and so is every other entry), so
``frequency_range()``'s direct name lookup always wins first and the
union path is dead code for the shipped registry. It only becomes
observable once a project registers multiple bands under one shared
method that is not itself a band name, as above.

Registering with an uppercase name also matters more than it looks.
``frequency_range()`` upper-cases its query before checking
``REGISTRY`` for a matching key, which works seamlessly for the
built-in bands because they are all uppercase already. A band
registered with a lowercase or mixed-case name is invisible to
*name*-based lookup through ``frequency_range()``, under either
casing of the query -- even though it is genuinely present in
``REGISTRY`` and fully reachable through :func:`band_for_frequency`
or direct dict access:

.. code-block:: pycon

   >>> lower = FrequencyBand(
   ...     name="mobilemt_mid", label="a lowercase-keyed band",
   ...     f_min=1000.0, f_max=2000.0, method="AIRBORNE",
   ... )
   >>> register_band(lower)
   >>> print("mobilemt_mid" in REGISTRY)
   True
   >>> try:
   ...     frequency_range("mobilemt_mid")
   ... except KeyError as exc:
   ...     print(type(exc).__name__)
   KeyError

Follow the built-in convention -- register custom bands with an
uppercase ``name`` -- unless the only access pattern needed is
``REGISTRY[name]`` or :func:`band_for_frequency`, neither of which
uppercases anything.

Choosing the Right Function
---------------------------

.. list-table::
   :header-rows: 1
   :widths: 32 30 38

   * - Need
     - Function/method
     - Notes
   * - Look up a named band's range
     - ``REGISTRY[name]`` or :func:`frequency_range`
     - ``frequency_range`` needs an uppercase name to match by name.
   * - Which bands contain a frequency
     - :func:`band_for_frequency`
     - Sorted narrowest to widest; scans by range, no name-key lookup.
   * - Rough depth of investigation
     - :meth:`~FrequencyBand.doi_range_m` / :func:`doi_estimate`
     - Skin-depth estimate only; not a substitute for real forward modelling.
   * - Filter an array to one band
     - :meth:`~FrequencyBand.clip_frequencies`
     - A range filter, not a data-quality tool -- see :doc:`processing_and_quality`.
   * - Add a project-specific band
     - :func:`register_band`
     - Use an uppercase ``name`` to match the built-in convention.

Next Steps
----------

This closes the current pass through :doc:`index`. The remaining
:mod:`pycsamt.metadata` modules -- ``geology`` (resistivity priors for
synthetic layered-earth generation) and ``rocks`` -- are forward-
modelling scaffolding rather than EMTF-adjacent metadata, and are
deliberately not covered here; see the module docstrings directly if
needed. From here:

* :doc:`site_and_survey` through :doc:`processing_and_quality` cover
  the format-neutral metadata that populates an
  :class:`~pycsamt.emtf.document.EMTF` document;
* the next step in the broader documentation roadmap is a deep pass
  through ``user_guide/emtf/``, which can now cross-reference this
  section instead of re-explaining these classes.
