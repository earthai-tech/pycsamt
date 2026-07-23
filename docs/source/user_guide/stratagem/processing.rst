.. _user_guide_stratagem_processing:

Static Shift and Noise Removal
==============================

:doc:`quality_control` cleared out frequencies that never deserved
trust in the first place. What survives is not automatically finished:
two more systematic effects sit in the impedance tensor of almost any
AMT survey, and neither is a bad-frequency problem that masking can
fix. :term:`Static shift` is a frequency-independent vertical offset of
the apparent-resistivity curve caused by small near-surface
conductivity heterogeneities near a station -- it does not go away by
dropping frequencies, because it affects every frequency at that
station by the same multiplicative amount.
:class:`~pycsamt.stratagem.process.StaticShiftCorrector` estimates and
removes it station by station. Powerline hum and isolated spectral
spikes are a separate, narrower-band problem that
:class:`~pycsamt.stratagem.process.NoiseRemover` handles instead. Both
classes delegate their real algorithms to :mod:`pycsamt.emtools.ss` and
:mod:`pycsamt.emtools.remove_noise` and modify ``EDIFile.Z.z``
in-place, the same ``fit()``/``out()`` shape as every other class in
this section.

AMA Spatial Filter
------------------

Static shift is galvanic and therefore has no signature in frequency,
but it does have one in *space*: a station sitting on an anomalous
patch of near-surface resistivity looks shifted relative to its
neighbours along the line, not relative to itself across frequency.
:class:`~pycsamt.stratagem.process.StaticShiftCorrector` exploits
exactly that, estimating each station's shift as its deviation from a
spatial moving average of log-resistivity over ``half_window``
neighbours on each side, ordered by ``sort_by`` -- ``'lon'`` for a
roughly E-W profile like K2's, ``'lat'`` for a N-S one. ``weights``
controls how neighbour influence falls off with distance (``'tri'``,
``'gauss'``, or ``'uniform'``), and an optional ``pband`` restricts the
estimate to a period range presumed free of near-surface distortion.
Continuing directly from :doc:`quality_control`'s injected 86-station
batch:

.. code-block:: pycon

   >>> import pandas as pd
   >>> from pathlib import Path
   >>> from tempfile import TemporaryDirectory
   >>> from pycsamt.stratagem import EDIBatch, CoordinateInjector
   >>> from pycsamt.stratagem.process import StaticShiftCorrector

   >>> aligned = pd.read_csv("data/stratagem/K2/k2-gps-aligned.csv")
   >>> survey_coords = aligned[aligned["use_for_survey"]]
   >>> with TemporaryDirectory() as tmp:
   ...     coord_csv = Path(tmp) / "coords.csv"
   ...     survey_coords.to_csv(coord_csv, index=False)
   ...     batch = EDIBatch("data/stratagem/K2/k2-edi").fit()
   ...     edi_objects = [e for i, e in enumerate(batch.edi_objects_) if i != 0]
   ...     injector = CoordinateInjector(epsg=32649, order="forward").fit(
   ...         edi_objects, coord_csv,
   ...         easting_col="easting", northing_col="northing",
   ...         elev_col="elev", station_col="edi_file",
   ...     )

   >>> sc = StaticShiftCorrector(sort_by="lon").fit(injector.edi_objects_)
   >>> sc.factors_.head(3)
      station  delta_log10_rho   fac_rho     fac_z  n_used
   0  Z2HX002         0.438309  0.364495  0.603734       7
   1  Z2HX003         0.366937  0.429599  0.655438       9
   2  Z2HX004        -0.045029  1.109249  1.053209       4

``fac_z`` is what actually multiplies each station's impedance; a
station reading 0.60 was overestimating resistivity by roughly a
factor of 1/0.60² once converted back through :math:`\rho_a \propto
|Z|^2`, relative to its neighbours' trend. ``n_used`` -- the number of
neighbour frequencies that actually contributed -- is the first thing
worth reading alongside a factor, not after it.

Skew-Based Station Exclusion
----------------------------

``sc.factors_`` has fewer rows than stations went in:

.. code-block:: pycon

   >>> len(injector.edi_objects_), len(sc.factors_)
   (86, 74)

The 12 missing stations are not sitting in the table with a null or
unit factor -- they are absent entirely. ``max_skew=6.0`` (the same
default seen in :doc:`quality_control`) excludes individual
*frequencies* whose phase-tensor skew exceeds 6°, not whole stations;
a station only disappears here once every one of its frequencies fails
that per-frequency screen. That is a different computation from
:class:`~pycsamt.stratagem.qc.QualityController`'s ``skew_med`` column,
which is a single per-station median across the whole frequency range
-- the two share a threshold value and a physical quantity, phase-tensor
skew, but not the same aggregation, so a station flagged ``high_skew``
in the QC report does not automatically vanish from ``factors_`` here,
and vice versa. Silent disappearance rather than a printed warning is
exactly why ``len(sc.factors_)`` is worth checking against the input
count before assuming every station was corrected.

Validating Correction Factors
-----------------------------

Even among the 74 stations that do get a factor, not all factors carry
equal weight:

.. code-block:: pycon

   >>> sc.factors_["fac_z"].describe().round(2)
   count    74.00
   mean      2.83
   std      10.97
   min       0.22
   25%       0.70
   50%       0.95
   75%       1.23
   max      93.03
   Name: fac_z, dtype: float64
   >>> sc.factors_.sort_values("fac_z", ascending=False).head(3)
       station  delta_log10_rho      fac_rho      fac_z  n_used
   52  Z2HX061        -3.937243  8654.511051  93.029625       1
   68  Z2HX082        -2.543103   349.222879  18.687506       1
   51  Z2HX060        -2.054870   113.467160  10.652097       1

The median correction is a modest 0.95 -- most of the line needed
almost no adjustment -- but station 61's factor is 93, built from
``n_used=1``: a single neighbour frequency, not a spatial trend. A
2-order-of-magnitude impedance correction resting on one data point is
not a static-shift estimate worth trusting blindly; it is a station to
inspect directly (its own resistivity curve, its neighbours', the
terrain notes in the field report) before deciding whether to accept,
discard, or manually override that factor. ``fac_z`` values close to 1
with a healthy ``n_used`` are the ones to trust without a second look;
everything else on this list is exactly what
:attr:`~pycsamt.stratagem.process.StaticShiftCorrector.factors_`
exists to let you find before the correction is baked into the
exported EDI.

Processing Order Dependency
---------------------------

The fit above ran on the coordinate-injected batch directly, before
any frequency masking. Running the same correction *after*
:doc:`quality_control`'s :class:`~pycsamt.stratagem.qc.FrequencyFilter`
instead -- on data that already has ``NaN`` gaps punched into it --
breaks the AMA estimation outright:

.. code-block:: pycon

   >>> from pycsamt.stratagem.qc import FrequencyFilter
   >>> from pycsamt.stratagem import StratagemRawReader
   >>> from copy import deepcopy

   >>> rdr = StratagemRawReader("data/stratagem/K2/k2-HX", component="X").fit()
   >>> edis_wrong = [deepcopy(e) for e in injector.edi_objects_]
   >>> filt = FrequencyFilter(fmin=10.0, fmax=10000.0).fit(edis_wrong, raw_reader=rdr)
   >>> sc_wrong = StaticShiftCorrector(sort_by="lon").fit(filt.edi_objects_)
   [StaticShiftCorrector] WARNING — AMA estimation failed (IndexError: index 23 is out of bounds for axis 0 with size 23).
     Tip: run StaticShiftCorrector BEFORE FrequencyFilter so that Z data is complete during spatial averaging.
     No static-shift correction applied.
   >>> sc_wrong.factors_.head(3)
      station  delta_log10_rho  fac_rho  fac_z  n_used
   0  Z2HX002              0.0      1.0    1.0       0
   1  Z2HX003              0.0      1.0    1.0       0
   2  Z2HX004              0.0      1.0    1.0       0
   >>> sc_wrong.factors_["n_used"].unique().tolist()
   [0]

Nothing raises past ``fit()`` -- the exception is caught internally, a
warning prints, and every station gets ``fac_z=1.0``, ``n_used=0``: a
correction table that looks complete but did nothing at all. Compare
that to the real result two sections above, where the same 86 stations
produced a genuine spread of factors and a median close to but not
exactly 1. A silently-uncorrected batch and a genuinely-uncorrected-by-
almost-nothing batch look identical unless ``n_used`` is checked, which
is the real, concrete reason
:class:`~pycsamt.stratagem.process.StaticShiftCorrector` has to run
*before* :class:`~pycsamt.stratagem.qc.FrequencyFilter` in the
pipeline, not the other way around -- exactly the ordering
:class:`~pycsamt.stratagem.survey.StratagemSurvey` enforces in
:doc:`pipeline`.

Noise Removal Stages
--------------------

:class:`~pycsamt.stratagem.process.NoiseRemover` runs three stages in
order on the statically-shift-corrected batch: a powerline notch at
``mains_hz`` and its harmonics, a Hampel median-absolute-deviation
outlier filter, and an optional log-frequency smoothing pass.
Continuing with ``sc.edi_objects_`` (the correctly-ordered, real
correction from the first section):

.. code-block:: pycon

   >>> import numpy as np
   >>> from pycsamt.stratagem.process import NoiseRemover

   >>> z_before = deepcopy(sc.edi_objects_[0].Z.z)
   >>> nr = NoiseRemover(mains_hz=50.0).fit(sc.edi_objects_)
   >>> np.allclose(np.abs(z_before), np.abs(nr.edi_objects_[0].Z.z), equal_nan=True)
   False
   >>> float(np.nansum(np.abs(np.abs(z_before) - np.abs(nr.edi_objects_[0].Z.z))))
   5097.653483673466

Station 2's impedance amplitudes genuinely changed under the notch and
Hampel stages -- ``mains_hz=50.0`` is the right mains frequency for
this survey region; North American data would need ``60.0`` instead.
The third, optional stage is where the class's own docstring carries a
documented caveat: ``smooth_win`` values above 4 can hit a shape
mismatch against a short frequency vector, and K2's stations run only
33-39 frequencies after band selection -- short enough to hit it in
practice, not just in theory:

.. code-block:: pycon

   >>> nr2 = NoiseRemover(mains_hz=50.0, smooth=True, smooth_win=5, verbose=1)
   >>> nr2.fit(deepcopy(sc.edi_objects_))
   [NoiseRemover] smooth_logfreq skipped (win=5): could not broadcast input array from shape (9,) into shape (8,)
   [NoiseRemover] processed 86 stations (notch + hampel + smooth)
   NoiseRemover(mains_hz=50.0, n_harm=30, hampel_win=3, smooth=True, n_stations_=86)

``fit()`` still completes -- the mismatch is caught around the
smoothing call specifically, notch and Hampel filtering both already
ran and stay applied, and only the smoothing stage is skipped, with a
message rather than a silent no-op. ``smooth_win=3`` (the class
default) is small enough to stay clear of this on K2's frequency
counts; raise it only after checking how many frequencies a station
actually has left post-filtering. :doc:`export_rename` picks up from
here, writing and renaming this corrected, denoised batch to its final
output form.
