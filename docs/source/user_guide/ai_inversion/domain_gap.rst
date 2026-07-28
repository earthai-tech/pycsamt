.. _ai_inversion_domain_gap:

Domain-gap and noise simulation
================================

A geological prior generator plus a verified forward solver, as
described in :doc:`geology_priors` and :doc:`dataset2d`, produces
clean, noiseless responses. Field data is never that well behaved.
:mod:`pycsamt.ai.domain_gap` exists to close some of that
:term:`domain gap` deliberately, by turning a clean
:class:`~pycsamt.ai.data.contracts.SurveyData` into one with realistic
:term:`error floor`, heteroscedastic noise, :term:`static shift`,
:term:`galvanic distortion`, missing stations or frequencies, and
undetected outliers — instead of hoping that a generic noise term
covers what a real survey actually looks like.

Two families of corruption
------------------------------

:mod:`~pycsamt.ai.domain_gap.simulator` splits its corruption
functions into two families that behave differently on purpose.
Systematic effects multiply the true impedance by a per-station real
matrix and leave the shape of the declared error to scale with it.
:func:`~pycsamt.ai.domain_gap.simulator.apply_static_shift` draws a
station's shift factor log-normally,

.. math::
   :label: eq-ai-domaingap-static-shift

   c_s = 10^{\,X_s}, \qquad X_s \sim \mathcal N(0, \sigma_{\log_{10}}^2),

so that :math:`\sigma_{\log_{10}}` (``static_shift_log10_sigma``) is a
decade-scale spread rather than a linear one, matching how static
shift is normally reported:

.. code-block:: pycon

   >>> import numpy as np
   >>> from pycsamt.ai.data.contracts import SurveyData
   >>> from pycsamt.ai.domain_gap.simulator import apply_static_shift

   >>> z = np.ones((3, 2, 2), dtype=complex) * (100 + 100j)
   >>> survey = SurveyData(
   ...     z, [10.0, 1.0], ["A", "B", "C"], ["xy", "yx"],
   ...     [[0, 0], [1, 0], [2, 0]],
   ... )
   >>> shifted, info = apply_static_shift(
   ...     survey, log10_sigma=0.15, rng=np.random.default_rng(2),
   ...     return_info=True,
   ... )
   >>> np.round(info["static_shift_factor"], 4)
   array([1.0675, 0.8348, 0.867 ])

:func:`~pycsamt.ai.domain_gap.simulator.apply_galvanic_distortion`
follows the same idea with a full Groom-Bailey-style real 2x2 matrix
per station,

.. math::
   :label: eq-ai-domaingap-distortion-matrix

   D_s = g_s\,R(\theta_s)
   \begin{bmatrix} 1+e_s & s_s \\ s_s & 1-e_s \end{bmatrix},

combining a gain :math:`g_s` (log-normal, like the static-shift
factor), a rotation :math:`R(\theta_s)` by a twist angle
:math:`\theta_s`, a shear :math:`s_s`, and an anisotropy :math:`e_s` —
each independently normal with its own configured standard deviation
— and applies :math:`D_s` to every frequency of that station's
impedance at once. Systematic effects change the signal; they do not,
by themselves, say anything about how much to trust it.

Random effects are the opposite: they update
``impedance_error`` and the validity mask directly, so a corrupted
survey's declared uncertainty stays honest about what was actually
injected.
:func:`~pycsamt.ai.domain_gap.simulator.add_heteroscedastic_noise`
adds complex Gaussian noise per ``(station, frequency)`` pair scaled
by a relative level drawn uniformly from a configured range, and
combines it with any pre-existing declared error in quadrature rather
than overwriting it. :func:`~pycsamt.ai.domain_gap.simulator.apply_error_floor`
clamps the declared error to a minimum fraction of :math:`|Z|`.
:func:`~pycsamt.ai.domain_gap.simulator.apply_dropout` marks whole
stations, whole frequencies, or individual observations missing by
setting impedance to ``NaN``, which
:class:`~pycsamt.ai.data.contracts.SurveyData` construction then
turns into the authoritative validity mask automatically — no separate
bookkeeping needed. :func:`~pycsamt.ai.domain_gap.simulator.inject_outliers`
does the opposite of dropout: it perturbs a fraction of otherwise
*valid* observations by a large factor without invalidating them,
simulating a bad reading that quality control failed to flag — exactly
the case a robust downstream inverter has to tolerate, not just
detect.

Composing one reproducible pass
------------------------------------

A real survey was not corrupted by one effect at a time, and neither
should a training example be.
:func:`~pycsamt.ai.domain_gap.simulator.apply_corruption_suite` runs
every step above in one fixed, physically motivated order — static
shift and distortion first, then noise, error floor, dropout, and
outliers — from a single seed that spawns an independent child
generator per step, so adding a new step later does not change what
every earlier step draws. It returns both the corrupted survey and a
:class:`~pycsamt.ai.domain_gap.simulator.CorruptionRecord` of exactly
what was sampled, not just what was requested. Four named severities
in :data:`~pycsamt.ai.domain_gap.simulator.SEVERITY_PRESETS` turn a
whole parameter block into one word:

.. list-table::
   :header-rows: 1
   :widths: 22 20 16 21 21

   * - Preset
     - Noise level
     - Static shift :math:`\sigma_{\log_{10}}`
     - Station dropout
     - Outlier rate
   * - ``clean``
     - 0
     - 0
     - 0
     - 0
   * - ``in_distribution``
     - 1-3%
     - 0.05
     - 0
     - 0
   * - ``severe``
     - 5-15%
     - 0.15
     - 5%
     - 2%
   * - ``held_out_corruption``
     - 10-25%
     - 0.30
     - 15%
     - 8%

``clean`` is the required control set every :doc:`losses`/
:doc:`scientific_validation` comparison needs a baseline from;
``held_out_corruption`` is deliberately outside the range a model
should have been trained to expect, for out-of-distribution checks.

Fitting corruption from a real survey
-------------------------------------------

Guessed corruption ranges are still a guess. The more useful path is
fitting them from a real survey's own quality-control diagnostics.
:func:`~pycsamt.ai.domain_gap.survey_fit.survey_data_from_sites`
bridges anything :func:`~pycsamt.emtools._core.ensure_sites` accepts
into a canonical :class:`~pycsamt.ai.data.contracts.SurveyData`,
converting impedance from the EDI-native ``[mV/km]/[nT]`` convention
:class:`~pycsamt.site.base.Site.z` stores to the SI convention
:class:`SurveyData` declares — a real unit mismatch found and fixed
while building this bridge, confirmed against
:attr:`~pycsamt.site.base.Site.rho`'s own field-convention apparent
resistivity to floating-point precision. From there,
:func:`~pycsamt.ai.domain_gap.survey_fit.fit_corruption_config` reads
noise and dropout directly from the survey's declared
``impedance_error`` and coverage, while
:func:`~pycsamt.ai.domain_gap.survey_fit.fit_distortion_priors_from_sites`
calls the existing :func:`pycsamt.emtools.gb.groom_bailey_table` and
:func:`pycsamt.emtools.ss.estimate_ss_ama` diagnostics directly on the
sites to get an empirical static-shift and distortion spread — the
one part of this package that genuinely depends on a specific survey
rather than on literature defaults:

.. code-block:: pycon

   >>> from pycsamt.emtools._core import ensure_sites
   >>> from pycsamt.ai.domain_gap import (
   ...     survey_data_from_sites,
   ...     fit_corruption_config,
   ...     fit_distortion_priors_from_sites,
   ... )

   >>> sites = ensure_sites(
   ...     "data/AMT/WILLY_data/L18PLT", recursive=True, verbose=0
   ... )
   >>> field = survey_data_from_sites(sites, recursive=False, verbose=0)
   >>> field.shape
   (28, 53, 4)

   >>> config = fit_corruption_config(field)
   >>> tuple(round(v, 4) for v in config.noise_level_range)
   (0.0509, 0.1582)
   >>> round(config.error_floor_fraction, 4)
   0.0297

The fitted noise range (5-16% relative) and error floor (about 3% of
:math:`|Z|`) come entirely from L18's own declared errors, not from a
default — a survey with tighter QC would fit a narrower range without
anyone having to retune a preset by hand.

Comparing simulated and field distributions
--------------------------------------------------

Fitting a corruption config is not, by itself, evidence that it
closes the domain gap. :mod:`~pycsamt.ai.domain_gap.report` computes
that evidence directly: per-feature summary statistics and a
two-sample Kolmogorov-Smirnov test between a simulated and a field
survey, on ``log_impedance_magnitude``, ``phase_deg``, and
``error_to_magnitude_ratio`` — features derived from the shared
:class:`SurveyData` contract so the comparison never depends on how
either survey was produced. The real M3 workflow applies a
field-fitted config to a *clean synthetic* survey, not back onto the
field survey itself, so the comparison actually tests whether a
generic geological prior looks like the target survey:

.. code-block:: pycon

   >>> from pycsamt.ai.geology import GeologyGrid
   >>> from pycsamt.ai.training.dataset2d import (
   ...     Maxwell2DDatasetConfig,
   ...     generate_2d_maxwell_dataset,
   ... )
   >>> from pycsamt.ai.domain_gap import (
   ...     apply_corruption_suite,
   ...     compare_survey_distributions,
   ... )

   >>> grid = GeologyGrid.regular_2d(nx=10, nz=8, dx_m=300, dz_m=150)
   >>> ds_config = Maxwell2DDatasetConfig(
   ...     dataset_id="willy-like-clean",
   ...     grid=grid,
   ...     correlation_length_x_m=(600.0, 1200.0),
   ...     correlation_length_z_m=(200.0, 400.0),
   ...     frequencies_hz=np.logspace(0, 4, 20),
   ...     station_x_m=list(np.linspace(500.0, 2500.0, 10)),
   ...     n_realizations=1,
   ...     seed=1,
   ...     log_resistivity_mean=2.0,
   ...     log_resistivity_std=0.4,
   ...     validation_fraction=0.0,
   ...     test_fraction=0.0,
   ... )
   >>> dataset = generate_2d_maxwell_dataset(ds_config)
   >>> clean_synthetic = dataset.samples[0].survey

   >>> corrupted, record = apply_corruption_suite(
   ...     clean_synthetic, config, seed=3
   ... )
   >>> report = compare_survey_distributions(corrupted, field)
   >>> magnitude = report.comparisons["log_impedance_magnitude"]
   >>> round(magnitude.ks_statistic, 4), round(magnitude.mean_difference, 4)
   (0.1376, 0.0445)
   >>> phase = report.comparisons["phase_deg"]
   >>> round(phase.ks_statistic, 4), round(phase.mean_difference, 2)
   (0.2908, -13.13)

The near-zero magnitude mean difference confirms the corruption
config's noise and error floor are pulling the synthetic magnitude
distribution onto the field survey's actual scale — the thing
:func:`fit_corruption_config` was built to do. The phase comparison
is the more honest result: a KS statistic of 0.29 and a 13-degree mean
offset show that a generic correlated-field prior with only
kilometre-scale correlation lengths does not reproduce L18's real
phase behaviour, which reflects genuine 2-D and 3-D structure this
small demonstration grid was never built to represent. That gap is
exactly what :doc:`geology_priors` and :doc:`../forward/index`'s
survey-specific tuning — not the noise model — would need to close;
:mod:`~pycsamt.ai.domain_gap.report` is what makes that distinction
visible instead of hiding it behind one aggregate "looks about right."

Accounting for every station, not just corruption
--------------------------------------------------------

A different, related question — is *this* survey internally coherent
enough to trust at all, station by station — is
:func:`~pycsamt.ai.domain_gap.audit.audit_survey`'s job, covered in
full in :ref:`ai_inversion_data_contracts_audit`. It lives in this
package for the same reason the functions above do: it depends on the
pandas-based dimensionality and distortion diagnostics
:mod:`pycsamt.ai.data` deliberately excludes.
