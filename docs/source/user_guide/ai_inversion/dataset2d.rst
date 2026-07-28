.. _ai_inversion_dataset2d:

2-D Maxwell training-data generation
=====================================

:doc:`forward_physics` solves one 2-D Maxwell problem at a time.
Training a network needs many, spanning a distribution of plausible
geology, converted into the same shape every later stage agrees on.
:func:`~pycsamt.ai.training.dataset2d.generate_2d_maxwell_dataset`
is that bridge: it draws a correlated resistivity field from
:mod:`~pycsamt.ai.geology` for each realization, solves it with
:class:`~pycsamt.forward.maxwell.mt2d.MT2DAdapter` through
:func:`~pycsamt.forward.maxwell.batch.solve_batch` (cached and
retried rather than resolved from scratch on every call), and
packages every realization that converges into a
:class:`~pycsamt.ai.training.dataset2d.Maxwell2DDataset` with a
realization-level train/validation/test split and a
:class:`~pycsamt.ai.data.manifest.DatasetManifest`. A realization
whose solve does not converge is recorded and excluded, never
silently included.

Configuring one coherent dataset
------------------------------------

A :class:`~pycsamt.ai.training.dataset2d.Maxwell2DDatasetConfig` fixes
everything that must stay identical across realizations for the
dataset to be one coherent, versioned thing: the shared ``(z, x)``
:class:`~pycsamt.ai.geology.GeologyGrid`, the horizontal and vertical
correlation-length ranges realizations are drawn from, the shared
frequency grid and receiver ``x`` positions, and a root
:term:`Random seed` from which every realization's field and solver
seeds are derived deterministically:

.. code-block:: pycon

   >>> from pycsamt.ai.geology import GeologyGrid
   >>> from pycsamt.ai.training.dataset2d import (
   ...     Maxwell2DDatasetConfig,
   ...     generate_2d_maxwell_dataset,
   ... )

   >>> grid = GeologyGrid.regular_2d(nx=6, nz=4, dx_m=300, dz_m=150)
   >>> config = Maxwell2DDatasetConfig(
   ...     dataset_id="demo-2d-v1",
   ...     grid=grid,
   ...     correlation_length_x_m=(600.0, 900.0),
   ...     correlation_length_z_m=(150.0, 300.0),
   ...     frequencies_hz=[10.0, 3.0],
   ...     station_x_m=[600.0, 900.0, 1200.0],
   ...     n_realizations=2,
   ...     seed=0,
   ...     validation_fraction=0.0,
   ...     test_fraction=0.0,
   ... )
   >>> config.components
   ('zxy', 'zyx')

   >>> dataset = generate_2d_maxwell_dataset(config)
   >>> len(dataset.samples), dataset.rejected
   (2, ())

Every accepted realization comes back as a
:class:`~pycsamt.ai.training.dataset2d.Maxwell2DSample` pairing the
true resistivity grid (the training target) with its simulated
:class:`~pycsamt.ai.data.contracts.SurveyData` response (the training
input) and the solver's own mesh size and residual, so a downstream
consumer never has to guess how well that one realization actually
converged:

.. code-block:: pycon

   >>> sample = dataset.samples[0]
   >>> sample.survey.shape
   (3, 2, 2)
   >>> sample.resistivity_ohm_m.shape
   (4, 6)
   >>> sample.mesh_cells, round(sample.relative_residual, 6)
   (11492, 0.0)

Two realizations, three stations, two frequencies, and both impedance
components finish in under a second here because the mesh and
realization count are kept deliberately small for the documentation
build; the mesh alone uses 11,492 padded cells even at this scale,
which is why :doc:`forward_physics`'s cache and batch runner, not a
plain loop, are what make a production-sized version of this call
tractable. The ``relative_residual`` of ``0.0`` and the empty
``rejected`` tuple are the honest signal to check before trusting a
generated dataset at all — a nonzero, unexcluded residual would mean
a realization entered training on an unconverged solve.

.. figure:: ../../images/user_guide/ai_inversion/dataset2d_realization_gallery.png
   :alt: Three independent resistivity-field realizations from one Maxwell2DDatasetConfig
   :align: center
   :width: 96%

   Three realizations drawn from the same
   :class:`Maxwell2DDatasetConfig` — same correlation-length ranges,
   same station layout (white triangles), independent random fields.
   Each is a genuinely 2-D section with structure varying in both
   ``x`` and ``z``, not a set of independent 1-D columns tiled
   side by side; that lateral continuity is the entire reason this
   module exists rather than reusing :doc:`../forward/index`'s
   station-wise 1-D batch generator.

Two impedance modes, and why one needed more care than the other
------------------------------------------------------------------------

A 2-D MT problem separates into two independent modes, and they are
not equally easy to discretize. With :math:`\mu=\mu_0` constant
everywhere, :term:`TE mode` solves for the electric field
:math:`E_y` under

.. math::
   :label: eq-ai-dataset2d-te

   \nabla^2 E_y = i\omega\mu_0\sigma\,E_y,

a *constant-coefficient* Laplacian: resistivity structure enters only
through the reaction term on the right, multiplying the field itself
rather than its derivatives. :term:`TM mode` solves for the magnetic
field :math:`H_y` under

.. math::
   :label: eq-ai-dataset2d-tm

   \nabla\cdot(\rho\,\nabla H_y) = i\omega\mu_0 H_y,

where resistivity sits *inside* the derivative — the diffusion
coefficient itself, not just a multiplier. Equation
:eq:`eq-ai-dataset2d-tm` is a genuinely harder numerical problem than
:eq:`eq-ai-dataset2d-te`: every bit of 2-D resistivity structure a
finite-difference scheme needs to resolve has to survive being
differentiated through, not just multiplied in afterward. In
practice this means TM mode needs finer *joint* horizontal and
vertical resolution than TE mode to reach comparable accuracy — and
that requirement, not a defect in the solver, was the real story
behind an earlier, mistaken finding that TM mode was unusable for
laterally fine structure. That investigation changed horizontal and
vertical resolution one at a time instead of together, and compared
independent random-field realizations across different
resolutions instead of refining one fixed field — both of which
produce exactly the appearance of a solver that never converges, even
when it does. Once resolution is refined in both directions together
on one fixed field, TM mode converges cleanly, changing by under 1-2%
between realistic production resolution and twice that.

A second, real but smaller issue was fixed alongside that finding:
:mod:`pycsamt.forward.em2d`'s TM-mode assembly combined resistivity
values at cell interfaces with a plain arithmetic mean, which is the
correct combination for two cells in series along the direction being
differentiated but not for the two cells actually involved here,
which act as parallel current paths. The physically correct
combination — a thickness-weighted harmonic mean of resistivity — is
what the surface-response extraction already used; the interior
assembly now matches it. Both impedance modes are therefore validated
and requested by default, and a near-uniform model confirms both
against the analytic half-space limit at once:

.. code-block:: pycon

   >>> import numpy as np
   >>> from pycsamt.forward.maxwell.benchmarks import half_space_impedance

   >>> grid = GeologyGrid.regular_2d(nx=4, nz=4, dx_m=200, dz_m=100)
   >>> config = Maxwell2DDatasetConfig(
   ...     dataset_id="halfspace-check",
   ...     grid=grid,
   ...     correlation_length_x_m=(1.0, 2.0),
   ...     correlation_length_z_m=(1.0, 2.0),
   ...     frequencies_hz=[10.0, 1.0],
   ...     station_x_m=[400.0],
   ...     n_realizations=1,
   ...     seed=0,
   ...     log_resistivity_mean=2.0,
   ...     log_resistivity_std=1e-6,
   ...     validation_fraction=0.0,
   ...     test_fraction=0.0,
   ... )
   >>> dataset = generate_2d_maxwell_dataset(config)
   >>> sample = dataset.samples[0]
   >>> analytic = half_space_impedance(100.0, config.frequencies_hz)
   >>> zxy = sample.survey.impedance[0, :, 0]
   >>> zyx = sample.survey.impedance[0, :, 1]
   >>> te_error = np.abs(zxy - analytic) / np.abs(analytic)
   >>> round(float(np.max(te_error)), 4)
   0.0424

   >>> # Zyx = -Zxy for an isotropic half-space (standard MT convention).
   >>> tm_error = np.abs(zyx - (-analytic)) / np.abs(analytic)
   >>> round(float(np.max(tm_error)), 4)
   0.0411

Both modes agree with the analytic reference to about 4%, using this
module's own default mesh construction rather than a specially tuned
one — the same order of accuracy for both, not TE working and TM
silently failing next to it.

Handing this to noise, or to training
-----------------------------------------

This module deliberately stops at a clean, forward-consistent
response. Turning it into training data that resembles a specific
field survey — adding heteroscedastic noise, dropout, static shift,
and distortion fitted from that survey's own quality control — is
:doc:`domain_gap`'s job, kept separate so that a clean control dataset
always exists independently of whatever noise model is layered on
top of it. From there, a :class:`Maxwell2DDataset`'s ``train``/
``validation``/``test`` partitions and its
:class:`~pycsamt.ai.data.manifest.DatasetManifest` are exactly the
inputs :doc:`losses` trains against and :doc:`scientific_validation`
reports evidence from — the realization-level split guarding against
the same leakage :doc:`data_contracts` already enforces for any other
kind of split.
