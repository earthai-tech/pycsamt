.. _ai_inversion_pinn_2d:

Physics-informed 2-D inversion
==============================

:class:`pycsamt.ai.inversion.PINNInverter2D` is a
:term:`physics-informed inversion` workflow that jointly optimizes layered
models for all stations on a profile. It uses a differentiable 1-D MT
:term:`forward operator` at each station and couples adjacent station models
with lateral smoothness. The result is a pseudo-2-D resistivity section: a set
of jointly regularized 1-D models, not a full finite-element or
finite-difference 2-D EM forward inversion.

This distinction is central to scientific use. The workflow can provide a
useful profile estimate when locally layered physics and lateral smoothness are
reasonable approximations. It should not be presented as equivalent to
Occam2D, MARE2DEM, ModEM, or another solver that explicitly models 2-D/3-D EM
fields.

.. admonition:: Meaning of “2-D”
   :class: important

   “2-D” here describes the jointly optimized station–layer section and its
   lateral regularization. Data fit is calculated with the differentiable
   Wait-style MT 1-D recursion independently at each station. Strong lateral
   induction, off-profile structure, anisotropy, topography, and other missing
   physics can therefore produce a smooth but biased section.

When to use this workflow
-------------------------

Consider PINN 2-D when:

* stations form an ordered survey profile;
* :term:`dimensionality` evidence suggests locally 1-D or gently varying 2-D
  structure, with :term:`phase tensor` and strike diagnostics reviewed;
* a common layered parameterization is scientifically acceptable as a
  :term:`model prior`;
* no labelled model dataset is available or desired;
* direct response-space optimization is preferred to supervised prediction;
* a fast pseudo-2-D baseline or starting model is useful;
* the result will be compared with classical inversion,
  :term:`response reconstruction`, and independent data.

Prefer a full 2-D or 3-D classical solver when:

* lateral induction is important to the target;
* strong contacts, topography, coast effects, anisotropy, or 3-D structure are
  expected;
* tensor components cannot be represented by per-station 1-D responses;
* model appraisal requires the corresponding numerical forward physics;
* regulatory or project standards require a recognized 2-D/3-D solver.

Workflow
--------

#. QC and order the profile data;
#. review strike and dimensionality evidence;
#. choose TE, TM, or the current averaged ``both`` mode;
#. define the common frequency grid and supported depth;
#. choose layer count, smoothness weights, optimizer, and epochs;
#. run baseline and sensitivity configurations;
#. inspect convergence and :term:`response-space metric` residuals;
#. compare section, thicknesses, and reconstructed responses across scenarios;
#. validate against classical inversion and independent evidence;
#. report the pseudo-2-D physics and all limitations explicitly.

1. Understand the model parameterization
----------------------------------------

For ``S`` stations and ``L`` layers, the optimizer uses:

.. code-block:: text

   log_rho.shape   = (S, L)
   log_thick.shape = (S, L - 1)

Both arrays are optimized jointly. Resistivity is represented as
:math:`\log_{10}(\rho/\Omega\mathrm{m})`; thickness is represented internally
in log10 metres and converted back for output.

The optimized station model is therefore

.. math::

   \mathbf{m}_s =
   \left[
      u_{s1},\ldots,u_{sL},
      q_{s1},\ldots,q_{s,L-1}
   \right],

where :math:`u_{s\ell}=\log_{10}\rho_{s\ell}` and
:math:`q_{s\ell}=\log_{10}h_{s\ell}`. The pseudo-2-D section is the matrix
:math:`U=[u_{s\ell}]`, but the forward response at station :math:`s` still
depends only on :math:`\mathbf{m}_s`. Lateral coupling enters through
regularization, not through a 2-D electromagnetic field calculation.

The public section methods transpose the internal station-major layout:

.. code-block:: text

   resistivity_section.shape = (L, S)
   thickness_section.shape   = (L - 1, S)

Each station has its own interface thicknesses. Consequently, layer index is a
parameter correspondence, not automatically a continuous geological horizon.
Inspect cumulative interface depths before drawing boundaries between
stations.

For each station, cumulative interface depth is

.. math::

   z_{sk}=\sum_{\ell=1}^{k}10^{q_{s\ell}}.

If thickness varies laterally, the same layer index can sit at different
physical depths from station to station. Plotting only layer index is useful
for debugging, but it is not a physical-depth interpretation.

``depth_max`` guides the layered depth parameterization and initialization; it
does not establish a geophysical depth of investigation. The defensible
interpretation depth must come from bandwidth, sensitivity, response fit, and
scenario stability.

2. Understand the loss function
-------------------------------

The implemented objective combines data fit, vertical smoothing, and lateral
smoothing:

.. math::

   \mathcal L = \mathcal L_{data}
   + \lambda_z\mathcal R_z
   + \lambda_x\mathcal R_x.

Data term
~~~~~~~~~

For valid station–frequency cells, the data term combines squared differences
in log10 apparent resistivity and normalized phase:

.. math::

   \mathcal L_{data} = \frac{1}{N_v}\sum_{s,f}
   \left[
   (\log_{10}\rho_{a,sf}^{pred}-\log_{10}\rho_{a,sf}^{obs})^2
   + \left(\frac{\phi_{sf}^{pred}-\phi_{sf}^{obs}}{90}\right)^2
   \right].

``N_v`` is the number of finite observations after common-grid interpolation.
Missing cells are excluded from the data term.

Equivalently, with finite-data mask :math:`M_{sf}`,

.. math::

   N_v = \sum_{s,f} M_{sf}.

A station with fewer valid frequencies contributes fewer residual terms. It
may still look smooth in the final section because lateral regularization
borrows structure from neighboring stations. Inspect valid-cell counts before
interpreting weakly constrained stations.

The 90-degree phase normalization makes the two residual blocks numerically
comparable but is not an observational covariance model. Formal errors and
station-specific uncertainty are not explicit weights in this objective.

If project errors are available, compare the final section with an external
normalized :term:`RMS misfit` calculation in :doc:`validation`. The optimizer
loss and a formal error-weighted inversion RMS are related diagnostics, not the
same statistic.

Vertical regularization
~~~~~~~~~~~~~~~~~~~~~~~

``smoothness_weight`` is :math:`\lambda_z`. It penalizes changes in log
resistivity between adjacent layers at each station. Larger values favor
vertically smooth models; smaller values permit sharper changes and greater
instability.

In station-layer notation, a representative vertical penalty is

.. math::

   \mathcal R_z =
   \sum_{s=1}^{S}\sum_{\ell=1}^{L-1}
   \left(u_{s,\ell+1}-u_{s,\ell}\right)^2.

Lateral regularization
~~~~~~~~~~~~~~~~~~~~~~

``lateral_weight`` is :math:`\lambda_x`. It penalizes log-resistivity
differences between adjacent stations at the same layer index. It assumes the
input station order represents profile adjacency.

The corresponding chain penalty is

.. math::

   \mathcal R_x =
   \sum_{s=1}^{S-1}\sum_{\ell=1}^{L}
   \left(u_{s+1,\ell}-u_{s,\ell}\right)^2.

The lateral term does not currently scale differences by physical station
spacing. Closely and widely spaced station pairs are treated as adjacent
indices. Resample carefully or interpret regularization strength in light of
irregular chainage.

Thickness regularization and constraints
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Resistivity smoothness is explicit in the documented objective. Layer
thicknesses are optimized in log space, but the same lateral/vertical
regularization description should not be assumed to apply to them unless
confirmed in the backend kernel. Audit thickness stability separately.

3. Prepare and order profile data
---------------------------------

Load the profile canonically and verify station order before constructing the
inverter:

.. code-block:: pycon

   >>> from pycsamt.emtools import ensure_sites
   >>>
   >>> sites = ensure_sites(  # doctest: +SKIP
   ...     "data/AMT/WILLY_DATA/L18",
   ...     recursive=True,
   ...     verbose=0,
   ... )

The inverter uses :func:`pycsamt.ai.inversion.sites_to_obs_2d` internally and
preserves the order returned by the input. It does not sort stations by
chainage. Provide a :class:`~pycsamt.site.Sites` object already ordered along
the intended profile direction.

Retain a table containing:

* station name;
* profile distance in metres;
* easting, northing, and elevation;
* original and inversion order;
* valid frequency range;
* TE/TM component completeness;
* QC exclusions and corrections.

Do not infer physical x coordinates from array column number in a final plot.

4. Select polarization mode
---------------------------

``mode`` accepts ``"te"``, ``"tm"``, or ``"both"``.

``te``
   Uses the component named by ``comp_te``, default ``"xy"``.

``tm``
   Uses the component named by ``comp_tm``, default ``"yx"``. The public
   observation bridge stores TM phase as magnitude.

``both``
   In the current implementation, averages TE and TM apparent resistivity and
   averages TE and TM phase before optimization. It does not form two
   independent data-misfit blocks or jointly fit both modes separately.

.. warning::

   Averaging TE and TM can suppress genuine mode differences caused by 2-D/3-D
   structure. Use ``mode="both"`` only as an explicitly reviewed approximation.
   For scientific comparison, run TE and TM separately and inspect their
   disagreement.

Mode labels depend on profile orientation and strike convention. Verify these
before accepting the defaults.

5. Choose the common frequency grid
-----------------------------------

The constructor extracts valid observations and builds a log-spaced common
grid with ``n_freqs`` points spanning the combined available range:

.. code-block:: pycon

   >>> from pycsamt.ai.inversion import PINNInverter2D
   >>>
   >>> inversion = PINNInverter2D(  # doctest: +SKIP
   ...     sites,
   ...     n_layers=10,
   ...     depth_max=2000.0,
   ...     n_freqs=32,
   ...     mode="te",
   ... )

Each station is interpolated to this grid in log-frequency space. Values
outside a station's own range become NaN and do not contribute to the data
loss.

Because the public constructor does not expose explicit ``freq_min`` and
``freq_max`` arguments, the common endpoints follow the extracted survey
coverage. To enforce a project-specific range, prefilter the reviewed data or
extend the API through a tested project change rather than mutating private
attributes.

Inspect coverage per station. A station constrained by only a small subset of
the common grid can borrow visual smoothness from neighbors without having the
same depth sensitivity.

6. Choose layer count and depth
-------------------------------

Start with the simplest parameterization able to represent expected electrical
units. Test several defensible combinations, for example:

.. code-block:: pycon

   >>> candidates = [
   ...     {"n_layers": 6, "depth_max": 1500.0},
   ...     {"n_layers": 10, "depth_max": 2000.0},
   ...     {"n_layers": 14, "depth_max": 3000.0},
   ... ]
   >>> len(candidates)
   3
   >>> candidates[1]
   {'n_layers': 10, 'depth_max': 2000.0}

Increasing layers can reduce response residuals while adding unstable
interfaces. Increasing ``depth_max`` can allocate parameters below meaningful
sensitivity. Compare:

* reconstructed response fit;
* interface-depth stability;
* deep-layer variation across runs;
* dependence on smoothing;
* classical inversion and sensitivity depth;
* independent borehole or geological evidence.

Do not interpret the bottom of the configured model as the investigation
depth.

7. Select regularization weights
--------------------------------

The defaults are:

.. code-block:: text

   smoothness_weight = 0.01
   lateral_weight    = 0.005

Treat them as starting values, not universal choices. Run a grid of plausible
weights:

.. code-block:: pycon

   >>> scenarios = [
   ...     (0.001, 0.0005),
   ...     (0.01, 0.005),
   ...     (0.05, 0.02),
   ... ]
   >>> len(scenarios)
   3
   >>> scenarios[0]
   (0.001, 0.0005)

For each scenario, retain total loss, response residuals, roughness measures,
section differences, and target stability. The reported convergence curve
contains total loss, not separate data and regularization histories, so
additional diagnostics may require calculation from outputs or backend
instrumentation.

Avoid choosing weights solely for a visually smooth section. Excessive lateral
smoothing can smear faults, isolated conductors, or station-specific problems.
Insufficient smoothing can fit noise with alternating layers.

8. Select optimizer controls
----------------------------

``epochs``
   Maximum Adam iterations, default 300.

``lr``
   Adam learning rate, default ``1e-2``.

``device``
   PyTorch/TensorFlow device resolved through the active backend.

The current high-level ``fit()`` does not expose early stopping, validation
data, learning-rate scheduling, or gradient clipping. It runs the requested
optimization iterations through the backend kernel. Monitor the curve and run
repeat configurations rather than assuming the final epoch is optimal.

A large learning rate can oscillate or diverge; a very small rate may appear
stable without reaching an adequate fit. Compare final sections across epochs
and learning rates.

9. Run a baseline inversion
----------------------------

.. code-block:: pycon

   >>> from pycsamt.ai.inversion import PINNInverter2D
   >>>
   >>> inversion = PINNInverter2D(  # doctest: +SKIP
   ...     sites,
   ...     n_layers=10,
   ...     depth_max=2000.0,
   ...     n_freqs=32,
   ...     mode="te",
   ...     smoothness_weight=0.01,
   ...     lateral_weight=0.005,
   ...     epochs=300,
   ...     lr=1e-2,
   ...     comp_te="xy",
   ...     comp_tm="yx",
   ...     device=None,
   ...     verbose=0,
   ... )
   >>> inversion.fit(verbose=True, log_every=25)  # doctest: +SKIP
   >>> print(inversion)  # doctest: +SKIP
   >>> print(inversion.stations)  # doctest: +SKIP
   >>> print(inversion.n_sites)  # doctest: +SKIP

The constructor performs data extraction immediately and can raise before
``fit()`` when no valid observations survive. The deep-learning backend is
required when fitting.

10. Extract resistivity and thickness
-------------------------------------

.. code-block:: pycon

   >>> import numpy as np
   >>> log10_rho = np.ones((4, 3)) * 2.0
   >>> rho_ohm_m = 10.0 ** log10_rho
   >>> thickness_m = np.array([
   ...     [50.0, 60.0, 55.0],
   ...     [100.0, 120.0, 110.0],
   ...     [200.0, 240.0, 220.0],
   ... ])
   >>> print(log10_rho.shape)
   (4, 3)
   >>> print(thickness_m.shape)
   (3, 3)

For a fitted inverter, use the public section methods:

.. code-block:: pycon

   >>> log10_rho = inversion.resistivity_section(as_log10=True)  # doctest: +SKIP
   >>> rho_ohm_m = inversion.resistivity_section(as_log10=False)  # doctest: +SKIP
   >>> thickness_m = inversion.thickness_section()  # doctest: +SKIP

Build interface depths per station from cumulative thickness:

.. code-block:: pycon

   >>> interface_depth_m = np.cumsum(thickness_m, axis=0)
   >>> print(interface_depth_m.shape)
   (3, 3)
   >>> print(interface_depth_m[:, 0].astype(int).tolist())
   [50, 150, 350]
   >>> assert np.all(np.isfinite(rho_ohm_m))
   >>> assert np.all(rho_ohm_m > 0)
   >>> assert np.all(np.isfinite(thickness_m))
   >>> assert np.all(thickness_m > 0)

Check predictions against explicit parameter bounds used in project scenarios.
The optimizer's positivity through log representation does not guarantee
geological plausibility.

11. Review convergence
----------------------

.. code-block:: pycon

   >>> curve = inversion.convergence_curve()  # doctest: +SKIP
   >>> print(curve.tail())  # doctest: +SKIP

The returned pandas DataFrame contains ``epoch`` and total ``loss``.

Review:

* initial reduction rate;
* late oscillation or divergence;
* whether loss is still decreasing strongly at the final epoch;
* sensitivity to learning rate;
* repeatability across backend/device and repeated runs;
* whether lower loss corresponds to better response residuals and stable
  structure.

There is no separate validation loss because the optimization is
observation-specific rather than supervised dataset training. Overfitting can
still occur by fitting noise with an overly flexible or weakly regularized
model.

12. Review residuals
--------------------

.. code-block:: pycon

   >>> residuals = inversion.residuals()  # doctest: +SKIP
   >>> residuals["log_rho_residual"] = (  # doctest: +SKIP
   ...     np.log10(residuals["rho_pred"])
   ...     - np.log10(residuals["rho_obs"])
   ... )
   >>> residuals["phase_residual_deg"] = (  # doctest: +SKIP
   ...     residuals["phase_pred"] - residuals["phase_obs"]
   ... )
   >>> summary = residuals.groupby("station").agg(  # doctest: +SKIP
   ...     rho_rms=("log_rho_residual",
   ...              lambda values: float(np.sqrt(np.nanmean(values**2)))),
   ...     phase_rms_deg=("phase_residual_deg",
   ...                    lambda values: float(np.sqrt(np.nanmean(values**2)))),
   ... )
   >>> print(summary)  # doctest: +SKIP

The residual table uses each station's original frequencies, not only the
common optimization grid.

.. warning::

   For ``mode="both"``, the optimizer fits averaged TE/TM observations, while
   ``residuals()`` selects TE observations for its reported observed columns.
   The diagnostic therefore does not directly report residuals against the
   averaged data used in optimization. Run TE and TM separately or construct a
   custom audited residual table for the averaged mode.

Also verify the shapes/components returned by the forward response in the
installed version. A single aggregate RMS can hide station-, frequency-, or
mode-dependent mismatch.

13. Plot the section correctly
------------------------------

The section is defined by varying thicknesses, so ``imshow`` on layer index can
misrepresent depth. For a quick layer-index diagnostic:

.. code-block:: pycon

   >>> import matplotlib.pyplot as plt  # doctest: +SKIP
   >>>
   >>> fig, ax = plt.subplots(figsize=(10, 5))  # doctest: +SKIP
   >>> image = ax.imshow(  # doctest: +SKIP
   ...     log10_rho,
   ...     aspect="auto",
   ...     origin="upper",
   ...     interpolation="nearest",
   ... )
   >>> ax.set_xlabel("Station index")  # doctest: +SKIP
   >>> ax.set_ylabel("Layer index")  # doctest: +SKIP
   >>> fig.colorbar(image, ax=ax, label="log10 resistivity (ohm m)")  # doctest: +SKIP
   >>> fig.savefig(  # doctest: +SKIP
   ...     "review/pinn2d_layer_index.png",
   ...     dpi=200,
   ...     bbox_inches="tight",
   ... )

Label this explicitly as layer index. For a physical-depth section, resample
each layered model onto a common depth grid using its cumulative thicknesses,
then plot against reviewed profile distances. Preserve the original layered
outputs alongside any resampled image.

Do not connect layer boundaries as geological horizons without checking their
meaning and stability.

14. Run TE/TM and regularization scenarios
------------------------------------------

A defensible study includes separate modes:

.. code-block:: pycon

   >>> results = {}  # doctest: +SKIP
   >>> for mode in ("te", "tm"):  # doctest: +SKIP
   ...     inv = PINNInverter2D(
   ...         sites,
   ...         n_layers=10,
   ...         depth_max=2000.0,
   ...         n_freqs=32,
   ...         mode=mode,
   ...         smoothness_weight=0.01,
   ...         lateral_weight=0.005,
   ...         epochs=300,
   ...         lr=1e-2,
   ...     ).fit(verbose=False)
   ...     results[mode] = {
   ...         "log10_rho": inv.resistivity_section(),
   ...         "thickness_m": inv.thickness_section(),
   ...         "residuals": inv.residuals(),
   ...     }

Compare:

* structures common to TE and TM;
* mode-specific conductive or resistive features;
* station and frequency residual patterns;
* boundary-depth differences;
* sensitivity to :math:`\lambda_z` and :math:`\lambda_x`;
* changes when suspicious stations or bands are excluded.

Large disagreement may indicate dimensionality, distortion, component quality,
or parameterization problems rather than a need to average modes.

15. Compare with 1-D and classical 2-D results
----------------------------------------------

The lateral weight should add demonstrable value over independent station
models. Compare against:

* :class:`pycsamt.ai.inversion.PINNInverter1D` at each station;
* a no/very-low lateral-weight scenario;
* a conventional layered 1-D inversion;
* Occam2D when a 2-D line assumption is appropriate;
* ModEM or MARE2DEM when geometry and resources justify them;
* boreholes and mapped geology.

Compare response fit using compatible error models. The PINN objective's phase
normalization and missing-value mask need not match the classical solver's RMS.
Do not rank methods by incomparable scalar values.

16. Use the PINN agent for orchestration
----------------------------------------

:class:`pycsamt.agents.PINNInversionAgent` wraps the workflow and returns an
:class:`pycsamt.agents.AgentResult`:

.. code-block:: pycon

   >>> from pycsamt.agents import PINNInversionAgent
   >>>
   >>> agent = PINNInversionAgent(  # doctest: +SKIP
   ...     dim=2,
   ...     n_layers=10,
   ...     depth_max=2000.0,
   ...     smoothness_weight=0.01,
   ...     lateral_weight=0.005,
   ...     epochs=300,
   ...     lr=1e-2,
   ...     solver="mt1d",
   ... )
   >>> result = agent.execute({  # doctest: +SKIP
   ...     "sites": sites,
   ...     "output_dir": "outputs/pinn2d/L18",
   ... })
   >>> if result.status == "failed":  # doctest: +SKIP
   ...     raise RuntimeError(result.error)
   >>> section = result["section"]  # doctest: +SKIP
   >>> loss_table = result.get("loss_df")  # doctest: +SKIP
   >>> residual_table = result.get("residuals_df")  # doctest: +SKIP

The agent is convenient for figures and standardized results, but inspect its
constructor mapping: the generic agent interface does not expose every
``PINNInverter2D`` option such as mode and component selection. Use the lower-
level class when those controls matter.

17. Consider hybrid 2-D refinement
----------------------------------

:class:`pycsamt.ai.inversion.HybridInverter2D` uses a fitted supervised
:class:`pycsamt.ai.inversion.EMInverter2D` to initialize the same joint physics
refinement:

.. code-block:: pycon

   >>> from pycsamt.ai.inversion import HybridInverter2D
   >>>
   >>> hybrid = HybridInverter2D(  # doctest: +SKIP
   ...     sites,
   ...     ai_inverter=fitted_ai_2d,
   ...     n_layers=10,
   ...     depth_max=2000.0,
   ...     n_freqs=32,
   ...     mode="te",
   ...     smoothness_weight=0.005,
   ...     lateral_weight=0.003,
   ...     epochs=150,
   ...     lr=5e-3,
   ... ).fit()
   >>> stage1 = hybrid.stage1_section()  # doctest: +SKIP
   >>> stage2 = hybrid.resistivity_section()  # doctest: +SKIP
   >>> residuals_stage1 = hybrid.residuals(stage=1)  # doctest: +SKIP
   >>> residuals_stage2 = hybrid.residuals(stage=2)  # doctest: +SKIP

The hybrid is justified only if Stage 2 improves response fit or scientific
stability without introducing unsupported structure. Preserve both stages.

.. warning::

   ``HybridInverter2D`` accepts an inverter object or a checkpoint path, but
   the current ``EMInverter2D`` class does not expose the explicit public
   ``save()/load()`` pair shown by some older examples/docstrings. Prefer a
   fitted in-memory object or a project-tested persistence mechanism; do not
   assume arbitrary ``.npz`` restoration works without verification.

18. Assess uncertainty
----------------------

``PINNInverter2D`` returns one optimized solution and a total loss history. It
does not itself return a posterior or calibrated interval.

Assess uncertainty with scenario ensembles:

* TE versus TM;
* layer count and ``depth_max``;
* vertical and lateral smoothness weights;
* learning rate and epochs;
* data masks and station exclusions;
* alternative initialization, including hybrid Stage 1;
* acceptable classical inversion configurations;
* data perturbations consistent with measurement errors.

Summarize model spread on a common physical depth grid and distinguish
optimizer variability from structural and data uncertainty. See
:doc:`uncertainty` for the wider framework.

For scenario set :math:`\mathcal{C}`, report spread in the same quantity used
for interpretation. For example, if :math:`u_c(z,x)` is the log-resistivity
section for scenario :math:`c`, a simple scenario standard deviation is

.. math::

   \sigma_u(z,x)
   =
   \sqrt{
      \frac{1}{|\mathcal{C}|-1}
      \sum_{c\in\mathcal{C}}
      \left(u_c(z,x)-\bar{u}(z,x)\right)^2
   }.

This is not a posterior distribution. It is a sensitivity diagnostic showing
where the interpretation depends on mode, regularization, layer count, data
masks, or optimizer settings.

19. Preserve a run record
-------------------------

.. code-block:: text

   pinn2d/L18_te_r001/
   ├── manifest.yml
   ├── input/
   │   ├── station_order.csv
   │   ├── qc_reference.yml
   │   └── frequency_coverage.csv
   ├── configuration/
   │   └── pinn2d.yml
   ├── outputs/
   │   ├── log10_resistivity.npy
   │   ├── thickness_m.npy
   │   ├── interface_depth_m.npy
   │   ├── convergence.csv
   │   └── residuals.csv
   ├── figures/
   └── review/
       ├── scenario_comparison.csv
       └── scientific_review.md

Record backend, device, software version, station order, modes/components,
frequency grid, layer/depth parameterization, weights, optimizer, epochs,
residual definitions, scenario set, classical comparison, reviewer, and
status.

Complete example
----------------

.. code-block:: pycon

   >>> from pathlib import Path  # doctest: +SKIP
   >>> import json  # doctest: +SKIP
   >>> import numpy as np  # doctest: +SKIP
   >>> from pycsamt.ai.inversion import PINNInverter2D  # doctest: +SKIP
   >>> from pycsamt.emtools import ensure_sites  # doctest: +SKIP
   >>>
   >>> output = Path("pinn2d/L18_te_r001")  # doctest: +SKIP
   >>> (output / "outputs").mkdir(parents=True, exist_ok=True)  # doctest: +SKIP
   >>> sites = ensure_sites(  # doctest: +SKIP
   ...     "data/AMT/WILLY_DATA/L18",
   ...     recursive=True,
   ...     verbose=0,
   ... )
   >>> inversion = PINNInverter2D(  # doctest: +SKIP
   ...     sites,
   ...     n_layers=10,
   ...     depth_max=2000.0,
   ...     n_freqs=32,
   ...     mode="te",
   ...     smoothness_weight=0.01,
   ...     lateral_weight=0.005,
   ...     epochs=300,
   ...     lr=1e-2,
   ...     comp_te="xy",
   ...     comp_tm="yx",
   ...     verbose=0,
   ... ).fit(verbose=True, log_every=25)
   >>> log10_rho = inversion.resistivity_section()  # doctest: +SKIP
   >>> thickness_m = inversion.thickness_section()  # doctest: +SKIP
   >>> interface_depth_m = np.cumsum(thickness_m, axis=0)  # doctest: +SKIP
   >>> convergence = inversion.convergence_curve()  # doctest: +SKIP
   >>> residuals = inversion.residuals()  # doctest: +SKIP
   >>> np.save(output / "outputs" / "log10_resistivity.npy", log10_rho)  # doctest: +SKIP
   >>> np.save(output / "outputs" / "thickness_m.npy", thickness_m)  # doctest: +SKIP
   >>> np.save(output / "outputs" / "interface_depth_m.npy", interface_depth_m)  # doctest: +SKIP
   >>> convergence.to_csv(output / "outputs" / "convergence.csv", index=False)  # doctest: +SKIP
   >>> residuals.to_csv(output / "outputs" / "residuals.csv", index=False)  # doctest: +SKIP
   >>> manifest = {  # doctest: +SKIP
   ...     "workflow": "PINNInverter2D",
   ...     "physics": "per_station_mt1d_with_lateral_smoothing",
   ...     "mode": "te",
   ...     "component": "xy",
   ...     "stations": inversion.stations,
   ...     "n_layers": 10,
   ...     "depth_max_m": 2000.0,
   ...     "n_freqs": 32,
   ...     "smoothness_weight": 0.01,
   ...     "lateral_weight": 0.005,
   ...     "epochs": 300,
   ...     "learning_rate": 1e-2,
   ... }
   >>> (output / "manifest.json").write_text(  # doctest: +SKIP
   ...     json.dumps(manifest, indent=2),
   ...     encoding="utf-8",
   ... )

Review checklist
----------------

.. list-table::
   :header-rows: 1
   :widths: 31 69

   * - Check
     - Required evidence
   * - Pseudo-2-D scope is explicit
     - Per-station MT1D physics, lateral coupling, and missing full 2-D effects.
   * - Profile order is verified
     - Station names, chainage, coordinates, direction, gaps, and duplicates.
   * - Modes are defensible
     - Strike convention, TE/TM components, separate results, and any averaging.
   * - Frequency coverage is adequate
     - Valid cells by station, common grid, missingness, and supported depth.
   * - Parameterization is tested
     - Layers, depth, thickness stability, and interface-depth scenarios.
   * - Regularization is appraised
     - Vertical/lateral weight grid, roughness, response fit, and stable targets.
   * - Optimization converges
     - Loss curve, learning-rate/epoch sensitivity, repeats, and no divergence.
   * - Residuals are reviewed
     - Apparent resistivity and phase by station/frequency, not only one RMS.
   * - Baselines are included
     - Independent 1-D, classical 2-D/3-D where appropriate, and geology.
   * - Uncertainty is conditional
     - Data perturbation, configuration scenarios, initialization, and omitted
       physics.
   * - Run is reproducible
     - Input/QC IDs, order, configuration, arrays, backend, residuals, reviewer,
       and status.

Common mistakes
---------------

Avoid these errors:

* describing the workflow as a full numerical 2-D EM inversion;
* allowing arbitrary input order to define lateral neighbors;
* treating layer index as physical depth;
* calling ``depth_max`` the depth of investigation;
* selecting ``both`` mode without recognizing TE/TM averaging;
* assuming lateral smoothing accounts for irregular station spacing;
* selecting weights for visual smoothness alone;
* interpreting total loss as pure data misfit;
* reporting the ``both`` residual table as residuals against averaged inputs;
* comparing PINN loss directly with a differently weighted classical RMS;
* assuming a single optimized solution provides uncertainty;
* discarding thickness outputs when plotting resistivity;
* relying on outdated Hybrid/EMInverter2D checkpoint examples without testing
  persistence.

Next steps
----------

Continue with:

* :doc:`validation` for acceptance and classical comparison;
* :doc:`uncertainty` for scenario and data-perturbation analysis;
* :doc:`inference` for distinctions between optimized and surrogate workflows;
* :doc:`agents` for PINN orchestration;
* :doc:`reporting` for pseudo-2-D provenance and deliverables;
* :doc:`../models/occam2d` for a classical smooth 2-D integration.
