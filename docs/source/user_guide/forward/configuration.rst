.. _forward_configuration:

Forward Configuration
=====================

Forward modelling should be reproducible. A synthetic response, training
dataset, or grid experiment is only useful when another user can reconstruct
the solver, model prior, sampled axis, station layout, noise model, and
random seed.

The forward package provides three configuration dataclasses:

* :class:`pycsamt.forward.ForwardConfig` for 1-D MT, CSAMT, and TEM dataset
  generation;
* :class:`pycsamt.forward.ForwardConfig2D` for 2-D MT finite-difference runs;
* :class:`pycsamt.forward.ForwardConfig3D` for quasi-3-D MT forward runs.

All three classes follow the same workflow:

#. write an annotated template;
#. edit the template as the project source of truth;
#. load the template back into Python;
#. validate before computation;
#. build arrays, grids, solvers, or dataset keyword arguments from the config;
#. archive the config next to datasets, figures, and interpretation notes.

Configuration Classes
---------------------

.. list-table::
   :header-rows: 1
   :widths: 24 34 42

   * - Class
     - Main purpose
     - Main helpers
   * - ``ForwardConfig``
     - 1-D MT/CSAMT/TEM synthetic dataset generation.
     - ``freq_grid()``, ``time_grid()``, ``to_dataset_kwargs()``,
       ``write_template()``, ``from_file()``, ``summary()``.
   * - ``ForwardConfig2D``
     - 2-D MT finite-difference grid and solver setup.
     - ``freq_grid()``, ``to_grid()``, ``to_solver_kwargs()``,
       ``write_template()``, ``from_file()``, ``summary()``.
   * - ``ForwardConfig3D``
     - Quasi-3-D MT grid and solver setup.
     - ``freq_grid()``, ``to_grid()``, ``to_solver_kwargs()``,
       ``write_template()``, ``from_file()``, ``summary()``.

Template Files
--------------

Forward configuration files can be written as Python, JSON, YML, or YAML. The
file extension selects the format automatically. Templates include parameter
comments, which makes them good review artifacts for scientific projects.

.. code-block:: python
   :linenos:

   from pycsamt.forward import ForwardConfig, ForwardConfig2D, ForwardConfig3D

   ForwardConfig.write_template("configs/forward_1d.yml")
   ForwardConfig2D.write_template("configs/forward_2d.yml")
   ForwardConfig3D.write_template("configs/forward_3d.yml")

Load the edited files with strict validation of keys:

.. code-block:: python
   :linenos:

   from pycsamt.forward import ForwardConfig

   cfg = ForwardConfig.from_file("configs/forward_1d.yml", strict=True)
   cfg.validate()
   print(cfg.summary())

``strict=True`` is recommended for production runs because misspelled keys
raise an error instead of being ignored. Use ``strict=False`` only when
loading a file that intentionally contains extra metadata.

1-D Configuration
-----------------

``ForwardConfig`` controls the 1-D solver, frequency or time grid, random
layered model prior, noise model, output path, and parallel worker count.

.. list-table::
   :header-rows: 1
   :widths: 28 72

   * - Setting group
     - Fields
   * - Solver
     - ``solver``, ``freq_min``, ``freq_max``, ``n_freqs``, ``time_min``,
       ``time_max``, ``n_times``, ``loop_radius``.
   * - Earth model prior
     - ``n_layers_min``, ``n_layers_max``, ``rho_min``, ``rho_max``,
       ``depth_max``, ``geology``.
   * - Dataset
     - ``n_samples``, ``noise_level``, ``noise_type``, ``include_phase``,
       ``seed``, ``n_jobs``.
   * - Output
     - ``output_dir``, ``output_name``, ``verbose``.

The solver controls which sampled axis is active:

.. list-table::
   :header-rows: 1
   :widths: 24 38 38

   * - ``solver``
     - Active axis
     - Main settings
   * - ``"mt1d"``
     - Frequency
     - ``freq_min``, ``freq_max``, ``n_freqs``, ``include_phase``.
   * - ``"csamt1d"``
     - Frequency
     - ``freq_min``, ``freq_max``, ``n_freqs``, ``include_phase``.
   * - ``"tem1d"``
     - Time gates
     - ``time_min``, ``time_max``, ``n_times``, ``loop_radius``.

Generate a 1-D MT training dataset:

.. code-block:: python
   :linenos:

   from pycsamt.forward import ForwardConfig, generate_dataset

   cfg = ForwardConfig(
       solver="mt1d",
       freq_min=1e-3,
       freq_max=1e4,
       n_freqs=40,
       n_layers_min=3,
       n_layers_max=7,
       rho_min=1.0,
       rho_max=10_000.0,
       depth_max=3000.0,
       n_samples=5000,
       noise_level=0.05,
       noise_type="field",
       include_phase=True,
       seed=42,
       n_jobs=1,
       output_dir="runs/forward",
       output_name="mt1d_training",
   )

   cfg.validate()
   dataset = generate_dataset(**cfg.to_dataset_kwargs())

``to_dataset_kwargs()`` builds the correct frequency or time array and
constructs the output path. This reduces the risk of the config and the
actual dataset generation call drifting apart.

Geological Priors
-----------------

``ForwardConfig.geology`` can replace broad ``rho_min``/``rho_max`` and
``depth_max`` sampling with a named geological prior. This is useful when an
AI dataset should represent a known target class instead of arbitrary
layered-earth variation.

Common prior names include:

.. code-block:: text
   :linenos:

   sedimentary
   crystalline
   geothermal
   marine
   permafrost

.. code-block:: python
   :linenos:

   from pycsamt.forward import ForwardConfig, generate_dataset

   cfg = ForwardConfig(
       solver="mt1d",
       geology="geothermal",
       n_layers_min=4,
       n_layers_max=7,
       n_samples=20_000,
       noise_type="field",
       seed=12,
       output_dir="runs/forward",
       output_name="geothermal_mt1d",
   )

   cfg.validate()
   dataset = generate_dataset(**cfg.to_dataset_kwargs())

When ``geology`` is set, the geological prior should be documented in the
project notes. A neural network trained on a narrow prior may not generalize
outside that geological setting.

TEM Configuration
-----------------

TEM uses time gates rather than frequency samples.

.. code-block:: python
   :linenos:

   from pycsamt.forward import ForwardConfig, generate_dataset

   cfg = ForwardConfig(
       solver="tem1d",
       time_min=1e-6,
       time_max=1e-2,
       n_times=25,
       loop_radius=50.0,
       n_layers_min=3,
       n_layers_max=6,
       n_samples=1000,
       noise_type="gaussian",
       noise_level=0.03,
       seed=5,
       output_dir="runs/forward",
       output_name="tem1d_training",
   )

   cfg.validate()
   dataset = generate_dataset(**cfg.to_dataset_kwargs())

TEM generation can be slower than MT1D because the current implementation
uses numerical integration for the step-off response. Start with a small
``n_samples`` value, inspect the response, then scale up.

2-D Configuration
-----------------

``ForwardConfig2D`` creates a :class:`pycsamt.forward.Grid2D` and solver
keyword arguments for :class:`pycsamt.forward.MT2DForward`.

.. list-table::
   :header-rows: 1
   :widths: 28 72

   * - Setting group
     - Fields
   * - Solver
     - ``freq_min``, ``freq_max``, ``n_freqs``.
   * - Grid
     - ``nx``, ``nz``, ``x_max``, ``z_max``, ``n_pad``, ``pad_factor``.
   * - Earth model
     - ``bg_rho``, ``model_type``, ``anomaly_rho``, ``anomaly_x_lo``,
       ``anomaly_x_hi``, ``anomaly_z_lo``, ``anomaly_z_hi``.
   * - Stations
     - ``n_stations``, ``station_x_max``.
   * - Output
     - ``verbose``.

The supported ``model_type`` values are:

.. list-table::
   :header-rows: 1
   :widths: 28 72

   * - ``model_type``
     - Meaning
   * - ``"halfspace"``
     - Uniform background resistivity.
   * - ``"anomaly"``
     - Rectangular anomaly inside a background.
   * - ``"random"``
     - Random 2-D grid model generated by ``Grid2D.random``.

Example 2-D anomaly run:

.. code-block:: python
   :linenos:

   from pycsamt.forward import ForwardConfig2D, MT2DForward

   cfg = ForwardConfig2D(
       freq_min=1e-2,
       freq_max=1e3,
       n_freqs=25,
       nx=50,
       nz=35,
       x_max=10_000.0,
       z_max=6000.0,
       n_pad=8,
       pad_factor=1.3,
       bg_rho=300.0,
       model_type="anomaly",
       anomaly_rho=10.0,
       anomaly_x_lo=2500.0,
       anomaly_x_hi=6500.0,
       anomaly_z_lo=400.0,
       anomaly_z_hi=1800.0,
       n_stations=16,
       verbose=True,
   )

   cfg.validate()
   print(cfg.summary())

   grid = cfg.to_grid()
   solver = MT2DForward(grid=grid, **cfg.to_solver_kwargs())
   response = solver.run()

``to_grid()`` accepts a ``seed`` argument for random models:

.. code-block:: python
   :linenos:

   cfg = ForwardConfig2D(model_type="random")
   cfg.validate()
   grid = cfg.to_grid(seed=42)

2-D Grid Tuning
---------------

2-D forward runs are sensitive to grid design.

.. list-table::
   :header-rows: 1
   :widths: 28 72

   * - Parameter
     - Practical guidance
   * - ``x_max``
     - Should include all stations and enough side room that boundaries do
       not dominate the response.
   * - ``z_max``
     - Should extend below the expected maximum investigation depth.
   * - ``nx`` and ``nz``
     - Increase when anomalies are small or response gradients are sharp.
   * - ``n_pad``
     - Add side and bottom padding to reduce boundary artifacts.
   * - ``pad_factor``
     - Values around ``1.2`` to ``1.5`` are typical; ``1.3`` is a stable
       starting point.
   * - ``n_stations``
     - Should resolve the lateral scale of the target.

3-D Configuration
-----------------

``ForwardConfig3D`` builds a :class:`pycsamt.forward.Grid3D` and keyword
arguments for :class:`pycsamt.forward.MT3DForward`.

.. list-table::
   :header-rows: 1
   :widths: 28 72

   * - Setting group
     - Fields
   * - Solver
     - ``freq_min``, ``freq_max``, ``n_freqs``, ``method``.
   * - Grid
     - ``nx``, ``ny``, ``nz``, ``x_max``, ``y_max``, ``z_max``,
       ``n_pad``, ``pad_factor``.
   * - Earth model
     - ``bg_rho``, ``model_type``, anomaly bounds, ``n_layers``,
       ``lateral_variation``, ``corr_length``.
   * - Stations
     - ``nx_stations``, ``ny_stations``.
   * - Output
     - ``verbose``.

The supported ``model_type`` values are:

.. list-table::
   :header-rows: 1
   :widths: 28 72

   * - ``model_type``
     - Meaning
   * - ``"halfspace"``
     - Uniform 3-D background.
   * - ``"block_anomaly"``
     - Rectangular 3-D anomaly inside a background.
   * - ``"random_layered"``
     - Random horizontal layers with optional Gaussian-random-field lateral
       variation.

Example quasi-3-D block anomaly:

.. code-block:: python
   :linenos:

   from pycsamt.forward import ForwardConfig3D, MT3DForward

   cfg = ForwardConfig3D(
       freq_min=1e-2,
       freq_max=1e3,
       n_freqs=15,
       method="quasi3d",
       nx=24,
       ny=24,
       nz=18,
       x_max=9000.0,
       y_max=9000.0,
       z_max=5000.0,
       n_pad=8,
       bg_rho=500.0,
       model_type="block_anomaly",
       anomaly_rho=20.0,
       anomaly_x_lo=2500.0,
       anomaly_x_hi=6500.0,
       anomaly_y_lo=2500.0,
       anomaly_y_hi=6500.0,
       anomaly_z_lo=500.0,
       anomaly_z_hi=2000.0,
       nx_stations=6,
       ny_stations=6,
       verbose=True,
   )

   cfg.validate()
   print(cfg.summary())

   grid = cfg.to_grid()
   response = MT3DForward(grid=grid, **cfg.to_solver_kwargs()).run()

Random layered 3-D model:

.. code-block:: python
   :linenos:

   cfg = ForwardConfig3D(
       model_type="random_layered",
       n_layers=5,
       lateral_variation=True,
       corr_length=2500.0,
       nx_stations=7,
       ny_stations=7,
   )

   cfg.validate()
   grid = cfg.to_grid(seed=42)

3-D Configuration Notes
----------------------

``ForwardConfig3D.method`` currently validates to ``"quasi3d"``. The
underlying solver code contains experimental hooks for fuller 3-D methods,
but the documented configuration path is quasi-3-D. Treat quasi-3-D outputs
as survey-scale synthetic responses, not as final production 3-D inversion
results.

For production 3-D modelling or inversion, continue to:

* :doc:`../models/modem`
* :doc:`../models/mare2dem`

Validation
----------

Call ``validate()`` before any expensive run. Validation catches basic range
errors such as negative frequencies, invalid model types, impossible anomaly
bounds, and invalid station counts.

.. code-block:: python
   :linenos:

   from pycsamt.forward import ForwardConfig2D

   cfg = ForwardConfig2D(freq_min=100.0, freq_max=10.0)

   try:
       cfg.validate()
   except ValueError as exc:
       print(f"Configuration problem: {exc}")

Validation does not prove that the model is geologically meaningful. It only
checks that the parameters are internally acceptable for the builder and
solver.

Summaries And Provenance
------------------------

Each config object provides ``summary()`` and ``repr`` output that can be
saved in logs or reports.

.. code-block:: python
   :linenos:

   from pathlib import Path

   cfg = ForwardConfig2D(model_type="anomaly")
   cfg.validate()

   run_dir = Path("runs/forward/2d_anomaly")
   run_dir.mkdir(parents=True, exist_ok=True)
   (run_dir / "summary.txt").write_text(cfg.summary())
   cfg.to_template(run_dir / "forward_config_2d.yml")

What To Record
--------------

For reproducibility, record:

* config file path and format;
* solver type and dimensionality;
* frequency grid or time-gate grid;
* model prior or grid constructor;
* resistivity bounds, anomaly bounds, or geology prior;
* station layout;
* noise type and noise level;
* random seed;
* output dataset or figure paths;
* pyCSAMT version or commit;
* any manual edits made after grid creation.

Recommended Run Layout
----------------------

.. code-block:: text
   :linenos:

   runs/
     forward/
       mt1d_training/
         forward_config.yml
         summary.txt
         mt1d_training.npz
         sample_responses.png
       profile_2d_anomaly/
         forward_config_2d.yml
         summary.txt
         model.png
         pseudosection_te.png
       survey_quasi3d/
         forward_config_3d.yml
         summary.txt
         response_maps/

Common Mistakes
---------------

``The generated dataset cannot be reproduced.``
   Record ``seed``, configuration file, pyCSAMT version, and output path.
   Avoid manually changing keyword arguments outside the config.

``The run validates but the response is not useful.``
   Validation only checks parameter ranges. Review frequency range, grid
   extents, anomaly size, station spacing, and plots.

``TEM generation is unexpectedly slow.``
   Start with fewer samples or fewer time gates. TEM1D currently uses
   numerical integration for the step-off response.

``The 2-D model has boundary artifacts.``
   Increase ``x_max``, ``z_max``, ``n_pad``, or adjust ``pad_factor``.

``The 3-D response is treated as a production inversion result.``
   ``ForwardConfig3D`` documents quasi-3-D synthetic modelling. Use ModEM or
   MARE2DEM for production external-engine workflows.

Next Steps
----------

* :doc:`solvers_and_grids` explains how configs create models and solvers.
* :doc:`synthetic_datasets` explains dataset generation and train/validation
  splits.
* :doc:`plotting` shows how to inspect configured runs.
* :doc:`forward_to_inversion` explains synthetic recovery workflows.
