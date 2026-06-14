.. _forward_synthetic_datasets:

Synthetic Datasets And Noise
============================

Synthetic datasets are central to AI inversion, regression tests, solver
benchmarks, and method demonstrations. pyCSAMT forward datasets store the
computed response as ``X`` and the model parameters as ``y``.

1-D Dataset Generation
----------------------

Use :func:`pycsamt.forward.generate_dataset` directly or through
:class:`pycsamt.forward.ForwardConfig`.

.. code-block:: python
   :linenos:

   import numpy as np

   from pycsamt.forward import generate_dataset

   dataset = generate_dataset(
       solver="mt1d",
       n_samples=1000,
       freqs=np.logspace(-3, 4, 40),
       n_layers=(3, 7),
       rho_range=(1.0, 10_000.0),
       depth_max=3000.0,
       noise_level=0.05,
       noise_type="field",
       include_phase=True,
       seed=42,
       n_jobs=1,
       output="runs/forward/mt1d_dataset.npz",
   )

   print(dataset.X.shape)
   print(dataset.y.shape)

``ForwardDataset`` can be saved, loaded, and split:

.. code-block:: python
   :linenos:

   from pycsamt.forward import ForwardDataset

   dataset = ForwardDataset.load("runs/forward/mt1d_dataset.npz")
   train, val, test = dataset.split(
       val_frac=0.1,
       test_frac=0.1,
       seed=0,
   )

Feature And Target Layout
-------------------------

For MT and CSAMT, a 1-D response feature vector usually contains
``log10(rho_a)`` and optionally phase. For TEM, the feature vector contains a
log-scaled transient response.

Targets are flattened model parameters. For a 1-D model, the layout is:

.. code-block:: text
   :linenos:

   [log10(rho_0), log10(rho_1), ..., log10(rho_n),
    thickness_0, thickness_1, ..., thickness_n_minus_1]

When layer counts vary, samples with fewer layers are padded by the dataset
generator. Keep this in mind when training neural networks or computing loss
functions.

Noise Models
------------

Noise makes synthetic data more realistic and prevents AI models from
learning only idealized solver curves.

.. list-table::
   :header-rows: 1
   :widths: 28 72

   * - Noise model
     - Use
   * - ``GaussianNoise``
     - Independent additive perturbations scaled by ``noise_level``.
   * - ``MultiplicativeNoise``
     - Pointwise scaling by a random factor, useful for relative response
       uncertainty.
   * - ``FieldRealisticNoise``
     - Frequency-dependent noise and occasional spikes intended to resemble
       field behaviour more closely.

.. code-block:: python
   :linenos:

   from pycsamt.forward import FieldRealisticNoise

   noisy_response = FieldRealisticNoise(
       base_level=0.05,
   ).apply(response, seed=42)

Pseudo-3-D Survey Datasets
--------------------------

:func:`pycsamt.forward.generate_dataset_3d` creates multi-station survey
datasets for graph-style or spatial AI models. Each survey has a fixed
station layout and spatially correlated layer resistivities.

.. code-block:: python
   :linenos:

   from pycsamt.forward import generate_dataset_3d

   surveys = generate_dataset_3d(
       n_surveys=500,
       n_stations=25,
       n_layers=4,
       extent=10_000.0,
       corr_length=2000.0,
       noise_level=0.03,
       include_phase=True,
       seed=7,
       output="runs/forward/survey3d_dataset.npz",
   )

   print(surveys.X.shape)
   print(surveys.y.shape)
   print(surveys.coords.shape)

Dataset QA
----------

Before using a synthetic dataset for training or benchmarking:

* plot random samples;
* check for NaN or infinite values;
* confirm units and log transforms;
* verify feature length matches the model architecture;
* compare noisy and noise-free examples;
* keep train, validation, and test splits deterministic;
* save the configuration used to create the dataset.
