Forward Commands
================

``pycsamt forward`` runs 1-D electromagnetic forward models and builds
synthetic datasets for examples, QA, and AI inversion experiments. The
commands are thin CLI wrappers around the real ``pycsamt.forward`` API:
``LayeredModel``, ``MT1DForward``, ``TEM1DForward``, ``CSAMT1DForward``,
``GEOLOGY_PRIORS``, and ``generate_dataset``.

Use this group when you need a quick response table, a reproducible
geology-sampled model, or a batch of model/response pairs without writing
Python.

Command Map
-----------

.. list-table::
   :header-rows: 1
   :widths: 34 66

   * - Command
     - Purpose
   * - ``pycsamt forward model geology``
     - List all built-in geological priors or inspect one prior in detail.
   * - ``pycsamt forward run``
     - Run one explicit or geology-sampled 1-D model and print the response.
   * - ``pycsamt forward generate``
     - Batch-generate synthetic training data as ``.npz`` or feature-only CSV.

Global Conventions
------------------

All forward subcommands support the common CLI formatting/logging options
where they make sense:

``--format {text,json,csv}``
    Output format for commands that print tables or records.

``-v`` / ``--verbose``
    Increase progress logging. For ``run`` and ``generate``, verbose mode
    writes model/run metadata to stderr.

``--no-color``
    Disable colored log output.

Frequencies are specified with ``FMIN:FMAX`` in Hz. The command creates
``--n-freqs`` log-spaced samples between the two endpoints.

.. code-block:: console

   pycsamt forward run --geology sedimentary --freqs 0.001:10000 --n-freqs 30

``forward model geology``
-------------------------

List or inspect the built-in geological priors used to draw random
``LayeredModel`` objects. These priors live in
``pycsamt.forward.synthetic.GEOLOGY_PRIORS``.

.. code-block:: console

   pycsamt forward model geology
   pycsamt forward model geology --format json
   pycsamt forward model geology --format csv > geology_priors.csv

Show one scenario in detail:

.. code-block:: console

   pycsamt forward model geology --name sedimentary
   pycsamt forward model geology --name geothermal --format json

Built-in scenario names are:

``sedimentary``
    Layered basin-style resistivity ranges.

``crystalline``
    Resistive basement and fractured-rock style ranges.

``geothermal``
    Hot-fluid and alteration-zone style ranges.

``marine``
    Conductive seawater/sediment style ranges.

``permafrost``
    Frozen-ground style ranges.

The list view includes ``scenario``, ``n_layers``, ``log10_rho_range``,
``depth_max_m``, and a short ``description``. The detail view prints the
full prior fields for the selected scenario.

Options
~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 26 22 52

   * - Option
     - Default
     - Meaning
   * - ``--name``, ``-n``
     - none
     - Scenario to inspect. When omitted, all scenarios are listed.
   * - ``--format``
     - ``text``
     - Print a text table, JSON records, or CSV.
   * - ``-v``, ``--verbose``
     - ``0``
     - Increase logging verbosity.
   * - ``--no-color``
     - false
     - Disable colored log output.

``forward run``
---------------

Run one 1-D model and print the forward response. You must define the
earth model in exactly one of two ways:

1. explicit ``--resistivities`` plus ``--thicknesses``;
2. ``--geology`` to draw a random model from a built-in prior.

Explicit Models
~~~~~~~~~~~~~~~

Layer resistivities are in Ohm.m, ordered from top to bottom.
Thicknesses are in metres. For ``N`` resistivity values, pass exactly
``N - 1`` thicknesses because the bottom half-space has no finite
thickness.

.. code-block:: console

   pycsamt forward run \
       --resistivities 100,10,500 \
       --thicknesses 300,800

CSV output is useful for spreadsheets and plotting scripts:

.. code-block:: console

   pycsamt forward run \
       --resistivities 50,500,20,800 \
       --thicknesses 200,600,1500 \
       --freqs 0.01:1000 \
       --n-freqs 40 \
       --format csv > response.csv

Geology-Sampled Models
~~~~~~~~~~~~~~~~~~~~~~

``--geology`` draws a model from one of the built-in priors. Use
``--seed`` to make the draw reproducible.

.. code-block:: console

   pycsamt forward run --geology sedimentary --seed 42
   pycsamt forward run --geology geothermal --solver mt --freqs 0.01:1000

Do not combine ``--geology`` with ``--resistivities`` or
``--thicknesses``.

Solvers
~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 20 40 40

   * - ``--solver``
     - Backend class
     - Output columns
   * - ``mt``
     - ``MT1DForward``
     - ``freq_Hz``, ``period_s``, ``rho_a_Ohm_m``, ``phase_deg``
   * - ``csamt``
     - ``CSAMT1DForward``
     - ``freq_Hz``, ``period_s``, ``rho_a_Ohm_m``, ``phase_deg``
   * - ``tem``
     - ``TEM1DForward``
     - ``time_s``, ``dBz_dt_T_per_s``

Example TEM run:

.. code-block:: console

   pycsamt forward run \
       --resistivities 200,5,1000 \
       --thicknesses 100,500 \
       --solver tem \
       --format json

Options
~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 28 20 52

   * - Option
     - Default
     - Meaning
   * - ``--resistivities``, ``-r``
     - none
     - Comma-separated layer resistivities in Ohm.m.
   * - ``--thicknesses``, ``-t``
     - none
     - Comma-separated layer thicknesses in metres.
   * - ``--geology``, ``-g``
     - none
     - Draw a random model from a built-in geological scenario.
   * - ``--seed``
     - none
     - Random seed for reproducible geology-sampled models.
   * - ``--solver``
     - ``mt``
     - One of ``mt``, ``tem``, or ``csamt``.
   * - ``--freqs``
     - ``0.001:10000``
     - Frequency range in Hz.
   * - ``--n-freqs``
     - ``30``
     - Number of log-spaced frequencies. Minimum is ``2``.
   * - ``--format``
     - ``text``
     - Print text, JSON, or CSV.
   * - ``-v``, ``--verbose``
     - ``0``
     - Print run metadata to stderr.
   * - ``--no-color``
     - false
     - Disable colored logs.

Validation Rules
~~~~~~~~~~~~~~~~

``forward run`` fails fast when:

* neither explicit layers nor ``--geology`` is provided;
* ``--resistivities`` is provided without ``--thicknesses``;
* the number of thicknesses is not one fewer than the number of
  resistivities;
* ``--geology`` is combined with explicit layer values;
* the geology scenario name is unknown;
* ``--n-freqs`` is less than ``2``.

``forward generate``
--------------------

Generate many random 1-D models and their synthetic responses. This is
the CLI entry point for training-data creation and wraps
``pycsamt.forward.generate_dataset``.

Quick MT dataset:

.. code-block:: console

   pycsamt forward generate --n-samples 5000 --output mt1d_train.npz

Large reproducible field-noise dataset:

.. code-block:: console

   pycsamt forward generate \
       --solver mt1d \
       --n-samples 50000 \
       --n-layers 3:8 \
       --freqs 0.001:10000 \
       --n-freqs 30 \
       --noise-type field \
       --noise-level 0.03 \
       --seed 42 \
       --jobs 4 \
       --output mt1d_field.npz

TEM dataset:

.. code-block:: console

   pycsamt forward generate \
       --solver tem1d \
       --n-samples 2000 \
       --noise-type multiplicative \
       --output tem1d_train.npz

Solvers and Noise
~~~~~~~~~~~~~~~~~

``--solver`` accepts ``mt1d``, ``tem1d``, and ``csamt1d``.

``--noise-type`` accepts:

``none``
    Disable noise. Internally the command passes a zero noise level.

``gaussian``
    Add relative Gaussian noise.

``field``
    Add field-style response noise.

``multiplicative``
    Add multiplicative relative noise.

``--noise-level`` is a relative amplitude, interpreted as a fraction.
For example, ``0.03`` means roughly three percent noise for supported
noise models.

Layer Counts
~~~~~~~~~~~~

``--n-layers`` accepts either a fixed integer or a uniform integer range:

.. code-block:: console

   pycsamt forward generate --n-layers 5
   pycsamt forward generate --n-layers 3:7

Malformed values, such as ``--n-layers bad``, are rejected before
dataset generation starts.

Options
~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 28 20 52

   * - Option
     - Default
     - Meaning
   * - ``--solver``
     - ``mt1d``
     - One of ``mt1d``, ``tem1d``, or ``csamt1d``.
   * - ``--n-samples``, ``-n``
     - ``1000``
     - Number of model/response pairs to generate. Minimum is ``1``.
   * - ``--n-layers``
     - ``3:7``
     - Fixed layer count or ``MIN:MAX`` range.
   * - ``--freqs``
     - ``0.001:10000``
     - Frequency range in Hz.
   * - ``--n-freqs``
     - ``30``
     - Number of log-spaced frequencies. Minimum is ``2``.
   * - ``--noise-type``
     - ``gaussian``
     - One of ``none``, ``gaussian``, ``field``, ``multiplicative``.
   * - ``--noise-level``
     - ``0.05``
     - Relative noise amplitude.
   * - ``--seed``
     - none
     - Global random seed for reproducible datasets.
   * - ``--output``, ``-o``
     - ``forward_dataset.npz``
     - Output file. ``.npz`` writes the full dataset; other suffixes write CSV.
   * - ``--jobs``
     - project default
     - Number of worker jobs passed to ``generate_dataset``.
   * - ``-v``, ``--verbose``
     - ``0``
     - Print generation metadata to stderr.
   * - ``--no-color``
     - false
     - Disable colored logs.

Dataset Files
-------------

When ``--output`` ends in ``.npz``, the command writes a compressed NumPy
archive from ``generate_dataset``. The expected arrays are:

``X``
    Feature matrix. For MT/CSAMT, features are based on
    ``log10(apparent resistivity)`` and phase. For TEM, features are
    based on ``log10(|dBz/dt|)``.

``y``
    Target matrix containing model parameters. The CLI help describes it
    as ``log10(rho)`` followed by thickness values, padded to the maximum
    sampled layer count.

``meta``
    Per-sample metadata, including layer counts and noise information
    when supplied by the generator.

When ``--output`` does not end in ``.npz``, the command writes only
``ds.X`` as CSV with columns named ``feature_0``, ``feature_1``, and so
on. Use ``.npz`` when you need the target matrix ``y`` for training.

Typical Workflows
-----------------

Inspect a prior, sample one response, then generate training data from
the same conceptual setting:

.. code-block:: console

   pycsamt forward model geology
   pycsamt forward model geology --name sedimentary
   pycsamt forward run --geology sedimentary --seed 7 --format csv > sedimentary_response.csv
   pycsamt forward generate --solver mt1d --n-samples 10000 --seed 7 --output train.npz

Create a small smoke-test dataset before launching a larger run:

.. code-block:: console

   pycsamt forward generate --n-samples 25 --n-freqs 8 --output smoke.npz
   pycsamt forward generate --n-samples 50000 --jobs 4 --output mt1d_train.npz

Troubleshooting
---------------

``Error: Provide --resistivities + --thicknesses, or --geology.``
    Add either explicit layer values or ``--geology SCENARIO``.

``--geology cannot be combined with --resistivities / --thicknesses.``
    Choose one model source. Use explicit layers for deterministic
    manual models, or ``--geology`` plus ``--seed`` for reproducible
    random models.

``Expected N thickness value(s) for M layers``
    The bottom layer is a half-space. Pass one fewer thickness than
    resistivity.

``Expected an integer (e.g. 5) or a range (e.g. 3:7).``
    Fix ``--n-layers`` for ``forward generate``.

API Equivalents
---------------

The CLI examples above correspond to direct Python calls such as:

.. code-block:: python

   import numpy as np
   from pycsamt.forward import LayeredModel, MT1DForward

   model = LayeredModel(
       resistivity=np.array([100.0, 10.0, 500.0]),
       thickness=np.array([300.0, 800.0]),
   )
   freqs = np.logspace(-3, 4, 30)
   response = MT1DForward(freqs).run(model)

For dataset generation:

.. code-block:: python

   import numpy as np
   from pycsamt.forward import generate_dataset

   freqs = np.logspace(-3, 4, 30)
   dataset = generate_dataset(
       solver="mt1d",
       n_samples=1000,
       freqs=freqs,
       n_layers=(3, 7),
       noise_level=0.05,
       noise_type="gaussian",
       seed=42,
   )

