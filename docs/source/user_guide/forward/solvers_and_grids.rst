.. _forward_solvers_and_grids:

Solvers And Grids
=================

The forward package separates three responsibilities:

``model containers``
   Objects such as :class:`pycsamt.forward.LayeredModel`,
   :class:`pycsamt.forward.Grid2D`, and :class:`pycsamt.forward.Grid3D` store
   resistivity, geometry, station positions, and metadata.

``solvers``
   Objects such as :class:`pycsamt.forward.MT1DForward`,
   :class:`pycsamt.forward.TEM1DForward`, :class:`pycsamt.forward.MT2DForward`,
   and :class:`pycsamt.forward.MT3DForward` compute synthetic responses from
   the model containers.

``response containers``
   Objects such as :class:`pycsamt.forward.ForwardResponse`,
   :class:`pycsamt.forward.ForwardResponse2D`, and
   :class:`pycsamt.forward.ForwardResponse3D` hold the predicted fields,
   apparent resistivities, phases, time-domain decays, station coordinates, and
   feature-array helpers.

This separation matters because survey design, synthetic dataset generation,
machine-learning workflows, and inversion handoff all need the same response
arrays but may use different model constructors and different solver settings.

Solver Map
----------

The currently documented forward solver families are:

.. list-table::
   :header-rows: 1
   :widths: 20 22 28 30

   * - Solver
     - Model container
     - Primary axis
     - Main outputs
   * - ``MT1DForward``
     - ``LayeredModel``
     - Frequencies in Hz
     - Impedance ``z``, apparent resistivity ``rho_a``, phase.
   * - ``CSAMT1DForward``
     - ``LayeredModel``
     - Frequencies in Hz
     - MT-like response with optional controlled-source correction.
   * - ``TEM1DForward``
     - ``LayeredModel``
     - Time gates in seconds
     - Step-off ``dBz_dt`` response and frequency-domain helper values.
   * - ``MT2DForward``
     - ``Grid2D``
     - Frequencies in Hz
     - TE/TM impedances, apparent resistivities, phases at profile stations.
   * - ``MT3DForward``
     - ``Grid3D``
     - Frequencies in Hz
     - Approximate tensor components on a station grid.

The forward solvers are deterministic for a fixed model, axis, and solver
configuration. Noise is applied later by the noise utilities, which keeps the
physical response and the observation model separate.

1-D Layered Models
------------------

:class:`pycsamt.forward.LayeredModel` represents a stack of horizontal layers.
The final resistivity is the halfspace, so the ``thickness`` array has one
fewer entry than ``resistivity``.

.. code-block:: python
   :linenos:

   import numpy as np

   from pycsamt.forward import LayeredModel, MT1DForward

   model = LayeredModel(
       resistivity=[100.0, 10.0, 500.0],
       thickness=[300.0, 800.0],
       name="conductive_middle_layer",
   )

   freqs = np.logspace(-3, 4, 40)
   response = MT1DForward(freqs=freqs).run(model)

   print(response.rho_a.shape)
   print(response.phase.shape)

Layered models are useful for:

* quick sanity checks against textbook MT behaviour;
* training small synthetic catalogues;
* building starting models for inversion;
* extending a 1-D model laterally into a simple 2-D grid;
* testing whether the chosen frequency or time range can sense the expected
  target depth.

The important depth scale for frequency-domain EM is the skin depth. A common
engineering estimate is:

.. math::

   \delta \approx 503 \sqrt{\frac{\rho}{f}},

where :math:`\delta` is in metres, :math:`\rho` is resistivity in
:math:`\Omega\,m`, and :math:`f` is frequency in Hz. This is not a substitute
for a forward model, but it is a useful first check: high frequencies see
shallower structure, while low frequencies sample deeper structure.

1-D MT And CSAMT Solvers
------------------------

``MT1DForward`` computes a plane-wave natural-source response using layered
earth impedance recursion. ``CSAMT1DForward`` uses the same layered-earth base
and can apply a controlled-source correction when source geometry is supplied.

.. code-block:: python
   :linenos:

   import numpy as np

   from pycsamt.forward import CSAMT1DForward, LayeredModel, MT1DForward

   model = LayeredModel(
       resistivity=[80.0, 25.0, 600.0],
       thickness=[250.0, 900.0],
   )

   freqs = np.logspace(-1, 4, 32)

   mt_response = MT1DForward(freqs=freqs).run(model)

   csamt_response = CSAMT1DForward(
       freqs=freqs,
       source_offset=5000.0,
       dipole_length=1000.0,
   ).run(model)

   mt_features = mt_response.to_array(log_rho=True, include_phase=True)
   csamt_features = csamt_response.to_array(log_rho=True, include_phase=True)

For MT and CSAMT 1-D responses:

* ``response.freqs`` has shape ``(n_freqs,)``;
* ``response.z`` has shape ``(n_freqs,)`` and is complex;
* ``response.rho_a`` has shape ``(n_freqs,)``;
* ``response.phase`` has shape ``(n_freqs,)`` in degrees;
* ``response.to_array()`` returns a one-dimensional feature vector.

When the model is one-dimensional, every surface station would see the same
response. Use 2-D or 3-D grids when lateral geometry is part of the experiment.

1-D TDEM Solver
---------------

``TEM1DForward`` computes a central-loop step-off response for a layered earth.
It uses time gates rather than frequencies as the primary output axis.

.. code-block:: python
   :linenos:

   import numpy as np

   from pycsamt.forward import LayeredModel, TEM1DForward

   model = LayeredModel(
       resistivity=[60.0, 250.0, 900.0],
       thickness=[120.0, 700.0],
   )

   times = np.logspace(-6, -3, 30)

   response = TEM1DForward(
       times=times,
       loop_radius=50.0,
       n_freqs=32,
       n_lam=80,
   ).run(model)

   print(response.times.shape)
   print(response.dBz_dt.shape)

For TDEM responses:

* ``response.times`` has shape ``(n_times,)``;
* ``response.dBz_dt`` has shape ``(n_times,)``;
* ``response.to_array()`` returns a log-scaled decay feature vector by default.

Use early gates to test shallow sensitivity and later gates to test deeper
sensitivity. If the time range is too narrow, the inversion may fit the decay
curve but remain insensitive to the target interval.

2-D Profile Grid Concepts
-------------------------

:class:`pycsamt.forward.Grid2D` stores a finite-difference grid for a profile.
It contains:

* horizontal cell widths ``dx``;
* vertical cell heights ``dz``;
* cell resistivity matrix ``resistivity`` with shape ``(nz, nx)``;
* surface station x-positions ``x_stations``;
* padding count ``n_pad``;
* node and cell-centre coordinate helpers.

The grid uses depth-positive ``z`` coordinates. Resistivity is stored
top-to-bottom and left-to-right. Stations must lie inside the grid extent.

The core region is the scientific model of interest. Padding cells are added on
the left, right, and bottom to reduce boundary influence. The ``n_pad`` value
records how many padding cells were added so plotting and interpretation tools
can distinguish core cells from numerical buffer cells.

.. code-block:: python
   :linenos:

   from pycsamt.forward import Grid2D

   grid = Grid2D.halfspace(
       rho=100.0,
       nx=40,
       nz=28,
       x_max=8000.0,
       z_max=5000.0,
       n_pad=8,
       pad_factor=1.3,
       n_stations=12,
   )

   print(grid.resistivity.shape)
   print(grid.x_stations)
   print(grid.x_nodes[0], grid.x_nodes[-1])

2-D Profile Grid Constructors
-----------------------------

Use the constructor that matches the experiment you want to run.

.. list-table::
   :header-rows: 1
   :widths: 25 35 40

   * - Constructor
     - Purpose
     - Typical use
   * - ``Grid2D.halfspace``
     - Uniform background model.
     - Baseline checks, boundary tests, and solver sanity checks.
   * - ``Grid2D.with_anomaly``
     - Background with one rectangular resistive or conductive block.
     - Target detectability, station-spacing tests, and pseudosection
       training examples.
   * - ``Grid2D.from_1d_layers``
     - Laterally extends a layered model across a profile.
     - Compare 1-D and 2-D solvers, or build a simple starting model.
   * - ``Grid2D.random``
     - Randomized 2-D resistivity model.
     - Synthetic catalogues and stress tests.

A compact block-anomaly grid looks like this:

.. code-block:: python
   :linenos:

   from pycsamt.forward import Grid2D, MT2DForward

   grid = Grid2D.with_anomaly(
       bg_rho=500.0,
       anomaly_rho=5.0,
       anomaly_bounds=(2000.0, 6000.0, 300.0, 1500.0),
       nx=50,
       nz=35,
       x_max=10000.0,
       z_max=6000.0,
       n_pad=8,
       n_stations=16,
   )

   response = MT2DForward(
       freqs=[1.0, 10.0, 100.0],
       grid=grid,
       verbose=True,
   ).run()

``anomaly_bounds`` are specified in core model coordinates as
``(x_lo, x_hi, z_lo, z_hi)``. The constructor internally accounts for padding
when placing the anomaly into the full numerical grid.

2-D MT Solver
-------------

``MT2DForward`` solves the 2-D MT finite-difference problem for TE and TM modes.
It returns a :class:`pycsamt.forward.ForwardResponse2D`.

.. code-block:: python
   :linenos:

   station_0 = response.station_response(0)

   print(station_0["rho_a_te"])
   print(station_0["phase_te"])

   features_te = response.to_feature_array(
       mode="te",
       log_rho=True,
       include_phase=True,
   )

   features_both = response.to_feature_array(
       mode="both",
       log_rho=True,
       include_phase=True,
   )

Response array shapes are always frequency first:

.. list-table::
   :header-rows: 1
   :widths: 30 35 35

   * - Attribute
     - Shape
     - Meaning
   * - ``response.freqs``
     - ``(n_freqs,)``
     - Frequency axis in Hz.
   * - ``response.stations_x``
     - ``(n_stations,)``
     - Surface station positions.
   * - ``response.zxy``
     - ``(n_freqs, n_stations)``
     - TE impedance component.
   * - ``response.zyx``
     - ``(n_freqs, n_stations)``
     - TM impedance component.
   * - ``response.rho_a_te``
     - ``(n_freqs, n_stations)``
     - TE apparent resistivity.
   * - ``response.phase_te``
     - ``(n_freqs, n_stations)``
     - TE phase in degrees.
   * - ``response.rho_a_tm``
     - ``(n_freqs, n_stations)``
     - TM apparent resistivity.
   * - ``response.phase_tm``
     - ``(n_freqs, n_stations)``
     - TM phase in degrees.
   * - ``response.to_feature_array()``
     - ``(n_stations, n_features)``
     - Station-first feature matrix for ML or downstream processing.

This difference is deliberate: physical response arrays are frequency-first,
while feature matrices are station-first because each row is one training or
analysis sample.

When passing a 2-D forward response to the inversion profile API, transpose the
selected response arrays:

.. code-block:: python
   :linenos:

   inversion_data = {
       "freqs": response.freqs,
       "rho_a": response.rho_a_te.T,
       "phase": response.phase_te.T,
       "station_x": response.stations_x,
   }

3-D Volume Grid Concepts
------------------------

:class:`pycsamt.forward.Grid3D` stores a 3-D resistivity volume and a 2-D
station layout. It contains:

* cell widths ``dx``, ``dy``, and ``dz``;
* resistivity with shape ``(nz, ny, nx)``;
* station coordinates ``stations_xy`` with shape ``(n_stations, 2)``;
* padding in x, y, and z;
* helpers that extract XZ and YZ slices for the quasi-3-D solver.

.. code-block:: python
   :linenos:

   from pycsamt.forward import Grid3D

   grid = Grid3D.halfspace(
       rho=100.0,
       nx=20,
       ny=20,
       nz=15,
       x_max=8000.0,
       y_max=8000.0,
       z_max=4000.0,
       n_pad=6,
       nx_stations=5,
       ny_stations=5,
   )

   print(grid.resistivity.shape)
   print(grid.stations_xy.shape)

The regular station grid is created over the core model and shifted into the
padded numerical coordinate system. As in 2-D, stations must lie inside the full
grid extent.

3-D Volume Grid Constructors
----------------------------

The main constructors are:

.. list-table::
   :header-rows: 1
   :widths: 25 35 40

   * - Constructor
     - Purpose
     - Typical use
   * - ``Grid3D.halfspace``
     - Uniform 3-D background.
     - Baseline tensor and station-grid checks.
   * - ``Grid3D.block_anomaly``
     - Background with one 3-D rectangular block.
     - Survey design and quasi-3-D anomaly tests.
   * - ``Grid3D.random_layered``
     - Random laterally variable model.
     - Synthetic training data and robustness tests.

.. code-block:: python
   :linenos:

   from pycsamt.forward import Grid3D, MT3DForward

   grid = Grid3D.block_anomaly(
       bg_rho=500.0,
       anomaly_rho=20.0,
       bounds=(2000.0, 6000.0, 2000.0, 6000.0, 300.0, 1500.0),
       nx=20,
       ny=20,
       nz=15,
       x_max=8000.0,
       y_max=8000.0,
       z_max=4000.0,
       n_pad=6,
       nx_stations=5,
       ny_stations=5,
   )

   response = MT3DForward(
       freqs=[1.0, 10.0, 100.0],
       grid=grid,
   ).run()

Quasi-3-D Solver
----------------

``MT3DForward`` is a quasi-3-D MT solver. It does not claim to replace a full
production 3-D modelling engine. Instead, it approximates tensor responses by
extracting orthogonal 2-D slices from the 3-D grid and running ``MT2DForward``
on those slices.

Conceptually:

* XZ slices are grouped by station y-row;
* YZ slices are grouped by station x-column;
* TE and TM responses from both slice families are combined;
* off-diagonal tensor components ``Z_xy`` and ``Z_yx`` carry the main response;
* diagonal components ``Z_xx`` and ``Z_yy`` are represented as approximate
  tensor outputs and should be interpreted carefully.

The output is :class:`pycsamt.forward.ForwardResponse3D`.

.. code-block:: python
   :linenos:

   print(response.rho_a_xy.shape)
   print(response.phase_xy.shape)
   print(response.stations_xy.shape)

   features = response.to_feature_array(
       components="xy_yx",
       log_rho=True,
       include_phase=True,
   )

   all_components = response.to_feature_array(
       components="all",
       log_rho=True,
       include_phase=False,
   )

As with the 2-D response, physical arrays are frequency-first:
``(n_freqs, n_stations)``. Feature arrays are station-first:
``(n_stations, n_features)``.

Grid Design Rules
-----------------

Forward results are only as useful as the grid that produced them. Before
trusting a response, check these points.

``Positive resistivity``
   Every model cell must have strictly positive resistivity. The grid
   constructors validate this, but it is still worth checking after manual
   edits.

``Station coverage``
   Stations must sit inside the grid. For profile experiments, station spacing
   should be small enough to sample the expected lateral anomaly signature.

``Padding``
   Padding reduces boundary effects. Increase ``n_pad`` or ``pad_factor`` if
   the response changes noticeably when the model extent is increased.

``Near-surface resolution``
   Shallow cells should resolve high-frequency or early-time sensitivity. A
   very coarse top layer can hide near-surface targets or create numerical
   artefacts.

``Depth extent``
   The bottom of the grid should be deeper than the target and deeper than the
   main sensitivity range of the lowest frequencies or latest gates.

``Frequency or time range``
   Do not rely on the solver to compensate for missing physics in the survey
   design. If the chosen axis cannot see the target, the inversion will not
   recover it reliably.

``Dimensionality``
   Use 1-D models for layered checks, 2-D grids for profile structure, and
   quasi-3-D grids for survey design or synthetic AI catalogues. Move to native
   external backends when production 2-D or 3-D inversion files are required.

A compact grid-check helper can be useful in notebooks:

.. code-block:: python
   :linenos:

   import numpy as np

   def describe_grid2d(grid):
       print(f"cells: nz={grid.nz}, nx={grid.nx}")
       print(f"stations: {grid.n_stations}")
       print(f"x extent: {grid.x_nodes[0]:.1f} to {grid.x_nodes[-1]:.1f} m")
       print(f"z extent: {grid.z_nodes[0]:.1f} to {grid.z_nodes[-1]:.1f} m")
       print(f"rho range: {np.nanmin(grid.resistivity):.3g} to "
             f"{np.nanmax(grid.resistivity):.3g} ohm m")
       print(f"padding cells: {grid.n_pad}")

   describe_grid2d(grid)

Response Containers And Feature Arrays
--------------------------------------

Response containers preserve the physical arrays and also expose feature helpers
for machine-learning and batch processing.

.. list-table::
   :header-rows: 1
   :widths: 26 34 40

   * - Response object
     - Physical shape
     - Feature shape
   * - ``ForwardResponse``
     - 1-D arrays over frequency or time.
     - ``(n_features,)``.
   * - ``ForwardResponse2D``
     - ``(n_freqs, n_stations)`` arrays.
     - ``(n_stations, n_features)``.
   * - ``ForwardResponse3D``
     - ``(n_freqs, n_stations)`` arrays.
     - ``(n_stations, n_features)``.

Use physical arrays when preparing inversion data, plotting scientific
responses, or checking units. Use feature arrays when training models, comparing
many station responses, or building tabular downstream datasets.

Debugging Solver Experiments
----------------------------

When a forward response looks wrong, debug in this order:

#. Run a halfspace model first. A halfspace should produce smooth, stable
   responses.
#. Reduce the number of frequencies and stations until the experiment is easy
   to inspect.
#. Plot the grid with stations visible and confirm that the target is where you
   think it is.
#. Check whether the anomaly is inside the sensitivity range of the selected
   frequencies or time gates.
#. Increase padding and model extent. If the response changes strongly, the
   original grid was too small.
#. Compare a 2-D halfspace result with a 1-D MT response for the same
   resistivity.
#. Check array orientation before passing data to inversion or ML routines.

Next Pages
----------

* :doc:`configuration` explains how configuration objects create grids and
  solvers.
* :doc:`synthetic_datasets` explains how to generate many forward responses.
* :doc:`plotting` shows how to inspect models and responses.
* :doc:`forward_to_inversion` explains how to pass synthetic responses to
  inversion workflows.
