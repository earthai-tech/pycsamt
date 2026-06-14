.. _forward_solvers_and_grids:

Solvers And Grids
=================

The forward package separates model containers from solvers. A model object
stores resistivity and geometry. A solver computes the response at requested
frequencies or time gates.

1-D Layered Models
------------------

:class:`pycsamt.forward.LayeredModel` represents a layered earth. The last
layer is the halfspace, so ``thickness`` has one fewer value than
``resistivity``.

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
   response = MT1DForward(freqs).run(model)

   print(response.rho_a.shape)
   print(response.phase.shape)

1-D Solvers
-----------

.. list-table::
   :header-rows: 1
   :widths: 25 35 40

   * - Solver
     - Input axis
     - Output
   * - ``MT1DForward``
     - Frequencies in Hz
     - Complex impedance, apparent resistivity, phase.
   * - ``CSAMT1DForward``
     - Frequencies in Hz and optional source offset
     - MT-like response with first-order controlled-source correction when
       ``source_offset`` is supplied.
   * - ``TEM1DForward``
     - Time gates in seconds and loop radius
     - Step-off :math:`\partial B_z / \partial t` response.

.. code-block:: python
   :linenos:

   import numpy as np

   from pycsamt.forward import CSAMT1DForward, TEM1DForward

   csamt = CSAMT1DForward(
       freqs=np.logspace(0, 4, 30),
       source_offset=5000.0,
       dipole_length=1000.0,
   ).run(model)

   tem = TEM1DForward(
       times=np.logspace(-6, -2, 25),
       loop_radius=50.0,
   ).run(model)

2-D Grids
---------

:class:`pycsamt.forward.Grid2D` stores a non-uniform 2-D resistivity model and
surface station positions. It can be built as a halfspace, a block anomaly, a
random model, or a horizontal extension of a 1-D layered model.

.. code-block:: python
   :linenos:

   from pycsamt.forward import Grid2D, MT2DForward

   grid = Grid2D.with_anomaly(
       bg_rho=500.0,
       anomaly_rho=5.0,
       anomaly_bounds=(2000.0, 6000.0, 300.0, 1500.0),
       nx=50,
       nz=35,
       x_max=10_000.0,
       z_max=6000.0,
       n_pad=8,
       n_stations=16,
   )

   solver = MT2DForward(freqs=[1.0, 10.0, 100.0], grid=grid, verbose=True)
   response = solver.run()

   station_0 = response.station_response(0)
   features = response.to_feature_array(mode="both")

The 2-D response arrays have shape ``(n_freqs, n_stations)`` for TE and TM
apparent resistivity and phase.

3-D Grids
---------

:class:`pycsamt.forward.Grid3D` stores a non-uniform 3-D resistivity volume
and a 2-D station layout. :class:`pycsamt.forward.MT3DForward` returns full
tensor-style response arrays for ``xy``, ``yx``, ``xx``, and ``yy``
components.

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
       nx_stations=5,
       ny_stations=5,
   )

   response = MT3DForward(freqs=[1.0, 10.0, 100.0], grid=grid).run()
   features = response.to_feature_array(components="xy_yx")

Grid Design Checks
------------------

Before trusting a forward response, check:

* all resistivity values are positive;
* station positions lie inside the model extent;
* padding is large enough to reduce boundary influence;
* shallow cells resolve the highest frequencies or earliest times;
* target depths are covered by the model;
* the response quantity matches the intended inversion feature vector.
