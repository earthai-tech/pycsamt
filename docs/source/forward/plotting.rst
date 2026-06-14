.. _forward_plotting:

Forward Plotting
================

Forward plots are quality-control tools. They help users see whether a model,
grid, response, or dataset looks physically plausible before the result is
used for inversion or AI training.

1-D Response And Model Plots
----------------------------

.. code-block:: python
   :linenos:

   import numpy as np

   from pycsamt.forward import (
       LayeredModel,
       MT1DForward,
       plot_response_1d,
       plot_model_1d,
       plot_response_and_model_1d,
   )

   model = LayeredModel([100.0, 10.0, 500.0], [300.0, 800.0])
   response = MT1DForward(np.logspace(-3, 4, 40)).run(model)

   axes = plot_response_1d(response, title="MT1D response")
   ax = plot_model_1d(model, title="Layered model")
   fig = plot_response_and_model_1d(response, model=model)

2-D Model And Pseudo-Sections
-----------------------------

.. code-block:: python
   :linenos:

   from pycsamt.forward import (
       Grid2D,
       MT2DForward,
       plot_model_2d,
       plot_pseudosection_2d,
       plot_response_profiles,
   )

   grid = Grid2D.with_anomaly(
       bg_rho=500.0,
       anomaly_rho=5.0,
       anomaly_bounds=(2000.0, 6000.0, 300.0, 1500.0),
       nx=50,
       nz=35,
       x_max=10_000.0,
       z_max=6000.0,
       n_stations=16,
   )
   response = MT2DForward(freqs=[1.0, 10.0, 100.0], grid=grid).run()

   plot_model_2d(grid, title="2-D anomaly model")
   plot_pseudosection_2d(response, mode="te", quantity="rho_a")
   plot_response_profiles(response, mode="both")

3-D Maps And Tensor Components
------------------------------

.. code-block:: python
   :linenos:

   from pycsamt.forward import (
       Grid3D,
       MT3DForward,
       plot_model_3d,
       plot_response_map_3d,
       plot_response_section_3d,
       plot_tensor_components_3d,
   )

   grid = Grid3D.halfspace(
       rho=100.0,
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

   plot_model_3d(grid)
   plot_response_map_3d(response, freq_idx=0, component="xy")
   plot_response_section_3d(response, component="xy")
   plot_tensor_components_3d(response, freq_idx=0)

Plotting Checklist
------------------

Use plots to check:

* apparent resistivity and phase are smooth for simple 1-D models;
* conductive anomalies produce plausible low-resistivity response zones;
* station positions align with the model;
* padding cells are not dominating the displayed interpretation region;
* noisy synthetic samples still look like possible field data;
* tensor components and TE/TM modes are labelled consistently.
