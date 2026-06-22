.. _user_guide_inversion:

Inversion
=========

``pycsamt.inversion`` provides one backend-neutral API for EM inversion across
MT, AMT, CSAMT, EMAP, and TDEM workflows. The package separates four concerns:

* ``InversionConfig`` describes the survey, backend, model, errors, and solver.
* ``InversionWorkflow`` runs the selected backend.
* ``InversionResult`` stores model, diagnostics, uncertainty, history, and files.
* ``pycsamt.inversion.export`` writes common products from any result.

Backends
--------

The current backend names are:

``builtin``
   Dependency-light local inversions for MT/AMT/CSAMT layered 1-D, TDEM 1-D,
   stitched 2-D profiles, and opt-in finite-difference 2-D profile inversion.

``simpeg``
   Optional SimPEG adapter for natural-source EM inversion paths. It is imported
   only when selected.

``pygimli``
   Optional pyGIMLi adapter for 1-D EM modelling/inversion and stitched profiles.

``occam2d``
   Occam2D lifecycle adapter. It prepares and validates native files, assembles
   runner commands, and can launch the executable when ``run_external=True``.

``modem``
   ModEM lifecycle adapter for 3-D and 2-D external workflows. Execution is
   opt-in through ``run_external=True``.

Minimal MT 1-D
--------------

.. code-block:: python

   import numpy as np
   from pycsamt.forward import LayeredModel, MT1DForward
   from pycsamt.inversion import InversionConfig, InversionWorkflow, StartingModel

   freqs = np.logspace(-2, 2, 12)
   truth = LayeredModel([80.0, 25.0, 600.0], [250.0, 900.0])
   response = MT1DForward(freqs=freqs).run(truth)

   cfg = InversionConfig(
       method="mt",
       dimension="1d",
       backend="builtin",
       data={"freqs": freqs, "rho_a": response.rho_a, "phase": response.phase},
       starting_model=StartingModel([100.0, 50.0, 500.0], [300.0, 1000.0]),
       max_iter=12,
   )
   result = InversionWorkflow(cfg).run()
   print(result.summary())

Minimal TDEM 1-D
----------------

.. code-block:: python

   import numpy as np
   from pycsamt.forward import LayeredModel, TEM1DForward
   from pycsamt.inversion import InversionConfig, InversionWorkflow, StartingModel

   times = np.logspace(-5, -3, 7)
   truth = LayeredModel([60.0, 250.0, 900.0], [120.0, 700.0])
   forward_options = {"loop_radius": 25.0, "n_freqs": 10, "n_lam": 16}
   response = TEM1DForward(times=times, **forward_options).run(truth)

   cfg = InversionConfig(
       method="tdem",
       dimension="1d",
       backend="builtin",
       data={"times": times, "values": response.dBz_dt},
       starting_model=StartingModel([80.0, 200.0, 700.0], [150.0, 800.0]),
       backend_options=forward_options,
       max_iter=8,
   )
   result = InversionWorkflow(cfg).run()

Stitched 2-D Profile
--------------------

For a profile where each station has a frequency response, pass ``rho_a`` and
``phase`` as arrays with shape ``(n_stations, n_frequencies)``.

.. code-block:: python

   cfg = InversionConfig(
       method="mt",
       dimension="2d",
       backend="builtin",
       data={
           "freqs": freqs,
           "rho_a": rho_by_station,
           "phase": phase_by_station,
           "station_x": [0.0, 400.0, 800.0],
           "station_names": ["S00", "S01", "S02"],
       },
       starting_model=StartingModel([100.0, 50.0, 600.0], [300.0, 900.0]),
       max_iter=10,
   )
   result = InversionWorkflow(cfg).run()
   model = result.to_resistivity_model()

Use the finite-difference 2-D built-in path when you want real 2-D physics
instead of stitched station inversions:

.. code-block:: python

   cfg.backend_options.update({
       "profile_mode": "fd2d",
       "nx": 12,
       "n_pad": 2,
       "components": ("te", "tm"),
   })

Exports
-------

All common export functions accept an ``InversionResult``:

.. code-block:: python

   from pycsamt.inversion import export, plot

   export.to_csv(result, "model.csv")
   export.to_npz(result, "model.npz")
   export.to_geojson(result, "model.geojson")
   export.to_vtk(result, "model.vtk")
   export.to_archive(result, "snapshot.zip")

   # Requires rasterio.
   export.to_geotiff(result, "model.tif")

   plot.plot_model(result)
   plot.plot_rms(result)

Hydrogeophysical Handoff
------------------------

Any inversion result that can be converted to ``ResistivityModel`` can enter the
``pycsamt.interp`` hydrogeophysical workflow:

.. code-block:: python

   from pycsamt.interp import EMHydroModel, PetrophysicalConfig
   from pycsamt.interp.petrophysics import ArchieModel

   resistivity_model = result.to_resistivity_model()
   hydro_cfg = PetrophysicalConfig(
       petro=ArchieModel(m=1.8, n=2.0),
       rho_w=20.0,
       porosity_prior=0.28,
   )
   hydro = EMHydroModel(resistivity_model, hydro_cfg, method_tag="MT").fit()
   print(hydro.water_table)

Runnable Demo
-------------

The complete synthetic demo lives at ``examples/demo/demo_inversion.py``:

.. code-block:: bash

   python examples/demo/demo_inversion.py

It runs MT 1-D, TDEM 1-D, stitched 2-D, writes CSV/NPZ/GeoJSON/VTK/archive
products, saves quick-look plots, and passes the recovered 2-D section into the
hydrogeophysical interpretation API.
