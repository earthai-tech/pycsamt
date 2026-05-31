.. _quickstart:

Quick start
===========

This page demonstrates the core v2 workflow end to end.

Processing EDI data
-------------------

.. code-block:: python

   from pycsamt.io import EDICollection
   from pycsamt.emtools import impedance

   # Load a folder of EDI files
   coll = EDICollection.from_dir("data/edi/")
   Z = impedance.apparent_resistivity(coll)

2-D inversion with Occam2D
--------------------------

.. code-block:: python

   from pycsamt.models.occam2d import OccamData, OccamRunner, InversionResult

   data = OccamData.from_edi("data/edi/", mode="TM")
   data.write("occam_run/OccamDataFile.dat")

   runner = OccamRunner(workdir="occam_run/")
   runner.run()

   result = InversionResult("occam_run/")
   result.plot_model()

Geological interpretation
--------------------------

.. code-block:: python

   from pycsamt.interp import ResistivityModel, ModelCalibrator
   from pycsamt.interp import export

   model = ResistivityModel.from_occam2d(result)
   cal   = ModelCalibrator(ptol=0.10).fit(model)
   logs  = cal.stratigraphic_logs()

   export.to_oasis_montaj_xyz(logs, "profile_K1.xyz")
   export.to_las(logs[0], "S17.las")

AI-based 1-D inversion
-----------------------

.. code-block:: python

   from pycsamt.ai.inversion import EMInverter1D

   inv = EMInverter1D(n_layers=5, arch="fcn")
   inv.fit(X_train, y_train, epochs=50)
   y_pred = inv.predict(X_test)
