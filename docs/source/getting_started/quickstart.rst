.. _quickstart:

Quick Start
===========

This page is the shortest practical tour of pyCSAMT v2. It starts with a
folder of EDI files, loads it into the public survey model, builds quality
tables and figures, runs an optional static-shift correction, records a
workflow session, and points toward inversion and AI workflows.

The examples use real v2 entry points. Replace paths such as
``data/AMT/WILLY_DATA/L18`` with your own survey directory.

Install And Verify
------------------

Install the package, then confirm that Python and the command-line interface
can see it.

.. code-block:: console

   python -m pip install pycsamt
   python -c "import pycsamt; print(pycsamt.__version__)"
   pycsamt --help

For source checkouts or optional GUI/AI dependencies, use
:doc:`installation` before continuing.

Load A Survey
-------------

For most users, the public reader is the best first entry point:

.. code-block:: python

   from pycsamt.api import read_edis

   edi_dir = "data/AMT/WILLY_DATA/L18"
   survey = read_edis(
       edi_dir,
       recursive=True,
       strict=False,
       on_dup="replace",
       verbose=0,
   )

   print(survey)
   print(survey.summary().df.head())

``read_edis`` returns an ``APISurvey`` view, which is convenient for display,
summaries, and API-style table outputs.

When a workflow needs the science ``Sites`` container directly, use
``ensure_sites`` from :mod:`pycsamt.emtools`:

.. code-block:: python

   from pycsamt.emtools import ensure_sites

   sites = ensure_sites(edi_dir, recursive=True, strict=False, verbose=0)
   print(len(list(sites)))

``ensure_sites`` accepts EDI directories, EDI files, EDI-like objects,
``EDICollection`` instances, and existing ``Sites`` objects.

Inspect Data Quality
--------------------

Quality control starts with compact tables. Use these before applying
corrections or preparing inversions.

.. code-block:: python

   from pycsamt.emtools.qc import (
       build_qc_table,
       station_confidence_table,
       frequency_confidence_table,
       qc_flags,
   )

   qc = build_qc_table(sites)
   stations = station_confidence_table(sites)
   freqs = frequency_confidence_table(sites)
   flags = qc_flags(sites, min_frac_ok=0.6, min_snr_med=2.0)

   print(qc.head())
   print(stations.head())
   print(freqs.head())
   print(flags.head())

Add quick-look figures when you want to inspect the profile visually:

.. code-block:: python

   from pathlib import Path
   from pycsamt.emtools.qc import (
       plot_frequency_confidence_psection,
       plot_station_confidence_dashboard,
   )

   out = Path("results/quickstart")
   out.mkdir(parents=True, exist_ok=True)

   fig = plot_frequency_confidence_psection(sites).get_figure()
   fig.savefig(out / "frequency_confidence.png", dpi=150)

   fig = plot_station_confidence_dashboard(sites).get_figure()
   fig.savefig(out / "station_confidence.png", dpi=150)

Run A Static-Shift Pass
-----------------------

The AMA static-shift workflow estimates station factors, optionally applies
them, and produces before/after figures. Factors must be finite and positive;
an empty factor table can be valid when the data do not satisfy the method's
assumptions.

.. code-block:: python

   from pycsamt.emtools.ss import (
       estimate_ss_ama,
       correct_ss_ama,
   )

   ss_table = estimate_ss_ama(
       sites,
       half_window=3,
       sort_by="lon",
       verbose=0,
   )
   print(ss_table)

   corrected = correct_ss_ama(
       sites,
       half_window=3,
       sort_by="lon",
       inplace=False,
       verbose=0,
   )

   print(len(list(corrected)))

Use ``inplace=False`` for exploration. Switch to ``inplace=True`` only when
you explicitly want to mutate the loaded objects in memory. For the full
static-shift figure workflow, continue with
:doc:`../tutorials/correct_static_shift`.

Normalize Mixed Inputs
----------------------

If your workflow may receive AVG, Jones J, EDI files, an ``EDICollection``, or
EDI-like objects, normalize at the boundary and then continue with
``ensure_sites``.

.. code-block:: python

   from pycsamt.session import normalize_session
   from pycsamt.emtools import ensure_sites

   with normalize_session("work/quickstart") as nz:
       edi_like = nz.load("raw_or_converted_input")

   sites = ensure_sites(edi_like, recursive=True, strict=False, verbose=0)

Use ``topo_src=...`` with ``normalize_session`` when AVG sources should receive
topography before conversion.

Record A Workflow Session
-------------------------

Use :class:`pycsamt.session.Session` when conversion results should be recorded
in a manifest. This is workflow bookkeeping, not the science data model.

.. code-block:: python

   from pycsamt.session import work_session

   with work_session("work/quickstart") as ses:
       ses.reg.add_object(sites, tags=["raw", "quickstart"])
       records = ses.reg.list()

   print(records)

The manifest is saved on context exit. You can import the same helpers from
the package root as ``pycsamt.work_session`` and ``pycsamt.normalize_session``.

Run A Pipeline Preset
---------------------

Pipelines chain common processing steps and write structured outputs. The
``basic_qc`` preset is a good first run; ``publication_ready`` is more
aggressive and should be reviewed carefully.

.. code-block:: python

   from pycsamt.pipeline import Pipeline, preset_catalogue

   print(preset_catalogue())

   pipe = Pipeline.from_preset("basic_qc")
   result = pipe.run(sites, outdir="results/quickstart_pipeline")

   print(result.summary())

Once this succeeds, inspect the written report and figures before moving to a
stronger preset.

Prepare Occam2D Inputs
----------------------

Occam2D preparation is separate from running the external Occam2D executable.
The first step is to build an ``OccamDataFile.dat`` from loaded sites.

.. code-block:: python

   from pathlib import Path
   from pycsamt.models.occam2d import OccamData

   occam_dir = Path("results/occam2d_l18")
   occam_dir.mkdir(parents=True, exist_ok=True)

   data = OccamData.from_edi(
       sites,
       modes=["RhoTM", "PhsTM"],
       title="L18 quickstart",
   )
   data.write(occam_dir / "OccamDataFile.dat")

To build the full Occam2D working directory, use the higher-level builder:

.. code-block:: python

   from pycsamt.models.occam2d import InputBuilder

   InputBuilder(sites, workdir=occam_dir).build()

Run Occam2D only when the executable is installed and configured:

.. code-block:: python

   from pycsamt.models.occam2d import OccamRunner, InversionResult

   runner = OccamRunner(occam_dir)
   runner.run(target_misfit=1.0)

   occam_result = InversionResult(occam_dir)
   occam_result.plot_model()

If the executable is unavailable, the preparation step is still useful because
it validates station ordering, offsets, modes, and data/error content.

AI-Based 1-D Inversion
----------------------

AI inversion needs a training dataset or a pretrained checkpoint. For a small
synthetic demonstration:

.. code-block:: python

   import numpy as np
   from pycsamt.forward.batch import generate_dataset
   from pycsamt.forward import MT1DForward, LayeredModel
   from pycsamt.ai.inversion.inv1d import EMInverter1D

   ds = generate_dataset(n_samples=2_000, seed=0, n_layers=5)

   inv = EMInverter1D(arch="resnet", n_layers=5, solver="mt1d")
   inv.fit(ds, epochs=30, batch_size=128, verbose=True)

   model = LayeredModel.random(n_layers=5, seed=99)
   response = MT1DForward(np.logspace(-3, 4, 30)).run(model)
   predicted_model = inv.predict_response(response)

For real studies, train on a synthetic distribution that reflects your local
geology, frequency range, noise level, and survey type.

CLI Equivalents
---------------

Use the CLI for quick checks before writing a script:

.. code-block:: console

   pycsamt --help
   pycsamt edi info data/AMT/WILLY_DATA/L18
   pycsamt edi validate data/AMT/WILLY_DATA/L18

The Python API is still preferred for reproducible notebooks, processing
pipelines, and inversion preparation.

Where To Go Next
----------------

* :doc:`first_survey` for a slower walkthrough of loading and QC;
* :doc:`data_formats` for AVG, Jones J, spectra, TEM/TDEM, and EDI details;
* :doc:`../tutorials/inspect_and_qc_survey` for richer diagnostics;
* :doc:`../tutorials/correct_static_shift` for the full static-shift workflow;
* :doc:`../pipeline/index` for custom processing chains;
* :doc:`../tutorials/prepare_occam2d_inversion` for inversion preparation.
