.. _user_guide_iot_data_bridge:

Bridging IoT Acquisition and the Data Model
===========================================

The IoT layer processes raw arrays and telemetry with numpy only and stays
independent of the heavier science API, which keeps ``import pycsamt.iot``
cheap on a field gateway. The :mod:`pycsamt.iot.bridge` module is the one
place the two worlds meet. Every science-API import (``seg``, ``site``) is
lazy, so pulling in the bridge does not force the geospatial stack on a
node that only needs telemetry.

The bridge works in both directions.

Forward: acquisition to the data model
--------------------------------------

The edge already forms per-window impedance estimates (see
:func:`~pycsamt.iot.assess_impedance_stability`).
:func:`~pycsamt.iot.impedance_to_z` turns them into a
:class:`pycsamt.z.z.Z` -- the impedance-tensor container the processing and
inversion flow consumes -- deriving an absolute error from the spread
across windows:

.. code-block:: python

   import numpy as np
   from pycsamt.iot import impedance_to_z, z_to_edi

   freq = np.logspace(4, 0, 12)
   z_windows = ...            # shape (n_windows, n_freq) complex, from the edge
   z = impedance_to_z(z_windows, freq, station="S01", method="amt")

:func:`~pycsamt.iot.z_to_edi` (and the in-memory
:func:`~pycsamt.iot.build_edifile`) writes a preliminary EDI from that
tensor, so a field node can emit an ``.edi`` on site:

.. code-block:: python

   path = z_to_edi(z, station="S01", lat=6.5, lon=3.4, elevation=120.0,
                   savepath="edi_out", method="amt")

For a whole session, :meth:`~pycsamt.iot.FieldSession.to_edifiles` builds
one EDI per station, enriched with the station geometry the session
recorded, and :meth:`~pycsamt.iot.FieldSession.to_sites_collection` returns
a :class:`pycsamt.site.base.Sites` collection ready to feed straight to
:meth:`pycsamt.pipeline.Pipeline.run`:

.. code-block:: python

   sites = session.to_sites_collection({"S01": z1, "S02": z2})
   # sites can now be handed to pycsamt.pipeline.Pipeline.run(sites)

When the optional geospatial stack is available,
:func:`~pycsamt.iot.z_to_site` returns a single EDI-backed
:class:`pycsamt.site.base.Site`.

Consistent QC via emtools
-------------------------

Field-side edge QC and downstream processing QC should not be two
independent notions of a "good" station. :func:`~pycsamt.iot.emtools_qc`
routes the sites the IoT layer produced through the *same*
:mod:`pycsamt.emtools` coherence/skew/SNR quality control the processing
flow uses:

.. code-block:: python

   from pycsamt.iot import emtools_qc

   table = emtools_qc(session, {"S01": z1, "S02": z2})       # full QC table
   flags = emtools_qc(session, {"S01": z1}, flags=True)      # pass/flag verdicts

It also accepts an already-built ``Sites`` collection or any EDI source, so
the same QC can be re-run on an archived survey.

.. note::

   For a *raw time series* rather than pre-computed impedance, use
   :func:`pycsamt.ts.ts_to_z` / :func:`pycsamt.ts.ts_to_edi`, which run the
   spectral estimation for you. The bridge starts one step later, from the
   impedance the edge has already estimated.

Reverse: seeding a re-occupation from an EDI survey
---------------------------------------------------

The bridge also reads an existing survey to plan a follow-up campaign.
:func:`~pycsamt.iot.field_session_from_edis` returns a
:class:`~pycsamt.iot.FieldSession` with every station's recorded geometry
and channels, plus a sensor node per station, ready to re-occupy:

.. code-block:: python

   from pycsamt.iot import field_session_from_edis, edi_survey_table

   session = field_session_from_edis("survey_2024/", survey_id="REOCCUPY",
                                     method="amt")
   session.n_stations, session.n_devices

   edi_survey_table("survey_2024/")   # station geometry + frequency coverage

:func:`~pycsamt.iot.deployment_from_edis` is a lighter counterpart that
returns just a ``DeploymentConfig`` inventory, and
:func:`~pycsamt.iot.read_edi_survey` yields the raw per-station summary
records.

Sources may be a directory of ``.edi`` files, a glob pattern, a single
file or ``EDIFile``, or any iterable mixing these.
