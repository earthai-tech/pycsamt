pycsamt.session
===============

Workflow sessions and normalization helpers for reproducible conversion work.

The :mod:`pycsamt.session` module provides a small provenance layer around
common file-conversion workflows. It is useful when a script, notebook, or
agent workflow produces intermediate objects that should be recorded in a
manifest without turning the session itself into the science API.

Use sessions for workflow bookkeeping. Use the science modules for scientific
work:

* load EDI survey data with :func:`pycsamt.emtools.ensure_sites`;
* run quality control with :func:`pycsamt.emtools.qc.build_qc_table`;
* estimate and apply static shift with :mod:`pycsamt.emtools.ss`;
* use :mod:`pycsamt.agents` when an agent should orchestrate a workflow.

Public Objects
--------------

.. automodule:: pycsamt.session
   :members:
   :show-inheritance:

Concepts
--------

``Session``
   Context manager for recording workflow outputs in a registry manifest.
   It can automatically capture selected conversion results from
   ``to_edi`` and the AVG/J-to-EDI transformers.

``work_session``
   Convenience constructor for :class:`pycsamt.session.Session`.

``Normalize``
   Context manager that accepts heterogeneous sources such as AVG, Jones J,
   EDI files, EDI-like objects, and collections, then returns an EDI-like
   collection whenever possible.

``normalize_session``
   Convenience constructor for :class:`pycsamt.session.Normalize`.

Recommended Pattern
-------------------

Use :class:`~pycsamt.session.Normalize` at the boundary of a workflow, then
pass the normalized EDI output to the real processing API.

.. code-block:: python

   from pycsamt.session import normalize_session
   from pycsamt.emtools import ensure_sites
   from pycsamt.emtools.qc import build_qc_table

   input_path = "data/AMT/WILLY_DATA/L18"
   work_dir = "work/willy_l18"

   with normalize_session(work_dir) as nz:
       edi_collection = nz.load(input_path)

   sites = ensure_sites(edi_collection, recursive=True, verbose=0)
   qc = build_qc_table(sites)

This keeps responsibilities clear: normalization and provenance live in the
session layer, while survey loading and QC live in the EM tools.

Record Conversion Outputs
-------------------------

Use :class:`~pycsamt.session.Session` when you want conversion outputs to be
registered automatically.

.. code-block:: python

   from pycsamt.session import work_session
   from pycsamt.zonge.avg import AVG
   from pycsamt.transformers.jedi import AVGtoEDI

   with work_session("work/avg_conversion") as ses:
       avg = AVG.from_file("raw/S01.AVG")
       edi = AVGtoEDI().transform(avg)

       records = ses.reg.find(tag="AVGtoEDI.transform")

The registry manifest is saved when the context exits. The captured result is
the transformer output, not the original input. If the raw input also needs to
be tracked, register it explicitly through ``ses.reg``.

Explicit Registry Entries
-------------------------

Automatic capture is helpful, but explicit registration is clearer when the
object is important to the workflow.

.. code-block:: python

   from pycsamt.session import Session
   from pycsamt.seg.edi import EDIFile

   with Session("work/manual_registry", auto_capture=False) as ses:
       edi = EDIFile.from_file("data/site001.edi")
       ses.reg.add_object(edi, tags=["raw", "edi"])
       ses.reg.save()

Disable automatic capture when a script should control exactly what appears in
the manifest.

Capture Children From Collections
---------------------------------

By default a collection result is registered as one object. Set
``capture_children=True`` to also register iterable children up to
``max_children``.

.. code-block:: python

   from pycsamt.session import work_session
   from pycsamt.transformers.jedi import JtoEDI
   from pycsamt.jones.collection import JCollection

   with work_session(
       "work/j_batch",
       capture_children=True,
       max_children=32,
   ) as ses:
       j_files = JCollection("raw/jones")
       edi_collection = JtoEDI().transform(j_files)

       batch_records = ses.reg.find(tag="JtoEDI.transform")

Use this mode when downstream audit needs station-level conversion records.
Keep ``max_children`` bounded for large surveys.

Topography During Normalization
-------------------------------

When converting AVG sources, :class:`~pycsamt.session.Normalize` can attach a
topography source before conversion.

.. code-block:: python

   from pycsamt.session import normalize_session

   topo = "data/topography.csv"

   with normalize_session("work/with_topography", topo_src=topo) as nz:
       edi_collection = nz.load("raw/line18.avg")

If the AVG object supports ``add_topography``, that method is used. Otherwise,
the normalizer attempts a conservative fallback for objects exposing a
``frame`` attribute.

Package Root Shortcuts
----------------------

The session helpers are also available lazily from the package root:

.. code-block:: python

   import pycsamt as pcs

   with pcs.work_session("work/root_shortcut") as ses:
       records = ses.reg.list()

This import style is convenient in notebooks, while direct imports from
``pycsamt.session`` are preferred in library code.

When Not To Use It
------------------

Do not use :class:`~pycsamt.session.Session` as a replacement for:

* a ``Sites`` object;
* an EDI collection;
* the desktop or web application session state;
* project memory in :mod:`pycsamt.assistant.memory`;
* workflow agents in :mod:`pycsamt.agents`.

The session module is intentionally small. It records and normalizes workflow
artifacts; it does not own electromagnetic interpretation.
