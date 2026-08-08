.. _forward_maxwell_overview:

Maxwell Adapter Layer
=====================

:doc:`solvers_and_grids` covers ``MT1DForward``, ``MT2DForward``, and
``MT3DForward``: fast, in-process solvers over ``LayeredModel``, ``Grid2D``,
and ``Grid3D`` with no formal validation harness of their own beyond the
tests that ship with pycsamt. :mod:`pycsamt.forward.maxwell` sits beside
that layer, not above it -- it wraps some of those exact same solvers (and
two production external ones) behind a solver-neutral contract that checks
a problem's compatibility before solving, validates a backend's result
after solving, and can prove the whole chain against a closed-form answer.
The six pages under this heading cover that layer in depth:
:doc:`maxwell_contracts`, :doc:`maxwell_meshing`, :doc:`maxwell_backends`,
:doc:`maxwell_adapters`, :doc:`maxwell_benchmarks`, and
:doc:`maxwell_caching_and_batch`. This page is the map between them, and the
answer to the first real question: when does any of this matter more than
calling ``MT2DForward`` directly?

Legacy And Validated Engines
----------------------------

.. list-table::
   :header-rows: 1
   :widths: 26 38 36

   * - Layer
     - Main responsibility
     - Typical use
   * - ``pycsamt.forward.em1d``/``em2d``/``em3d``
     - Fast, in-process 1-D/2-D/quasi-3-D solvers over ``LayeredModel``,
       ``Grid2D``, and ``Grid3D``. No capability declaration, benchmark
       record, or cache of their own.
     - Prototyping, survey design, and generating a
       :term:`synthetic dataset` at scale through
       :doc:`synthetic_datasets`'s ``ForwardDataset``.
   * - :mod:`pycsamt.forward.maxwell`
     - Solver-neutral problem/result contracts, capability-checked
       adapters over both in-repository research solvers and production
       external ones (:term:`ModEM`, :term:`MARE2DEM`), analytic
       benchmarks, content-addressed caching, and resumable batch solving.
     - A response whose accuracy is a checked claim rather than an
       assumption; AI-inversion training data; a production-grade forward
       check before committing to a full external inversion project.

Both layers can answer the same physical question for the same model --
:doc:`maxwell_adapters`' ``MT2DAdapter`` literally wraps ``MT2DForward``.
What ``pycsamt.forward.maxwell`` adds is not new physics; it is a validation
boundary, a benchmark record, and machinery for running that boundary at
scale. Reach for the legacy solvers directly when a quick, throwaway check
is enough. Reach for this layer when the result needs to be defensible on
its own -- reused in a report, fed into training data, or compared across
solver versions.

Package Map
-----------

.. list-table::
   :header-rows: 1
   :widths: 26 34 20

   * - Concern
     - Modules
     - Page
   * - Problem, mesh, and result contracts
     - ``contracts.py``, ``contracts_tri.py``
     - :doc:`maxwell_contracts`
   * - Structured and triangular mesh building
     - ``mesh.py``, ``tri_mesh_gen.py``
     - :doc:`maxwell_meshing`
   * - Backend capability and registry
     - ``backends.py``
     - :doc:`maxwell_backends`
   * - Adapter validation boundary and five concrete adapters
     - ``adapters.py``, ``external.py``, ``mt2d.py``, ``mt3d.py``,
       ``tri_fem2d.py``, ``modem3d.py``, ``mare2dem.py``
     - :doc:`maxwell_adapters`
   * - Analytic benchmarks
     - ``benchmarks.py``
     - :doc:`maxwell_benchmarks`
   * - Result caching and batch solving
     - ``cache.py``, ``batch.py``
     - :doc:`maxwell_caching_and_batch`

Every module above is reachable from :mod:`pycsamt.forward.maxwell`
directly except the five concrete adapter modules and ``external.py``,
which are imported from their own submodules -- importing the package
itself never pulls in an optional solver or an external executable, only
the contracts and registry machinery that describe them.

Core Workflow
-------------

The clearest way to see that this layer adds a boundary rather than a
second solver is to run the identical model both ways and compare the raw
impedance. ``MT2DAdapter`` builds the exact same
:class:`~pycsamt.forward.grid2d.Grid2D` internally that a direct
``MT2DForward`` call would use, so the two paths have nothing left to
disagree about:

.. code-block:: pycon

   >>> import numpy as np
   >>> from pycsamt.forward.grid2d import Grid2D
   >>> from pycsamt.forward.em2d import MT2DForward
   >>> from pycsamt.forward.maxwell import MaxwellMesh, MaxwellProblem, ReceiverSet
   >>> from pycsamt.forward.maxwell.mt2d import MT2DAdapter

   >>> x_edges = np.linspace(0, 10_000, 41)
   >>> z_edges = np.linspace(0, 5_000, 31)
   >>> resistivity = np.full((30, 40), 100.0)
   >>> resistivity[8:14, 15:26] = 8.0
   >>> frequencies = np.geomspace(300, 1, 8)
   >>> station_x = np.linspace(2_000, 8_000, 5)

   >>> grid = Grid2D(
   ...     dx=np.diff(x_edges), dz=np.diff(z_edges), resistivity=resistivity,
   ...     x_stations=station_x, n_pad=0, name="overview-demo",
   ... )
   >>> legacy = MT2DForward(frequencies, grid, verbose=False).run()

   >>> mesh = MaxwellMesh(x_edges, z_edges)
   >>> receivers = ReceiverSet(
   ...     [[x, 0.0] for x in station_x], [f"S{i:02d}" for i in range(5)],
   ... )
   >>> problem = MaxwellProblem(
   ...     mesh, 1.0 / resistivity, frequencies, receivers, ("zxy", "zyx"),
   ... )
   >>> result = MT2DAdapter(verbose=False).solve(problem)
   >>> zxy_maxwell = result.impedance_v_a[:, :, 0].T
   >>> zyx_maxwell = result.impedance_v_a[:, :, 1].T
   >>> float(np.max(np.abs(legacy.zxy - zxy_maxwell)))
   0.0
   >>> float(np.max(np.abs(legacy.zyx - zyx_maxwell)))
   0.0

Zero, not merely small -- the two arrays are computed by the same finite-
difference code underneath. What ``result`` carries that ``legacy`` does
not is a :term:`problem hash` tying it to this exact input, a
:term:`solver diagnostics` record, and a backend identity that passed
through :term:`preflight assessment` and :term:`postflight validation`
before ever reaching this comparison. The executed apparent-resistivity and
phase figure for this same shallow-conductor setup already lives in
:doc:`maxwell_adapters`; nothing about rerunning it here would show
anything new, which is itself evidence for the point being made -- the
physics did not change, only what is checked around it.

Choosing A Path
---------------

.. list-table::
   :header-rows: 1
   :widths: 34 26 30

   * - Goal
     - Start here
     - Then read
   * - Fast, approximate response for prototyping or survey design
     - :doc:`solvers_and_grids`
     - :doc:`synthetic_datasets`, :doc:`plotting`
   * - A response with a checkable accuracy claim
     - :doc:`maxwell_adapters`
     - :doc:`maxwell_benchmarks`
   * - Understand the problem/mesh/result data model
     - :doc:`maxwell_contracts`
     - :doc:`maxwell_meshing`
   * - Build a mesh from real geology or topography
     - :doc:`maxwell_meshing`
     - :doc:`maxwell_adapters`
   * - Pick or register a backend by capability, not by import
     - :doc:`maxwell_backends`
     - :doc:`maxwell_adapters`
   * - Solve many problems reliably, or reuse results across runs
     - :doc:`maxwell_caching_and_batch`
     - :doc:`maxwell_benchmarks`
   * - Generate AI-inversion training data from this layer
     - :doc:`/user_guide/ai_inversion/forward_physics`
     - :doc:`/user_guide/ai_inversion/dataset2d`,
       :doc:`/user_guide/ai_inversion/roadmap`
   * - Run a full external inversion project, not just a forward check
     - :doc:`/user_guide/models/overview`
     - :doc:`/user_guide/models/modem`, :doc:`/user_guide/models/mare2dem`

Relationship To Theory
----------------------

:doc:`/theory/maxwell_forward` derives the physics and the contract
decisions this whole layer is built on -- why the mesh scale follows skin
depth, what a capability declaration is actually claiming, and why
postflight validation has to check axis identity rather than trust a
returned array's shape. Read it when a page here states a rule without
re-deriving it, which is deliberate: these pages point to that derivation
rather than repeating it. :doc:`/user_guide/ai_inversion/forward_physics`
and :doc:`/user_guide/ai_inversion/roadmap` cover the same modules from the
opposite direction -- not "how do I get one validated response" but "how do
these guarantees let a training pipeline trust ten thousand of them
unattended."
