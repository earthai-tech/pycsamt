.. _forward_maxwell_backends:

Maxwell Backend Registry
========================

:doc:`maxwell_adapters` walked through five concrete adapters and used
``report = adapter.assess(problem)`` to reject a buried receiver before any
solver ran. That check has two layers: a generic assessment every adapter
inherits from its declared :class:`~pycsamt.forward.maxwell.BackendCapabilities`,
and the solver-specific extras each adapter's own ``assess`` adds on top --
``MT2DAdapter``'s surface-receiver rule, ``Mare2DEMAdapter``'s per-region
conductivity rule, and so on. This page stays below that adapter-specific
layer and covers the two things every backend shares regardless of which
solver it wraps: exactly what the generic assessment checks, and how a
backend gets found by capability instead of by remembering which Python
module implements it.

Capability Assessment
---------------------

:meth:`~pycsamt.forward.maxwell.BackendCapabilities.assess` is the method
every adapter's own ``assess`` calls first, before adding its extra checks.
It never touches the numerical solver -- every check reads only the problem's
shape and the capability declaration -- and it does not stop at the first
failure:

.. code-block:: pycon

   >>> import numpy as np
   >>> from pycsamt.forward.maxwell import (
   ...     BackendCapabilities, MaxwellMesh, MaxwellProblem, ReceiverSet,
   ... )

   >>> mesh = MaxwellMesh(np.linspace(0, 1_000, 6), np.linspace(0, 1_000, 6))
   >>> problem = MaxwellProblem(
   ...     mesh, np.full(mesh.shape, 0.01), [1.0],
   ...     ReceiverSet([[500.0, 0.0]], ["S00"]), ("zxy", "zyx"),
   ... )
   >>> demo3d = BackendCapabilities(
   ...     "demo3d", "1", (3,), ("zxx", "zyy"), time_conventions=("exp(-iwt)",)
   ... )
   >>> report = demo3d.assess(problem)
   >>> report.compatible
   False
   >>> report.errors
   ('2-D problems are unsupported', "unsupported impedance components: ['zxy', 'zyx']", "time convention 'exp(+iwt)' is unsupported")

A 2-D problem asking for ``zxy``/``zyx`` at the default
``exp(+iwt)`` convention fails all three of ``demo3d``'s declared axes at
once, and every failure is collected rather than short-circuited on the
first one -- useful when a script is deciding what to fix rather than just
whether to proceed. :meth:`~pycsamt.forward.maxwell.CompatibilityReport.require`
turns the same errors into one exception when the caller just wants to fail
fast:

.. code-block:: pycon

   >>> report.require()
   Traceback (most recent call last):
   ...
   ValueError: backend 'demo3d' is incompatible: 2-D problems are unsupported; unsupported impedance components: ['zxy', 'zyx']; time convention 'exp(+iwt)' is unsupported

Mesh structure is checked the same way. A capability that declares
``supports_nonuniform_mesh=False`` rejects a mesh with variable cell widths
before ever reaching a solver that assumes a constant grid spacing:

.. code-block:: pycon

   >>> nonuniform_mesh = MaxwellMesh([0, 100, 300, 1_000], [0, 200, 1_000])
   >>> nonuniform_problem = MaxwellProblem(
   ...     nonuniform_mesh, np.full(nonuniform_mesh.shape, 0.01), [1.0],
   ...     ReceiverSet([[150.0, 0.0]], ["S00"]), ("zxy", "zyx"),
   ... )
   >>> uniform_only = BackendCapabilities(
   ...     "demo-uniform", "1", (2,), ("zxy", "zyx"), supports_nonuniform_mesh=False,
   ... )
   >>> uniform_only.assess(nonuniform_problem).errors
   ('nonuniform meshes are unsupported',)

The inactive-cell check is the subtlest one, because it reads two flags
together rather than one. ``supports_inactive_cells`` alone governs whether
:attr:`~pycsamt.forward.maxwell.MaxwellProblem.active_cells` may contain any
``False`` entries at all:

.. code-block:: pycon

   >>> active = np.ones(mesh.shape, dtype=bool)
   >>> active[0, :] = False
   >>> air_problem = MaxwellProblem(
   ...     mesh, np.full(mesh.shape, 0.01), [1.0],
   ...     ReceiverSet([[500.0, 0.0]], ["S00"]), ("zxy", "zyx"), active_cells=active,
   ... )
   >>> no_inactive = BackendCapabilities(
   ...     "demo-no-inactive", "1", (2,), ("zxy", "zyx"), supports_inactive_cells=False,
   ... )
   >>> no_inactive.assess(air_problem).errors
   ('inactive cells are unsupported',)

But a whole flat top row of inactive cells -- the same air mask under every
column -- is a materially simpler claim than an inactive mask that varies
laterally, because a laterally uniform mask is just an air-layer count,
while a laterally varying one is a non-flat air/earth boundary in disguise.
The assessment treats those as different questions: the second one is only
accepted once the capability also declares ``supports_topography=True``,
never merely ``supports_inactive_cells=True``:

.. code-block:: pycon

   >>> lateral_active = np.ones(mesh.shape, dtype=bool)
   >>> lateral_active[0, :2] = False
   >>> ridge_problem = MaxwellProblem(
   ...     mesh, np.full(mesh.shape, 0.01), [1.0],
   ...     ReceiverSet([[500.0, 0.0]], ["S00"]), ("zxy", "zyx"), active_cells=lateral_active,
   ... )
   >>> flat_inactive_ok = BackendCapabilities(
   ...     "demo-inactive-flat", "1", (2,), ("zxy", "zyx"),
   ...     supports_inactive_cells=True, supports_topography=False,
   ... )
   >>> flat_inactive_ok.assess(ridge_problem).errors
   ('laterally varying inactive cells require topography support',)
   >>> topo_capable = BackendCapabilities(
   ...     "demo-topo", "1", (2,), ("zxy", "zyx"),
   ...     supports_inactive_cells=True, supports_topography=True,
   ... )
   >>> topo_capable.assess(ridge_problem).compatible
   True

``TriFEM2DAdapter`` and ``Mare2DEMAdapter`` are exactly the ``supports_topography=True``
case; :doc:`maxwell_adapters` shows what a real laterally-varying terrain
mask looks like once it comes from an actual topography polyline rather
than two flipped array entries. Cell and frequency ceilings are simpler
counts, checked directly against the mesh and axis size:

.. code-block:: pycon

   >>> big_mesh = MaxwellMesh(np.linspace(0, 1_000, 9), np.linspace(0, 1_000, 9))
   >>> big_problem = MaxwellProblem(
   ...     big_mesh, np.full(big_mesh.shape, 0.01), [10.0, 1.0],
   ...     ReceiverSet([[500.0, 0.0]], ["S00"]), ("zxy", "zyx"),
   ... )
   >>> capped = BackendCapabilities(
   ...     "demo-small", "1", (2,), ("zxy", "zyx"), maximum_cells=10, maximum_frequencies=1,
   ... )
   >>> capped.assess(big_problem).errors
   ('cell count 64 exceeds limit 10', 'frequency count 2 exceeds limit 1')

This is exactly the check ``MT3DAdapter``'s default 6,000-cell ceiling uses
in :doc:`maxwell_adapters`. Finally, one check never produces a hard error at
all -- it is advisory by construction, appended as a warning rather than an
error whenever a capability's :term:`verified benchmark` list is empty:

.. code-block:: pycon

   >>> unverified = BackendCapabilities("demo-plain", "1", (2,), ("zxy", "zyx"))
   >>> unverified.assess(problem).compatible, unverified.assess(problem).warnings
   (True, ('backend declares no verified benchmarks',))

This is the same warning ``Mare2DEMAdapter`` carries in its current
capability declaration, shown from the adapter's own side in
:doc:`maxwell_adapters`: a problem can be perfectly *compatible* with a
backend's declared shape while that backend has never been *checked*
against a known answer. Compatibility and verification are independent
questions, and only ``verified_benchmarks`` answers the second one.

The Backend Protocol
--------------------

:class:`~pycsamt.forward.maxwell.MaxwellBackend` is a
:func:`~typing.runtime_checkable` :class:`~typing.Protocol`: any object with
a ``capabilities`` property and a ``solve`` method satisfies it structurally,
without inheriting from any particular base class. Every concrete adapter
from :doc:`maxwell_adapters` conforms; most ordinary objects do not:

.. code-block:: pycon

   >>> from pycsamt.forward.maxwell import MaxwellBackend
   >>> from pycsamt.forward.maxwell.mt2d import MT2DAdapter
   >>> isinstance(MT2DAdapter(), MaxwellBackend)
   True
   >>> isinstance(object(), MaxwellBackend)
   False

   >>> class NoSolve:
   ...     capabilities = BackendCapabilities("nosolve", "1", (2,), ("zxy",))
   >>> isinstance(NoSolve(), MaxwellBackend)
   False

``NoSolve`` declares the right ``capabilities`` attribute but has no
``solve`` method, so it fails the protocol check the same way it would fail
at the first call site that tried to use it as a backend -- structural
typing catches the omission immediately rather than deep inside a batch run.
This is what lets the registry below store arbitrary factories and still
guarantee, before ever calling one, that whatever it builds can be used the
same way regardless of which module defined it.

Registering Backends
--------------------

A :class:`~pycsamt.forward.maxwell.BackendRegistration` bundles three
things: a capability declaration, a zero-argument-or-keyword factory that
builds the backend, and an optional :term:`availability probe`. The last
argument is what lets ``list_backends()`` report whether an external
executable can currently be resolved, without importing or constructing
anything:

.. code-block:: pycon

   >>> from pycsamt.forward.maxwell import BackendRegistration
   >>> demo_cap = BackendCapabilities("registration-demo", "1", (2,), ("zxy",))
   >>> registration = BackendRegistration(demo_cap, lambda **kw: MT2DAdapter())
   >>> registration.availability()
   (True, None)

:meth:`~pycsamt.forward.maxwell.BackendRegistration.create` does one more
check that neither ``assess`` nor the :class:`MaxwellBackend` protocol
covers: it confirms the factory actually built something whose
``capabilities`` match the registration it was filed under, byte for byte.
A factory that quietly returns a different version than it registered is
rejected rather than silently served:

.. code-block:: pycon

   >>> from pycsamt.forward.maxwell import BackendRegistry
   >>> registered_cap = BackendCapabilities("consistency-demo", "1.0", (2,), ("zxy", "zyx"))
   >>> drifted_cap = BackendCapabilities("consistency-demo", "2.0", (2,), ("zxy", "zyx"))
   >>> class DriftedBackend:
   ...     capabilities = drifted_cap
   ...     def solve(self, problem):
   ...         raise NotImplementedError
   >>> drifted_registration = BackendRegistration(registered_cap, lambda **kw: DriftedBackend())
   >>> drifted_registration.create()
   Traceback (most recent call last):
   ...
   RuntimeError: backend instance capabilities differ from its registration.

That failure mode matters in practice whenever a factory closes over a
mutable configuration object: bump the wrapped solver's version after
registering but before the first ``create_backend`` call, and this check is
what turns a silent version mismatch into an immediate, loud error instead
of a mislabeled result three steps later.

The four in-repository and external adapters that ship a
``register_*_backend`` entry point all register into the same process-wide
:data:`~pycsamt.forward.maxwell.backend_registry`, so a caller can ask what
is available without importing any adapter module by name:

.. code-block:: pycon

   >>> from pycsamt.forward.maxwell import (
   ...     register_mt2d_backend, list_backends,
   ... )
   >>> from pycsamt.forward.maxwell.tri_fem2d import register_trifem2d_backend
   >>> from pycsamt.forward.maxwell.modem3d import register_modem3d_backend
   >>> from pycsamt.forward.maxwell.mare2dem import register_mare2dem_backend
   >>> register_mt2d_backend(replace=True)
   >>> register_trifem2d_backend(replace=True)
   >>> register_modem3d_backend(replace=True)
   >>> register_mare2dem_backend(replace=True)
   >>> for name, info in sorted(list_backends().items()):
   ...     cap = info["capabilities"]
   ...     print(name, "dims=", cap["dimensions"],
   ...           "verified=", cap["verified_benchmarks"], "available=", info["available"])
   mare2dem dims= [2] verified= [] available= False
   modem3d dims= [3] verified= ['half-space', 'layered-earth'] available= True
   mt2d dims= [2] verified= ['half-space', 'layered-earth'] available= True
   trifem2d dims= [2] verified= ['half-space', 'layered-earth'] available= True

``MT3DAdapter`` never appears here, on any machine: it has no
``register_mt3d_backend`` function, a deliberate omission
:doc:`maxwell_adapters` already explains -- research-only and cell-capped,
it is meant to be constructed directly rather than discovered by a
production dataset generator scanning the registry. ``modem3d`` reads
``available=True`` on the machine that built this page because a compiled
``Mod3DMT`` sits in pycsamt's vendored source tree, exactly as
:doc:`maxwell_adapters` found through :func:`~pycsamt.forward.maxwell.external.resolve_executable`
directly; ``mare2dem`` reads ``False`` for the same ``mpirun`` reason shown
there. ``register_mt2d_backend()`` a second time without ``replace=True``
raises rather than silently swapping the factory, since registration is
process-wide and a silent solver swap could change numerical behaviour
somewhere else in a long-running session:

.. code-block:: pycon

   >>> register_mt2d_backend()
   Traceback (most recent call last):
   ...
   ValueError: backend 'mt2d' is already registered.

:meth:`~pycsamt.forward.maxwell.BackendRegistry.names` can filter to only
what actually resolves right now, which is the cheap pre-check a script
should run before attempting to build every registered backend in a loop:

.. code-block:: pycon

   >>> from pycsamt.forward.maxwell import backend_registry
   >>> sorted(backend_registry.names())
   ['mare2dem', 'modem3d', 'mt2d', 'trifem2d']
   >>> backend_registry.names(available_only=True)
   ('modem3d', 'mt2d', 'trifem2d')

Private Registries
------------------

``register_backend``/``create_backend``/``list_backends`` all operate on one
process-wide :term:`backend registry`, which is exactly right for adapter
packages announcing themselves at import time. A test or an isolated
application that wants its own throwaway backend, without touching or being
affected by whatever else has already registered into the global registry
in the same process, can build a private :class:`~pycsamt.forward.maxwell.BackendRegistry`
instead:

.. code-block:: pycon

   >>> private_registry = BackendRegistry()
   >>> private_registry.register(
   ...     BackendRegistration(
   ...         BackendCapabilities("demo", "1", (2,), ("zxy",)), lambda **kw: None,
   ...     )
   ... )
   >>> private_registry.names()
   ('demo',)
   >>> "demo" in list_backends()
   False

``"demo"`` exists only inside ``private_registry``; the process-wide
registry queried through ``list_backends()`` never saw it. This is the
pattern worth reaching for in a test suite that registers a fake or
minimal backend to exercise adapter-selection logic -- a private registry
cannot leak a stray registration into a later, unrelated test the way a
``replace=True`` call against the global one eventually would.

Common Mistakes
---------------

``Reading assess() as a solver dry run``
   ``BackendCapabilities.assess`` and every adapter's own ``assess`` read
   shapes and declared flags only; they never touch conductivity values,
   receiver geometry beyond coordinate bounds, or numerical stability. A
   problem passing ``assess`` is compatible in shape, not guaranteed to
   solve cleanly.

``Treating a warning as a pass``
   ``report.compatible`` can be ``True`` while ``report.warnings`` is
   non-empty -- the "backend declares no verified benchmarks" warning fires
   on every solve through an unverified capability without ever blocking
   it. Read warnings explicitly rather than trusting ``compatible`` alone.

``Assuming supports_inactive_cells alone permits topography``
   A capability needs both ``supports_inactive_cells=True`` and
   ``supports_topography=True`` to accept a laterally varying air mask; the
   first flag alone only covers a flat, uniform air-layer count.

``Comparing capabilities objects instead of their name``
   ``BackendRegistration.create`` compares the *whole* returned
   ``capabilities`` object, not just its name -- a factory that bumps a
   version string or flips one capability flag after registration fails
   loudly rather than serving a mismatched instance.

``Polluting the global registry from a test``
   Calling ``register_mt2d_backend(replace=True)`` from inside a test
   mutates process-wide state that other tests or an interactive session
   may depend on. Prefer a private :class:`BackendRegistry` unless a test
   is specifically exercising the global registration functions themselves.

Next Pages
----------

:doc:`maxwell_adapters` is the natural companion to this page: it shows
what each concrete adapter's own ``assess`` adds on top of the generic
checks covered here. :doc:`maxwell_benchmarks` covers how the
``verified_benchmarks`` values read from the registry above are actually
produced and scored.
