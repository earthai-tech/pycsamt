.. _theory_maxwell_forward:

Maxwell Forward Modelling and Solver Contracts
==============================================

Electromagnetic inversion is meaningful only through the forward problem that
connects an earth model to observable fields.  The
:mod:`pycsamt.forward.maxwell` package is the boundary where pyCSAMT makes that
connection explicit.  It does not present every numerical solver as equivalent.
It defines immutable mesh, problem, receiver, result, diagnostic, capability,
benchmark, cache, and batch contracts; a numerical backend enters this namespace
only through an adapter that declares and enforces what it can actually solve.

That distinction matters in routine work.  A conductivity array is not yet a
Maxwell problem.  It needs a coordinate system, cell geometry, air and earth
regions, frequencies, receiver locations, tensor components, time convention,
boundary treatment, and an identified solver.  Likewise, a complex response
array is not yet scientific evidence.  Its axes must match the problem exactly,
and its convergence and validity diagnostics must pass a declared policy.

This page develops the physics and numerical reasoning behind those contracts,
then shows how to construct and assess them with the public API.  The practical
AI workflow in :doc:`../user_guide/ai_inversion/forward_physics` builds on the
same foundation.

From Maxwell's equations to an MT response
------------------------------------------

For harmonic fields with the time dependence
:math:`\exp(+i\omega t)`, Maxwell's equations in an isotropic conductive earth
can be written

.. math::
   :label: eq-maxwell-faraday

   \nabla\times\mathbf E=-i\omega\mu\mathbf H,

.. math::
   :label: eq-maxwell-ampere

   \nabla\times\mathbf H=
   \mathbf J_s+(\sigma+i\omega\varepsilon)\mathbf E,

where :math:`\mathbf E` is electric field, :math:`\mathbf H` magnetic field,
:math:`\mathbf J_s` an imposed source current, :math:`\sigma` conductivity,
:math:`\varepsilon` permittivity, and :math:`\mu` permeability.  In the
magnetotelluric frequency range and a sufficiently conductive earth, displacement
current is commonly neglected when

.. math::
   :label: eq-maxwell-quasistatic-condition

   \frac{\omega\varepsilon}{\sigma}\ll1.

Eliminating :math:`\mathbf H` then gives a curl--curl equation,

.. math::
   :label: eq-maxwell-curl-curl

   \nabla\times\mu^{-1}\nabla\times\mathbf E
   +i\omega\sigma\mathbf E
   =-i\omega\mathbf J_s.

Signs change consistently under :math:`\exp(-i\omega t)`.  This is why
``MaxwellProblem.time_dependence`` is a physical input and why a backend must
declare supported conventions.  Mixing conventions conjugates phase behavior;
it cannot be repaired by relabeling a plot.

At a surface receiver, the horizontal impedance tensor relates electric and
magnetic fields,

.. math::
   :label: eq-maxwell-impedance-tensor

   \begin{bmatrix}E_x\\E_y\end{bmatrix}
   =
   \begin{bmatrix}Z_{xx}&Z_{xy}\\Z_{yx}&Z_{yy}\end{bmatrix}
   \begin{bmatrix}H_x\\H_y\end{bmatrix}.

The canonical component names in pyCSAMT are ``zxx``, ``zxy``, ``zyx``, and
``zyy``.  A 2-D earth invariant along one horizontal direction separates into
TE and TM modes and supports the off-diagonal components ``zxy`` and ``zyx``.
Asking a 2-D problem for diagonal components is therefore rejected by the
contract rather than filled with arbitrary zeros.

Apparent resistivity and phase are derived responses,

.. math::
   :label: eq-maxwell-rhoa-phase

   \rho_a=\frac{|Z|^2}{\mu_0\omega},
   \qquad
   \phi=\operatorname{atan2}(\Im Z,\Re Z).

The solver contract retains complex impedance because its real and imaginary
parts preserve the quantity directly returned by the field ratio.  Derived
features can always be calculated later with the same unit and sign convention.

Why diffusion controls the mesh scale
-------------------------------------

In a uniform conductor, equations :eq:`eq-maxwell-faraday` and
:eq:`eq-maxwell-ampere` lead to a complex propagation constant.  The amplitude
decay scale is the skin depth

.. math::
   :label: eq-maxwell-skin-depth

   \delta=\sqrt{\frac{2}{\omega\mu\sigma}}
   =\sqrt{\frac{\rho}{\pi\mu f}}.

For :math:`\mu=\mu_0`, this is approximately

.. math::
   :label: eq-maxwell-skin-depth-rule

   \delta\approx503
   \sqrt{\frac{\rho\,[\Omega\,\mathrm m]}{f\,[\mathrm{Hz}]}}\ \mathrm m.

.. code-block:: pycon

   >>> from pycsamt.forward.maxwell import skin_depth_m
   >>> round(float(skin_depth_m(100, 1)))
   5033
   >>> round(float(skin_depth_m(10, 1000)))
   50

The smallest skin depth usually comes from the most conductive material at the
highest frequency.  It sets a demanding local cell scale.  The largest skin
depth, often associated with low frequency and resistive material, motivates
the extent of the earth and padding.  These two requirements compete: a mesh
that is fine enough everywhere and wide enough everywhere can become enormous,
especially in 3-D.

Skin depth is only a planning scale.  Thin layers, sharp contacts, receiver
interpolation, topography, anisotropy, and discretization error can require
finer cells.  Conversely, padding far from receivers may grow because it mainly
moves an artificial boundary away from the sensitive region.

.. figure:: /images/theory/maxwell_mesh_anatomy.png
   :alt: Padded topographic Maxwell mesh and skin depth curves
   :width: 100%

The left panel separates the geological core with dashed lines from air,
horizontal padding, and bottom padding.  Padding conductivity is extended from
the nearest earth cells, while the terrain mask replaces cells above the local
surface with numerical air.  The right panel shows that frequency and
resistivity jointly control penetration.  Its horizontal line is only a
four-cell scale.  In the executed example the minimum skin depth contains
``0.38`` core cells, so the mesh is deliberately too coarse for the most
conductive, highest-frequency case.  The warning is actionable numerical
evidence, not an exception to hide.

The canonical mesh contract
---------------------------

:class:`~pycsamt.forward.maxwell.MaxwellMesh` stores strictly increasing cell
edges, not only an image shape.  Its canonical array order is ``(z, x)`` in 2-D
and ``(z, y, x)`` in 3-D.  Depth increases downward.  Horizontal coordinates
and the optional CRS remain attached to the mesh.

.. code-block:: pycon

   >>> from pycsamt.forward.maxwell import MaxwellMesh
   >>> mesh = MaxwellMesh(
   ...     x_edges_m=[0, 100, 250, 500],
   ...     z_edges_m=[0, 50, 150, 350],
   ...     crs="EPSG:32630",
   ... )
   >>> mesh.dimension, mesh.shape
   (2, (3, 3))
   >>> mesh.cell_widths_m["x"].tolist()
   [100.0, 150.0, 250.0]
   >>> mesh.cell_centres_m["z"].tolist()
   [25.0, 100.0, 250.0]

Edges make cell width, volume, adjacency ratio, receiver containment, and
interpolation unambiguous.  A tuple such as ``(30, 40)`` cannot reveal whether
cells are 10 m or 10 km wide, uniform or geometrically stretched.

The arrays returned by the contracts are read-only.  This prevents a problem
hash from referring to one model while an in-place mutation silently changes
the conductivity later.  To alter a physical input, create a new contract and
therefore a new hash.

Building a solver mesh from geology
-----------------------------------

The geological grid describes the region where the earth model is interpreted.
A solver mesh normally needs more: air layers, horizontal padding, bottom
padding, extended material values, and explicit region masks.
:func:`~pycsamt.forward.maxwell.build_solver_mesh` performs that transformation
without choosing a numerical backend.

.. code-block:: pycon

   >>> import numpy as np
   >>> from pycsamt.ai.geology import GeologyGrid
   >>> from pycsamt.forward.maxwell import MeshDesign, build_solver_mesh
   >>> grid = GeologyGrid.regular_2d(
   ...     nx=20, nz=12, dx_m=150, dz_m=100,
   ... )
   >>> rho = np.full(grid.shape, 300.0)
   >>> rho[:3] = 40.0
   >>> design = MeshDesign(
   ...     horizontal_padding_cells=(4, 4),
   ...     bottom_padding_cells=5,
   ...     air_layers=4,
   ...     padding_expansion=1.3,
   ...     air_expansion=1.2,
   ... )
   >>> solver_model = build_solver_mesh(
   ...     grid,
   ...     resistivity_ohm_m=rho,
   ...     frequencies_hz=[100, 10, 1],
   ...     design=design,
   ... )
   >>> solver_model.mesh.shape
   (21, 28)
   >>> solver_model.core_slices
   (slice(4, 16, None), slice(4, 24, None))
   >>> solver_model.quality.cell_count
   588

Exactly one of ``conductivity_s_m`` and ``resistivity_ohm_m`` must be supplied.
The returned model always stores conductivity in SI units.  This explicit
conversion avoids a particularly destructive error: passing resistivity values
to a solver that expects conductivity.

The core slices map the original geological model into the padded array.  They
are essential when a response-based optimization runs on the solver mesh but
the recovered model must be exported on the geological grid.  Do not reconstruct
them from padding counts after serialization; use the stored contract.

What ``MeshDesign`` controls
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``horizontal_padding_cells`` accepts a symmetric integer or a low/high pair.
``bottom_padding_cells`` extends the deepest model value.  ``air_layers`` adds
cells above depth zero.  ``padding_expansion`` and ``air_expansion`` are
geometric growth factors moving away from the core.

For core boundary cell width :math:`h_0` and growth factor :math:`g`, padding
width :math:`k` is

.. math::
   :label: eq-maxwell-padding-width

   h_k=h_0g^k,

and the total extent of :math:`n` cells is

.. math::
   :label: eq-maxwell-padding-extent

   L_{pad}=h_0\sum_{k=1}^{n}g^k
   =h_0g\frac{g^n-1}{g-1}.

Large growth moves the boundary away with few cells but creates abrupt width
changes.  Small growth is smoother but increases cell count.  The correct
choice is established by mesh convergence, not appearance.

``air_conductivity_s_m`` is a positive numerical conductivity, not geological
resistivity.  Zero conductivity can make discrete systems singular or violate
backend assumptions.  Whether air should also be marked inactive depends on
the backend capability.

Topography is a physical region mask
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

An aligned :class:`~pycsamt.ai.geology.TopographicSurface` defines the local
depth of earth.  The builder interpolates that surface onto padded horizontal
cell centres and constructs complementary ``earth_mask`` and ``air_mask``
arrays.  It does not merely warp a display.

.. code-block:: pycon

   >>> from pycsamt.ai.geology import TopographicSurface
   >>> elevation = 420 + 25 * np.sin(
   ...     2 * np.pi * grid.x_m / np.ptp(grid.x_m)
   ... )
   >>> surface = TopographicSurface(
   ...     grid, elevation, float(elevation.max()),
   ...     source="station elevations",
   ... )
   >>> topo_model = build_solver_mesh(
   ...     grid, resistivity_ohm_m=rho,
   ...     frequencies_hz=[100, 10, 1],
   ...     topography=surface, design=design,
   ... )
   >>> topo_model.earth_mask.shape == topo_model.mesh.shape
   True
   >>> np.array_equal(topo_model.air_mask, ~topo_model.earth_mask)
   True

Topography and the grid must share dimension, coordinates, and CRS.  Failing
early is preferable to aligning two arrays by shape when they represent
different locations.

Reading mesh-quality diagnostics
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:class:`~pycsamt.forward.maxwell.MeshQuality` records total cells, width
extremes, global aspect ratio, worst adjacent-width ratio, minimum skin depth,
cells per minimum skin depth, and advisory warnings.  Its ``acceptable``
property means only that configured advisory limits were not violated.  It is
not proof of forward accuracy.

The reported global width ratio is

.. math::
   :label: eq-maxwell-global-width-ratio

   R_h=\frac{\max_{a,k}h_{a,k}}{\min_{a,k}h_{a,k}},

while the adjacent ratio is

.. math::
   :label: eq-maxwell-adjacent-ratio

   R_{adj}=\max_{a,k}
   \left(\frac{h_{a,k+1}}{h_{a,k}},
         \frac{h_{a,k}}{h_{a,k+1}}\right).

The first describes overall scale separation; the second detects abrupt local
growth.  A large global ratio can be acceptable across many gradual padding
cells, whereas a large adjacent ratio can damage numerical accuracy locally.

If the skin-depth warning fires, possible corrections are to refine core cells,
remove unreliable high frequencies, limit unrealistically conductive prior
values, or use local refinement supported by another mesh builder/backend.  Do
not simply lower ``minimum_cells_per_skin_depth`` to make the warning disappear.

The solver-neutral Maxwell problem
----------------------------------

:class:`~pycsamt.forward.maxwell.MaxwellProblem` binds the complete physical
input.  Its required relationships are validated at construction:

* conductivity is positive, finite, and shaped like the mesh;
* receiver dimension equals mesh dimension;
* frequencies are positive, finite, and unique;
* component names are canonical and unique;
* 2-D components are restricted to ``zxy`` and ``zyx``;
* active cells match the mesh and include at least one cell;
* permeability is positive and the time convention is supported by the
  contract;
* metadata is finite and JSON-compatible.

Receivers are a separate contract because names, coordinates, and orientation
are part of output identity.

.. code-block:: pycon

   >>> from pycsamt.forward.maxwell import ReceiverSet
   >>> receivers = ReceiverSet(
   ...     [[450, 0], [1050, 0], [1650, 0], [2250, 0]],
   ...     ["S01", "S02", "S03", "S04"],
   ... )
   >>> topo_model.assess_receivers(receivers)
   ()
   >>> problem = topo_model.to_problem(
   ...     [100, 10, 1], receivers,
   ...     components=("zxy", "zyx"),
   ...     mark_air_inactive=False,
   ...     metadata={"survey": "teaching-line"},
   ... )
   >>> problem.mesh.dimension, problem.receivers.count
   (2, 4)
   >>> len(problem.problem_hash)
   64

``assess_receivers`` reports positions outside the mesh or below the local
discretized terrain without silently snapping them.  The example retains air
as active low-conductivity cells because the in-repository ``mt2d`` adapter does
not support inactive cells or topographic masks.  That also means a truly
topographic problem is not within that adapter's stated capability merely
because a padded conductivity array can be formed.

Problem hashes and provenance
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The deterministic SHA-256 ``problem_hash`` covers conductivity, frequencies,
active cells, mesh, receivers, components, convention, permeability, and
metadata.  It is the identity used to validate results and key caches.

.. math::
   :label: eq-maxwell-problem-identity

   h_p=\operatorname{SHA256}
   (\mathbf\sigma,\mathbf f,\mathbf A,
   \mathcal M,\mathcal R,\mathbf c,\tau,\mu,\mathcal P).

A changed station order produces a changed problem even if coordinates are the
same, because output axis order changes.  A changed metadata field also changes
the hash by design; metadata should therefore contain stable physical
provenance, not timestamps that differ on every identical run.

Backend capabilities are scientific claims
-------------------------------------------

:class:`~pycsamt.forward.maxwell.BackendCapabilities` declares dimension,
components, time conventions, nonuniform-mesh support, inactive-cell support,
topography, anisotropy, cell/frequency limits, and verified benchmarks.

.. code-block:: pycon

   >>> from pycsamt.forward.maxwell import BackendCapabilities
   >>> capability = BackendCapabilities(
   ...     name="teaching-backend", version="1.0",
   ...     dimensions=(2,), components=("zxy", "zyx"),
   ...     supports_nonuniform_mesh=True,
   ...     supports_inactive_cells=False,
   ...     supports_topography=False,
   ...     maximum_cells=10000,
   ...     verified_benchmarks=("half-space",),
   ... )
   >>> report = capability.assess(problem)
   >>> report.compatible
   True
   >>> report.errors
   ()

A capability declaration does not certify itself.  ``verified_benchmarks``
must correspond to stored outcomes for the exact adapter and solver version.
Changing translation logic or numerical dependencies can require the adapter
version and evidence to be updated together.

Compatibility assessment is preflight, not execution.  It gathers hard errors
such as wrong dimension, unsupported components, inactive cells, nonuniform
mesh, topography, or size limits before an expensive solve begins.  Warnings are
advisory validation concerns.  Users can inspect both; adapters call the same
assessment automatically.

The adapter lifecycle
---------------------

An adapter translates between the canonical contracts and one backend's native
input/output.  It is deliberately more than a thin wrapper.

.. figure:: /images/theory/maxwell_adapter_lifecycle.png
   :alt: Maxwell adapter preflight, backend execution, and postflight lifecycle
   :width: 100%

Preflight rejects problems outside declared capabilities.  The backend then
solves only compatible input.  Postflight verifies result type, problem hash,
frequency order, station order, components, backend identity, convergence,
residual threshold, and response validity.  A backend cannot silently return a
smaller frequency set or reordered stations and still pass.

:class:`~pycsamt.forward.maxwell.AdapterPolicy` controls result acceptance:

.. code-block:: pycon

   >>> from pycsamt.forward.maxwell import AdapterPolicy
   >>> policy = AdapterPolicy(
   ...     require_convergence=True,
   ...     maximum_relative_residual=1e-6,
   ...     require_all_valid=True,
   ... )
   >>> policy.maximum_relative_residual
   1e-06

The exceptions distinguish failure stages:

* ``IncompatibleProblemError`` means preflight rejected the physics or
  numerical scope;
* ``BackendExecutionError`` wraps an ordinary solver failure;
* ``InvalidBackendResultError`` means output violated the canonical contract;
* ``SolverConvergenceError`` means diagnostics or validity failed policy.

Treating all four as a generic empty result would destroy the evidence needed
to correct the workflow.

Adapting a trusted callable
~~~~~~~~~~~~~~~~~~~~~~~~~~~

:class:`~pycsamt.forward.maxwell.CallableMaxwellAdapter` is useful when a
trusted Python callable already accepts ``MaxwellProblem`` and returns a fully
formed ``ForwardResult``.  It does not convert arbitrary arrays automatically;
the callback remains responsible for canonical output and diagnostics.

.. code-block:: pycon

   >>> from pycsamt.forward.maxwell import (
   ...     BackendCapabilities, CallableMaxwellAdapter,
   ...     ForwardResult, SolverDiagnostics, half_space_impedance,
   ... )
   >>> demo_capability = BackendCapabilities(
   ...     "analytic-demo", "1", (2,), ("zxy", "zyx"),
   ... )
   >>> def analytic_solver(value):
   ...     z = np.zeros((value.receivers.count,
   ...                   len(value.frequencies_hz),
   ...                   len(value.components)), dtype=complex)
   ...     reference = half_space_impedance(
   ...         100.0, value.frequencies_hz,
   ...     )
   ...     z[:] = reference[None, :, None]
   ...     diagnostics = SolverDiagnostics(
   ...         np.ones((len(value.frequencies_hz), 1), bool),
   ...         np.zeros((len(value.frequencies_hz), 1), int),
   ...         np.zeros((len(value.frequencies_hz), 1)), 0.0,
   ...     )
   ...     return ForwardResult(
   ...         value.problem_hash, value.frequencies_hz,
   ...         value.receivers.names, value.components,
   ...         z, None, "analytic-demo", "1", diagnostics,
   ...     )
   >>> adapter = CallableMaxwellAdapter(
   ...     demo_capability, analytic_solver, policy,
   ... )
   >>> result = adapter.solve(problem)
   >>> result.shape, result.success
   ((4, 3, 2), True)

This example is an analytic half-space response used to demonstrate the
contract.  It is not a multidimensional solution of the heterogeneous
``problem`` conductivity.  Naming it ``analytic-demo`` prevents that distinction
from being hidden behind a plausible response array.

In-repository and external adapters
-----------------------------------

``MT2DAdapter`` translates compatible problems to the in-repository 2-D
finite-difference implementation.  Its declared scope is precise:

* two-dimensional ``zxy`` and ``zyx`` responses;
* surface receivers at ``z=0``;
* :math:`\exp(+i\omega t)`;
* vacuum permeability;
* nonuniform meshes;
* no inactive cells, topography, or anisotropy.

It declares passing half-space and layered-earth benchmarks.  These statements
are stronger and more useful than saying merely that “2-D is supported.”

``MT3DAdapter`` is explicitly research-only.  It supports uniform 3-D meshes
and all four impedance components but defaults to a 6,000-cell safety ceiling,
has no inactive cells or topography, and is not a production-scale 3-D solver.
Increasing the ceiling changes computational risk, not numerical maturity.

``ModEm3DAdapter`` provides an external-process integration.  External adapters
must additionally resolve an executable, isolate a run directory, write native
files, execute under a timeout policy, parse output, convert units, and retain
stdout/stderr provenance.  Availability of an executable is not the same as
compatibility of a problem or validation of its result.

The external run policy controls timeout, environment, retained workspaces, and
process behavior.  A reproducible record should preserve executable path and
version, command arguments, native input hashes, return code, and parsed-output
identity.  Never infer success solely because a results file exists; it may be
stale or incomplete.

Lazy backend registration
-------------------------

The registry stores factories rather than importing every optional solver when
``pycsamt.forward.maxwell`` is imported.  This keeps the core contracts usable
without external executables or optional solver libraries.

.. code-block:: pycon

   >>> from pycsamt.forward.maxwell import (
   ...     create_backend, list_backends, register_mt2d_backend,
   ... )
   >>> register_mt2d_backend(replace=True)
   >>> "mt2d" in list_backends()
   True
   >>> mt2d = create_backend("mt2d", verbose=False)
   >>> mt2d.capabilities.dimensions
   (2,)
   >>> mt2d.capabilities.verified_benchmarks
   ('half-space', 'layered-earth')

Registration is process-wide.  ``replace=True`` should be explicit because
silently replacing a factory could change numerical behavior elsewhere in a
long-running experiment.  List capability metadata before creating or running
a backend when building user-facing selection interfaces.

Canonical results and diagnostics
---------------------------------

:class:`~pycsamt.forward.maxwell.ForwardResult` stores impedance with shape
``(station, frequency, component)``.  It includes all three axis labels, a
validity mask, backend identity, problem hash, diagnostics, and metadata.

:class:`~pycsamt.forward.maxwell.SolverDiagnostics` records convergence,
iterations, and relative residual for each frequency/source solve plus total
runtime.  For a linear system

.. math::
   :label: eq-maxwell-linear-system

   \mathbf A\mathbf u=\mathbf b,

a common reported relative residual is

.. math::
   :label: eq-maxwell-relative-residual

   r_{rel}=\frac{\|\mathbf A\widehat{\mathbf u}-\mathbf b\|_2}
                  {\|\mathbf b\|_2}.

A small algebraic residual says the discrete system was solved accurately.  It
does not prove that the mesh, boundary condition, source representation, or
earth model accurately represents the continuous problem.  Solver convergence
and mesh convergence answer different questions.

.. code-block:: pycon

   >>> result.backend_name, result.backend_version
   ('analytic-demo', '1')
   >>> result.receiver_names
   ('S01', 'S02', 'S03', 'S04')
   >>> result.components
   ('zxy', 'zyx')
   >>> result.diagnostics.success
   True
   >>> result.validate_against(problem) is None
   True

The validity mask separates missing or unusable responses from legitimate
complex zeros.  Downstream misfit code should combine it with observed-data
validity rather than replacing invalid entries by zero.

Analytic benchmarks before geological complexity
------------------------------------------------

A half-space has a closed-form impedance under the package convention,

.. math::
   :label: eq-maxwell-half-space-impedance

   Z(\omega)=\sqrt{i\omega\mu\rho},

with the square-root branch consistent with the chosen time dependence.  Its
apparent resistivity is constant and its phase is 45 degrees in the ideal
quasi-static case.  A layered-earth recurrence provides a stronger frequency-
dependent reference.

.. code-block:: pycon

   >>> from pycsamt.forward.maxwell import (
   ...     half_space_impedance, layered_earth_impedance,
   ... )
   >>> frequencies = np.array([100.0, 10.0, 1.0])
   >>> z_half = half_space_impedance(100.0, frequencies)
   >>> np.round(np.rad2deg(np.angle(z_half)), 1).tolist()
   [45.0, 45.0, 45.0]
   >>> z_layered = layered_earth_impedance(
   ...     [30.0, 600.0, 10.0], [250.0, 700.0], frequencies,
   ... )
   >>> z_layered.shape
   (3,)

.. figure:: /images/theory/maxwell_analytic_benchmarks.png
   :alt: Half-space and layered-earth apparent resistivity and phase responses
   :width: 94%
   :align: center

The half-space curve remains flat in apparent resistivity and constant in phase.
The layered model changes with frequency because the relative influence of its
conductive cover, resistive middle layer, and conductive basement changes with
penetration.  A solver that cannot reproduce these responses within declared
amplitude and phase tolerances should not be trusted on a more complicated
geology merely because its image looks smooth.

Executable benchmark objects
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``half_space_benchmark`` and ``layered_earth_benchmark`` package a problem,
analytic reference, and identity.  ``run_benchmarks`` evaluates an adapter and
returns metrics against :class:`~pycsamt.forward.maxwell.BenchmarkThresholds`.
The default criteria include normalized complex RMS, maximum relative amplitude
error, maximum phase error, valid fraction, and convergence.

For prediction :math:`Z_j` and reference :math:`Z_j^*`, normalized RMS is

.. math::
   :label: eq-maxwell-benchmark-nrms

   \operatorname{NRMS}=
   \sqrt{\frac{\sum_j|Z_j-Z_j^*|^2}
                {\sum_j|Z_j^*|^2}}.

Circular phase error must wrap differences into a principal interval so values
near :math:`-180^\circ` and :math:`180^\circ` are recognized as neighbors.
Amplitude, phase, validity, and convergence gates are retained separately;
passing one cannot compensate for failing another.

Mesh-convergence evidence
-------------------------

Analytic agreement on canonical cases is necessary but not sufficient.  For a
representative model, solve a sequence of meshes with characteristic core sizes
:math:`h`, :math:`h/r`, and :math:`h/r^2`.  Compare responses after exact axis
alignment.  A relative change measure is

.. math::
   :label: eq-maxwell-mesh-convergence

   C_h=\frac{\|\mathbf Z_h-\mathbf Z_{h/r}\|_2}
             {\|\mathbf Z_{h/r}\|_2}.

The result is credible only when changes decline toward the application
tolerance and the finest mesh is itself within capability limits.  Perform
similar tests for padding extent and boundary placement.  A single fine-looking
mesh supplies no convergence trend.

Topography convergence should refine both horizontal terrain sampling and
vertical cells near the surface.  If moving receivers between terrain and a
flat reference changes the physical question, do not mix that geometry change
with mesh refinement in the same comparison.

Caching without confusing identity
----------------------------------

:class:`~pycsamt.forward.maxwell.MaxwellResultCache` stores results under stable
keys with file digests and locking.  A cache hit is valid only when the exact
problem and expected backend identity match.  Cache statistics distinguish hits,
misses, writes, and corrupt entries.

Caching is especially valuable for AI dataset generation and repeated loss
evaluation because Maxwell solves dominate cost.  It must not conceal a solver
upgrade: include backend version in the cache expectation or use a new cache
namespace.  Corruption errors should trigger quarantine or recomputation, never
silent acceptance of a partially written archive.

Concurrent writers require locking.  A lock timeout means another process may
be producing the same key or a stale lock needs investigation.  Deleting broad
cache directories automatically is not an appropriate response; preserve the
failure and resolve the exact entry.

Batch solving and failure manifests
-----------------------------------

``solve_batch`` applies one adapter to multiple problems with a
:class:`~pycsamt.forward.maxwell.BatchPolicy`.  Batch behavior must distinguish
fail-fast scientific experiments from large synthetic campaigns where a small
number of documented failures can be retained and retried.

A :class:`~pycsamt.forward.maxwell.FailureManifest` stores the problem hash,
failure type, message, and position.  This is essential because dropping failed
realizations without recording why can bias a training distribution.  For
example, if high-contrast conductive bodies fail more often, silently removing
them narrows the prior and makes the eventual inverse model weakest on the very
targets of interest.

Batch parallelism should respect backend thread safety, external executable
licenses, memory per solve, and cache locking.  More workers can reduce
throughput when sparse factorizations compete for memory.  Benchmark sustained
throughput and failure rate, not only the fastest individual solve.

Choosing a backend responsibly
------------------------------

Backend selection begins with the physical question:

#. Is the earth approximation 1-D, 2-D, or genuinely 3-D?
#. Which impedance components and time convention are required?
#. Are cells nonuniform, inactive, anisotropic, or topographic?
#. Are receivers on the surface or buried?
#. How many cells and frequencies does the problem require?
#. Which canonical benchmarks has the exact version passed?
#. What mesh-convergence evidence exists at comparable contrasts and geometry?

Only after those questions should runtime and convenience decide among
compatible backends.  A fast solver outside its declared dimension or terrain
support is not an approximation with known error; it is a different problem.

Common failure modes
--------------------

**The mesh builds but quality is unacceptable.**  The builder validates data
structure and reports advisory numerical risks.  Refine the core, adjust the
frequency band, examine extreme conductivity, or redesign padding.  Do not
equate successful construction with adequate discretization.

**A topographic model is rejected by ``mt2d``.**  The mesh builder can represent
topographic air/earth regions, but the adapter does not claim inactive-cell or
topography support.  Use a validated capable backend or explicitly solve a flat-
surface approximation and state that limitation.

**Receivers are below local terrain.**  Check the vertical datum, the meaning of
depth zero, station/topography alignment, and coordinate order.  Automatic
snapping would conceal these errors, so ``assess_receivers`` reports them.

**The solver converges but benchmarks fail.**  Algebraic convergence applies to
the implemented discrete system.  Investigate field normalization, boundary
conditions, sign convention, receiver interpolation, units, and mesh extent.

**The result has the right shape but validation fails.**  A reordered frequency,
station, or component axis is scientifically different despite equal shape.
Return the axes exactly as supplied by ``MaxwellProblem``.

**A 3-D solve exceeds the cell ceiling.**  The research adapter's direct sparse
method is not production-scalable.  Coarsening until it runs may destroy the
physics.  Prefer a validated external production backend and establish mesh
convergence within available resources.

**Responses look identical across heterogeneous models.**  Confirm conductivity
was not passed as resistivity, the model was mapped into ``core_slices``, the
backend actually consumes the supplied array, and cached results use different
problem hashes.

**Phase has the opposite sign.**  Check time dependence, impedance definition,
component orientation, and unit conversion before applying any manual sign
change.  The convention must be consistent through analytic reference, backend,
observations, and plots.

What reproducibility requires
-----------------------------

A defensible forward result should retain:

* geological model coordinates, units, CRS, and property parameterization;
* the complete ``MeshDesign``, solver mesh edges, core slices, region masks,
  model hash, and quality warnings;
* receiver names, coordinates, orientation, and topography/datum source;
* frequencies, components, time convention, permeability, and problem hash;
* backend and adapter name/version, executable version where applicable, and
  capability report;
* boundary conditions and native solver settings not represented by the common
  contract;
* convergence, iterations, residuals, validity, runtime, and backend messages;
* benchmark outcomes and mesh-convergence evidence for the declared use;
* cache identity or external run-directory provenance when used.

The ``provenance()``, ``to_dict()``, and archive methods on the contracts provide
the machine-readable foundation.  They do not replace the scientific narrative
explaining why the mesh and backend are suitable.

Relationship to inversion and AI datasets
------------------------------------------

In classical inversion, the forward operator is evaluated repeatedly as the
earth model changes.  Numerical error must remain below the data-error scale;
otherwise the optimizer can fit discretization artifacts.  Mesh and solver
identity should remain fixed during one objective comparison unless the change
is explicitly accounted for.

In supervised AI inversion, each label-response pair is

.. math::
   :label: eq-maxwell-ai-pair

   \mathbf m_i\sim p(\mathbf m),
   \qquad
   \mathbf d_i=\mathcal F_{h,b}(\mathbf m_i)+\boldsymbol\epsilon_i,

where :math:`h` denotes mesh/discretization and :math:`b` backend.  A network can
learn systematic error in :math:`\mathcal F_{h,b}` as readily as earth physics.
Dataset manifests must therefore retain problem hashes, backend versions,
diagnostics, and failures.  Cross-solver validation on a subset estimates how
strongly the learned distribution depends on one implementation.

For response-consistency losses, predictions must be mapped from the geological
grid to the same solver mesh and evaluated through a compatible adapter.  A
skin-depth proxy or analytic half-space response is useful for teaching and unit
tests but is not a substitute for heterogeneous Maxwell evidence.

Reproduce the figures
---------------------

.. code-dropdown:: ../../scripts/generate_maxwell_theory_figures.py
   :language: python
   :linenos:
   :title: View and copy the complete Maxwell theory figure generator

Executed output:

.. code-block:: text

   solver mesh cells: 1064
   cells per minimum skin depth: 0.38
   mesh quality warnings: 1
   layered impedance magnitude range: 0.000937496 1.53906
   figures generated: 3

Continue with :doc:`impedance_tensor` for tensor interpretation,
:doc:`dimensionality` for selecting an earth approximation,
:doc:`inversion_concepts` for the inverse objective, and
:doc:`ai_inversion` for the role of Maxwell simulations in learned inversion.
