.. _forward_maxwell_contracts:

Maxwell Problem, Mesh, And Result Contracts
===========================================

Every page in this section has already used
:class:`~pycsamt.forward.maxwell.MaxwellMesh`,
:class:`~pycsamt.forward.maxwell.MaxwellProblem`, and
:class:`~pycsamt.forward.maxwell.ForwardResult` without stopping to explain
what they actually enforce. :doc:`/theory/maxwell_forward` derives why that
enforcement has to exist, and
:doc:`/user_guide/ai_inversion/forward_physics` shows the same contracts
from a training-pipeline's point of view. This page stays with the objects
themselves: what each one stores, why its arrays cannot be edited in place,
the two different axis orders the package uses for two different purposes,
and how a problem or result survives being written to disk without pickle.
:class:`~pycsamt.forward.maxwell.contracts_tri.TriMesh` and
:class:`~pycsamt.forward.maxwell.contracts_tri.TriProblem`, the unstructured
sibling pair :doc:`maxwell_adapters` already used for ``TriFEM2DAdapter``
and ``Mare2DEMAdapter``, close the page.

Two Canonical Axis Orders
-------------------------

The package uses two different canonical orders, and conflating them is the
single easiest mistake to make against this contract.
:class:`~pycsamt.forward.maxwell.MaxwellMesh` cell arrays -- conductivity,
active cells -- follow the mesh's own ``(z, x)`` or ``(z, y, x)`` order,
depth first. :class:`~pycsamt.forward.maxwell.ForwardResult` impedance
follows a completely different order, ``(station, frequency, component)``,
because the output is organized around what a user actually indexes into: a
receiver's whole spectrum, not a depth slice.

.. code-block:: pycon

   >>> import numpy as np
   >>> from pycsamt.forward.maxwell import MaxwellMesh, ReceiverSet, MaxwellProblem

   >>> mesh = MaxwellMesh(np.linspace(0, 10_000, 41), np.linspace(0, 5_000, 31))
   >>> mesh.shape
   (30, 40)
   >>> conductivity = np.full(mesh.shape, 0.01)
   >>> receivers = ReceiverSet([[5_000.0, 0.0], [7_000.0, 0.0]], ["S00", "S01"])
   >>> problem = MaxwellProblem(mesh, conductivity, [10.0, 1.0], receivers, ("zxy", "zyx"))
   >>> problem.conductivity_s_m.shape
   (30, 40)

``mesh.shape`` is ``(nz, nx)``, and the conductivity array supplied to
``MaxwellProblem`` must match it exactly -- not transposed, even though a
transposed array is often still a valid-looking ``(nx, nz)`` shape for a
square-ish mesh and would otherwise pass silently into the wrong physical
orientation. A genuinely mismatched shape is rejected immediately:

.. code-block:: pycon

   >>> MaxwellProblem(mesh, conductivity.T, [10.0, 1.0], receivers, ("zxy", "zyx"))
   Traceback (most recent call last):
   ...
   ValueError: conductivity_s_m must be positive, finite, and shaped (30, 40).

Every array the contracts hand back is also read-only, on purpose: a
:term:`problem hash` describes exactly one conductivity array, and an
in-place edit after construction would silently invalidate that identity
without changing the hash a cache or benchmark had already recorded against
it.

.. code-block:: pycon

   >>> problem.conductivity_s_m[0, 0] = 5.0
   Traceback (most recent call last):
   ...
   ValueError: assignment destination is read-only

Changing a physical input means building a new ``MaxwellProblem``, which
produces a new, honest hash -- never mutating one in place.

Receivers And Immutability
--------------------------

:class:`~pycsamt.forward.maxwell.ReceiverSet` is a contract in its own
right, not just a coordinate array, because station names and orientation
are part of a result's identity: two receiver sets with identical
coordinates but different name order produce different
:attr:`~pycsamt.forward.maxwell.ForwardResult.receiver_names` axes.
``orientation_deg`` normalizes into ``[0, 360)`` on construction rather than
preserving whatever raw value was passed:

.. code-block:: pycon

   >>> rotated = ReceiverSet([[0.0, 0.0]], ["S00"], orientation_deg=-30)
   >>> rotated.orientation_deg
   330.0
   >>> rotated.dimension, rotated.count
   (2, 1)

Like every other contract in this module, it round-trips through a
versioned, JSON-compatible dictionary rather than pickle:

.. code-block:: pycon

   >>> state = rotated.to_dict()
   >>> state["orientation_deg"]
   330.0
   >>> ReceiverSet.from_dict(state).orientation_deg == rotated.orientation_deg
   True

Problem Identity And Provenance
-------------------------------

:attr:`~pycsamt.forward.maxwell.MaxwellProblem.problem_hash` is
deterministic in physical content, not object identity -- building the
exact same mesh, conductivity, frequencies, and receivers a second time
reproduces the identical digest, which is what lets
:doc:`maxwell_caching_and_batch`'s cache and :doc:`maxwell_benchmarks`'s
provenance record key on it at all:

.. code-block:: pycon

   >>> same_again = MaxwellProblem(mesh, conductivity, [10.0, 1.0], receivers, ("zxy", "zyx"))
   >>> same_again.problem_hash == problem.problem_hash
   True

The hash also covers ``metadata``, deliberately -- two otherwise-identical
problems that differ only in a metadata field are, by this contract's
definition, different problems:

.. code-block:: pycon

   >>> tagged_a = MaxwellProblem(
   ...     mesh, conductivity, [10.0, 1.0], receivers, ("zxy", "zyx"), metadata={"run": "a"},
   ... )
   >>> tagged_b = MaxwellProblem(
   ...     mesh, conductivity, [10.0, 1.0], receivers, ("zxy", "zyx"), metadata={"run": "b"},
   ... )
   >>> tagged_a.problem_hash != tagged_b.problem_hash
   True

That is exactly why :doc:`/theory/maxwell_forward` warns against putting a
generation timestamp in ``metadata``: a field meant only as a label would
silently turn every fresh run of otherwise-identical physics into a
guaranteed cache miss. :meth:`~pycsamt.forward.maxwell.MaxwellProblem.provenance`
returns the same identity information without the large numeric arrays,
which is what a stored record should hold instead of re-deriving the hash
from a code diff later:

.. code-block:: pycon

   >>> sorted(problem.provenance().keys())
   ['components', 'magnetic_permeability_h_m', 'mesh', 'metadata', 'receivers', 'schema_version', 'time_dependence']

Diagnostics And Result Validity
-------------------------------

:class:`~pycsamt.forward.maxwell.SolverDiagnostics` records convergence,
iteration count, and relative residual per frequency and per source
polarization -- one row per frequency, columns for however many independent
solves that frequency needed:

.. code-block:: pycon

   >>> from pycsamt.forward.maxwell import SolverDiagnostics
   >>> diagnostics = SolverDiagnostics([[True, True]], [[0, 0]], [[1e-9, 2e-9]], 0.05)
   >>> diagnostics.success, diagnostics.maximum_relative_residual
   (True, 2e-09)

``ForwardResult.success`` is not the same question as
``diagnostics.success``, and the difference matters: the first also
requires every entry of the :term:`validity mask` to be true, so a result
can have a fully converged linear solve and still read as unsuccessful
because a backend explicitly marked one output entry unusable for a reason
diagnostics never sees -- an extrapolated station, a degenerate geometry, a
known solver blind spot:

.. code-block:: pycon

   >>> from pycsamt.forward.maxwell import ForwardResult
   >>> single_freq_problem = MaxwellProblem(
   ...     mesh, conductivity, [10.0], ReceiverSet([[5_000.0, 0.0]], ["S00"]), ("zxy", "zyx"),
   ... )
   >>> impedance = np.array([[[1 + 1j, 1 + 1j]]])
   >>> partly_invalid = np.array([[[True, False]]])
   >>> converged = SolverDiagnostics([[True]], [[0]], [[1e-9]], 0.01)
   >>> result = ForwardResult(
   ...     single_freq_problem.problem_hash, [10.0], ["S00"], ("zxy", "zyx"),
   ...     impedance, partly_invalid, "demo", "1", converged,
   ... )
   >>> result.diagnostics.success
   True
   >>> result.success
   False

:meth:`~pycsamt.forward.maxwell.ForwardResult.validate_against` is the
check every adapter and cache runs before trusting a result belongs to a
given problem -- hash, frequency order, receiver order, and component set
must all match exactly:

.. code-block:: pycon

   >>> result.validate_against(single_freq_problem)
   >>> other_problem = MaxwellProblem(
   ...     mesh, conductivity, [10.0], ReceiverSet([[5_000.0, 0.0]], ["S00"]),
   ...     ("zxy", "zyx"), metadata={"run": "other"},
   ... )
   >>> result.validate_against(other_problem)
   Traceback (most recent call last):
   ...
   ValueError: result problem_hash does not match the problem.

Persistence Without Pickle
--------------------------

Both ``MaxwellProblem`` and ``ForwardResult`` write themselves to a
compressed, versioned ``.npz`` archive plus a JSON provenance block --
never pickle, which would tie the archive to a specific Python and package
version. This is exactly the mechanism
:class:`~pycsamt.forward.maxwell.MaxwellResultCache` wraps in
:doc:`maxwell_caching_and_batch`; here it runs directly, with nothing else
involved:

.. code-block:: pycon

   >>> from pathlib import Path
   >>> from tempfile import TemporaryDirectory
   >>> directory = TemporaryDirectory()
   >>> problem_path = Path(directory.name) / "problem.npz"
   >>> _ = problem.to_npz(problem_path)
   >>> restored_problem = MaxwellProblem.from_npz(problem_path)
   >>> restored_problem.problem_hash == problem.problem_hash
   True

   >>> result_path = Path(directory.name) / "result.npz"
   >>> _ = result.to_npz(result_path)
   >>> restored_result = ForwardResult.from_npz(result_path)
   >>> np.array_equal(restored_result.impedance_v_a, result.impedance_v_a)
   True
   >>> restored_result.backend_name
   'demo'

A round trip restores the exact hash and array content, not merely
something that looks similar -- ``from_npz`` reconstructs each contract
through its own validated constructor, so a corrupted or hand-edited
archive fails the same checks a freshly built object would.

The Triangular Mesh Contract
----------------------------

:class:`~pycsamt.forward.maxwell.contracts_tri.TriMesh` and
:class:`~pycsamt.forward.maxwell.contracts_tri.TriProblem` describe an
unstructured triangulation instead of a rectilinear grid, but they are
built from the same design: frozen, validated, JSON- and ``.npz``-
serializable, and identified by the same kind of :term:`problem hash`.

.. list-table::
   :header-rows: 1
   :widths: 22 32 32

   * -
     - ``MaxwellMesh`` / ``MaxwellProblem``
     - ``TriMesh`` / ``TriProblem``
   * - Geometry
     - Strictly increasing cell edges per axis
     - Node coordinates plus triangle connectivity
   * - Per-cell array shape
     - ``mesh.shape`` -- ``(nz, nx)`` or ``(nz, ny, nx)``
     - ``mesh.shape`` -- ``(n_triangles,)``
   * - Dimension
     - 2 or 3
     - Always 2
   * - Components
     - ``zxy``/``zyx`` in 2-D, all four in 3-D
     - ``zxy``/``zyx`` only

A triangle's own centroid and area are derived properties, not stored
redundantly, and every triangle is guaranteed to have strictly positive
area at construction:

.. code-block:: pycon

   >>> from pycsamt.forward.maxwell.contracts_tri import TriMesh, TriProblem
   >>> nodes = [[0, 0], [1_000, 0], [1_000, 500], [0, 500], [500, 250]]
   >>> triangles = [[0, 1, 4], [1, 2, 4], [2, 3, 4], [3, 0, 4]]
   >>> tri_mesh = TriMesh(nodes, triangles, boundary_segments=[[0, 1], [1, 2], [2, 3], [3, 0]])
   >>> tri_mesh.dimension, tri_mesh.n_nodes, tri_mesh.n_triangles, tri_mesh.shape
   (2, 5, 4, (4,))
   >>> tri_mesh.triangle_areas_m2.tolist()
   [125000.0, 125000.0, 125000.0, 125000.0]
   >>> tri_mesh.region_ids.tolist()
   [0, 0, 0, 0]

   >>> TriMesh([[0, 0], [1, 0], [2, 0]], [[0, 1, 2]])
   Traceback (most recent call last):
   ...
   ValueError: triangles must have strictly positive area (degenerate or collinear triangle detected).

``TriProblem``'s conductivity array is one value per triangle, matching
``mesh.shape`` exactly the same way ``MaxwellProblem``'s matches its
rectilinear mesh, and the same 2-D-only component restriction applies:

.. code-block:: pycon

   >>> tri_receivers = ReceiverSet([[500.0, 0.0]], ["S00"])
   >>> tri_problem = TriProblem(tri_mesh, np.full(tri_mesh.shape, 0.02), [10.0, 1.0], tri_receivers)
   >>> tri_problem.conductivity_s_m.shape
   (4,)
   >>> TriProblem(tri_mesh, np.full(tri_mesh.shape, 0.02), [10.0], tri_receivers, components=("zxx",))
   Traceback (most recent call last):
   ...
   ValueError: 2-D triangular problems support only zxy and zyx impedance components.

Common Mistakes
---------------

``Assuming a transposed conductivity array will just work``
   A conductivity array shaped ``(nx, nz)`` instead of ``(nz, nx)`` is only
   caught when the mesh is not square -- on a square mesh it would pass
   ``MaxwellProblem``'s shape check while describing a rotated model. Build
   conductivity arrays directly from ``mesh.shape``, never assume an axis
   order from a plotting convention.

``Editing metadata after the fact``
   Because ``metadata`` is part of the hash, "just adding a label" to an
   already-built problem means constructing a new ``MaxwellProblem`` -- there
   is no in-place update, by design.

``Treating diagnostics.success as the whole answer``
   A converged solve is necessary, not sufficient. Always check
   ``result.success`` (or the :term:`validity mask` directly) rather than
   ``result.diagnostics.success`` alone when deciding whether a result is
   trustworthy.

``Comparing TriProblem conductivity shape to a node count``
   ``TriProblem.conductivity_s_m`` has one entry per *triangle*
   (``mesh.n_triangles``), not per node. Building it from
   ``mesh.n_nodes`` is a shape mismatch that fails construction rather than
   producing a silently wrong per-node interpretation.

Next Pages
----------

:doc:`maxwell_meshing` covers how a real ``MaxwellMesh`` or ``TriMesh`` --
padded, air-layered, topography-aware, or graded around receivers -- gets
built in the first place, rather than the small ``np.linspace`` grids used
to demonstrate the contract on this page.
