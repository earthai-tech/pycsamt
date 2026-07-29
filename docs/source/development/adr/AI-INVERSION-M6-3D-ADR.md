# M6 — 3-D Maxwell Solver Feasibility and Architecture Decision Record

Status: decided (gate evaluated 2026-07-27)
Scope: `pycsamt.forward.maxwell` M6 milestone per `AI-INVERSION-IMPLEMENTATION-PLAN.md`
Rule: `mt3d.py` (an in-house production 3-D adapter) may not be implemented until this record's gate passes; per the plan, a failed gate means 3-D production training uses a trusted external backend and any in-house solver stays research-only.

## 1. Governing equation

The frequency-domain equation this milestone must eventually support, with one documented time convention:

`curl(mu^-1 curl E) + i*omega*sigma*E = source`

with complex fields, anisotropic conductivity when supported, divergence/null-space handling, and impedance/tipper extraction at receivers.

## 2. What was actually checked (not assumed)

- **Installed Python 3-D EM packages**: `SimPEG`, `simpeg`, `discretize`, `pymatsolver`, `pygimli`, `emg3d` were all probed by direct `importlib.import_module` in this environment on 2026-07-27 and are **not installed**. Only `mtpy` (1.1.5) is present, which is a data-processing/plotting package, not a forward solver.
- **ModEM**: `pycsamt/models/modem/_source/{2D,3D}/*.f90` bundles ModEM v6.2.6 Fortran source (`Mod3DMT.f90`, Copyright 2004-2014 Oregon State University; Egbert, Kelbert, Meqbel) with build instructions in `pycsamt/models/modem/_source/README.txt`. No compiled `Mod3DMT`/`Mod2DMT` binary exists in this environment or on `PATH` — it requires a separate `gfortran`+LAPACK/BLAS(+MPI) build step the user must perform. `pycsamt.models.modem.runner.ModEmRunner` already exists as an application-level subprocess wrapper (resolves the binary, builds the command line, launches it, loads results) but has not been wrapped as a `pycsamt.forward.maxwell.backends.MaxwellBackend`.
- **No explicit ModEM license file** was found alongside the vendored source (only the copyright header naming the OSU authors). ModEM's exact redistribution/production-use terms should be confirmed with the original authors before any commercial or redistributed production deployment; this is called out as an open item below, not asserted as resolved.
- **Existing "3-D" forward code**: `pycsamt/forward/em3d.py::MT3DForward` is explicitly documented as **quasi-3-D** — it averages two independent 2-D `MT2DForward` solves (XZ and YZ cross-sections) per station, sets `Zxx = Zyy = 0` unconditionally, and its own docstring states the approximation "breaks down for strongly 3-D structures where galvanic coupling between adjacent columns is significant." It is not a candidate for "genuine 3-D" per the plan's own definition and is out of scope for M6/M7.
- **Mesh infrastructure reuse**: `pycsamt/forward/maxwell/mesh.py` (`MaxwellMesh`, `MeshDesign`, `build_solver_mesh`, `SolverMeshModel`) is already dimension-agnostic — its 3-D branch (padding on `y`, 3-D earth/air masking via topography, `core_slices`) is implemented and covered by `pycsamt/forward/tests/test_maxwell_mesh.py`. Whichever 3-D backend is chosen below can reuse this mesh layer unchanged; it does not need to be rebuilt.

## 3. Feasibility benchmark: direct vs. iterative solves at 3-D scale

**Method.** A complex-valued, diagonally-dominant 7-point-stencil sparse operator (`6.0 + 0.5j` on the diagonal, `-1` on the six axis neighbors) was built at increasing grid sizes `n x n x n` as a *computational proxy* for the DOF count and sparsity pattern of a real edge-based 3-D EM assembly — this is **not** a physical Maxwell solve, and the matrix is deliberately well-conditioned (diagonally dominant), unlike a real indefinite Maxwell system. It is used only to measure how direct vs. iterative solve **cost scales with problem size** on this machine, which is a property of the linear-algebra, not the physics. Measured with `scipy 1.10.0`, `scipy.sparse.linalg.spsolve` (direct, SuperLU) and `bicgstab` (iterative, no preconditioner), wall-clock via `time.perf_counter`:

| n (cells/side) | N = n³ (DOF) | nnz | direct `spsolve` | iterative `bicgstab` (unpreconditioned) |
|---:|---:|---:|---:|---:|
| 10 | 1,000 | 6,400 | 0.012 s | 0.011 s (23 iters) |
| 16 | 4,096 | 27,136 | 0.100 s | 0.006 s (28 iters) |
| 22 | 10,648 | 71,632 | 1.052 s | 0.013 s (29 iters) |
| 28 | 21,952 | 148,960 | 4.628 s | 0.030 s (28 iters) |

Direct-solve time grew ~4.4x for a ~2.06x increase in DOF (n: 22->28), consistent with the well-known O(N^2) (or worse) fill-in scaling of 3-D sparse LU factorization (Davis, *Direct Methods for Sparse Linear Systems*, SIAM, 2006) — 3-D sparse matrices lack the good nested-dissection separators that make 2-D direct solves (used successfully by this project's own `mt2d.py`) practical. Note: peak-memory was also measured via Python `tracemalloc` but is **not reported** above because it only captures Python-level allocations, not the C-level SuperLU fill-in buffers where the real memory cost lives; the timing trend alone is sufficient to establish the scaling conclusion, and inflating it with a known-to-undercount memory number would overstate this benchmark's precision.

**Extrapolation.** A single small production-scale 3-D MT mesh (e.g. 40-60 core cells per axis plus 6-8 padding cells and 6-8 air layers per `MeshDesign`'s defaults) reaches at least n~60-80 per axis, i.e. N = 216,000-512,000 cells before accounting for the ~3x multiplier of a real edge-element DOF count (three field components per cell). Extrapolating the measured O(N^~2.1) direct-solve trend from n=28 (4.63 s) to n=80 (23.3x more DOF) gives an estimated **~35-45 minutes for a single direct solve of one frequency/polarization** — before multiplying by polarizations (2), frequencies (10-30 typical), and training realizations (hundreds to thousands). This is not tractable for dataset generation. The unpreconditioned iterative solve stayed fast in this proxy benchmark, but this **cannot** be read as evidence that iterative 3-D EM solves are easy in general: real Maxwell systems are complex-symmetric/indefinite and famously ill-conditioned near DC, which is exactly why the published 3-D EM literature (ModEM, SimPEG, emg3d) universally pairs Krylov solvers with problem-specific preconditioners (multigrid, shifted-Laplacian, or similar) rather than using them unpreconditioned. Designing, implementing, and validating such a preconditioner in-house is itself a substantial, specialized research effort.

## 4. Decisions

### 4.1 Preferred verified backend

**Decision: ModEM, via a to-be-built `BaseExternalMaxwellAdapter` subclass (a future `modem3d.py` in `pycsamt/forward/maxwell/`) — not an in-house Python 3-D solver.**

Rationale: no verified Python 3-D EM package is available in this environment (§2); a genuine in-house FD/FE 3-D solver would additionally need a problem-specific preconditioner, which is out of scope for near-term production per §3's cost estimate; ModEM source is already vendored with a working application-level runner (`pycsamt.models.modem.runner.ModEmRunner`) and is the specific example the plan itself names ("an existing pycsamt-compatible ModEM workflow"). The adapter foundation this ADR can build on (`pycsamt.forward.maxwell.external.BaseExternalMaxwellAdapter`, `ExternalRunPolicy`, `resolve_executable`, `make_availability_probe`) already exists and was purpose-built for exactly this kind of file-based external-executable integration (M4 PR4).

Open item: confirm ModEM's exact redistribution/production-use terms with the original OSU authors before any commercial or redistributed production deployment; no explicit license file accompanies the vendored source (§2).

### 4.2 Discretization, sparse solver, and preconditioner

**Decision: if an in-house research-only prototype is pursued later (never as the production path — see §4.1), it should use an edge-based finite-volume/finite-difference formulation on a Yee-type staggered grid** (the same family ModEM itself uses; Egbert & Kelbert, 2012), **solved with a preconditioned Krylov method (BiCGSTAB or QMR)**, never a direct sparse solve, per the scaling evidence in §3. Suitable preconditioners (multigrid, e.g. Newman & Alumbaugh, 2000; or a shifted-Laplacian/ILU-based approach) are a distinct, non-trivial research task and are explicitly **not** committed to by this ADR — they are the reason the in-house path remains research-only rather than a scheduled deliverable.

### 4.3 Field formulation, boundary conditions, air treatment, topography, mesh padding/refinement

**Decision: secondary/scattered-field (primary/secondary) formulation** for better accuracy near strong conductivity contrasts and thin features, **Dirichlet (zero scattered-field) conditions at padded outer boundaries**, **explicit air layers with small numerical conductivity**, all built on the **existing, already-3-D-capable `pycsamt.forward.maxwell.mesh` module** (`MaxwellMesh`, `MeshDesign`, `build_solver_mesh`, `SolverMeshModel` — confirmed in §2 to already handle 3-D padding, air/earth masking, and topography). No new mesh infrastructure is needed regardless of which backend (ModEM adapter now, or an in-house solver later) consumes it.

### 4.4 CPU/GPU, distributed batch strategy, licensing, checkpointing, expected cost per realization

- **CPU/GPU**: ModEM is CPU-only Fortran (no GPU path); this is acceptable since the adapter route (§4.1) does not require new numerical development.
- **Distributed batch strategy**: ModEM supports an MPI build (`make mpi`, per the vendored `README.txt`) for a single large 3-D solve; at the *dataset-generation* level (many independent realizations), `pycsamt.forward.maxwell.batch.solve_batch`'s existing `BatchPolicy.max_workers` thread pool (M4 PR4, already built) is the right layer for parallelizing across realizations, each of which can itself optionally use an MPI-parallel ModEM run.
- **Licensing**: open item, see §4.1.
- **Checkpointing**: already available for free at the realization level via `pycsamt.forward.maxwell.cache.MaxwellResultCache` (content-addressed, multi-process-safe) plus `pycsamt.forward.maxwell.batch.FailureManifest` (resumable failure tracking) — both already built in M4 PR4 and backend-agnostic, requiring no 3-D-specific work.
- **Expected cost per realization**: dominated by ModEM's own solve time (external to this codebase, not benchmarked here) rather than by anything this ADR controls; the in-house-solver cost estimate in §3 is the reason that path is rejected for production, not a cost commitment for the chosen backend.

## 5. Gate verdict

**Gate result: PASS for the external-backend path (§4.1); FAIL for an in-house production 3-D solver.** This is exactly the plan's own documented fallback ("If the gate fails, 3-D production training uses a trusted external backend; a new in-house solver remains research-only") — not a deviation from the plan.

## 6. Consequences / what this unblocks

- Per this ADR's original verdict, `pycsamt/forward/maxwell/mt3d.py` was not scheduled as a production deliverable. The user subsequently asked, explicitly and with the gate's conclusion already known, for both an in-house research prototype and an external adapter to be built, so users can choose ("sure why not build the two mt3d and modem3d.py so user will specify"). `mt3d.py` was therefore implemented on 2026-07-28 as a **research-only** module — `BackendCapabilities.maximum_cells` (default 6,000) enforces the small-grid restriction in code, not just in this document — and is not a reversal of the gate verdict above, which stands for the *production* recommendation.
- `mt3d.py` implements the governing equation from §1 as a genuine (not quasi-3-D) discrete curl-curl solve on a uniform Yee grid, and was validated in two stages: (1) the discrete curl operators independently, via exact linear-algebra identities (uniform/linear test fields, and `div(curl(E))=0` as an exact matrix product); (2) the assembled physics against `pycsamt.forward.maxwell.benchmarks.half_space_benchmark`, which it passes with real margin (~3% normalized RMS, ~1.5 degree phase error). It does **not** pass `layered_earth_benchmark` (30-45% error, confirmed by mesh refinement to be a boundary-condition bias, not discretization error) — a known, documented limitation of its single-layer boundary approximation (§4.3's scattered-field recommendation was not followed in `mt3d.py`, by design, to keep the research prototype's implementation risk low; see the module docstring for the total-field approximation actually used).
- `pycsamt/forward/maxwell/modem3d.py` (`ModEm3DAdapter`) was implemented on 2026-07-28, reusing `pycsamt.models.modem`'s existing `ModEmConfig`/`ModEmData`/`ModEmModel3D` file I/O via the `external.py` `BaseExternalMaxwellAdapter` foundation rather than a second parser. **This is input/output plumbing only, not physics-validated**: no compiled `Mod3DMT` binary exists in the environment it was built in, so — unlike `mt3d.py` — its generated files and response parsing could not be checked end-to-end against a live run or an analytic benchmark. `_VERIFIED_BENCHMARKS` is empty for exactly this reason; do not treat this adapter as numerically trustworthy until it has been run against a real compiled binary and checked against `half_space_benchmark`/`layered_earth_benchmark`. Its tests instead validate everything reachable without a live binary: real round trips through `ModEmData.read`/`ModEmModel3D.read`, and a full pipeline run against a scripted dependency-free stand-in executable. This work also surfaced and fixed a real, previously-uncovered bug in `pycsamt/models/modem/model3d.py::_parse_model3d` (a trailing-footer-detection heuristic that silently truncated the resistivity grid for any model with `nx<=3`, invisible before because the package's own tests for this path are all skipped in this environment).
- M7 ("Verified 3-D forward backend") should be re-scoped, if pursued, around validating `ModEm3DAdapter` against a real compiled `Mod3DMT` binary and the project's own analytic benchmarks (half-space/1-D-layered/2-D-extrusion/isolated-3-D-body) rather than around extending `mt3d.py`.
- M8 (correlated 3-D training) remains blocked on M7 regardless of path chosen here.
