# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Solver-neutral contracts and adapters for multidimensional Maxwell physics.

This package is the integration boundary for verified 2-D and 3-D EM solver
backends.  It will own common problem/result contracts, backend capability
reporting, convergence diagnostics, and canonical benchmark definitions.

Existing :mod:`pycsamt.forward.em2d` and :mod:`pycsamt.forward.em3d` APIs are
not re-exported here.  A backend is exposed through this namespace only after
its dimensional coupling and numerical validation gates are documented and
tested.  Importing the package must not import optional solver dependencies.
"""

from .adapters import (
    AdapterPolicy,
    BackendExecutionError,
    BaseMaxwellAdapter,
    CallableMaxwellAdapter,
    IncompatibleProblemError,
    InvalidBackendResultError,
    MaxwellAdapterError,
    SolverConvergenceError,
)
from .backends import (
    BackendCapabilities,
    BackendRegistration,
    BackendRegistry,
    CompatibilityReport,
    MaxwellBackend,
    backend_registry,
    create_backend,
    list_backends,
    register_backend,
    unregister_backend,
)
from .batch import (
    BatchAbortedError,
    BatchPolicy,
    BatchReport,
    FailureManifest,
    ProblemFailure,
    solve_batch,
)
from .benchmarks import (
    BenchmarkMetrics,
    BenchmarkOutcome,
    BenchmarkReport,
    BenchmarkThresholds,
    MaxwellBenchmark,
    half_space_benchmark,
    half_space_impedance,
    layered_earth_benchmark,
    layered_earth_impedance,
    run_benchmarks,
)
from .cache import (
    CacheCorruptionError,
    CacheEntry,
    CacheLockTimeoutError,
    CacheStatistics,
    MaxwellResultCache,
)
from .contracts import (
    ForwardResult,
    MaxwellMesh,
    MaxwellProblem,
    ReceiverSet,
    SolverDiagnostics,
)
from .contracts_tri import TriMesh, TriProblem
from .external import (
    BaseExternalMaxwellAdapter,
    ExecutableNotFoundError,
    ExternalProcessError,
    ExternalRunPolicy,
    ExternalRunResult,
    make_availability_probe,
    probe_executable_version,
    resolve_executable,
)
from .mesh import (
    MeshDesign,
    MeshQuality,
    SolverMeshModel,
    build_solver_mesh,
    skin_depth_m,
)
from .mare2dem import Mare2DEMAdapter, register_mare2dem_backend
from .modem3d import ModEm3DAdapter, register_modem3d_backend
from .mt2d import MT2DAdapter, register_mt2d_backend
from .mt3d import MT3DAdapter
from .tri_fem2d import TriFEM2DAdapter, register_trifem2d_backend
from .tri_mesh_gen import build_graded_tri_mesh

__all__ = [
    "MaxwellMesh",
    "ReceiverSet",
    "MaxwellProblem",
    "SolverDiagnostics",
    "ForwardResult",
    "TriMesh",
    "TriProblem",
    "BackendCapabilities",
    "CompatibilityReport",
    "MaxwellBackend",
    "BackendRegistration",
    "BackendRegistry",
    "backend_registry",
    "register_backend",
    "unregister_backend",
    "create_backend",
    "list_backends",
    "MaxwellAdapterError",
    "IncompatibleProblemError",
    "BackendExecutionError",
    "InvalidBackendResultError",
    "SolverConvergenceError",
    "AdapterPolicy",
    "BaseMaxwellAdapter",
    "CallableMaxwellAdapter",
    "MeshDesign",
    "MeshQuality",
    "SolverMeshModel",
    "skin_depth_m",
    "build_solver_mesh",
    "CacheCorruptionError",
    "CacheLockTimeoutError",
    "CacheEntry",
    "CacheStatistics",
    "MaxwellResultCache",
    "BenchmarkThresholds",
    "BenchmarkMetrics",
    "BenchmarkOutcome",
    "MaxwellBenchmark",
    "BenchmarkReport",
    "half_space_impedance",
    "layered_earth_impedance",
    "half_space_benchmark",
    "layered_earth_benchmark",
    "run_benchmarks",
    "ExecutableNotFoundError",
    "ExternalProcessError",
    "ExternalRunPolicy",
    "ExternalRunResult",
    "BaseExternalMaxwellAdapter",
    "resolve_executable",
    "probe_executable_version",
    "make_availability_probe",
    "MT2DAdapter",
    "register_mt2d_backend",
    "MT3DAdapter",
    "ModEm3DAdapter",
    "register_modem3d_backend",
    "Mare2DEMAdapter",
    "register_mare2dem_backend",
    "TriFEM2DAdapter",
    "register_trifem2d_backend",
    "build_graded_tri_mesh",
    "BatchAbortedError",
    "BatchPolicy",
    "BatchReport",
    "FailureManifest",
    "ProblemFailure",
    "solve_batch",
]
