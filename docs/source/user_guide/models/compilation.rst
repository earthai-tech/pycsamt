.. _models_compilation:

Compiling the External Solvers
==============================

:term:`Occam2D`, :term:`ModEM`, and :term:`MARE2DEM` are external Fortran
programs, not Python packages -- pyCSAMT builds their input files, resolves
their executable, and launches them as a subprocess (see
:doc:`configuration_and_io`), but it does not ship a pre-compiled binary for
any of them. This page covers the one step that has to happen before any of
that: getting a working executable in the first place.

``pycsamt/models/_solver_build/`` holds a small set of bash scripts written
for exactly this. They work identically in intent on Windows, Linux, and
macOS, detect (and, if asked, install) a Fortran toolchain, and print exactly
what they are about to do before doing it -- nothing here silently
apt-installs or conda-creates anything.

.. list-table::
   :header-rows: 1
   :widths: 20 30 50

   * - Interface
     - Command
     - When to use it
   * - ``pycsamt build`` CLI
     - ``pycsamt build modem3d --auto-install``
     - **Recommended.** Works from any directory once pyCSAMT is
       installed -- no need to know or type the path below. Every
       ``build.sh``/``<solver>.sh`` option is accepted after the
       sub-command; ``pycsamt build <solver> --help`` shows the full
       list.
   * - ``make`` (repo checkout only)
     - ``make modem3d ARGS="--auto-install"``
     - A convenience for contributors already sitting in a git checkout
       with ``make`` on ``PATH``. Forwards to the exact same script; not
       shipped in the installed package.
   * - Scripts directly
     - ``pycsamt/models/_solver_build/modem3d.sh --auto-install``
     - What both interfaces above actually call. Useful before pyCSAMT
       itself is installed (e.g. bootstrapping a toolchain ahead of
       ``pip install -e .``), or to run a script from a source checkout
       whose path you already have open.

The three are equivalent -- the CLI and ``make`` both resolve to the same
script with the same arguments, so anything documented below in terms of
the script's own flags applies to all three.

.. list-table::
   :header-rows: 1
   :widths: 16 22 20 42

   * - Script
     - Builds
     - Source
     - Toolchain
   * - ``modem2d.sh``
     - ModEM ``Mod2DMT`` (2-D)
     - vendored, committed
     - gfortran + make (+ LAPACK/BLAS)
   * - ``modem3d.sh``
     - ModEM ``Mod3DMT`` (3-D)
     - vendored, committed
     - gfortran + make (+ LAPACK/BLAS)
   * - ``occam2d.sh``
     - Occam2DMT
     - vendored, committed
     - gfortran + make
   * - ``mare2dem.sh``
     - MARE2DEM
     - downloaded on demand
     - Intel ``mpiifort``/``mpiicc`` + MKL, Linux/macOS/WSL only

A ``build.sh`` dispatcher forwards to whichever of these you name -- it is
what ``pycsamt build <solver>`` itself calls into:

.. code-block:: bash

   pycsamt/models/_solver_build/build.sh modem3d --auto-install
   # is equivalent to
   pycsamt/models/_solver_build/modem3d.sh --auto-install
   # is equivalent to
   pycsamt build modem3d --auto-install

Everything below was run for real while writing this page -- the terminal
output is not invented.

Quick Start
-----------

.. code-block:: bash

   # ModEM 3-D, installing a toolchain automatically if one is missing
   pycsamt build modem3d --auto-install

   # ModEM 2-D, assuming a toolchain is already present
   pycsamt build modem2d

   # Occam2D
   pycsamt build occam2d --auto-install -y

   # MARE2DEM: status only -- nothing is downloaded or built without
   # --auto-install (see "MARE2DEM" below for why this one is different)
   pycsamt build mare2dem

Every script prints ``--help`` for its own flags (``--auto-install``,
``-y``/``--yes``, ``--clean``, ``--mpi``, ``--intel``, ``--prefix``), which
are consistent across the three ModEM/Occam2D scripts, and
``pycsamt build <solver> --help`` forwards to that same help text rather
than duplicating it:

.. code-block:: console

   $ pycsamt build occam2d --help
   Author: LKouadio <etanoyau@gmail.com>
   License: LGPL-3.0

   Build Occam2DMT (deGroot-Hedlin & Constable) from the vendored
   source at pycsamt/models/occam2d/_source/.

   Usage:
     pycsamt/models/_solver_build/occam2d.sh [options]

   Options:
     --auto-install     Install a missing Fortran toolchain (see
                         lib/common.sh for exactly what that means per
                         OS). Off by default. Occam2D itself does not
                         need LAPACK/BLAS (its own linear algebra is
                         self-contained; a commented-out fast-LAPACK
                         path exists in Occam.f90 for anyone who wants
                         to enable it by hand), but the shared
                         auto-installer provisions the same toolchain
                         as the ModEM scripts for consistency.
     -y, --yes          Skip the confirmation prompt --auto-install
                         would otherwise show.
     --clean            Remove previous object/module files and the
                         binary before building.
     --prefix DIR        After a successful build, also copy the
                         binary into DIR.
     -h, --help          Show this help and exit.

   The vendored Makefile's own default compiler ("f90") is a legacy
   alias unlikely to exist on a modern system; this script always
   overrides it to the detected gfortran.

   Exit status: 0 on a successful build, 1 otherwise.

Auto-Install Behavior
---------------------

Nothing runs silently: every install command is printed first and needs a
``y`` confirmation unless ``-y``/``--yes`` is also given.

* **Linux** -- ``apt-get``/``dnf``/``yum install gfortran`` plus LAPACK/BLAS
  development packages, whichever package manager is found first. Runs
  without ``sudo`` automatically if already root (e.g. inside a minimal
  container) or if ``sudo`` itself is not installed.
* **macOS** -- ``brew install gcc openblas`` (Homebrew's ``gcc`` formula
  includes ``gfortran``). Toolchain *detection* also searches versioned
  binaries (``gfortran-15`` down to ``gfortran-9``), since Homebrew does not
  always leave an unversioned ``gfortran`` symlink on ``PATH`` when more than
  one GCC version is installed side by side.
* **Windows** -- gfortran and ``make`` are not reliably available outside a
  dedicated MinGW-w64 toolchain, so an **isolated** conda environment
  (``pycsamt-fortran`` by default, override with ``PYCSAMT_FORTRAN_ENV``) is
  created with ``m2w64-gcc-fortran``, ``m2w64-openblas``, and ``make`` from
  conda-forge. It does not touch your existing pycsamt environment, and you
  never need to ``conda activate`` it -- the scripts resolve its binaries by
  absolute path.

.. _modem_compilation:

ModEM
-----

ModEM's source is vendored under ``pycsamt/models/modem/_source/{2D,3D}/``.
Both Makefiles had two real, platform-independent bugs (they would have
failed identically on Linux or macOS, not just here) that these scripts'
first real build run surfaced and fixed in place:

#. A missing build rule for ``FIELDS/FiniteDiff3D/sg_spherical.f90``
   (3-D only) -- it is ``use``\ d directly by ``GridCalc.f90`` but was never
   in the object list at all.
#. ``MPI/Declaration_MPI.f90``, ``{2D,3D}_MT/Sub_MPI.f90``, and
   ``MPI/Main_MPI.f90`` were never compiled, despite being unconditionally
   ``use``\ d by several files (``INVcore.f90``, ``DCG.f90``,
   ``Mod{2,3}DMT.f90``, ``SymmetryTest.f90``, ``ModelSpace.f90``)
   *regardless* of whether MPI is enabled -- their entire module bodies are
   wrapped in ``#ifdef MPI`` / ``#endif``, and the serial build never passed
   ``-cpp`` to strip those guards.

Because the vendored object lists are not a strict topological sort of
Fortran module dependencies (a pre-existing property of the Makefiles as
adapted from upstream, not something these scripts try to reorder), a
single-pass ``make`` can still fail with "module file not found" purely
because of build order. ``modem2d.sh``/``modem3d.sh`` retry with
``make -k`` (keep-going) for up to six passes -- never ``-j`` (parallel),
since that would make the already-imperfect ordering non-deterministic
instead of merely sub-optimal.

ModEM 3-D
~~~~~~~~~

.. code-block:: console

   $ pycsamt build modem3d --clean
   ==> Building ModEM Mod3DMT from .../pycsamt/models/modem/_source/3D
   [info]  Using F90=.../envs/pycsamt-fortran/Library/mingw-w64/bin/gfortran.exe  make=.../Library/bin/make.exe
   ==> make clean
   rm -rf objs Mod3DMT
   ==> Compiling (this can take a minute; a few make -k passes is normal)
   [info]  make pass 1/6: 64 error line(s) remaining
   [info]  make pass 2/6: 52 error line(s) remaining
   [info]  make pass 3/6: 24 error line(s) remaining
   [info]  make pass 4/6: 0 error line(s) remaining
   [info]  Copied 5 runtime DLL(s) into .../3D for a self-contained build.
   ==> Built: .../pycsamt/models/modem/_source/3D/Mod3DMT.exe
   [info]  Point pycsamt at it via ModEmConfig(binary_3d="..."), or leave the default
   [info]  "Mod3DMT" and rely on this directory being pycsamt's own search-path fallback.
   [warn]  Remember: models need >= 10 earth z-cells (ModEM's own hardcoded
   [warn]  air-layer count) or the solver will crash -- see this script's header.

Four ``make -k`` passes to converge is normal and expected, not a sign
something is wrong -- each pass produces ``.mod`` files a later file in the
list depends on. Point :class:`~pycsamt.models.modem.ModEmConfig` at the
result the same way :doc:`modem` already shows:

.. code-block:: pycon

   >>> from pycsamt.models.modem import ModEmConfig, ModEmRunner
   >>> cfg = ModEmConfig(mode="3d", binary_3d="Mod3DMT")
   >>> runner = ModEmRunner("runs/modem_3d_v01/native", config=cfg)

The final warning above is worth taking seriously: ModEM's own Fortran (not
pyCSAMT's code) hardcodes 10 air layers for the 3-D solver, and its default
"mirror" air-layer sizing method reads that many earth-layer widths from the
model with no bounds check against the actual earth cell count. A model with
fewer than 10 earth z-cells will build and launch fine, then crash mid-solve
with ``Error: b in QMR contains NaNs``.
:class:`~pycsamt.forward.maxwell.modem3d.ModEm3DAdapter` already guards
against this on the AI-inversion side, but it applies equally if you drive
``Mod3DMT`` directly.

ModEM 2-D
~~~~~~~~~

The 2-D build needs the same MPI-stub fix and nothing else -- its object
list was already in a working order once the missing files were added, so it
converges in a single pass:

.. code-block:: console

   $ pycsamt build modem2d --clean
   ==> Building ModEM Mod2DMT from .../pycsamt/models/modem/_source/2D
   [info]  Using F90=.../gfortran.exe  make=.../make.exe
   ==> make clean
   rm -rf objs Mod2DMT
   ==> Compiling (this can take a minute; a few make -k passes is normal)
   [info]  make pass 1/6: 0 error line(s) remaining
   [info]  Copied 5 runtime DLL(s) into .../2D for a self-contained build.
   ==> Built: .../pycsamt/models/modem/_source/2D/Mod2DMT.exe
   [info]  Point pycsamt at it via ModEmConfig(binary_2d="..."), or leave the default
   [info]  "Mod2DMT" and rely on this directory being pycsamt's own search-path fallback.

MPI and Intel builds
~~~~~~~~~~~~~~~~~~~~

Both ModEM scripts accept ``--mpi`` (needs ``mpif90`` already on ``PATH`` --
not something ``--auto-install`` provisions on any platform) and ``--intel``
(needs ``ifort`` already on ``PATH`` -- a commercial compiler,
``--auto-install`` never installs one):

.. code-block:: bash

   pycsamt build modem3d --mpi
   pycsamt build modem3d --intel

.. _occam2d_compilation:

Occam2D
-------

Occam2D's vendored Makefile (``pycsamt/models/occam2d/_source/Makefile``)
needed no source fix at all -- only a compiler override, since its own
default (``FC90 = f90``) is a legacy alias unlikely to exist on a modern
system. Its four source files compile and link in a single ``gfortran``
invocation:

.. code-block:: console

   $ pycsamt build occam2d --clean
   ==> Building Occam2D from .../pycsamt/models/occam2d/_source
   [info]  Using F90=.../gfortran.exe  make=.../make.exe
   ==> make clean
   rm -f *.o *~ core *.mod
   rm -f Occam2D
   ==> Compiling
   [info]  make pass 1/3: 0 error line(s) remaining
   [info]  Copied 5 runtime DLL(s) into .../_source for a self-contained build.
   ==> Built: .../pycsamt/models/occam2d/_source/Occam2D.exe
   [info]  Point pycsamt.models.occam2d.runner at this binary's path.

Occam2D itself does not need LAPACK/BLAS (a commented-out fast-LAPACK path
exists in ``Occam.f90`` for anyone who wants to enable it by hand), but
``occam2d.sh`` still provisions the same toolchain as the ModEM scripts for
consistency, rather than adding a second, narrower detection path.

.. _mare2dem_compilation:

MARE2DEM
--------

MARE2DEM is different in kind, not just in build recipe, so
``mare2dem.sh`` is a thin wrapper around the pre-existing
:class:`~pycsamt.models.mare2dem.SourceManager` rather than a
reimplementation:

* its source is **not vendored** (``pycsamt/models/mare2dem/_source/`` ships
  only a ``.gitkeep`` -- the real tree is roughly 49 MB and under its own
  license), so it has to be downloaded first;
* it genuinely needs an Intel MPI Fortran/C toolchain and the Intel MKL for
  ScaLAPACK/BLACS -- there is no generic gfortran/conda path for it the way
  there is for ModEM and Occam2D;
* it **cannot be built natively on Windows at all** --
  ``SourceManager.build()`` raises immediately on ``sys.platform ==
  "win32"``.

Because of this, ``mare2dem.sh`` always reports status first and only
downloads or builds when ``--auto-install`` is given -- downloading ~49 MB
and compiling a large MPI codebase is a heavier action than the other three
scripts perform by default. On Windows it stops immediately with WSL
guidance, which is exactly what happened while writing this page:

.. code-block:: console

   $ pycsamt build mare2dem
   [error] MARE2DEM cannot be built natively on Windows (Intel MKL/ScaLAPACK
   toolchain is Linux/macOS-oriented). Run this script from inside WSL2
   instead: wsl bash pycsamt/models/_solver_build/mare2dem.sh

On Linux, macOS, or inside WSL2, with Intel oneAPI already sourced
(``source /opt/intel/oneapi/setvars.sh``):

.. code-block:: bash

   pycsamt build mare2dem --auto-install -y

which runs the same two calls documented on
:class:`~pycsamt.models.mare2dem.SourceManager`:

.. code-block:: pycon

   >>> from pycsamt.models.mare2dem import SourceManager
   >>> sm = SourceManager(verbose=1)
   >>> sm.download()          # clones from Bitbucket, or falls back to a tar.gz
   >>> sm.build()              # generates a Make include, compiles, locates the binary

**This has now genuinely been done, inside WSL2 Ubuntu with Intel's free
oneAPI toolchain (2026.1)** -- everything below is real output from that
build, not a description of a remaining risk. It surfaced several real bugs
in ``SourceManager`` itself, all fixed in place; the sections after this one
walk through what they were.

Setting up the toolchain
~~~~~~~~~~~~~~~~~~~~~~~~

Intel oneAPI is not preinstalled on a plain Ubuntu (or WSL2 Ubuntu) image.
Installed from Intel's own apt repository:

.. code-block:: bash

   wget -qO- https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB \
       | gpg --dearmor | sudo tee /usr/share/keyrings/oneapi-archive-keyring.gpg > /dev/null
   echo "deb [signed-by=/usr/share/keyrings/oneapi-archive-keyring.gpg] https://apt.repos.intel.com/oneapi all main" \
       | sudo tee /etc/apt/sources.list.d/oneAPI.list
   sudo apt-get update
   sudo apt-get install -y intel-oneapi-compiler-fortran intel-oneapi-mpi-devel intel-oneapi-mkl-devel

   source /opt/intel/oneapi/setvars.sh

``setvars.sh`` has to be sourced in every new shell that will build or run
MARE2DEM -- it puts ``mpiifx``/``mpiicx``/``ifx``/``icx`` and the MKL
runtime libraries on ``PATH``/``LD_LIBRARY_PATH``. ``SourceManager`` does
not source it for you (it cannot -- that has to happen in the calling
shell before Python even starts), which is why ``print_status()``/``build()``
report an unusable toolchain if you forget this step.

Real status, real build
~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: pycon

   >>> from pycsamt.models.mare2dem import Mare2DEMConfig, SourceManager
   >>> cfg = Mare2DEMConfig(source_dir="/root/build/mare2dem")
   >>> sm = SourceManager(config=cfg, verbose=1)
   >>> sm.print_status()
   MARE2DEM SourceManager status
   ────────────────────────────────────────
     source_dir  : /root/build/mare2dem
     downloaded  : True
     built       : False
     binary      : (not found)
     FC compiler : mpiifx
     CC compiler : mpiicx
     MKLROOT     : /opt/intel/oneapi/mkl/2026.1

``FC compiler``/``CC compiler`` read ``mpiifx``/``mpiicx`` here, not the
``mpiifort``/``mpiicc`` this page used to show -- see `A real toolchain bug,
found and fixed`_ below for why that distinction matters.

.. code-block:: pycon

   >>> binary = sm.build()

.. code-block:: console

   Building MARE2DEM (this may take several minutes) …
     FC = mpiifx
     CC = mpiicx
     MKLROOT = /opt/intel/oneapi/mkl/2026.1
     Include file = /root/build/mare2dem/_pycsamt_build/pycsamt_auto.inc

   mpiicx -O2 -fPIC -std=gnu89  -c -o triangle.o -DTRILIBRARY -O2 -fPIC triangle.c
   triangle.c:12871:6: warning: a function definition without a prototype is
   deprecated in all versions of C and is not supported in C23 [-Wdeprecated-non-prototype]
   ...
   115 warnings generated.
   mpiifx -cxxlib -O2 -fpp -fPIC -c -o mt1d.o mt1d.f90
   mpiifx -cxxlib -O2 -fpp -fPIC -c -o kx_io.o kx_io.f90
   mpiifx -cxxlib -O2 -fpp -fPIC -c -o em2dkx.o em2dkx.f90
   mpiifx -cxxlib -O2 -fpp -fPIC -c -o dc2dkx.o dc2dkx.f90
   mpiifx -cxxlib -O2 -fpp -fPIC -c -o mare2dem_scalapack.o mare2dem_scalapack.f90
   mpiifx -cxxlib -O2 -fpp -fPIC -c -o occam.o occam.f90
   mpiicx -O2 -fPIC -std=gnu89  -c -o c_fortran_triangle.o c_fortran_triangle.c
   mpiifx -cxxlib -O2 -fpp -fPIC -c -o filtermodules.o filtermodules.f90
   mpiifx -cxxlib -O2 -fpp -fPIC -c -o mare2dem_common.o mare2dem_common.f90
   mpiifx -cxxlib -O2 -fpp -fPIC -c -o spline_kx_module.o spline_kx_module.f90
   mpiifx -cxxlib -O2 -fpp -fPIC -c -o mare2dem_worker.o mare2dem_worker.f90
   mpiifx -cxxlib -O2 -fpp -fPIC -c -o mare2dem_mpi.o mare2dem_mpi.f90
   mpiifx -cxxlib -O2 -fpp -fPIC -c -o mare2dem_penaltymatrix.o mare2dem_penaltymatrix.f90
   mpiifx -cxxlib -O2 -fpp -fPIC -c -o mare2dem_io.o mare2dem_io.f90
   mpiifx -cxxlib -O2 -fpp -fPIC -c -o em2d.o em2d.f90
   mpiifx -cxxlib -O2 -fpp -fPIC -c -o runmare2dem.o runmare2dem.f90
   mpiifx -cxxlib -O2 -fpp -fPIC em_constants.o kdtree2.o ... runmare2dem.o \
       ./libraries/scalapack-2.2.0/libscalapack.a \
       -L/opt/intel/oneapi/mkl/2026.1/lib -Wl,-rpath,/opt/intel/oneapi/mkl/2026.1/lib \
       -I/opt/intel/oneapi/mkl/2026.1/include -lmkl_intel_lp64 -lmkl_sequential \
       -lmkl_core -lpthread -lm -ldl -o MARE2DEM
   #
   # All done!
   #

   MARE2DEM binary ready: /root/build/mare2dem/MARE2DEM

The ``triangle.c`` warnings (``-Wdeprecated-non-prototype``, 115 of them) are
expected and harmless -- Triangle's own C source predates modern function-
prototype conventions by decades; ``icx`` warns about it but still compiles
it correctly. ScaLAPACK's BLACS layer (``./libraries/scalapack-2.2.0/``, not
shown above since it only needs building once and was already cached here)
produces a similar, equally harmless volume of ``-std=gnu89`` compatibility
output the first time it is built. Point
:class:`~pycsamt.models.mare2dem.Mare2DEMConfig` at the result the same way
:doc:`mare2dem` already shows:

.. code-block:: pycon

   >>> cfg = Mare2DEMConfig(binary=str(binary), use_mpi=True, n_procs=4)

A real physics check, not just a compile check
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A binary that compiles is not the same claim as a binary that produces
correct physics. ``pycsamt/forward/tests/test_maxwell_mare2dem.py`` has a
``requires_real_mare2dem``-gated test that only runs when
``PYCSAMT_MARE2DEM_BINARY`` points at a real executable (and ``mpirun`` is
on ``PATH``) -- it solves a real half-space forward problem through
:class:`~pycsamt.forward.maxwell.mare2dem.Mare2DEMAdapter` and checks the
result against the analytic
:func:`~pycsamt.forward.maxwell.benchmarks.half_space_impedance` formula:

.. code-block:: console

   $ export PYCSAMT_MARE2DEM_BINARY=/root/build/mare2dem/MARE2DEM
   $ python -m pytest pycsamt/forward/tests/test_maxwell_mare2dem.py -k real -v
   pycsamt/forward/tests/test_maxwell_mare2dem.py::test_real_mare2dem_passes_half_space_benchmark PASSED
   ======================= 1 passed in 1.13s =======================

This is the adapter-level counterpart to the compile-level bugs below --
getting MARE2DEM to *build* surfaced the toolchain issues on this page;
getting a real solve to match known physics surfaced four further, genuinely
different bugs in how :class:`~pycsamt.forward.maxwell.mare2dem.Mare2DEMAdapter`
talks to the compiled binary (mesh input format, resistivity units, required
filename conventions). None of those are build issues -- see that module's
own docstring for the full account, since they belong there, not here.

A real toolchain bug, found and fixed
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Before this real build, ``SourceManager`` detected and generated Make
include flags for the *classic* Intel toolchain
(``mpiifort``/``mpiicc``/``xiar``, the same names the vendored example
include files ``habanero.inc``/``macos.inc`` use). That toolchain is gone or
broken on a current oneAPI installation (2024.0+):

* ``mpiifort`` is either absent, or present but internally tries to exec
  the now-removed classic ``ifort`` and fails with ``eval: ifort: not
  found``.
* ``xiar`` (the classic Intel archiver) no longer exists at all.
* The vendored ``TRICOPTS`` example (``-fp-model precise -fp-model
  source``) is classic-``icc``-only syntax; ``icx`` rejects it outright
  (``unsupported argument 'source' to option '-ffp-model='``).

None of this is a documentation gap -- it is exactly the kind of thing that
only a real build run finds, and it would have failed identically for
anyone following the previous version of this page on a current oneAPI
install. Fixed in ``pycsamt/models/mare2dem/source.py``:

#. Compiler detection now tries ``mpiifx``/``mpiicx`` (the modern wrappers
   for ``ifx``/``icx``) first, falling back to the classic
   ``mpiifort``/``mpiicc`` only for older oneAPI installations where they
   still genuinely work.
#. The generated ``ARCH`` line is always plain ``ar`` now, not conditionally
   ``xiar`` -- GNU ``ar`` links ``ifx``-compiled object files on Linux/macOS
   fine, confirmed by the build above.
#. ``TRICOPTS`` is simplified to ``-O2 -fPIC`` unconditionally; the
   fp-model flags were never required for correctness, just inherited from
   the vendored example.
#. Two flags this page's previous version did not mention at all, both
   required for a successful link/compile and both visible in the real
   transcript above: ``-cxxlib`` in ``FFLAGS`` (``ifx``-specific; the final
   link fails without it) and ``-std=gnu89`` in ``CFLAGS`` (the vendored
   BLACS/ScaLAPACK C sources rely on 1990s-style implicit function
   declarations that both ``icx`` and modern GCC reject under their C99+
   defaults otherwise -- ``error: call to undeclared function
   'BI_imvcopy'``, among others, without it).

Two further, smaller fixes worth knowing if you build on a different
system: ``call_triangle.f90``'s ``#IF DEFINED(...)``-style preprocessor
directives are uppercase, which only Intel's ``-fpp`` accepts case-
insensitively -- GNU-style ``-cpp`` does not, so a gfortran-based attempt at
this build fails on that file specifically, which is also why there is no
gfortran/conda fallback path for MARE2DEM the way there is for ModEM and
Occam2D. And ``MKLLIB`` now always adds ``-Wl,-rpath,<mklroot>/lib``, so the
built binary finds its MKL runtime libraries without needing
``setvars.sh``/``LD_LIBRARY_PATH`` sourced again at *run* time, only at
*build* time.

.. _solver_windows_binaries:

Self-Contained Windows Binaries
-------------------------------

A MinGW-w64-built ``.exe`` needs its own runtime DLLs either on ``PATH`` or
sitting next to itself -- Windows checks the executable's own directory
first. All four scripts copy the handful of DLLs the ``pycsamt-fortran``
conda environment provides (``libgcc_s_seh-1.dll``, ``libgfortran-3.dll``,
``libopenblas.dll``, ``libquadmath-0.dll``, ``libwinpthread-1.dll``)
alongside the binary automatically, so it runs standalone without activating
anything -- confirmed by running each built ``.exe`` from a plain shell with
none of the conda environment's directories on ``PATH``.

.. _solver_compilation_troubleshooting:

Troubleshooting
---------------

**"internal compiler error: Segmentation fault", inconsistently, always on
the same files (Windows)** -- this is a real bug that was caught and fixed
while first writing these scripts, not a description of a remaining risk:
``conda env list`` prints Windows-style backslash paths
(``C:\Users\...\envs\pycsamt-fortran``) regardless of shell, and splicing
that literally into a colon-separated bash ``PATH`` string is invalid -- the
colon right after the drive letter is read as a second ``PATH`` entry
separator, corrupting ``PATH`` silently. gfortran internally re-execs its
own ``cc1``/``collect2``/``as``/``ld`` helpers by searching ``PATH``, so a
corrupted ``PATH`` let it pick up mismatched helpers from elsewhere on the
system -- on specific, larger source files, not randomly. ``lib/common.sh``
now converts every Windows path it resolves to POSIX form
(``/c/Users/...``) before it ever reaches ``PATH``. If you see this error
again, suspect a raw backslash in ``$PATH`` before suspecting the Fortran
source.

**"module file not found" during a ModEM build that never resolves** -- this
should not happen after the Makefile fixes above, but if it does, check that
``make`` is actually being invoked without ``-j``/parallel jobs (these
scripts never add one, but a wrapping build system might) -- the object list
order depends on a serial, retried build.

**"b in QMR contains NaNs" while running a compiled ``Mod3DMT``, not while
building it** -- see the "Remember" warning in the ModEM 3-D section above.
This is a model-configuration issue (fewer than 10 earth z-cells), not a
build problem.

**gfortran not found on macOS after ``brew install gcc``** -- Homebrew often
installs a versioned binary (``gfortran-14``) without an unversioned
``gfortran`` symlink when multiple GCC versions coexist. Re-run the build
script; its own detection already searches versioned names, so this
typically resolves itself on retry once ``brew`` has finished linking.

**``eval: ifort: not found`` while building MARE2DEM** -- ``mpiifort`` is
present on ``PATH`` but is a wrapper for the classic ``ifort``, which
current oneAPI releases no longer ship. ``SourceManager`` now prefers
``mpiifx`` automatically (see `A real toolchain bug, found and fixed`_
above); if you still see this, your oneAPI installation predates the fix
being able to find ``mpiifx`` at all -- install ``intel-oneapi-compiler-
fortran`` from a current oneAPI release, or pass an explicit
``inc_file=`` to :meth:`~pycsamt.models.mare2dem.SourceManager.build` for
your cluster's own toolchain.

**``error: call to undeclared function 'BI_imvcopy'`` (or similar) while
building MARE2DEM's ScaLAPACK dependency** -- a pre-2026 ``SourceManager``
did not pass ``-std=gnu89`` to the C compiler; the vendored BLACS sources
rely on implicit function declarations that both ``icx`` and modern GCC
reject by default under C99+. Already fixed in ``_generate_inc`` -- if you
see this with a current pycsamt, you are likely using a hand-written
``inc_file=`` that needs the same flag added.

**``unsupported argument 'source' to option '-ffp-model='`` while building
MARE2DEM** -- the vendored example include files' ``TRICOPTS`` value
(``-fp-model precise -fp-model source``) is classic-``icc``-only syntax;
already dropped from the auto-generated include file (see above). Only
relevant if you are editing a Make include file by hand.

Related Reading
---------------

* :doc:`modem`, :doc:`occam2d`, :doc:`mare2dem` for using each engine once a
  binary exists: building native files, configuring the runner, and loading
  results.
* :doc:`configuration_and_io` for how each ``*Config``/``*Runner`` pair
  resolves an executable from ``PATH``, the run directory, and each engine's
  own ``_source`` folder.
* ``pycsamt/models/_solver_build/README.md`` for the scripts' own reference
  documentation (flags, environment variables, and the full bug list) kept
  alongside the code they document.
* ``pycsamt build --help`` for the packaged CLI's own command listing, and
  the repo-root ``Makefile`` for the ``make <solver>`` convenience targets
  used by contributors working from a checkout.
