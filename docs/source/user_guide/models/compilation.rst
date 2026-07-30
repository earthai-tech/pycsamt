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
* it genuinely needs Intel ``mpiifort``/``mpiicc`` and the Intel MKL for
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

There is no compiled ``MARE2DEM`` binary in the environment this page was
written in (Intel oneAPI is not installed here), so the download/build call
above is not shown with fabricated output -- run it yourself once the
prerequisites are in place, and ``sm.print_status()`` will confirm what
happened.

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
