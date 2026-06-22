ModEM Fortran Source — Build Instructions
==========================================

The ModEM v6.2.6 Fortran sources are already bundled in this directory:

  _source/2D/   — Mod2DMT (2-D MT serial)
  _source/3D/   — Mod3DMT (3-D MT serial + MPI)

You do NOT need the external ModEMv626/ folder to compile.
Just install the prerequisites and run make.

Step 1 — Install prerequisites
--------------------------------
Ubuntu/Debian:
  sudo apt install gfortran liblapack-dev libblas-dev

CentOS/RHEL:
  sudo yum install gcc-gfortran lapack-devel blas-devel

For MPI builds (3D only — recommended for production):
  sudo apt install libopenmpi-dev   # or: module load mpi (HPC clusters)

Step 2 — Compile
-----------------
Serial (2D):
  cd pycsamt/models/modem/_source/2D/
  make

Serial (3D):
  cd pycsamt/models/modem/_source/3D/
  make

MPI-parallel (3D — recommended for production runs):
  cd pycsamt/models/modem/_source/3D/
  make mpi

Intel compiler:
  make intel          # serial (2D or 3D)
  make intel_mpi      # MPI, 3D only

After compilation, copy or symlink the binary to a directory on PATH:
  ln -s $(pwd)/Mod3DMT ~/bin/Mod3DMT
  ln -s $(pwd)/Mod2DMT ~/bin/Mod2DMT   # 2D build

The pycsamt runner locates the binary via PATH or the absolute path
stored in ModEmConfig.binary_3d / ModEmConfig.binary_2d.

Updating sources (advanced)
-----------------------------
If you obtain a newer ModEM release, run install_sources.sh to refresh:

  bash pycsamt/models/modem/_source/install_sources.sh /path/to/newer/ModEM

This overwrites the bundled sources with the new version.

References
----------
Egbert, G.D. & Kelbert, A. (2012). Computational recipes for electromagnetic
inverse problems. Geophysical Journal International, 189(1), 251-267.

Kelbert, A., Meqbel, N., Egbert, G.D., Tandon, K. (2014). ModEM: A modular
system for inversion of electromagnetic geophysical data.
Computers & Geosciences, 66, 40-53.
