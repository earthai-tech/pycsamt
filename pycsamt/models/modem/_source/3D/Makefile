# Makefile for Mod3DMT (ModEM 3D MT inversion)
# Adapted from ModEM v6.2.6 distribution (Egbert & Kelbert, 2012)
#
# Usage:
#   make            — serial build with gfortran (default)
#   make intel      — serial build with ifort
#   make mpi        — MPI-parallel build with mpif90  ← recommended for 3D
#   make clean      — remove object files and binary
#
# Prerequisites:
#   gfortran >= 6.0 (or ifort / mpif90 for MPI)
#   LAPACK + BLAS (system packages: liblapack-dev libblas-dev)
#   For MPI build: OpenMPI or IntelMPI headers and libraries
#
# Output:  ./Mod3DMT

# ---- compiler settings ----
F90      ?= gfortran
FFLAGS   ?= -O3 -ffree-line-length-none
MPIFLAGS ?= -x f95-cpp-input
LIBS     ?= -llapack -lblas

# ---- directories ----
SRCDIR   = .
OBJDIR   = objs
SRC3D    = $(SRCDIR)/3D_MT
UTILS    = $(SRCDIR)/UTILS
INV      = $(SRCDIR)/INV
SENS     = $(SRCDIR)/SENS
FIELDS   = $(SRCDIR)/FIELDS
LAPACK   = $(SRCDIR)/LAPACK
MPI      = $(SRCDIR)/MPI

MODULE   = -J$(OBJDIR)
OUT      = Mod3DMT

# ---- object list ----
OBJ = \
  $(OBJDIR)/math_constants.o \
  $(OBJDIR)/utilities.o \
  $(OBJDIR)/file_units.o \
  $(OBJDIR)/polpak.o \
  $(OBJDIR)/elements.o \
  $(OBJDIR)/sg_scalar.o \
  $(OBJDIR)/sg_sparse_vector.o \
  $(OBJDIR)/sg_vector.o \
  $(OBJDIR)/sg_boundary.o \
  $(OBJDIR)/sg_diff_oper.o \
  $(OBJDIR)/DataSpace.o \
  $(OBJDIR)/GridDef.o \
  $(OBJDIR)/GridCalc.o \
  $(OBJDIR)/ModelSpace.o \
  $(OBJDIR)/SensMatrix.o \
  $(OBJDIR)/EMfieldInterp.o \
  $(OBJDIR)/transmitters.o \
  $(OBJDIR)/SolnSpace.o \
  $(OBJDIR)/receivers.o \
  $(OBJDIR)/dataTypes.o \
  $(OBJDIR)/DataFunc.o \
  $(OBJDIR)/DataSens.o \
  $(OBJDIR)/SolverSens.o \
  $(OBJDIR)/WSfwd2Dmod.o \
  $(OBJDIR)/WSfwd1Dmod.o \
  $(OBJDIR)/WSfwd2Dpar.o \
  $(OBJDIR)/WSutils.o \
  $(OBJDIR)/FwdTEmod.o \
  $(OBJDIR)/FwdTMmod.o \
  $(OBJDIR)/boundary_ws.o \
  $(OBJDIR)/modelOperator3D.o \
  $(OBJDIR)/solver.o \
  $(OBJDIR)/EMsolve3D.o \
  $(OBJDIR)/nestedEM.o \
  $(OBJDIR)/ForwardSolver.o \
  $(OBJDIR)/SensComp.o \
  $(OBJDIR)/SymmetryTest.o \
  $(OBJDIR)/UserCtrl.o \
  $(OBJDIR)/ioAscii.o \
  $(OBJDIR)/DataIO.o \
  $(OBJDIR)/Main.o \
  $(OBJDIR)/INVcore.o \
  $(OBJDIR)/NLCG.o \
  $(OBJDIR)/DCG.o \
  $(OBJDIR)/Mod3DMT.o

all: $(OBJDIR) $(OUT)

$(OBJDIR):
	mkdir -p $(OBJDIR)

$(OUT): $(OBJ)
	$(F90) $(FFLAGS) -o $@ $(OBJ) $(LIBS)

# ---- UTILS ----
$(OBJDIR)/math_constants.o: $(UTILS)/math_constants.f90
	$(F90) $(FFLAGS) $(MODULE) -c $< -o $@

$(OBJDIR)/utilities.o: $(UTILS)/utilities.f90
	$(F90) $(FFLAGS) $(MODULE) -c $< -o $@

$(OBJDIR)/file_units.o: $(UTILS)/file_units.f90
	$(F90) $(FFLAGS) $(MODULE) -c $< -o $@

$(OBJDIR)/polpak.o: $(UTILS)/polpak.f90
	$(F90) $(FFLAGS) $(MODULE) -c $< -o $@

$(OBJDIR)/elements.o: $(UTILS)/elements.f90
	$(F90) $(FFLAGS) $(MODULE) -c $< -o $@

# ---- FIELDS ----
$(OBJDIR)/sg_scalar.o: $(FIELDS)/FiniteDiff3D/sg_scalar.f90
	$(F90) $(FFLAGS) $(MODULE) -c $< -o $@

$(OBJDIR)/sg_sparse_vector.o: $(FIELDS)/FiniteDiff3D/sg_sparse_vector.f90
	$(F90) $(FFLAGS) $(MODULE) -c $< -o $@

$(OBJDIR)/sg_vector.o: $(FIELDS)/FiniteDiff3D/sg_vector.f90
	$(F90) $(FFLAGS) $(MODULE) -c $< -o $@

$(OBJDIR)/sg_boundary.o: $(FIELDS)/FiniteDiff3D/sg_boundary.f90
	$(F90) $(FFLAGS) $(MODULE) -c $< -o $@

$(OBJDIR)/sg_diff_oper.o: $(FIELDS)/FiniteDiff3D/sg_diff_oper.f90
	$(F90) $(FFLAGS) $(MODULE) -c $< -o $@

# ---- SENS ----
$(OBJDIR)/DataSpace.o: $(SENS)/DataSpace.f90
	$(F90) $(FFLAGS) $(MODULE) -c $< -o $@

$(OBJDIR)/SensMatrix.o: $(SENS)/SensMatrix.f90
	$(F90) $(FFLAGS) $(MODULE) -c $< -o $@

$(OBJDIR)/SensComp.o: $(SENS)/SensComp.f90
	$(F90) $(FFLAGS) $(MODULE) -c $< -o $@

$(OBJDIR)/SymmetryTest.o: $(SENS)/SymmetryTest.f90
	$(F90) $(FFLAGS) $(MODULE) -c $< -o $@

$(OBJDIR)/DataSens.o: $(SENS)/DataSens.f90
	$(F90) $(FFLAGS) $(MODULE) -c $< -o $@

# ---- INV ----
$(OBJDIR)/INVcore.o: $(INV)/INVcore.f90
	$(F90) $(FFLAGS) $(MODULE) -c $< -o $@

$(OBJDIR)/NLCG.o: $(INV)/NLCG.f90
	$(F90) $(FFLAGS) $(MODULE) -c $< -o $@

$(OBJDIR)/DCG.o: $(INV)/DCG.f90
	$(F90) $(FFLAGS) $(MODULE) -c $< -o $@

# ---- 3D_MT ----
$(OBJDIR)/GridDef.o: $(SRC3D)/GridDef.f90
	$(F90) $(FFLAGS) $(MODULE) -c $< -o $@

$(OBJDIR)/GridCalc.o: $(SRC3D)/GridCalc.f90
	$(F90) $(FFLAGS) $(MODULE) -c $< -o $@

$(OBJDIR)/EMfieldInterp.o: $(SRC3D)/EMfieldInterp.f90
	$(F90) $(FFLAGS) $(MODULE) -c $< -o $@

$(OBJDIR)/ModelSpace.o: $(SRC3D)/modelParam/ModelSpace.f90
	$(F90) $(FFLAGS) $(MODULE) -c $< -o $@

$(OBJDIR)/transmitters.o: $(SRC3D)/DICT/transmitters.f90
	$(F90) $(FFLAGS) $(MODULE) -c $< -o $@

$(OBJDIR)/receivers.o: $(SRC3D)/DICT/receivers.f90
	$(F90) $(FFLAGS) $(MODULE) -c $< -o $@

$(OBJDIR)/dataTypes.o: $(SRC3D)/DICT/dataTypes.f90
	$(F90) $(FFLAGS) $(MODULE) -c $< -o $@

$(OBJDIR)/SolnSpace.o: $(SRC3D)/SolnSpace.f90
	$(F90) $(FFLAGS) $(MODULE) -c $< -o $@

$(OBJDIR)/SolverSens.o: $(SRC3D)/SolverSens.f90
	$(F90) $(FFLAGS) $(MODULE) -c $< -o $@

$(OBJDIR)/DataFunc.o: $(SRC3D)/DataFunc.f90
	$(F90) $(FFLAGS) $(MODULE) -c $< -o $@

$(OBJDIR)/DataIO.o: $(SRC3D)/DataIO.f90
	$(F90) $(FFLAGS) $(MODULE) -c $< -o $@

$(OBJDIR)/ForwardSolver.o: $(SRC3D)/ForwardSolver.f90
	$(F90) $(FFLAGS) $(MODULE) -c $< -o $@

$(OBJDIR)/WSfwd2Dmod.o: $(SRC3D)/FWD/Mod2d/WSfwd2Dmod.f90
	$(F90) $(FFLAGS) $(MODULE) -c $< -o $@

$(OBJDIR)/WSfwd1Dmod.o: $(SRC3D)/FWD/Mod2d/WSfwd1Dmod.f90
	$(F90) $(FFLAGS) $(MODULE) -c $< -o $@

$(OBJDIR)/WSfwd2Dpar.o: $(SRC3D)/FWD/Mod2d/WSfwd2Dpar.f90
	$(F90) $(FFLAGS) $(MODULE) -c $< -o $@

$(OBJDIR)/WSutils.o: $(SRC3D)/FWD/Mod2d/WSutils.f90
	$(F90) $(FFLAGS) $(MODULE) -c $< -o $@

$(OBJDIR)/FwdTEmod.o: $(SRC3D)/FWD/Mod2d/FwdTEmod.f90
	$(F90) $(FFLAGS) $(MODULE) -c $< -o $@

$(OBJDIR)/FwdTMmod.o: $(SRC3D)/FWD/Mod2d/FwdTMmod.f90
	$(F90) $(FFLAGS) $(MODULE) -c $< -o $@

$(OBJDIR)/boundary_ws.o: $(SRC3D)/FWD/boundary_ws.f90
	$(F90) $(FFLAGS) $(MODULE) -c $< -o $@

$(OBJDIR)/modelOperator3D.o: $(SRC3D)/FWD/modelOperator3D.f90
	$(F90) $(FFLAGS) $(MODULE) -c $< -o $@

$(OBJDIR)/solver.o: $(SRC3D)/FWD/solver.f90
	$(F90) $(FFLAGS) $(MODULE) -c $< -o $@

$(OBJDIR)/EMsolve3D.o: $(SRC3D)/FWD/EMsolve3D.f90
	$(F90) $(FFLAGS) $(MODULE) -c $< -o $@

$(OBJDIR)/nestedEM.o: $(SRC3D)/FWD/nestedEM.f90
	$(F90) $(FFLAGS) $(MODULE) -c $< -o $@

$(OBJDIR)/ioAscii.o: $(SRC3D)/ioMod/ioAscii.f90
	$(F90) $(FFLAGS) $(MODULE) -c $< -o $@

$(OBJDIR)/Main.o: $(SRC3D)/Main.f90
	$(F90) $(FFLAGS) $(MODULE) -c $< -o $@

# ---- top-level ----
$(OBJDIR)/UserCtrl.o: $(SRCDIR)/UserCtrl.f90
	$(F90) $(FFLAGS) $(MODULE) -c $< -o $@

$(OBJDIR)/Mod3DMT.o: $(SRCDIR)/Mod3DMT.f90
	$(F90) $(FFLAGS) $(MODULE) -c $< -o $@

# ---- MPI target ----
mpi: F90 = mpif90
mpi: FFLAGS += -x f95-cpp-input -DMPI
mpi: $(OBJDIR) $(OUT)

# ---- Intel target ----
intel: F90 = ifort
intel: FFLAGS = -O3 -heap-arrays -parallel
intel: MODULE = -module $(OBJDIR)
intel: $(OBJDIR) $(OUT)

# ---- Intel MPI target ----
intel_mpi: F90 = mpiifort
intel_mpi: FFLAGS = -O3 -heap-arrays -DMPI
intel_mpi: MODULE = -module $(OBJDIR)
intel_mpi: $(OBJDIR) $(OUT)

clean:
	rm -rf $(OBJDIR) $(OUT)

.PHONY: all mpi intel intel_mpi clean
