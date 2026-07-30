# Makefile for Mod2DMT (ModEM 2D MT inversion)
# Adapted from ModEM v6.2.6 distribution (Egbert & Kelbert, 2012)
#
# Usage:
#   make            — serial build with gfortran (default)
#   make intel      — serial build with ifort
#   make mpi        — MPI build with mpif90
#   make clean      — remove object files and binary
#
# Prerequisites:
#   gfortran >= 6.0 (or ifort/mpif90)
#   LAPACK + BLAS (system packages: liblapack-dev libblas-dev)
#
# Output:  ./Mod2DMT

# ---- compiler settings (override on command line) ----
# -cpp is required even for the serial build: Declaration_MPI.f90,
# Sub_MPI.f90, and Main_MPI.f90 are unconditionally `use`d by several
# files below (INVcore.f90, DCG.f90, Mod2DMT.f90, SymmetryTest.f90,
# ModelSpace.f90) but their entire module body is wrapped in
# #ifdef MPI / #endif; without -cpp, gfortran does not strip that and
# fails on the raw preprocessor directives / the unconditional
# `include 'mpif.h'` inside (same issue as, and fixed identically to,
# the 3D Makefile).
F90      ?= gfortran
FFLAGS   ?= -O3 -ffree-line-length-none -cpp
MPIFLAGS ?= -x f95-cpp-input
LIBS     ?= -llapack -lblas

# ---- directories ----
SRCDIR   = .
OBJDIR   = objs
SRC2D    = $(SRCDIR)/2D_MT
UTILS    = $(SRCDIR)/UTILS
INV      = $(SRCDIR)/INV
SENS     = $(SRCDIR)/SENS
LAPACK   = $(SRCDIR)/LAPACK
MPI      = $(SRCDIR)/MPI

MODULE   = -J$(OBJDIR)
OUT      = Mod2DMT

# ---- object list (order follows module dependencies) ----
OBJ = \
  $(OBJDIR)/math_constants.o \
  $(OBJDIR)/utilities.o \
  $(OBJDIR)/file_units.o \
  $(OBJDIR)/declaration_mpi.o \
  $(OBJDIR)/sub_mpi.o \
  $(OBJDIR)/main_mpi.o \
  $(OBJDIR)/DataSpace.o \
  $(OBJDIR)/GridDef.o \
  $(OBJDIR)/EMfield.o \
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
  $(OBJDIR)/WSfwd2D.o \
  $(OBJDIR)/FwdTEmod.o \
  $(OBJDIR)/FwdTMmod.o \
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
  $(OBJDIR)/Mod2DMT.o

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

# ---- MPI stubs (module body is entirely #ifdef MPI / #endif; -cpp
# strips it to an empty module for the serial build, matching how
# INVcore.f90/DCG.f90/Mod2DMT.f90/SymmetryTest.f90/ModelSpace.f90
# unconditionally `use` these regardless of MPI) ----
$(OBJDIR)/declaration_mpi.o: $(MPI)/Declaration_MPI.f90
	$(F90) $(FFLAGS) $(MODULE) -c $< -o $@

$(OBJDIR)/sub_mpi.o: $(SRC2D)/Sub_MPI.f90
	$(F90) $(FFLAGS) $(MODULE) -c $< -o $@

$(OBJDIR)/main_mpi.o: $(MPI)/Main_MPI.f90
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

# ---- 2D_MT ----
$(OBJDIR)/GridDef.o: $(SRC2D)/GridDef.f90
	$(F90) $(FFLAGS) $(MODULE) -c $< -o $@

$(OBJDIR)/EMfield.o: $(SRC2D)/EMfield.f90
	$(F90) $(FFLAGS) $(MODULE) -c $< -o $@

$(OBJDIR)/EMfieldInterp.o: $(SRC2D)/EMfieldInterp.f90
	$(F90) $(FFLAGS) $(MODULE) -c $< -o $@

$(OBJDIR)/ModelSpace.o: $(SRC2D)/modelParam/ModelSpace.f90
	$(F90) $(FFLAGS) $(MODULE) -c $< -o $@

$(OBJDIR)/transmitters.o: $(SRC2D)/DICT/transmitters.f90
	$(F90) $(FFLAGS) $(MODULE) -c $< -o $@

$(OBJDIR)/receivers.o: $(SRC2D)/DICT/receivers.f90
	$(F90) $(FFLAGS) $(MODULE) -c $< -o $@

$(OBJDIR)/dataTypes.o: $(SRC2D)/DICT/dataTypes.f90
	$(F90) $(FFLAGS) $(MODULE) -c $< -o $@

$(OBJDIR)/SolnSpace.o: $(SRC2D)/SolnSpace.f90
	$(F90) $(FFLAGS) $(MODULE) -c $< -o $@

$(OBJDIR)/SolverSens.o: $(SRC2D)/SolverSens.f90
	$(F90) $(FFLAGS) $(MODULE) -c $< -o $@

$(OBJDIR)/DataFunc.o: $(SRC2D)/DataFunc.f90
	$(F90) $(FFLAGS) $(MODULE) -c $< -o $@

$(OBJDIR)/DataIO.o: $(SRC2D)/DataIO.f90
	$(F90) $(FFLAGS) $(MODULE) -c $< -o $@

$(OBJDIR)/ForwardSolver.o: $(SRC2D)/ForwardSolver.f90
	$(F90) $(FFLAGS) $(MODULE) -c $< -o $@

$(OBJDIR)/WSfwd2D.o: $(SRC2D)/FWD/WSfwd2D.f90
	$(F90) $(FFLAGS) $(MODULE) -c $< -o $@

$(OBJDIR)/FwdTEmod.o: $(SRC2D)/FWD/FwdTEmod.f90
	$(F90) $(FFLAGS) $(MODULE) -c $< -o $@

$(OBJDIR)/FwdTMmod.o: $(SRC2D)/FWD/FwdTMmod.f90
	$(F90) $(FFLAGS) $(MODULE) -c $< -o $@

$(OBJDIR)/ioAscii.o: $(SRC2D)/ioMod/ioAscii.f90
	$(F90) $(FFLAGS) $(MODULE) -c $< -o $@

$(OBJDIR)/Main.o: $(SRC2D)/Main.f90
	$(F90) $(FFLAGS) $(MODULE) -c $< -o $@

# ---- top-level ----
$(OBJDIR)/UserCtrl.o: $(SRCDIR)/UserCtrl.f90
	$(F90) $(FFLAGS) $(MODULE) -c $< -o $@

$(OBJDIR)/Mod2DMT.o: $(SRCDIR)/Mod2DMT.f90
	$(F90) $(FFLAGS) $(MODULE) -c $< -o $@

# ---- MPI target ----
mpi: F90 = mpif90
mpi: MPIFLAGS += -DMPI
mpi: FFLAGS += $(MPIFLAGS)
mpi: $(OBJDIR) $(OUT)

# ---- Intel target ----
intel: F90 = ifort
intel: FFLAGS = -O3 -heap-arrays
intel: MODULE = -module $(OBJDIR)
intel: $(OBJDIR) $(OUT)

clean:
	rm -rf $(OBJDIR) $(OUT)

.PHONY: all mpi intel clean
