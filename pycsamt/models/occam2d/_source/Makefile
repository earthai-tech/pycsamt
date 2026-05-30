#----------------------------------------------------
# Makefile for Occam 2D MT Version 3.0
#----------------------------------------------------

#----------------------------------------------------
# Default compiler
#----------------------------------------------------
 FC90 = f90
 FCFLAGS= -O2

#----------------------------------------------------
# GFortran Compiler, Debian LINUX
#----------------------------------------------------
# FC90 = gfortran
# FCFLAGS= -O2 -Wall

#----------------------------------------------------
# IBM XLF Compiler, Mac OS X G5
#----------------------------------------------------
#FC90 =/opt/ibmcmp/xlf/8.1/bin/f90
#FCFLAGS= -O3 -qstrict
# FCFLAGS= -O3 -qstrict -qarch=g5  -Wl,-framework,vecLib
# use the line above for vector library with fast Lapack, 
# make sure to first uncomment LAPACK sgesv call in Occam.f90


#----------------------------------------------------
# Intel iFort Compiler, Mac OS X, using vector library with fast Lapack
#----------------------------------------------------
# FC90 = f90
# FCFLAGS= -fast  -Wl,-framework,vecLib -warn uncalled
# use the line above for vector library with fast Lapack, 
# make sure to first uncomment LAPACK sgesv call in Occam.f90


#----------------------------------------------------
# Absoft Compiler, Mac OS X , using vector library with fast Lapack
# Note, we have found a segmentation fault using Absoft and
# have not been able to account for it.  Beware of Absoft.
#
#----------------------------------------------------
# FC90 =/Applications/Absoft/bin/f90
# FCFLAGS= -O3 -cpu:g5 
# LIBS = -L/Applications/Absoft/lib -lU77  
# FCFLAGS= -O3 -cpu:g5 -Wl,-framework,vecLib    
# use above for vector library with fast Lapack, 
# make sure to uncomment LAPACK sgesv call in Occam.f90

#----------------------------------------------------
TARGETS= clean  Occam2D

OBJSOC= OccamModules.f90 MT2DModules.f90 Occam.f90 MT2D.f90
        
all:  $(TARGETS)

clean:
		rm -f *.o *~ core *.mod
		rm -f Occam2D

Occam2D:$(OBJSOC)
		$(FC90) $(FCFLAGS) -o $@ $(OBJSOC)  $(LIBS)

# General compile rules

.SUFFIXES: .f90 .o 

.f90.o:	
	$(FC90) $(FCFLAGS)	-c -o $@ $< 

