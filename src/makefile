#!/bin/bash

# This makefile file is part of EQquasi.
# Machines supported include 
#     ls6
#     ubuntu
#     utig (with locally installed mumps through scivision/mumps) 
#     conda-linux (everything from an activated conda-forge env; see environment.yml)

SYSTEM = ${MACHINE}

# Fail loudly on an unset MACHINE. Every branch below sets FFLAGS; with no
# branch taken FFLAGS is empty, gfortran falls back to its 132-column limit and
# the build dies deep in globalvar.f90 with "Line truncated" errors that read
# like source corruption rather than a missing variable.
ifeq ($(strip $(SYSTEM)),)
$(error MACHINE is not set. Build with e.g. `EQQUASIROOT=/path/to/repo MACHINE=utig make` -- see README. Supported: utig, ubuntu, local, conda-linux, ls6)
endif

# Some hosts ship a broken mpif90.openmpi wrapper (missing wrapper-data file)
# while plain mpif90 works. Prefer the .openmpi wrapper when it actually runs,
# otherwise fall back. Override with `make FC=...` at any time.
MPIFC := $(shell command -v mpif90.openmpi >/dev/null 2>&1 && mpif90.openmpi -show >/dev/null 2>&1 && echo mpif90.openmpi || echo mpif90)

# Row 8: eqquasi.f90 unconditionally references solveTimeLoopPETSc regardless
# of MACHINE, but only conda-linux compiles and links that object (PETSc via
# conda is not available on utig/ubuntu/local/ls6). Left as a plain Fortran
# CALL, that is an undefined-reference link failure on every other machine.
# PETSC_DEFINE is empty by default and set to -DHAVE_PETSC only in the
# conda-linux branch below; eqquasi.o is preprocessed so the call site
# itself compiles out on machines without it, replaced by a loud runtime stop
# (rule 2) instead of a link error. CPPFLAG is the preprocessor flag itself,
# gfortran-style (-cpp) by default; ls6's mpiifort wants -fpp instead
# (overridden in that branch below) -- not exercised on any host yet
# (rule 9), so this is a plausible-not-verified guess, flagged as such.
PETSC_DEFINE =
CPPFLAG = -cpp

ifeq ($(SYSTEM), utig)
    FFLAGS     = -fopenmp -ffree-line-length-none -O3
    FC         = ${MPIFC}
    
    LIB        = /usr/lib/x86_64-linux-gnu
    INC        = /usr/include
    NETCDF_LIB = -L${LIB} -lnetcdf -lnetcdff
    NETCDF_INC = -I${INC}
    
    MUMPS_ROOT = ${EQQUASIROOT}/mumps/build
    MUMPS_INC  = ${MUMPS_ROOT}/_deps/mumps-src/include
    MUMPS_LIB  = -L${MUMPS_ROOT} -ldmumps -lmumps_common
    PORD_LIB   = -L${MUMPS_ROOT} -lpord

    BLAS_LIB   = -L${LIB} -lblas
    # MUMPS here is built with parallel-root support, so it needs ScaLAPACK and
    # BLACS symbols. Debian's libscalapack-openmpi is not installed on every
    # utig host; AMD AOCL ships a usable one. Prefer the system package, fall
    # back to AOCL with an embedded rpath so the binary runs without callers
    # having to set LD_LIBRARY_PATH. AOCL's static libscalapack.a is not built
    # -fPIC and cannot go into a PIE executable, so link the shared one.
    SCALAPACK_LIB := $(shell \
        if ls ${LIB}/libscalapack-openmpi.so* >/dev/null 2>&1; then \
            echo "-L${LIB} -lscalapack-openmpi"; \
        elif ls /opt/AMD/aocl/aocl-linux-gcc-3.0-6/lib/libscalapack.so* >/dev/null 2>&1; then \
            echo "-L/opt/AMD/aocl/aocl-linux-gcc-3.0-6/lib -lscalapack -Wl,-rpath,/opt/AMD/aocl/aocl-linux-gcc-3.0-6/lib"; \
        else \
            echo "-L${LIB} -lscalapack-openmpi"; \
        fi)
    LAPACK_LIB = -L${LIB} -llapack

else ifeq ($(SYSTEM), ubuntu)
    FFLAGS     = -fopenmp -ffree-line-length-none -O3 
    FC         = ${MPIFC}
    LIB        = /usr/lib/x86_64-linux-gnu
    INC        = /usr/include
    
    NETCDF_LIB = -L${LIB} -lnetcdf -lnetcdff
    NETCDF_INC = -I${INC}
    
    MUMPS_LIB  = -L${LIB} -ldmumps -lmumps_common 
    MUMPS_INC  = ${INC}
    PORD_LIB   = -L${LIB} -lpord
    
    BLAS_LIB   = -L${LIB} -lblas
    SCALAPACK_LIB = -L${LIB} -lscalapack-openmpi 
    LAPACK_LIB = -L${LIB} -llapack

else ifeq ($(SYSTEM), local)
    FFLAGS     = -fopenmp -ffree-line-length-none -O3
    FC         = ${MPIFC}
    
    LIB        = /usr/lib/x86_64-linux-gnu
    INC        = /usr/include
    NETCDF_LIB = -L${LIB} -lnetcdf -lnetcdff
    NETCDF_INC = -I${INC}
    
    MUMPS_ROOT = ../mumps/build/local
    MUMPS_INC  = ${MUMPS_ROOT}/include
    MUMPS_LIB  = -L${MUMPS_ROOT}/lib -ldmumps -lmumps_common -lpord
    PORD_LIB   =

    BLAS_LIB   = -L${LIB} -lblas
    # install.mumps.sh builds ScaLAPACK and BLACS alongside MUMPS, so take them
    # from there. The system -lscalapack-openmpi is not present on hosts where
    # MUMPS had to be built locally in the first place.
    SCALAPACK_LIB = -L${MUMPS_ROOT}/lib -lscalapack -lblacs
    LAPACK_LIB = -L${LIB} -llapack
    
else ifeq ($(SYSTEM), conda-linux)
    # Compilers, MPI, MUMPS, ScaLAPACK, BLAS and netCDF all come from one
    # activated conda-forge env (environment.yml), so the build is the same on
    # every host and needs no admin. Refuse rather than silently mixing in the
    # system MPI, which is the classic ABI-mismatch crash.
    ifeq ($(strip $(CONDA_PREFIX)),)
    $(error MACHINE=conda-linux needs an activated env: conda env create -f environment.yml && conda activate eqquasi-petsc)
    endif
    FFLAGS     = -fopenmp -ffree-line-length-none -O3
    FC         = ${CONDA_PREFIX}/bin/mpif90
    LIB        = ${CONDA_PREFIX}/lib
    INC        = ${CONDA_PREFIX}/include
    # rpath so bin/eqquasi-<version> runs without LD_LIBRARY_PATH.
    RPATH      = -Wl,-rpath,${LIB}
    NETCDF_LIB = -L${LIB} -lnetcdff -lnetcdf ${RPATH}
    NETCDF_INC = -I${INC}
    MUMPS_INC  = ${INC}
    MUMPS_LIB  = -L${LIB} -ldmumps -lmumps_common
    PORD_LIB   = -L${LIB} -lpord
    BLAS_LIB   = -L${LIB} -lblas
    SCALAPACK_LIB = -L${LIB} -lscalapack
    LAPACK_LIB = -L${LIB} -llapack
    # PETSc (row 8, PATHWAY_FORWARD): conda-forge ships petsc.pc under
    # $CONDA_PREFIX/lib/pkgconfig reporting `Cflags: -I${includedir}` and
    # `Libs: -L${libdir} -lpetsc` -- both already covered by this branch's
    # INC/LIB, so only the -lpetsc flag itself is new. Only wired up here:
    # PETSc via conda is not available on the utig/ubuntu/local/ls6 branches.
    PETSC_INC  = -I${INC}
    PETSC_LIB  = -L${LIB} -lpetsc ${RPATH}
    OBJ_PETSC  = solveTimeLoopPETSc.o
    PETSC_DEFINE = -DHAVE_PETSC

else ifeq ($(SYSTEM), ls6)
    FFLAGS     = -qopenmp
    FC         = mpiifort
    CPPFLAG    = -fpp
    LIB        = ${TACC_MKL_LIB}
    MUMPS_INC  = ${TACC_MUMPS_INC}
    NETCDF_INC = -I${TACC_NETCDF_INC}
    NETCDF_LIB = -L${TACC_NETCDF_LIB} -lnetcdf -L${TACC_NETCDF_LIB} -lnetcdff 
    MUMPS_LIB  = -L${TACC_MUMPS_LIB} -ldmumps -lmumps_common -lpord -lptscotch \
                                    -lptscotcherr -lptesmumps -lesmumps -lscotch \
                                    -lscotcherr -lmetis -lparmetis
    BLAS_LIB   = -L${LIB} -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lpthread
    SCALAPACK_LIB = -L${LIB} -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64
    LAPACK_LIB = -L${LIB} -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core
endif 

OPT  = -c ${FFLAGS} 
OPT1 = ${FFLAGS}
OBJ  =  eqquasi.o   globalvar.o mesh4num.o  meshgen.o   \
        solveTimeLoopMUMPS.o      c8g2.o      shl.o       faulting.o  \
        fric.o               library_output.o \
        createMatrixHolderInCRSFormat.o       elemAssembleInCRS.o   read_input.o            \
        func_lib.o  netcdf_io.o porepressure.o ${OBJ_PETSC}\

clean : eqquasi
	rm -r *.o

# order: MUMPS, SCALAPACK, BLAS, LAPACK, NETCDF
eqquasi: $(OBJ) 
	#$(FC) $(OPT1) $(OBJ) -o eqquasi ${BLAS_LIB} ${SCALAPACK_LIB} \
	#${LAPACK_LIB} ${NETCDF_LIB} ${MUMPS_LIB} ${PORD_LIB}
	$(FC) $(OPT1) $(OBJ) -o eqquasi ${MUMPS_LIB} ${PORD_LIB} \
	${SCALAPACK_LIB} ${BLAS_LIB} ${LAPACK_LIB} ${NETCDF_LIB} ${PETSC_LIB}

eqquasi.o:   eqquasi.f90 solveTimeLoopMUMPS.o      globalvar.o mesh4num.o  \
             meshgen.o   faulting.o  fric.o      \
             library_output.o        createMatrixHolderInCRSFormat.o       elemAssembleInCRS.o   \
             read_input.o            func_lib.o  \
             netcdf_io.o
	$(FC) $(OPT) ${CPPFLAG} ${PETSC_DEFINE} eqquasi.f90
globalvar.o: globalvar.f90
	$(FC) $(OPT) globalvar.f90     
mesh4num.o:  mesh4num.f90  globalvar.o
	$(FC) $(OPT) mesh4num.f90
meshgen.o:   meshgen.f90 globalvar.o func_lib.o
	$(FC) $(OPT) meshgen.f90
solveTimeLoopMUMPS.o:      solveTimeLoopMUMPS.f90    globalvar.o c8g2.o      shl.o       \
             faulting.o  fric.o                \
             library_output.o        createMatrixHolderInCRSFormat.o       elemAssembleInCRS.o   \
             netcdf_io.o
	$(FC) $(OPT) solveTimeLoopMUMPS.f90 -I${MUMPS_INC}        
solveTimeLoopPETSc.o:      solveTimeLoopPETSc.f90    globalvar.o c8g2.o      shl.o       \
             faulting.o  fric.o                \
             library_output.o        createMatrixHolderInCRSFormat.o       elemAssembleInCRS.o   \
             netcdf_io.o
	$(FC) $(OPT) ${CPPFLAG} solveTimeLoopPETSc.f90 ${PETSC_INC}
c8g2.o:      c8g2.f90    globalvar.o shl.o
	$(FC) $(OPT) c8g2.f90    
shl.o:       shl.f90     globalvar.o
	$(FC) $(OPT) shl.f90
faulting.o:  faulting.f90 fric.o globalvar.o library_output.o func_lib.o
	$(FC) $(OPT) faulting.f90        
fric.o:      fric.f90 globalvar.o
	$(FC) $(OPT) fric.f90        
library_output.o: library_output.f90  globalvar.o
	$(FC) $(OPT) library_output.f90
createMatrixHolderInCRSFormat.o :      createMatrixHolderInCRSFormat.f90  globalvar.o
	$(FC) $(OPT) createMatrixHolderInCRSFormat.f90
elemAssembleInCRS.o :  elemAssembleInCRS.f90  globalvar.o
	$(FC) $(OPT) elemAssembleInCRS.f90    
read_input.o: read_input.f90 globalvar.o
	$(FC) $(OPT) read_input.f90
func_lib.o:  func_lib.f90
	$(FC) $(OPT) func_lib.f90
netcdf_io.o: netcdf_io.f90
	$(FC) $(OPT) netcdf_io.f90 ${NETCDF_INC}
porepressure.o: porepressure.f90 globalvar.o
	$(FC) $(OPT) porepressure.f90    
