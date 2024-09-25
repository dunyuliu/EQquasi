#!/bin/bash

# This makefile file is part of EQquasi.
# Machines supported include 
#     ls6
#     ubuntu
#     utig (with locally installed mumps through scivision/mumps) 

SYSTEM = ${MACHINE}

ifeq ($(SYSTEM), utig)
    FFLAGS     = -fopenmp -ffree-line-length-none -O3
    FC         = mpif90.openmpi
    
    LIB        = /usr/lib/x86_64-linux-gnu
    INC        = /usr/include
    NETCDF_LIB = -L${LIB} -lnetcdf -lnetcdff
    NETCDF_INC = -I${INC}
    
    MUMPS_ROOT = ${EQQUASIROOT}/mumps/build
    MUMPS_INC  = ${MUMPS_ROOT}/_deps/mumps-src/include
    MUMPS_LIB  = -L${MUMPS_ROOT} -ldmumps -lmumps_common
    PORD_LIB   = -L${MUMPS_ROOT} -lpord
    
    BLAS_LIB   = -L${LIB} -lblas
    SCALAPACK_LIB = -L${LIB} -lscalapack-openmpi 
    LAPACK_LIB = -L${LIB} -llapack

else ifeq ($(SYSTEM), ubuntu)
    FFLAGS     = -fopenmp -ffree-line-length-none -O3 
    FC         = mpif90.openmpi
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
    FC         = mpif90.openmpi
    
    LIB        = /usr/lib/x86_64-linux-gnu
    INC        = /usr/include
    NETCDF_LIB = -L${LIB} -lnetcdf -lnetcdff
    NETCDF_INC = -I${INC}
    
    MUMPS_ROOT = ../mumps/build/local
    MUMPS_INC  = ${MUMPS_ROOT}/include
    MUMPS_LIB  = -L${MUMPS_ROOT}/lib -ldmumps -lmumps_common -lpord
    PORD_LIB   = 
    
    BLAS_LIB   = -L${LIB} -lblas
    SCALAPACK_LIB = -L${LIB} -lscalapack-openmpi 
    LAPACK_LIB = -L${LIB} -llapack
    
else ifeq ($(SYSTEM), ls6)
    FFLAGS     = -qopenmp
    FC         = mpiifort
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
        msr.o       func_lib.o  netcdf_io.o\
        
clean : eqquasi
	rm -r *.o

# order: MUMPS, SCALAPACK, BLAS, LAPACK, NETCDF
eqquasi: $(OBJ) 
	#$(FC) $(OPT1) $(OBJ) -o eqquasi ${BLAS_LIB} ${SCALAPACK_LIB} \
	#${LAPACK_LIB} ${NETCDF_LIB} ${MUMPS_LIB} ${PORD_LIB}
	$(FC) $(OPT1) $(OBJ) -o eqquasi ${MUMPS_LIB} ${PORD_LIB} \
	${SCALAPACK_LIB} ${BLAS_LIB} ${LAPACK_LIB} ${NETCDF_LIB}

eqquasi.o:   eqquasi.f90 solveTimeLoopMUMPS.o      globalvar.o mesh4num.o  \
             meshgen.o   faulting.o  fric.o      \
             library_output.o        createMatrixHolderInCRSFormat.o       elemAssembleInCRS.o   \
             read_input.o            msr.o       func_lib.o  \
             netcdf_io.o
	$(FC) $(OPT) eqquasi.f90 
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
#main_aztec.o: main_aztec.f90 globalvar.o c8g2.o shl.o faulting.o fric.o   library_output.o msr.o elemcal_aztec.o 
#    $(FC) $(OPT) main_aztec.f90 -I${inc_dir_aztec}        
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
#elemcal_aztec.o : elemcal_aztec.f90  globalvar.o
#    $(FC) $(OPT) elemcal_aztec.f90        
msr.o :      msr.f90  globalvar.o
	$(FC) $(OPT) msr.f90        
read_input.o: read_input.f90 globalvar.o
	$(FC) $(OPT) read_input.f90
func_lib.o:  func_lib.f90
	$(FC) $(OPT) func_lib.f90
netcdf_io.o: netcdf_io.f90 
	$(FC) $(OPT) netcdf_io.f90 ${NETCDF_INC}    
