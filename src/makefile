#!/bin/bash

# This makefile file is part of software EQquasi.

# This makefile is to compile EQquasi v1.3.0 with MUMPS, a high-performance, 
# for Ubuntu.

# machine is defined in install-eqquasi.shl
# currently, it supports: 
#     ls6
#     ubuntu
SYSTEM = ${MACHINE}

ifeq ($(SYSTEM), ubuntu)
    FFLAGS     = -fopenmp -ffree-line-length-none
    FC         = mpif90
    LIB        = /usr/lib/x86_64-linux-gnu
    INC        = /usr/include
    NETCDF_LIB = -L${LIB} -lnetcdf -L{LIB} -lnetcdff
    NETCDF_INC = -I${INC}
    MUMPS_LIB  = -L${LIB} -lmumps_common -L${LIB} -ldmumps
    PORD_LIB   = -L${LIB} -lpord
    BLAS_LIB   = -L${LIB} -lblas
    SCALAPACK_LIB = -L${LIB} -lscalapack-openmpi 
    LAPACK_LIB = -L${LIB} -llapack

else ifeq ($(SYSTEM), ls6)
    FFLAGS     = -qopenmp
    FC         = mpiifort
    LIB        = ${TACC_MKL_LIB}
    INC        = ${TACC_MUMPS_INC}
    NETCDF_INC = -I${TACC_NETCDF_INC}
    NETCDF_LIB = -L${TACC_NETCDF_LIB} -lnetcdf -L${TACC_NETCDF_LIB} -lnetcdff 
    MUMPS_LIB  = ${TACC_MUMPS_LIB}/libmumps_common.a ${TACC_MUMPS_LIB}/libdmumps.a
    PORD_LIB   = -L${TACC_MUMPS_LIB} -lpord
    BLAS_LIB   = -L${LIB} -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lpthread
    SCALAPACK_LIB = -L${LIB} -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64
    LAPACK_LIB = -L${LIB} -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core
endif 

OPT  = -c ${FFLAGS} 
OPT1 = ${FFLAGS}
OBJ  =  eqquasi.o   globalvar.o mesh4num.o  meshgen.o   \
        main.o      c8g2.o      shl.o       faulting.o  \
        fric.o      library_bound.o         library_output.o \
        crs.o       elemcal.o   read_input.o            \
        msr.o       func_lib.o  netcdf_io.o\
        
clean : eqquasi
	rm -r *.o
    
eqquasi: $(OBJ) 
	$(FC) $(OPT1) $(OBJ) -o eqquasi ${BLAS_LIB} ${SCALAPACK_LIB} \
            ${LAPACK_LIB} ${NETCDF_LIB} ${MUMPS_LIB} ${PORD_LIB}
eqquasi.o:   eqquasi.f90 main.o      globalvar.o mesh4num.o  \
             meshgen.o   faulting.o  fric.o      library_bound.o \
             library_output.o        crs.o       elemcal.o   \
             read_input.o            msr.o       func_lib.o  \
             netcdf_io.o
	$(FC) $(OPT) eqquasi.f90 
globalvar.o: globalvar.f90
	$(FC) $(OPT) globalvar.f90     
mesh4num.o:  mesh4num.f90  globalvar.o
	$(FC) $(OPT) mesh4num.f90
meshgen.o:   meshgen.f90 globalvar.o func_lib.o
	$(FC) $(OPT) meshgen.f90
main.o:      main.f90    globalvar.o c8g2.o      shl.o       \
             faulting.o  fric.o      library_bound.o         \
             library_output.o        crs.o       elemcal.o   \
             netcdf_io.o
	$(FC) $(OPT) main.f90 -I${INC}        
#main_aztec.o: main_aztec.f90 globalvar.o c8g2.o shl.o faulting.o fric.o library_bound.o library_output.o msr.o elemcal_aztec.o 
#    $(FC) $(OPT) main_aztec.f90 -I${inc_dir_aztec}        
c8g2.o:      c8g2.f90    globalvar.o shl.o
	$(FC) $(OPT) c8g2.f90    
shl.o:       shl.f90     globalvar.o
	$(FC) $(OPT) shl.f90
faulting.o:  faulting.f90 fric.o globalvar.o library_output.o func_lib.o
	$(FC) $(OPT) faulting.f90        
fric.o:      fric.f90 globalvar.o
	$(FC) $(OPT) fric.f90        
library_bound.o: library_bound.f90  globalvar.o
	$(FC) $(OPT) library_bound.f90
library_output.o: library_output.f90  globalvar.o
	$(FC) $(OPT) library_output.f90
crs.o :      crs.f90  globalvar.o
	$(FC) $(OPT) crs.f90
elemcal.o :  elemcal.f90  globalvar.o
	$(FC) $(OPT) elemcal.f90    
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
