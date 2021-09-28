#!/bin/bash

# This makefile file is part of software EQquasi.

# This makefile is to compile EQquasi v1.2.1 with MUMPS, a high-performance, 
# parallel direct solver to sparse linear system of equations. It is modified from
# the Makefile.inc that comes with MUMPS5.1.2 package. 

# Created in 2018. 
# Last updated on 09/28/2021. 
# Dunyu Liu (dliu@ig.utexas.edu)
# 
HPC=1 

# For HPC Frontera with MUMPS installed with PETSC
#ml mumps petcs
#topdir = ${TACC_PETSC_DIR}
#libdir = ${TACC_MUMPS_LIB}
#incdir = ${TACC_MUMPS_INC}
#topdirpord = /work2/07931/dunyuliu/0.Dunyu/P25-SEAS-BP5/MUMPS_5.1.2

# Locally installed mumps
topdir = /work2/07931/dunyuliu/0.Dunyu/P25-SEAS-BP5/MUMPS_5.1.2
libdir = ${topdir}/lib
incdir = ${topdir}/include
topdirpord = ${topdir}

MKLROOT=/opt/intel/mkl/lib/intel64

LPORDDIR = $(topdirpord)/PORD/lib/
IPORD    = -I$(topdirpord)/PORD/include/
LPORD    = -L$(LPORDDIR) -lpord

PLAT    =
LIBEXT  = .a
OUTC    = -o 
OUTF    = -o 
RM = /bin/rm -f
CC = mpiicc
FC = mpiifort
FL = mpiifort
AR = ar vr 
#RANLIB = ranlib
RANLIB  = echo
# Make this variable point to the path where the Intel MKL library is
# installed. It is set to the default install directory for Intel MKL.

LAPACK = -L$(MKLROOT) -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core
SCALAP = -L$(MKLROOT) -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64

LIBPAR = $(SCALAP) $(LAPACK)

INCSEQ = -I$(topdir)/libseq
LIBSEQ  = $(LAPACK) -L$(topdir)/libseq -lmpiseq

LIBBLAS = -L$(MKLROOT) -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core 
LIBOTHERS = -lpthread

#Preprocessor defs for calling Fortran from C (-DAdd_ or -DAdd__ or -DUPPER)
CDEFS   = -DAdd_

#Begin Optimized options
OPTF    = -O -nofor_main -DBLR_MT -qopenmp # or -openmp for old compilers
OPTL    = -O -nofor_main -qopenmp
OPTC    = -O -qopenmp
#End Optimized options
 
INCS = $(INCPAR)
LIBS = $(LIBPAR)
LIBSEQNEEDED =

LIBMUMPS_COMMON = ${libdir}/libmumps_common${PLAT}${LIBEXT}
LIBDMUMPS = ${libdir}/libdmumps${PLAT}${LIBEXT} ${LIBMUMPS_COMMON}

OPT = -c -qopenmp 
OPT1 = -qopenmp
OBJ = eqquasi.o globalvar.o mesh4num.o meshgen.o \
		main.o c8g2.o shl.o faulting.o fric.o library_bound.o library_output.o crs.o elemcal.o Read_Input_Files.o\
		
clean : eqquasi
	rm -r *.o
	
F90 = ${FC}
eqquasi: ${OBJ} 
	${F90} ${OPT1} ${OBJ} -o equasi ${LIBDMUMPS} ${LORDERINGS} ${LIBS} ${LIBBLAS} ${LIBOTHERS} ${LPORD}

##${LORDERINGS} ${LIBS} ${LIBBLAS} ${LIBOTHERS}
	
eqquasi.o: eqquasi.f90 main.o globalvar.o mesh4num.o meshgen.o faulting.o fric.o library_bound.o library_output.o crs.o elemcal.o Read_Input_Files.o
	${F90} ${OPT} eqquasi.f90 
globalvar.o:  globalvar.f90
	${F90} ${OPT} globalvar.f90 	
mesh4num.o:  mesh4num.f90  globalvar.o
	${F90} ${OPT} mesh4num.f90
meshgen.o:  meshgen.f90 globalvar.o
	${F90} ${OPT} meshgen.f90
main.o: main.f90 globalvar.o c8g2.o shl.o faulting.o fric.o library_bound.o library_output.o crs.o elemcal.o 
	${F90} ${OPT} main.f90 -I${incdir}			
c8g2.o:  c8g2.f90 globalvar.o shl.o
	${F90} ${OPT} c8g2.f90	
shl.o:  shl.f90 globalvar.o
	${F90} ${OPT} shl.f90
faulting.o:  faulting.f90 fric.o globalvar.o library_output.o
	${F90} ${OPT} faulting.f90		
fric.o: fric.f90 globalvar.o
	${F90} ${OPT} fric.f90		
library_bound.o: library_bound.f90  globalvar.o
	${F90} ${OPT} library_bound.f90
library_output.o: library_output.f90  globalvar.o
	${F90} ${OPT} library_output.f90
crs.o : crs.f90  globalvar.o
	${F90} ${OPT} crs.f90
elemcal.o : elemcal.f90  globalvar.o
	${F90} ${OPT} elemcal.f90	
Read_Input_Files.o: Read_Input_Files.f90 globalvar.o
	${F90} ${OPT} Read_Input_Files.f90
