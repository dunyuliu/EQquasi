#
#  This file is for EQquasi integrated with MUMPS 5.1.2, released
#  on Thur March 22 10:27:00 CDT 2018
#

topdir = /scratch/user/dunyuliu/0LIBRARY/MUMPS_5.1.2
libdir = $(topdir)/lib

include $(topdir)/Makefile.inc

LIBMUMPS_COMMON = $(libdir)/libmumps_common$(PLAT)$(LIBEXT)

LIBDMUMPS = $(libdir)/libdmumps$(PLAT)$(LIBEXT) $(LIBMUMPS_COMMON)

F90 = mpiifort
OPT = -c -qopenmp 
OPT1 = -qopenmp
OBJ = EQQuasi3d.o Quasi-globalvar.o Quasi-mesh4num.o Quasi-meshgen.o \
		Quasi-CRSelem.o Quasi-c8g2.o Quasi-SHL.o Quasi-faulting.o fric.o\
		
eqquasi3drun: $(OBJ) 
	$(F90) $(OPT1) $(OBJ) -o eqquasi3drun $(LIBDMUMPS) $(LORDERINGS) $(LIBS) $(LIBBLAS) $(LIBOTHERS)

EQQuasi3d.o: EQQuasi3d.f90 Quasi-CRSelem.o Quasi-globalvar.o Quasi-mesh4num.o Quasi-meshgen.o Quasi-faulting.o fric.o
	$(F90) $(OPT) EQQuasi3d.f90 
Quasi-globalvar.o: Quasi-globalvar.f90
	$(F90) $(OPT) Quasi-globalvar.f90 	
Quasi-mesh4num.o: Quasi-mesh4num.f90 Quasi-globalvar.o
	$(F90) $(OPT) Quasi-mesh4num.f90
Quasi-meshgen.o: Quasi-meshgen.f90 Quasi-globalvar.o
	$(F90) $(OPT) Quasi-meshgen.f90
Quasi-CRSelem.o: Quasi-CRSelem.f90 Quasi-globalvar.o Quasi-c8g2.o Quasi-SHL.o Quasi-faulting.o fric.o
	$(F90) $(OPT) Quasi-CRSelem.f90			
Quasi-c8g2.o: Quasi-c8g2.f90 Quasi-globalvar.o Quasi-SHL.o
	$(F90) $(OPT) Quasi-c8g2.f90	
Quasi-SHL.o: Quasi-SHL.f90 Quasi-globalvar.o
	$(F90) $(OPT) Quasi-SHL.f90
Quasi-faulting.o: Quasi-faulting.f90 fric.o Quasi-globalvar.o
	$(F90) $(OPT) Quasi-faulting.f90		
fric.o: fric.f90 Quasi-globalvar.o
	$(F90) $(OPT) fric.f90		
clear: eqquasi3drun
	rm *.o
