#
#  This file is for EQquasi integrated with MUMPS 5.1.2, released
#  on Thur March 22 10:27:00 CDT 2018
#

topdir = /work2/07931/dunyuliu/0.Dunyu/P25-SEAS-BP5/MUMPS_5.1.2
libdir = $(topdir)/lib
incdir = /work2/07931/dunyuliu/0.Dunyu/P25-SEAS-BP5/MUMPS_5.1.2/include
include $(topdir)/Makefile.inc

LIBMUMPS_COMMON = $(libdir)/libmumps_common$(PLAT)$(LIBEXT)

LIBDMUMPS = $(libdir)/libdmumps$(PLAT)$(LIBEXT) $(LIBMUMPS_COMMON)

F90 = mpif90
OPT = -c -qopenmp 
OPT1 = -qopenmp
OBJ = eqquasi.o globalvar.o mesh4num.o meshgen.o \
		main.o c8g2.o shl.o faulting.o fric.o library_bound.o library_output.o crs.o elemcal.o Read_Input_Files.o\
		
clean : eqquasi
	rm -r *.o
	
eqquasi: $(OBJ) 
	$(F90) $(OPT1) $(OBJ) -o eqquasi $(LIBDMUMPS) $(LORDERINGS) $(LIBS) $(LIBBLAS) $(LIBOTHERS)
	
eqquasi.o: eqquasi.f90 main.o globalvar.o mesh4num.o meshgen.o faulting.o fric.o library_bound.o library_output.o crs.o elemcal.o Read_Input_Files.o
	$(F90) $(OPT) eqquasi.f90 
globalvar.o:  globalvar.f90
	$(F90) $(OPT) globalvar.f90 	
mesh4num.o:  mesh4num.f90  globalvar.o
	$(F90) $(OPT) mesh4num.f90
meshgen.o:  meshgen.f90 globalvar.o
	$(F90) $(OPT) meshgen.f90
main.o: main.f90 globalvar.o c8g2.o shl.o faulting.o fric.o library_bound.o library_output.o crs.o elemcal.o 
	$(F90) $(OPT) main.f90			
c8g2.o:  c8g2.f90 globalvar.o shl.o
	$(F90) $(OPT) c8g2.f90	
shl.o:  shl.f90 globalvar.o
	$(F90) $(OPT) shl.f90
faulting.o:  faulting.f90 fric.o globalvar.o library_output.o
	$(F90) $(OPT) faulting.f90		
fric.o: fric.f90 globalvar.o
	$(F90) $(OPT) fric.f90		
library_bound.o: library_bound.f90  globalvar.o
	$(F90) $(OPT) library_bound.f90
library_output.o: library_output.f90  globalvar.o
	$(F90) $(OPT) library_output.f90
crs.o : crs.f90  globalvar.o
	$(F90) $(OPT) crs.f90
elemcal.o : elemcal.f90  globalvar.o
	$(F90) $(OPT) elemcal.f90	
Read_Input_Files.o: Read_Input_Files.f90 globalvar.o
	$(F90) $(OPT) Read_Input_Files.f90
