# Version 1.2.0
*modified:   LICENSE
*renamed:    batch/AutoEQcycRSF.txt -> batch/ADArun.txt
*deleted:    code/Q0/EQQuasi3d.f90
*deleted:    code/Q0/Quasi-CRSelem.f90
*deleted:    code/Q0/Quasi-c8g2.f90
*deleted:    code/Q0/Quasi-faulting.f90
*deleted:    code/Q0/Quasi-globalvar.f90
*deleted:    code/Q0/Quasi-hrglss.f90
*deleted:    code/Q0/Quasi-mesh4num.f90
*deleted:    code/Q0/Quasi-meshgen.f90
*deleted:    code/Q0/fric.f90
*deleted:    code/Q0/makefile
*deleted:    code/Qi/EQQuasi3d.f90
*deleted:    code/Qi/Quasi-CRSelem.f90
*deleted:    code/Qi/Quasi-SHL.f90
*deleted:    code/Qi/Quasi-c8g2.f90
*deleted:    code/Qi/Quasi-faulting.f90
*deleted:    code/Qi/Quasi-globalvar.f90
*deleted:    code/Qi/Quasi-hrglss.f90
*deleted:    code/Qi/Quasi-mesh4num.f90
*deleted:    code/Qi/Quasi-meshgen.f90
*deleted:    code/Qi/makefile
*new file:   code/Read_Input_Files.f90
*new file:   code/c8g2.f90
*new file:   code/crs.f90
*new file:   code/elemcal.f90
*new file:   code/eqquasi.f90
*new file:   code/faulting.f90
*renamed:    code/Qi/fric.f90 -> code/fric.f90
*new file:   code/globalvar.f90
*new file:   code/library_bound.f90
*new file:   code/library_output.f90
*new file:   code/main.f90
*new file:   code/makefile
*new file:   code/mesh4num.f90
*new file:   code/meshgen.f90
*renamed:    code/Q0/Quasi-SHL.f90 -> code/shl.f90
*new file:   input/FE_Stations.txt
*new file:   script/CombineGlobaldat.m
*new file:   script/Convert_slip_stress_depth_prof.m
*new file:   script/Convert_slip_stress_strike_prof.m
*new file:   script/Make_mov_plot_Quasi.m
*new file:   script/Plot_Globaldat.m
*new file:   script/Plot_Single_onfault_st.m
*new file:   script/Plot_rpt.m

#Changes:
New: A license is created.
New: Use one code to simulate all earthquake cycles. For the first cycle, initial conditions 
	are assigned. For later ones, initial conditions are imported from cplot_EQquasi.txt.
New: Significant simplification of the code structure. No variable is transferred between major subroutines. 
New: files are renamed.
New: library_bound.f90 is created. It contains two subroutines bound_load and bound_ft_ku to 
	handle far_field loading and calculate KU on faults for use in subroutine faulting. 
New: library_output.f90 is created. It contains several subroutines output_onfault_st, 
	output_offfault_st, output_onfault_transfer, output_timedy, output_globaldat, and output_prof. 
	Future outputs could be implemented here.
New: Read_Input_Files.f90 is copied from EQdyna 3D version 5.2.0 to import station locations and to 
	read in icstart from currentcycle.txt, which is created by the batch script. Other systematic 
	parameters should be implemented in a later time. 
New: crs.f90 is created to handle the CRS format of stiffness structure construction.
New: elemcal.f90 is created to handle the construction of the stiffness matrix and mass array in 
	the CRS format. 
New: ./script is created for MATLAB scripts for post processing. 
Improvement: on-screen information printing is improved. 
[Verification] SCEC SEAS benchmark BP-5-QD.

# Version 1.1.0
*modified:   REAME.md
*modified:   code/Q0/EQQuasi3d.f90
*modified:   code/Q0/Quasi-CRSelem.f90
*modified:   code/Q0/Quasi-faulting.f90
*modified:   code/Q0/Quasi-globalvar.f90
*modified:   code/Q0/Quasi-meshgen.f90
*modified:   code/Q0/makefile
*modified:   code/Qi/EQQuasi3d.f90
*modified:   code/Qi/Quasi-CRSelem.f90
*modified:   code/Qi/Quasi-faulting.f90
*modified:   code/Qi/Quasi-globalvar.f90
*modified:   code/Qi/Quasi-mesh4num.f90
*modified:   code/Qi/Quasi-meshgen.f90
*modified:   code/Qi/makefile

# Version 1.0.0
*new file:   Changelog.md
*new file:   LICENSE
*new file:   REAME.md
*new file:   batch/AutoEQcycRSF.txt
*new file:   code/Q0/EQQuasi3d.f90
*new file:   code/Q0/Quasi-CRSelem.f90
*new file:   code/Q0/Quasi-SHL.f90
*new file:   code/Q0/Quasi-c8g2.f90
*new file:   code/Q0/Quasi-faulting.f90
*new file:   code/Q0/Quasi-globalvar.f90
*new file:   code/Q0/Quasi-hrglss.f90
*new file:   code/Q0/Quasi-mesh4num.f90
*new file:   code/Q0/Quasi-meshgen.f90
*new file:   code/Q0/fric.f90
*new file:   code/Q0/makefile
*new file:   code/Qi/EQQuasi3d.f90
*new file:   code/Qi/Quasi-CRSelem.f90
*new file:   code/Qi/Quasi-SHL.f90
*new file:   code/Qi/Quasi-c8g2.f90
*new file:   code/Qi/Quasi-faulting.f90
*new file:   code/Qi/Quasi-globalvar.f90
*new file:   code/Qi/Quasi-hrglss.f90
*new file:   code/Qi/Quasi-mesh4num.f90
*new file:   code/Qi/Quasi-meshgen.f90
*new file:   code/Qi/fric.f90
*new file:   code/Qi/makefile
*new file:   library/MUMPS/include/cmumps_c.h
*new file:   library/MUMPS/include/cmumps_root.h
*new file:   library/MUMPS/include/cmumps_struc.h
*new file:   library/MUMPS/include/dmumps_c.h
*new file:   library/MUMPS/include/dmumps_root.h
*new file:   library/MUMPS/include/dmumps_struc.h
*new file:   library/MUMPS/include/mumps_c_types.h
*new file:   library/MUMPS/include/mumps_compat.h
*new file:   library/MUMPS/include/smumps_c.h
*new file:   library/MUMPS/include/smumps_root.h
*new file:   library/MUMPS/include/smumps_struc.h
*new file:   library/MUMPS/include/zmumps_c.h
*new file:   library/MUMPS/include/zmumps_root.h
*new file:   library/MUMPS/include/zmumps_struc.h
*new file:   library/MUMPS/lib/libcmumps.a
*new file:   library/MUMPS/lib/libdmumps.a
*new file:   library/MUMPS/lib/libmumps_common.a
*new file:   library/MUMPS/lib/libpord.a
*new file:   library/MUMPS/lib/libsmumps.a
*new file:   library/MUMPS/lib/libzmumps.a
#Detailed changes:
none