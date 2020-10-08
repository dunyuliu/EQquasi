/* Copyright 2018-2020, Dunyu Liu, dunyuliu@tamu.edu.
* All rights reserved. This file is part of EQquasi, see the LICENSE attached.*/ 

*Introduction to EQquasi.
EQquasi is a paralell finite element software to simulate quasi-static deformation caused by
earthquake fault slips. In particular, EQquasi is part of fully dynamic earthquake simulator
EQsimu (Liu et al., 2020) to simulate deformation during the inter-seismic, nucleation, and 
post-seismic phases of an earthquake cycle. It adopts parallel solvers MUMPS and AZTEC to 
handle the heavy computing loads. Rate- and state- friction law with various forms are imple-
mented. 

### Author:  Dunyu Liu
### Date:    08/2020-10/08/2020
### Contact: dunyuliu@tamu.edu
## Version 1.2.0; Git tag v1.2.0; Parent 1.1.0.
# Major changes:
* significant simplification of the code structure.
* files are renamed. 
* New files library_bound.f90, library_output.f90, Read_Input_Files.f90, crs.f90, elemcal.f9 
*	are created. 
* Subroutines bound_load, bound_ft_ku, output_onfault_st, output_offfault_st, output_onfault_transfer, 
*	output_timedy, output_globaldat, and output_prof  crs, elemcal, readcurrentcycle, 
*	readstations1 and readstations2 are created. 
* ./script/ is created to hold MATLAB post processing scripts. 
* The license is improved. 
* On-sreen printing is improved. 
* [Verification] It is applied to SCEC SEAS BP5-QD. 
* Detailed changes please refer to the Changelog.md/the google doc.

### Author:  Dunyu Liu
### Date:    11/25/2019-01/06/2020
### Contact: dunyuliu@tamu.edu
## Version 1.1.0; Git tag v1.1.0; Parent 1.0.0.
# Major changes:
* Change fault geometry, frictional parameters for SEAS BP4. 
* Change output formats in accordance to BP4.


### Author:  Dunyu Liu
### Date:    03/20/2019
### Contact: dunyuliu@tamu.edu
## Version 1.0.0; Git tag v1.0.0; Root
# Major changes:
* EQquasi is developed to simulate fully dynamic earthquake cycles on a bent fault.
* MUMPS, a parallel direct solver, is implemented to handle the heavy computing loads.
* The non-co-seismic deformation in the first and later cycles are handled in two 
*	codes due to initial conditions are set in the Q0/ and loaded from input files 
*	in Qi/.
* libs and headfiles of MUMPS are put in ./library/
* License, REAME.md, Chagelog.md are created. 




