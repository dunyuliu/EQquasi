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




