EQquasi
=======

EQquasi is a parallel Finite Element software to simulate quasi-static/quasi-dynamic deformation induced by earthquake fault slips. It is part of the fully dynamic earthquake cycle simulator EQsimu [Liu et al., 2020, GJI](https://doi.org/10.1093/gji/ggz475) to simulate deformation during the inter-seismic, nucleation, and post-seismic phases of an earthquake cycle. It adopts parallel solvers MUMPS and AZTEC to handle the computing loads. Rate- and state- friction law with various forms govern the fault dynamics.

This repository provides the source code, Matlab postprocessing scripts, batch scripts for various HPC systems, and input files different models.

Software requirements
---------------------
* a Unix-like operating system
* build tools cmake
* Fortran compilers
* a functioning MPI environment
* MUMPS/AZTEC installed and paths correctly set up
* Matlab installed for postprocessing

Note
----
EQquasi is still under development and comes without any guaranteed functionality. At the moment, we cannot provide much support for general users.

If you are interested in using EQquasi, please contact Dunyu Liu (dliu@ig.utexas.edu). 
