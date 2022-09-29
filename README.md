*EQquasi*
=======
*```EQquasi```* is a parallel finite-element software to simulate quasi-static/quasi-dynamic earthquake cycle deformation induced by fault slips governed by rate- and state- friction. It is part of the fully dynamic earthquake cycle simulator *```EQsimu```* [(*Liu et al.*, 2020, *GJI*)](https://doi.org/10.1093/gji/ggz475) to simulate deformation during the inter-seismic, nucleation, post-seismic, and dynamic rupture phases of earthquake cycles. It relies on parallel solvers [*MUMPS*](http://mumps-solver.org) or [*AZTEC*](https://trilinos.github.io/aztecoo.html#aztec-21-foundation-for-aztecoo) to handle the computing loads. It is written in FORTRAN90 with post-process scripts written in [*MATLAB*](https://www.mathworks.com/products/matlab.html).

This repository hosts the source code, compiling instruction, post-processing scripts, batch scripts for some High Perforamnce Computing systems, and input files of different benchmark models.

*MUMPS* is distributed under the [CeCILL-C license](http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html) and for proper ackowledgement, please read the LICENCE for *MUMPS*. The newest version of MUMPS can be downloaded through this [link](http://mumps-tech.com/mumps-2/). <br/>
*AZTEC* now comes with [*Trilinos*](https://github.com/trilinos/Trilinos) in the name of *AZTECOO*, but the current *EQquasi* still uses the standalone *AZTEC2.1*. (To-do-list: need to update its license.)  <br/>

Dependence
---------------------
The code requires netCDF libraries, python, FORTRAN and MPI compilers. Currently, it is tested to work on Lonestar6 at TACC.

Installation on Lonestar6
---------------------
```
git clone https://github.com/dunyuliu/EQquasi.git
cd EQquasi
source install-eqquasi.sh
```
Quick Start Guide
---------------------
After the installation, you need three steps to run a new case. <br/>
First, create a new case with the utility *create_newcase*. <br/> 
*create_newcase* takes in two parameters: <br/> 
  (a) case_dir - the case directory where you want to run the job, and
  (b) compset  - the predefined case (bp5, bp1001)
```
create_newcase case_dir compset
```
Second, cd into the case_dir, modifiy the *user_defined_params.py* and then run the utility *case.setup*.
```
case.setup
```
Third, run the utility *case.submit*, which will submit the job to the HPC system.
```
case.submit
```

Example
---------------------
A good starting example would be compset==bp5-qd-2000 (benchmark problem 5, quasi-dynamic, 2000 m on-fault resolution). <br/>
The case can be created by the following command:
```
create_newcase case_dir bp5-qd-2000
```
With the default user_defined_params.py, the code should take about 27 minutes to finish the 1st earthquake cycle. For 1st + 2nd earthquake cycles, 1 hour 33 minutes will be expected on Lonestar6. <br/>

Key progress
---------------------
*v1.2.1* with *MUMPS* is benchmarked in [SEAS BP5](https://strike.scec.org/cvws/seas/benchmark_descriptions.html) and results are published in [*Jiang et al. (2022, JGR)*](https://doi.org/10.1029/2021JB023519).

Operating environment
---------------------
* a Unix-like operating system
* build tools cmake
* Fortran compilers
* a functioning MPI environment
* MUMPS/AZTEC installed and paths correctly set up
* Matlab installed for postprocessing

Note
----
*EQquasi* is still under development and comes without any guaranteed functionality.

If you are interested in using *EQquasi* or have any questions or comments, please contact Dunyu Liu (dliu@ig.utexas.edu). 
