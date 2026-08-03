# News in 2023
* 20230327 *```EQquasi```* works on Ubuntu.

*EQquasi*
=======
*```EQquasi```* is a parallel finite-element software to simulate quasi-static/quasi-dynamic earthquake cycle deformation induced by fault slips governed by rate- and state- friction. It is part of the fully dynamic earthquake cycle simulator *```EQsimu```* [(*Liu et al.*, 2020, *GJI*)](https://doi.org/10.1093/gji/ggz475) to simulate deformation during the inter-seismic, nucleation, post-seismic, and dynamic rupture phases of earthquake cycles. It relies on parallel solvers [*MUMPS*](http://mumps-solver.org) or [*AZTEC*](https://trilinos.github.io/aztecoo.html#aztec-21-foundation-for-aztecoo) to handle the computing loads. It is written in FORTRAN90 with pre-staging and post-processing scripts in Python3.

*MUMPS* is distributed under the [CeCILL-C license](http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html) and for proper ackowledgement, please read the LICENCE for *MUMPS*. The newest version of MUMPS can be downloaded through this [link](http://mumps-tech.com/mumps-2/). <br/>
*AZTEC* now comes with [*Trilinos*](https://github.com/trilinos/Trilinos) in the name of *AZTECOO*, but the current *EQquasi* still uses the standalone *AZTEC2.1*. (To-do-list: need to update its license.)  <br/>

Setup of computing environment
---------------------
*```EQquasi```* relies on the following packages for pre-staging and computing. <br/>
  - FORTRAN compilers and MPI <br/>
  - python3 <br/>
  - netCDF <br/>
  - pip <br/>
    - numpy=1.26.4 (or older, due to a change of dtype size in later versions.)
    - netCDF4 (```scripts/case.setup``` imports it to write on_fault_vars_input.nc)

For post-processing, additional Python packages are needed:
  - xarray
  - imageio
  - pdf2image
      
```ubuntu.env.setup.sh``` is a bash shell script to install necessary packages for Ubuntu system.

Simply, 
```
sudo bash ubuntu.env.setup.sh
```

Installation
---------------------
To install and test *```EQquasi```*  on *Ubuntu*, 
```
git clone https://github.com/dunyuliu/EQquasi.git
cd EQquasi
bash make.scripts.executable.sh
python3 testAll.py
```
To install *```EQquasi```* without testing, try
```
./install.eqquasi.sh -m ubuntu
```
instead of ```python3 testAll.py```. In this case, *MUMPS* should have been installed via Ubuntu's apt-get. <br/>

To install *```EQquasi```* on TACC's HPC Lonestar6, try
```
./install.eqquasi.sh -m ls6
```
In this case, *MUMPS* has been installed by TACC admin. <br/>

To install *```EQquasi```* with local installation of *MUMPS* on Ubuntu, try
```
./install.eqquasi.sh -m local
```

To activate bash environment variables $EQQUASIROOT and add executable scripts to $PATH,
```
source install.eqquasi.sh
```
or manually add the paths to .bashrc:
```
export EQQUASIROOT=$(pwd)
export PATH=$(pwd)/bin:$PATH
export PATH=$(pwd)/scripts:$PATH
```
where $(pwd) is the root path for your *```EQquasi```* installation. <br/>
<br/>

How-to-use
---------------------
Only three steps are required to setup and run a new case <br/>
```
create.newcase directoryForYourCase compset
cd directoryForYourCase
# modify user_defined_params.py
./case.setup
bash run.sh
```
Here, ```compset``` stands for predefined cases with each defiend via a single parameter file ```user_defined_params.py``` under /case_input. <br/>
Currently supported compsets are listed in ```case_input/compsets.txt```:
  - bp5.qdc.2000
  - bp7.qdc.a.10
  - bp1001.fdc.250
  - bp1001.fdc.rough.250
  - bp1001.qdc.rough.250
  - liu2020.fdc.planar
  - liu2020.fdc.rough.250

In addition, ```test.bp5.qdc```, ```test.bp5.qdc.dip90``` and
```test.bp7.qdc``` are small, fast compsets used by the regression suite; they
are deliberately not listed in ```compsets.txt```.

Where things run
---------------------
Every scratch artifact -- generated cases, simulation output, build products --
belongs under ```work/``` at the repository root, which is gitignored. Nothing
scratch is written to the repo root itself.

```
work/                    # gitignored; create your cases here
work/test/               # created and wiped by testAll.py
test.reference.results/  # committed oracles; NOT scratch, never wiped
bin/                     # gitignored build product
```

```testAll.py``` begins by deleting its scratch directory, so keeping that
directory under ```work/``` is what stops a stray path from removing tracked
files, and it keeps ```git status``` clean after a run.

Example
---------------------
A good starting example would be compset==bp5.qdc.2000 (benchmark problem 5, quasi-dynamic, 2000 m on-fault resolution). <br/>
The case can be created by the following command:
```
create.newcase caseDir bp5.qdc.2000
```
With the default user_defined_params.py, it should take about 27 minutes to finish the 1st earthquake cycle. For 1st + 2nd earthquake cycles, 1 hour 33 minutes will be expected on Lonestar6. <br/>

Key progress
---------------------
*v1.2.1* with *MUMPS* is benchmarked in [SEAS BP5](https://strike.scec.org/cvws/seas/benchmark_descriptions.html) and results are published in [*Jiang et al. (2022, JGR)*](https://doi.org/10.1029/2021JB023519).

Verification
---------------------

### Regression suite

```
python3 -m pytest tests/   # static/structural guards, ~1 s, no MPI or MUMPS
python3 testAll.py         # physics regression: builds, runs, compares to references
```

`testAll.py` runs `test.bp5.qdc`, `test.bp5.qdc.dip90` and `test.bp7.qdc` and
compares `fault.00101.nc` against `test.reference.results/`. Both BP5 cases
reproduce their references with **max absolute difference 0.0** across all 17
variables.

The committed references are *regression locks*, not physics validations. The
test compsets are deliberately coarse (`test.bp5.qdc` uses `dx = 4000 m` against
BP5's 2000 m; `test.bp7.qdc` uses `dx = 50 m`, which spans BP7's
velocity-weakening disc with only about four cells). They catch unintended
changes; they do not establish that a benchmark is reproduced correctly. That
requires the production compset at its native resolution.

### BP8 pore fluid diffusion

The along-fault pore pressure solver used by `bp == 8` is verified against the
closed-form solution in eq. (22) of the SEAS BP8 benchmark description, which
gives `p(x2, x3, t)` for a Gaussian injection source in an unbounded domain.

Relative error in pore pressure at the injection point, `BP8-QD-GS`:

| t (s)   | `dx` = 50 m | `dx` = 25 m |
|---------|-------------|-------------|
| 25 000  | 5.11 %      | 1.51 %      |
| 50 000  | 5.75 %      | 1.50 %      |
| 75 000  | 5.74 %      | 1.44 %      |
| 100 000 | 5.60 %      | 1.37 %      |

Halving `dx` reduces the error by a factor of about 3.9, i.e. **second-order
convergence**, as expected for the FTCS Laplacian combined with a discretely
normalized Gaussian source. Extrapolating to the production resolution
`dx = 10 m` gives roughly 0.2 %.

The remaining error is source under-resolution, not a modelling error: the
Gaussian has a characteristic size `L_gauss = 50 m`, so at `dx = 50 m` it is
resolved by a single cell and at `dx = 25 m` by two.

Two caveats on the comparison itself:

  - The analytic solution assumes an unbounded fault, whereas the solver imposes
    zero fluid flux on the boundary of the frictional domain at +/-400 m. The
    two only agree while the diffusion length `sqrt(4*alpha*t)` stays well
    inside that boundary; it reaches 400 m at around `t = 8e5 s`, i.e. roughly
    9 days into the 30-day benchmark.
  - Eq. (20) of the benchmark description writes the Gaussian source as
    `(1/L^2) exp(-r^2 / 2L^2)`, which integrates to `2*pi` rather than 1, and is
    therefore inconsistent with the closed form in eq. (22) by that factor. The
    implementation follows eq. (22)/(25) and normalizes the source weights
    discretely so that they sum to exactly 1.

The friction solve was checked independently of the fluid: with
`tau0 = 14.6 MPa` on `sigma_bar = 25 MPa` and `theta = Dc/V_init`, the
regularized rate-and-state law gives `V = 6.54e-11 m/s`, and the code reports
`6.55e-11 m/s` at the first output step.

Note
----
*```EQquasi```* is still under development and comes without any guaranteed functionality.

If you are interested in using *```EQquasi```*, please contact Dunyu Liu (dliu@ig.utexas.edu). 
