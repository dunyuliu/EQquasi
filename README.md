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
python3 `pytest -m e2e`
```
To install *```EQquasi```* without testing, try
```
./install.eqquasi.sh -m ubuntu
```
instead of ```python3 `pytest -m e2e```. In this case, *MUMPS* should have been installed via Ubuntu's apt-get. <br/>

To install *```EQquasi```* on TACC's HPC Lonestar6, try
```
./install.eqquasi.sh -m ls6
```
In this case, *MUMPS* has been installed by TACC admin. <br/>

To install *```EQquasi```* with local installation of *MUMPS* on Ubuntu, try
```
./install.eqquasi.sh -m local
```

### Rebuilding just the Fortran

`cd src && make` on its own does **not** work: the makefile selects all of its
flags from `MACHINE`, so with it unset `FFLAGS` is empty, gfortran applies its
default 132-column limit and the build fails with `Line truncated` errors in
`globalvar.f90` that look like corrupted source. Both variables are required:

```
EQQUASIROOT=/path/to/EQquasi MACHINE=utig make -C src
mv src/eqquasi bin/
```

The makefile now stops with an explicit message when `MACHINE` is unset.

Two host-specific traps on the utig workstation, both handled automatically but
worth knowing when a build breaks:

  - **`dmumps_struc.h` not found.** The `utig` target takes MUMPS headers from
    `${EQQUASIROOT}/mumps/build/_deps/mumps-src/include`. That directory holds
    fetched *source*, and pruning it leaves the built `.a` libraries in place
    while making the tree unbuildable. The header must match the libraries: check
    the version in `mumps/build/CMakeCache.txt` against the `MUMPS x.y.z` string
    in the header, since a mismatched `DMUMPS_STRUC` layout corrupts memory
    silently rather than failing to link.
  - **`cannot find -lscalapack-openmpi`.** This MUMPS is built with parallel-root
    support and needs ScaLAPACK/BLACS, but the Debian package is not installed on
    every host. The makefile falls back to AMD AOCL with an embedded rpath. Link
    the AOCL *shared* library; its `libscalapack.a` is not built `-fPIC` and
    cannot go into a PIE executable.

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

In addition, ```test.*``` compsets are small, fast versions used by the test
suite and deliberately not listed in ```compsets.txt```.

Where things run
---------------------
Every scratch artifact -- generated cases, simulation output, build products --
belongs under ```work/``` at the repository root, which is gitignored. Nothing
scratch is written to the repo root itself.

```
work/          # gitignored; create your cases here
reference/     # committed reference results; never wiped
bin/           # gitignored build product
```

A case run through ```run.sh``` keeps each cycle's output in ```cycle0/```,
```cycle1/```, ... The post-processing utilities discover those automatically:
from inside a case, ```plotRuptureTime.py``` with no arguments processes every
cycle.

Example
---------------------
A good starting example would be compset==bp5.qdc.2000 (benchmark problem 5, quasi-dynamic, 2000 m on-fault resolution). <br/>
The case can be created by the following command:
```
create.newcase caseDir bp5.qdc.2000
```
With the default ```user_defined_params.py``` the first cycle takes about
27 minutes on Lonestar6, or ~51 minutes on 4 MPI ranks of a shared
workstation. <br/>

Key progress
---------------------
*v1.2.1* with *MUMPS* is benchmarked in [SEAS BP5](https://strike.scec.org/cvws/seas/benchmark_descriptions.html) and results are published in [*Jiang et al. (2022, JGR)*](https://doi.org/10.1029/2021JB023519).

Verification
---------------------

### Test suite

Five tiers, selected by marker (`pytest.ini`). The first three read source and
reference files only -- no MPI, no MUMPS, seconds. The `e2e` tiers build the
code and run benchmarks.

```
python3 -m pytest tests/              # unit + contract + regression, ~4 s
python3 -m pytest tests/ -m e2e_fast  # what CI runs on every push, ~20 min
python3 -m pytest tests/ -m e2e       # adds the full BP5 cycle, ~75 min
```

| tier | checks | cost |
|---|---|---|
| `unit` | numerics in isolation: the pore-pressure operator against its analytic solution, initial conditions against the benchmark equations, physical invariants | ms |
| `contract` | file formats and cross-file agreement: BP8 section 4 conformance, `model.txt`'s positional contract, compset registration, the CI gate, and that every reference file has a reader | ~1 s |
| `regression` | one guard per defect that actually occurred here | ~2 s |
| `e2e_fast` | BP5 at 101 steps, BP8 over 30 days, and a clean `install.eqquasi.sh` build | ~20 min |
| `e2e` | the above plus BP5's full first cycle | ~75 min |

### References

`reference/<benchmark>/` holds frozen results. A reference is a **run**, not a
file: `reference/bp5/cycle0/` is a full earthquake cycle,
`reference/bp5/cycle0-step101-fast/` the same case stopped at step 101.

`tests/e2e/cases.py` holds the case table and the single runner; a new
benchmark is a row there plus a reference directory. What gets compared is
decided by what the reference contains -- fault snapshots per fault at
max\|diff\| 0, station and profile series in full, and the scalars each
`summary.json` names. Results are reported per file:

```
SUCCESS fltst_strk000dp000.txt: 4483 rows x 9 cols, max|diff| = 0.0e+00
FAIL    global.dat: max|diff| = 3.1e-04 at column 6
```

References are **regression locks, not validations**. They detect unintended
change; they do not establish that a benchmark is reproduced correctly. See
`reference/bp8/README.md` for what BP8's does and does not establish.

### BP8

BP8-specific findings -- pore-pressure solver convergence, the domain-size
sweep, the time-step study, the initial condition, and how to run and package a
submission -- are in `reference/bp8/README.md` and `reference/bp8/`. They are
kept there because they belong with the reference they describe, and because
they change as the benchmark does.

Note
----
*```EQquasi```* is still under development and comes without any guaranteed functionality.

If you are interested in using *```EQquasi```*, please contact Dunyu Liu (dliu@ig.utexas.edu). 
