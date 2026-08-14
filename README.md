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
    - netCDF4 (```script/case.setup``` imports it to write on_fault_vars_input.nc)

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
python3 -m pytest -m e2e
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
export PATH=$(pwd)/script:$PATH
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
Here, ```compset``` stands for predefined cases with each defiend via a single parameter file ```user_defined_params.py``` under /compsets. <br/>
Currently supported compsets, with what gates each one, are in
```compsets/README.md```:
  - bp5.qdc.2000
  - bp7.qdc.a.10
  - bp8.qdc.gs.10
  - bp1001.fdc.250
  - bp1001.fdc.rough.250
  - bp1001.qdc.rough.250
  - bp1002.qdc.2500 (two-fault step-over; the only `ntotft > 1` compset)
  - liu2020.qdc.kink.300
  - liu2020.fdc.planar.300
  - liu2020.fdc.rough.250
  - das.qdc.10

In addition, ```test.*``` compsets (```test.bp5.qdc```,
```test.bp5.qdc.dip90```, ```test.bp7.qdc```, ```test.bp8.qdc```,
```test.stepover.qdc```) are small, fast versions used by the test
suite and deliberately outside the register.

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

A case created by ```create.newcase``` holds two files you touch and five
folders:

```
<case>/
    user_defined_params.py   the only file you edit
    case.setup  case.submit  run.sh    what you run
    input/                   everything a cycle consumes; frozen once a cycle has run
    result/                  cycle0/ ... cycleN/, one per completed cycle
    scratch/                 the cycle currently running, and nothing else
    tool/                    machinery the case imports, not for editing
    log/                     setup.log, run.log
```

The solver reads ```input/``` and writes ```scratch/```; ```run.sh``` moves
```scratch/``` into ```result/cycleN``` only once the cycle has succeeded, so a
killed cycle can never be mistaken for a finished one and recovery is
```rm -rf scratch```. The post-processing utilities discover cycles
automatically, under ```result/cycle*``` or, for cases predating v1.13.0,
```cycle*``` or ```Q*```.

Post-processing
---------------------
Five utilities cover a run. From inside a case, each takes cycle directories
as arguments, or none to process every cycle:

| tool | writes | scope |
|---|---|---|
| ```plotRuptureTime.py``` | ```rupture_time.png``` | per cycle |
| ```plotPeakSliprateTime.py``` | ```peak_slip_rate_vs_time.png``` | all cycles, concatenated |
| ```plotOnFaultVars``` | ```fault.NNNNN.nc[.fN].png``` + gif | per cycle |
| ```plotAccumulated``` | ```accumulatedSlip.{horizontal,vertical}.png``` | all cycles, stacked |
| ```plotStations.py``` | on-/off-fault station time series | per cycle |

```
plotRuptureTime.py result/cycle2
plotPeakSliprateTime.py result/cycle0 result/cycle1 result/cycle2
plotAccumulated --fault 1 --depth-km -10 --ylim 0 18 result/cycle0 result/cycle2
plotStations.py result/cycle2
```

**Multi-fault runs.** ```plotAccumulated``` takes ```--fault N``` and draws one
fault per figure; pass the same ```--ylim``` to both so the two are comparable,
since autoscaling otherwise invents a difference. ```plotOnFaultVars``` loops
the faults itself and tags output ```.f0```, ```.f1```. ```plotStations.py```
finds the ```fltst_ft<N>_*``` files faults 2+ write. ```plotRuptureTime.py```
reads ```cplot_ft<N>_EQquasi.txt```, so it needs output from v1.13.0 or later;
earlier runs wrote fault 1 only and the tool will report "nothing ruptured"
for a cycle whose event was on another fault.
```plotPeakSliprateTime.py``` is global by construction -- ```global.dat```
records one peak slip rate for the whole model, not one per fault.

**Reading accumulated slip.** The profile crosses the velocity-weakening zone,
the velocity-strengthening margins and the non-rate-and-state creep region, so
a plane-wide maximum is not the earthquake: on BP1002 the largest number in
the file is the imposed creep at |x| > 50 km. Restrict to VW before quoting a
coseismic slip. Each cycle's curves start from the sum of all previous cycles,
so a cycle that adds nothing traces the previous total exactly -- that
coincidence is the inherited initial state, not a duplicate.

```plotAccumulated``` takes each cycle's end state from ```fault.r.nc```, so
its totals do not depend on ```--ninterval```; thinning changes how many
intermediate lines are drawn and nothing else.

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
python3 -m pytest testsys/              # unit + contract + regression, ~1 min
python3 -m pytest testsys/ -m e2e_fast  # what CI runs on every push, ~20 min
python3 -m pytest testsys/ -m e2e       # adds the full BP5 cycle, ~75 min
```

| tier | checks | cost |
|---|---|---|
| `unit` | numerics in isolation: the pore-pressure operator against its analytic solution, initial conditions against the benchmark equations, physical invariants | ms |
| `contract` | file formats and cross-file agreement: BP8 section 4 conformance, `model.txt`'s positional contract, compset registration, the CI gate, and that every reference file has a reader | ~1 s |
| `regression` | one guard per defect that actually occurred here | ~2 s |
| `e2e_fast` | BP5 and BP7 at 101 steps, BP8 over 30 days, BP1002's complete step-over event, and a clean `install.eqquasi.sh` build | ~20 min |
| `e2e` | the above plus BP5's full first cycle | ~75 min |

### References

`reference/<benchmark>/` holds frozen results. A reference is a **run**, not a
file: `reference/bp5/cycle0/` is a full earthquake cycle,
`reference/bp5/cycle0-step101-fast/` the same case stopped at step 101.

`testsys/e2e/cases.py` holds the case table and the single runner; a new
benchmark is a row there plus a reference directory. What gets compared is
decided by what the reference contains -- fault snapshots per fault, station
and profile series in full, and the scalars each `summary.json` names.

Comparison is relative, never bit-for-bit. Two runs of the same case, same
binary, same host differ by ~2e-14 from MPI reduction ordering alone, so exact
equality would fail on noise and could not distinguish it from a regression.
Each value is judged against the size of the quantity it belongs to -- its
column, for a table; its vector, for a fault variable -- so a component that
vanishes by symmetry (dip shear on a strike-slip fault, fault opening under
the contact condition) is measured against the component that does not.
Columns the SEAS format stores as log10 are linearised first, since -30 and
-24 are both zero. A 1 mPa change in dip shear still fails. Results are
reported per file:

```
SUCCESS fltst_strk000dp000.txt: 4483 rows x 9 cols, worst diff = 0.0e+00 of the file's scale
FAIL    global.dat: 12 of 8964 entries outside rtol=1e-09; worst at row 40, column 6
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
