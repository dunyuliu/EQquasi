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

### Elastic domain size does not matter here

The whole space is truncated by a finite box, so the obvious worry is that the
box is too small. Measured, at `dx = 50 m`, extending the half-width equally in
x, y and z:

| half-width | elements | peak log10 V     | peak moment rate | slip     |
|------------|----------|------------------|------------------|----------|
| 500 m      |    8 000 | -6.377 @ 0.84 d  | 1.03e9 N/s       | 42.76 mm |
| 1000 m     |   64 000 | -6.377 @ 0.84 d  | 9.32e8 N/s       | 42.33 mm |

Eight times the elements leaves peak slip rate unchanged to three decimals and
slip within 1 %. The fault is locked outside `Omega_f` at +/-400 m anyway, so
the locked collar -- not the mesh boundary -- sets the compliance, and there is
no benefit in paying for a larger box at this resolution.

### The time-step factor xi can be loosened

`xi` is the Lapusta et al. (2009) safety factor in `dtev = xi * Dc / V_max`
(`par.xi`). BP5 uses 0.015; BP7 and BP8 use 0.05. At `dx = 50 m`, 4000 steps
each:

| xi   | wall   | simulated time reached | peak log10 V | peak moment rate | slip     |
|------|--------|------------------------|--------------|------------------|----------|
| 0.05 | 9m22s  |  6.03 d                | -6.377       | 1.029e9          | 42.76 mm |
| 0.1  | 9m22s  | 17.00 d                | -6.376       | 1.032e9          | 42.74 mm |
| 0.2  | 9m27s  | 22.33 d                | -6.374       | 1.032e9          | 42.74 mm |

Same wall clock, 3.7x more simulated time, and the answers agree to 0.003 log
units in peak slip rate, 0.3 % in moment rate and 0.05 % in slip. For BP8 at
this resolution `xi = 0.2` is essentially free, cutting a 30-day run from about
20 000 steps to about 5 500.

Note this was measured at `dx = 50 m`. A finer mesh resolves faster transients,
so re-check before assuming it carries to `dx = 10 m`.

### Submitting to the CRESCENT DET platform

```
checkBP8Submission <result_dir> --zip auto --both-station-names
```

validates a result directory against the benchmark description and writes the
upload zip. Naming convention:

```
<modeler>_<code>-<version>-<resolution>m.zip
e.g. dliu_eqquasi-1-4-3-10m.zip
```

Version and resolution are read from `runInfo.json`, not typed, so the dataset
label on the platform always matches the binary and mesh that produced the data.
Version dots become dashes so everything after the underscore is one
dash-separated token, which is what the platform displays.

All BP8 output files use the `.dat` extension, **including the station time
series** (`fltst_strk+000dp+000.dat`). The platform routes purely on filename
and processes `.dat` for every file type; a station file written as `.txt` is
silently ignored, so the section 4.3 profiles render and the section 4.1 time
series do not, with no error reported anywhere. Section 4.1 of the description
lists station names bare, without any extension, which is what led us to `.txt`
in the first place -- the platform's own `processed_files` listing is the
authority.

BP5 and BP7 keep `.txt`, which is their existing convention.

### Full 30-day run and discrete mass conservation

`test.bp8.qdc` run to the benchmark's full `t_f = 30 days` (5185 steps at
`dx = 50 m`, 17 min on 2 ranks) exits on the time criterion as intended and
gives an independent check that does not rely on the Green's function at all.

After injection stops at `t_off = 100 h`, the zero-flux condition on the
boundary of `Omega_f` traps every injected drop inside the frictional domain,
so the pore pressure must relax towards a uniform value fixed only by the total
volume injected, `V = Q0 * t_off = 1080 m^3`:

```
p_uniform = V / (beta * phi * A * L_fwid)
```

| quantity                                | value       |
|-----------------------------------------|-------------|
| predicted, discrete area (17x17 nodes)  | 1.4948 MPa  |
| **code, centre station at t = 30 d**    | **1.5013 MPa** |
| error                                   | **+0.43 %** |

This exercises the diffusion operator, the source normalization and the
boundary condition together.

**A discretization detail this exposes.** The agreement holds against the
*discrete* area, `(17 * 50)^2 = 722 500 m^2`, not the continuum area of
`Omega_f`, `800^2 = 640 000 m^2`. Pore pressure lives at nodes and every node
carries a full `dx^2` cell, including those on the boundary, so the effective
storage area overshoots by one cell width on each side. At `dx = 50 m` that is
a 13 % area overestimate and therefore an 11 % *underestimate* of the
long-time uniform pressure. The error is first order in `dx` and shrinks to
about 2.5 % at the production `dx = 10 m`, but it is a real bias in late-time,
boundary-influenced results and is not removed by the second-order convergence
reported above, which is measured while the pressure front is still far from
the boundary.

Physical behaviour over the full run, at the centre station:

From `work/bp8.prod10`, dx = 10 m, `xi` = 0.2, literal specification initial
condition (see "The initial condition is over-determined" below):

| quantity                        | value                              |
|---------------------------------|------------------------------------|
| peak pore pressure              | 13.051 MPa at `t_off` = 4.17 d     |
| minimum effective normal stress | 11.949 MPa (52 % reduction)        |
| peak slip rate                  | 3.2e-7 m/s (log10 -6.489) at 0.80 d |
| secondary slip-rate peak        | log10 -6.804 at 4.37 d             |
| total slip at 30 d              | 45.04 mm                           |

Effective normal stress never reaches zero.

An earlier revision of this section reported a peak at 8.6 d and read it as "the
expected signature of a slip front that keeps propagating once injection ends".
That reading was wrong. The dominant peak is at **0.80 d**, well *before*
shutoff, and it is an artefact of the over-determined initial condition: the
patch starts 1.67 MPa overstressed and sheds that stress as an early aseismic
transient. The physical response to shutoff is the *secondary* peak at 4.37 d.
Removing the overstress removes the early peak entirely and leaves a maximum of
about log10 -6.80, which is close to the amplitude of the secondary peak here.

The friction solve was checked independently of the fluid: with
`tau0 = 14.6 MPa` on `sigma_bar = 25 MPa` and `theta = Dc/V_init`, the
regularized rate-and-state law gives `V = 6.54e-11 m/s`, and the code reports
`6.55e-11 m/s` at the first output step.

That check was narrower than it looked: it confirms the solver inverts the
friction law correctly, but not that `6.54e-11 m/s` is the value the benchmark
wants. It is not what eq. (28) asks for, and that leads to the following.

### The initial condition is over-determined

BP8 specifies three quantities that one equation already ties together. Eq. (13)
relates `tau`, `V` and `theta`, so a case may prescribe **two** of them. The
benchmark gives all three:

| source    | quantity                                    |
|-----------|---------------------------------------------|
| Table 1   | `tau_init` = 14.6 MPa, `sigma_bar_0` = 25 MPa |
| eq. (28)  | `V(t=0)` = `V_init` = 1e-12 m/s             |
| eq. (30)  | `theta_0` = `D_RS/V_init` = 5.0e8 s         |

They are mutually inconsistent. With eq. (30)'s `theta_0`, eq. (13) gives
`f = 0.51707`, hence `tau = 12.93 MPa`, not the 14.6 MPa of Table 1. Equivalently,
holding Table 1's `tau_init` and eq. (30)'s `theta_0` forces `V(0) = 6.54e-11 m/s`,
65x the `V_init` that eq. (28) prescribes.

Three readings are possible, each sacrificing exactly one constraint:

| reading | keeps | sacrifices |
|---------|-------|------------|
| **A** literal | Table 1 `tau`, eq. (30) `theta` | eq. (28) `V` |
| **B** equilibrated | eq. (28) `V`, eq. (30) `theta`; `tau_0` -> 12.9277 MPa | Table 1 `tau_init` |
| **C** derived state | eq. (28) `V`, Table 1 `tau`; `theta_0` -> 4.02e11 s | eq. (30) `theta` |

Measured in EQquasi, dx = 50 m, `xi` = 0.2, 4000 steps:

| reading | slip at centre | peak `log10 Vmax` | `V` at first solved step |
|---------|----------------|-------------------|--------------------------|
| A (`work/c3.xi0.2`)     | 42.74 mm at 22.3 d | -6.374 @ 0.85 d | 6.542e-11 m/s |
| B (`work/bp8.equil50`)  | 23.97 mm at 23.2 d | -6.748 @ 1.64 d | 1.000e-12 m/s |
| C (`work/bp8.icfix50`)  | 37.30 mm at 22.5 d | -6.395 @ 2.37 d | 1.000e-12 m/s |

The independent whole-space BEM predicted 24.0 mm for reading B; the FEM gives
23.97 mm, agreeing to 0.1 %. Only reading B starts at `V_init` as eq. (28)
requires -- A begins 65x faster.

The choice changes the answer by **1.8x**, so it is not a detail. Reading A is
what the cases ship, because it takes Table 1 and section 3's eq. (30) at face
value. Reading B satisfies both equations of section 3 and adjusts only a Table 1
number, and it is the one that matches the `taehoKim_ref` comparison (~21 mm,
with the slip-rate maximum near `t_off` rather than at 0.8 d). This has been
raised with the benchmark authors; until they resolve it, `tests/
test_initial_conditions.py` pins reading A and records the size of the
inconsistency so it cannot drift unnoticed.

Note also that **Table 1 omits Poisson's ratio**, which the whole-space elastic
kernel depends on. That is a second specification gap worth raising.

### What the discrepancy is *not*

An independent whole-space spectral boundary-integral solver, written from the
PDF and carrying no boundaries at all, reproduces the FEM under the identical
initial condition: **43.7 mm against 45.0 mm**, peak `log10 Vmax` -6.450 against
-6.489. Domain truncation therefore does **not** explain the gap, despite
section 6 of the description warning that it might. Also excluded, with
measurements: the radiation damping term (`faulting.f90:242` gives
`vs*rho/2 = 4.6244e6 Pa s/m` = `mu/(2 cs)`, applied per eq. 8), slip-magnitude
versus `slip_2` confusion, and spatial resolution (the BEM is converged to 0.5 %
across dx = 50/25/10 m).

Computational performance
---------------------

All numbers below were measured on a single shared workstation, **not** on HPC:

  - 2 x AMD EPYC 7532, 32 cores each (64 cores total), 2.4 GHz, 512 MiB L3
  - 1 TB RAM
  - Ubuntu 22.04.5, kernel 6.8.0
  - gfortran 11.4.0, Open MPI 4.1.1, netCDF 4.8.1, MUMPS 5.8.2 built locally
    via [scivision/mumps](https://github.com/scivision/mumps)
  - built with `-O3 -fopenmp`, run on **2 MPI ranks**

| compset            | `dx`   | elements | steps | solver time | wall  | peak RSS |
|--------------------|--------|----------|-------|-------------|-------|----------|
| `test.bp5.qdc`     | 4000 m | 13 500   | 102   | 15.2 s      | 16.8 s| 217 MB   |
| `test.bp7.qdc`     | 50 m   | 12 800   | 102   | 20.0 s      | 21.7 s| 233 MB   |
| `test.bp8.qdc`     | 50 m   | 12 800   | 201   | 39.5 s      | 41.2 s| 234 MB   |
| `test.bp8.qdc`     | 25 m   | 64 000   | 201   | 162.8 s     | 4m56s | 1.23 GB  |

Peak RSS is the maximum resident set size reported by `/usr/bin/time -v` for
the `mpirun` process tree, so it is roughly the per-rank peak rather than the
total across ranks.

### Parallel scaling

Measured on the 8000-element `test.bp8.qdc` case, 200 steps, `OMP_NUM_THREADS=1`,
otherwise idle machine:

| MPI ranks | wall  | s / step | speedup |
|-----------|-------|----------|---------|
| 1         | 29.4 s| 0.147    | 1.00x   |
| 2         | 25.6 s| 0.128    | 1.15x   |
| 4         | 24.1 s| 0.121    | 1.22x   |
| 8         | 23.8 s| 0.119    | 1.24x   |

Eight times the ranks buys 24 %, i.e. about 15 % parallel efficiency. At this
problem size the direct solve is dominated by work that does not distribute, so
**one rank is the sensible default**; add ranks only when the mesh is large
enough that memory, not time, is the constraint.

**Always set `OMP_NUM_THREADS`.** The build uses `-fopenmp`, so if it is unset
each MPI rank may spawn as many threads as it likes; several ranks on one node
then oversubscribe badly. `runInfo.json` records `omp_threads_per_rank = 0` when
it was unset, precisely so this is visible after the fact.

### Scaling with mesh size

Comparing two `test.bp8.qdc` runs that differ only in resolution:

| `dx`   | elements | s / step | peak RSS |
|--------|----------|----------|----------|
| 50 m   | 12 800   | 0.196    | 234 MB   |
| 25 m   | 64 000   | 0.810    | 1.23 GB  |

5x the elements costs about 4.1x the time per step and 5.2x the memory, so over
this range both are close to linear in element count. Treat any extrapolation
beyond it as an order-of-magnitude estimate only: *```EQquasi```* solves with a
direct sparse factorization (*MUMPS*), whose cost eventually grows faster than
the number of unknowns.

**Measure on an idle machine.** An earlier timing of the `dx = 25 m` case on
this same hardware gave 1549 s instead of 162.8 s -- a 10x error -- purely
because the node was under load average ~33. Check `uptime` first.

### Reading the numbers

  - "solver time" is *```EQquasi```*'s own reported time for the time loop;
    "wall" additionally includes mesh generation, matrix assembly, the initial
    *MUMPS* analysis and factorization, and netCDF output.
  - `test.bp8.qdc` runs 201 steps against the others' 102, so compare the
    per-step cost (0.15, 0.20, 0.20 s respectively), not the totals.
  - BP7 and BP8 cost slightly more per step than BP5 at a comparable element
    count because their meshes are more strongly graded along `y`.

Note
----
*```EQquasi```* is still under development and comes without any guaranteed functionality.

If you are interested in using *```EQquasi```*, please contact Dunyu Liu (dliu@ig.utexas.edu). 
