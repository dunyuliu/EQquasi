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

| quantity                        | value                              |
|---------------------------------|------------------------------------|
| pore pressure at shutoff        | 13.54 MPa (`sigma_bar` 25 -> 11.5) |
| peak pore pressure              | 13.63 MPa at 4.16 d, i.e. at `t_off` |
| peak slip rate                  | 6.5e-9 m/s at 8.6 d                |
| total slip at 30 d              | 5.58 mm                            |

Effective normal stress is reduced by 54 % at its lowest but never reaches
zero. Note the slip rate peaks about 4.5 days *after* injection stops. That is
the expected signature of a slip front that keeps propagating once injection
ends, with the growing slipping patch loading the centre elastically, but it
has **not** been verified against independent results -- BP8 was released in
July 2026 and no community comparison exists yet.

The friction solve was checked independently of the fluid: with
`tau0 = 14.6 MPa` on `sigma_bar = 25 MPa` and `theta = Dc/V_init`, the
regularized rate-and-state law gives `V = 6.54e-11 m/s`, and the code reports
`6.55e-11 m/s` at the first output step.

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

### Scaling with mesh size

Comparing the two `test.bp8.qdc` rows, which differ only in resolution:

| `dx`   | elements | s / step | peak RSS |
|--------|----------|----------|----------|
| 50 m   | 12 800   | 0.196    | 234 MB   |
| 25 m   | 64 000   | 0.810    | 1.23 GB  |

5x the elements costs about 4.1x the time per step and 5.2x the memory, so over
this range both are close to linear in element count.

Extrapolating to `bp8.qdc.gs.10` at `dx = 10 m` -- roughly 640 000 elements,
10x the `dx = 25 m` case -- suggests order 8 s per step and order 12 GB on 2
ranks, i.e. something like 12 hours for the ~5200 steps of the 30-day
benchmark. Treat that as an order-of-magnitude estimate from two data points,
not a measurement: *```EQquasi```* solves with a direct sparse factorization
(*MUMPS*), whose cost does eventually grow faster than the number of unknowns,
and the trend above may not hold another decade out. Use more MPI ranks, and
more nodes if memory rather than time becomes the binding constraint.

**Measure on an idle machine.** An earlier timing of the `dx = 25 m` case on
this same hardware gave 1549 s instead of 162.8 s -- a 10x error -- purely
because the node was under load average ~33 at the time. Check `uptime` before
trusting any number here.

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
