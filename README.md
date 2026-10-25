# EQquasi

EQquasi is a parallel finite-element code for simulating quasi-static and
quasi-dynamic earthquake cycle deformation on faults governed by rate- and
state-dependent friction. It is the quasi-dynamic half of the fully dynamic
earthquake cycle simulator EQsimu ([Liu, Duan and Luo, 2020,
*GJI*](https://doi.org/10.1093/gji/ggz475)), which pairs it with the dynamic
rupture code EQdyna. The solver is Fortran 90; case setup and
post-processing are Python 3. It solves its linear systems with
[MUMPS](http://mumps-solver.org) (CeCILL-C license) or, on conda builds,
PETSc.

EQquasi is released under the [MIT License](LICENSE) and is still under
development, without guaranteed functionality.

Full documentation -- install options, the parameter reference, output
files, benchmarks, performance and troubleshooting -- is published at
https://dunyuliu.github.io/EQquasi/, built from `docs/user/`.

## Install

Requires a Fortran compiler, MPI, netCDF (with Fortran bindings), and
Python 3 with numpy 1.26 (or older) and netCDF4 (`scripts/case.setup`
imports it to write `on_fault_vars_input.nc`).

On any Linux host, with no admin rights, everything (compilers and MPI
included) comes from conda-forge:

```
conda env create -f environment.yml
conda activate eqquasi-petsc
bash install.eqquasi.sh -m conda-linux
```

Other targets for `install.eqquasi.sh -m <machine>`: `ubuntu` (system
MUMPS via `sudo bash ubuntu.env.setup.sh` first), `local` (builds a private
MUMPS), `ls6` (TACC Lonestar6, cluster modules), `utig` (UTIG
workstations). Details and troubleshooting for each are in
[Getting Started](https://dunyuliu.github.io/EQquasi/getting-started/).

Each new shell needs `EQQUASIROOT` and `bin/`/`scripts/` on its `PATH`:
`source install.eqquasi.sh`, or add the equivalent `export` lines to your
shell startup file.

## Run an example

```
create.newcase mycase bp5.qdc.2000
cd mycase
./case.setup
bash run.sh
```

`bp5.qdc.2000` is SEAS benchmark BP5 at 2000 m on-fault resolution; the
default `user_defined_params.py` runs one full earthquake cycle in about
27 minutes on 4 MPI ranks. The full compset list is in
[Benchmarks](https://dunyuliu.github.io/EQquasi/benchmarks/); see
[Running a Case](https://dunyuliu.github.io/EQquasi/running-a-case/) for
what each step does and how to configure, resume or plot a case.

## Output

Each completed cycle's output lands in `result/cycle<N>/` inside the case
directory; `log/cycle<N>.log` holds its screen output. Plot it from inside
the case with `plotRuptureTime.py`, `plotPeakSliprateTime.py`,
`plotOnFaultVars`, `plotAccumulated` and `plotStations.py`. File formats
and columns are documented in [Output Files](https://dunyuliu.github.io/EQquasi/outputs/).

## Citing and help

If you use EQquasi in published work, please cite
[Liu, Duan and Luo (2020, *GJI*)](https://doi.org/10.1093/gji/ggz475) and,
for SEAS BP5 results, [Jiang et al. (2022, *JGR*)](https://doi.org/10.1029/2021JB023519);
see [Citing](https://dunyuliu.github.io/EQquasi/citing/) for the full
reference list and MUMPS acknowledgement.

For questions, bug reports or feature requests, open a GitHub issue.
