# Getting Started

## Requirements

* A Fortran compiler and MPI (the builds below use gfortran and Open MPI, or Intel MPI on Lonestar6)
* MUMPS, the parallel direct solver
* netCDF with its Fortran bindings
* Python 3 with numpy 1.26 (or older), netCDF4, xarray and matplotlib; imageio for the animated plots

numpy is pinned at 1.26 because later versions change a dtype size that the
netCDF stack on these systems was built against.

## Install

Clone the repository and make the scripts executable:

```
git clone https://github.com/dunyuliu/EQquasi.git
cd EQquasi
bash make.scripts.executable.sh
```

Then build for your machine. `install.eqquasi.sh -m <machine>` compiles the
solver and installs it as `bin/eqquasi-<version>`, next to a file
`bin/eqquasi-<version>.cfg` that records which MPI launcher to use with it.
There is deliberately no plain `bin/eqquasi`: each case's `run.sh` names the
exact version it runs.

### Any Linux host, no admin rights: conda

Everything, compilers and MPI included, comes from conda-forge
(`environment.yml`):

```
conda env create -f environment.yml
conda activate eqquasi-petsc
bash install.eqquasi.sh -m conda-linux
```

This is also the only build that includes the PETSc solver (`par.solver = 2`).

### Ubuntu 22.04

Install the system packages and Python modules (needs root), then build:

```
sudo bash ubuntu.env.setup.sh
bash install.eqquasi.sh -m ubuntu
```

MUMPS comes from Ubuntu's `libmumps-dev` package. To build a private copy of
MUMPS instead, use `bash install.eqquasi.sh -m local`.

### UTIG workstations

```
bash install.eqquasi.sh -m utig
```

This build uses the MUMPS built under the repository's `mumps/` folder. See
[Troubleshooting](troubleshooting.md#build-problems) if it cannot find the
MUMPS header or ScaLAPACK.

### TACC Lonestar6

```
./install.eqquasi.sh -m ls6
```

This loads the cluster's `netcdf` and `mumps` modules and records `ibrun`
as the launcher.

### Every session

Each new shell needs `EQQUASIROOT` set and the `bin/` and `scripts/` folders
on its path. From the repository root:

```
source install.eqquasi.sh
```

or add the equivalent lines to your shell startup file:

```
export EQQUASIROOT=/path/to/EQquasi
export PATH=$EQQUASIROOT/bin:$EQQUASIROOT/scripts:$PATH
```

Do not build with a bare `make` in `src/`. The makefile takes all of its
flags from `MACHINE`, which the installer sets, and the installer also
writes the launcher record that `case.setup` needs.

## Check the install

The fast test tiers read the source and reference files only. They need
`pytest` but no MPI run, and take about one to two minutes:

```
python3 -m pytest tests/
```

To also build and run small benchmark cases end to end (about 20 minutes):

```
python3 -m pytest tests/ -m e2e_fast
```

## Run a first case

`test.bp5.qdc.2000` is SEAS benchmark BP5 at 2000 m resolution, cut to 101
time steps. It runs in about a minute and a half on 2 MPI ranks. From the
repository root:

```
create.newcase scratch/myfirst test.bp5.qdc.2000
cd scratch/myfirst
./case.setup
bash run.sh
```

The cycle's output lands in `result/cycle0/`. Plot it from inside the case:

```
plotRuptureTime.py
plotPeakSliprateTime.py
plotOnFaultVars
plotStations.py
```

A 101-step run stops while the first earthquake is still nucleating. For a
full cycle, create the case from `bp5.qdc.2000` instead: about 27 minutes on
Lonestar6, or about 51 minutes on 4 MPI ranks of a shared workstation. See
[Running a case](running-a-case.md) for what each step does.
