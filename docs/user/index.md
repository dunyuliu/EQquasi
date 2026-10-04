# EQquasi User Guide

EQquasi is a parallel finite element code for quasi-static and
quasi-dynamic earthquake cycle simulation on faults governed by rate- and
state-dependent friction. It is written in Fortran 90, solves with the
direct sparse solver MUMPS (PETSc is an option on conda builds), and uses
Python 3 to set up cases and plot results. EQquasi is the quasi-dynamic half
of the fully dynamic earthquake cycle simulator EQsimu, which pairs it with
the dynamic rupture code EQdyna.

This site covers installing EQquasi, setting up and running a case, the
parameter reference, the output files, the benchmark cases, measured
performance, and what to do when a run stops.

## Where to go

* [Getting started](getting-started.md) -- install on a workstation, Ubuntu, or TACC Lonestar6, and run a first case.
* [Running a case](running-a-case.md) -- the case layout, the earthquake-cycle loop, launcher settings, and plotting.
* [Parameters](parameters.md) -- every default in `user_defined_params.py`, generated from the code.
* [Output files](outputs.md) -- what each cycle writes and what the columns mean.
* [Benchmarks](benchmarks.md) -- the SEAS benchmark cases, the compset list, and how results are checked.
* [Performance](performance.md) -- measured run times and scaling, with the hardware they were measured on.
* [Troubleshooting](troubleshooting.md) -- the solver's stop codes and the setup and build errors you may meet.
* [Citing](citing.md) -- the papers to cite.

## Getting the code

EQquasi is released under the MIT License and hosted on GitHub at
https://github.com/dunyuliu/EQquasi. MUMPS and AZTEC are distributed under
their own licenses.

EQquasi is still under development and comes without any guaranteed
functionality.
