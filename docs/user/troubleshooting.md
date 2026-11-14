# Troubleshooting

## When a cycle fails

`run.sh` stops at the first failed cycle. It prints the exit status (and
the signal, for a status above 128), then the last 20 lines of
`log/cycle<N>.log`, and leaves `scratch/` as the solver left it. The full
screen output is in `log/cycle<N>.log`. Fix the cause, then:

```
rm -rf scratch
bash run.sh
```

A status above 128 means the process was killed by a signal, most often by
the system running out of memory or by a job time limit, not by the solver.

## Solver stop codes

When the solver refuses an input or meets a state it cannot continue from,
it prints a boxed message naming the problem and stops with one of these
codes. The launcher normally passes the code through as the exit status.

| code | meaning | what to do |
|---:|---|---|
| 2 | Mesh consistency check failed: node, element, fault-node or equation counts disagree between the counting and building passes. | Report it, with the case's `user_defined_params.py`. |
| 3 | No fault-adjacent elements were tagged, so every fault node would have zero mass. The model is too thin in y for this `dx`. | Reduce `dx`, or widen the model along y (`fymin`, `fymax`). |
| 4 | Either the rate-and-state box (`xminc`, `xmaxc`, `zminc`) falls outside the model, or a fault's y offset is not a whole multiple of `dy`. | Fix the box or the model bounds; or choose `dy` (the message lists valid values) or move the fault. |
| 5 | A fault's x or z bound is not a whole multiple of `dx` from the mesh origin, so its edge falls between node lines. | Fix `faultgeom` or `dx`. |
| 6 | A required input is missing: the solver was not started from a prepared case (`input/` and `scratch/`), `on_fault_vars_input.nc` or the restart file `fault.r.nc` is absent, or a declared fault ended up with zero nodes. | Run through `case.setup` and `run.sh`; for a zero-node fault, check `faultgeom` against `dx` and the domain. |
| 7 | Normal-stress caps are on and the initial normal stress lies outside them. | Widen `min_norm`/`max_norm`, move `init_norm` into range, or set `C_normal_stress_caps = 0`. |
| 8 | `model.txt` has a malformed normal-stress-caps line. | Regenerate the inputs with the current `case.setup`. |
| 9 | BP8 Peaceman well: no fault node lies within half a cell of the injection point. | Centre the fault mesh on the origin and choose a `dx` that divides the domain evenly. |
| 10 | BP8 Peaceman well: the equivalent radius is not larger than `fluid_rwell`. | Increase `dx` or decrease `fluid_rwell`. |
| 11 | `model.txt` has a malformed `enlarging_ratio_xz` line, or the value is below 1. | Set `enlarging_ratio_xz >= 1` and rerun `case.setup`. |
| 12 | BP8: the fault's node count is not a whole number of columns, so the pore-pressure stencil would use the wrong neighbours. | Check that the fault's z extent is a whole multiple of `dx`. |
| 13 | `model.txt` has a malformed `xi_gamma` line, or `xi_gamma > 0` and a fault node has no stable time step: `a*sigma <= 0`, or its cell is larger than h* = gamma*mu*Dc/((b-a)*sigma) (Lapusta et al. 2000, eq. 10). | Refine `dx`, or set `xi_gamma = 0` to use the constant `xi`. |
| 501 | Zero mass at a fault node. The message prints the node's coordinates. | Check the fault geometry and mesh near that point. |
| 502 | A negative time step was computed from a negative trial slip rate: the run has become unstable at that node. | The message prints the location and the slip rate; inspect the fault state there. |
| 508 | The effective normal stress at a node is no longer compressive, so rate-and-state friction has no solution. The message prints the time, location, and stresses. | Common at the tips of step-overs, bends and rough faults. Normal-stress caps (`C_normal_stress_caps = 1`) are the regularization `case.setup` suggests. |

### Other failures from the solver

* Exit status 1 after `MUMPS factorization (JOB=4) failed` or
  `MUMPS solve (JOB=3) failed`: MUMPS reported an error (printed as `INFOG(1)` and
  `INFOG(2)`), for example a singular matrix or running out of memory.
* `STOP Stopped` after a netCDF message: a netCDF read or write failed. The
  message names the file. A message that `on_fault_vars_input.nc` has no
  `nid_fault` dimension means it was written by an old `case.setup`;
  regenerate it.
* `par.solver=2 (PETSc) needs MACHINE=conda-linux` (or the same for
  `par.solver=3`): this binary was built without PETSc. Use `solver = 1`
  (MUMPS) or the conda build.
* `par.solver=3: more blocks than cells in a direction`: too many MPI ranks
  for the mesh. Lower `HPC_ncpu`.
* `par.solver=3: fault N is not one rectangle of nodes on a single interior
  y line, or shares nodes with another fault`: the distributed loop numbers
  each rank's part of the mesh in closed form, which needs every fault to be
  a full rectangle of nodes on one y line inside the model, and no two faults
  to touch. Use `solver = 1` or `2`.
* `par.solver=3: an off-fault station without exactly one owning rank`: a
  node matched to an off-fault station is owned by no rank, or by more than
  one. This is a fault in the block bookkeeping, not in the case: report it
  and use `solver = 1` or `2`.
* `par.solver=3: rank N MPI-IO failed on <file>`: a rank could not open or
  write one of the text files every rank writes its own lines of (`cplot_*`,
  `p1output.txt`, `p2output.txt`). Check that the output directory is
  writable and on a file system all ranks see.
* `par.solver=3: rank N owns M nodes of fault F, not one rectangle` or
  `par.solver=3: rank N cannot tell which rank owns each fault node for the
  pore-pressure exchange`: a rank's part of a fault is not the rectangle of
  nodes the fault snapshot (`fault.*.nc`) and the BP8 pore-pressure exchange
  between ranks assume. This is a fault in the block bookkeeping, not in the
  case: report it and use `solver = 1` or `2`.
* `ON-FAULT STATIONS REQUESTED BUT NOT FOUND`: a warning, not a stop. Some
  requested stations do not lie on a fault node; see [Output
  files](outputs.md).

## Messages from case.setup and create.newcase

* `EQQUASIROOT is not set` -- run `source install.eqquasi.sh` in the repository first.
* `no compset <name>` -- the second argument to `create.newcase` must be a
  directory under `compset/`.
* `no MPIRUN in .../bin/eqquasi-<version>.cfg` -- no binary is installed for
  the current source version. Install with
  `install.eqquasi.sh -m <machine>`.
* `this case has already run` -- `case.setup` will not regenerate the inputs
  of a case that has results. Create a new case, or pass `--force` if the
  inputs are certain to be identical.
* `generateFaultInterface failed` -- the fault-surface generator for
  `insertFaultType > 0` did not finish. Its own error is printed above; a
  missing Python module (such as matplotlib) is a common cause.

`case.setup` also prints warnings that do not stop it: a cohesive zone
narrower than one cell, a maximum time step that can overload the fault in
one step, and a non-planar or multi-fault geometry running without
normal-stress caps. Read them before starting a long run.

## Build problems

* **`Line truncated` errors in `globalvar.f90`.** The build ran with
  `MACHINE` unset, so no compiler flags were set. Build through
  `install.eqquasi.sh -m <machine>`, never a bare `make` in `src/`.
* **`dmumps_struc.h` not found (utig).** The `utig` build reads MUMPS
  headers from `mumps/build/_deps/mumps-src/include`. That folder holds
  fetched source; if it was removed, the built libraries remain but the tree
  cannot be compiled. The header must match the libraries: check the version
  in `mumps/build/CMakeCache.txt` against the `MUMPS x.y.z` string in the
  header.
* **`cannot find -lscalapack-openmpi` (utig).** The makefile falls back to
  AMD AOCL's ScaLAPACK; link its shared library, since its static
  `libscalapack.a` is not built position-independent.
* **`MACHINE=conda-linux needs an activated env`.** Run `conda activate eqquasi-petsc` first.
* **`numpy.dtype size changed`.** Use numpy 1.26 or older with the system netCDF4 package.
