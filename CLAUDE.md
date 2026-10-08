# CLAUDE.md

What an agent must know before touching EQquasi. Rules live in
`PROJECT_RULES.md` (read its card first, every time); open work lives in
`PATHWAY_FORWARD.md`; user documentation lives in `docs/user/`. This file
does not repeat them.

## What this is

EQquasi is a 3D quasi-dynamic finite-element code for earthquake cycles on
rate-and-state faults: Fortran 90 solver (`src/`), Python case tooling and
post-processing (`script/`), case templates (`compset/`), read-only gold
results (`reference/`), and the test gate (`testsys/`). SEAS benchmarks it
targets include BP5, BP7, BP8 (GS and PW) and the BP1002 step-over.

## Build

```
bash install.eqquasi.sh -m <machine>     # conda-linux | utig | ubuntu | local | ls6
```

The binary is `bin/eqquasi-<version>`, never a bare `eqquasi`; the version is
`EQQUASI_VERSION` in `src/globalvar.f90`. `bin/eqquasi-<version>.cfg` records
the MPI launcher (`MPIRUN`) and its flags (`MPIRUN_ARGS`) that every `run.sh`
uses.

- On theo4, build with `-m conda-linux`; `utig` fails there (no `netcdf.mod`).
  `conda` is a zsh function, so from a bash script activate with
  `eval "$(/home/staff/dliu/anaconda_knox/bin/conda shell.bash hook)"; conda activate eqquasi-petsc`.
- Any `src/` change bumps the version in the same PR; never rebuild the
  binary a running case is using.

## Run a case

```
create.newcase <case_dir> <compset>
cd <case_dir>            # edit user_defined_params.py
python3 ./case.setup     # the case's own copy, not script/case.setup
bash run.sh
```

Never call the binary directly. A cycle cannot resume mid-way: a killed run
restarts from step 0.

## Test

```
EQQUASIROOT=<checkout> python3 -m pytest -q -m "not e2e" testsys   # fast tier
EQQUASIROOT=<checkout> python3 -m pytest -q -m e2e_fast testsys    # local sweep, ~20-40 min
```

CI runs the fast tier, the build and two smokes (`-m e2e_ci`). The full
`-m e2e` tier is periodic, not per-PR.

## Traps that cost real time

- `gh` is not on PATH: `/home/staff/dliu/bin/gh`.
- The interactive shell is zsh: `for v in "a b"; set -- $v` does not word-split.
  Put multi-step launches in a bash script.
- Concurrent MPI runs on one host: check `taskset -cp <pid>` and live CPU
  (`top`, not `ps`'s lifetime %) right after launch; a run that goes quiet is
  caught only by watching its `log/cycle*.log` mtime.
- Timing comparisons need the host to themselves; a killed OpenMPI job can
  leave a PRRTE daemon running.
- `git worktree remove` deletes gitignored files: move run output out first.
- Slides for this project live in the separate manuscripts repo, never here.
