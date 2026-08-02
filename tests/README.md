# Static regression suite

Fast, dependency-light guards. No MPI, no MUMPS, no simulation output — the
whole suite runs in about a second.

```
python3 -m pytest tests/ -v
```

This is **not** the physics regression. That one lives in `testAll.py` +
`check.test.py` + `test.reference.results/`, needs a full build and a real
run, and is the release gate named by rule 3 of `PROJECT_RULES.md`. The suite
here exists to catch the cheap, structural failures that the physics gate
cannot see and that were all hit for real while following the README.

## Note on the directory name

It is `tests/`, plural, deliberately. `testAll.py:6` does `rm -rf test`, so a
directory named `test/` at the repo root gets deleted on every run of the
physics suite.

## What each file guards

| file | defect it exists to catch |
|---|---|
| `test_build_scripts.py` | `install.eqquasi.sh`'s `local` branch was dead code (`[ $MACHINE == "local"]`, `[! -f ...]` — both missing a space, so bash never selected it); `mkdir bin` failed on reinstall; a failing `make` fell through silently to `mv`; the `local` makefile target linked a system `-lscalapack-openmpi` that cannot exist on a host that needed a local MUMPS build; `mpif90.openmpi` was hardcoded and is broken on some hosts |
| `test_release_gate.py` | CI's grep pattern (`error\|failed\|fatal\|abort`) did not match `check.test.py`'s literal `FAIL ` output, and `testAll.py` ran under `\|\| true` — **a BP5 regression printed FAIL and CI still went green**; `check.test.py` also silently skipped every reference file that did not exist, and never set a non-zero exit status |
| `test_compsets.py` | `case_input/bp7.qdc.a.10` sets `par.xmin/...` but `case.setup` reads `par.fxmin/...`, so BP7 silently inherited BP5's 120x100x60 km domain and was OOM-killed at `dx = 10 m`; compsets that override `par.dx` leave `par.dy`/`par.dz` at the class-level default; `compsets.txt` drifted from `case_input/`; the README named a compset (`bp5.qd.2000`) and a script (`create_newcase`) that do not exist, and omitted `netCDF4` from the dependency list |
| `test_io_contracts.py` | `model.txt` is positional and unlabelled — `case.setup`'s writes and `read_input.f90`'s reads must stay in lockstep or every later parameter shifts silently; the on-fault netCDF variable names must agree between writer and reader; the `fric()`, `fltsta()` and `globaldat()` magic indices must stay inside their allocations |

## Strict xfails

Six cases in `test_compsets.py` are marked `xfail(strict=True)` against the
compsets known to carry the shadow-attribute defect. They are owned by branch
`bp7-domain-fix`. Because the marks are strict, the moment that fix lands
these turn into **XPASS failures** — at which point delete the marks and the
`KNOWN_BROKEN_DOMAIN` list. A fixed bug does not get a permanent excuse.
