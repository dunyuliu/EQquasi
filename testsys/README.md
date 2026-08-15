# Test tiers

```
python3 -m pytest testsys/            # default: everything except e2e, ~4 s
python3 -m pytest -m unit           # numerics only
python3 -m pytest -m contract       # formats, compset, build, release gates
python3 -m pytest -m regression     # one guard per defect that actually happened
python3 -m pytest -m e2e testsys/     # builds and runs; ~20 min for BP8, longer with BP5/BP7
```

Tiers are directories; `conftest.py` marks each test by the directory it sits in,
so `-m <tier>` works without decorating anything.

| tier | what it is | cost | needs |
|------|------------|------|-------|
| `unit/` | numerics checkable in-process: the pore-pressure operator, the initial-condition inversion | ms | numpy |
| `contract/` | output formats against sections 4.1-4.3, compset/README consistency, build invocation, CI and release gates | ms | nothing |
| `regression/` | one guard per defect this project actually hit | ms | nothing |
| `e2e/` | builds the solver and runs a benchmark, diffed against a frozen result | minutes | gfortran, MPI, MUMPS |

`e2e` is excluded by default via `addopts` in `pytest.ini`, because it compiles
and runs. Everything else is static and fast enough to run on every save.

## Why the e2e tier matters

The other three tiers read source and compset files; they never execute the
solver. They cannot catch a change that alters the physics while leaving the
text intact, and exactly that happened twice here -- the pore-pressure boundary
treatment and the initial condition each moved the answer by tens of percent
without tripping a single guard.

Oracles:

  - BP5, BP7 -- `test.reference.results/`, via `pytest -m e2e` and `testsys/e2e/test_benchmarks.py`.
    These predate this suite; `testsys/e2e/test_bp5_bp7_regression.py` only wraps
    them so all three benchmarks report through one command.
  - BP8 -- `reference/bp8/summary.json`, frozen from v1.4.7.

Note the difference in strength. The BP5/BP7 references are *regression locks*
on deliberately coarse compset: they detect change, but they were never checked
against another code. The BP8 gold is also a regression lock, but it additionally
agrees with an independent implementation (`taehoKim_ref`) to -2 % on slip at the
injection point, +5 % at 200 m, and 0.7 % on late-time pore pressure. Only the
second kind tells you the physics is right.

## Adding a guard when something breaks

The rule this suite is built on: **every defect gets a test, and the test is
verified to fail on the defect before it is allowed to pass.** A guard written
against already-fixed code proves nothing -- several here were caught doing
exactly that during development and had to be re-checked by reintroducing the
bug.

Put it in `regression/` if it guards one specific past failure, `unit/` if it
checks a numerical property that should hold generally, `contract/` if it pins an
external format or a build/release invariant.

State the failure in the docstring, with numbers. The value of these files is as
much the record of what went wrong as the assertion itself.
