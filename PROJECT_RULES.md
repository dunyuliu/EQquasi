# EQquasi Project Rules

Twenty rules, each earned by an incident in this repository. Rule numbers are
cited from `testsys/`, `script/defaultParameters.py`, `src/read_input.f90` and
commit messages, so they are never reused or renumbered.

## Before you act

Read this card. It is the whole rule book at the moment of action; everything
below it is the record of why.

| | |
|---|---|
| 1 | New file, name, or parameter? Say so and get agreement **first**. |
| 2 | Missing input fails loudly. No default, no skip, no warn-and-proceed. |
| 3 | Before merging `src/*.f90` or `case.setup`: run the fast tiers **and** `-m e2e`. |
| 4 | `model.txt` / `stations.txt` / `fric.txt`: **append** a field, never insert. |
| 5 | New `fric()` slot: constant in globalvar.f90 AND defaultParameters.py + table row. |
| 6 | New benchmark = new `if (bp == N)` block. Never edit an existing one. |
| 7 | New compset: `compsets/<bench>.<mode>.<res>/` + a row in `compsets/README.md`. |
| 8 | `reference/` is read-only and only grows. Every gold file needs a reader. |
| 9 | Build claims are hypotheses until run **on the target host**. `MACHINE` is required. |
| 10 | Docs move with the code, in the same commit. |
| 11 | Comparing two runs? Check both `runInfo.json` versions first. |
| 12 | The gate must keep an `ntotft > 1` case. |
| 13 | No `if (ntotft == 1)` branching. `nint` is shadowed. Watch implicit interfaces. |
| 14 | Work in a worktree. Stage explicit paths — **never** `git add -A`. |
| 15 | Shared 64-core box: check `uptime`, ≤2 runs, wait above load 56. |
| 16 | A subagent's or audit's finding is a hypothesis. Check the source. |
| 17 | CI is the gate. Red CI stops work — do not tag on top of one. |
| 18 | Every compset must create, run, **and** post-process. |
| 19 | Change params → `case.setup` → `run.sh`. Never call the binary directly. |
| 20 | "Plot" = five figures, fixed order. Per fault, or say which fault. |

Two habits that would have caught most violations of the above: read the file
list a commit prints before pushing it, and re-read this card before staging.

## Index

1. [Minimal changes; no new files until necessary](#1-minimal-changes-no-new-files-until-necessary)
2. [No silent fallbacks, swallowed errors, or placeholder data](#2-no-silent-fallbacks-swallowed-errors-or-placeholder-data)
3. [The regression gate is the release gate](#3-the-regression-gate-is-the-release-gate)
4. [model.txt is a positional contract](#4-modeltxt-is-a-positional-contract)
5. [fric()/on_fault_vars magic indices have one authoritative source](#5-friconfaultvars-magic-indices-have-one-authoritative-source)
6. [New bp values are additive branches](#6-new-bp-values-are-additive-branches)
7. [New compsets: naming, placement, and registration](#7-new-compsets-naming-placement-and-registration)
8. [Reference results are read-only, only grow, and every file has a reader](#8-reference-results-are-read-only-only-grow-and-every-file-has-a-reader)
9. [Build/environment claims must be verified on the target host](#9-buildenvironment-claims-must-be-verified-on-the-target-host)
10. [Docs move with the code, in the same change](#10-docs-move-with-the-code-in-the-same-change)
11. [Compare like with like](#11-compare-like-with-like)
12. [The regression gate must include a multi-fault case](#12-the-regression-gate-must-include-a-multi-fault-case)
13. [Known Fortran landmines must not recur](#13-known-fortran-landmines-must-not-recur)
14. [Shared-checkout discipline](#14-shared-checkout-discipline)
15. [Operational safety on the shared compute host](#15-operational-safety-on-the-shared-compute-host)
16. [A finding is a hypothesis until checked against source](#16-a-finding-is-a-hypothesis-until-checked-against-source)
17. [CI is the gate, not the local suite](#17-ci-is-the-gate-not-the-local-suite)
18. [The whole workflow must work for every example](#18-the-whole-workflow-must-work-for-every-example)
19. [Drive runs through the workflow, never the binary directly](#19-drive-runs-through-the-workflow-never-the-binary-directly)
20. [“Plot” means a fixed set of five figures](#20-plot-means-a-fixed-set-of-five-figures)

## How these rules overlap

Rules 3, 12, 17, 18 and 19 are one family — what counts as verified — cited
separately because they were learned separately. 3 is what the gate consists
of; 12 that it must contain a multi-fault case; 17 that CI, not the local
suite, is the verdict; 18 that the whole workflow must run; 19 that a run must
go through it. Rules 4 and 5 are the other pair: both are "a contract with one
authoritative source", 4 for a file format, 5 for magic indices.

---

## 1. Minimal changes; no new files until necessary

Smallest edit that solves the problem. Fold new logic into the file that owns
the concern — a new `bp == 8` block goes in the existing chain in
`src/faulting.f90`, not a new source file. Never refactor unrelated code in the
same change.

The same applies to **names**. A new variable, parameter or constant is a new
thing to keep consistent: reuse one, or inline the value, before inventing a
name. `par` carries a schema (`script/defaultParameters.py`); a key not in
that schema is not a parameter however well it reads.

*Incidents.* A pore-pressure feature that also reorganized `globalvar.f90`
made the diff impossible to bisect against the gate. The kink compset was
written with `par.kinkAngleDeg`/`par.kinkX` hung off `par`, unflagged; they
worked only because the generator read them with `getattr(..., default)` — the
schema never knew them. Violated four times on 2026-08-14 alone.

**Apply**: before adding a file OR a name, say so and get agreement. Never
introduce either as a side effect of doing something else.

---

## 2. No silent fallbacks, swallowed errors, or placeholder data

Missing input, binary, or config fails loudly. No substituted default, no
skipped step, no warn-and-proceed.

**2a. Silent degradation is the same violation, not a lesser one.** A
multi-fault mesh whose fault plane fell between mesh lines got zero nodes and
ran on. A guard that clamps and continues is worse than one that stops: the
run produces numbers.

*Incident.* `install.eqquasi.sh`'s malformed `[ $MACHINE == "local"]` test and
a workflow line that discarded `pytest -m e2e`'s exit code both made a broken
branch look green. Both fixed; `testsys/contract/test_release_gate.py` now
guards the exit-code path.

---

## 3. The regression gate is the release gate

`reference/bp5/`, `reference/bp5.dip90/`, `reference/bp7/`, `reference/bp8/`
and `reference/bp1002/` are the safety net. Before merging any change to
`src/*.f90` or `script/defaultParameters.py`/`case.setup`:

```
python3 -m pytest testsys/          # unit + contract + regression, ~1 min
python3 -m pytest -m e2e testsys/   # builds and runs the benchmarks, ~3 h
```

A gate only counts if it **stops** the commit. Chaining the test run with `;`
instead of `&&` let a red suite through on 2026-08-14.

---

## 4. model.txt is a positional contract

`case.setup`'s `create_model_input_file` writes fixed-order lines;
`read_input.f90`'s `readmodel` reads them back positionally with unlabeled
`read(1002,*)`. `stations.txt` and `fric.txt` follow the same convention.

A new field is **appended after the last existing write/read pair on both
sides**, never inserted — an insertion shifts every field after it and both
sides "succeed" while reading garbage into the wrong variable.

Optional trailing blocks (per-fault `faultgeom`) are read with `iostat`, so
append scalars **before** them: 13 of 16 compsets write no faultgeom, and their
faultgeom read would swallow a scalar placed after it.

**Check**: `testsys/contract/test_io_contracts.py` counts writes against reads.

---

## 5. fric()/on_fault_vars slots are named, and the registry is the authority

`fric(slot, node, fault)` slots carry names now, defined twice with the SAME
names and numbers: `src/globalvar.f90` (Fortran `integer, parameter ::
FR_RSF_A = 9`, ...) and `script/defaultParameters.py` (Python module level),
which also holds the full slot -> meaning -> unit -> writer table. Slots 1-5
stay numeric: their meaning depends on `friclaw`.

- **Adding a slot**: pick an unused number from the table, add the constant to
  BOTH files and a table row, and use the name at every new site.
- **Never repurpose a slot number** -- it silently corrupts every compset
  still writing the old meaning.
- **Compsets may use raw integers** (they are user-facing configuration); the
  table is what a raw integer is checked against. Do not "helpfully" rename
  inside compsets as a side effect of other work (rule 1).
- The netCDF variable-name mapping in `case.setup`'s
  `netcdf_write_on_fault_vars` and `src/netcdf_io.f90` remains the third
  surface that must agree; it now reads in names on both sides.

*History*: until 2026-08-15 every slot was a bare integer at ~370 call
sites; the revamp renamed them in one pass, gated on the full e2e tier
(references compared at their usual tolerances) precisely because a partial
renaming would have created a fourth source of truth.

**Route**: a slot added to one file but not the other -> the compiler or
ImportError catches the named form; a mismatched NUMBER between the two
constant lists is the dangerous one -> `lars-eriksson`.

**Tier**: judgment for additions; the named form makes drift loud where the
bare integers were silent.
---

## 6. New bp values are additive branches

A new benchmark adds a new `if (bp == N) then ... endif` adjacent to the
existing ones — never modify the body of an existing `bp == 5`/`7`/`8` block,
and never change the default path.

This is what makes rule 3's gate meaningful: if new work can only add branches,
a passing regression run is proof the two are independent rather than
coincidence. Rule 13a makes the same argument for `ntotft`.

**Apply**: `git diff` should show 0 modified lines inside existing `bp`
blocks.

---

## 7. New compsets: naming, placement, and registration

- **Production**: `compsets/<benchmark>.<mode>.<res>/`, plus a row in the
  register table in `compsets/README.md`.
- **CI regression**: a smaller `test.<benchmark>.<mode>/`, added to
  `testsys/e2e/cases.py`, with its reference under `reference/<benchmark>/`.

`create.newcase` does **not** validate the name — it does `os.listdir` and
`shutil.copy`. The register is documentation, not an allowlist, so keep it
accurate rather than relying on it to reject typos. It records what gates each
compset, its oracle and where it was published; a compset with no gate is
unverified and must say so.

---

## 8. Reference results are read-only, only grow, and every file has a reader

- Never edit or regenerate a gold file as a side effect of unrelated work.
- A new gold file is committed only after you have personally read the
  comparison verdict. The first reference for a benchmark is unverified by
  definition — say so in the commit rather than implying zero-diff history.
- Re-blessing is **additive**: compare every pre-existing entry first, and
  only add what is new. If an existing file differs, stop and investigate.

**8a. Every reference artifact needs a test that reads it.** Gold nothing
asserts on is dead weight masquerading as a safety net.
`testsys/contract/test_reference_gold_is_referenced.py` enforces this.

---

## 9. Build/environment claims must be verified on the target host

- `src/makefile` hardcodes `mpif90.openmpi` and `-lscalapack-openmpi`, which
  are Debian `update-alternatives` names. On a host without them only bare
  `mpif90` exists. Run `which mpif90.openmpi mpirun.openmpi` first.
- **`MACHINE` is required**; bare `make` fails by design. Do not export a
  default to make a build succeed on a host you have not checked.
- `install.eqquasi.sh` has no `utig` branch — `MACHINE=utig` falls through
  every case, setting no `LD_LIBRARY_PATH`. It works on cotopaxi by accident
  and fails on knox with a missing `libscalapack.so`.
- **Never rebuild while a multi-cycle run is in flight.** `run.sh` invokes the
  versioned binary afresh each cycle, so a mid-run rebuild splits the dataset
  across two builds under one version number. Cost a restart on 2026-08-14.

---

## 10. Docs move with the code, in the same change

When a compset, script, or file is added, moved, or renamed, `README.md`'s
matching section changes in the same commit — and a pre-existing typo one line
above your edit gets fixed while you are there.

---

## 11. Compare like with like

Three ways a comparison here has produced a wrong conclusion:

- **Different binaries.** A BP8 "resolution regression" was chased for hours
  before `runInfo.json` showed v1.4.5 against v1.4.8. Read both runs'
  `runInfo.json` before trusting any diff.
- **More than one variable.** The step-over A/B experiment was read as
  isolating segment length; shortening the segments also moved VW from 36% to
  49–56% and removed the non-RSF ends.
- **Hashes where tolerances belong.** Two runs of the same binary on identical
  inputs are never byte-identical — MPI reduction ordering gives 1e-15 to
  1e-21. Only a numeric comparator can answer "did this change anything".

---

## 12. The regression gate must include a multi-fault case

Every `test.bp*` compset has `ntotft = 1`. Three multi-fault bugs were found by
hand; none by the gate, because nothing in it exercised `ntotft > 1`.
`bp1002.qdc.2500` now carries this, with a reference, an e2e row, and a
physical invariant asserting the seed lands on fault 0 only.

It sits in the **full** tier, not every-push, because the event is 3821 steps
and ~2600 s on 3 ranks. So every-push CI still has no multi-fault case — a
known, costed gap recorded beside the row in `testsys/e2e/cases.py`.

---

## 13. Known Fortran landmines must not recur

Guarded by `testsys/regression/test_code_convention_landmines.py`:

- **13a. No `if (ntotft == 1) ... else ...` branching.** The fault-node engine
  was made `ntotft`-neutral so that `ntotft = 1` exercises the identical path
  `ntotft > 1` does. A branch here leaves the multi-fault path unexercised —
  the same failure rule 12 guards from the gate side.
- **13b. `nint` is not the intrinsic.** `globalvar.f90` declares
  `integer, parameter :: nint = 8`, shadowing it. Use `int(x + 0.5d0)`.
- **13c. Implicit interfaces bite character arguments.** Passing
  `trim(dir)//"name.nc"` to an external subroutine whose dummy is
  `character(len=50)` hands over a temporary of the expression's length; the
  callee still reads 50 bytes and appends adjacent memory to the path. Assign
  into a fixed-length local first.
- **13d. Hardcoded physics constants.** `faulting.f90` carried normal-stress
  caps as literals that silently overrode any compset reaching that branch.
  Physics thresholds belong in `model.txt`, not in source.

---

## 14. Shared-checkout discipline

- **Work in your own `git worktree`**, never `git checkout` a branch in the
  shared checkout. Branching in place has put two commits on another agent's
  branch. Worktrees are siblings of the main checkout.
- **Stage explicit paths, never `git add -A` / `git add .`.** `git status`
  first, then `git add <path> <path>`.
- **Read the file list the commit prints.** A path you did not intend to touch
  is the signal.
- **Machine-local state gets a `.gitignore` entry the first time it appears.**
  The `-A` rule is judgment and will be broken again; an ignore entry cannot
  be.

*Incident.* On 2026-08-03 `git add -A` swept `.claude/scheduled_tasks.lock` —
a session id, pid and process start time — into a commit about the README, and
from there into a public repo with forks, where a history rewrite cannot
retract it.

---

## 15. Operational safety on the shared compute host

Work runs on shared 64-core boxes (`cotopaxi`, `knox`). All simulation
artifacts belong under `work/` (gitignored). Check `uptime` before launching;
at most two of your runs at once at 3 ranks; **wait if load is above ~56**
rather than adding to it. A long run is not more urgent than someone else's
interactive session.

Three shell footguns, each of which has cost a run:

- **`pkill -f <pattern>` can match the issuing shell** and kill it. Kill by
  PID, found with `readlink /proc/<pid>/cwd`.
- **Background commands inherit a stale cwd.** Always `cd` to an absolute path
  inside the command.
- **Killed runs leave NFS `.nfs*` files.** Kill the ranks, wait, then `rm -rf`.

---

## 16. A finding is a hypothesis until checked against source

An audit reported a "critical" friction-law drift where the code was in fact
correct; an agent's nucleation-length estimate was off by 33×. A finding —
including one in this file, or from the process this file describes — is
evidence to go check, not a conclusion to act on. Verify against the cited
source before it drives a fix or a rule change.

This cuts both ways: a subagent caught two claims in a release brief that were
mine and wrong, both belonging to the previous release.

---

## 17. CI is the gate, not the local suite

A red CI run stops work. Do not commit, tag or release on top of one. Check
after pushing: `gh run list --limit 5`, and `gh run view <id> --log-failed`
when red. `gh` is not on PATH on these hosts — use its full path.

*Incident, 2026-08-12.* CI went red on the `mfault` merge and stayed red
through five tagged releases because the local suite was green and nobody
looked. The merge added a `nid_fault` dimension to every on-fault netCDF
write, so runs produced `(1, nz, nx)` against gold's `(nz, nx)`.

*Incident, 2026-08-14.* CI was red from v1.12.0 onward because the Build step
asserted `test -x bin/eqquasi`, an unversioned name that stopped existing when
binaries became versioned. It died before the e2e tier ran at all.

---

## 18. The whole workflow must work for every example

Create the case, run it, post-process it. All three must work for **every**
compset in `compsets/`, and above all for the runs frozen under `reference/`.
A gold reference nobody can regenerate or plot is a file, not a reference.

- `create.newcase <dir> <compset>` then `./case.setup` produces a runnable case
- the solver runs to its own exit criterion, **not a step cap** — a cycle that
  ends at `nstep` is truncated, and its slip is a lower bound
- every post-processing utility runs and produces a figure, or says clearly
  why not (BP8 is aseismic; an empty rupture-time plot is the right answer)

**Check**: `testsys/contract/test_utilities_run_on_every_reference.py`.

---

## 19. Drive runs through the workflow, never the binary directly

Change `user_defined_params.py`, then `./case.setup`, then `bash run.sh`. Do
not invoke `mpirun ... bin/eqquasi` by hand.

`run.sh` owns the cycle loop: it writes `currentcycle.txt`, fetches the restart
pair from `result/cycleN-1`, runs the solver at the case root, and packages
`scratch/` into `result/cycleN` only on a zero exit. Calling the solver
directly skips all of it.

*Incident, 2026-08-12.* The BP5 full-cycle gold was produced by running the
binary directly, so its output landed flat in the case root and no utility
could read it.

---

## 20. “Plot” means a fixed set of five figures

When the user says **"plot"** with no qualification, produce these five, in
this order, and show them:

| # | tool | scope |
|---|---|---|
| 1 | `plotRuptureTime.py` | per cycle |
| 2 | `plotPeakSliprateTime.py` | all cycles |
| 3 | `plotOnFaultVars` | per cycle |
| 4 | `plotAccumulated` | all cycles |
| 5 | `plotStations.py` | per cycle |

Every tool taking `--fault` runs once per fault and both figures are shown: a
single-fault figure from a two-fault model claims half the model did nothing,
and that must be stated, not implied by omission.

Choosing which figures to show after seeing the numbers is how a result gets
talked into existence. The set is fixed in advance so the same five appear
whether they support the story or not.

**Reading accumulated slip**: the profile crosses VW, VS and the non-RSF creep
region, so a plane-wide maximum is not the earthquake — on BP1002 the largest
number in the file is imposed creep at |x| > 50 km. Restrict to VW before
quoting coseismic slip.

---

*History and the full incident record live in git. `PATHWAY_FORWARD.md` holds
open items and the work queue.*
