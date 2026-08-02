# EQquasi Project Rules

Index — read this list first; open a rule only when it is load-bearing.

1. [Minimal changes; no new files until necessary](#1-minimal-changes-no-new-files-until-necessary)
2. [No silent fallbacks, swallowed errors, or placeholder data](#2-no-silent-fallbacks-swallowed-errors-or-placeholder-data)
3. [The BP5 oracle is the release gate — it must actually run and actually fail loud](#3-the-bp5-oracle-is-the-release-gate--it-must-actually-run-and-actually-fail-loud)
4. [model.txt is a positional contract — never insert, always append](#4-modeltxt-is-a-positional-contract--never-insert-always-append)
5. [fric()/on_fault_vars magic indices have one authoritative source per direction](#5-fricon_fault_vars-magic-indices-have-one-authoritative-source-per-direction)
6. [New bp values are additive branches, never edits to existing bp branches](#6-new-bp-values-are-additive-branches-never-edits-to-existing-bp-branches)
7. [New compsets: naming, placement, and registration](#7-new-compsets-naming-placement-and-registration)
8. [Reference results are read-only and only grow, never shrink or overwrite](#8-reference-results-are-read-only-and-only-grow-never-shrink-or-overwrite)
9. [Build/environment claims must be verified on the target host, not assumed](#9-buildenvironment-claims-must-be-verified-on-the-target-host-not-assumed)
10. [Docs move with the code, in the same change](#10-docs-move-with-the-code-in-the-same-change)

---

## 1. Minimal changes; no new files until necessary

Smallest edit that solves the problem. Fold new logic into the file it
belongs with (e.g., a new `bp==8` block goes inside the existing `if (bp ==
N)` chain in `src/faulting.f90`, `src/library_output.f90`,
`src/solveTimeLoopMUMPS.f90` — not a new source file). Never refactor
unrelated code in the same change.

**Rationale**: every other rule constrains work you were asked to do; this
one bounds how much you do at all. A pore-pressure-diffusion feature that
also reorganizes `globalvar.f90` makes the diff impossible to bisect against
the BP5 oracle (rule 3).

**How to apply**: before adding a file, grep for whether an existing file
already owns that concern (`src/faulting.f90` owns fault-node physics,
`src/library_output.f90` owns benchmark-format output, `scripts/lib.py` owns
shared Python helpers).

**Tier**: judgment (2).

---

## 2. No silent fallbacks, swallowed errors, or placeholder data

Missing input, binary, or config fails loudly. No substituted default, no
skipped step. Two existing violations to not repeat:

- `install.eqquasi.sh:67` — `elif [ $MACHINE == "local"]` (missing space
  before `]`) and `install.eqquasi.sh:72` — `if [! -f "$MUMPS_LIB_DIR/$file"
  ]` (missing space after `[`) are both malformed bash tests. The `local`
  branch has never executed correctly; nothing has ever caught this because
  nothing exercises `-m local` in CI.
- `.github/workflows/test.yml` line `python3 testAll.py >> testRunLog.txt ||
  true` — the `|| true` discards `testAll.py`'s exit code. The workflow's
  only actual gate is a keyword grep over the log (rule 3), so a crash that
  doesn't print one of `error|exception|failed|fatal|abort` (case-insensitive)
  passes CI silently.

**How to apply**: a new pore-pressure input file, a new `bp==8` compset
parameter, or a new solver path must `stop` with a descriptive message on
missing/malformed input — following the existing pattern at
`src/read_input.f90:23-34` (`INQUIRE` + explicit `stop`) — never fall back to
a default value or an empty array.

**Tier**: mechanical for the two cited syntax bugs (reproduce with `bash -n
install.eqquasi.sh` or `MACHINE=local bash -x install.eqquasi.sh -m local`);
judgment for new code.

---

## 3. The BP5 oracle is the release gate — it must actually run and actually fail loud

`test.reference.results/test.bp5.qdc/Q0/fault.00101.nc` and
`test.reference.results/test.bp5.qdc.dip90/Q0/fault.00101.nc` are the
project's safety net. Any change to `src/*.f90` or to
`scripts/defaultParameters.py`/`scripts/case.setup` must be followed by:

```
python3 testAll.py
```

and `check.test.py`'s printed verdicts for both `test.bp5.qdc` and
`test.bp5.qdc.dip90` must read `SUCCESS` for every compared file — never
`FAIL`, never skipped.

Two gaps in the current gate mean a human must read the log, not just trust
green CI:

- **`check.test.py`'s pass criterion is `np.allclose(..., rtol=1e-3,
  atol=1e-3)`** (`check.test.py`, `compare_nc_files`), not bit-exact. If the
  intended invariant is exact reproduction, that must be verified by eye
  (`diff` the raw values, or drop the tolerance to something far tighter) —
  the stated tolerance alone does not prove "zero difference on every
  variable."
- **CI's success/failure detection does not look for `check.test.py`'s own
  output.** `check.test.py` prints the literal string `'FAIL '+fn1+' '+fn2`
  on divergence, but `.github/workflows/test.yml`'s grep pattern is
  `'error|Error|ERROR|exception|failed|fatal|abort'` — none of these match
  the substring `FAIL ` (`failed` requires the letters `-ed`, which `FAIL `
  does not have). **A regression that changes `fault.00101.nc` values will
  print `FAIL` to the log and CI will still go green.** This is the single
  highest-priority mechanical gap for BP8 work, since BP8 will exercise new
  I/O paths that could quietly perturb BP5/BP7 shared code
  (`src/faulting.f90`, `src/library_output.f90`, `src/netcdf_io.f90`).
- `check.test.py`'s `fileNameList = ['fault.00001.nc','fault.00101.nc',
  'global.dat', 'tdyna.txt']` is guarded by `if os.path.exists(refPath)`
  (`check.test.py`, main loop) — only `fault.00101.nc` exists under either
  reference dir today, so `fault.00001.nc`, `global.dat`, `tdyna.txt` are
  silently never compared. If BP8 needs new ASCII output files compared,
  they must be added to a reference dir *and* actually produced, or this
  silent-skip will hide their absence too.

**Route**: the grep-pattern fix and the `|| true` are release-gate defects →
`haruto-nakamura`. Missing regression coverage for the tolerance question →
`iris-vermeulen`.

**Tier**: mechanical (run `python3 testAll.py`, inspect `check.test.py`
stdout for literal `FAIL`; reproduce the grep gap with `echo 'FAIL x y' |
grep -qiE 'error|Error|ERROR|exception|failed|fatal|abort'; echo $?` → exits
1, confirming the miss).

---

## 4. model.txt is a positional contract — never insert, always append

`scripts/case.setup`'s `create_model_input_file` (lines 13-38) writes 26
lines in a fixed order; `src/read_input.f90`'s `readmodel` (lines 38-63)
reads them back with 26 unlabeled `read(1002,*)` statements in the same
order, matched positionally, not by name. `stations.txt` (`case.setup`
`create_station_input_file` / `src/read_input.f90` `readstations1`/`2`) and
`fric.txt` (`src/read_input.f90` `readfric`, lines 134-205, whose format
branches on `friclaw`) follow the same convention.

Any new parameter for BP8 (e.g., a pore-pressure diffusivity, an injection
rate, a diffusion-grid spacing) **must be appended after the last existing
`f.write(...)`/`read(1002,*)` pair on both sides**, never inserted in the
middle — an insertion silently shifts every field after it and both sides
will "succeed" (no crash) while reading garbage into the wrong variable.

**Incident risk**: this is exactly the kind of change that would pass a
naive smoke test (`case.setup` runs, `eqquasi` starts) while producing wrong
physics, because Fortran's list-directed `read(unit,*)` has no field names to
catch a shift.

**How to apply**: when adding a BP8 line, edit both files in the same
change, add it as the new last line, and rerun the BP5 regression (rule 3)
to confirm the *unrelated* existing fields still parse correctly (BP5's
`model.txt` will now have one more trailing line than the reference case's
was generated with — confirm `readmodel` doesn't choke on end-of-file for
existing shorter model.txt files stored anywhere, e.g. already-created case
dirs under `test/`).

**Tier**: judgment (2) — requires reading the diff of both files side by
side; not mechanically checkable without a schema, which does not exist (see
Unenforceable rules).

---

## 5. fric()/on_fault_vars magic indices have one authoritative source per direction

`fric(1:100, node, fault)` (allocated at `src/eqquasi.f90:163`) is indexed by
bare integers with no named constants. The same indices are duplicated in
three places that must never disagree:

| Index | Meaning | Python writer | netCDF variable name (case.setup:59-67) | Fortran reader |
|---|---|---|---|---|
| 9 | `a` (RSF) | `scripts/defaultParameters.py:88-94` | `a` | `src/netcdf_io.f90:363` |
| 10 | `b` (RSF) | `scripts/defaultParameters.py:95` | `b` | `src/netcdf_io.f90:364` |
| 11 | `Dc` | `scripts/defaultParameters.py:96-98` | `Dc` | `src/netcdf_io.f90:365` |
| 12 | `v0` | `scripts/defaultParameters.py:99` | `v0` | `src/netcdf_io.f90:366` |
| 13 | `r0` | `scripts/defaultParameters.py:100` | `r0` | `src/netcdf_io.f90:367` |
| 46 | init slip rate | `scripts/defaultParameters.py:102-104` | `init_slip_rate` | `src/netcdf_io.f90:368` |
| 8 | init shear | `scripts/defaultParameters.py:107` | `init_shear_stress` | `src/netcdf_io.f90:369` |
| 7 | init normal | `scripts/defaultParameters.py:106` | `init_normal_stress` | `src/netcdf_io.f90:370` |
| 20 | init state | `scripts/defaultParameters.py:105` | `init_state` | `src/netcdf_io.f90:371` |

A BP8 on-fault pore-pressure state variable needs a new index. Rules:

- **Pick an index not already in use.** Grep `fric(<N>,` across `src/*.f90`
  first (indices 22, 23, 44, 47, 71, 72, 28, 29 are already occupied — see
  `src/faulting.f90:319`, `src/netcdf_io.f90:373`, `src/library_output.f90:164`).
- **Write it in exactly three places in the same change**:
  `scripts/defaultParameters.py` (Python-side array), `scripts/case.setup`
  (`netcdf_write_on_fault_vars`, both the `createVariable` call and the
  `var[:,:] = par.on_fault_vars[:,:,N]` assignment), and the corresponding
  Fortran read in `src/netcdf_io.f90` (mirror the pattern at lines 363-373).
- **Never repurpose an existing index** for a new meaning — it will silently
  corrupt every compset that still writes to it under the old meaning.

**Route**: because there is no single source of truth, any future addition
that misses one of the three files is a doc-vs-code drift bug →
`sophia-okafor`; a mismatched index value across the three (data corruption,
not a doc issue) → `lars-eriksson`.

**Tier**: judgment (2) for new additions; mechanical only in the narrow
sense that `grep -n "fric(<N>," src/*.f90` before assigning a new index is a
checkable precondition.

---

## 6. New bp values are additive branches, never edits to existing bp branches

Existing benchmark selection follows one pattern throughout `src/*.f90`:

```fortran
if (bp == 5) then
    ...
endif
if (bp == 7) then
    ...
endif
```

(see `src/faulting.f90:139,325`, `src/library_output.f90:14,57,161,178`,
`src/solveTimeLoopMUMPS.f90:61`). BP8 must add a new `if (bp == 8) then
... endif` block adjacent to the existing ones in each of these files —
**never modify the body of an existing `bp == 5` or `bp == 7` block** to
accommodate BP8 logic, and never change the default (no-`bp`-match) code
path that BP5/BP7 currently fall through to.

**Rationale**: this is what makes rule 3's oracle meaningful — if BP8 work
can only add branches, a passing BP5 regression after a BP8 change is proof
the two are actually independent, not proof by coincidence.

**How to apply**: `git diff` before committing BP8 work should show 0
modified lines inside any pre-existing `if (bp == 5)` / `if (bp == 7)` block
in `src/faulting.f90`, `src/library_output.f90`, `src/solveTimeLoopMUMPS.f90`.

**Tier**: mechanical — `git diff -U0 src/faulting.f90 src/library_output.f90
src/solveTimeLoopMUMPS.f90` can be scripted to reject any `-`/`+` line pair
falling inside the byte range of an existing `if (bp == 5)`/`if (bp == 7)`
block (requires a small script; not implemented today).

---

## 7. New compsets: naming, placement, and registration

Two parallel registries exist and BP8 needs entries in both, for different
reasons:

- **Production compset** → new directory `case_input/bp8.qdc.<res>/`
  (matching the `<benchmark>.<mode>.<resolution-or-variant>` convention of
  `bp5.qdc.2000`, `bp7.qdc.a.10`, `bp1001.qdc.rough.250`) containing a
  `user_defined_params.py`, **and** an entry appended to
  `case_input/compsets.txt` (currently 7 lines, one compset name per line —
  this file is what the README's "Currently supported compset includes"
  list should be generated from, see rule 10).
- **CI regression compset** → a second, smaller/faster directory named
  `test.bp8.qdc` (mirroring `test.bp5.qdc` vs `bp5.qdc.2000`: the `test.*`
  variant cuts `nstep` from thousands to ~101 and coarsens `dx`, compare
  `case_input/test.bp5.qdc/user_defined_params.py:43` `nstep=101` against
  `case_input/bp5.qdc.2000/user_defined_params.py` `nstep=10000`), added to
  `nameList`/`coreNumList` in `testNameList.py`, and its reference output
  placed at `test.reference.results/test.bp8.qdc/Q0/`.

Note `scripts/create.newcase` (lines 13-27) does **not** validate `compset`
against `compsets.txt` — it just does `os.listdir` + `shutil.copy` from
`case_input/<compset>/`. `compsets.txt` is documentation-only today, not an
enforced allowlist; treat it as such and keep it accurate rather than
relying on it to reject typos.

**Route**: `compsets.txt` drifting from the actual `case_input/` directories
present → `sophia-okafor`.

**Tier**: mechanical (`diff <(ls case_input | grep -v -E 'compsets.txt|user_defined_params.py|^test\.') <(sort case_input/compsets.txt)` — not implemented today, listed under Unenforceable).

---

## 8. Reference results are read-only and only grow, never shrink or overwrite

`test.reference.results/<name>/Q0/*.nc` are oracles (rule 3). Adding BP8:

- Never edit or regenerate an existing `test.bp5.qdc` or `test.bp5.qdc.dip90`
  reference file as a side effect of BP8 work.
- A new `test.reference.results/test.bp8.qdc/Q0/fault.00101.nc` is only
  committed after a run whose `check.test.py` verdict you have personally
  read as `SUCCESS` against itself is meaningless — the first BP8 reference
  is, by definition, unverified against any prior oracle. Say so explicitly
  in the commit/PR: "first reference for BP8, not yet independently
  verified" rather than implying it carries the same evidentiary weight as
  the BP5 zero-diff history.

**Tier**: mechanical (checksum/mtime diff on `test.reference.results/` before
vs. after a BP8-only change should show only additions).

---

## 9. Build/environment claims must be verified on the target host, not assumed

Known landmines, do not rediscover them the hard way:

- `src/makefile` hardcodes `FC = mpif90.openmpi` (lines 13, 31, 48 — `utig`,
  `ubuntu`, `local` targets) and `SCALAPACK_LIB = -L${LIB}
  -lscalapack-openmpi` (lines 26, 40, 56). These names are Debian/Ubuntu
  `update-alternatives` symlinks; on a host without the `openmpi-bin`
  alternatives system (or with MPICH/other MPI installed), only bare
  `mpif90`/`mpirun` exist and the build/run will fail. Before assuming a new
  host works, run `which mpif90.openmpi mpirun.openmpi` and fall back to
  editing the makefile's `FC`/`MPIRUN` for that host only — do not silently
  `ln -s` a fake `mpif90.openmpi` and call it fixed without recording that
  choice in a per-host note.
- `local` target's `SCALAPACK_LIB` (`src/makefile:56`) points at a
  system-wide `-lscalapack-openmpi`, which does not exist if MUMPS was built
  locally with its own bundled ScaLAPACK — check `./mumps/build/local/lib`
  before trusting the `local` target to link.
- `install.eqquasi.sh -m local` cannot currently be validated end-to-end
  because of the syntax bugs in rule 2 — until those are fixed (route to
  `lars-eriksson`), do not treat "`install.eqquasi.sh -m local` completed
  without error" as evidence of anything, since the branch never runs.
- `README.md`'s dependency list omits `netCDF4` even though `scripts/
  case.setup:6` does `import netCDF4 as nc` — `.github/workflows/test.yml`'s
  `pip install` line does include `netCDF4`, so CI works, but a human
  following just the README's "Setup of computing environment" section will
  not have it. This is a doc-vs-code drift, not a rule violation to fix here
  — route to `sophia-okafor`.

**Tier**: mechanical for the `which mpif90.openmpi` check and the
`ldconfig -p | grep scalapack-openmpi` check; judgment for whether a given
host's fix belongs in the shared makefile vs. a local override.

---

## 10. Docs move with the code, in the same change

`README.md` currently has two verified inaccuracies that any BP8-driven edit
to the same sections must not compound:

- "`create.newcase caseDir bp5.qd.2000`" — real compset name is
  `bp5.qdc.2000` (`case_input/compsets.txt` line 1).
- "`create_newcase directoryForYourCase compset`" — real script name (and
  the one used two lines later and in the Installation section) is
  `create.newcase` (`scripts/create.newcase`).

When BP8 is added, `README.md`'s "Currently supported compset includes"
list must gain `bp8.qdc.<res>` in the same commit that adds
`case_input/bp8.qdc.<res>/`, and the two pre-existing typos above should be
fixed if that section is touched anyway (adjacent-line rule: don't leave a
typo one line above your own edit).

**Tier**: mechanical for the two cited typos (`grep -n 'bp5.qd.2000\|create_newcase' README.md`); judgment for "is the compset list still complete."

---

## Unenforceable rules (as written)

| Rule | Why it can't be checked today | What would fix it |
|---|---|---|
| 3 (bit-exact oracle) | `check.test.py` uses `rtol=atol=1e-3`; nothing in the repo asserts the stronger "zero difference" property the project actually relies on for confidence | Add a second, tighter comparison mode (e.g. `atol=0` byte-for-byte or `atol=1e-12`) invoked explicitly for release/tag builds, distinct from the day-to-day 1e-3 smoke check |
| 3 (CI gate reads FAIL) | `.github/workflows/test.yml`'s grep pattern does not match `check.test.py`'s literal `FAIL ` output string | Change the grep pattern to include `FAIL` (case-sensitive, word-boundary), or have `check.test.py` call `sys.exit(1)` on any failure and drop the `\|\| true` in the workflow |
| 4 (model.txt positional contract) | No schema file exists pairing `case.setup` writes with `read_input.f90` reads — verifying order-correctness today means manually diffing two files | A generated/shared schema (e.g. a YAML or Python list of `(name, type)` consumed by both a `case.setup` codegen step and a Fortran read-order check) would make this mechanical |
| 6 (additive-only bp branches) | No script currently maps byte ranges of `if (bp == N)` blocks to a diff and rejects edits inside them | A small AST/regex-based pre-commit check scoped to `src/faulting.f90`, `src/library_output.f90`, `src/solveTimeLoopMUMPS.f90` |
| 7 (compsets.txt registration) | `compsets.txt` is not read by any script (`create.newcase` doesn't validate against it) — it can drift from `case_input/` silently and nothing notices | A CI step: `diff` the compset directory listing against `compsets.txt` and fail on mismatch |

---

## Rule-book health

No prior rule book, stub, or duplicate existed at seed time — this is the
canonical file. Starter-set rules not carried over: rule 9 (cheap targeted
check before expensive run) and rule 11-numbered provenance/benchmark rule
were folded into rule 3 rather than kept generic, since this repo's one
expensive-run hazard (regenerating BP5 references) already has a named gate.
