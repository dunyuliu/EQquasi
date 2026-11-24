# EQquasi v1.24.11 — 2026-11-23

Tag target: `71f6c8d778e86e67c87f157310c92c797d079166` ("Tests: floor comparators agree on non-finite, in-plane floor cases bite, ICNTL(14) source-order contract (#249)"). `EQQUASI_VERSION = '1.24.11'` (`src/globalvar.f90:10`, read at that SHA). Previous tag: v1.24.9 (`f4639024036dacb8a68fabd213e138bc41b416d7`). This notes file is committed on master AFTER the tag target: master commits after 71f6c8d are board text and this file only (`git diff --stat 71f6c8d origin/master` at the time of writing: `PATHWAY_FORWARD.md` only, 17 insertions), so they are not in the tagged tree and the gates below did not run them.

**There is no v1.24.10 tag.** Version 1.24.10 was set by #242 (`b4d8051`, the `-eq_s2t_plant_nan` hook) and 1.24.11 by #244 (`6149bce`, ICNTL(14)). The release gate never ran on `b4d8051` (rule 21), so no tag was cut there; v1.24.11 covers both bumps.

## Summary of scope

`git log v1.24.9..71f6c8d` is 19 commits; `git diff --shortstat v1.24.9 71f6c8d` is 33 files, +27621 / -27025 (almost all of it the regenerated liu2020 gold); `src/` is 3 files, +29 / -1 (`git diff --stat v1.24.9 71f6c8d -- src`); tracked `*.md *.sh *.py *.f90` is 14 files, +636 / -48.

Behaviour a user can see:

- New test hook `-eq_s2t_plant_nan` (#242, `src/solveTimeLoopPETScDD.f90`): plants a NaN in one row of the slip-to-traction operator product so the stop path can be exercised. Inert unless the option is given; an out-of-range row prints a message and stops all ranks with `MPI_ABORT`. Documented in `docs/user/troubleshooting.md`.
- MUMPS workspace: `mumps_par%ICNTL(14) = 100` on the native MUMPS path (`src/solveTimeLoopMUMPS.f90`, #244), after the PETSc and PETScDD loops already set it (v1.24.9), so all three solvers now set it; a user-set value wins on the PETSc paths. One case measured: MUMPS-allocated memory on the largest rank rose from 721 to 979 MB, peak RSS unchanged (4 ranks, one case, host load about 16; UNAUDITED, PR 244).
- Comparator floors for out-of-plane station columns (#232) and the comparators' agreement on non-finite entries (#249): test-side only, no solver change.
- Rule 22 now caps tracked files at 10 MB (#241, `PROJECT_RULES.md`, `tests/contract/test_root_layout.py`); the largest tracked file is `data/liu2020.qdc.kink.600/cycle0/disp.10000.nc`, 6,346,312 bytes.
- `data/liu2020.qdc.kink.600` regenerated from the v1.24.8 run (#245, owner-approved); the reference is a regression lock, no reproduction claim (`data/liu2020.qdc.kink.600/README.md`).
- README: clone step, `eqquasi-petsc` environment, fast-test command (#246).

No performance claim is made by this release.

## Files added / removed / renamed / cleaned up

Added: `tests/contract/test_mumps_workspace_setting.py`, `data/liu2020.qdc.kink.600/cycle0/disp.10000.nc`. Removed, renamed: none. Everything else in `git diff --name-status v1.24.9 71f6c8d` is a modification (liu2020 gold, `src/` 3 files, `tests/` 4 files, `README.md`, `PROJECT_RULES.md`, `PATHWAY_FORWARD.md`, `docs/user/troubleshooting.md`). This release adds `docs/dev/release_notes_v1.24.11.md` after the tag target.

## Content updates to master documents

`README.md` (#246), `PROJECT_RULES.md` (rule 22 cap and largest-file text, #241 and #248), `docs/user/troubleshooting.md` (`-eq_s2t_plant_nan`), `PATHWAY_FORWARD.md` through the conductor 40 gate result (#253).

## Audit findings and fixes

Code audit of `v1.24.9..71f6c8d`, verbatim from the read-only code auditor (2026-11-23), complete, nothing dropped:

> Code audit report, range v1.24.9..71f6c8d (read-only; verbatim from the code auditor, 2026-11-23)
>
> Verdict for v1.24.9..71f6c8d: 1 blocker, 0 major, 2 minor, 2 low. The only blocker is the one the previous audit found: no tag exists. No code blocker.
>
> Scope: `git diff --stat ae974dd 71f6c8d -- src` is empty, confirming that `src/` did not change after ae974dd. The `src/` delta for the whole range is 3 files and +29/-1 lines: the version bump to 1.24.11, `mumps_par%ICNTL(14) = 100` on the native MUMPS path, and the `-eq_s2t_plant_nan` test hook. The checks below were run at 71f6c8d.
>
> Blocker
> - `src/globalvar.f90:10`, blocker: `EQQUASI_VERSION = '1.24.11'` and the comment above it says the string "must match a tag that exists", but no `v1.24.10` or `v1.24.11` tag exists locally or on origin. This release closes it, and it is the only blocker.
>
> Minor
> - `docs/user/troubleshooting.md:77-82`, minor: the new entry says the operator check's NaN flag stops a `-eq_s2t_plant_nan` run. The code comment (`src/solveTimeLoopPETScDD.f90:634-637`) and the test's expected message (the `nan_row0` case in `tests/e2e/test_benchmarks.py`) say a different check stops it first: the finite check on the operator product S g. That check's abort message, "EQ_S2T: NaN or Inf in the operator product S g, step" (`src/solveTimeLoopPETScDD.f90:1615`), is not documented anywhere in `docs/user/troubleshooting.md`.
> - `tests/e2e/cases.py:663-666` and `tests/e2e/test_benchmarks.py:523-524`, minor (carried over; the earlier audit rated it major): the 9-column station floors still apply to columns 2, 4 and 6 (the dip columns) based only on the column count. Those columns are out-of-plane only for a vertical fault. PR 249 pinned this behaviour with `test_compare_series_fltst_nine_column_layout_gets_the_floors`; the code is unchanged. So this was not fixed in PR 249. Downgraded because every gold case with a 9-column layout (BP5/BP7) has a vertical fault, so no gold comparison is affected. A future case with a dipping fault would compare in-plane dip columns with these floors.
>
> Low
> - `data/liu2020.qdc.kink.600/cycle0/runInfo.json:6` and `data/liu2020.qdc.kink.600/README.md:9-21`, low: the regenerated reference records the build host's machine name as provenance. The superseded set also did this. The file is written by the solver; no filesystem path or personal name in it.
> - `data/liu2020.qdc.kink.600/cycle0/user_defined_params.py:2,88`, low: the header says a cycle takes "~5 h", the README says ~6 h (measured 6.19 h). The dx comment gives Lambda0/dx = 1.26, the README gives 0.94.
>
> Acknowledged, not findings
> - `tests/contract/test_mumps_workspace_setting.py` checks only the order of statements in the source, not a run, and its docstring says so. Run-time evidence is in PR 244: MUMPS-allocated memory on the largest rank rose from 721 to 979 MB (not audited).
> - No e2e test reaches the operator check's NaN flag `v(7)`; the board records this.
> - The 13 e2e skips ("this reference holds no restart/offfault/perfault/profile files", `tests/e2e/test_benchmarks.py:97`) are a data gap that existed at v1.24.9. No gold file that decides a skip changed in this range.
> - The liu2020 reference was produced by v1.24.8, before the ICNTL(14) change; ICNTL(14) only sets MUMPS workspace size, so it does not change the numbers. The README claims a regression lock only, no reproduction claim.
>
> Checked and clean
> - ICNTL(14) on all three solvers: MUMPS `JOB=-1` (line 17), then `call DMUMPS`, then `ICNTL(14)=100` (line 22), then `JOB=4` (line 53); factorization error check at line 64. PETSc (`:430-439`) and PETScDD (`:577-583`): set between `PCFactorSetMatSolverType` and `KSPSetFromOptions`; a user-set value wins.
> - NaN hook (`src/solveTimeLoopPETScDD.f90:634-655`): inert unless the option is set; an out-of-range row prints a message and stops all ranks with `MPI_ABORT` (a test case covers this); only the owning rank writes the NaN, with the same local indexing as `-eq_s2t_plant_row`. PR 242 records 4 plant cases passing and a fast tier of 384 passed; not re-run by the auditor (needs a PETSc build).
> - Comparator floors: `station_floor_exempt` and `_dd_station_series` refuse a floor on a non-finite entry; `_dd_station_series` rejects any non-finite reference value (`:533-536`); the solver clamps log10 slip rate at 1e-30 (`src/library_output.f90:138-139`). In-plane columns and `srfst_` files get no floor, and tests check that a difference above a floor still fails. Run by the auditor: 29 comparator tests (`-m e2e_fast -k "station_series or compare_series"`) passed.
> - Rule 22 cap: rule text and `MAX_TRACKED_FILE_BYTES = 10*1024*1024` agree; largest tracked file `disp.10000.nc` 6,346,312 bytes. liu2020 gold is in `REGENERATED_8BYTE`. Contract tests for the MUMPS setting, root layout and gold netCDF shape: 104 passed, 0 skipped.
> - No fallbacks, placeholders or silent failures added in `src/`, `tests/`, `scripts/`. No machine-local paths or personal names in added lines, apart from the liu2020 hostname above; "Liu, Duan & Luo 2020" is a literature citation.
> - README additions (clone step, `eqquasi-petsc` environment, fast-test command) match `CLAUDE.md` and `pytest.ini`.
>
> Sign-off rests with the human reviewer.

Reading of that audit (release engineer, not part of the quote):

- The blocker (no tag) is closed by this tag.
- The 9-column floors finding was NOT fixed in PR 249: that PR only pinned the behaviour with a test. It is now minor and stays open.
- No fix was applied in this release; it tags an existing commit, so the only new file is this one.

Rules audit (from the board, `PATHWAY_FORWARD.md` conductor entry, quoted): "rules audit (5 findings: PR 244 had no release audit, PR 242/244 bodies lack a full `-m e2e_fast` line, PR 245 replaced files where rule 8 says re-blessing is additive, rule 22 text stale, 3a listing)". None is a code defect. The rule 22 text was corrected in #248.

Full-tier gate on the tag target, `runs/20261123_c40_gate_71f6c8d/` (fresh clone of 71f6c8d under `env -i`): `status.txt` reads `build rc=0`, `fast rc=0` (`fast.log`: 401 passed, 190 deselected, 94.24 s), `e2e_full rc=1` at 21:47:29; `e2e_full.log` tail: `1 failed, 143 passed, 46 skipped, 401 deselected, 3 warnings in 13875.88s (3:51:15)`.

- The one failure, `test_install_script_produces_a_working_binary` ("EQquasi build FAILED."), was a harness defect: `MACHINE` was unset in the gate script's environment, so the install defaulted to `ubuntu`. Re-run alone on the same clone and SHA with `MACHINE=conda-linux`: `1 passed, 590 deselected in 19.76s`, rc 0 (`install_rerun.log`; the conductor's own re-run took 19.59 s). The gate script now exports `MACHINE`. The full tier was not re-run end to end after that fix; the re-run is the single test.
- Stranger clone of 71f6c8d under `env -i` following the README, `runs/20261123_c40_stranger/`: clone rc 0, steps rc 0 (build rc 0, `create.newcase`, `case.setup`), `run.sh` stopped by a 150 s timeout after printing `Welcome to EQquasi 1.24.11`. README gap found: `conda env create -f environment.yml` fails with `CondaValueError: prefix already exists` on a host that already has `eqquasi-petsc`, and that run used the existing environment, not a fresh one.

## Remaining open issues or pending items

1. Minor: `docs/user/troubleshooting.md:77-82` names the wrong check for a `-eq_s2t_plant_nan` stop; the abort message "EQ_S2T: NaN or Inf in the operator product S g, step" is undocumented.
2. Minor: 9-column station floors apply to the dip columns by column count (`tests/e2e/cases.py:663-666`, `tests/e2e/test_benchmarks.py:523-524`); wrong for a dipping fault; no gold case affected. Not fixed in #249.
3. Low: `data/liu2020.qdc.kink.600/cycle0/runInfo.json:6` records the build host's machine name.
4. Low: `data/liu2020.qdc.kink.600/cycle0/user_defined_params.py:2,88` says "~5 h" and Lambda0/dx = 1.26; the README says ~6 h and 0.94.
5. 46 skips in the full tier, all absent reference data, not covered by rule 3a's binary clause (from `e2e_full.log`, `grep '^SKIPPED'`): 11 at `tests/e2e/test_benchmarks.py:80` (no reference for `data/test.bp1002.qdc.zone.2500/cycle0`); at `:97` no restart files 1, no offfault files 4, no perfault files 8, no profile files 9; at `:126` no `summary.json` 4; at `:176` no section 4.3 profiles 9. The owner's earlier figure of 13 skips (and the audit's "13 e2e skips") differs from 46; the tier selection differs and the two are not reconciled. The v1.24.9 gate had the same 46.
6. Rules-audit items above (5 findings), as open.
7. Unaudited: the ICNTL(14) memory figures (721 to 979 MB largest-rank MUMPS allocation; peak RSS unchanged) rest on one case, 4 ranks, host load about 16.
8. README gap: `conda env create` fails where the environment already exists.
9. The full tier was not re-run end to end after the gate-script `MACHINE` fix (see above).
10. No `tests/release_gate.sh` exists in the repository; the five Release gate rows below were walked by hand.

## Totals or cost changes

Gate on the tag target: fast tier 401 passed; full tier 143 passed, 46 skipped, 1 harness-defect failure re-run green. Wall: fast 94 s, full 3:51:15, install re-run 19.76 s. No solver speed or memory change is claimed beyond the single ICNTL(14) measurement above.

## Assumptions used

- The owner-decided regeneration of the liu2020 gold (#245) stands.
- The conductor's decisions for this release: audit quoted as above; no v1.24.10 tag; annotated tag on 71f6c8d.
- The 46 skips are data absence, as at v1.24.9.

## CI run this release was gated on

Run 45283376909, workflow "Automatic Testing of EQquasi", event push, conclusion `success`, SHA `71f6c8d778e86e67c87f157310c92c797d079166` (`gh run view 45283376909 --json conclusion,headSha`, re-read 2026-11-23); it is the only run `gh run list --commit 71f6c8d...` returns.

## Trend since v1.24.9

The project has no `tests/check.sh`, `evals/run.sh` or board VERIFIED/BROKEN columns, so three of the five measures do not exist here and are not invented.

- Gate assertions: not measured (no `tests/check.sh`). Full tier at v1.24.9: 119 passed, 46 skipped, 0 failed (v1.24.9 Release body); now 143 passed, 46 skipped, 1 harness failure re-run green. Skips unchanged.
- Fixture verdicts: not measured (no `evals/run.sh`).
- Tracked text lines, `git diff --shortstat v1.24.9 71f6c8d -- '*.md' '*.sh' '*.py' '*.f90'`: 14 files, +636 / -48. Gate numbers improved (24 more passes, same skips), so the growth is matched by tests and docs, not unexplained.
- Board currency: not measured (no VERIFIED/BROKEN columns).
- CI green-on-first-try, `gh run list --limit 200` since 2026-11-19T22:50Z: 81 runs, 58 success, 22 cancelled (superseded by the workflow's concurrency group), 0 failure, 1 not yet concluded, 0 re-runs (attempt above 1). Of the 58 concluded and not cancelled, all 58 were green at attempt 1.

## Work record

- audit: code audit of `v1.24.9..71f6c8d`, quoted in full above: 1 blocker (no tag; closed by this tag), 0 major, 2 minor, 2 low.
- correctness: no code blocker; NaN hook and ICNTL(14) order checked by the auditor at file:line; gate on the tag target green apart from the harness-defect test, re-run green.
- conciseness: not assessed by a dedicated pass; the release changes no code, and the `src/` delta is +29 / -1.
- fixes: none applied; deferred findings are open issues 1 to 4 above.
- docs: this file checked against the filesystem and `git diff --name-status v1.24.9 71f6c8d`; open README gap recorded; no release-notes file existed under `docs/dev/` before this one.
- refactor: not run; no refactor pass was dispatched, since this release tags an existing commit and changes no code.
- rules: rules audit from the board (5 findings, none a code defect); no separate rules-enforcer run in this release.

## Release gate

- tree: release branch cut from `origin/master` in its own worktree; tag target 71f6c8d is an ancestor of master (`git merge-base --is-ancestor`, exit 0); only `PATHWAY_FORWARD.md` and this file follow it.
- ci: run 45283376909 success on 71f6c8d778e86e67c87f157310c92c797d079166.
- publish: note v1.24.11, annotated tag v1.24.11 on 71f6c8d, remote origin; the tag exists only after this note merges, so its `git ls-remote --tags origin` line is in the release report, not here.
- release: created with `gh release create v1.24.11 --verify-tag --latest --notes-file docs/dev/release_notes_v1.24.11.md` after this note merges; `gh release view` is in the release report, not here.
- clone: PASS 71f6c8d778e86e67c87f157310c92c797d079166 before the tag: clone of the pushed SHA into an empty directory under `env -i HOME=<empty dir>` (`runs/20261123_c41_release_clone/`), README steps: `make.scripts.executable.sh`; `conda env create -f environment.yml` run with `-p <dir>` added, because the default `eqquasi-petsc` already exists on the host, so the environment was built fresh from the manifest (22:02:21, rc 0); `install.eqquasi.sh -m conda-linux` rc 0; fast tier `401 passed, 190 deselected in 99.91s`; `create.newcase` + `case.setup` rc 0; `run.sh` printed `Welcome to EQquasi 1.24.11` and ran to the 150 s timeout (no full cycle). Not run: the full `-m e2e` tier from this clone (it ran from the earlier gate clone, above).
