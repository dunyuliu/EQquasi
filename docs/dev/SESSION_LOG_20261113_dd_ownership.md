# Session log 2026-11-13: board row 1, distributed PETSc loop (conductor handoff at ~80 tool calls)

Times CST. Evidence paths are under `runs/` of the main checkout (gitignored).

## Landed
| PR | Squash SHA | CI run | What |
|---|---|---|---|
| #75 | 5896868 | 43520841041 success | `test.yml` pins Python 3.12: `'3.x'` resolved to 3.15, no wheels for pinned numpy, `pip` built sdists and `Install dependencies` stalled on every PR |
| #73 | 14c3033 | 43521617143 success | board row 1 spec, steps (e)-(h) |
| #74 | 9bfc743 | 43522459530 success | board row 6, xi above the stability limit |
| #72 | df21bc6 | 43534844531/-37 success | v1.23.0 per-node xi (eq. 15); audit MAJOR (no test ran the path) closed by `test_xi_gamma_per_node_step_runs_and_matches_eq15`, my own run 1 passed rc 0 (`runs/20261113_gate-72/xi.log`); own fast tier 287 passed; own `-m e2e_fast` 43 passed (install test needs `MACHINE=conda-linux`), 21 skipped, all "reference holds no <profile/perfault/offfault/restart/summary.json> files" |

## Open PRs (not merged)
- #76 `dd-ownership` e5d74f8, ready, v1.22.0 source, solver 3 loop (`src/solveTimeLoopPETScDD.f90`). Audit: no blocker; MAJOR = version bump after #72 (now due: 1.24.0). Own parity check: bp5 4 ranks worst rel 9.3e-14 vs solver 2 (`runs/20261113_gate-76/`). Holds: owner has NOT confirmed `par.solver = 3` and the new file name; CI cannot reach solver 3 (CI builds `-m ubuntu`, no PETSc), only local `-m e2e_fast` does.
- #77 `dd-mesh` 39ec587, draft, base `dd-ownership`, stage 2a per-rank mesh. Audit: no blocker; MAJORs: behind master (1.23.0; trial merge conflicts in `tests/e2e/test_benchmarks.py`), rank-local `stop 7` at `src/eqquasi.f90:159-161,185`, comparison omits `srfst_*`, mesh files, `disp.*.nc`. A fix agent is on the last two (see still-owns).

## bp1002 baseline, `test.bp1002.qdc.zone.2500`, 101 steps, v1.22.0 + #76 head e5d74f8 (NOT host-exclusive: load 7-30 during runs, `runs/20261113_baseline/load.log`)
rank-0 cpu_time, seconds:
| ranks | solver 1 time loop | solver 3 time loop | solver 1 factor | solver 3 factor |
|---|---|---|---|---|
| 8 | 62.85 | 13.73 | 3.19 | 3.21 |
| 16 | 62.25 | 13.54 | 2.60 | 2.39 |
| 32 | 60.91 | 13.10 | 2.28 | 2.13 |
Time loop is flat in rank count for both (the triangular solves do not scale at this size); solver 3 is about 4.5x lower per step at every rank count. Contended host, one size: indicative only. BP8 10 m wide box not run: solver 3 refuses bp == 8 (stage 2c), native MUMPS there is not started.

## bp7 dip component (item 4)
`fault.00001.nc`, `test.bp7.qdc.a.10`, solver 2, 1/4/8/16 ranks: `shear_dip` max 1.0e-18 Pa vs `shear_strike` 1.6e7 Pa; differences between rank counts 3.1e-19..3.3e-19 (same size as the values); all other fields bit-identical. Symmetry noise present in solver 2 itself. Solver 3 on bp7 is in my gate job below.

## Still owns (as of this entry)
- (closed 19:40) Gate job on #77 39ec587 finished, `runs/20261113_gate-77/rc` all 0: build ok, fast tier ok, `-m e2e_fast` 45 passed / 21 skipped (same data-absence skip reasons, native-path regression included), bp7 solver 3 at 4 ranks vs solver 2 at 4 ranks (`runs/20261113_baseline/`): worst rel 3.8e-13 (fltst), 9.9e-15 (p?output), 1.6e-14 (global.dat), 1.2e-10 on the 7-digit `srfst_*` files; float32 netCDF worst 0.448 normalised by field max, read as the dip-component symmetry noise (~1e-18 Pa, same size as solver 2 rank-to-rank) but the field was not identified: check it before claiming bp7 covered.
- Fix agent senior-researcher agent task `a1f095e6caffd6cc1`, worktree `.claude/worktrees/agent-a1f095e6caffd6cc1`, branch `dd-mesh-fix`, pushes a fast-forward to `dd-mesh`; queued on the heavy lock (parity run PID 873305). Do not dispatch onto `dd-mesh` meanwhile.
- Worktrees: `board-dd`, `ci-pin`, `owner-relays` (all merged, reap after confirming no live child), `gate-77`, `log-handoff`, `baseline-76` (reap; evidence is in `runs/`), `agent-a5ec6fe5d7e4cc36c` (copy `runs/` and `scratch/dd-mesh/` out before removal: its runs and notes live inside it).
- Heavy-run lock: `scratch/locks/heavy.lock`, `flock`, timeout INSIDE the flock (a timeout around flock expired while waiting and produced a false failure once).

## Next, in order (fresh conductor)
1. Owner confirms `par.solver = 3` and `src/solveTimeLoopPETScDD.f90` (not decided; ask).
2. Rebase `dd-ownership` on master, bump 1.24.0, rebase `dd-mesh` onto it, retarget #77, re-gate each (fast tier, own `-m e2e_fast` with `MACHINE=conda-linux`, auditor on the delta only), merge #76 then #77 serially.
3. Stage 2b (per-rank output, no rank-0 gather), then 2c (bp == 8 pore pressure), each a stacked PR, one agent at a time.
4. Follow-up MINORs, none blocking: #72: test asserts nonzero exit not exactly 13; `parameters.md` says "each node" but creeping nodes keep constant xi; `case.setup` has no upper bound while `read_input.f90:232` caps at 1e6; mu vs mu/(1-nu) in eq. 15 unchecked against the paper. #76: off-fault `srfst_*` untested, `if (j <= 0) j = 1` silent fallback at `eqquasi.f90:74,86` still live for solvers 1/2, test does not assert rank count. #77: `check()` exit at `netcdf_io.f90:771`, `ddPeakRSSkB` -1 line, 501/502 text.
5. Rule-book lesson for the rules agent: no floating `python-version: '3.x'` in workflows (`docs.yml:32` still has it); a wrapped `timeout` around `flock` is not a test timeout.

## Owner rows (unchanged, not decided here)
`par.solver` meaning at (h); PETSc hard dependency on ls6/utig/ubuntu; retiring `solveTimeLoopMUMPS.f90` (major bump); row 6 sub-items (a)-(c); carried-over list in the brief.

## Entry 2, evening of 2026-11-13 (second conductor, tool calls 1-34 of 100)
- #76 rebased on master 0f06403 (one test-file add/add conflict, both blocks kept), v1.24.0, head 0974fbd. Gate (conda-linux build, `eqquasi-petsc` env): fast tier 287 passed, `-m e2e_fast` 46 passed / 21 skipped, all skips data absence (no summary.json, no section 4.3 profiles, no offfault/perfault/profile/restart files in the references). `runs/20261113_gate-76b/`, CI green.
- #77 fix agent finished (ba0cb49: rank-local stops abort every rank under solver 3; test compares srfst, mesh, disp). Audit of that delta: PASS, one MINOR (test comment claims double precision at 1e-10; fltst/srfst columns are 7 digits, effective ~2e-6), folded into stage 2b. `dd-mesh` rebased on #76, v1.24.1, head 910ec98, base `dd-ownership`. Gate: fast 287 passed, `-m e2e_fast` 46 passed / 21 skipped, same skips. `runs/20261113_gate-77b/`, CI green.
- bp7 float32 open point closed: the fields differing by 0.39-0.41 of their own maximum at solver 3 np4 vs solver 2 np4 are `shear_dip` (max 1.2e-18 Pa), `vzm`, `vzs` (4e-35), `slipd` (8e-38) in `fault.00001.nc`; solver 2 np1 vs np4 differs by 0.31 on the same four. Dip component is symmetry noise in both solvers. All other netCDF fields agree to <= 1.9e-12 of their maximum (`disp.00101.nc`). bp7 is covered.
- The NumPy 2.2.6 failure of `plotMagnitudeTime.py` in the contract tests occurs with `/usr/bin/python3` on master as well; gates run in the `eqquasi-petsc` env.
- Stage 2b dispatched (senior-researcher agent, branch `dd-2b` off `dd-mesh`, draft PR base `dd-mesh`, v1.24.2): per-rank output, no field gather onto rank 0.
- Reaped: worktrees board-dd, ci-pin, owner-relays, log-handoff, gate-77, baseline-76, both stage agents' worktrees; their evidence is in `runs/20261113_agent-a5ec/` and `runs/20261113_dd-mesh-fix/`.
- Owner decision still open: `par.solver = 3` as selector and `src/solveTimeLoopPETScDD.f90`; #76 and #77 are held on it.

## Entry 3, 2026-11-14 (second conductor, tool calls 35-77 of 100)
Stack, all draft or open, none merged (owner decision on `par.solver = 3` and `src/solveTimeLoopPETScDD.f90` still open): #76 `dd-ownership` 0974fbd (v1.24.0) <- #77 `dd-mesh` 910ec98 (1.24.1) <- #82 `dd-2b` 77e04b1 (1.24.2) <- #83 `dd-2c` ae7e3b1 (1.24.3).
- #82 stage 2b (per-rank output: MPI-IO for fixed-width text, owner-written station files, rank 0 as single netCDF writer one rank part at a time). Code audit PASS, 4 MINORs folded into 2c. Own gate: fast 287 passed, `-m e2e_fast` 46 passed / 21 skipped, rc 0 (`runs/20261113_gate-82/`), CI green.
- #83 stage 2c (BP8 pore pressure on owned fault nodes, `ddPfHalo`, profile-line slots, bp7 and bp8 added to the parity test, `stop 508` aborts under solver 3). Code audit PASS (physics read as identical to solver 2; coverage gap: only blocks (2,2,1), so z-halo, in-plane corners, depth-line slots across ranks and the `fluid_src = 2` well path have no committed test; expected-kinds check wants one file per prefix). CI went RED at bfa0fcf: `-m ubuntu` link, `ddblocks_` defined only in the PETSc-only object; the local conda-linux gate links it. Fix `ae7e3b1` moves `ddBlocks` into `porepressure.f90`; CI green (run 43574544861); own gate at ae7e3b1: fast 287 passed, `-m e2e_fast` 48 passed / 21 skipped (data absence, same reasons), rc 0, `runs/20261114_gate-83b/`.
- Lesson for the rule book (for the rules agent): any `src/` change is gated by a link of the non-PETSc object set as well, or by CI on that head, before it is called gated; the conda-linux local gate cannot see a symbol that only the PETSc objects define.
- Open, workable: intermittent `KSPConvergedReason = -11` at the first solve, 1 of 11 identical 8-rank bp8 solver 3 runs (author's report, failing log not kept). Control loop: 20 solver 3 and 20 solver 2 runs, 8 ranks, 3 steps, `PETSC_OPTIONS=-mat_mumps_icntl_4 3`: 40 of 40 rc 0, no -11 (`runs/20261114_ksp-stress/results.txt`). Not reproduced, not attributed. `KSPSetUp` at `solveTimeLoopPETScDD.f90:567` does not read the MUMPS failure code (solver 1 reads INFOG at `solveTimeLoopMUMPS.f90:60,153`): keep the failing run's stdout if it recurs.
- Open, workable, not started: stage (f) memory evidence (peak RSS per rank, solver 3 vs 2, bp1002 and 10 m BP8 wide at 8/16/32 ranks; needs a quiet host), (g) `MATIS` + `PCBDDC`, tests for the uncovered halo paths above, `ddBlocks` still includes `mpif.h` beside `use mpi`.
- Over cap: the first stage 2b agent used 86 tool calls against a stop line of about 80; the second stage 2c run was cut by a model rate limit and resumed.
