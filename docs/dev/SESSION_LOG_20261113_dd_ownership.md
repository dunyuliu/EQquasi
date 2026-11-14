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
