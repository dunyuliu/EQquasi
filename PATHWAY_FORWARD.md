# Pathway forward

Open work only. Updated 2026-11-13, at v1.22.0.

Close an item by deleting it, and record what was learned where it will be
read again: `PROJECT_RULES.md` for a rule, `data/<bench>/README.md` for a
benchmark finding, a commit message for a fix. Closed rows and the findings
behind them (BP8 domain and boundary study, step-over science, caps, known
limits) are in `docs/dev/board_history.md`. Rules are in `PROJECT_RULES.md`;
build, run and traps are in `CLAUDE.md`.

---

## Open

1. [ ] **P1 -- PETSc scaling in EQquasi.** Owner, 2026-10-09. Known today:
   native MUMPS scales 1.9x over 1->8 ranks (BP5); 10 m BP8 wide box ~15 s/step
   at 8 ranks, ~79 GB of factors; CG+GAMG scales better but is slower at every
   size tried (1.4x bp1002, 7.7x kink300); faulting and BP8 pore pressure run
   on rank 0 only. Steps, in order:
   (a) profile per-step time (solve / faulting+pressure / communication / I/O)
   at 8, 16, 32 ranks on 10 m BP8 wide and bp1002, host to ourselves;
   (b) MUMPS block-low-rank compression, ParMETIS ordering, threaded BLAS per
   rank; (c) hypre BoomerAMG vs tuned GAMG with the setup reused across steps
   (K is constant); (d) distribute the rank-0 work only if (a) says it matters.
   Done when: 10 m BP8 wide <= 5 s/step at 32 ranks, or a measured reason it
   cannot be, with parity to the gold at every change.
   Evidence: the per-step table from (a) and the final run's runInfo.json.
   **(a), re-verified 2026-10-12** (binary eqquasi-1.21.2, host theo4,
   host-exclusive, 10-step profiling runs, native MUMPS `par.solver=1`;
   `runs/20261009_petsc-profile/{bp1002,bp8wide}-r{8,16,32}/result/cycle0/runInfo.json`):
   bp1002.qdc.2500 (127k eq) -- 8r 0.436 s/step, 16r 0.409, 32r 0.394 s/step
   (flat; rank-0-only faulting/pressure work dominates at this size). BP8 10 m
   wide (3.52M eq, scratch-only `fxmin/fxmax/fzmin/fzmax=-1500/1500`,
   `enlarging_ratio_xz=1.2`, NOT the committed `compset/bp8.qdc.gs.10` which
   is a smaller -500/500 box) -- 8r 11.01 s/step, 16r 8.01, 32r 6.91 s/step.
   Correction to a prior conductor's note: an earlier pass flagged these as
   possibly contaminated by a concurrent external job and claimed a rerun was
   needed; on re-running r16 and r32 clean (host confirmed idle, `uptime`
   load 0.35, no other eqquasi/python procs), the new numbers (8.01, 6.91)
   match the flagged-contaminated ones (7.98, 6.76) to ~2-4%, i.e. within
   normal run-to-run noise -- the contamination was not material here.
   (b)/(c)/(d) not started. Finer solve / faulting+pressure / communication /
   I/O split (as (a) asks) is NOT yet instrumented -- only aggregate
   time-loop and factorization time exist today (`solveTimeLoopMUMPS.f90`'s
   `cpu_time` calls). Recommendation (not acted on -- rule 1, new-compset
   naming needs agreement): register the -1500/1500 wide box as its own
   compset (e.g. `bp8.qdc.gs.10.wide`) instead of leaving it only in `scratch/`.
   **Owner direction, 2026-11-13**: "i think we can domain decompose from the
   beginning with petcs?"; on the three changes below, "all three changes
   sounds great". Why: the time loop is two triangular solves per step and
   those do not scale (above; Wang 2026 Case 3a, 3.61M eq, 24 ranks, v1.22.0:
   factorization 137 s, time loop 15093 s per cycle, ~6.5 s/step). The PETSc
   path today is not a domain decomposition: rows are split by `PETSC_DECIDE`,
   rank 0 builds the whole RHS and receives the whole solution
   (`VecScatterCreateToZero`), and the default solve is still MUMPS. The
   three changes, in the Fortran PETSc path, replacing (c)/(d):
   (e) each rank owns a block of the mesh: it generates only its own
   elements (the generator is analytic, no partitioner), assembles its own
   rows and RHS, and runs friction and pore pressure on its own fault nodes;
   (f) memory scales out: no rank holds the full mesh or the global factors
   (79 GB at 10 m BP8 wide, 141 GB in-core estimate at Case 3);
   (g) subdomain matrices kept unassembled (`MATIS`) and solved with
   `PCBDDC`, the setup built once and reused every step (K is constant).
   Owner, 2026-11-14, on the solver under (g): "This is a really simple system.
   It's symmetric. You can use maybe just a regular CG or GRMS." (verbatim;
   read as: try plain Krylov, CG or GMRES, on the distributed matrix, not
   only BDDC with direct local solves.)
   Owner, 2026-11-14: "If you cut the domain into subdomains, and each
   subdomain you can use mumps with much less memory cost, then it also
   sounds a good strategy." (verbatim; read as: a subdomain-wise direct
   solve under an outer Krylov iteration is also a candidate.)
   Owner, 2026-11-14: "you can just factorize once and just solve for each
   subdomain and then piece them together. But you just factorize once."
   (verbatim.)
   Owner, 2026-11-14: "was a really good strategy and put it on board as
   well. So in the end, we may have several solver options. And I think the
   last that you proposed, this 3-step one, is pretty good. Give it a try."
   (verbatim.) The 3-step strategy, exact substructuring with no outer
   iteration: (1) factorize each subdomain interior once; (2) build and
   factorize once the Schur complement that couples only the interface
   unknowns; (3) each step is one solve per subdomain, one interface solve,
   one more solve per subdomain. Known cost: the interface matrix is dense,
   its size not yet computed. Status: TO TRY as a spike, beside plain CG and
   subdomain-direct solves under CG; none measured yet.
   (h) end state, owner 2026-11-13: "in the end, mumps is a choice o petcs".
   One PETSc time loop; MUMPS, BDDC and AMG are solver choices inside it and
   `solveTimeLoopMUMPS.f90` retires. Before it retires: PETSc-with-MUMPS
   matches the gold on every gated compset and is no slower per step than
   native MUMPS; PETSc builds on every supported machine (today only
   `conda-linux` links it); one assembly path serves both `MATIS` and AIJ;
   what `par.solver` then means is agreed first (rule 1).
   Two design principles for all of (e)-(h), owner 2026-11-13: "things
   parallelize from beginning, and thing that can armotize once should be
   done that way".
   Parallel from the beginning: no rank-0 phase anywhere -- mesh, numbering,
   assembly, initial conditions, friction, pore pressure, station and fault
   output are all per-rank.
   Amortize once: whatever does not change is built once and reused -- the
   solver setup across steps (K is constant); mesh, partition and setup
   across the cycles of a run (each cycle rebuilds and refactorizes today,
   137 s per cycle at Case 3); and, as a candidate to measure, the fault
   slip-to-traction operator (row 4's note), which would amortize the volume
   solve itself across steps and across runs on one mesh and modulus.
   Design choices shown to the owner 2026-11-13, answer "go":
   a new distributed time loop beside the two existing ones, which retire
   only at (h), so the gold gate stays meaningful throughout; a regular
   x-y-z block decomposition from the analytic generator, not `DMPlex`;
   split fault nodes may cross subdomain boundaries, with a ghost exchange
   per step for friction; the fault slip-to-traction operator is a measured
   side experiment, outside (e)-(h); an iterative solve is held to parity
   within a stated tolerance, shown not to move event times. Still to
   agree at (h): what `par.solver` means, and PETSc as a hard dependency.
   Evidence: parity to the gold at each of (e)-(h), and the per-step table
   against native MUMPS at 8/16/32 ranks on bp1002 and 10 m BP8 wide.
   Unmeasured: whether BDDC beats MUMPS at these sizes (CG+GAMG did not).
   Host, owner 2026-11-13: "do it on cotopaxi" (64 cores, the
   `eqquasi-petsc` env is on the shared home). Timings there are not
   comparable with the theo4 numbers above: re-measure native MUMPS on
   cotopaxi as the baseline first.
   **State, 2026-11-14** (MEASURED, not VERIFIED; no independent oracle).
   (e) stack, nothing merged, held on the owner's answer for `par.solver = 3`
   and `src/solveTimeLoopPETScDD.f90`: #76 (v1.24.0) <- #77 <- #82 <- #83
   <- #85 (v1.24.4: coverage of z and corner blocks, the well path, an abort;
   `nsmp` padding zeroed). #85 own gate at 81629ad from a fresh clone: fast
   287 passed, `-m e2e_fast` 53 passed / 21 skipped (data absence in the
   reference), CI run 43583429026 success. (g) `MATIS` + `PCBDDC`, spike only (draft PR #87, not for merge): see below.
   Evidence, binary eqquasi-1.24.3 (`dd-2c` ae7e3b1), host cotopaxi, NOT
   exclusive (1-min load 4.1-32.9 before each run, all users, rank-0
   `cpu_time`; `runs/20261114_row1-evidence/table.csv`, `matrix.out`,
   `load.log`). bp1002 101 steps, time loop s; wide = 10 m BP8 wide, 10 steps,
   s/step; solver 1 native MUMPS, solver 3 distributed loop:
   | case | ranks | s1 | s3 | max RSS per rank GB, s1 / s3 |
   |---|---|---|---|---|
   | bp1002 | 8 | 66.5 | 20.2 | 0.64 / 0.49 |
   | bp1002 | 16 | 65.6 | 18.6 | 0.56 / 0.33 |
   | bp1002 | 32 | 64.2 | 14.5 | 0.44 / 0.21 |
   | wide | 8 | 14.80 | 11.01 | 25.7 / 21.0 |
   | wide | 16 | 11.58 | 9.01 | 16.7 / 11.9 |
   | wide | 32 | 9.43 | 6.36 | 11.2 / 6.9 |
   Single rank, solver 3: bp1002 56.9 s loop, 1.95 GB; wide 101.3 GB, 62 s/step
   over 3 steps, factorization 3269 s. (f): the largest rank holds 21 / 12 /
   7 % of the single-rank total on wide at 8 / 16 / 32 ranks (25 / 17 / 11 %
   on bp1002 at 0.49 / 0.33 / 0.21 of 1.95 GB), so no rank holds the whole
   problem; the sum over ranks is 1.5-1.8x the single-rank total (150 / 176
   / 185 GB), the factors' fill grows with ranks. Peak RSS is the process
   total (MUMPS factors dominate), the mesh-only share is not separated.
   Not at the target: wide solver 3 at 32 ranks is 6.36 s/step against
   <= 5 s, on a loaded host. bp1002 solver 1 is flat in rank count (66.5 /
   65.6 / 64.2); solver 3 falls 28 % over 4x ranks but a conductor's earlier
   run of the same case measured 13.7 / 13.5 / 13.1 at lower load (`runs/20261113_baseline/`, binary 1.22.0), so the
   rank-count trend is not established. What dominates bp1002 at 8 ranks:
   KSPSolve is 93.7 % of the solver 3 time loop (`runs/20261114_prof/results.txt`; 95 % was the wide 8-rank share). A
   host-exclusive repeat is needed before any scaling claim.
   (g) BDDC spike (CG + `PCBDDC`, MUMPS local and coarse solves, `-eq_mat_is`,
   branch `dd-3g-spike2` 7b77221; an earlier run was void: the `MATIS` map was
   on `PETSC_COMM_SELF`, BDDC saw no interface and stopped at iteration 1).
   Repaired probe, control: relative residual 7.5e-09 / 5.8e-13 at rtol 1e-8 /
   1e-12, solution within 5.3e-10 / 3.7e-14 of the MUMPS solve at the same
   step. bp1002, 8 ranks, 101 steps: 27 / 41 iterations per solve (flat over
   all 202 solves), 3.17 / 3.71 s/step against 0.15-0.18 for MUMPS (17-25x
   slower), setup 1.0 s against 3.3 s, parity to MUMPS 5.6e-07 (rtol 1e-8) and
   1.0e-11 (1e-12) on the fault, station and cplot files, step count and
   peak-slip-rate step identical, no event inside 101 steps. 10 m BP8 wide,
   16 ranks, 3 steps, rtol 1e-8 (`runs/20261114_bddc-wide/`): BDDC 74
   iterations per solve, 339.7 s/step against 7.18 s/step for MUMPS (~47x
   slower), setup 60 s against 359.6 s, peak RSS of the largest rank 6.9 against
   11.9 GB; parity on wide not checked. Host not exclusive (load 7-22). So
   BDDC with these options does not beat MUMPS at either size: its gain is
   setup time and memory, its cost is ~0.06 s (bp1002) to ~2.3 s (wide, corrected
   2026-11-14: 149 iterations per step, two solves) per CG iteration; untested: tuned BDDC options (coarse space, local solvers,
   vertex/edge primal sets), more ranks.
   Mesh-only memory (stage 2 per-rank mesh), solver 3 on wide at 16 ranks:
   peak RSS after mesh setup 90.8 MB for the largest rank, 1.44 GB summed,
   against 11.9 GB per rank at the end: factors are >99 % of the footprint.
   **2026-11-14, fourth conductor** (MEASURED, indicative: host load 4-43 on
   every wide run; binary 1.24.4 of `dd-3g-spike2`, solver 3, no code change,
   runtime PETSc options only; no `par.` parameter).
   Where a step goes (`-eq_profile`, branch `dd-prof` eabe30a, spike, not for
   merge; `runs/20261114_prof/`): KSPSolve (the MUMPS forward/back solve) is
   87-95 % of the loop on bp1002 and on wide at 8/16/32 ranks; output 2-9 %,
   RHS vector assembly 1-2 %, friction + fault ghost exchange <= 2 % (imbalance
   up to 10x, still small), pore pressure ~0. Solve speed-up 8 -> 32 ranks is
   only 1.85x (wide MatSolve 6.06 / 4.34 / 3.27 s). BDDC, wide, 16 ranks: 2.29 s
   per CG iteration, of which local Neumann solve 1.07 s, local Dirichlet 0.88 s
   (MUMPS triangular solves), coarse solve 0.001 s, MatMult 0.08 s.
   MUMPS block-low-rank compression inside the PETSc direct solve
   (`-mat_mumps_icntl_35 1 -mat_mumps_cntl_7 <eps>`; `runs/20261114_blr-wide*`,
   wide, 3 steps, s/step / factorization s / peak RSS largest rank GB):
   | ranks | variant | s/step | factor s | GB |
   |---|---|---|---|---|
   | 16 | none | 8.32 | 344 | 12.2 |
   | 16 | BLR 1e-8 | 3.59 | 156 | 9.5 |
   | 16 | BLR 1e-6 | 2.69 | 129 | 9.7 |
   | 16 | BLR 1e-4 | 2.58 | 120 | 8.5 |
   | 32 | none | 7.55 | 242 | 6.9 |
   | 32 | BLR 1e-6 | 3.25 | 151 | 5.2 |
   | 32 | BLR 1e-4 | 2.47 | 135 | 5.0 |
   Ordering and threads do not help: METIS 8.11 s/step (= none), PORD 13.0 and
   26.5 GB, 8 ranks x 2 threads 12.8 s/step. Parity against the uncompressed
   run on the same binary: wide, 3 steps, worst relative difference over
   `cplot_EQquasi.txt` 1.4e-13 (1e-8), 4.6e-11 (1e-6), but 3 aseismic steps barely
   test it (the `fltst` log10 slip-rate column differs 2.1e-3 / 3.9e-3, its floor
   value is -20 against -30: a sentinel, excluded; audit 2026-11-14). `test.bp5.qdc.2000`, 8 ranks, 4483 steps to the stop tag (an event
   is included; `runs/20261114_blr-small/`): 1e-8 worst 4.5e-8, step count and
   time column identical; 1e-6 worst 5.7e-6 (`global.dat`), same step count and
   time column; 1e-4 worst 3.4e-3 and 4486 steps against 4483: REJECTED. bp1002
   101 steps: 1e-8 3e-11, 1e-6 4.6e-8, 1e-4 5.7e-6. So at the 5 s target: wide
   at 16 and at 32 ranks is below 5 s with BLR 1e-8 or 1e-6 at 16 ranks and 1e-6 at 32 (no 32-rank 1e-8 run; indicative, loaded
   host), the loop is still ~flat from 16 to 32 ranks. Not done: parity of BLR
   against the gold compsets through `-m e2e_fast` (needs a solver 3 gate case
   with these options), a wide run long enough to include slip-rate change, a
   host-exclusive repeat. Reaching the options in production needs the owner's
   answer on `par.solver = 3` (the options act on the PETSc loops only; native
   `par.solver = 1` has its own MUMPS controls).
   **BLR against the gold, 2026-11-14 09:31** (`runs/20261114_gate-blr/`, the
   `e2e_fast` reference tests unchanged, `par.solver = 3` forced in a throwaway
   branch, binary 1.24.4 of `dd-2d` 81629ad, 4-rank compsets, host load 3-12):
   exact solver 3: 39 passed, 21 skipped (data absence: no summary.json,
   section 4.3 profiles, offfault, perfault in the references), rc 0.
   BLR 1e-8: 21 failed / 18 passed, rc 1; BLR 1e-6: 32 failed / 7 passed,
   rc 1. Failures are at the gate's rtol 1e-9 in cplot, global, onfault,
   offfault, perfault, profile, restart and snapshot categories (e.g. a fltst column 1.115424e-05
   against 1.1154279e-05), i.e. the compressed solve differs from the gold at
   the level of its own tolerance, not at round-off. So BLR is NOT verified:
   it does not pass the gold comparator, and loosening that tolerance is not
   mine to do. Open, owner's call: what parity tolerance an approximate solve
   is held to (the row's own text: "within a stated tolerance, shown not to
   move event times"); the bp5 event run above shows event steps and times
   unmoved at 1e-8 and 1e-6.
   Long wide run (`runs/20261114_blr-widelong/`, 16 ranks, 200 steps = 1e5 s,
   slip rate rises from 1.3e-12 at step 3 to 1.479e-07 m/s; load before each
   run 3.7-5.1, own 16 ranks bring it to ~19, so near-quiet): none 5.52 s/step
   (setup 336 s), BLR 1e-8 2.87 (140 s), BLR 1e-6 2.48 (133 s); max slip rate
   at the end identical to 7 digits (1e-8) and 6 digits (1e-6, 1.47935e-07);
   worst difference over `cplot_EQquasi.txt` against none 3.3e-9 (1e-8), 6.8e-7
   (1e-6); `fltst` floor sentinel 1.1e-4 / 1.3e-4. Speedup 1.9-2.2x at 200 steps, smaller than the 2.3x / 3.1x (1e-8 / 1e-6) of 3 steps (the
   first steps' solves were cheaper to compress or the load differed).
   **BDDC tuning, 2026-11-14 fifth conductor** (MEASURED, indicative; spike
   binary 1.24.4 of `dd-3g-spike2`, runtime PETSc options only;
   `runs/20261114_bddc-tune/`, total 1-min load 6-13 before each bp1002 run).
   bp1002, 8 ranks, 11 steps = 22 solves, rtol 1e-8, exact MUMPS local and
   coarse solves, all converged on reason 2; MUMPS reference 0.15-0.18 s/step
   = 1.7-2.0 s per 11 steps:
   | variant | iterations per solve | time loop s | setup s |
   |---|---|---|---|
   | default | 27 | 29.5 | 1.03 |
   | + rigid-body near-null space (`-eq_near_null`) | 16 | 17.5 | 1.08 |
   | + faces, edges, vertices primal | 25 | 26.1 | 1.05 |
   | + deluxe scaling | 27 | 39.9 | 2.45 |
   | near-null + deluxe | 16 | 22.4 | 2.54 |
   | near-null + faces/edges/vertices | 15 | 15.7 | 1.09 |
   | default, 16 ranks | 37 | 23.2 | 0.45 |
   | default, 32 ranks | KSP reason -8 (indefinite PC) at step 1, run aborts | | |
   | near-null + adaptive threshold 5 | `PCSetUp_BDDC` error (`bddc.c:1528`), rc 56 | | |
   Best variant is ~8x slower than MUMPS on bp1002. Wide, 16 ranks, 3 steps,
   near-null + faces/edges/vertices: 31 iterations per solve (74 with default
   options), 125.9 s/step, setup 61.9 s, 2.02 s per iteration, against 339.7
   (default BDDC) and 7.18 (MUMPS), i.e. 17.5x slower. Why it cannot win with
   exact local solves: one iteration applies a local Neumann and a local
   Dirichlet triangular solve (~2 s on wide at 16 ranks) and one MUMPS solve of
   the whole system costs ~3.6 s, so BDDC would need fewer than two iterations
   per solve. Verdict: BDDC with direct local solves is measured negative for
   time per step at both sizes; its gains are setup time (62 s against 360 s)
   and memory. Remaining levers (inexact local solves, more ranks above the
   32-rank breakdown) are untested and unlikely to change this bound.
   Setup is built once and reused: `PCSetUp` 1, `KSPSetUp` 1, `PCApply` 151
   over 2 solves (`-log_view`, `runs/20261114_prof/bp8.qdc.gs.10.wide_s3_np16_bddc`;
   `KSPSetUp` is called once before the loop, `solveTimeLoopPETScDD.f90:680` at
   7b77221). Independent re-derivation of every figure above in this row
   against `runs/`, 2026-11-14: corrections folded in place (a mixed-up 95 %,
   17x as 17-25x, 362 s as 359.6 s, no 32-rank BLR 1e-8 run, parity metric
   named, load ranges, 3x as 2.3x / 3.1x, gate categories and time).
   **`KSPConvergedReason = -11` (`PCSETUP_FAILED`) at the first solve, 8-rank
   BP8**: 1 of 11 runs reported, 0 of 40 in `runs/20261114_ksp-stress/`. Not
   reproduced. Closes when `KSPSetUp` reads the MUMPS INFOG(1)/(2) on failure
   (the native loop does, `solveTimeLoopMUMPS.f90:60,153`) so a recurrence names
   its cause, or on 60 clean runs; until then it stays an open observation.
   **Slip-to-traction operator, sizing only (2026-11-14, arithmetic from the
   run logs, NOT a measurement)**: fault elements 7,440 (bp1002) and 53,868
   (10 m BP8 wide), so about 15 k and 108 k slip unknowns (two components per
   fault node). Dense operator 15 k^2 x 8 B = 1.8 GB (bp1002), 108 k^2 x 8 B =
   93 GB (wide; ~6 GB per rank at 16 ranks, half that symmetric or in float32).
   Build cost is one volume solve per unknown: 108 k x ~3.6 s (a MUMPS solve on
   wide, 16 ranks) = ~108 h one rank-group, unless batched (`MatMatSolve`),
   whose per-right-hand-side throughput is the one number not yet measured.
   Per step: one matvec streaming the matrix, ~0.6 s at an assumed 160 GB/s
   against 7.18 s now, i.e. ~10x if the build is paid. bp1002 gains ~5x at
   most (0.03 s matvec against 0.15 s). Feasible in memory, open in build time;
   **Correction, 2026-11-14 (independent audit; `src/meshgen.f90:637`)**: "Fault
   elements (et=2)" are the volume hexes within 2*dy of the fault, not fault
   quads, so 7,440 and 53,868 are not fault sizes. Slip unknowns are two per
   fault node: wide has 10,201 split nodes (6,561 "fault nodes in Omega_f"),
   i.e. ~20.4 k (13.1 k) unknowns, bp1002 594 nodes, ~1.2 k. Re-derived (arithmetic,
   not measured): dense operator 3.3 GB on wide (not 93 GB), build at 0.231 s per
   right-hand side ~4,700 s (not 6.9 h), break-even ~650 steps (not 3,800). The
   15 k / 108 k / 93 GB / 6.9 h / 3,800-step figures in this paragraph and the
   next rest on the element count and are void. The 7.18 s/step baseline has no
   run behind it (measured wide MUMPS at 16 ranks: 6.9-7.0 s/step).
   **Measured 2026-11-14** (`runs/20261114_matmat/`, spike branch `spike-mm`, not
   merge-bound, binary 1.24.4 plus a timer; PETSc `MatMatSolve` on the MUMPS
   Cholesky factor, dense random right-hand sides, `-eq_spike_log -eq_matmat 32`,
   host load 15 after launch, indicative). Wide, 16 ranks: 1 column 4.41 s, 8
   columns 3.41 s (0.426 s each), 32 columns 7.40 s (**0.231 s each**); the
   1-column call is the first solve and includes warm-up. bp1002, 8 ranks: 0.0636
   s, 0.0127 s each (8), 0.0083 s each (32). So the build for wide is 108 k x
   0.231 s = ~25,000 s = ~6.9 h at 16 ranks (bp1002: 15 k x 0.0083 = ~2 min),
   not the ~108 h of the arithmetic above; sparse fault-only right-hand sides
   and more columns per call are untested and may lower it. Break-even against
   7.18 s/step with a ~0.6 s matvec (assumed, not measured): ~3,800 steps. The
   matvec, the dense memory (93 GB) and whether the operator is valid under the
   moving-load terms are still not measured.


   **Plain Krylov, subdomain-direct, exact substructuring (2026-11-14;
   `runs/20261114_krylov*`, spike branch `spike-mm` ddd36b2, not merge-bound,
   indicative: host load 13-21 during runs).** Owner lines: plain CG/GMRES; direct
   solves per subdomain under an outer iteration; factorize once. MEASURED
   (same-script re-run), not VERIFIED. K is symmetric to round-off (`-eq_check_sym`).
   - bp1002, 8 ranks: CG+GAMG with rigid-body near-null space is the best plain
     Krylov, 2.6-3.1x slower than MUMPS (loop 4.4-5.2 s against 1.7 s), 0.32 GB against 0.50 GB peak, setup 0.97 s
     against 3.23 s. Jacobi and block-Jacobi+ICC are slower still. CG breaks
     down only in the BDDC inexact variants (KSP -8), where GMRES with right
     preconditioning was used instead.
   - GAMG+near-null+previous-step guess, bp1002 11 steps (22 solves; "3 steps" was wrong): 10.2 iterations per
     solve (max 11), relative residual 9e-9, at 16 ranks (loop 2.46 s, setup
     0.54 s, peak RSS 200 MB per rank; unit confirmed: `rss.*` is GNU `time -f %M`, kB, and agrees with the in-code line, 478,760 kB against 470.5 MB, 2026-11-14) and 32 ranks (loop 1.33
     s, setup 0.34 s, 140 MB). Worst parity against exact at 16 and 32 ranks 1.9e-8
     (netCDF `disp` and `fault` fields, relative to the field maximum; 5e-11 was
     wrong), 2.5e-10 over the text files (`srfst_*`).
   - **Wide, 16 ranks, 3 steps:** exact MUMPS setup 337 s, loop 20.8 s (6.9
     s/step), peak RSS 12.1 GB per rank; GAMG CG setup 10.6 s, loop 97.4 s
     (32 s/step, 4.7x slower), 63 iterations per solve (max 67), 2.5 GB per rank
     (4.8x less), relative residual 9e-9; worst-case parity against exact over the 33 compared files of
     `result/cycle0` (per file max|a-b| / max|a|; log10 slip-rate column of the `fltst`
     files excluded, floor sentinel; 3 aseismic steps): text and `.dat` files
     <= 4.1e-15, `fltst` <= 6.7e-19, netCDF fields 6.6e-9 to 1.1e-7 (`fault.00001.nc`
     `shear_dip`, `disp.00001.nc` `disp`; the worst `fault` file is 1.75e-8, `fault.00003.nc` and `fault.r.nc`), i.e. 1e-7 against the gold gate's 1e-9, so
     GAMG CG at rtol 1e-8 would not pass the comparator either.
     Conclusion: the iterative route trades 4.7x time for 4.8x memory and a
     30x shorter setup; it does not beat MUMPS on time at 3.5 M equations.
   - Subdomain-direct (CG, block Jacobi / additive Schwarz, exact MUMPS blocks,
     bp1002 8 ranks): factorized once (`-log_view`: setup 1, apply per
     iteration); a single block pass leaves residual 3.6e-2 / 3.1e-1 / 7.1e-2
     (per variant); one-level variants need a coarse space (iteration counts grow
     with block count). Not measured: wide, and a slip-rate-change case (the
     11-step bp1002 runs are aseismic).
   - Inexact BDDC local solves (GMRES, BLR 1e-4, ICC, ILU, GAMG): none beats
     exact local solves; BLR 1e-4 15 iterations in 22 s, ICC/ILU 69 iterations in
     16-18 s; np32 fails (-8) and adaptive threshold errors at `bddc.c:1528`.
   - **Exact 3-step substructuring, MEASURED 2026-11-14** (spike branch
     `spike-sub` 7d469dc, one flag, not merge-bound; host load 6-14, indicative).
     Runtime route: `-eq_mat_is -eq_all_primal` makes every interface dof a BDDC
     primal vertex (`PCBDDCSetPrimalVerticesLocalIS`), so the coarse problem is the
     whole interface Schur complement, factorized once by MUMPS, with exact local
     solves (MUMPS) and CG: 1 iteration per solve, relative residual 1.2-1.3e-15
     in all 8 runs. bp1002, 11 steps, exact solver 3 beside it (same binary 1.24.4,
     `runs/20261114_substruct/`); setup s / time loop s / peak RSS GB per rank:
     | ranks | interface dofs | exact | substructuring | worst diff from exact |
     |---|---|---|---|---|
     | 4 | 6,600 | 5.40 / 2.51 / 0.73 | 19.9 / 7.19 / 2.29 | 1.3e-15 |
     | 8 | not logged (<= 3,492 per rank) | 3.07 / 1.43 / 0.48 | 10.9 / 2.89 / 1.98 | 1.3e-15 |
     | 16 | 17,793 | 2.22 / 0.97 / 0.33 | 6.20 / 2.12 / 1.09 | 2.0e-19 |
     | 32 | 26,433 | 1.84 / 1.02 / 0.20 | 3.87 / 1.76 / 0.76 | 7.5e-16 |
     (worst diff: per-file max|a-b| / max|a| over `result/cycle0`; identical digits at 4 and 8
     ranks are the same last-printed-digit flip of `global.dat`, the two files differ.)
     So it costs 2-3.7x setup, 1.7-2.9x time loop and 3.1-4.1x memory against exact
     MUMPS on bp1002, and it is exact. **Gold gate** (`e2e_fast` reference tests,
     comparator unchanged, `par.solver = 3` forced in throwaway branch `gate-sub`
     b57a6a1, `runs/20261114_gate-sub/`, `runInfo.json` shows the spike path ran): 39
     passed, 21 skipped (data absence, six reasons: no summary.json, section 4.3
     profiles, offfault, perfault, profile, restart files), rc 0, 894 s: it passes at
     1e-9, where BLR failed 21 and 32 and exact solver 3 has 0 failures. Only 4-rank
     compsets. Wide not run: interface ~183-192 k dofs (audit estimates, not the 170 k
     above), dense Schur 270-295 GB plus its factor, and the build needs one local
     solve per interface dof (~170 k x 0.23 s ~ 11 h at 16 ranks); not an
     unattended run. Mid-size feasibility (`bp8.qdc.gs.10`, 16 ranks) is queued: see
     the session log.
   - **Substructuring on an event-resolving run (bp5, 4483 steps, 8 ranks), MEASURED 2026-11-14** (`runs/20261114_substruct-bp5/`, spike binary, same script): exact MUMPS setup 3.4 s, loop 577 s, kspsolve 516 s; `-eq_mat_is -eq_all_primal` setup 8.0 s, loop 1202 s (2.1x slower), kspsolve 1139 s. Per-file worst max|a-b|/max|a| over the 120 compared files, no variable excluded: 3.9e-11 (`fault.03001.nc` `slipn`; 4.6e-14 under the consolidated table's metric, which drops variables with max < 1e-16), so the iterative-in-name solve still carries ~1e-11 drift over 4483 steps, far below the 1e-9 comparator. Peak RSS (GNU time, largest of the 8 ranks, kB): exact 520,616, all-primal 1,236,568 (2.4x; corrected 2026-11-14: the first-posted 491,348 and 1,091,204, ratio 2.2x, were rank 7's, not the maximum). GAMG CG on bp5 (harvested 2026-11-14 12:20, same dir, same script): setup 1.0 s, loop 1500.8 s (2.6x exact at load 5.8, 1.9x exact at load 25.7, 799 s; GAMG ran at load 19.2), 6.8 KSP iterations per solve (max 12), peak RSS 0.35 GB per rank (1.5x less than exact), worst per-file difference 1.5e-7 (`fault.02601.nc` `vzm`, float32 floor), 116 of 120 files fail the comparator: inexact, as on bp1002 and wide. The subdomain-direct runs are still queued (sub4). MEASURED (same-script re-run), not VERIFIED.
   - **gs10 (`bp8.qdc.gs.10`, 1.96 M equations), 16 ranks, 3 steps, exact vs
     all-primal substructuring, MEASURED 2026-11-14 12:00** (`runs/20261114_substruct-gs10/`,
     spike binary, load 11-12 before each run, both rc 0, `DONE` present): setup 80.8 s
     against 268.4 s (3.3x), loop 8.1 s against 13.3 s (2.71 against 4.45 s/step, 1.6x), peak RSS
     of the largest rank 4.47 against 20.68 GB (4.6x), summed over ranks (`rss.*`) 63.3 against
     261.2 GB (4.1x); 1 CG iteration per solve; interface N_Gamma = 76,992 dofs (sum of
     `EQ_PRIMAL` owned-shared over ranks), so a dense N_Gamma x N_Gamma matrix is 47 GB.
     The same-comparator check below passes 33 of 33 files (worst per-file difference
     4.3e-10). Substructuring is feasible at 2 M equations and loses on every column.
   - **Why all-primal substructuring costs more memory than global MUMPS** (diagnostic,
     `runs/20261114_substruct-mem/`, bp1002, 8 ranks, 3 steps, MUMPS `ICNTL(4)=2` on the
     BDDC solvers; host load 13-14, memory is not load-sensitive; MEASURED, partial
     attribution). Exact: one global factor, 2,914 MB in total (INFOG(17)) of 3,432 MB summed
     peak RSS (`rss.*`; the in-code `Peak RSS` line reads 3,377, 13,314, 1,931 and 489.5 MB, a different counter). All-primal, same case: 16 local factorizations, two per rank (read as the Neumann and
     Dirichlet solvers, each of local order 18.1-19.0 k (MUMPS `JOB, N, NNZ` lines of `run.log`: 18,135 to 18,972), 89-100 MB each, 1.55 GB in total),
     plus one coarse factorization of the whole interface Schur complement, order 12,924,
     3,412 MB in total (645 MB estimated on the largest rank): the coarse factor alone is
     larger than the entire global factor, and the local factors come on top of it. Dense
     local Schur blocks (sum of n_loc^2 x 8 B) add 0.7 GB here and 12.8 GB on gs10. Summed
     peak RSS is 13,634 MB (largest rank 1,977 MB against 498 MB); the factors above are
     about 5.0 GB, the dense local Schur 0.7 GB, and ~8.0 GB (58 %) is not attributed (MUMPS
     work arrays above the factor size, coarse-matrix assembly copies; not measured). Why it
     is structural (a reading; the control that would test it, wide-box Schur order against global fill, was not run): MUMPS' elimination tree already is a multilevel substructuring, with
     the interfaces eliminated sparsely as separators; the all-primal route flattens it
     into one level, factorizes the interior once per subdomain twice over and the
     interface as one dense-blocked problem whose order grows with rank count (bp1002:
     12.9 k at 8 ranks, 17.8 k at 16, 26.4 k at 32; memory ratio all-primal / exact
     3.1, 4.1, 3.3, 3.8x at 4, 8, 16, 32 ranks).
   - **Does any variant change the verdict at these sizes without new merge-bound code?
     No.** What the spike binary can reach by runtime options: fewer primal dofs (BDDC
     default, near-null, faces/edges/vertices): memory 0.40-0.41 GB against 0.50 on bp1002
     and 6.9 against 11.9 GB on wide, but 9-49x slower and not exact; GAMG CG: 4.8x less
     memory on wide, 2.6-4.7x slower, not exact; MUMPS BLR on the global solve: 2-3x faster,
     ~20 % less memory, fails the gold comparator; subdomain-direct under CG: 16-19x slower
     on bp1002 (table). BLR inside the all-primal local and coarse solves was not run; it
     would make the solve inexact, the property that is the only thing substructuring
     buys over GAMG or BDDC. What would change the memory side is a multilevel or
     compressed Schur factorization, which is new code (and re-implements what MUMPS does).
     Plainly: at 0.2 to 2 M equations and 4-32 ranks, exact substructuring is 1.6-2.9x
     slower in the loop, 2-3.7x slower in setup and 2.4-4.6x larger than exact global
     MUMPS, and is exact; wide is out of reach. Among exact options (passing the
     comparator against the exact run) MUMPS is the fastest and smallest; among inexact
     ones only BLR is faster than MUMPS, and it does not pass the gold comparator.
   - **Consolidated comparison, 2026-11-14 12:15** (MEASURED, same-script re-runs, not
     VERIFIED; written by script from `runInfo.json`, `rss.*` and `result/cycle0`, csv
     parse-checked; each row cites its run directory, "not run" means no run exists;
     the rows tagged running fill in when their `DONE` lands). Definitions: setup =
     `factorization_seconds`, loop = `time_loop_seconds` of `runInfo.json` (rows
     elsewhere in this row quote the in-code `EQ_WALL` line, up to 5 % different, 3.5-4.6 % on the 32-rank bp1002 loops); peak
     RSS = largest `rss.*` over ranks (GNU time `%M`, kB / 1e6); worst diff = per-file
     max|a-b| / max|a| against the exact run of the same family over text, `.dat` and
     netCDF files, log10 slip-rate columns of `fltst` and netCDF variables with max
     < 1e-16 (symmetry-noise dip components) excluded; netCDF is float32, so 6e-8 to
     1.2e-7 is one to two ulps (GAMG and both BDDC rows on wide read 1.1e-7 alike for
     that reason: a quantization floor, not a solver error figure); the earlier text-only
     figures in this row (BLR 1.4e-13, 4.6e-11 on `cplot`) are a different, weaker metric.
     Comparator vs exact = `compare_series` / `compare_netcdf` of `tests/e2e/cases.py` at
     1e-9 between the run and its exact run (all files except mesh and json). Gold
     comparator column = the `e2e_fast` reference tests with `par.solver = 3` forced, per
     option on the 4-rank compsets, not on the row's own case (`runs/20261114_gate-blr/`,
     `runs/20261114_gate-sub/`). Load = 1-minute load before the run: bp5 exact is 799 s
     at load 25.7 and 573 s at load 5.8, so timings compare within a case family only.
     Subdomain-direct setup reads ~0: the block factorizations run inside the first
     solve and are counted in the loop.

     | case | option | run dir (under runs/20261114_) | setup s | loop s | peak RSS GB per rank | KSP its/solve | worst diff from exact | comparator vs exact run | gold comparator (e2e_fast, solver 3) | load before |
     |---|---|---|---|---|---|---|---|---|---|---|
     | bp1002, 8 ranks, 11 steps | exact MUMPS | krylov/test.bp1002.qdc.zone.2500_s3_np8_exact | 3.23 | 1.7 (0.153/step, 11 st) | 0.50 | 1.0 | reference | reference | yes: 39 passed, gate-blr | 12.02 |
     | bp1002, 8 ranks, 11 steps | BLR 1e-8, 1e-6 | not run: at 11 steps; see the 101-step case | not run | not run | not run | not run | not run | not run | not run | not run |
     | bp1002, 8 ranks, 11 steps | all-primal substructuring | substruct/test.bp1002.qdc.zone.2500_s3_np8_allprimal | 10.90 | 2.9 (0.261/step, 11 st) | 1.98 | 1.0 | 1.5e-19 | 25/25 files pass | yes: 39 passed, gate-sub | 6.88 |
     | bp1002, 8 ranks, 11 steps | GAMG CG | krylov/test.bp1002.qdc.zone.2500_s3_np8_gamgg | 0.97 | 4.4 (0.401/step, 11 st) | 0.32 | 10.5 | 1.9e-08 | 12 of 25 files fail | not run | 13.27 |
     | bp1002, 8 ranks, 11 steps | subdomain-direct CG, 16 blocks | krylov2/test.bp1002.qdc.zone.2500_s3_np8_bj16 | 0.00 | 26.8 (2.435/step, 11 st) | 0.14 | 144.4 | 3.7e-08 | 15 of 25 files fail | not run | 13.45 |
     | bp1002, 8 ranks, 11 steps | subdomain-direct CG, 16 blocks (rerun) | subdirect/test.bp1002.qdc.zone.2500_s3_np8_bj16 | 0.01 | 31.3 (2.843/step, 11 st) | 0.23 | 66.5 | 1.9e-08 | 14 of 25 files fail | not run | 13.07 |
     | bp1002, 8 ranks, 11 steps | subdomain-direct CG, 64 blocks | subdirect/test.bp1002.qdc.zone.2500_s3_np8_bj64 | 0.01 | 27.5 (2.500/step, 11 st) | 0.17 | 107.7 | 7.5e-08 | 14 of 25 files fail | not run | 11.61 |
     | bp1002, 16 ranks, 11 steps | subdomain-direct CG, 1 block per rank | subdirect/test.bp1002.qdc.zone.2500_s3_np16_bj16r | 0.00 | 15.7 (1.425/step, 11 st) | 0.14 | 66.0 | 1.9e-08 | 14 of 25 files fail | not run | 11.77 |
     | bp1002, 8 ranks, 11 steps | BDDC default | bddc-tune/test.bp1002.qdc.zone.2500_s3_np8_base | 1.03 | 29.4 (2.675/step, 11 st) | 0.40 | 27.0 | 1.5e-08 | 13 of 25 files fail | not run | 2.82 |
     | bp1002, 8 ranks, 11 steps | BDDC near-null + faces/edges/vertices | bddc-tune/test.bp1002.qdc.zone.2500_s3_np8_nn_fev | 1.09 | 15.7 (1.428/step, 11 st) | 0.41 | 15.0 | 9.3e-09 | 12 of 25 files fail | not run | 10.83 |
     | bp1002, 8 ranks, 101 steps | exact MUMPS | blr-small/test.bp1002.qdc.zone.2500_s3_np8_base | 3.43 | 19.2 (0.191/step, 101 st) | 0.50 | 1.0 | reference | reference | yes: 39 passed, gate-blr | 20.05 |
     | bp1002, 8 ranks, 101 steps | BLR 1e-8 | blr-small/test.bp1002.qdc.zone.2500_s3_np8_blr8 | 3.76 | 20.5 (0.203/step, 101 st) | 0.48 | 1.0 | 1.5e-08 | 13 of 25 files fail | no: 21 failed, gate-blr | 22.42 |
     | bp1002, 8 ranks, 101 steps | BLR 1e-6 | blr-small/test.bp1002.qdc.zone.2500_s3_np8_blr6 | 3.06 | 13.7 (0.136/step, 101 st) | 0.47 | 1.0 | 1.0e-06 | 21 of 25 files fail | no: 32 failed, gate-blr | 24.29 |
     | bp1002, 8 ranks, 101 steps | BLR 1e-4 (rejected) | blr-small/test.bp1002.qdc.zone.2500_s3_np8_blr4 | 2.69 | 9.4 (0.093/step, 101 st) | 0.43 | 1.0 | 2.0e-04 | 23 of 25 files fail | not run | 25.13 |
     | bp1002, 8 ranks, 101 steps | all-primal, GAMG CG, subdomain-direct, BDDC | not run: not run at 101 steps | not run | not run | not run | not run | not run | not run | not run | not run |
     | bp5 event run, 8 ranks, 4483 steps | exact MUMPS (BLR family, load 25.7) | blr-small/test.bp5.qdc.2000_s3_np8_base | 3.72 | 799.4 (0.178/step, 4483 st) | 0.55 | 1.0 | reference | reference | yes: 39 passed, gate-blr | 25.69 |
     | bp5 event run, 8 ranks, 4483 steps | BLR 1e-8 | blr-small/test.bp5.qdc.2000_s3_np8_blr8 | 3.69 | 850.9 (0.190/step, 4483 st) | 0.55 | 1.0 | 1.2e-07 | 115 of 120 files fail | no: 21 failed, gate-blr | 22.20 |
     | bp5 event run, 8 ranks, 4483 steps | BLR 1e-6 | blr-small/test.bp5.qdc.2000_s3_np8_blr6 | 4.10 | 786.3 (0.175/step, 4483 st) | 0.54 | 1.0 | 4.6e-05 | 118 of 120 files fail | no: 32 failed, gate-blr | 25.49 |
     | bp5 event run, 8 ranks, 4483 steps | BLR 1e-4 (rejected: 4486 steps) | blr-small/test.bp5.qdc.2000_s3_np8_blr4 | 2.91 | 380.8 (0.085/step, 4486 st) | 0.52 | 1.0 | n/a (MISSING disp.04483.nc) | not comparable (4486 steps vs reference) | not run | 13.37 |
     | bp5 event run, 8 ranks, 4483 steps | exact MUMPS (substructuring family, load 5.8) | substruct-bp5/test.bp5.qdc.2000_s3_np8_exact | 3.44 | 573.0 (0.128/step, 4483 st) | 0.52 | 1.0 | reference | reference | yes: 39 passed, gate-blr | 5.82 |
     | bp5 event run, 8 ranks, 4483 steps | all-primal substructuring | substruct-bp5/test.bp5.qdc.2000_s3_np8_allprimal | 8.02 | 1198.2 (0.267/step, 4483 st) | 1.24 | 1.0 | 4.6e-14 | 120/120 files pass | yes: 39 passed, gate-sub | 17.43 |
     | bp5 event run, 8 ranks, 4483 steps | GAMG CG | substruct-bp5/test.bp5.qdc.2000_s3_np8_gamgg | 1.04 | 1500.8 (0.335/step, 4483 st) | 0.35 | 6.8 | 1.5e-07 | 116 of 120 files fail | not run | 19.22 |
     | bp5 event run, 8 ranks, 4483 steps | subdomain-direct, BDDC | not run: not run on bp5 | not run | not run | not run | not run | not run | not run | not run | not run |
     | gs10 (bp8.qdc.gs.10), 16 ranks, 3 steps | exact MUMPS | substruct-gs10/bp8.qdc.gs.10_s3_np16_exact | 80.78 | 8.1 (2.708/step, 3 st) | 4.47 | 1.0 | reference | reference | yes: 39 passed, gate-blr | 12.32 |
     | gs10 (bp8.qdc.gs.10), 16 ranks, 3 steps | all-primal substructuring | substruct-gs10/bp8.qdc.gs.10_s3_np16_allprimal | 268.41 | 13.3 (4.450/step, 3 st) | 20.68 | 1.0 | 4.3e-10 | 33/33 files pass | yes: 39 passed, gate-sub | 11.35 |
     | gs10 (bp8.qdc.gs.10), 16 ranks, 3 steps | BLR, GAMG CG, subdomain-direct, BDDC | not run: not run on gs10 | not run | not run | not run | not run | not run | not run | not run | not run |
     | wide (10 m BP8), 16 ranks, 3 steps | exact MUMPS | krylov4/bp8.qdc.gs.10.wide_s3_np16_exact | 331.26 | 20.6 (6.881/step, 3 st) | 12.07 | 1.0 | reference | reference | yes: 39 passed, gate-blr | 13.27 |
     | wide (10 m BP8), 16 ranks, 3 steps | BLR 1e-8 | blr-wide/bp8.qdc.gs.10.wide_s3_np16_blr8 | 156.38 | 10.8 (3.585/step, 3 st) | 9.52 | 1.0 | 4.4e-07 | 13 of 33 files fail | no: 21 failed, gate-blr | 16.74 |
     | wide (10 m BP8), 16 ranks, 3 steps | BLR 1e-6 | blr-wide/bp8.qdc.gs.10.wide_s3_np16_blr6 | 128.59 | 8.1 (2.690/step, 3 st) | 9.66 | 1.0 | 1.6e-04 | 15 of 33 files fail | no: 32 failed, gate-blr | 17.02 |
     | wide (10 m BP8), 16 ranks, 3 steps | all-primal substructuring | not run: not run: dense Schur 270-295 GB, ~11 h build (sizing) | not run | not run | not run | not run | not run | not run | not run | not run |
     | wide (10 m BP8), 16 ranks, 3 steps | GAMG CG | krylov4/bp8.qdc.gs.10.wide_s3_np16_gamgg | 10.61 | 97.2 (32.416/step, 3 st) | 2.54 | 63.0 | 1.1e-07 | 9 of 33 files fail | not run | 21.47 |
     | wide (10 m BP8), 16 ranks, 3 steps | subdomain-direct CG, 1 block per rank | subdirect/bp8.qdc.gs.10.wide_s3_np16_bj_wide | 0.00 | 1000.6 (333.5/step, 3 st) | 3.36 | 229.3 | 1.1e-07 | 9 of 33 files fail | not run | 13.53 |
     | wide (10 m BP8), 16 ranks, 3 steps | BDDC default | bddc-wide/bp8.qdc.gs.10.wide_s3_np16_bddc | 60.19 | 1019.2 (339.720/step, 3 st) | 6.95 | 74.2 | 1.1e-07 | 9 of 33 files fail | not run | 20.86 |
     | wide (10 m BP8), 16 ranks, 3 steps | BDDC near-null + faces/edges/vertices | bddc-tune/bp8.qdc.gs.10.wide_s3_np16_nnfev | 61.91 | 377.6 (125.883/step, 3 st) | 6.97 | 31.0 | 1.1e-07 | 9 of 33 files fail | not run | 9.82 |
     | wide (10 m BP8), 16 ranks, 200 steps | exact MUMPS | blr-widelong/bp8.qdc.gs.10.wide_s3_np16_base | 335.57 | 1103.4 (5.517/step, 200 st) | 11.66 | 1.0 | reference | reference | yes: 39 passed, gate-blr | 4.17 |
     | wide (10 m BP8), 16 ranks, 200 steps | BLR 1e-8 | blr-widelong/bp8.qdc.gs.10.wide_s3_np16_blr8 | 140.30 | 573.9 (2.869/step, 200 st) | 9.53 | 1.0 | 5.4e-07 | 15 of 33 files fail | no: 21 failed, gate-blr | 3.70 |
     | wide (10 m BP8), 16 ranks, 200 steps | BLR 1e-6 | blr-widelong/bp8.qdc.gs.10.wide_s3_np16_blr6 | 132.91 | 495.6 (2.478/step, 200 st) | 9.25 | 1.0 | 2.0e-04 | 19 of 33 files fail | no: 32 failed, gate-blr | 5.10 |
     | wide (10 m BP8), 16 ranks, 200 steps | all-primal, GAMG CG, subdomain-direct, BDDC | not run: not run at 200 steps | not run | not run | not run | not run | not run | not run | not run | not run |

   - Open (2026-11-14 12:15): `-11` 60-run batch running, 2 of 60 clean (`runs/20261114_ksp-stress2/`, binary 1.24.4,
     ~8-10 min per run behind the heavy lock; `rc` file, closes on `ALLDONE`); sub4 and the bp5 GAMG CG run are harvested (table rows; 12:35). Wide subdomain-direct, 1 block per rank (`runs/20261114_subdirect/`, MEASURED, host load 13.5, indicative): setup ~0, loop 1000.6 s against 20.6 s exact (48x slower), 229 CG iterations per solve (max 247), peak RSS 3.36 against 12.07 GB per rank (3.6x less), worst difference 1.1e-7 (`disp.00001.nc`, float32 floor), 9 of 33 files fail: it loses on time by a wider margin than GAMG CG (97 s) and is inexact. bp1002 at 16 ranks: loop 15.7 s against 0.94 s exact (16.7x), 66 iterations, 14 of 25 fail.
3. [ ] **P1 -- EQquasi BP8 discrepancy (bug hunt).** Owner, 2026-10-09: "hunt
   the bug". An independent FEniCSx linear-hex model with EQquasi's box and BCs
   converges to Kim within ~2 % at 25 m (`docs/dev/fenicsx_feasibility_2026-10-09.md`);
   EQquasi does not (depth station -18 % at 50 m, -13 % at 10 m). Pressure
   matches; domain and dtmax ruled out. The uploaded BP8 entries are suspect
   along depth until this closes. Done when: the cause is found and fixed, and
   EQquasi at 25 m agrees with Kim within a few % at all 9 stations.
   **Escalation, 2026-10-12**: 6 code-audit findings on this hunt
   rejected after a source/gate check (4 before today: nodal area/dx^2, GS
   z-symmetry, and two more per the 2026-10-09 brief; 2 today, PRs #58 and
   #59, closed unmerged -- see `docs/dev/SESSION_LOG_2026-10-12_board-sweep.md`
   for each). Today's two both targeted `src/faulting.f90`'s split-node
   traction formula (lines 112-128): #58 claimed "ruled out" without ever
   computing the required EQquasi-vs-FEniCSx ratio; #59 computed only the
   FEniCSx half and flagged the mass-weighted cross term as a "structural
   bug," which does not survive the rule-16a sanity check -- that formula
   sits outside every `bp ==` branch (the only ones in the file are at
   lines 159/162/359, pore-pressure/BP7-specific), so a bias large enough to
   explain BP8's -13%/-18% would also break BP5/BP1002's passing regression
   gold, and it doesn't. The split-node coupling itself is very likely not
   the bug (3 independent rejected attempts have now landed on it from
   different angles). **Needs owner diagnosis, not a 3rd same-session
   retry**: live lead is still unconfirmed-but-unrefuted (tau/pressure match
   Kim within ~0.5% but V amplifies a small tau gap exponentially); next
   untried candidate is `src/qdct*.f90`'s element-to-fault force
   contribution (assembly path, not yet audited by anyone).
4. [ ] **P2 -- FEniCSx backend, EQdyna pattern.** Owner, 2026-10-09: "like
   eqdyna, we will have a fortran and fenicsx backend now ... follow eqdyna's
   success". One compset and one gold set, a test matrix that runs each case
   on each backend, a backend counts only at parity with the gold. Order:
   after row 3 (the Fortran gold must be right); BP5 and step-over parity
   first, BP8 against Kim. The prototype's precomputed slip-to-traction matrix
   (quasi-static, linear) is also a scaling idea for row 1.
6. [ ] **P2 -- Time-step factor xi above the stability limit.** Owner,
   2026-11-13: "that ksi lesson, put on board". Lesson: a constant `par.xi`
   above the Lapusta et al. (2000) eq. 15 limit makes the time step
   oscillate between events. Reported by the session that ran the
   comparison, unaudited here: Wang 2026 Case 3 at 1 km with xi 0.2
   (`runs/20261106_wang2026_case3/{c3a,c3b}`, `result/cycle*/global.dat`),
   90-97 % of interseismic steps reverse the sign of the step-to-step change
   in log10 Vmax; with per-node xi (PR #72, `par.xi_gamma = pi/4`) 0-1 %.
   The eq. 15 limit with k = gamma mu / h, a 0.005, b 0.01, L 0.05 m,
   sigma 25 MPa, 36 GPa is 0.097 at 1 km and 0.046 at 500 m, so xi 0.2 at
   1 km is about twice too large. Event times at 1 km moved by at most
   0.8 yr over five events; the defect is between events, not in the cycles.
   Open, owner's call: (a) `case.setup` computes the eq. 15 limit from
   material, h and friction and warns or stops when `par.xi` exceeds it;
   (b) per-node xi as the default with `par.xi` as the cap; (c) PR #72's
   open audit items -- `xi_gamma = NaN` silently off, the stop-13 check also
   runs at creeping nodes, stop 13 exits without `MPI_Abort`.
   Done when: PR #72 is merged or closed with (c) resolved, and (a) and (b)
   are each decided and, if yes, landed with a test.
