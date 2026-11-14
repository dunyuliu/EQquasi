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
   next measurable: `MatMatSolve` seconds per right-hand side on wide at 16
   ranks, 32 columns.

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
