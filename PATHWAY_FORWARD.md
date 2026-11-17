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
   **ROW 1 SUMMARY FOR THE OWNER, 2026-11-14 (conductor 12; every figure is a table or paragraph below, none new; 8-16 ranks, MEASURED, hosts shared so loads vary run to run).**
   - *Passes the gold-style comparator (1e-9):* exact MUMPS (reference); all-primal substructuring (bp1002 25/25, bp5 120/120, gs10 33/33; wide not runnable, dense Schur 270-295 GB); GAMG CG at tight rtol (wide 1e-14: 0 of 33 fail; bp1002 1e-11/1e-12); BLR-preconditioned CG at tight rtol. *Fails it:* BLR 1e-8/1e-6/1e-4 as a solver (bp1002 13-23 of 25, bp5 115-118 of 120, wide 13-19 of 33), GAMG CG at rtol 1e-8 (wide 9 of 33, bp5 116 of 120), subdomain-direct CG and BDDC (9-15 of 25/33).
   - *Seconds per step (loop):* bp5 0.128-0.178 exact, 0.267 all-primal, 0.335 GAMG; gs10 16 ranks 2.71 exact, 4.45 all-primal; wide 16 ranks 6.9-8.3 exact (5.5 at 200 steps), 32.4 GAMG at rtol 1e-8, 49.9-53.4 GAMG at the passing rtol (6-8x slower), 2.9 BLR 1e-8 (1.9x faster, fails). Exact wide at 8/16/32 ranks: 12.7 / 6.9 / 7.1; the row's target (wide <= 5 s/step at 32 ranks) is not met by anything that passes.
   - *Setup and memory per rank (wide, 16 ranks):* exact 331-344 s and 12.1-12.3 GB; GAMG CG 10 s and 2.54 GB (4.8x less); BLR 128-156 s and 9.5 GB; all-primal on gs10 268 s and 20.7 GB against exact 81 s and 4.5 GB.
   - *Merge gate:* full `e2e` tier on the stack head, solver 3 forced: 90 passed, 45 skipped (10 pending bp1002 zone gold, 35 category-absent), 0 failed (`runs/20261114_gate-full/`). Audit of the stack diff: 0 BLOCKER, 0 MAJOR, 9 MINOR (all listed above). Stack #76, #77, #82, #83, #85 is MERGEABLE with CI green; #87 (BDDC spike) is not for merge.
   - *64-rank `-11` abort:* MUMPS -9, cured by `-mat_mumps_icntl_14 100` (3 of 3 pass against 5 of 6 fail; +20-30 % rank memory), passed through `PETSC_OPTIONS`; 60 of 60 clean at 8 ranks.
   - *Owner decisions open:* (1) merge the stack and in what order; what `par.solver` means at (h) and PETSc as a hard dependency; (2) default ICNTL(14) in source, and `KSPSetUp` reading INFOG(1)/(2); (3) move the comparator tests to `tests/unit` (new file) and document `PETSC_OPTIONS`; (4) audit items 2-6, 8, 9; (5) closing this row, or taking the measured reason: no tested solver passes at <= 5 s/step on wide at 32 ranks, and the solver survey is saturated (no new runs are queued).
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
   Owner, 2026-11-16, on the 3-step solve after its measured cost (dense
   interface, 3-4x the memory of exact): "how come? the interface shouldn't
   grow that crazy" and, on the untried ways to cut the dense interface cost
   (BLR compression of the Schur complement offered as the cheapest test on
   the existing spike; a sparse hierarchical interface factorization; fewer,
   larger subdomains): "keep pursuing that line" (verbatim.) Owner, 2026-11-16:
   "keep pursuing that line, optimize it with memory and speed." (verbatim.)
   Status: steps 1 and 2 measured 2026-11-16 (conductor 13, block below);
   the cheaper levers are in progress, each locked here before the next.
   Owner, 2026-11-16, on the parity tolerance for approximate solves (the
   question further down this row): "let's use 1e-9 to machine precesion for
   now" (verbatim.) Read as (a reading, not the owner's words): the 1e-9
   comparator stays the bar for every solver option for now, and iterative
   solves are run as tight as it takes to pass it, down to machine precision.
   Owner, 2026-11-16, after the survey of the other solver lines (GAMG CG
   best on memory and scaling; its untested lever, rigid-body modes per side
   of the fault's split nodes in the near-null space, offered as a second
   line): "so the winner is gamg cg. what's the near null space on the split
   node?" and then "keep proposing, keep exploring. 3 days budget. go"
   (verbatim.) Budget read as 2026-11-16 16:20 CST to 2026-11-19 16:20 CST;
   scope read as (a reading): the solver lines on this row, the 3-step solve
   and GAMG CG first, new proposals recorded here before they are run.
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
   7b77221). Independent re-derivation of the 2026-11-14 clean 64-rank pair, ICNTL(14) test and 200-step pair against `runs/`: setup of the clean pair is `runInfo.json` 7.4 / 169.6 s (8.1 / 182.1 were `EQ_WALL`), exact 64-rank loops 22.8-25.1 s, station-failure wording; the other figures recompute. Independent re-derivation of every figure above in this row
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
   - **CONDUCTOR 15, 2026-11-16 20:10 CST (autopilot; owner direction "keep proposing, keep exploring. 3 days budget. go", row 1 (g)). Budget read as wall clock 2026-11-16 16:20 to 2026-11-19 16:20 CST, bounded in practice by 100 tool calls per conductor; scope = the solver lines on this row; the 1e-9 comparator is the bar; proposals are written here before they run. Owner-held and not started: a merge-bound single-precision preconditioner (rule 1), `par.solver = 3` and the new source file, default ICNTL(14), everything under (h).**
    *Independent audit of #124-#127 figures (quantitative-claim auditor, 8 tool calls, read-only, 2026-11-16 20:09):* every ratio, memory figure, iteration count and file count recomputes from `runs/20261116_c14-*` (wide-alt 1.142 / 1.231, setup 0.536 / 0.426, 33 of 33 twice by `cmp_c11.py`; bp5 1.308, 120 of 120; wide-tight and 1e-12 rows; run order of the alternating batch is exa, spfa, exb, spfb by `runInfo.json` mtime). Corrected in the text above: the 6.9 s/step wide exact figure was unsourced (measured 6.06 at load 13, so the wide s2t break-even is 2,333 steps at 6.06, 2,049 at 6.9); the bp5 1.31x pair is sequential, not alternating; the 298-step solve-only break-even does not reproduce (293 / 389). Noted, not changed: wide-alt single slows 1.63x and exact 1.52x between the two load levels; bp5 `ex2` against `spf` is 2 of 120 (the same two files), only `ex` against `ex2` is 0 of 120; the s2t operator is exact for the fault, not for the six off-fault displacement files. `harvest.py` has no baseline name for `exa`/`exb`, so the wide-alt counts come from `cmp_c11.py` directly (re-run by the auditor, 33 of 33 twice). My re-run of the wide 1e-12 pair: `compared 33 fail 1` (`fltst_strk+000dp+000.dat`, 3 of 451 entries outside 2e-6), as the row states. No test pins any of these ratios (route to the test architect when a landing needs it).
    *Proposals, written before they run (all spike binary `spike-blr` 7d469dc, runtime options only, no code change, not merge-bound; each same-batch, alternating, one `heavy.lock` hold per run; host load recorded per run, so ratios are indicative):*
    1. **Single-precision factor at 32 and 64 ranks on wide** (`scratch/bench/c15/np.sh`, `runs/20261116_c15-np/`, 40 steps, rtol 1e-14, Fischer 2 / 10): 32 ranks exact / single pairs twice; 64 ranks single at default ICNTL(14) twice (does halved workspace lift the `-11` that hit 5 of 6 exact runs?), exact and single at ICNTL(14) = 100 twice; then 16 ranks single + BLR (ICNTL(35) = 1, CNTL(7) = 1e-8) against exact twice. Question: does rank scaling and the 64-rank failure change, and does BLR on the single-precision factor cut memory further without breaking 1e-9. Falsifier: any file failing 1e-9, or a `-11` at 64 ranks with single.
    2. **Fischer guess sweep** (`scratch/bench/c15/fsweep.sh`, `runs/20261116_c15-fsweep/`, bp1002, 8 ranks, 101 steps, single-precision MUMPS preconditioner, rtol 1e-14): no Fischer, model 1 sizes 2 / 4 / 10, model 2 sizes 4 / 10 / 20 / 40, exact control first, middle, last. Question: smallest guess size that holds iterations near 1 and the comparator at 0 of 25.
    3. **Slip-to-traction follow-up on `spike-s2t`** (6475482): an output-step solve for the six off-fault displacement files so all 25 bp1002 files compare, then a wide build with dense-operator memory and build time measured (replaces the arithmetic break-even). Dispatch to a specialist after batch 1 frees the host.
    4. Hand-built 3-step version stays last (owner order on this row).
    *Conductor 15 hand-off, 2026-11-16 20:40 CST (stopped after a context summary, 22 tool calls; times from `date`, state from git, PRs and `ps` only; a fresh conductor is requested).* Landed: #127 (69a9a21), #128 (cb8f879, audit corrections and proposals 1-4 above). CORRECTION by conductor 16, 2026-11-16 20:41 CST: the batch-A statement written here by conductor 15 (running) was false at merge; batch A was killed at 20:40 when conductor 15 stopped (it was a child of her session). Completed: `runs/20261116_c15-np/` e32a, e32b, s32a (see `matrix.out`; host load 38-74, only alternating pairs count). Relaunched DETACHED by conductor 16 at 20:41 CST: `scratch/bench/c15/driver.sh` PID 3479524 (ppid 1) running the remainder of `np.sh` (s32b, 64-rank, BLR 16) then `fsweep.sh` (batch B, `runs/20261116_c15-fsweep/`); `DONE` file per directory on finish; the scripts now live in the main `scratch/bench/c15/`. Budget read by conductor 16: 2026-11-16 16:20 to 2026-11-19 16:20 CST. Next step: when `DONE` exists in each directory, harvest with `cmp_c11.py` (baseline names `e32a`, `e64ia`, ... differ from `harvest.py`'s `_ex`), write the tables on this row, check the 64-rank default-ICNTL(14) runs for `-11` in `log/cycle0.log`. Proposal 3 branch `spike-s2t-2` 79c8fa6 (spike, not for merge; worktree `.claude/worktrees/agent-a1448e25e7231dcf8`; mira-volkov report, unaudited): 10 insertions in `src/solveTimeLoopPETScDD.f90`; `compared 25 fail 0` at 1e-9 on bp1002 8 ranks 101 steps (`runs/20261116_c15-s2t/*_v2s2t`, `*_v2ex`); but the stations record every step, so the full solve runs every step and the loop is 9.67 s against exact 16.41 s (`EQ_WALL`), not the 3.26 s of the 19-of-25 version; the recovery (station rows in the dense operator) and the wide build measurement are not done; the report also found `ddS2tBuild` leaving a nonzero `uloc` (fixed in that commit, step-1 station velocity). Re-check the comparator line and both `runInfo.json` before quoting any of it. Conductor 16, 2026-11-16: independent read-only audit of 79c8fa6 (lars-eriksson): comparator rerun `compared 25 fail 0 rtol 1e-09` is real, no fallback; the 9.67 s vs 16.41 s is one solve per step against two (`allocated(dout)` forces a full solve each step), build 40.2 s outside `EQ_WALL`. Follow-up `spike-s2t-3` 4f513b9 (porting specialist; spike, not for merge; worktree `.claude/worktrees/agent-ac1e90238e672c9bb`; UNAUDITED, audit dispatched): station rows added to the dense operator, version 1.24.6; reported comparator `compared 25 fail 0 rtol 1e-09`, 2 KSP solves (steps 1 and 101) against 202, `EQ_WALL` loop 5.66 s vs 13.07 s exact (back to back, load 67-89, indicative), build 49.0 s outside the loop figure, in-run operator-vs-solve worst relative difference about 2.5e-15 traction, 5e-16 displacement; non-station free-node `cons` stale between full steps (reported, not checked). Next: audit verdict, then the wide (bp8) build measurement (dense-operator memory and build time). AUDIT of 4f513b9 (lars-eriksson, read-only, 2026-11-16): comparator rerun `compared 25 fail 0 rtol 1e-09` against both the 1.24.5 and 1.24.6 exact runs, largest column-scaled difference 2.9e-13: the pass is real. BLOCKER F1 (`solveTimeLoopPETScDD.f90:887-889`, `:1057-1068`): a stop on a non-output step writes `disp`/`disp.r.nc`/`fault.r.nc` from stale free-node values (the gate never reaches it: both runs end at step 101, a full step); fix in `spike-s2t-4` 6449009 (`.or. stoptag == 1`, version 1.24.7). MAJOR: the 5.66 s loop figure excludes the 49.0 s build, so 101 steps cost about 58 s against about 16 s exact (break-even about 670 steps); the in-run check runs only at steps 1 and 101 and checks column assembly, not time stepping; suspected free-node `consv` zeroed on skipped steps (`:1030-1036`, no visible effect, fault pairs likely constrained). Gate gap: no test stops on a non-output step. F1 retest in `runs/20261116_s2t4-f1/` is VACUOUS (`par.mode = 2`, `exit_slip_rate = 1.02e-9` stopped at step 1, `EQ_WALL` steps = 1; comparator `compared 21 fail 0` means nothing for F1). Wide build at 32 ranks (`runs/20261116_s2t4-wide/`, mira's worktree): `rc=124` at the 7000 s runcase timeout with no `EQ_S2T built` line and no progress print, so the build takes more than 1.9 h at 32 ranks against 49 s for bp1002 (3564 columns), and the wide dense operator is unmeasured; the e32a baseline is 683 s wall, kspsolve 371.6 s over 40 steps. Next (dispatched, spike-s2t-5): a sensitive F1 retest and a per-batch progress line with a 12-batch wide probe to extrapolate build time and operator bytes per rank. WIDE PROBE, 2026-11-17 (conductor 16; binary 1.24.8 = 6449009 + a per-batch `EQ_S2T_PROG` line, uncommitted in the worktree `.claude/worktrees/agent-a63de7651de39cf07`, `runs/20261116_s2t5-wide/`; host load 90-100, indicative): the build is killed by the 2400 s timeout (`rc=124`) after 391 of 1913 batches; `ncol` = 61206, 32 columns per batch, 32 ranks (rank 0 owns no fault rows, `nr` = 0); 2718.8 s elapsed at batch 391, 3.3 s per batch over the last three, so the build is about 2.1 h at 32 ranks (about 6.9 s per batch averaged, 3.7 h if that rate held: the range is the load). Dense operator, ESTIMATE not measured: 6 x 6561 fault rows x 61206 columns x 8 B = about 19 GB over all ranks, against e32a's 181 GB summed peak RSS (+10%); the measured peak RSS is owed. Loop saving: e32a is 385 s over 40 steps, 9.6 s per step, so break-even is about 800 steps. PROPOSAL, written before running: a full wide s2t run at 32 ranks, 40 steps, `-eq_s2t_build`, to measure the build, peak RSS and the `compared N fail 0 rtol 1e-09` line against e32a; cost one 2.5-3 h slot on the heavy lock (`runs/20261116_s2t5-widefull/`, queued behind batch A, `timeout 21600`); the F1 stop-step retest is still open (mode-2 slip rate 3e-2 at step 1 and 3e-1 at step 21 means a stop at steps 40-60 may be unreachable on bp1002). CORRECTION 2026-11-17 02:00 CST (conductor 16): the probe was not killed by its 2400 s timeout (that killed only the shell; the MPI job kept the heavy lock) but ran on until `runcase.sh`'s own 7000 s timeout; the log ends at batch 986 of 1913 at 6785 s (`runs/20261116_s2t5-wide/`, mira's worktree; avg 6.9 s per batch, 5.4 s early, 11.5 s over batches 300-600, 4.5 s over 600-900; host load 90-100). So the build is about 3.6 h at 32 ranks (13,200 s at the run average), not 2.1 h; the 2.1 h figure above rested on a 3-batch rate and is withdrawn. No `EQ_S2T built` line yet, so no comparator, RSS or operator-size figure. The proposed full wide run is queued DETACHED (`scratch/c17/go_wide_long.sh` in that worktree, PID 3057725, ppid 1): it waits for the batch-A driver (3479524) and the requeued stop-step bisect (2258661; its first copy timed out waiting for the lock, 7200 s, and ran nothing) to exit, then runs `runs/20261117_s2t5-widelong/` with a 12 h limit (`runcase_long.sh`, 36000 s inside). BATCH A HARVEST, 2026-11-17 02:50 CST (conductor 16; `runs/20261116_c15-np/`, `DONE` present, `scratch/bench/c15/cmp_all.sh` re-prints the comparator lines; wide bp8, 40 steps, `EQ_WALL` loop s / kspsolve s / max rank RSS GB; host load 40-100, loops contaminated, same-config pairs differ up to 2x so only memory and the comparator count): 32 ranks: e32a 385.0/371.6/6.81; e32b INVALID (killed with conductor 15, log ends at `EQ_TOL`, no `result/`); s32a 365.7/357.8/3.58; s32b 713.2/701.3/3.65 (2x its twin: load); comparator vs e32a: s32a `compared 33 fail 0 rtol 1e-09`, s32b `compared 33 fail 0 rtol 1e-09`; memory 0.53x exact. 64 ranks: e64ia 419.5/407.2/4.56, e64ib 384.3/372.6/4.59 (exact vs exact: `compared 33 fail 0`); s64ia (ICNTL(14)=100) 1610/1590/2.47, comparator `fail 0`; s64ib (ICNTL(14)=100) 529.0/519.9/2.44 and s64db (default ICNTL(14)) 487.8/478.4/2.23: BOTH FAIL the comparator, 1 of 451 entries of `fltst_strk+000dp+000.dat` outside rtol 2e-6 (row 38 column 5, -29.9702 vs -30 and -29.9583 vs -30, against e64ia/e64ib); s64da (default ICNTL(14)) ABORTED (`rc=1`, `run.log` shows the prterun abort, no `-11` or `ICNTL` string in `run.log`, so the cause is not identified: the same options as s64db, which completed). So the single-precision preconditioner passes 4 of 6 wide runs at 1e-9 and fails the comparator in 2 of 3 completed 64-rank runs: not a bit-reliable candidate at 64 ranks, and no tolerance is bumped. 16 ranks, single-precision plus BLR (`-mat_mumps_icntl_35 1 -mat_mumps_cntl_7 1e-8`): e16a 175.0/168.4/11.90, b16a 120.8/115.5/5.28; e16b 178.4/170.1/11.97, b16b 140.1/134.7/5.35; comparator `compared 33 fail 0 rtol 1e-09` for both pairs; memory 0.44x exact, loop 0.69x and 0.79x (load-contaminated; MEASURED, not VERIFIED: no independent oracle). Stop-step bisect: its first requeue also ran nothing (variables not exported into the `bash -c`; renamed `BIS_DONE.empty_unexported_vars_20261117`); a fixed chain `scratch/c17/go_chain.sh` (PID 3133306, ppid 1; in `.claude/worktrees/agent-a63de7651de39cf07`) runs the bisect after the batch-B driver exits, then the long wide s2t run (supersedes the PID 3057725 waiter, stopped). Batch B (Fischer sweep) finished at 02:38, see below. BATCH B HARVEST AND F1 RETEST, 2026-11-17 (conductor 16; `scratch/bench/c15/cmp_fs.sh` re-prints the sweep): Fischer sweep, bp1002 8 ranks 101 steps, single-precision MUMPS PC, rtol 1e-14, `runs/20261116_c15-fsweep/` (`DONE`; host load 40-100): ALL 11 runs `compared 25 fail 0 rtol 1e-09` against ex1. Loop s / kspsolve s: exact ex1 14.04/12.42, ex2 11.71/10.40, ex3 11.72/10.38 (exact-vs-exact spread about 20%); single-precision without Fischer 29.07/27.75 (2.2x exact); Fischer model 1 size 2/4/10: 14.78/13.38, 14.31/13.14, 13.78/12.66; model 2 size 4/10/20/40: 15.67/14.23, 14.38/13.23, 14.85/13.61, 14.04/12.90. Reading: Fischer is what makes the single-precision PC viable on bp1002 (about 2x over no guess); model and size change the loop by less than the exact-vs-exact spread; no setting beats exact; the iteration-count column of `EQ_KSP` was not read, time only. F1 stop-step retest (bisect `runs/20261116_s2t5/` in the specialist worktree `.claude/worktrees/agent-a63de7651de39cf07`, mode 2 with `exit_slip_rate` 4e-2/6e-2/1e-1/2e-1 stops at step 2/2/2/11; 2e-1 stops at step 11, off an output step): s2t at 4f513b9 (1.24.6, built in `.claude/worktrees/agent-ac1e90238e672c9bb`) gives `disp.00011.nc`/`disp.r.nc` differing from the exact run in 163584 of 176607 entries, max 1.353e-02 against field scale 3.074e-02 (rel 4.402e-01), and `cmp_c11.py` reports the two files failing at 1e-9: F1 CONFIRMED as a wrong-output bug; s2t with the F1 fix (1.24.8, `spike-s2t-5` working tree, uncommitted) gives `compared 25 fail 0 rtol 1e-09` and 0 of 176607 differing entries in `disp.00011.nc` and `disp.r.nc` (nodal diff by `scratch/c17/nodal_diff.py`, which loads the netCDF values; `cmp_c11.py` does compare the `.nc` files). F1 is closed on the spike; the fix is a 1-line condition and would need the gate test (stop on a non-output step, `exit_slip_rate` 2e-1 on bp1002, 101 steps) if the operator ever becomes merge-bound (owner decision).
  - **CONDUCTOR 14, 2026-11-16 (autopilot; owner direction 'keep proposing, keep exploring', 3-day budget, board commit d54addf; tolerance bar 1e-9 per 0d9923b). Owner-readable paragraph first, then the measurements. Tool calls so far: ~45.**
     *Verdict on the 3-step (exact substructuring) solve, for the owner.* It is exact and passes the 1e-9 comparator (bp1002 25/25, bp5 120/120, gs10 33/33), and it does not beat exact MUMPS: on bp1002 at 8-32 ranks it needs 3.4-4.1x the memory (peak per rank 3.6-4.1x, summed 3.4-4.0x), 2.3-3.6x the setup and 1.7-2.3x the loop time; on gs10 at 16 ranks 4.6x peak memory (20.7 against 4.47 GB), 4.1x summed (261 against 63 GB), 3.3x setup (268 against 81 s), 1.6x loop (4.45 against 2.71 s per step). Every cheaper variant tried (BLR on the interface, deluxe scaling, two-level BDDC) closed negative or gained 5-40 % of the gap only. The one remaining lever is a hand-built version (MUMPS native Schur complement per rank plus a dense parallel factorization of the interface): by the estimate on this row it would sit at about 1.4x exact memory and 6x exact setup, so it could win only on loop-time scaling at 64 or more ranks where exact MUMPS stops scaling; cost: new source of 300-500 lines on a spike branch, about one specialist-day, and no memory gain. I do not recommend it before the two measured routes below are exhausted. Numbers re-derived from the raw run files by an independent audit this run; the corrections it found are applied in place above (marked 'Corrected').
     *Independent numeric audit of the conductor-13 blocks (2026-11-16, read-only, 24 tool calls):* all 36 setup/loop/peak/summed cells of the BLR table, the step-1 matrix sizes, the step-3 table and the failing-file counts re-derive from `table.csv` and the case logs; 11 text discrepancies, all corrected in the blocks above (BLR on bp1002 does compress the Schur factor by 13.5 %, 86.5 % kept, not 100 %; ratio ranges; `-log_view` percentages are of run time, not of setup; the 8.77 M entries are the largest rank; gs10 BLR figures are `EQ_WALL`, the controls `runInfo.json`; two-level BDDC '20-25x' is against all-primal, 31-37x against exact). The 1.4x / 6x / 86 GB / 8 min estimates are arithmetic from the sizes: interface 76,992 dofs, dense triangle 76,992^2 x 8 B / 2 = 23.7 GB (1.5 GB per rank at 16 ranks), Cholesky n^3/3 = 1.5e14 flop, 8 min assuming about 20 Gflop/s per rank (assumed, not measured). No test pins these figures to the raw files (gap noted).
     *GAMG levers on gs10, 16 ranks, 3 steps, rtol 1e-8 (`runs/20261116_c13-gamg-opts/`, harvest `scratch/bench/c13/harvest.py`, MEASURED, load 20-24):* base GAMG 43.5 iterations per solve, loop 40.4 s, peak 1.44 GB; two smoother sweeps 43.5 its, 47.3 s (no gain); strength threshold 0.01 40.2 its, 41.0 s (no gain); threshold 0.05 27.7 its but loop 73.2 s, setup 160 s, peak 4.13 GB (worse); aggressive coarsening off: aborts at step 1, `DIVERGED_PC_FAILED` / `SUBPC_ERROR` (KSP -11, `log/cycle0.log` line 41-43); hypre BoomerAMG with nodal coarsening 6 and rigid-body interpolation (variant 3): SEGV in all 16 ranks at step 1 (`log/cycle0.log` line 31). So none of those levers cuts the iteration count into a time gain; the iterations per solve are not what the next two routes attack.
     *Fischer guess and single-precision MUMPS preconditioner, bp1002, 8 ranks, 101 steps (202 solves; the exact run ends at max slip rate 5.66e-2 m/s), CG on the true residual, same-script, same binary (spike 7d469dc, 1.24.4), comparator `cmp_c11.py` 25 files at 1e-9, MEASURED; loop = `runInfo.json` time_loop_seconds, setup = `EQ_WALL`; peak GB largest rank / summed GB; host load at launch in brackets (64 cores, 28 users, so timings are indicative; each batch has its own same-load exact control):*
     | batch | variant | its per solve (max) | setup s | loop s | peak / summed GB | files failing |
     |---|---|---|---|---|---|---|
     | 1 (load 16-17), rtol 1e-10 | exact MUMPS (control) | 1 | 3.19 | 12.08 | 0.50 / 3.43 | 0 |
     | | GAMG, last-solution guess per tag | 9.63 (12) | 0.94 | 36.42 | 0.32 / 2.36 | 3 |
     | | GAMG, Fischer model 1 size 6 | 1.96 (12) | 0.97 | 9.92 | 0.32 / 2.36 | 10 |
     | | GAMG, Fischer model 2 size 10 | 2.01 (12) | 0.94 | 10.43 | 0.32 / 2.36 | 4 |
     | | same, no per-tag guess | 2.01 (12) | 0.96 | 12.48 | 0.32 / 2.36 | 4 |
     | 2 (load 49-51) | exact MUMPS (control) | 1 | 3.25 | 13.84 | 0.48 / 3.47 | 0 |
     | rtol 1e-12 | GAMG, last-solution guess | 12.21 (14) | 0.96 | 51.55 | 0.32 / 2.36 | 0 |
     | | GAMG, Fischer 1 / 6 | 3.48 (14) | 0.97 | 18.18 | 0.32 / 2.36 | 1 |
     | | GAMG, Fischer 2 / 10 | 3.38 (14) | 0.95 | 18.73 | 0.32 / 2.36 | 1 |
     | rtol 1e-14 | **GAMG, Fischer 1 / 6** | 5.23 (17) | 0.97 | 26.11 | 0.32 / 2.36 | **0** |
     | 3 (load 48-49), rtol 1e-12 | exact MUMPS (control) | 1 | 3.25 | 13.49 | 0.50 / 3.47 | 0 |
     | | **MUMPS single precision Cholesky as the CG preconditioner** (`-pc_precision single`) | 2.01 (3) | 2.55 | 21.25 | 0.33 / 2.25 | 0 |
     | | **same plus Fischer 2 / 10** | 1.06 (3) | 2.33 | 14.44 | 0.33 / 2.24 | **0** |
     *Reading (a reading from the table, one host, one case; not the owner's words).* (1) The previous-solution reuse the board had (guess per tag) leaves GAMG at 10-12 iterations; a Fischer projection onto the last 6-10 solutions cuts that to 2-5 at the same tolerance, because K is constant and consecutive right-hand sides are nearly collinear. At the tolerance that passes (1e-14) GAMG with Fischer is 1.9x exact in loop (26.1 against 13.8 s), 3.4x shorter in setup and 1.5x smaller in peak memory (0.32 against 0.48 GB), and without Fischer at 1e-12 it is 3.7x (51.6 s); so Fischer is worth about 2x of the gap and does not close it on bp1002 at 8 ranks. (2) A single-precision exact factor under CG with Fischer passes the comparator at 1e-12 (0 of 25, 1.06 iterations per solve), costs 1.07x exact in loop (14.44 against 13.49 s), is 0.72x in setup, and holds 0.66x the peak and 0.65x the summed memory (0.33 against 0.50 GB; 2.24 against 3.47 GB). It is the first variant on this row that matches exact in time and beats it in memory while passing the 1e-9 comparator. Its limit here: single-precision solves are not faster than double at this size (without Fischer it needs 2 iterations and 21.3 s), so the gain is memory only. PETSc 3.25.5 in `eqquasi-petsc` exposes it at run time (`PETSC_HAVE_MUMPS_MIXED_PRECISION`; `-pc_precision single`), no source change. (3) hypre BoomerAMG with default settings on bp1002 (8 ranks, rtol 1e-10, per-tag guess): 28-29 iterations per solve and about 9.5 s per solve (46 of 101 steps in 833 s, loop never finished; I stopped it by `kill` of its recorded timeout and prterun PIDs, 16:53 CST to free the lock); about 50x slower than GAMG there. (4) Not shown: whether (2) keeps its advantage where the factor memory is large (wide: exact 12.1 GB per rank) and over an event-resolving run (bp5, 4483 steps); both are queued, below. Exact control repeats: setup 3.19-3.25 s and loop 12.1-13.8 s across the three batches at load 16-51 (the loop varies 14 % with load; the setups agree to 2 % because they are in the same code path).
     *Event run and wide, single-precision MUMPS preconditioner + Fischer 2/10 at rtol 1e-12 (same driver and comparator, MEASURED 2026-11-16 17:37-17:51, load 46-56, indicative; `runs/20261116_c14-bp5-event/`, `runs/20261116_c14-wide-single/`).* bp5, 8 ranks, 5000 steps (event run), exact against single+Fischer: setup 3.69 / 2.59 s, loop 810 / 564 s (0.70x; loads 48 / 56), 1 / 0.83 iterations per solve (max 3), peak 0.57 / 0.41 GB (0.72x), summed 3.73 / 2.50 GB, **2 of 120 files fail** (`srfst_strk000st016dp-010.txt`, `srfst_strk000st032dp000.txt`: 4 and 5 of 31,381 entries outside the station tolerance, values near 1e-7 to 1e-8; cause not tested, tighter rtol not yet run). Wide, 16 ranks, 40 steps (aseismic: final max slip rate 9.3e-11 m/s, so no event), exact against single+Fischer: setup 479 / 227 s (0.47x), loop 461.5 / 574.1 s (1.24x; 11.5 against 14.4 s per step, loads 50 / 46 so both are slower than the 6.06 s/step measured at load 13 (wide-alt exact a; the 6.9 s/step quoted before was unsourced, audit 2026-11-16)), 1 / 1.38 iterations per solve (max 3), **peak 11.8 / 6.61 GB per rank (0.56x), summed 175.5 / 93.6 GB (0.53x)**, **1 of 33 files fails** (`fltst_strk+000dp+000.dat`, 3 of 451 entries, row 32 column 5: -29.20 against -30.x, the log10 slip-rate column at its floor value; earlier board comparisons excluded that column by hand, the comparator does not, and I have not exempted it). Reading (not the owner's words): the single-precision factor halves memory and setup at wide for a 1.24x loop, is faster than exact on the bp5 event run, and misses the 1e-9 comparator by 1-2 files in the two cases, in near-zero entries; whether a tighter rtol clears them is the next run (bp5, rtol 1e-14; wide is a 25-minute run). *bp5 event run, tighter tolerance and GAMG, MEASURED 2026-11-16 18:00-18:51 (`runs/20261116_c14-bp5-event/` second exact control and GAMG, `runs/20261116_c14-bp5-tight/`; load 47-56).* Single-precision MUMPS + Fischer 2/10 at **rtol 1e-14: 0 of 120 files fail**, 1.25 iterations per solve (max 3), setup 2.55 s against 3.49 s, loop 852.7 s against 652.0 s for the exact control run made after it at the same load (47-48; sequential, not alternating; 1.31x), peak 0.35 against 0.53 GB (0.66x), summed 2.49 against 3.76 GB (0.66x). So the 2 files that failed at 1e-12 clear at 1e-14; at 1e-12 the loop was 564 s against exact 810 s (first control) and 765 s (second control, 3.78 s setup, 0 of 120 fail), the exact loop alone varying 652-810 s across four runs at load 47-56 (24 %), so loop ratios from 0.70x to 1.31x across tolerance and control are inside the load noise only partly: a same-load pair at one tolerance (exact beside 1e-14) reads 1.31x, and I take that as the figure to quote. GAMG CG + Fischer 1/6 at 1e-14: 5.02 iterations per solve (max 17), setup 1.31 s, loop 1302.7 s (1.7-2.0x exact, 1.6-1.7x against the 765-810 s controls), peak 0.35 / summed 2.68 GB, 1 of 120 files fail. Reading (not the owner's words): on the event run the single-precision factor passes the comparator at 1e-14, halves nothing in time but cuts memory to 0.66x at 1.3x the loop time, and GAMG is both slower and still misses one file; the wide run at 1e-14 follows.
     *Wide, 16 ranks, 40 steps, single-precision MUMPS + Fischer 2/10 at rtol 1e-14, MEASURED 2026-11-16 18:55-19:16 (`runs/20261116_c14-wide-tight/`, `gamg8.sh`; the exact control ran after it, load 43.7 before the single run, 55.6 before the exact).* **0 of 33 files fail** (`cmp_c11.py`; the `fltst` floor-column miss at 1e-12 is gone), 1.66 iterations per solve (max 4), setup 206.7 s against 487.8 s (0.42x), peak 6.28 against 11.66 GB per rank (0.54x), summed 93.0 against 176.6 GB (0.53x), loop 247.2 s against 425.3 s. The loop figure is not usable as a ratio: the two runs had different loads, and the single-precision run at 1e-12 on the same case took 574 s at load 46 (and the exact 461.5 s at load 50), so the same configuration class ranges 247-574 s at similar load with 1.4-1.7 iterations per solve in both; I cannot separate the solver from host contention here, and a same-load repeat with the exact run beside it (alternating, not in sequence) is needed before any wide speed claim. What the table does support at wide: the exact memory halves (0.53-0.54x), the setup more than halves, and the 1e-9 comparator passes 33 of 33 at 1e-14 and 32 of 33 at 1e-12. Combined with bp5 (120 of 120 at 1e-14, memory 0.66x, loop 1.31x on a same-batch pair) this makes single-precision MUMPS under CG with Fischer the variant on this row that passes the bar at about half the exact memory at wide; it needs no source change in the solver library, but in EQquasi it is a spike flag (`-eq_plain_ksp`, `-eq_guess_per_tag`) behind `par.solver = 3`, so a merge-bound form needs the owner's agreement on `par.solver` (rule 1) and a test that triggers the path.
     *Wide, same-load repeat, alternating exact and single-precision + Fischer 2/10 at 1e-14, MEASURED 2026-11-16 19:20-20:05 (`runs/20261116_c14-wide-alt/`, `gamg9.sh`; 16 ranks, 40 steps; loop s / load before):* exact a 242.5 / 13.3, single a 276.9 / 25.3, exact b 367.6 / 41.4, single b 452.5 / 53.4. Setup 367.9 / 197.3 s and 459.9 / 196.1 s (single 0.54x and 0.43x), peak 12.08 / 6.53 and 12.30 / 6.28 GB (0.54x, 0.51x), summed 175.5 / 94.4 and 176.0 / 93.4 GB (0.54x, 0.53x), 1.66 iterations per solve (max 4), 0 of 33 files fail in both single runs. **Loop: single is 1.14x (a) and 1.23x (b) the exact run beside it**, each single run starting at a load 12 higher than its exact partner, so the ratio is biased against it; both solvers slow down 1.5x between load 13-25 and 41-53. The 0.58x reading in the previous paragraph was a load artifact and is withdrawn: single-precision MUMPS + Fischer at wide costs about 1.1-1.25x exact in loop time, 0.4-0.5x in setup and 0.5x in memory, and passes 33 of 33 at 1e-14 twice.
     *Slip-to-traction operator, spike (owner not asked beforehand; I run it as a spike on the reading that 'keep exploring' covers it, and tell the owner here), bp1002, 8 ranks, 101 steps, MEASURED, spike branch `spike-s2t` 6475482 on origin (+293 lines in `src/solveTimeLoopPETScDD.f90`, version 1.24.5 in its worktree, not merge-bound), `runs/20261116_c14-s2t/`, load 48-50.* Design: with `-eq_s2t_build` the dense operator S (traction on the fault from the constrained displacements) is built once by `MatMatSolve` on the MUMPS factor, 3,564 columns (6 per split-node pair, 594 pairs: three components on each of two nodes, not the 1.2 k the board's sizing assumed) plus one column for the far-field term (alpha times a fixed pattern; the path checks to 1e-12 each step that the far field really has that form, and never fired), distributed by rows; each step then does an `Allgatherv` of the 3,564 slip values and a local matmul instead of the KSP solve and `ddBoundFtKu`. The exact path is linear in the constrained displacements (tags 0 and 1 are the same operator), not affine. Parity against exact (`cmp_c11.py`, 25 files, 1e-9): the 19 fault-side files pass; the 6 that fail (`disp.00001.nc`, `disp.00101.nc`, `disp.r.nc`: 0.6 / 0.49 of the field scale; three `srfst_*` station files) fail because the s2t path never solves for free-node displacement (`cons(free)` stays 0): not a defect in S, a missing output solve; a merge-bound version would run one exact solve at each output step. Timing: build 35.69 s (0.0100 s per column, outside the loop); loop 3.256 s (0.0322 s/step; apply 0.0021 s/step) against exact 14.58 s (0.1443 s/step; KSPSolve 0.1240 s/step; the exact run carried `-eq_spike_log`, whose residual MatMult and norms the s2t path skips, so the solve-only figure is cleaner); **break-even 318 steps (the earlier 298 solve-only does not reproduce: 293 with the apply time, 389 with the s2t loop per step)** on bp1002. Wide, arithmetic only (mine, not measured): 6 x 10,201 split nodes = 61.2 k columns; at the measured 0.231 s per column (wide, 16 ranks, 32 columns per call, board above) the build is about 14,100 s (3.9 h), and the break-even against 6.9 s/step is about 2,000-2,400 steps for an apply of 0-1 s per step (the apply and the dense memory, 61 k columns by the number of traction rows, about 15 GB if three per pair, an assumption, are not measured). So the operator pays back inside one event-resolving run (bp5 4,483 steps) on bp1002 and would at wide only for runs of several thousand steps or repeated runs on one mesh and moduli. Why it matters against the solver lines: it is exact (the same factor), so no tolerance question arises.
     *Proposals, written before running, ranked by expected gain over cost (suggestions from the invoking session and my own sizing, not owner decisions; checked against what this row rules out: none repeats a rejected lever).* (1) Single-precision MUMPS preconditioner + Fischer at wide and on bp5: gain up to about 0.65x memory at exact-level time; cost two runs, already queued; ranked first because it needs no source. (2) Slip-to-traction operator with an output-step solve for free-node displacement and a wide build: gain 4.7x per step on bp1002 and about 3-7x projected at wide at a one-time build of hours; cost one specialist-day to add the output solve and a MatMatSolve build that stores only the fault rows; ranked second, exact, and the only variant that attacks the time-loop cost rather than the memory. (3) Fischer settings sweep (size, model) on the GAMG line and on (1) at a tolerance 1e-11 to 1e-13 across the 202 solves: gain below 2x on GAMG (floor about 2 iterations); cost minutes per run; ranked third. (4) Krylov recycling or deflation beyond Fischer: not tried; this build lists `KSPHPDDM` types and `PCDEFLATION`; expected gain over Fischer small (Fischer already reaches 1-2 iterations with the single-precision factor), so not run. (5) Hand-built exact 3-step with a dense interface factorization: new source, 300-500 lines, about one specialist-day, memory 1.4x exact (estimate), only a time-scaling gain at 64+ ranks; ranked last, held until (1)-(2) say the time loop at high rank counts is the remaining problem.
   - **3-step line, conductor 13, 2026-11-16 (autopilot, owner direction 7a115a5, tolerance 0d9923b); step 1: what the coarse problem is, step 2: BLR. MEASURED, same-script runs, spike `spike-sub` 7d469dc (binary eqquasi-1.24.4 built after the last source change, tree clean), not merge-bound; hosts shared, load 5-20, timings indicative; evidence `runs/20261116_c13-diag/`, `_c13-blr-bp1002/` (`table.csv`), `_c13-blr-gs10/`, harness `scratch/bench/c13/`.**
     *Step 1, what the interface problem is* (`-ksp_view`, MUMPS `ICNTL(4)=2`, `-log_view`, bp1002, 8 ranks, 3 steps). The coarse matrix is assembled sparse, not dense: type `mpiaij`, 12,924 rows, 70.1 M nonzeros against 167 M for a full matrix (42 % filled), factorized by parallel MUMPS Cholesky into 50.5 M entries (60 % of a dense triangle, 0.40 GB of factor entries; MUMPS internal data 3.4 GB summed (allocated total `INFOG(19)` 3.419 GB; `INFOG(17)` is the 3,412 MB estimate), 0.65 GB on the largest rank, `INFOG(18)`). It is block-dense: each subdomain's interface (2.8-3.5 k dofs) is one fully coupled block, and the global matrix is their union. On gs10 (16 ranks) the same: 76,992 rows, 1.26 G nonzeros against 5.93 G (21 %), factor 860 M entries (29 % of dense), MUMPS internal 70 GB summed. So "dense" is right for the factor (the fill is 29-60 % of a full triangle) and wrong for the assembled input (21-42 %); a sparse input does not make the factorization sparse here because the Schur complement has no sparsity inside a subdomain's interface. Local solves: Dirichlet and Neumann are two factorizations of the same interior matrix (largest rank: order 18,972, 8.77 M factor entries each; local orders 18.1-19.0 k and 7.9-8.8 M entries across ranks; interior 18.1-19.0 k, interface 2.8-3.5 k per rank, one subdomain per rank, all interface dofs primal). Setup is 11.8 s on bp1002 at 8 ranks (`-log_view`): local Schur by one multi-right-hand-side `MatMatSolve` 4.5 s, `PCBDDCCorr` (coarse basis) 5.1 s, coarse setup 3.9 s (coarse MUMPS factor 3.3 s, `MatConvert` MATIS to MPIAIJ 2.5 s); these three nest or overlap and sum to 13.5 s, more than the 11.8 s setup, and their `-log_view` percentages (29, 33, 24 %) are of total run time, so they are not shares of the setup (38, 43, 33 % of it). Memory: the coarse factor (and its MPIAIJ input) is not the whole of it. gs10: peak RSS 20.7 GB on the largest rank, 262.8 GB summed (exact 4.5 and 63.3), of which the coarse MUMPS factor is 70 GB; the rest (about 190 GB) is outside the coarse factor and consistent with the dense per-rank coarse-basis block, interior order 115-120 k times 7.3-12.1 k primal columns times 8 B (about 7-12 GB per rank), which exists whenever every interface dof is primal; that attribution is a reading (the buffers were not instrumented), the control would be a run with fewer primal dofs. Consequence for the options: (a) BLR can only touch the 70 GB of coarse factor; (b)/(c) are about the coarse factor and the primal count.
     *Step 2, BLR* (bp1002, 11 steps, CG rtol 1e-12 on the true residual; setup s / loop s / peak GB largest rank / summed GB; ex = exact MUMPS solver 3, ap = all-primal; comparator = `cmp_c11.py`, 25 files, 1e-9). np 8: ex 3.09 / 1.65 / 0.48 / 3.44; ap 11.13 / 2.76 / 1.98 / 13.63 (1 iteration, 0 fail); coarse BLR 1e-8 13.37 / 5.66 / 2.01 / 14.27 (2 iterations, 0 fail); coarse BLR 1e-4 12.86 / 11.56 / 1.97 / 14.15 (5 iterations, 0 fail). np 16: ex 2.27 / 0.99 / 0.33 / 4.57; ap 6.45 / 2.23 / 1.19 / 16.40; BLR 1e-8 7.82 / 4.01 / 1.38 / 16.07 (2 its, 0 fail); BLR 1e-4 7.18 / 8.11 / 1.24 / 16.57 (5 its, 1 of 25 fails, `disp.00001.nc` 3.8e-9: rtol 1e-12 is not tight enough there, the displacement amplification noted below). np 32: ex 1.81 / 0.98 / 0.20 / 5.94; ap 4.14 / 1.94 / 0.76 / 19.88; BLR 1e-8 5.11 / 2.53 / 0.84 / 21.44 (2 its); BLR 1e-4 4.95 / 4.49 / 0.96 / 20.25 (5 its, 0 fail). **Coarse BLR gives no memory: MUMPS reports effective entries of 86.5 % of the theoretical count on bp1002 (`INFOG(35)` 43.6 M against `INFOG(29)` 50.5 M) and 95.0 % on gs10 (816.6 M of 859.9 M), so the Schur factor compresses 5-13 % at 1e-8, and the summed RSS is unchanged or up (np 16 BLR 1e-8 is 2 % down, the other seven cells up); it costs 1.3-4.2x in loop time against all-primal (BLR solve phase plus the extra iterations).** (Corrected 2026-11-16 by the numeric audit: this paragraph first read 100.0 % for bp1002, which is the theoretical-count line, and 2-4x.) gs10, 16 ranks, 3 steps, coarse BLR 1e-8 (BLR timings are `EQ_WALL`, the ap and ex figures are `runInfo.json`; the two differ by under 1 %): setup 309.9 s (ap 268.4, ex 80.8), loop 20.3 s (6.8 s/step; ap 4.45, ex 2.71), 2 iterations per solve, peak 20.72 GB largest rank and 262.8 GB summed (ap 20.68 and 261.2, ex 4.47 and 63.3), 33 of 33 files pass (host load 11 before, 20 after). BLR on the local solves (Dirichlet and Neumann, as for the earlier non-all-primal BDDC): CG aborts with `KSPConvergedReason = -8` (indefinite preconditioner) at step 1 at BLR 1e-8 and 1e-10 on bp1002 (8 ranks; the same -8 at 16 and 32 ranks at 1e-8, with and without coarse BLR, all 8 failing logs read); 1e-12 runs (1 iteration, 0 fail) with no memory change (peak 2.00 GB, summed 13.79 against 13.63) and loop 3.26 s against 2.76 (local compression not read from MUMPS statistics); 1e-8 with GMRES runs (2 iterations, 0 fail) at 8.66 s loop against 2.76. Reading under the owner's tolerance: every BLR setting that passes needs 2-5 outer iterations and is slower and no smaller than all-primal. Step 2 closes negative: BLR on the interface does not close the gap to exact MUMPS (peak 3.6-4.1x and summed 3.4-4.0x memory, setup 2.3-3.6x, loop 1.7-2.3x for all-primal on this table (the 2026-11-14 11-step table below gives 3.1-4.1x, 2-3.7x, 1.7-2.9x); BLR leaves memory unchanged and adds loop time). Not verified: no independent audit yet, no tuning of the MUMPS BLR cluster parameters (`ICNTL(38)`, `CNTL(7)` only), one host.
   - **3-step line, conductor 13, 2026-11-16: step 3 (runtime-option levers) and the verdict so far; GAMG line, source reading and proposals written before running. MEASURED, same-script, bp1002 11 steps, host load 13-33 (indicative), `runs/20261116_c13-levers/table.csv`.**
     *Deluxe scaling* (`-pc_bddc_use_deluxe_scaling`, 8 ranks, same batch as its all-primal control): setup 10.34 against 10.89 s, loop 2.04 against 2.57 s, peak 1.86 against 1.98 GB, summed 12.74 against 13.64 GB, 1 iteration, 0 of 25 fail: a 5-20 % gain, no change of verdict. *Two-level BDDC* (`-pc_bddc_levels 1`, the interface problem solved by a second BDDC with MUMPS local solves, which is the run-time form of a hierarchical interface factorization; rtol 1e-12, 0 of 25 fail in all three): 16 ranks ratio 4: 25 iterations, setup 32.1 s, loop 33.8 s, summed 11.67 GB; 32 ranks ratio 4: 41 iterations, 25.1 s, 32.6 s, 12.63 GB; 32 ranks ratio 8: 24 iterations, 53.8 s, 33.4 s, 13.70 GB; against all-primal 16 / 32 ranks: setup 5.98 / 3.88 s, loop 1.69 / 1.32 s, summed 16.41 / 21.30 GB, and exact 2.23 / 1.77 s, 1.07 / 0.91 s, 4.56 / 5.95 GB. So it saves 29-41 % of memory against all-primal, stays 2.1-2.6x above exact in memory, and is 20-25x slower in loop than all-primal (31-37x than exact; it iterates on the interface). Not a route to the bar. *Dense ScaLAPACK Cholesky for the coarse problem* (`-pc_bddc_coarse_pc_factor_mat_solver_type scalapack`): aborts, "MatSolverType scalapack does not support matrix type mpiaij"; the coarse matrix type is not reachable by an option (`-pc_bddc_coarse_mat_type` is not read), so it needs source that builds the interface matrix itself (spike only, see below).
     *Where the line stands (a reading from the measured sizes, with an estimate where marked).* Exact MUMPS already eliminates the interfaces sparsely inside one nested-dissection tree. One-level substructuring replaces that by a factor of the full interface: on gs10 at 16 ranks that factor alone is 860 M entries (6.9 GB of entries, 70 GB of MUMPS memory summed) against 63 GB for all of exact MUMPS. A hand-built exact version (MUMPS native Schur complement per rank, `ICNTL(19)`, then a dense ScaLAPACK Cholesky of the 76,992-order interface) would hold about 1.2 GB local dense Schur plus 1.2 GB local factors plus 3 GB of the interface per rank, about 86 GB summed against 63 (an estimate, 1.4x exact), and the dense factorization is 1.5e14 flop (about 8 min at 16 ranks, 6x the exact setup; estimate). So at 8-32 ranks no variant tested or sized reaches exact MUMPS in memory or setup, and the per-step time stays 1.7-2.9x slower (all-primal) or worse. The line's possible gain is time scaling at rank counts where exact stops scaling (exact loop is flat from 16 ranks, board above), not memory; what that needs is a hand-built spike (new source, estimated 300-500 lines: local MUMPS Schur, interface assembly, interface factorization, three-phase solve, on a spike branch only) and a rank count of 64 or more to show a crossing; a merge-bound version would additionally need rule 1 agreement on `par.solver` and PETSc. Proposal, not run: build that spike only if the owner wants the time-scaling question answered; expected gain: loop time at 64+ ranks, no memory gain; cost: one specialist dispatch, about a day. Considered and rejected: more BLR (no compression at 1e-8, 95-100 % of entries kept); fewer subdomains than ranks (MATIS takes one subdomain per rank, source change in PETSc's use, and the interface only shrinks as the rank count falls toward exact MUMPS itself); sharing the Dirichlet and Neumann factorizations (identical matrices, 1.55 of 13.6 GB on bp1002 at 8 ranks, a 11 % ceiling, and PETSc has no option for it).
     *GAMG line, near-null space per fault side: the hypothesis does not apply as stated (source read by the mesh specialist, 8 calls, my spot check of the compset lines).* K has no term tying the two copies of a split node (the fault enters only through the right-hand side), but each fault is a finite planar patch: bp1002 faults A and B span z from 0 to -20 km in a domain to -60 km and have x tips inside the mesh (`compset/test.bp1002.qdc.zone.2500/user_defined_params.py:93-94,110-112`), and the nodes just outside a patch are single, so the sides are tied around the tips and below. No block moves rigidly on its own, K is SPD with the boundary-constrained dofs removed, and bp1002 has two faults at y = 0 and -5 km with overlapping x, so a side rule is ambiguous; a hard half-space cut would put a jump through intact material (a high-energy vector, not near-null). Not tested: a smooth-weighted per-side mode, which is a different transform from the one proposed and for which there is no evidence. The 63 iterations on gs10 and wide against 10 on bp1002 are the same on both (board above), which points at mesh grading and anisotropy (element size 10 m at the fault to km far away) rather than at the near-null space. Proposal, run next (options only, gs10 at 16 ranks as the proxy since its iteration count equals wide's, 3 steps, rtol 1e-8, same base as the earlier GAMG runs; `runs/20261116_c13-gamg-opts/`): strength threshold 0.01 and 0.05, aggressive coarsening off, two smoother sweeps, and hypre BoomerAMG with nodal coarsening and the rigid-body near-null space; expected gain: iterations 63 down toward 30-40 if grading is the cause (a guess); cost: one batch of about 7 short runs; then the best at wide 16 ranks at rtol 1e-14 against 58-74 iterations. Rejected: reuse of the previous step's solution (already in).
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

   - **GAMG CG tuning on bp1002, 16 ranks, 11 steps, MEASURED 2026-11-14 (runtime options only)**
     (`runs/20261114_scaling/`, spike binary, CG on the unpreconditioned norm, `-eq_near_null`;
     exact at 16 ranks: setup 2.19 s, loop 0.94 s, 0.33 GB per rank; comparator = the table's,
     25 files). Setup is `-ksp_view`-confirmed to be GAMG V-cycle, 4 levels, aggressive coarsening
     on level 0 and 1 smoothing step in the prolongator. Block size 3 and the rigid-body near-null space
     are not printed by `-ksp_view`: they are set at `src/solveTimeLoopPETSc.f90:199` and `:354`, and the
     effect is shown by the control without `-eq_near_null` (8 ranks): 78.0 iterations per solve
     against 10.5 with it. Timings below ran at load 4.4 to 9.2 rising through the sequence (one
     1-minute figure per run, 1-2 s loops), so iterations, the load-independent figure, carry the
     comparison; loops are indicative.

     | variant | KSP its/solve (max) | setup s | loop s | peak RSS GB per rank | worst diff from exact | files failing |
     |---|---|---|---|---|---|---|
     | reference (rtol 1e-8) | 10.2 (11) | 0.53 | 2.13 | 0.200 | 1.9e-08 | 13 of 25 |
     | rtol 1e-6 | 7.3 (8) | 0.55 | 1.42 | 0.200 | 5.4e-07 | 21 of 25 |
     | rtol 1e-10 | 13.2 (14) | 0.54 | 2.17 | 0.200 | 3.8e-09 | 3 of 25 |
     | aggregation threshold 0.02 | 10.2 (11) | 0.79 | 2.32 | 0.165 | 1.9e-08 | 13 of 25 |
     | unsmoothed aggregation (`-pc_gamg_agg_nsmooths 0`) | 22.0 (25) | 0.42 | 3.40 | 0.198 | 1.1e-07 | 16 of 25 |
     | Richardson + SOR smoothers | 8.3 (9) | 0.54 | 1.64 | 0.201 | 1.1e-07 | 16 of 25 |
     | 2 smoother iterations | 10.2 (11) | 0.53 | 1.68 | 0.202 | 1.9e-08 | 13 of 25 |
     | coarse equation limit 1000 | 10.2 (11) | 0.60 | 1.72 | 0.200 | 1.9e-08 | 13 of 25 |

     Reading: only the tolerance and the smoother change the iteration count; the threshold, a second
     smoothing step and the coarse limit leave it at 10.2. Even rtol 1e-10 (13.2 iterations) leaves
     3 of 25 files above the gold comparator's 1e-9, so passing it needs a tighter rtol than any run here
     (iterations rise 1.5 per decade of rtol, 7.3 / 10.2 / 13.2 at 1e-6 / 1e-8 / 1e-10; whether 1e-11 or 1e-12 passes the comparator is not run and not estimated): the cost of
     exactness looks small in iterations, and is not shown in wall time. The two 1.1e-7 rows are float32-floor `slip_rate` differences, not a
     convergence figure (the residual tolerance is the same as the reference). Exact 0.94 s loop
     against GAMG 1.4-3.4 s: GAMG is slower at this size in every variant. Wide tuning and the
     scaling curve are running (`runs/20261114_scaling/`).
   - **Strong scaling of GAMG CG against exact MUMPS on wide, 3 steps, MEASURED 2026-11-14 12:50**
     (`runs/20261114_scaling/`, spike binary, rtol 1e-8, same options at every rank count; exact
     from existing runs of the same case: 8 and 32 ranks `prof/`, 16 ranks `krylov4/`; host load
     shown, it was 10-20 and other sessions' jobs ran beside every run, so absolute times are
     indicative and the curve shape is the finding). s/step = loop / 3; s per iteration =
     kspsolve / (6 x iterations per solve): 2 CG solves per step, 6 `KSP Object` blocks per 3-step
     log (corrected 2026-11-14 audit; first posted with one solve per step, 2x too high); kspsolve
     is the in-code `EQ_WALL` line, loop and setup are `runInfo.json`.

     | ranks | option | setup s | loop s | s/step | its/solve (min-max) | s per iteration | peak RSS GB per rank | load before |
     |---|---|---|---|---|---|---|---|---|
     | 8 | GAMG CG | 18.5 | 138.7 | 46.2 | 63.0 (60-67) | 0.361 | 4.80 | 10.11 |
     | 8 | exact MUMPS | 640.4 | 38.2 | 12.7 | 1 | n/a | 20.91 | n/a |
     | 16 | GAMG CG | 10.6 | 98.4 | 32.8 | 63.0 (60-67) | 0.255 | 2.54 | 16.49 |
     | 16 | exact MUMPS | 331.3 | 20.6 | 6.9 | 1 | n/a | 12.07 | 13.27 |
     | 32 | GAMG CG | 5.5 | 55.4 | 18.5 | 63.2 (61-67) | 0.142 | 1.35 | 19.89 |
     | 32 | exact MUMPS | 255.6 | 21.2 | 7.1 | 1 | n/a | 6.79 | n/a |
     | 64 | both | not run: load1 was 30.8 at 12:47 (other sessions), 64 ranks would oversubscribe the 64 cores | | | | | | |

     Reading: GAMG CG scales and MUMPS does not. From 8 to 32 ranks GAMG's loop falls 2.5x (63 %
     parallel efficiency, ranks^-0.66; s per iteration halves per doubling, the iteration count stays
     at 63 so the cost is the work per iteration, not convergence), setup 3.4x, memory per rank 3.5x;
     exact's loop falls 38.2 to 20.6 s and then stays flat (21.2 s at 32), its setup
     falls 2.5x. The GAMG/exact loop ratio is 3.6x, 4.8x, 2.6x at 8, 16, 32 ranks. Extrapolating
     GAMG's ranks^-0.66 against exact's flat ~21 s puts the crossing near 140 ranks (an
     extrapolation from three noisy points, anchored at the 32-rank point: 139; a 3-point fit gives 147; not a measurement; the exact 16-rank baseline is 20.6 s in `krylov4/`, a second run `20261114_prof/` np16 gave 27.8 s, so "flat" rests on one of two 16-rank runs), beyond this node; GAMG's setup
     advantage (5.5 s against 255.6 s) is already 47x at 32 ranks and memory 5x, but the solve is
     inexact (9 of 33 files fail at 8 and 32 ranks, worst 1.1e-7, float32 floor) and over a 200-step
     run the loop, not setup, dominates (exact wide at 16 ranks: loop 1103 s for 200 steps), so the
     verdict at this size and rank count is unchanged. The 64-rank point, measured twice (matched-load block below),
     does not test the extrapolation: GAMG's loop does not fall past 32 ranks on this node.
   - **GAMG CG tuning on wide, 16 ranks, 3 steps, and the iteration-count control, MEASURED
     2026-11-14 13:05** (`runs/20261114_scaling/`, rtol 1e-8 unless stated, comparator = the table's,
     33 files; exact at 16 ranks: loop 20.6 s; load before each run 16-29 (13.7 for the gs10 control), so iterations carry the
     comparison and loops are indicative):

     | variant | its/solve (min-max) | setup s | loop s | peak RSS GB per rank | worst diff from exact | files failing |
     |---|---|---|---|---|---|---|
     | reference | 63.0 (60-67) | 10.6 | 98.4 | 2.54 | 1.1e-07 | 9 of 33 |
     | rtol 1e-6 | 42.2 (40-46) | 9.9 | 62.7 | 2.54 | 1.3e-06 | 15 of 33 |
     | rtol 1e-10 | 83.5 (81-87) | 9.9 | 122.6 | 2.54 | 3.4e-08 | 6 of 33 |
     | aggregation threshold 0.02 | 60.5 (58-64) | 42.9 | 133.5 | 2.34 | 1.1e-07 | 9 of 33 |
     | unsmoothed aggregation | 93.0 (90-97) | 7.2 | 125.8 | 2.54 | 1.1e-07 | 9 of 33 |
     | Richardson + SOR smoothers | 35.2 (34-37) | 10.0 | 80.8 | 2.54 | 1.1e-07 | 9 of 33 |

     Same ordering as on bp1002: the smoother (SOR, 1.8x fewer iterations) and the tolerance move the
     iteration count (about 10 iterations per decade of rtol here, 42.2 / 63.0 / 83.5, against 1.5 on bp1002); the threshold only
     raises setup (4x) and unsmoothed aggregation raises iterations 1.5x. No tolerance tested reaches
     the gold comparator (rtol 1e-10: 6 of 33 files still fail). SOR smoothers scale with ranks the same way:
     35.2, 35.2, 35.3 iterations at 8, 16, 32 ranks (rank-independent), loop 115.6 (load 28), 80.8, 50.7 s
     against the reference 138.7, 98.4, 55.4 s (SOR 8 to 32 ranks: 2.3x, ranks^-0.60; same ~140-rank
     crossing against exact's flat ~21 s, an extrapolation).
     **Matched-load pairs, wide, 3 steps, 2026-11-14 13:07-13:24** (MEASURED, indicative; `runs/20261114_matched/`,
     `scratch/bench/sub8.sh`, spike binary, GAMG CG with Richardson+SOR smoothers rtol 1e-8 run
     immediately before exact MUMPS, same case and rank count, one heavy-lock hold; loop = `runInfo.json`
     `time_loop_seconds`, RSS = maximum over rank files, load = 1-minute load before the run, which at 32
     and 16 ranks still contains the previous run of the pair; other jobs of the host owner, about 4
     cores, ran beside every run):

     | ranks | GAMG CG (SOR) setup / loop s | its/solve | RSS GB | load | exact setup / loop s | RSS GB | load | loop ratio |
     |---|---|---|---|---|---|---|---|---|
     | 16 | 9.7 / 60.4 | 35.2 (34-37) | 2.54 | 14.96 | 330.1 / 19.5 | 12.28 | 18.76 | 3.1x |
     | 32 | 5.4 / 40.5 | 35.3 (34-38) | 1.36 | 25.66 | 194.0 / 16.2 | 6.81 | 29.56 | 2.5x |
     | 64 | 8.3 / 52.2 (4.60); repeat `..._g64b` 8.3 / 47.1 (load 4.81; foreign jobs on the host not checked) | 37.2 (36-39) | 0.78 | 4.60 | first run aborted (-11, see below); rerun `e64b` 177.9 / 25.1 | 3.77 | 18.84 | 1.9x vs the repeat (47.1 / 25.1) |
     | 64, clean pair 2026-11-14 14:51-14:57 (`runs/20261114_clean64/`, `g64c`, `e64c`; load 4.77 before GAMG) | 7.4 / 50.7 | 37.2 (36-39) | 0.78 | 4.77 | 169.6 / 24.3 (ICNTL(14)=100, INFOG(1)=0) | 5.06 | 44.1 (contains the GAMG run) | 2.1x |

     Reading: 16 to 32 ranks GAMG's loop falls 1.49x (ranks^-0.58) against exact's 1.20x; at 64 ranks GAMG's
     loop does not fall: 52.2, 47.1 and 50.7 s in three runs (the last two started at load 4.8; all three above the
     32-rank 40.5 s, 1.16-1.29x). The 64-rank run has no spare core on this 64-core node (the host owner's other jobs, about 4
     cores, share cores with it, and a bulk-synchronous run waits for its slowest rank), so the 64-rank point is
     a node-saturated result, reproduced in three runs, not a clean scaling point; whether foreign cores were
     busy during `g64c` was not recorded. Exact's loop, same node: 19.5 s (16 ranks), 16.2 s (32), 22.8-25.1 s (64: e64c, e64b, w1-w3), so exact also slows from 32 to 64 ranks, by 1.5x; the GAMG to exact loop ratio at 64 ranks is 2.1x against 2.5x at 32. Earlier unmatched SOR 32-rank figure was 50.7 s at
     load 18-28; this run 40.5 s: run-to-run spread on this node is of order 25 %. Crossing re-derived from the
     matched 16/32 pair, an extrapolation: with exact flat at 16.2 s it is near 32 x (40.5/16.2)^(1/0.576) = 157
     ranks (the unmatched estimate was 139); if exact keeps falling at its measured -0.26 exponent the
     ratio shrinks only as ranks^-0.31 and the crossing is near 600 ranks. Either way it is beyond this 64-core
     node. A test would need a cluster run at 256 to 512 ranks (4 to 8 nodes of 64 cores, low-latency
     interconnect, exclusive nodes): GAMG needs about 0.15-0.26 GB per rank at those counts (measured per-rank RSS fell 1.87x 16 to 32 and
     1.74x 32 to 64 ranks; 1.74x per doubling extrapolated), exact MUMPS 6.8 GB per rank at 32 and 3.8 GB
     at 64 (e64b) with a setup of 178-330 s repeated; and the exact 64-rank runs at the default ICNTL(14) abort with `-11` in five of six tries (MUMPS -9 read in four of them; below). SOR accuracy against exact:
     see the 200-step block below.
     **`KSPConvergedReason = -11` reproduced at 64 ranks, exact MUMPS, wide; cause read: MUMPS INFOG(1) = -9**
     (2026-11-14; `runs/20261114_matched/..np64_e64` and `runs/20261114_minus11/..np64_r1` to `r4`; MUMPS
     ICNTL(4)=2 in r1-r4, default in e64). `log/cycle0.log` reads `par.solver=3: KSPSolve did NOT converge
     (KSPConvergedReason = -11) at step 1 iiTag 0`, then `MPI_ABORT`, rc 1, 117-147 s wall, no `runInfo.json`.
     Five `-11` aborts in six exact 64-rank runs at the default ICNTL(14) (the sixth, `..._e64b`, load 18.8, passed: setup 177.9 s,
     loop 25.1 s, RSS 3.77 GB). The MUMPS code was read in four of the five, r1-r4 only (ICNTL(4)=2): `run.log` prints, once
     per run, `** ERROR RETURN ** FROM DMUMPS INFO(1)= -9` with `On return from DMUMPS, INFOG(2)=` 4,799,232 / 8,804,430 /
     5,913,752 / 6,650,676 (and `INFO(2)` equal to that value on the failing rank), and on the other 63 ranks
     `INFO(1)= -1` with `INFO(2)` = 48 (r1, r3, r4) or 0 (r2). On those 63 ranks -1 means "an error occurred on another
     process" and `INFO(2)` is that process's rank (48, 0, 48, 48), not a quantity: it is constant over 63 ranks and
     small, which fits a rank id and cannot be a size. The MUMPS user guide could not be opened here (the conda env ships
     headers only; the guide URLs tried returned 404), so the meaning of -1 and of -9 (main internal work array too small,
     INFOG(2) the size short, cure a larger ICNTL(14)) is from the guide as recalled, not re-read; the log pattern above
     is the only local evidence. `..._e64` ran at the default print level: its logs hold no `INFO(1)` line, only the `-11`,
     so it is a fifth `-11` with its MUMPS code unread, not a fifth -9.
     So the failure is factorization workspace (ICNTL(14), default 20 % over the symbolic estimate) at 64-way
     parallelism, not a PETSc setup defect; whether the earlier single 8-rank sighting (1 of 11, then 0 of 100)
     had the same cause is unknown, its codes were never printed. Load before the failing runs was 35-68
     (their own predecessor included) and 18.8 before the one that passed: a reading, not a controlled
     comparison. Test of the cure, `-mat_mumps_icntl_14 100` (with ICNTL(4)=2) at 64 ranks, `runs/20261114_minus11/..np64_w1` to `w3`: all three passed (rc 0 in `rc.txt`; MEASURED 2026-11-14 14:46-14:55 CST, `DONE_w` present): loop 24.5 / 23.3 / 22.8 s (`runInfo.json`; `EQ_WALL` setup 191.9 / 179.9 / 183.5 s, e64b by the same field 189.6 s: the 177.9 and 179.6 s figures above came from another timer not re-identified), max-rank RSS 4.50 / 4.99 / 4.49 GB against e64b 3.77 GB (+19 to +32 %), MUMPS INFOG(1) = 0 in each. Against the default ICNTL(14): 0 of 3 failed with 100, 5 of 6 failed at default (4 of 4 with the same ICNTL(4)=2); Fisher exact one-sided p about 0.05, and load was 44-68 before the ICNTL(14)=100 runs against 35-68 before the failing ones, so it is a strong indication, not a controlled proof. Cause (with the guide caveat above) and runtime cure: factorization workspace at 64-way parallelism, `-mat_mumps_icntl_14 100`; cost about +20-30 % rank memory. The earlier 8-rank sighting is not shown to be the same cause: its code was never printed, and 0 of 100 later 8-rank runs failed at the default. The
     board's first closing condition (KSPSetUp reading INFOG(1)/(2) on failure, code in
     `solveTimeLoopMUMPS.f90`-style) is not done: the codes above came from ICNTL(4)=2 printing, not from code.
     **Slip-rate change, GAMG CG (SOR, rtol 1e-8, previous-step initial guess), wide, 32 ranks, 200 steps,
     MEASURED 2026-11-14** (`runs/20261114_event/bp8.qdc.gs.10.wide_s3_np32_ev32`, `-ksp_monitor`; rc 0,
     load 62.9 before (my own preceding runs) and 35.6 after, loop 2035.7 s = 10.2 s/step, setup 5.5 s).
     Iterations per solve counted from the monitor's iteration-0 resets: 400 solves = 2 per step x 200,
     mean 27.28 (agrees with `runInfo.json` `ksp_iterations_mean` 27.275), range 25-38. By 20-step block:
     32.15, 29.35, 28.3, 27.4, 27.0, 26.25, 26.0, 26.0, 25.3, 25.0 (re-counted from `run.log` 2026-11-14: 400 solves, mean 27.275, range 25-38; first block holds the maximum). The max
     slip rate of this case rises from 1e-12 to 1.479e-7 m/s in these 200 steps (final value identical to
     0.1479348E-06 in the exact 200-step run of `runs/20261114_blr-widelong/..np16_base`: different run
     directories, ranks (32 against 16) and solvers, so not a wrong-file read; the comparator below is the
     second statistic), so the 200 steps do cross five orders of magnitude of slip rate, and iterations per
     solve fall with the warm start rather than rising: there is no iteration growth through the change at
     this sampling (a per-step slip-rate trajectory was not recorded; `log/cycle0.log` prints it at
     output steps only; the run is aseismic-to-accelerating, not an event with rupture). Output against the
     exact baseline (repo comparator, `compare_series` / `compare_netcdf`, rtol 1e-9, 33 compared files):
     13 of 33 fail; worst per-file max|a-b|/max|a| 1.1e-4 (`fltst_strk+000dp+000.dat`, and 1.0e-7 to 3e-8
     for netCDF `disp`/`fault` fields: float32 floor). Most station-file failures are entries near zero
     (first message: row 200 column 7, `2.466e-10` against the reference; two files fail on non-zero values, see the matched-pair table below), not diagnosed
     further. So the SOR variant is inexact over 200 steps as it is over 3, and not by a margin that the
     tolerance of the gold gate covers.
     Slip-rate-change summary (200 steps, wide, 32 ranks, MEASURED 2026-11-14 15:11-15:41 CST, matched pair back to back,
     `runs/20261114_event2/`, `scratch/bench/sub13.sh`; the earlier unmatched `..._ev32` and 16-rank `np16_base` runs
     are superseded by it):

     | solver | its/solve through the change | setup / loop s | s/step (loop / 200) | host load before | accuracy against the exact run of the same ranks |
     |---|---|---|---|---|---|
     | GAMG CG (SOR), `ev32b` | 27.275 mean, 25-38 (identical to `ev32`'s count, a deterministic solve; loops differ, 2035.7 s there, so not a wrong-file read) | 5.6 / 1747.7 | 8.74 | 67.9 | 13 of 36 compared files fail (rtol 1e-9 netCDF, 2e-6 station files); worst netCDF 1.1e-7 of field scale (`disp.00001.nc`, float32 floor), station files fail mostly on near-zero entries (`fltst_strk+000dp+000`, row 200, column 7: 2.466e-10 against 4.19e-17; `srfst` rows 166 and 199, column 4; `cplot_EQquasi.txt` row 5099, column 6) and two at 5.7e-4 to 1.4e-3 relative on non-zero values (`fltst_strk+200dp+000` row 106, `fltst_strk-200dp+000` row 166, column 7) |
     | exact MUMPS (ICNTL(14)=100), `ex32` | 1 (preonly) | 191.8 / 725.8 | 3.63 | 36.5 (contains the GAMG run's tail) | reference |

     Loop ratio 2.4x (GAMG slower), against 2.5x in the 3-step matched pair at 32 ranks: no change in the ratio through the
     slip-rate rise (1e-12 to 1.5e-7 m/s), where GAMG's iterations per solve fall from 32.15 to 25.0 per 20-step block.
     Exact's setup is 191.8 s against GAMG's 5.6 s, but at 5.1 s/step less in the loop (8.74 - 3.63) exact is ahead
     after 37 steps, so GAMG does not win on total time past that. Comparator: `compare_series`/`compare_netcdf` of
     `tests/e2e/cases.py` on `result/cycle0`, one-off script, not committed.

     Spike-run reporting defect, 2026-11-14 (audit): `runInfo.json` of the GAMG runs reads `ksp_type` preonly, `pc_type`
     cholesky (g64b checked), not the solver used; the real iteration counts (`ksp_iterations_*`) are right. Merge-bound
     solver 3 path: `solveTimeLoopPETSc.f90:496-509` sets the names from the whitelist booleans (`isPreonly`,
     `isCholesky`) and the guard at `:472` rejects other combinations, so it names what it ran by construction; read from
     source only, the merge-bound binary was not run on a GAMG option. Not changed (source change, owner). Owner
     questions, not decided here: a test pinning the INFO(2) meaning (new test); `KSPSetUp` reading INFOG(1)/(2) on
     failure and any default ICNTL(14) (changes merge-bound solver behaviour).
     **Why wide needs 63 iterations and bp1002 10:** not wide-specific. The control `bp8.qdc.gs.10`
     (1.96 M equations, no enlarged box) at 16 ranks takes 64.0 iterations per solve (62-67), the same
     as wide's 63.0; GAMG CG on it: setup 5.6 s against 80.8 s exact, loop 49.3 s against 8.1 s (6.1x),
     RSS 1.44 against 4.47 GB, worst diff 1.1e-07, 9 of 33 files fail. The AMG hierarchy has 5 levels on
     wide against 4 on bp1002 with operator complexity ~1.15 on both, so the growth is not coarse-grid
     cost; what separates bp1002 (163 k equations, 10 iterations) from the BP8 meshes is size or the
     BP8 fault geometry, which these two points cannot separate (a mesh-refinement series of one geometry
     would; not run). Not examined: the near-null space quality on fault split nodes.
   - **Tight tolerances pass; BLR as a CG preconditioner does not win (2026-11-14, conductor 11; MEASURED, indicative, host load 6-12 before each bp1002 run; spike binary `spike-sub` 7d469dc, rebuilt in a fresh worktree, runtime options only: `-eq_plain_ksp -ksp_type cg -eq_guess_per_tag`, `-pc_type cholesky` with MUMPS and `-mat_mumps_icntl_35 1 -mat_mumps_cntl_7 <eps>` for BLR, `-pc_type gamg` with near-null and SOR for GAMG).** Comparator: the table's, run against an exact run of the same binary, rtol 1e-9 on every compared file; NOT the `e2e_fast` gold gate, which was not run on these options. bp1002, 8 ranks, 101 steps (`runs/20261114_blr-pc/`, `runs/20261114_gamg-tight/`; exact beside them: setup 3.21 s, loop 11.42 s, 0.494 GB per rank):
     | option | KSP its/solve (max) | setup s | loop s | peak RSS GB per rank | files failing (of 25) |
     |---|---|---|---|---|---|
     | BLR 1e-8 PC, CG rtol 1e-12 | 2.00 (2) | 3.32 | 21.4 | 0.501 | 0 |
     | BLR 1e-6 PC, CG rtol 1e-12 | 2.18 (3) | 2.81 | 23.5 | 0.463 | 0 |
     | BLR 1e-4 PC, CG rtol 1e-12 | 4.21 (5) | 2.47 | 31.8 | 0.445 | 0 |
     | BLR 1e-4 PC, CG rtol 1e-10 | 3.32 (4) | 2.53 | 24.9 | 0.440 | 5 (worst 3.8e-9, `disp`) |
     | GAMG CG (SOR) rtol 1e-11 | 11.12 (13) | 0.95 | 40.9 | 0.319 | 0 |
     | GAMG CG (SOR) rtol 1e-12 | 12.21 (14) | 0.94 | 44.2 | 0.319 | 0 |
     So the two cheap "not run" items on GAMG CG are answered at this size: rtol 1e-11 and 1e-12 pass the 1e-9 comparator (rtol 1e-10 failed 3 of 25 on 11 steps at 16 ranks, a different case); and BLR 1e-4, which fails as a direct solve, passes as a preconditioner under tight CG. Neither is faster than exact on bp1002 (loops 1.9-2.8x and 3.6-3.9x exact). Wide, 16 ranks, 10 steps (`runs/20261114_blr-pc-wide/`, host load 12-22, ICNTL(14)=100 throughout): exact setup 317.7 s, loop 53.2 s (5.32 s/step), 13.46 GB; BLR 1e-8 PC, CG rtol 1e-12: setup 139.6 s, loop 54.1 s (5.41 s/step), 2.00 its/solve, 10.99 GB, 3 of 33 files fail (worst 1.4e-8 of field scale, `disp.00001.nc`); BLR 1e-6: setup 116.0 s, loop 77.4 s (7.74 s/step), 3.70 its, 10.65 GB, 2 of 33 fail (worst 1.2e-8). So no loop speedup on wide: two iterations of a solve that is 1.9x faster than exact is parity, and it does not pass 1e-9 there (netCDF `disp` at 1e-8 of field scale, below the float32 one-ulp floor of 6e-8 but above the comparator). Gains are setup (2.3-2.7x) and 18-21 % memory. GAMG CG at tight rtol on wide, the near-null space on fault split nodes, and the gold gate on these options were not run.
     **Controls on the wide BLR-PC result, 2026-11-14 16:30 (conductor 11, MEASURED; coordinator's suggestions).** (1) Noise floor: exact against exact on wide passes the same comparator, 33 of 33 files at 1e-9, for a same-rank rerun (`krylov4` against `blr-wide`, 16 ranks, 3 steps) and across rank counts (16 against 32, `blr-wide` against `blr-wide2`); at 200 steps (16 against 32 ranks, `blr-widelong` base against `event2` ex32) the netCDF files pass and only three `fltst` `.dat` log10 slip-rate columns fail (the comparator applies rtol 2e-6 to series files; worst -26.095 against the floor sentinel -30 at row 132, column 5; the column the consolidated table excludes). So the comparator floor on wide is below 1e-9 and the 1.4e-8 `disp` differences are solver error. (2) Stopping test: the runs already used `-ksp_norm_type unpreconditioned`; `-ksp_monitor_true_residual` (`runs/20261114_blr-pc-wide2/`, wide, 16 ranks, 3 steps, BLR 1e-8 PC) shows the monitored and true residuals identical, and at CG rtol 1e-12 (2 iterations, ||r||/||b|| 4.0e-13 at step 1) `disp.00001.nc` still differs from exact by 1.388e-8 of field scale (1 of 33 files fail), i.e. a relative residual of 4e-13 leaves a displacement error ~3e4 times larger (an ill-conditioning amplification, read from these two numbers, not a measured condition number). At rtol 1e-15 (3 iterations, 1.0e-18) 33 of 33 files pass. So BLR-preconditioned CG on wide can pass the 1e-9 comparator, at one more iteration: loop 25.5 s (rtol 1e-15, 3 iterations) against 17.8 s (rtol 1e-12, 2 iterations) for the 3 steps (host load not controlled), i.e. 8.5 s/step against exact 5.3 s/step at 10 steps, so it still does not beat exact on time. This also changes how "GAMG CG at rtol 1e-12 passes on bp1002" reads for wide: a tolerance that passes on one case is not evidence for another, and GAMG CG on wide was not run below rtol 1e-10 (6 of 33 failing) until the next paragraph. (3) Per application: `EQ_WALL` kspsolve over the 10-step wide runs, 20 solves: exact 49.9 s = 2.50 s per application; BLR 1e-8 PC 51.3 s over 40 applications = 1.28 s per application (1.95x cheaper, matching the direct BLR 5.52 against 2.87 s/step, 1.92x). The tie of 5.41 against 5.32 s/step is therefore the iteration count (2), not a slower BLR application; the gain over exact is setup and memory only.
     **GAMG CG on wide at tight rtol, 2026-11-14 16:55 (conductor 11, MEASURED, indicative; `runs/20261114_gamg-wide-tight/`, spike binary, 16 ranks, 3 steps, near-null + Richardson/SOR + previous-step guess, `-ksp_norm_type unpreconditioned`; compared with the 3-step exact run `20261114_blr-wide/..np16_base` at the table's comparator, 33 files).** rtol 1e-12: 58.0 iterations per solve (max 60), setup 10.0 s, loop 122.4 s (load 5.2), 2 of 33 files fail (`disp.00003.nc`, `disp.r.nc`, worst 8.4e-9 of field scale); rtol 1e-14: 68.7 (70), loop 149.8 s (load 17.8), 0 of 33 fail; rtol 1e-15: 74.2 (76), loop 160.2 s (load 18.9), 0 of 33 fail; peak RSS 2.54 GB per rank in all three (exact 12.1-13.5 GB). So an iterative solve passes the 1e-9 comparator on wide at rtol 1e-14, at 41-50 s/step against exact 5.3-6.9 s/step (7-9x slower; two exact baselines exist, 6.0-6.4x against 8.32 s/step and 7.3-7.8x against 6.88, see conductor 12 below), setup 10 s against 318-336 s and 4.8x less memory; the two tolerances that pass differ by one decade from the one that fails, consistent with the displacement amplification above. Not verified: that the pass holds over 200 steps or through the slip-rate rise, the gold gate on these options, any rank count other than 16.
     Stack check, 2026-11-14 16:30: #76, #77, #82, #83, #85 show MERGEABLE and CLEAN with CI success (`gh pr view`); #87 is a draft spike with CI success. None needed a rebase.
   - Open (2026-11-14 12:15): `-11` 60-run batch FINISHED 2026-11-14 12:34, 60 of 60 clean (`runs/20261114_ksp-stress2/results.txt`: `rc=0`, `rc=0` in all 60 (an abort returns rc 1, as the 64-rank sighting did; the `-11` text goes to `log/cycle0.log`, which the batch's grep on the run log did not read), 60 of 60 logs read `INFOG(1) = 0` and none negative (re-checked 2026-11-14 with the pattern `INFOG\(1\) *= *-`; the batch script's own INFOG pattern had no space before `=` and could not match); `ALLDONE` in `rc`; binary 1.24.4, 8 ranks, `test.bp8.qdc.gs.10`, 3 steps, solver 3, MUMPS ICNTL(4)=3): the row's second closing condition (60 clean runs) is met; with the 0 of 40 earlier and 1 of 11 reported, the `-11` stays an unreproduced single observation (closing the row is the board owner's); sub4 and the bp5 GAMG CG run are harvested (table rows; 12:35). Wide subdomain-direct, 1 block per rank (`runs/20261114_subdirect/`, MEASURED, host load 13.5, indicative): setup ~0, loop 1000.6 s against 20.6 s exact (48x slower), 229 CG iterations per solve (max 247), peak RSS 3.36 against 12.07 GB per rank (3.6x less), worst difference 1.1e-7 (`disp.00001.nc`, float32 floor), 9 of 33 files fail: it loses on time by a wider margin than GAMG CG (97 s) and is inexact. bp1002 at 16 ranks: loop 15.7 s against 0.94 s exact (16.7x), 66 iterations, 14 of 25 fail.
     **Conductor 12, 2026-11-14 16:55 CST (MEASURED where marked; the rest read from source, not run).**
     *Full `e2e` tier on the stack head* (`dd-2d` c0a3160, `par.solver = 3` forced, `scratch/bench/gate_full.sh`): running since 16:40, result not in yet; it goes on the board in the next PR with counts and failures by name. Nothing is pushed to the stack branches until it ends.
     *Code audit of the full stack diff (`origin/master...origin/dd-2d`), desk read only, nothing built or run: no BLOCKER, no MAJOR, 9 MINOR (line numbers on `dd-2d`).* (1) `docs/user/parameters.md:131` and `scripts/defaultParameters.py:161` say "no bp 8" for solver 3, but bp 8 runs and is tested. (2) `faulting.f90` stops 501/502 and the solver-2-without-PETSc stop (`eqquasi.f90:~103`) became `MPI_ABORT`: exit status now comes from the launcher (behaviour change on solvers 1/2). (3) `nsmp = 0` (`eqquasi.f90:~303`) makes the padding columns of `eqquasi.mesh.nsmp.nc` 0 instead of uninitialised when faults differ in size (output change on solvers 1/2). (4) The solver-3 restart path (`netcdf_io.f90:~671/761/805`) is run by no test; every DD test starts at cycle 0. (5) Every solver-3 test skips on a binary built without PETSc. (6) `tests/e2e/test_benchmarks.py:~33-35` detects PETSc through `ldd`: a statically linked PETSc reads as absent and every DD test skips. (7) The pure-Python comparator tests (`test_benchmarks.py:~807-883`) are marked `e2e_fast`, so neither the fast tier nor CI runs them (an existing owner question: move them to `tests/unit`). (8) The 4-byte node-count guard at `solveTimeLoopPETScDD.f90:128` runs after the counts are already made in 4-byte integers. (9) `netcdf_io.f90:~356-365` zero-fills each `fault.*.nc` with many separate `nf90_put_var` calls per snapshot (cost, not correctness). Found clean: collectives and deadlocks, fault-node ownership and halo exchange, node numbering, integer overflow at 3.6 M equations (the first 4-byte limit is the MPI-IO write count past 4.17 M owned fault nodes on one rank), solvers 1/2 otherwise unchanged, build without PETSc, version 1.24.4. Not read: `faulting.f90` outside the hunks, the `exitCriteria` body, the profile-output routines, `solveTimeLoopPETSc.f90`. Plan after the tier result: fix 1, 7 (as far as the owner question allows: a marker change in an existing file, not a new file) and 8 on the owning branches in existing files, re-gate; the rest wait for the owner.
     *How a user passes `-mat_mumps_icntl_14 100` on the merge-bound solver 2/3 paths today (read from source, the merge-bound binary not run with it):* through the PETSc options database, with no source change. `PetscInitialize(PETSC_NULL_CHARACTER, ...)` (`solveTimeLoopPETSc.f90:126`, `solveTimeLoopPETScDD.f90:109`) reads `PETSC_OPTIONS`, `~/.petscrc` and `./.petscrc`; `KSPSetFromOptions` runs after `PCFactorSetMatSolverType(..., MATSOLVERMUMPS)` (`:428`, `:549`), so the MUMPS options are read at factorization. `export PETSC_OPTIONS="-mat_mumps_icntl_14 100"` before `bash run.sh` is enough on one node (the generated `run.sh` passes the launcher the inherited environment; for several nodes add `-x PETSC_OPTIONS` through `MPIRUN_ARGS`). The solver-2 whitelist guard (`solveTimeLoopPETSc.f90:~449-472`) rejects other KSP/PC type combinations, not MUMPS control parameters. Every ICNTL(14) run on this row used exactly that route (`PETSC_OPTIONS` set from `SPIKE_OPTS` in `scratch/bench/runcase_sub.sh:8`, spike binary). `docs/user/` does not mention it (docs gap, owner's call whether to document). A default ICNTL(14) in the source stays the existing owner question above.
     *Independent audit of #114's wide GAMG table (quantitative-claim auditor, 19 tool calls; logs read, comparator re-run on the 33 files):* iterations per solve (58.0/60, 68.7/70, 74.2/76), setup 10.0 s, loops 122.4/149.8/160.2 s, host loads, 2/0/0 of 33 failing (`disp.00003.nc`, `disp.r.nc`, 8.4e-9), RSS 2.54 GB and the 4.8x memory ratio all recompute. The EXACT baseline figures did not: against the stated baseline `20261114_blr-wide/..np16_base` (`runInfo.json` `time_loop_seconds` 24.96 s over 3 steps, 343.6 s setup, 12.25 GB per rank) the correct numbers are 8.32 s/step, setup 343.6 s, 12.25 GB; the board's exact figures came from a different exact run of the same case (`krylov4/..np16_exact`: 6.88 s/step, setup 331.3 s, 12.07 GB, table row above; conductor 12 first reported no matching file, wrongly, and corrected here on re-reading the table). Corrected reading: the two passing runs are 49.9 s/step (rtol 1e-14) and 53.4 s/step (1e-15), so 6.0x and 6.4x slower than exact on wide at 16 ranks (the failing rtol 1e-12 run, 40.8 s/step, is 4.9x), setup 10 s against 343.6 s, memory 4.8x less. Also: the 41-50 s/step range spanned 1e-12 to 1e-14 and included the failing run; the repo's `CATEGORIES` has no `disp.*` entry, so `pytest` alone would not compare the two files that fail at 1e-12 and the 33-file harness that produced the table is not in the repo; the GAMG `runInfo.json` records `preonly`/`cholesky` although the solves were CG (the log's `EQ_KSP` lines are the evidence of the solver used); RSS is `rss.*` kB/1000, the log's largest rank reads 2.42 GiB. Whether the exact wide 3-step runs elsewhere on this row (5.3-6.9 s/step) are other loads or other rank counts was not traced; the ratio is 6.0-6.4x against `blr-wide/..np16_base` (8.32 s/step) and 7.3-7.8x against `krylov4/..np16_exact` (6.88 s/step); two exact runs of one case differ by 21 % in loop time, so read it as 6-8x.
     **Full `e2e` tier on the stack head, 2026-11-14 18:22 CST (conductor 12, MEASURED; `runs/20261114_gate-full/{full.out,junit.xml,full.rc,log.txt}`, `scratch/bench/gate_full.sh`).** `dd-2d` c0a3160 (worktree detached there, binary built from it, host `cotopaxi`, `MACHINE=conda-linux`), `par.solver = 3` forced through a throwaway edit of `tests/e2e/test_benchmarks.py` in that worktree, `python3 -m pytest -q -rs -m e2e tests`: rc 0 (`full.rc`), 90 passed, 45 skipped, 0 failed, 0 errors, 287 deselected, 6098 s (16:40:35-18:22:14). Skips by reason (`-rs`): 10 x no reference at `data/test.bp1002.qdc.zone.2500/cycle0` (that gold does not exist anywhere: `tests/e2e/cases.py:96` records the case as pending its reference, so these ten are not exercised, not a provisioning artifact; `data/` in the main checkout and the worktree list the same nine directories); 1 no restart files, 4 no offfault, 8 no perfault, 9 no profile files, 4 no summary.json, 9 no section 4.3 profiles, each "this reference holds none" for the category in that gold; none is a missing PETSc or a missing binary. The tier ran with solver 3 forced, so it exercises the DD path on every case it has a gold for; solver 1/2 on the stack head are covered by CI's tier, not by this run. Docs-only commit 0bc146b (audit item 1: "no bp 8" removed from `docs/user/parameters.md` and `scripts/defaultParameters.py`) was then pushed to `dd-2d`; no `src/` change, so this result stands for it. Audit item 7 cannot be done by a marker change: `tests/conftest.py:111` adds the tier marker from the directory, so only moving the comparator tests under `tests/unit` (a new file, owner) takes them out of the e2e tier; item 8 needs a `src/` version bump and a restack of #76/#77/#82/#83/#85, not worth it before the owner decides on the stack. Both listed for the owner with items 2-6 and 9.
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
