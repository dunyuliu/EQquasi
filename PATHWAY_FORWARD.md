# Pathway forward

Open work only. Updated 2026-10-08, at v1.21.2.

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
   **(a) done, 2026-10-09** (conductor, binary eqquasi-1.21.2, HEAD 605eacf,
   theo4 host-exclusive, fresh builds, 10-step profiling runs, native MUMPS
   `par.solver=1`): bp1002.qdc.2500 (127k eq) -- 8r 0.436 s/step, 16r 0.409,
   32r 0.394 s/step (flat; rank-0-only faulting/pressure work dominates at
   this size). BP8 10 m wide (3.52M eq / 1.23M nodes, the scratch-only
   `fxmin/fxmax/fzmin/fzmax=-1500/1500`, `enlarging_ratio_xz=1.2` variant, NOT
   the committed `compset/bp8.qdc.gs.10` which is a smaller `-500/500` box,
   694k nodes) -- first sweep (8r 10.63 s/step, 16r 7.98, 32r 6.76 s/step)
   **discarded as contaminated**: row 2's agent ran heavy serial compute
   (up to 19 GB/11 min) on the same host in that exact window. Clean rerun
   in progress as of this entry; see the next evidence line or runInfo.json
   under `runs/20261009_petsc-profile/bp8wide-r{8,16,32}/` for the real
   numbers. (b)/(c)/(d) not started either way. Finer solve /
   faulting+pressure / communication / I/O split (as (a) asks) is NOT yet
   instrumented -- only aggregate time-loop and factorization time exist
   today (`solveTimeLoopMUMPS.f90`'s `cpu_time` calls). Runs + logs +
   runInfo.json in `runs/20261009_petsc-profile/{bp1002,bp8wide}-r{8,16,32}/`.
   Recommendation (not acted on -- rule 1, new-compset naming needs
   agreement): register the -1500/1500 wide box as its own compset (e.g.
   `bp8.qdc.gs.10.wide`) instead of leaving it only in `scratch/`.
2. [ ] **P2 -- FEniCSx feasibility spike.** Owner, 2026-10-09. Not a port: a
   one-fault quasi-static prototype (BP5 or BP8) outside src/, in runs/, to
   answer three questions before any port is considered: can split-node (or
   interface) fault slip with rate-and-state be done cleanly; do P2 elements
   close the BP8 gap to Kim that linear hexes leave (depth -18 % at 50 m,
   -13 % at 10 m); what does a step cost against EQquasi on the same problem.
   Uses the host's fenicsx conda env; no new file in the repo until the
   answer says so (rule 1). Done when: a one-page verdict with numbers, under
   docs/dev/, and a go / no-go for a port.
   **Q1/Q2 done, 2026-10-09** (dunyu-liu): verdict landed at
   `docs/dev/fenicsx_feasibility_2026-10-09.md` (PR #52, merged 45ab982).
   Split-node RSF works cleanly in dolfinx. P2 elements do NOT close the
   BP8-vs-Kim depth gap -- the gap is EQquasi-specific (an independent
   FEniCSx linear-hex model converges to Kim; EQquasi's own 50m/10m runs sit
   at -18%/-13% and get worse with refinement at the centre station) --
   **new finding, routed to lars-eriksson, not row 2's to fix**. Lean verdict
   so far: no-go on a port for the depth-gap motivation. **Q3 (step cost vs
   EQquasi) still open**, blocked on an idle host (theo4 had an unrelated
   external job running as of 13:52; row 1's own timing needs it too).
   Not done yet: row 2 isn't closed until Q3 has a number.
