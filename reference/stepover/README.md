# reference/stepover

First reference for `test.stepover.qdc` (two-fault step-over, ntotft = 2,
101 steps). Regression lock, not a verified benchmark.

Why this row matters more than its size: it is the only ntotft > 1 case cheap
enough for the fast tier, so it puts per-fault routing of on-fault input into
every-push CI for the first time — the gap PROJECT_RULES rule 12 records
(BP1002 covers it only in the full tier). Three multi-fault bugs in this
project's history were found by hand for lack of exactly this row; the
per-fault `cplot_ft2_*` files frozen here are the artifacts that catch the
next one.

- Frozen 2026-08-15 from `work/deorphan.stepover` (`create.newcase` +
  `case.setup` + `bash run.sh`, unmodified).
- Binary `eqquasi-1.16.0`, knox, 2 MPI ranks (`HPC_ncpu = 2`).
- cycle0: 101 steps (nstep cap, by design), interseismic throughout
  (peak Vmax 1.0e-9): the gate locks the creeping stress build-up state,
  including both faults' cplot output and all six station files (three per
  fault, fault 2 under the `ft2_` tag).
- The compset's original stations matched no node — both fault corners sit
  at x = 500 mod 1000 m, so integer-km x coordinates miss the lattice (the
  exact silent-station failure the eqquasi.f90 warning was written for; here
  it hit all six). Station x moved to half-km values in the same change.
