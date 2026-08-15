# reference/stepover.con

First reference for `test.stepover.qdc.con.1000` — the constraining twin of
`test.stepover.qdc.1000`. Regression lock, not a verified benchmark.

The ONE flip from the releasing twin is the sign of `par.far_vel_load`
(-4e-10: left-lateral). Identical geometry, so the step that dilates under
right-lateral slip is compressional here — over cycles the interior tips
clamp instead of unclamping (the stop-508 mode of the releasing case).

- Frozen 2026-08-15 from `work/deorphan.con` (`create.newcase` +
  `case.setup` + `bash run.sh`, unmodified). Binary `eqquasi-1.16.0`, knox,
  2 MPI ranks. cycle0: 101 steps, interseismic (V decays 1.0e-9 -> 6.0e-10).
- KNOWN FLAG (not fixed here): the compset's initial shear is
  magnitude-only — `shear_steady_state` does not carry the loading sign, so
  the left-lateral run starts slightly off steady state; the V decay above
  is that transient. Deterministic, so valid as a lock; flagged for the
  owner as the laterality-switch design question.
