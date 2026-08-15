# reference/stepover.con

First reference for `test.stepover.qdc.con.1000` — the constraining twin of
`test.stepover.qdc.1000`. Regression lock, not a verified benchmark.

The ONE flip from the releasing twin is the sign of `par.far_vel_load`
(-4e-10: left-lateral). Identical geometry, so over cycles the interior tips
are EXPECTED to clamp instead of unclamping (the stop-508 mode of the
releasing case). Not exhibited in this 101-step lock — see the flag below:
the window is not yet left-lateral, only dip components and stress
magnitudes respond.

- Frozen 2026-08-15 from `work/deorphan.con` (`create.newcase` +
  `case.setup` + `bash run.sh`, unmodified). Binary `eqquasi-1.16.0`, knox,
  2 MPI ranks. cycle0: 101 steps, interseismic (V decays 1.0e-9 -> 6.0e-10).
- KNOWN FLAG (not fixed here): the compset's initial shear is
  magnitude-only — `shear_steady_state` does not carry the loading sign, so
  the left-lateral run starts slightly off steady state; the V decay above
  is that transient. Deterministic, so valid as a lock; flagged for the
  owner as the laterality-switch design question.
