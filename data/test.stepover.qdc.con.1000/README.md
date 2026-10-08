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
- **RE-BLESSED 2026-09-24** (owner decision "go as recommended", the rule-8
  exception, not a precedent): the KNOWN FLAG above is fixed. Initial shear
  now follows `sign(par.far_vel_load)` in
  `compset/test.stepover.qdc.con.1000/user_defined_params.py`, so the
  left-lateral case starts close to steady state instead of off it. Refrozen
  from `create.newcase` + `case.setup` + `bash run.sh`, unmodified, binary
  `eqquasi-1.19.0`, theo4, 2 MPI ranks: cycle0 101 steps, V decays
  1.0e-9 -> 8.07e-10 (was 1.0e-9 -> 6.0e-10) — a smaller transient, consistent
  with starting nearer steady state, not perfectly flat (only dip components
  and stress magnitudes fully respond to the sign flip; the residual decay is
  the same kind of settling the releasing twin's positive-sign case also
  shows over its first few steps, not a reintroduction of the old bug).
- **Additive re-bless 2026-09-25** (row 6, rule 8): `peak_sliprate_per_fault.dat`
  added (time, then one peak-V column per fault; column count = 1 + ntotft).
  Every pre-existing oracle file was first compared against the fresh run
  (`create.newcase` + `case.setup` + `bash run.sh`, binary `eqquasi-1.19.1`,
  theo4, 2 MPI ranks) with the e2e comparators: all pass, worst 1.8e-13 of
  scale (cplot_ft2), 1.1e-15 (global.dat), 2.9e-15 (fault.*.nc). Nothing else
  was replaced. max over the two fault columns equals global.dat column 2
  exactly.
