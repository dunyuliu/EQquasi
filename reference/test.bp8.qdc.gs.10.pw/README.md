# reference/test.bp8.qdc.gs.10.pw

FIRST, UNVERIFIED reference for BP8 with Peaceman-well fluid loading
(`par.fluid_src = 2`) on the `test.bp8.qdc.gs.10` compset, gated only in the
full e2e tier (`-m e2e`), never fast/CI — a 30-day-simulated run is not a
per-push cost (PATHWAY_FORWARD.md row 15, owner decision 2026-09-24).

`name` != `compset` in `testsys/e2e/cases.py`'s CASES row on purpose: this
reference is reused off the GS compset (`fluid_src=2` override only —
`fluid_q0`/`fluid_rwell` already default to the 2026-08-13 spec values,
0.0015 and 0.05, in that compset's own `user_defined_params.py` since row 10,
PR #22/v1.18.3) rather than a new compset directory, so as not to rename
`bp8.qdc.gs.10` and orphan its own read-only GS reference.

- Frozen 2026-09-24/25 from `create.newcase` + `case.setup` + `bash run.sh`,
  unmodified. Binary `eqquasi-1.19.0`, theo4, 3 MPI ranks. Exit on the
  solver's own criterion (`fluid_tend`, not `nstep` — `nstep=8000` is only a
  safety cap above the observed exit step): 5196 steps, simulated 30.0026
  days, final max slip rate 3.8167e-10 m/s.
- Sanity-checked (not diffed byte-for-byte — rule 11, different binary) against
  a prior run at `scratch/bp8.pw.spec0813` (binary `eqquasi-1.18.2`, predates
  the v1.18.0-1.18.2 PETSc solver landing): that run also exited at exactly
  5196 steps with final max slip rate 3.82e-10 m/s, peak well-cell pressure
  14.26 MPa, slip 0.0201 m. Step count and final Vmax both agree to 3
  significant figures across the two binaries — no reason to suspect the
  PETSc-solver work changed this case's physics.
- Rule 8: this is a FIRST reference — no prior run of THIS binary/version
  exists to diff against, so "verified" here means the physical sanity check
  above, not a zero-diff regression history. The next run against it is the
  first real regression check.
- **OPEN QUESTION, found 2026-09-25, not decided here**: re-running this exact
  case against this exact reference (`pytest -m e2e -k gs.10.pw`, same binary,
  same 3 ranks, immediately after freezing it) FAILS the `onfault` category at
  4 of 9 stations — `fltst_strk+000dp+000/+200/-200.dat` and
  `fltst_strk+200dp+000.dat`/`-200dp+000.dat` (the stations on the strike=0 or
  dip=0 axes, i.e. nearest the well) diverge from their own just-frozen values
  by up to ~8% in column 5 (effective normal stress, MPa) by the run's end;
  the other 4 off-axis stations match to `0.0e+00`. Values on both sides stay
  physically plausible (~-24 to -30 MPa, near `par.init_norm = -25e6`), not
  clamped or wrong-looking — consistent with the kind of MPI-reduction-order
  -> RSF chaotic amplification this project has already documented elsewhere
  (`PATHWAY_FORWARD.md`'s bp1002caps finding: "binary version plus chaotic
  divergence, not the caps"), not an obviously new bug, but NOT independently
  confirmed to that specific mechanism here. This means the `onfault` category
  for THIS row is not a reliable regression gate at `rtol=2e-06` even against
  itself — a future CI/full-tier run may legitimately fail 4 of 9 station
  files with no code change at all. Loosening the tolerance or dropping those
  4 stations from the reference would change what "parity" means for this
  row (a methodology decision), so neither is done here; flagged for the
  project owner. The reference is committed as-is (unmodified from the frozen
  run) because it is still the correct record of what that run produced.
- Included files mirror `reference/test.bp8.qdc.gs.10`'s categories
  (`fault.*.nc`, `fltst_strk*`, `srfst_strk*`, `global.dat`, `*_strike.dat`/
  `*_depth.dat`, `cplot_EQquasi.txt`, `runInfo.json`); `disp.*.nc` and the
  `eqquasi.mesh.*.nc` mesh files were dropped — no comparison category reads
  them (`testsys/e2e/cases.py`'s `CATEGORIES`), so keeping them would be gold
  nothing asserts on (rule 8a).
