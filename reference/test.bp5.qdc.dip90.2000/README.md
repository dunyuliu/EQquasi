# reference/bp5.dip90

First reference for `test.bp5.qdc.dip90` (BP5 friction, dip-90 planar control
with the surface-kink generator geometry shipped statically). Regression lock,
not a verified benchmark: the numbers are checked against nothing external.

- Frozen 2026-08-15 from `work/deorphan.dip90` (created by `create.newcase`
  from the compset, `case.setup`, `bash run.sh`, unmodified).
- Binary `eqquasi-1.16.0`, knox, 2 MPI ranks (`HPC_ncpu = 2`).
- cycle0: 101 steps (nstep cap, by design for the gate), peak Vmax 0.0698 —
  the run ends mid-nucleation; that is the state the gate locks.
- This compset was orphaned (no e2e row) from 2026-08-12 to 2026-08-15; the
  `generateFaultInterface` input/-vs-root layout bug stayed latent exactly
  because nothing ran it (PATHWAY_FORWARD, queue item 4).
