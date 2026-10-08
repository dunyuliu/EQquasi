# bp5.kink gold — first reference, regression lock only

Cycle 0 of `compsets/bp5.qdc.kink.2000` — the BP5 dip90 test configuration
with **only the fault surface changed** to a 10° kink at x = 0 (user-designed
control) — produced by `eqquasi-1.14.0` on knox, 2026-08-15, from
`work/kink.dip90ctrl.dx2000`.

The result this locks: under BP5 friction the rupture **crosses the bend
freely** (nucleates ~−25 km, sweeps to +30 km, no contour bunching at x = 0),
and cycles 0–2 of the source run are all physical (peaks 1.168 / 1.384 /
0.083 m/s, 0 swings). Contrast pair to BP1002, where a 5 km step-over is an
absolute seismic barrier.

First reference: a regression lock, not validated against any independent
oracle. Cycle 0 ends on the exit criterion (4876 steps, final 3.84e-6). Runs
in the e2e full tier.
