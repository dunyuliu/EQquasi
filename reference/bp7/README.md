# BP7-QD

SEAS benchmark problem 7, quasi-dynamic. Description:
[SEAS_BP7_May2023rev.pdf](https://strike.scec.org/cvws/seas/download/SEAS_BP7_May2023rev.pdf)

## The nucleation is prescribed, off-centre, and lives in the Fortran

This is the thing to know before reading a BP7 result, and it cost time here.

The compset (`case_input/test.bp7.qdc/user_defined_params.py`) looks perfectly
symmetric: a centred velocity-weakening disc, `radii = sqrt(x^2 + z^2) <= rad`,
uniform `Dc`, uniform initial slip rate, no seeded patch. So a result that
nucleates at (-50, -40) m rather than the origin looks like the solver breaking
symmetry — which is what I concluded, wrongly, before reading the description.

**Section 3.1 of the benchmark prescribes exactly that asymmetry.** A shear
traction perturbation is imposed from `t = 0`, smooth in space and time,

```
dtau^0(x2, x3, t) = dtau0 * G1(r) * G2(t),     r = sqrt((x2-y2)^2 + (x3-y3)^2)

G1(r) = exp( r^2 / (r^2 - Rnuc^2) )   for r < Rnuc,   0 otherwise      (31)
G2(t) = exp( (t-T)^2 / (t(t-2T)) )    for 0 < t < T,  1 for t >= T     (32)
```

with Table 1 giving

| | | |
|---|---|---|
| `dtau0` | max amplitude | 1.75 MPa |
| `Rnuc` | radius | 150 m |
| `T` | duration | 1 s |
| `(y2, y3)` | **hypocenter** | **(-50 m, -50 m)** |

The hypocenter is deliberately off-centre. A BP7 rupture that nucleates at the
origin would be wrong.

### Where it is implemented

Not in the compset — in the solver, under `bp == 7`:

| | |
|---|---|
| `src/globalvar.f90:131-133` | `nucx = -50`, `nucz = -50`, `nucdtao0 = 1.75d6`, `nucr = 150.0d0`, `nuct` |
| `src/faulting.f90:556` | `nucleation1`, equations (31) and (32) |
| `src/faulting.f90:162` | applies it, `icstart == 1` only -- the *first* rupture |
| `src/solveTimeLoopMUMPS.f90:64` | holds `dtev1 = dt` while `time < nuct`, so the 1 s ramp is resolved |
| `src/library_output.f90:337` | BP7 output |

They are compile-time constants. That is fine while they match Table 1, but it
means a reader of `user_defined_params.py` cannot see them, which is how a
specified asymmetry gets mistaken for a defect. Same shape as `hang_time` in
`solveTimeLoopMUMPS.f90`, which sets how long a cycle runs and therefore the
length of every full-cycle reference here.

**These belong in the compset, and moving them is tracked as work to do.** Until
then, a benchmark-specific treatment that lives in the source gets recorded in
that benchmark's reference README -- this section is the BP7 instance. If you
add another, write it down in the same place.

## What the reference holds

| | |
|---|---|
| `cycle0-step101-fast/` | the compset's own 101 steps. Pre-nucleation and still exactly symmetric under x- and z-reflection, since the perturbation has barely ramped. Used by `e2e_fast`. |
| `cycle0/` | the first full cycle: 1358 steps, 1.17 d, 202 of 1681 fault nodes rupturing between 0.575 and 1.675 s. |

Both produced through `create.newcase` -> `case.setup` -> `run.sh`, at
`dx = 25 m` (the benchmark suggests 10 m; this is the test resolution).
`disp.*.nc` omitted -- no comparison reads them.

The rupture is fast: 1.1 s across the whole patch, against BP5's 117 s. That is
why `plotRuptureTime.py` falls back from its fixed 5 s contour interval here --
at 5 s there would be no contour line inside the data at all.

## Symmetry, as a diagnostic

Useful when something looks wrong:

```
step    x-flip asym    z-flip asym    max|V|
 1      0.000e+00      0.000e+00      1.0e-09    <- exactly symmetric
 1001   1.473e-05      1.867e-05      3.3e-05    <- perturbation has acted
 1358   4.353e-09      4.880e-09      1.7e-08
```

Exact symmetry at step 1 and asymmetry later is the *expected* signature, not a
bug: `G2(t)` ramps the perturbation in over 1 s, so the initial condition is
symmetric and the solution is not.

Like every reference here, this is a regression lock, not a validation.


## Provenance

Blessed at **v1.7.0** (`cycle0/runInfo.json`). Confirmed still reproducible by the full e2e
tier on 2026-08-14, at the v1.13.0 workflow-revamp commit set.

The version gap is not staleness. These numbers are what current code
produces, checked every full-tier run; what was missing was any record of
which version produced them, so a future divergence could not be dated.
