# BP8-QD-GS: what the CRESCENT DET comparison against `taehoKim_ref` shows

Frozen record of the platform comparison, read from the DET viewer on
**2026-08-07**. The screenshots themselves are not recoverable (image
attachments are not written to disk), so the numbers are transcribed here.
Values read off log/linear axes are approximate and marked `~`; values from our
own output files are exact.

**Both panels below are station `fltst_strk+000dp+000`, i.e. `(x2, x3) = (0, 0)`,
the injection point.** Getting this wrong is what derailed several days of work
-- see "The station mix-up" at the end.

Datasets compared:

| legend | what it is |
|--------|------------|
| `dliu_eqquasi14650m` | ours, v1.4.6, dx = 50 m, aging law, `tau^0` derived = 12.9277 MPa (reading B) |
| `taehoKim_ref`       | the reference, HBI (hierarchical boundary integral) |

## Panel-by-panel, station (0,0)

| quantity | `taehoKim_ref` | ours (reading B, dx = 50 m) |
|----------|----------------|------------------------------|
| `slip_2` | plateaus at **~0.038 m**, reached by ~0.4e6 s | 0.02397 m |
| `slip_3` | flat 0 | ~-2e-19 m (numerical noise) |
| `slip_rate_2` | rises to a peak ~-6.5, decays to ~-15 by 2.5e6 s | peak -6.748, decays to ~-16.5 |
| `slip_rate_3` | flat 0 | ~-25 to -30 (numerical noise) |
| **`shear_stress_2`** | **starts at 14.6 MPa**, minimum ~7.0, recovers to ~8.2 | **starts at 12.9277 MPa**, minimum ~6.9, stays ~7.0 |
| `shear_stress_3` | flat 0 | +-1e-15 MPa (numerical noise) |
| `pore_pressure` | peak ~13.5 MPa at ~0.36e6 s, decays to ~1.5 | 13.625 MPa, decays to ~1.5 |
| `darcy_vel_2` | dips to **~-4.5e-7 m/s** | **0** |
| `darcy_vel_3` | dips to **~-4.5e-7 m/s** | **0** |
| `state` | **flat at ~2.8** for the whole run | 8.699 -> 3.5 -> 6.35 |

## The four things this settles

**1. The reference uses Table 1's `tau_init` = 14.6 MPa.** Its `shear_stress_2`
starts there, ours at 12.9277 MPa. So the benchmark's own author resolves the
over-determined initial condition by keeping Table 1 and letting `V(0)` be
6.54e-11 m/s -- 65x `V_init` -- rather than honouring eq. (27)'s uniform
`V_init`. Reading A is the community-comparable choice; reading B is the one
consistent with section 3's prose. They differ by 1.8x in slip. See README,
"The initial condition is over-determined".

**2. There is no 1.8x discrepancy with the reference.** At matched `tau^0`
(reading A, dx = 10 m) we give 45.04 mm against the reference's ~38 mm, i.e.
**+18%**. That is the real open question and it is an ordinary-sized one.

**3. Our Darcy velocity at (0,0) is 0 and the reference's is not.** Ours is
correct: the Gaussian pressure field has its maximum at the injection point, a
maximum has zero gradient, and `q = -(k/eta) grad p` must vanish there. Verified
against the analytic solution at all nine stations -- exact at r = 200 m, 0.4% at
r = 283 m. The reference's `|q| ~ 6.4e-7` corresponds to the analytic value about
one cell from the origin, and its two components being roughly equal implies an
offset along the diagonal. Most likely a collocation convention (cell-centre vs
node), not a physics error -- its pressure field agrees with ours closely.

**4. The reference's state variable is flat at ~2.8 for 30 days.** The aging law
cannot do that: once slip decays `d(theta)/dt -> 1` and `theta` grows to ~1e6 s,
as ours does. Since the slip law was removed from BP8 on 2026-08-06, this is
**unexplained**. Candidates: the column is not `theta` in seconds; the reference
predates the scope change and was run with the slip law; or it is wrong. Asking
for the actual data files is the way to settle it.

## Cosmetic issue in our own output

`slip_3`, `slip_rate_3` and `shear_stress_3` are numerically zero for this
problem but we emit float noise (~1e-19 m, ~1e-15 MPa), which the platform plots
as a dense band filling the panel. The reference emits exact zeros and plots as a
clean line. Not wrong, but it reads badly next to the reference and is worth
flushing to zero below a threshold before the next upload.

## The station mix-up

An earlier comparison in this project was made at station `(-200, 0)`, where the
reference reads ~21 mm and ours ~32 mm. That 21 mm was then carried forward as
if it were the centre-station value and compared against our centre-station
45 mm, manufacturing an apparent 1.8x disagreement. It sent the investigation
after domain truncation, radiation damping, the evolution law, and a rewrite of
the initial condition -- none of which were the problem, and one of which
(switching the shipped default to reading B) made our submission *less*
comparable to the reference, not more.

Our own values at both stations, reading A, dx = 10 m (`work/bp8.prod10`):

    x =    0 m : 45.04 mm
    x = -200 m : 31.95 mm

The lesson is narrow and worth stating plainly: when a comparison number comes
off someone else's plot, record the station with it. A number without its
coordinates is not a measurement.
