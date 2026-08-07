# BP8-QD-GS: frozen record of the CRESCENT DET comparison against `taehoKim_ref`

Read from the DET viewer on **2026-08-07**, stations `(0, 0)` and `(-200, 0)`.

The screenshots themselves cannot be archived -- image attachments are not
written to disk -- so the reference curves are **transcribed by eye** and marked
`~`. Our own curves are exact and archived under `data/` (subsampled to 400 rows
so they can live in git).

| legend on the platform | what it is |
|---|---|
| `dliu_eqquasi14650m` | ours, v1.4.6, dx = 50 m, aging law, `tau^0` derived = 12.9277 MPa (reading B) |
| `taehoKim_ref` | the reference, HBI (hierarchical boundary integral) |

`data/` also carries our **reading A** (`tau^0` = 14.6 MPa, dx = 10 m) curves,
which are the ones actually comparable to the reference -- see below.

## FROZEN: what `taehoKim_ref` is, and how the three readings compare

Read directly off the DET plots, station (0,0) and (-200,0), 2026-08-07:

| reference quantity | value | read from |
|---|---|---|
| `tau^0` at t = 0 | **14.6 MPa** | `shear_stress_2` panel, both stations |
| `V` at t = 0 | **1e-12 m/s** (log10 = -12) | `slip_rate_2` panel, both stations |
| `slip_2` at 30 d, (0,0) | **~38 mm** | `slip_2` panel |
| `slip_2` at 30 d, (-200,0) | **~21 mm** | `slip_2` panel |
| edge/centre ratio | **0.553** | derived |
| `state` | **flat ~2.8** for 30 d | `state` panel, both stations |
| pore pressure peak, (0,0) | ~13.4 MPa | `pore_pressure` panel |
| late-time uniform pressure | **~1.7 MPa** | `pore_pressure` panel, (-200,0) |

`tau^0` = 14.6 MPa **and** `V(0)` = `V_init` together force
`theta_0` = 4.0188e11 s, i.e. Table 1 + eq. (27), with eq. (29) given up. That is
**reading C**, and it is what the compsets now ship.

| reading | `tau^0` | `V`(step 1) | slip(0,0) | slip(-200,0) | ratio |
|---|---|---|---|---|---|
| A  Table 1 + eq. (29) | 14.6000 | 6.54e-11 | 42.74 mm | 27.65 mm | 0.647 |
| **C  Table 1 + eq. (27)** | **14.6000** | **1.00e-12** | **37.30 mm** | **19.56 mm** | **0.524** |
| B  eq. (27) + eq. (29) | 12.9277 | 1.00e-12 | 23.97 mm | 5.29 mm | 0.221 |
| **`taehoKim_ref`** | **14.6000** | **1.0e-12** | **~38 mm** | **~21 mm** | **0.553** |

Ours at dx = 50 m, old fluid scheme, ~22 d; the reference at 30 d. C is within
2 % at the centre and 7 % at 200 m.

## Slip at 30 d, the headline comparison

| station | `taehoKim_ref` | ours reading A, dx = 10 | ours reading B, dx = 50 |
|---|---|---|---|
| `(0, 0)`    | ~38 mm | 45.04 mm (+18 %) | 23.97 mm (-37 %) |
| `(-200, 0)` | ~21 mm | 31.95 mm (+52 %) | 5.35 mm (-75 %) |

**Reading A is much the closer of the two**, at both stations. Reading B, which
the compsets currently ship, is the worse match -- see "What went wrong" below.

Our slip profile is also **broader** than the reference. Edge/centre ratio,
`slip(-200,0) / slip(0,0)`:

    taehoKim_ref     0.553
    ours reading A   0.709

We put relatively more slip out at 200 m than the reference does. That shape
difference, not the amplitude, is the substantive discrepancy -- and see
"Elastic domain truncation" in README: it shrinks as the computational box
grows, so it is at least partly the truncation that section 6 warns about.

## Station (0, 0)

| quantity | `taehoKim_ref` | ours (reading B, dx = 50) |
|---|---|---|
| `slip_2` | ~0.038 m, plateau by ~0.4e6 s | 0.02397 m |
| `slip_3` | flat 0 | ~-2e-19 m (float noise) |
| `slip_rate_2` | peak ~-6.5, decays to ~-15 | peak -6.748, decays to ~-16.5 |
| `slip_rate_3` | flat 0 | ~-25 to -30 (float noise) |
| **`shear_stress_2`** | **starts 14.6 MPa**, min ~7.0, recovers to ~8.3 | **starts 12.9277 MPa**, min ~6.9, flat ~7.15 |
| `shear_stress_3` | flat 0 | +-1e-15 MPa (float noise) |
| `pore_pressure` | peak ~13.4 MPa at `t_off` | 13.625 MPa |
| `darcy_vel_2` | ~-4.4e-7 m/s, roughly constant during injection | **0** |
| `darcy_vel_3` | ~-4.4e-7 m/s | **0** |
| `state` | **flat ~2.8** for 30 days | 8.699 -> 3.5 -> 6.35 |

## Station (-200, 0)

| quantity | `taehoKim_ref` | ours (reading B, dx = 50) |
|---|---|---|
| `slip_2` | ~0.021 m, still rising at 30 d | 0.00535 m, plateaued |
| `slip_rate_2` | peak ~-7.4 at ~0.35e6 s, decays to ~-8.9 | peak ~-7.8, decays to ~-10 |
| `shear_stress_2` | rises 14.6 -> **peak ~17.8** at ~0.3e6 s, then falls to ~13 | rises 12.93 -> peak ~15, falls to ~12.8 |
| `pore_pressure` | peak ~2.75 MPa, settles ~1.7 | 2.856 MPa, settles ~1.5 |
| `darcy_vel_2` | ~-1.37e-6 m/s at the dip | -1.471e-6 m/s (good agreement) |
| `darcy_vel_3` | ~-3.1e-8 m/s | **0** |
| `state` | flat ~2.75 | 8.699 -> 5.0 -> 6.2 |

Note the shear stress *rises* before falling here: the slipping patch loads the
surrounding fault before the front arrives. Both codes show it.

## The three things this settles

**1. The reference uses Table 1's `tau_init` = 14.6 MPa.** Its `shear_stress_2`
starts there at both stations. So the benchmark's own author resolves the
over-determined initial condition by keeping Table 1 and letting `V(0)` be
6.54e-11 m/s, 65x `V_init`, rather than honouring eq. (27). Reading A is the
community-comparable choice. See README, "The initial condition is
over-determined".

**2. The reference's state variable is flat at ~2.8 for 30 days**, at both
stations. The aging law cannot do that -- once slip decays `d(theta)/dt -> 1`
and `theta` grows to ~1e6 s, as ours does. With the slip law removed from BP8 on
2026-08-06 this is **unexplained**.

**3. The reference's Darcy velocity is inconsistent with its own pressure
field**, and the earlier "half-cell collocation offset" explanation does **not**
survive measurement:

  - At `(-200, 0)`, `q_2` agrees well (-1.37e-6 vs our -1.471e-6), but `q_3` is
    -3.1e-8 where symmetry on the `x3 = 0` axis requires exactly 0. The ratio
    `q_3/q_2 = 0.023` implies an offset of ~4.5 m at r = 200 m, i.e. about half
    a cell. Consistent with a collocation convention.
  - At `(0, 0)`, both components are ~-4.4e-7, so `|q| = 6.22e-7 m/s`. Our field
    -- verified against the analytic solution at all nine stations, exact at
    r = 200 m and 0.4 % at r = 283 m -- reaches that magnitude only at
    **r ~ 278 m**. That is 28 cells, not half of one. Meanwhile its pore
    pressure at the same station peaks at ~13.4 MPa, which is unambiguously the
    injection point.

  So within a single file the pressure says "centre" and the Darcy velocity says
  "r ~ 280 m". A half-cell offset explains `(-200, 0)` but cannot explain
  `(0, 0)`. Our value at `(0, 0)` is 0 and is correct on principle: the pressure
  maximum sits there, a maximum has zero gradient, and `q = -(k/eta) grad p`.

## A cosmetic defect of ours

`slip_3`, `slip_rate_3` and `shear_stress_3` are identically zero for this
problem, but we emit float noise (~1e-19 m, ~1e-15 to 1e-16 MPa) that the
platform renders as a dense band filling the panel, where the reference plots a
clean line. Not wrong, but it reads badly. Worth flushing to zero below a
threshold before the next upload.

## What went wrong, and the lesson

An earlier comparison was made at `(-200, 0)`, where the reference reads ~21 mm.
That number was then carried forward as though it were the **centre**-station
value and set against our centre-station 45 mm, manufacturing an apparent 1.8x
disagreement that did not exist. The real figure at matched `tau^0` is +18 %.

The invented gap drove an investigation into domain truncation, radiation
damping, the state evolution law and a rewrite of the initial condition. None
explained a 1.8x that was never there, and it led to switching the shipped
default from reading A to reading B, which made our submission *less* comparable
at both stations, not more.

The irony is that domain truncation was dismissed too early and *is* a genuine
contributor -- just to the shape, not the amplitude. It was ruled out on the
centre-station slip and the peak slip rate, which are precisely the quantities
least sensitive to it. Doubling the box moves the centre by 2.6 % and the edge by
10.8 %. Ruling a cause out requires measuring the quantity it would actually
affect.

Two things would have prevented it. Record the station alongside any number read
off someone else's plot -- a value without its coordinates is not a measurement.
And when a discrepancy appears, check it at more than one station before
theorising about its cause; the second station here would have exposed the error
immediately.

## Files

    README.md   this record
    data/readingA-dx10_fltst_strk+000dp+000.csv
    data/readingA-dx10_fltst_strk-200dp+000.csv
    data/readingB-dx50_fltst_strk+000dp+000.csv
    data/readingB-dx50_fltst_strk-200dp+000.csv

400 rows each, subsampled from the full time series, 11 columns in the section
4.1 field order. Reading A is `work/bp8.prod10`; reading B is `work/bp8.test50`.
