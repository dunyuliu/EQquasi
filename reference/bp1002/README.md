# BP1002 — two-fault step-over

BP5 parameters, geometry changed only: two right-lateral segments offset
across strike, so this is the gate's **only `ntotft > 1` case**.

| | |
|---|---|
| Segment A | `x ∈ [-42.5, -2.5]` km at `y = 0` km |
| Segment B | `x ∈ [-7.5, 32.5]` km at `y = -5` km |
| Both | `z ∈ [-20, 0]` km, `dx = 2500` m |

The segments **overlap** along strike by 5 km: B's near end sits alongside A,
not beyond its tip. Cycle 0 nucleates on A and, after a delay, jumps — B
nucleates around step 3000 and peaks near step 6000 as A arrests, and both
faults end up slipping their full length (9.49 m and 9.41 m, 153/153 nodes
each, global peak 1.19 m/s).

An earlier version of this reference reported the opposite — B locked in a
stress shadow, never rupturing. That was an artifact of `meshgen.f90` /
`mesh4num.f90` tagging **every** node on a fault's y-plane as a fault
boundary, on the y coordinate alone. Nodes outside the fault's own x-z extent
got no equation number and were never split, so 1222 of the 1528 nodes on each
plane sat clamped at exactly zero displacement while the material around them
moved metres. The rigid sheet manufactured +25 MPa of normal-stress change at
the fault edges against a ~7 MPa real signal, and killed the run in later
cycles with STOP 508. Single-fault benchmarks never exposed it, because their
fault spans the whole model and the plane is fault everywhere.

## Why this case is in the gate

Every other benchmark here has `ntotft = 1`. Nothing in the gate exercised
per-fault routing of on-fault input, and three multi-fault bugs in this
project's history were found by hand for exactly that reason
(`PROJECT_RULES.md`, rule 7). This case is the fix: `case.setup` writes a
`(ntotft, nfz, nfx, nvar)` array, `netcdf_read_on_fault` reads it back with
an explicit fault dimension, and the reference locks the result.

The seed is on **fault 0 only** — fault 1 starts flat at the creep rate and
has to be brought to failure by the elastic interaction alone. A routing bug
that swapped or aliased the faults changes which fault carries the seed, which
is the point.

## Why dx = 2500 m and not BP5's 2000 m

The fault-normal offset between the segments is 5000 m, and the mesher
(`src/func_lib.f90`, `build_yline_belt`) requires that offset to be an integer
multiple of `dy = dx`. It refuses to mesh otherwise rather than silently
producing a fault with zero nodes: `5000/2000 = 2.5` is refused,
`5000/2500 = 2` meshes.

## The initial condition this case got wrong once

The compset builds each node's initial shear as the steady state for **that
node's own** initial slip rate, `on_fault_vars[..., 46]`. It once passed
`creep_slip_rate` for every node including the nucleation patch, whose initial
slip rate is 0.03 m/s. The patch then started 1.86 MPa short of the stress its
own slip rate needs; the Newton solve sided with the stress, `V` collapsed to
the creep rate at step 1, and `dtev = xi*Dc/V` jumped `dt` from 0.065 s to
~2e6 s, flattening both faults. See the commit "Seed the step-over with shear
its own slip rate can sustain".

## Tier

`cycle0/` — the complete event, in the **full** tier only. It is not in the
every-push fast tier, which leaves the gate with no `ntotft > 1` case; that
gap is deliberate and recorded in `tests/e2e/cases.py`, not an oversight.
The event is 7417 steps and ~700 s on 3 ranks.

**Caveat.** First reference for this benchmark — a regression lock, not a
validation. It has not been checked against any prior oracle.
