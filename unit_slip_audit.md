# Row 3 Audit: EQquasi Traction Formula vs FEniCSx Green's Function

**Date**: 2026-10-12  
**Node tested**: Interior patch at (x, z) = (-300.0, -300.0) m, mesh h=100 m, p=1 (linear elements)  
**Setup**: 81 fault nodes, lumped area per node = 1.0e4 m²

## FEniCSx Green's Function Result

For 1 m unit strike-slip (s_x = 1 m, s_z = 0), the traction is:
- tau_x = -3.227414e+08 Pa  (strike direction)
- tau_z = -2.99e+02 Pa  (dip direction, essentially zero as expected)

This gives an implicit shear stiffness of 3.227e+08 Pa for the fault, or 0.010× the material shear modulus μ = 3.204e+10 Pa. This is physically reasonable for a fault element at this resolution.

## EQquasi Traction Formula Analysis

**Formula from src/faulting.f90 lines 125-126:**

```fortran
tstk = (mslav*mmast*(fvd(5,2,2)-fvd(5,1,2))/dt + mslav*fvd(5,2,1) &
    - mmast*fvd(5,1,1)) / mtotl + fsfault
```

where `mtotl = (mslav + mmast) * arn(i,ift)`.

**In quasi-static limit (accel → 0):**

```
tstk_static = (mslav*fvd(5,2,1) - mmast*fvd(5,1,1)) / mtotl + fsfault
```

## Critical Issue: Cross-Term Structure

The formula uses:
```
mslav * F_master - mmast * F_slave
```

This **crosses master/slave labels with mass weightings**, which is structurally suspicious.

A standard mass-weighted average would be:
```
mslav * F_slave + mmast * F_master   (natural weighting)
```

**Analysis**: If forces satisfy equilibrium (F_master ≈ -F_slave), then:
- **Current formula**: mslav*(-F_slave) - mmast*F_slave = -(mslav + mmast)*F_slave
- **Correct formula**: mslav*F_slave + mmast*(-F_slave) = (mslav - mmast)*F_slave

For balanced masses (mslav ≈ mmast), the ratio of magnitudes is:
```
|EQquasi| / |Correct| ≈ (mslav + mmast) / (mslav - mmast) ≈ ∞  (singularity)
```

More realistically, if masses are nearly equal, the current formula can amplify or sign-flip the traction incorrectly.

## Discrepancy Consistency with Board Finding

The board (PATHWAY_FORWARD.md row 3) reports:
- **FEniCSx at 25 m**: converges to Kim within ~2% (linear hexes verified)
- **EQquasi at 50 m / 10 m**: -18% / -13% vs Kim
- **Pressure, domain, BC**: ruled out as causes

The -13% / -18% magnitude error is consistent with a systematic **weighting or scaling error** in how EQquasi combines master/slave forces to compute traction — not a sign flip (which would be ±100%), but a moderate systematic bias.

## Deliverable

**Real numeric measurement (signed):**
- Node 10, interior fault
- FEniCSx tau_x for 1 m unit slip: **-3.227e+08 Pa** (negative = compressive, physically expected)
- Ratio vs reference (Kim): **not measured** (would require EQquasi run on same case)

**Structural finding:** The cross-weighting in faulting.f90:125-126 (mslav×F_master - mmast×F_slave) is a red flag for a master/slave index swap or incorrect force combination logic.

**Next step for owner:** Audit whether the formula should instead use (mmast×F_master + mslav×F_slave) or (mslav×F_slave + mmast×F_master), and verify against first-principles elasticity or a reference like Freund (1979) / Melosh & Raefsky (1981).
