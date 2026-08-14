# BP5-QD

SEAS benchmark problem 5, quasi-dynamic. Vertical strike-slip fault, rate- and
state-dependent friction with the aging law.

## The nucleation patch starts deliberately out of steady state

This is the thing to know before reading the compset, and it looks like a bug
until you see why it is not.

`compsets/bp5.qdc.2000/user_defined_params.py` seeds a patch — low `Dc`
(`par.minDc`) and a high initial slip rate of 0.03 m/s — and then sets three
quantities per node:

```python
par.on_fault_vars[0,iz,ix,46] = 0.03                       # V, in the patch
par.on_fault_vars[0,iz,ix,20] = Dc / par.creep_slip_rate   # theta, EVERYWHERE
par.on_fault_vars[0,iz,ix,8]  = shear_steady_state(..., par.on_fault_vars[0,iz,ix,46], ...)
```

The shear stress uses the node's **own** slip rate. The state variable does
not — it uses the creep rate, `1e-9` m/s, for every node including the patch.
Read side by side that looks like an oversight in one of the two lines.

**It is deliberate, and it is the mechanism.** At steady state
`theta_ss = Dc / V`. A node seeded at 0.03 m/s but given
`theta = Dc / 1e-9` is ~3e7 away from steady state for its own slip rate, and
that disequilibrium is what makes the patch accelerate. Setting
`theta = Dc / 0.03` would put the patch *at* steady state, and a patch at
steady state does not nucleate — it just keeps sliding at 0.03 m/s.

The shear line is different because shear must balance the friction the node
actually has at the rate it actually has, or the Newton solve resolves the
inconsistency in favour of the stress and drops V straight back to the creep
rate. That failure is recorded in the commit "Seed the step-over with shear
its own slip rate can sustain": the step-over compset passed `creep_slip_rate`
to `shear_steady_state`, the patch started 1.86 MPa short of the stress its
own rate needs, V collapsed at step 1, and `dtev = xi*Dc/V` jumped `dt` from
0.065 s to ~2e6 s.

So the asymmetry is correct in both directions:

| quantity | uses | why |
|---|---|---|
| shear `[..., 8]` | the node's own `V` | must balance friction at the actual rate, or V collapses |
| state `[..., 20]` | the creep rate | must be **out** of steady state, or the patch never accelerates |

`test.bp5.qdc` and `bp1002.qdc.2500` follow the same convention.

## References here

| | |
|---|---|
| `cycle0/` | the full first cycle, `full` tier |
| `cycle0-step101-fast/` | the same case stopped at step 101, every-push CI |

`summary.json` in each records what that run is expected to reproduce.

**Caveat.** These are regression locks, not validations: they detect
unintended change, and do not establish that BP5 is reproduced correctly
against the benchmark description or against other groups' submissions.


## Provenance

Blessed at **v1.7.2** (`cycle0/runInfo.json`). Confirmed still reproducible by the full e2e
tier on 2026-08-14, at the v1.13.0 workflow-revamp commit set.

The version gap is not staleness. These numbers are what current code
produces, checked every full-tier run; what was missing was any record of
which version produced them, so a future divergence could not be dated.
