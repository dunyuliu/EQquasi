# EQquasi v1.24.14 — 2026-11-25

Tag target: `e11190e1bb120a3d76b6ac77a466d13241753079` ("Release prep v1.24.14: dip comment versions, gate-script test can fail (#283)"). `EQQUASI_VERSION = '1.24.14'` (`src/globalvar.f90`). Previous tag: v1.24.11.

**Versions covered (rule 21c).** 1.24.12 (#266, `c392cc4`) and 1.24.14 (#278, `97405bf`). There is no 1.24.13 on master: that label belonged to the station-velocity divisor fix (draft PR #277, not merged), so no binary from master ever reported it. Neither 1.24.12 nor 1.24.14 was tagged separately; this release covers both.

## Behaviour a user can see

- **An unconverged Newton solve now aborts (1.24.12, #266).** The nodal Newton solve in `src/faulting.f90` that reaches its cap (`ivmax = 100`) prints `NEWTON-RAPHSON DID NOT CONVERGE: iteration cap reached`, with the time step, predictor/corrector pass, fault index, node (x, z), `v_trial` and residual, then calls `MPI_ABORT` with code 509 (exit status 253 under OpenMPI). Before, it printed a warning and continued; a run that used to finish with an unconverged node now stops. Rule 2b codifies this.
- **Planar faults now dip at `par.dip` (#275, closes #271; script only, no version bump).** `scripts/generateFaultInterface` offset each row by `dx*cos(dip)`, which built a plane dipping atan(1/cos(dip)) (50 degrees became 57.3). It now offsets by `dx/tan(dip)`. Dip 90 is unchanged (both forms write 0.000000). **Any case with `par.dip` other than 90 meshes differently from before and is not comparable with a run from an earlier version.** No committed reference changed.
- **Solver 4 far-field proportionality check (1.24.14, #278, closes #273).** In `ddS2tApply` (`src/solveTimeLoopPETScDD.f90:1582-1583`) the check on `constmp/w` over the constrained far-field dofs was `amax - amin > 1.0d-12*max(abs(amin),abs(amax))` (and the same for the zero-spread test). The tolerance is now `max(1.0d-12, 4*epsilon(1.0d0)*n)` where `n` (`s2t_far_nacc`) counts time steps in the corrector pass and is reset with the far-field running sum. The floor is reached up to about 1126 steps; the tolerance is 2.1e-11 at 24000 steps. Reason: the sum is a running sum, so its round-off spread grows with the step count (arithmetic: 1.339e-12 at 24000 steps against the old 1e-12). The restart message for a `fault.r.nc` lacking `dip_node_spacing` now names both causes (a coupled restart from EQdyna, a file written by an older EQquasi) and a remedy.
- Tests and tooling: e2e reference tests are selected at collection (#261); `tests/release_gate.sh` is the checked-in full-tier gate (#263, #267); rules 2b and others amended (#264, #279).

## #270 status: open, not reproduced by a run

Issue 270 (solver 4 aborts long runs at 1e-12) is **not closed**. A coarse BP5 solver-4 case with a non-zero far-field rate ran 30000 steps on the v1.24.12 build (302 s) and the v1.24.14 build (289 s), both exit 0, no `EQ_S2T: far-field` abort (`runs/20261125_c44_s2t_270/`; an earlier run at `nstep 30000` stopped at step 4966 because the post-rupture exit criterion fired, `src/onFaultHelpers.f90:29-50`). The tolerance change rests on the round-off arithmetic and a source-text test only.

## Gate evidence

- Full tier (`tests/release_gate.sh`): 695d92c (#278), `runs/20261124_c43_gate_695d92c/`, 145 passed, 11 skipped; 7bc744f (#275), `runs/20261124_c43_gate_7bc744f/`, PASS; 6f48623 (#266), `runs/20261124_c42_gate_6f48623/`, 145 passed, 11 skipped. All 11 skips are one reason: no reference at `data/test.bp1002.qdc.zone.2500/cycle0` (`tests/e2e/test_benchmarks.py:121`; open owner question).
- **The gate was not rerun on the tag SHA.** Decision, by rule 21 (which requires CI on the exact SHA, not a gate on it): `git diff 695d92c e11190e -- src tests/e2e` is empty; the only code-side differences from the gated 695d92c tree are the `generateFaultInterface` change (gated at 7bc744f), its unit test, a comment in `defaultParameters.py` and a tightened contract test (#283). On the tag SHA: CI and the stranger clone, next section.
- No performance claim is made by this release.

## Tag-SHA evidence (CI and stranger clone)

- CI on `e11190e`: run 45578623599, conclusion success (`build` success, `deploy` skipped).
- Stranger clone (`runs/20261125_c44_stranger/`, bare clone under `env -i`, empty HOME, existing `eqquasi-petsc` conda env reused instead of `conda env create`): README clone and `make.scripts.executable.sh` rc 0; `bash install.eqquasi.sh -m conda-linux` rc 0 (`bin/eqquasi-1.24.14`); `pytest -q -m "not e2e" tests` 413 passed, 156 deselected, 0 skipped (98 s); `create.newcase mycase bp5.qdc.2000`, `./case.setup`, `bash run.sh` rc 0, 57 min (02:50 to 03:47 CST) on a host at load about 47. The conductor's own script died once at `source install.eqquasi.sh` on its `set -u`, not on the README (`stranger.part1.out`).

## Audit

Code audit of `v1.24.11..origin/master` (head b09594b), verbatim from the read-only code auditor, 2026-11-25:

> 0 blocker, 0 major, 6 minor, 3 advisory. Code read only; executed: a dip-90 float check and git queries.
>
> Minor
> - `scripts/generateFaultInterface:108,114-117`: only `dip % 180 == 0` is refused; a dip outside 0 < dip <= 90 is accepted. `np.sign(dip)/tan(dip)` is even in dip, so a negative dip builds the same plane as its positive value (the old `cos` form mirrored it). No compset or doc uses a negative dip. The guard runs after line 108 has computed 0/0 = NaN.
> - `tests/e2e/test_benchmarks.py:91-92`: a reference lacking a category removes that test at collection instead of skipping; if a category is deleted, its check disappears with no skip line, so the rule 3a skip-count check cannot catch it.
> - `tests/unit/test_s2t_far_tolerance.py:410-414`: the arithmetic test exercises only numpy, never the Fortran check; the only solver-side guard is a source-text match (417-424). #270 has no run yet that aborts on the old version and passes on the new.
> - `tests/contract/test_release_gate.py:136`: `"MACHINE=" in src` is satisfied by the usage comment at `tests/release_gate.sh:4`; the real guard is the `die` at line 17.
> - `scripts/defaultParameters.py:119-123`: the comment says the `dx/tan(dip)` offset applies "from 1.24.12 on", but c392cc4 and 2b13241 also report 1.24.12 and still use `cos`; 1.24.14 is the first version with the fix.
> - `src/netcdf_io.f90:780-781`: the remedy command opens the fixed name `fault.r.nc` in append mode; the message says to stamp "a copy", so run as printed it edits the original in place.
>
> Advisory: version labels collide across branches (1.24.14 on unmerged spike d78d822, 1.24.13 on c245aef and 3e8fb71); the Newton cap test (`tests/e2e/test_benchmarks.py:377-378`) is `e2e_fast`, so not in CI's `e2e_ci`; `src/faulting.f90:598-603` `newton_solver` is an empty stub nothing calls (predates the range).
>
> 1. Newton cap: aborts on every path. The only Newton loop in `src/` is `src/faulting.f90:269-316`; `newton_converged` is set only at 304-306, so any other loop end, a NaN residual included, reaches 321-331: message, flush, `MPI_ABORT` 509. The separate `v_trial > 20` m/s warning (333-337) still warns and continues; it predates the range.
> 2. Far-field check: old `if (amax - amin > 1.0d-12*max(abs(amin), abs(amax)) .or. zmax > 1.0d-12*scale)`; new `tol = max(1.0d-12, 4.0d0*epsilon(1.0d0)*dble(s2t_far_nacc))` in the same test (1582-1583). `s2t_far_nacc` rises by one per step in the corrector pass (line 783), zeroed with `cons` at `src/eqquasi.f90:290` and `src/onFaultHelpers.f90:236`. The tolerance is 1e-12 up to about 1126 steps, 2.1e-11 at 24000, 8.9e-10 at 1e6; a defect would have to put under about 1e-9 relative error into the far field to hide.
> 3. Dip: old `sign(dip)*(nz-1-j)*dx*cos(dip)`, new `.../tan(dip)` (line 108); `tan` is never zero for 0 < dip <= 90; dip 0 stops the run (114-117); dip 90 is bit-identical (6.123233995736766e-15 for dx = 100, written as 0.000000 by `%f`, line 142).
> 4. Tests: no skip added; one runtime skip remains (`PENDING_REFERENCE`, bp1002 zone); four former runtime skips became collection-time deselection plus asserts; no tolerance loosened, no comparison stubbed.
> 5. Restart message (#273) names both causes and a remedy (`src/netcdf_io.f90:776-783`); the e2e test asserts four key phrases.

## Audit findings and fixes

Fixed in #283 (`e11190e`): the `defaultParameters.py` dip comment now says 1.24.14 / 1.24.12 (and `docs/user/parameters.md` regenerated); the `MACHINE` contract assertion now needs a non-comment `export MACHINE` line and the unset-`MACHINE` refusal.

Open issues, not fixed here (judgment calls or `src/` changes): negative or out-of-range dip accepted by `generateFaultInterface`; reference-less categories vanish at collection instead of skipping (owner decision D1 chose collection-time selection); the arithmetic-only #270 test; the restart-message remedy edits the original file (needs a `src/netcdf_io.f90` change and its own version); #270 itself.

## Files added / removed

Added: `tests/release_gate.sh`, `tests/unit/test_fault_interface_dip.py`, `tests/unit/test_s2t_far_tolerance.py`, this file. Removed: none.
