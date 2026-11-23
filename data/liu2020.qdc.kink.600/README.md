# liu2020.kink gold — regression lock, regenerated with v1.24.8

Cycle 0 of `compset/liu2020.qdc.kink.600` (paper parameters of Liu, Duan &
Luo 2020 at dx = 600 m, Λ₀/dx = 0.94). Regenerated 2026-11-20 on the owner's
decision (board row 1), replacing the `eqquasi-1.14.0` set of 2026-08-15.

## Provenance (from the run's own `runInfo.json` and log)

| Item | Value |
|---|---|
| Binary | `eqquasi-1.24.8` (banner of `log/cycle0.log`) |
| Commit that built the binary | not recorded in the run; the version string was introduced by the NetCDF 8-byte-reals change (#203), nothing stronger is claimed |
| Host | cotopaxi, AMD EPYC 7532, 64 logical cores |
| MPI ranks / threads | 3 ranks, 1 OpenMP thread per rank (`par.HPC_ncpu = 3`; the compset's 24 is the only difference in `user_defined_params.py`) |
| Solver | `solver` 1 (native), `friclaw` 3 |
| Mesh | dx = 600 m, 210000 elements, 226644 nodes, 5151 fault nodes, 618120 equations |
| Steps | 10000 (`exit_reason: nstep_reached`) |
| Simulated time | 5.871558e8 s |
| Time loop wall | 22285.9 s (6.19 h), 2.2286 s/step; factorization 29.6 s |
| Final max slip rate | 2.161455e-3 m/s |

The run is a clean start at step 0 (the first attempt of the same case was
killed at step 2821 and is not used).

## What is and is not claimed

- A regression lock for this configuration, **not** a validation against the
  paper (that needs dx = 300 m and the paper's ruptures are EQdyna's; this is
  quasi-dynamic).
- Cycle 0 ends at the `nstep = 10000` cap mid-decay, not on the exit
  criterion: a truncated cycle.
- Measured against the superseded v1.14.0 set (3 ranks, same case): step
  count, simulated time and final max slip rate agree at the printed
  precision of `runInfo.json`; the largest difference in `cplot_EQquasi.txt`
  is 1.5e-7 and in `global.dat` 3.9e-6, each as a fraction of its column's
  maximum. That is one comparison of two versions, not a tolerance.
- **Not claimed:** bit-for-bit or any tolerance of reproduction of this set on
  another rank count, host or build. It was produced once; no repeat run was
  made, so run-to-run and rank-count spread are unknown. The earlier text
  calling it "deterministically reproducible" was unsupported and is withdrawn.
- No e2e row: one cycle costs ~6 h. The utilities contract
  (`test_utilities_run_on_every_reference`) is the reader that keeps these
  files live (rule 8a).

`disp.10000.nc` is the last-step displacement (6.3 MB, under the rule 22
10 MB per-file cap). `model.txt` is the run's `input/model.txt`.

Sequence context (from the superseded run): cycles 0–2 were coherent (0
oscillation swings) — the dx = 1200 oscillation is absent at 600 m.
