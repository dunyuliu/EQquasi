# liu2020.kink gold — first reference, regression lock only

Cycle 0 of `compsets/liu2020.qdc.kink.600` (paper parameters of Liu, Duan &
Luo 2020 at dx = 600 m, Λ₀/dx = 0.94), produced by `eqquasi-1.14.0` on knox,
2026-08-15, from `work/liu2020.kink.10cyc`.

**Caveats, all deliberate:**
- First reference for this configuration — a regression lock, **not** a
  validation against the paper (that needs dx = 300 m / 64-bit MUMPS, and the
  paper's ruptures are EQdyna's; this is quasi-dynamic).
- Cycle 0 ends at the `nstep = 10000` cap mid-decay (final V 2.16e-3), not on
  the exit criterion. Deterministically reproducible, but a truncated cycle.
- No e2e row: one cycle costs ~5 h. The utilities contract
  (`test_utilities_run_on_every_reference`) is the reader that keeps these
  files live (rule 8a).

Sequence context: cycles 0–2 of the source run are coherent (0 oscillation
swings) — the dx = 1200 oscillation is absent at 600 m.
