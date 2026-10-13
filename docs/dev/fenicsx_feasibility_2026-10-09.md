# FEniCSx feasibility spike (board row 2, P2) — 2026-10-09

Destination: `docs/dev/fenicsx_feasibility_2026-10-09.md` (approved name). Prototype:
`runs/20261009_fenicsx_spike/` (untracked). dolfinx 0.10.0, petsc 3.25.0 + MUMPS, env
`/home/staff/dliu/anaconda_knox/envs/fenicsx` unmodified. Host theo4, serial, 1 thread.

## Verdict

| question | answer |
|---|---|
| Split-node fault slip with rate-and-state in dolfinx | **Yes, cleanly** (below) |
| Do P2 elements close the BP8 gap to Kim | **The gap is not an element-order problem.** An independent *linear*-hex FEniCSx model converges to Kim; EQquasi does not. |
| Step cost vs EQquasi | **~225-1600x cheaper per RSF step** (below), but not an apples comparison -- see caveats. |
| Port | **No-go**, on both the depth-gap motivation and the step-cost gap (the latter is algorithmic, not language/library; see Q3). |

## Prototype

`bp8fe.py`: BP8-QD-GS, aging law, quasi-dynamic (`eta = rho*cs/2`), `tau0 = 14.6 MPa`,
`sigma0 = 25 MPa`, `V_init = 1e-12`, Omega_f = ±400 m, the rest of the fault locked. Elastic box
as EQquasi: x, z ±1500 m traction-free, y ±2000 m clamped; uniform `h` within ±500 m (x, z) and
0–400 m (y), graded 1.2 beyond. Half-space by antisymmetry: `u_t = s/2` on y = 0, `sigma_yy = 0`.
Fault traction = reaction / lumped nodal area (the same as EQquasi's split-node nodal force / area).
The slip-to-traction map is assembled exactly as a dense matrix G from MUMPS solves, so element
order is the only variable. Pore pressure comes from the exact Neumann cosine series (120² modes,
zero flux on Omega_f). Time integration: RK45, rtol 1e-6.

## Q1 — split node

Full space, conforming mesh. The fault is interior facets, slip enters as a discontinuous (DG_p)
lifting of ±s/2 on the cells either side (Melosh–Raefsky), and the fault force is assembled over
the + side cells only. No `dolfinx_mpc` needed (and it is not in the env).

| check (`split p h`) | Q1 h=100 | Q2 h=100 | Q1 h=50 |
|---|---|---|---|
| full-space vs half-space fault force, max rel. diff | 5.5e-16 | 1.6e-15 | 1.7e-15 |
| fault equilibrium `max|R+ + R-|` (rel.) | 2.1e-16 | 8.4e-16 | 6.6e-16 |

The first attempt differed by 0.25 %. Cause: DG dof coordinates at x = 400 differ by ulps between
cells, and a tolerance-free slip cut-off gave the two sides different slip. Any port must locate
fault dofs topologically, not by coordinate tests. Operator check against the analytic circular
crack (a = 300 m), centre traction / stress drop: Q1 0.9998 (h = 25 m), Q2 0.9992 (h = 50 m).
Not exercised: non-planar faults (topological fault-dof lookup) and tetrahedra (P2 triangle
vertices have zero lumped area, so they need a consistent face-mass solve).

## Q2 — the depth gap

slip_2 at 30 d, % vs Kim (`taehoKim_ref`, `scratch/crescent/portal_bp8.json`):

| station | Kim mm | EQq 50 m | EQq 10 m | Q1 100 | Q1 50 | Q1 25 | Q2 100 | Q2 50 |
|---|---|---|---|---|---|---|---|---|
| (0,0) | 38.112 | −2.18 | −6.16 | +5.32 | +1.10 | +0.49 | −2.04 | +0.03 |
| (−200,0) | 21.597 | −5.16 | −5.91 | −14.38 | +0.76 | +2.06 | +1.89 | +2.28 |
| (0,+200) | 19.736 | **−18.00** | **−13.10** | −28.01 | −6.02 | +0.52 | −1.17 | +1.62 |
| (+200,0) | 22.084 | −7.24 | −7.98 | −16.27 | −1.46 | −0.18 | −0.35 | +0.03 |
| (0,−200) | 19.971 | −18.96 | −14.12 | −28.85 | −7.13 | −0.66 | −2.33 | +0.43 |
| (±200,±200) range | 16.23–16.59 | −19.4…−21.1 | −13.1…−15.0 | −92.7…−92.9 | −5.4…−7.5 | +0.2…+2.4 | −0.9…+1.3 | +1.4…+3.7 |

- unknowns, half model: Q1 100/50/25 = 20,631 / 91,260 / 392,925; Q2 100/50 = 151,875 / 693,693.
  EQquasi: 50 m 141,398 nodes (full box, v1.20.1); 10 m 1,231,276 nodes (v1.21.0).
- **Rule 11:** the two EQquasi runs are different versions (1.20.1 and 1.21.0) with different
  meshes (50 m: uniform x/z ±1500; 10 m: ±500 box graded 1.2). The board's −18 %/−13 % pair is not
  a same-binary sequence.
- Prototype convergence: the two finest (Q1 25, Q2 50) agree within 1.2 % at every station. Q1
  observed order at (0,+200) = 1.75 (100/50/25), Richardson limit +3.3 % vs Kim — indicative only,
  since h = 100 is pre-asymptotic (it loses the diagonal stations). RK45 rtol 1e-8
  changes Q1 50 by < 1 µm.
- Pressure (series vs Kim), centre: 7.13/7.10 MPa at 1 d, 13.06/12.99 at t_off, 2.59/2.59 at
  10 d. EQquasi 10 m vs Kim at t_off: 12.94/12.94. Pressure does not explain the gap.
- EQquasi integrates with 2×2×2 Gauss points (`nint=8`, `src/globalvar.f90`), as the prototype
  does, and its 50 m mesh is at least as fine as the prototype's. Yet it gives −18 % where
  prototype Q1 50 gives −6 %, and its centre station moves *away* from Kim from 50 → 10 m. The
  residual is EQquasi-specific (fault update, time loop or setup), not linear-hex discretization.
  Routed to lars-eriksson; not chased here.
- What P2 buys in general: at matched fault-node spacing (50 m), Q2 100 gives −1.2 % at (0,+200)
  against Q1 50's −6.0 %, and at 100 m Q1 loses the diagonal stations (−93 %) where Q2 is within
  1.3 %. Cost: memory — Q2 50 peaked at 19.3 GB, Q1 25 at 7.6 GB, for a similar result.
- Kim is not exact either: pressure at (0,+200) vs (+200,0) differs by 12 % at t_off
  (2.43 vs 2.71 MPa) in a radially symmetric diffusion problem; mirror stations differ by up to 2.2 %.

## Q3 — step cost (closed, 2026-10-12)

Run 2026-10-12 22:46-23:03, `bash pipeline.sh` (`green.log`, `rsf.log`), host load 0.15-1.19
(`uptime` before/during) -- **contended, not host-exclusive** (CLAUDE.md: timing comparisons need
the host to themselves); the fields below are the pipeline's own `wall_contended_s`, same label
the prior session used. No exclusive rerun was done: the gap is wide enough (2-3 orders of
magnitude) that ordinary background load cannot be the explanation.

| | nfault (Omega_f) | per-step (RSF stepper) |
|---|---|---|
| FEniCSx Q1 h=100 | 81 | 0.71 ms (1.3 s / 1840 steps) |
| FEniCSx Q1 h=50 | 289 | 1.25 ms (2.0 s / 1594 steps) |
| FEniCSx Q1 h=25 | 1089 | 4.98 ms (7.1 s / 1426 steps) |
| FEniCSx Q2 h=100 | 289 | 1.22 ms (2.0 s / 1642 steps) |
| FEniCSx Q2 h=50 | 1089 | 4.70 ms (7.8 s / 1660 steps) |
| EQquasi 50 m (8 ranks) | 3,721 | 1,123 ms (`scratch/bp8.gs.50m.xz1500.y2000`) |
| EQquasi 10 m wide (8 ranks) | 10,201 | 13,632 ms (`scratch/bp8.gs.10m.wide.v1.21.0`) |

Not an apples comparison, and the gap does not transfer to a port decision:

- **Different problem per step.** The FEniCSx RSF step is a dense matvec against a precomputed
  fault-only Green's matrix (81-1089 dof, planar fault, linear elasticity, quasi-static). EQquasi
  solves the full 3D volume FEM system every step (380k-3.5M equations) and supports nonlinear and
  dynamic physics the dense-G approach does not touch.
- **The precompute is excluded from the per-step number and does not scale well**: building +
  factoring G cost 1.8 s (Q1 100) up to 668 s (Q2 50, 19.3 GB peak RSS) -- it is the dominant cost
  at this problem size and would dominate further for a non-planar fault or a larger Omega_f.
- **Rank mismatch.** EQquasi's numbers are at 8 MPI ranks; FEniCSx ran single-threaded serial.
  Neither EQquasi single-rank nor a parallelized dense solve was measured.

The precomputed slip-to-traction matrix is a scaling idea for row 1 (already noted there and in
row 4), not evidence for a port: it only applies where the fault problem is quasi-static-linear
and planar, which is not EQquasi's general case.

## Provenance and caveat

Commands: `bash run.sh "split 1 100" "split 2 100" "split 1 50"`; `bash pipeline.sh`;
`python compare.py`, run 2026-10-09 12:55–13:21. The originating worktree
(`.claude/worktrees/agent-aed14af1bc8baa76b`) was removed at about 13:19 while the pipeline was
running. The scripts here are reconstructed verbatim from that session. The raw logs, G matrices
and `rsf_*.npz` were lost. The numbers above are transcribed from that session's console output
and can be reproduced in about 25 min serial (do not run during timing work).

Q3's numbers (`green.log`, `rsf.log`) are from a fresh `bash pipeline.sh` run 2026-10-12
22:46-23:03, `runs/20261009_fenicsx_spike/`, not host-exclusive (see Q3).
