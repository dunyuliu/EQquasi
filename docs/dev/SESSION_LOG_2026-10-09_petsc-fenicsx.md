# Session log 2026-10-09 — row 1 (PETSc scaling) + row 2 (FEniCSx spike)

Conductor: wei-lin (autopilot). Budget read as milestone-based (rule 19/Phase
3a): run until row 1 and row 2 each reach their board "Done when" criterion or
a genuine blocker. Merge grant: rule 21 (conductor may merge own PRs to master
on the rule-3 gate, cut patch/minor tags unattended; no major/force-tag/publish).

## Orientation
- Master 605eacf, clean, no open PRs, one worktree (main). theo4: 48 cores,
  idle at session start (load 0.26).
- `bin/eqquasi-1.21.2` did not exist in either checkout despite being named as
  "today's binary" — built fresh (rule 9) via `install.eqquasi.sh -m
  conda-linux`, matches `src/globalvar.f90`'s `EQQUASI_VERSION` at HEAD.

## Row 1 — PETSc scaling, step (a): per-step profiling

10-step profiling runs (not gold/parity runs), native MUMPS (`par.solver=1`),
`OMP_NUM_THREADS=1`, binary eqquasi-1.21.2, HEAD 605eacf.

**bp1002.qdc.2500** (127,221 eq, 2-fault, clean throughout):
8r 0.436 s/step, 16r 0.409, 32r 0.394 s/step — flat; rank-0-only
faulting/pressure work dominates at this size (board's known gap, item (d)).

**BP8 "10 m wide" box** — finding: the committed `compset/bp8.qdc.gs.10`
(`fxmin/fxmax/fzmin/fzmax = -500/500`, 694k nodes) is NOT the box the board's
done-when criterion refers to. The actual "10 m BP8 wide" (1.23M nodes,
`-1500/1500`, `enlarging_ratio_xz=1.2`) only exists as a hand-edited
`scratch/bp8.gs.10m.wide.v1.21.0/user_defined_params.py`, never registered as
its own compset (rule 7). Recommendation left on the board row, not acted on
(rule 1 — new-compset naming needs agreement): register it as e.g.
`bp8.qdc.gs.10.wide`.

First sweep on the correct wide box (12:55–13:17) — **discarded, contaminated**:
dunyu-liu's row-2 Q2 convergence runs (up to 19 GB/11 min, serial but heavy)
ran on the same host in that exact window despite being briefed to hold only
*comparative timing* against EQquasi. Lesson: "no timing work while row 1 is
timing" must mean no heavy compute of any kind, not just a head-to-head
comparison — logged back to dunyu-liu, and the host's 1-min load average must
be back near baseline (~0.2-0.3) before trusting a rerun as clean, not just
`ps` showing no matching process (5/15-min load averages lag ~10 min behind
a finished heavy job).

Clean rerun launched 13:23 (host settled, load 0.97 1-min, no competing
procs): results pending — see board row 1 for final numbers once all three
rank counts complete.

## Row 2 — FEniCSx feasibility spike (dunyu-liu)

Dispatched in worktree `.claude/worktrees/agent-aed14af1bc8baa76b` /
`worktree-agent-aed14af1bc8baa76b`.

**Verdict so far (Q1, Q2 answered; Q3 pending)**:
- Q1 (split-node RSF in dolfinx): clean, no extra package needed, verified
  against an independent half-space formulation to 1e-15.
- Q2 (does P2 close the BP8-vs-Kim depth gap): **no** — the gap is
  EQquasi-specific. An independent FEniCSx linear-hex model on the same
  problem converges to Kim (+0.5% at 25 m); EQquasi's own runs sit at -18%
  (50 m, v1.20.1) / -13% (10 m, v1.21.0) and the error at the centre station
  gets *worse* with refinement (-2.2% → -6.2%), consistent with a real bug,
  not a discretization limit. **New finding, routed to lars-eriksson, not
  row 2's to fix**: the two reference runs compared are different EQquasi
  versions (1.20.1 vs 1.21.0) per rule 11 — flagged in the doc.
- Q3 (step cost vs EQquasi): not run yet; held for host-exclusive time after
  row 1's timing is done.
- Lean verdict: no-go on a port for the depth-gap motivation, pending Q3 and
  the bug hunt.

**Incident — worktree auto-reaped mid-task, lost raw data.** The dispatch
brief told dunyu-liu to keep its scratch work *inside its own worktree*
(`scratch/fenicsx_spike/`) to avoid colliding with other runs/ directories.
That was the wrong call: `isolation: "worktree"` auto-cleans a worktree with
no *tracked* diff, and gitignored/untracked scratch data does not count —
this is the exact failure rule 14 already names ("row 7 and row 16 runs were
lost this way... a run left inside a worktree goes with it"). The worktree
was removed around 13:19, mid-pipeline-rerun, taking the raw logs, fault
matrices (`G_*.npz`) and RSF results (`rsf_*.npz`) with it. dunyu-liu rebuilt
the verdict doc and scripts from its own console transcript (same numbers
it had already reported to the conductor) and placed them in the MAIN
checkout's `runs/20261009_fenicsx_spike/` (gitignored, no collision with
row 1's `runs/20261009_petsc-profile/`) — consistent with rule 14's actual
guidance that run *data* belongs in the main checkout's `runs/`, not inside a
worktree. **Lesson, generalizes beyond this project (candidate for consilium
inbox)**: brief every worktree-isolated agent whose deliverable is run output
(not a tracked code diff) to write that output directly to the shared `runs/`
location in the main checkout from the start, never to worktree-local scratch
— the auto-reap-on-no-tracked-diff behaviour makes "move it out before reap"
too late if the agent's own turn ends before the conductor can intervene.

Pending: land `runs/20261009_fenicsx_spike/fenicsx_feasibility_2026-10-09.md`
as `docs/dev/fenicsx_feasibility_2026-10-09.md` via a docs-only PR (fast lane,
no audit) once content is reviewed; prototype scripts stay in `runs/`
(untracked, per row 2's own "no new file until the answer says so").

## Roster of live/dispatched children
- dunyu-liu (agent id `aed14af1bc8baa76b`): worktree removed (auto-reap,
  see incident above); currently idle, awaiting go-ahead for Q3 timing.
- Row-1 profiling: run directly by the conductor (Bash, no subagent — no
  roster specialist fits per-rank MUMPS/PETSc timing; dunyu-liu's budget
  for this session is scoped to row 2 only, per owner note).

## Update 13:50 — host contention, step (b)/(c) scoping

Two more clean-rerun attempts of BP8-wide r8 (13:23-13:36, 13:36-13:50) both
landed at ~10.6-11.0 s/step, suspiciously matching the first "contaminated"
run. Root cause found: an unrelated 8-process job (`run_inversion.py`,
`work3d/3d_exp086_modelerr_twin_fig9data_hashimalyr`, fenicsx env, ~90% CPU
each, started ~13:38) is running on theo4 — not dunyu-liu's, not row 1's,
likely the human owner's own unrelated work on this shared host. Per rule 15,
not ours to kill. BP8-wide step (a) numbers stay **blocked** until the host
is genuinely idle again; bp1002's numbers (collected 12:43-12:45, before any
contamination) stand.

Scoping check for (b)/(c), no host time needed:
- **(b) MUMPS BLR**: not enabled anywhere in `src/solveTimeLoopMUMPS.f90` —
  only `ICNTL(3)` (output stream) is set today. `ICNTL(35)`/`CNTL(7)` are
  untouched. ParMETIS: `libparmetis.so` is present in the eqquasi-petsc conda
  env (no install needed), but MUMPS's `ICNTL(28)`/`ICNTL(29)` (parallel
  analysis / parmetis ordering) are likewise unset. Real code change, real
  work.
- **(c) hypre BoomerAMG vs GAMG**: `libHYPRE.so` is present, but
  `solveTimeLoopPETSc.f90:476-486` is a hard correctness whitelist — only
  `KSPPREONLY+PCCHOLESKY` (MUMPS) and `KSPCG+PCGAMG` at `ksp_rtol<=1e-12`
  are permitted; anything else (including `-pc_type hypre`) hits
  `MPI_ABORT` by design (rule 2, no silent unverified answers). Testing
  BoomerAMG is therefore not a runtime-flag experiment — it needs its own
  reference-parity verification (same rigor as row 8b's CG+GAMG gate) before
  the whitelist can be widened, i.e. a real `src/*.f90` change through the
  rule-3 gate, not a profiling run.

Neither (b) nor (c) is startable as a quick measurement; both are scoped
engineering tasks for a dedicated session. Not started this session given
budget and the live host contention.
