# Pathway forward

Open work only. Updated 2026-11-13, at v1.22.0.

Close an item by deleting it, and record what was learned where it will be
read again: `PROJECT_RULES.md` for a rule, `data/<bench>/README.md` for a
benchmark finding, a commit message for a fix. Closed rows and the findings
behind them (BP8 domain and boundary study, step-over science, caps, known
limits) are in `docs/dev/board_history.md`. Rules are in `PROJECT_RULES.md`;
build, run and traps are in `CLAUDE.md`.

---

## Open

1. [ ] **P1 -- PETSc scaling in EQquasi.** Owner, 2026-10-09. Known today:
   native MUMPS scales 1.9x over 1->8 ranks (BP5); 10 m BP8 wide box ~15 s/step
   at 8 ranks, ~79 GB of factors; CG+GAMG scales better but is slower at every
   size tried (1.4x bp1002, 7.7x kink300); faulting and BP8 pore pressure run
   on rank 0 only. Steps, in order:
   (a) profile per-step time (solve / faulting+pressure / communication / I/O)
   at 8, 16, 32 ranks on 10 m BP8 wide and bp1002, host to ourselves;
   (b) MUMPS block-low-rank compression, ParMETIS ordering, threaded BLAS per
   rank; (c) hypre BoomerAMG vs tuned GAMG with the setup reused across steps
   (K is constant); (d) distribute the rank-0 work only if (a) says it matters.
   Done when: 10 m BP8 wide <= 5 s/step at 32 ranks, or a measured reason it
   cannot be, with parity to the gold at every change.
   Evidence: the per-step table from (a) and the final run's runInfo.json.
   **ROW 1 SUMMARY FOR THE OWNER, 2026-11-14 (conductor 12; every figure is a table or paragraph below, none new; 8-16 ranks, MEASURED, hosts shared so loads vary run to run).**
   (Row 1 history referred to as "below" here: (archived: board_history.md, section "Row 1 measurement history, 2026-11-14 to 2026-11-18 (archived 2026-11-18)"))
   - *Passes the gold-style comparator (1e-9):* exact MUMPS (reference); all-primal substructuring (bp1002 25/25, bp5 120/120, gs10 33/33; wide not runnable, dense Schur 270-295 GB); GAMG CG at tight rtol (wide 1e-14: 0 of 33 fail; bp1002 1e-11/1e-12); BLR-preconditioned CG at tight rtol. *Fails it:* BLR 1e-8/1e-6/1e-4 as a solver (bp1002 13-23 of 25, bp5 115-118 of 120, wide 13-19 of 33), GAMG CG at rtol 1e-8 (wide 9 of 33, bp5 116 of 120), subdomain-direct CG and BDDC (9-15 of 25/33).
   - *Seconds per step (loop):* bp5 0.128-0.178 exact, 0.267 all-primal, 0.335 GAMG; gs10 16 ranks 2.71 exact, 4.45 all-primal; wide 16 ranks 6.9-8.3 exact (5.5 at 200 steps), 32.4 GAMG at rtol 1e-8, 49.9-53.4 GAMG at the passing rtol (6-8x slower), 2.9 BLR 1e-8 (1.9x faster, fails). Exact wide at 8/16/32 ranks: 12.7 / 6.9 / 7.1; the row's target (wide <= 5 s/step at 32 ranks) is not met by anything that passes.
   - *Setup and memory per rank (wide, 16 ranks):* exact 331-344 s and 12.1-12.3 GB; GAMG CG 10 s and 2.54 GB (4.8x less); BLR 128-156 s and 9.5 GB; all-primal on gs10 268 s and 20.7 GB against exact 81 s and 4.5 GB.
   - *Merge gate:* full `e2e` tier on the stack head, solver 3 forced: 90 passed, 45 skipped (10 pending bp1002 zone gold, 35 category-absent), 0 failed (`runs/20261114_gate-full/`). Audit of the stack diff: 0 BLOCKER, 0 MAJOR, 9 MINOR (all listed above). Stack #76, #77, #82, #83, #85 is MERGEABLE with CI green; #87 (BDDC spike) is not for merge.
   - *64-rank `-11` abort:* MUMPS -9, cured by `-mat_mumps_icntl_14 100` (3 of 3 pass against 5 of 6 fail; +20-30 % rank memory), passed through `PETSC_OPTIONS`; 60 of 60 clean at 8 ranks.
   - *Owner decisions open:* (1) merge the stack and in what order; what `par.solver` means at (h) and PETSc as a hard dependency; (2) default ICNTL(14) in source, and `KSPSetUp` reading INFOG(1)/(2); (3) move the comparator tests to `tests/unit` (new file) and document `PETSC_OPTIONS`; (4) audit items 2-6, 8, 9; (5) closing this row, or taking the measured reason: no tested solver passes at <= 5 s/step on wide at 32 ranks, and the solver survey is saturated (no new runs are queued).
   **(a), re-verified 2026-10-12** (binary eqquasi-1.21.2, host theo4,
   host-exclusive, 10-step profiling runs, native MUMPS `par.solver=1`;
   `runs/20261009_petsc-profile/{bp1002,bp8wide}-r{8,16,32}/result/cycle0/runInfo.json`):
   bp1002.qdc.2500 (127k eq) -- 8r 0.436 s/step, 16r 0.409, 32r 0.394 s/step
   (flat; rank-0-only faulting/pressure work dominates at this size). BP8 10 m
   wide (3.52M eq, scratch-only `fxmin/fxmax/fzmin/fzmax=-1500/1500`,
   `enlarging_ratio_xz=1.2`, NOT the committed `compset/bp8.qdc.gs.10` which
   is a smaller -500/500 box) -- 8r 11.01 s/step, 16r 8.01, 32r 6.91 s/step.
   Correction to a prior conductor's note: an earlier pass flagged these as
   possibly contaminated by a concurrent external job and claimed a rerun was
   needed; on re-running r16 and r32 clean (host confirmed idle, `uptime`
   load 0.35, no other eqquasi/python procs), the new numbers (8.01, 6.91)
   match the flagged-contaminated ones (7.98, 6.76) to ~2-4%, i.e. within
   normal run-to-run noise -- the contamination was not material here.
   (b)/(c)/(d) not started. Finer solve / faulting+pressure / communication /
   I/O split (as (a) asks) is NOT yet instrumented -- only aggregate
   time-loop and factorization time exist today (`solveTimeLoopMUMPS.f90`'s
   `cpu_time` calls). Recommendation (not acted on -- rule 1, new-compset
   naming needs agreement): register the -1500/1500 wide box as its own
   compset (e.g. `bp8.qdc.gs.10.wide`) instead of leaving it only in `scratch/`.
   **Owner direction, 2026-11-13**: "i think we can domain decompose from the
   beginning with petcs?"; on the three changes below, "all three changes
   sounds great". Why: the time loop is two triangular solves per step and
   those do not scale (above; Wang 2026 Case 3a, 3.61M eq, 24 ranks, v1.22.0:
   factorization 137 s, time loop 15093 s per cycle, ~6.5 s/step). The PETSc
   path today is not a domain decomposition: rows are split by `PETSC_DECIDE`,
   rank 0 builds the whole RHS and receives the whole solution
   (`VecScatterCreateToZero`), and the default solve is still MUMPS. The
   three changes, in the Fortran PETSc path, replacing (c)/(d):
   (e) each rank owns a block of the mesh: it generates only its own
   elements (the generator is analytic, no partitioner), assembles its own
   rows and RHS, and runs friction and pore pressure on its own fault nodes;
   (f) memory scales out: no rank holds the full mesh or the global factors
   (79 GB at 10 m BP8 wide, 141 GB in-core estimate at Case 3);
   (g) subdomain matrices kept unassembled (`MATIS`) and solved with
   `PCBDDC`, the setup built once and reused every step (K is constant).
   Owner, 2026-11-14, on the solver under (g): "This is a really simple system.
   It's symmetric. You can use maybe just a regular CG or GRMS." (verbatim;
   read as: try plain Krylov, CG or GMRES, on the distributed matrix, not
   only BDDC with direct local solves.)
   Owner, 2026-11-14: "If you cut the domain into subdomains, and each
   subdomain you can use mumps with much less memory cost, then it also
   sounds a good strategy." (verbatim; read as: a subdomain-wise direct
   solve under an outer Krylov iteration is also a candidate.)
   Owner, 2026-11-14: "you can just factorize once and just solve for each
   subdomain and then piece them together. But you just factorize once."
   (verbatim.)
   Owner, 2026-11-14: "was a really good strategy and put it on board as
   well. So in the end, we may have several solver options. And I think the
   last that you proposed, this 3-step one, is pretty good. Give it a try."
   (verbatim.) The 3-step strategy, exact substructuring with no outer
   iteration: (1) factorize each subdomain interior once; (2) build and
   factorize once the Schur complement that couples only the interface
   unknowns; (3) each step is one solve per subdomain, one interface solve,
   one more solve per subdomain. Known cost: the interface matrix is dense,
   its size not yet computed. Status: TO TRY as a spike, beside plain CG and
   subdomain-direct solves under CG; none measured yet.
   Owner, 2026-11-16, on the 3-step solve after its measured cost (dense
   interface, 3-4x the memory of exact): "how come? the interface shouldn't
   grow that crazy" and, on the untried ways to cut the dense interface cost
   (BLR compression of the Schur complement offered as the cheapest test on
   the existing spike; a sparse hierarchical interface factorization; fewer,
   larger subdomains): "keep pursuing that line" (verbatim.) Owner, 2026-11-16:
   "keep pursuing that line, optimize it with memory and speed." (verbatim.)
   Status: steps 1 and 2 measured 2026-11-16 (conductor 13, block below);
   the cheaper levers are in progress, each locked here before the next.
    (The "block below" is archived: (archived: board_history.md, section "Row 1 measurement history, 2026-11-14 to 2026-11-18 (archived 2026-11-18)"))
   Owner, 2026-11-16, on the parity tolerance for approximate solves (the
   question further down this row): "let's use 1e-9 to machine precesion for
   now" (verbatim.) Read as (a reading, not the owner's words): the 1e-9
   comparator stays the bar for every solver option for now, and iterative
   solves are run as tight as it takes to pass it, down to machine precision.
   Owner, 2026-11-16, after the survey of the other solver lines (GAMG CG
   best on memory and scaling; its untested lever, rigid-body modes per side
   of the fault's split nodes in the near-null space, offered as a second
   line): "so the winner is gamg cg. what's the near null space on the split
   node?" and then "keep proposing, keep exploring. 3 days budget. go"
   (verbatim.) Budget read as 2026-11-16 16:20 CST to 2026-11-19 16:20 CST;
   scope read as (a reading): the solver lines on this row, the 3-step solve
   and GAMG CG first, new proposals recorded here before they are run.
   Owner, 2026-11-17 (asked which clusters and allocations are reachable;
   two TACC systems are): the scale test is to run on Vista under an
   allocation the owner named (marked paraphrase; the allocation code is
   left out of this public file), followed by "just need more cores to
   scale?" (verbatim.) Read as (a reading): the cluster test under (g) is
   released on that one allocation; the build there is spike work until a
   machine target is agreed; the question is answered with measurements.
   Owner, 2026-11-17, on which solver options move from spike to a
   merge-bound branch: "let me decide promotion" (verbatim.) Read as (a
   reading): no agent promotes, proposes a merge-bound rewrite of, or
   retires a solver option; rankings are offered, the owner picks.
   Owner, 2026-11-17, shown a ranking whose item 1 was the distributed
   PETSc loop with exact MUMPS (`par.solver = 3`, the stack #76, #77, #82,
   #83, #85): "cool. promote 1" (verbatim.) Read as (a reading):
   `par.solver = 3` and the new file `src/solveTimeLoopPETScDD.f90` are
   agreed, and the stack lands through gated merges in stack order. Not
   read into it: a change of default solver, retiring solver 1 or 2,
   anything under (h), or promotion of any other option.
   Owner, 2026-11-17, asked for a better name for that option and offered
   three: "distributed direct" (verbatim.) Read as (a reading): the
   option is called distributed direct (distributed assembly, one global
   direct factorization) in docs, the board and PR titles from here on;
   renaming the source file or any parameter is not decided.
   Owner, 2026-11-17: "and let's /refactor the code base to make eqquasi
   naturally adopts more solver options on the line?" (verbatim.) Read as
   (a reading): after the stack lands, a behaviour-preserving refactor of
   the solver layer so a new solver option plugs into the distributed
   direct loop without a new time loop; gold results identical; the
   design goes on this row before code; no new user-facing parameter,
   no default change and no retired solver without the owner's word.
   Owner, 2026-11-17: "/release to lock this promotion. other lines still
   pending my decision and test later" (verbatim.) Read as (a reading):
   once the stack is on master, cut a release (audit first, tag gated on
   green CI for that commit) so the distributed direct option is locked;
   the release comes before the refactor; every other solver line stays
   spike-only and undecided.
   Owner, same exchange: "and this should go through full test" and
   "with a full cycle?" (verbatim.) Read as (a reading): the release is
   gated on the full `e2e` tier on the exact commit with the distributed
   direct option exercised, not only the fast tier; "full cycle" is read
   two ways and both are done unless the owner narrows it: the full
   milestone cycle (rules audit, code audit, fixes, release), and a run
   through a complete earthquake cycle compared against the gold result.
   Owner, same exchange: "also, always measure perfs" (verbatim.) Read
   as (a reading): every gate, release and refactor step on this row
   records performance beside correctness (setup, time per step, peak
   memory, rank count, host load), before and after the change.
   The owner then named the three: "memory cost, scaling with more cores,
   and compute speed per step" (verbatim.)
   Owner, 2026-11-17: "and whatever metrics you'd better collect, say
   machine used, hardware used, i think eqdyna has some setups you can
   borrow" (verbatim.) Read as (a reading): each performance record also
   names the machine and its hardware (CPU model, cores, memory, compiler,
   MPI and library versions), plus whatever else makes the number
   reproducible; look at how EQdyna records its runs and borrow that
   layout rather than invent one. A new file or script for it is still a
   rule 1 proposal on this row first.
   Owner, 2026-11-18 13:55 CST, answering three numbered questions put by
   the invoking session: "3A. 2A" (verbatim.) The options as they were put:
   question 2, MUMPS workspace setting ICNTL(14), "A) default it to 100,
   since Vista fails at the default from 64 ranks (recommended), B) leave
   at 35"; question 3, solver-layer refactor design on this row, "A) go,
   gated at 1e-9 (recommended, since byte-identical only holds at 1 rank),
   B) hold". Read as (a reading): (2A) ICNTL(14) defaults to 100 in source
   on the MUMPS-through-PETSc paths, a user option still overriding it; a
   `src/` change, so a version bump and a gated PR, with setup time, time
   per step and memory recorded before and after; whether the native MUMPS
   loop (solver 1) gets the same default, and `KSPSetUp` reading
   INFOG(1)/(2), are NOT decided. (3A) the refactor design on this row
   (#166) may start; its gate is the 1e-9 comparator against gold and
   against the pre-refactor binary, not byte identity; it stays
   behaviour-preserving, with no new user-facing parameter, file rename or
   default change without the owner's word. Question 1 (an absolute floor
   for the out-of-plane columns in the comparator) was NOT answered and
   stays open: strict 1e-9 holds until it is.
   Owner, 2026-11-18 13:57 CST, answering question 1: "1A" (verbatim.) The
   option as it was put: "Comparator bar for out-of-plane columns: A) add
   an absolute floor so roundoff passes (recommended), B) keep strict
   1e-9". Read as (a reading): the comparator keeps rtol 1e-9 on every
   column that carries signal and gains an absolute floor on the
   out-of-plane columns (`Slip_3`, `Slip_rate_3`, `Shear_stress_3` and
   their log10 forms in the `fltst` files), so values at roundoff level
   pass. The floor value was NOT given: it is proposed on this row with
   the measured exact-against-exact differences behind it before any
   result is re-read as a pass; the cross-architecture float32 `disp`
   files are NOT covered by this answer.
   Owner, 2026-11-18 13:59 CST, after being told what stands between the
   slip-to-traction option and regular use (and, one message earlier, the
   Vista limits of distributed direct): "fix those vista issue, and run
   full cycles, a few, and verify. if we can pass that and the CI, let's
   promote it" (verbatim.) Read as (a reading, the owner may correct it):
   "it" is the slip-to-traction option; promotion is CONDITIONAL on
   (i) the Vista issues fixed: ICNTL(14) (2A above), the NUL holes in
   output written to NFS scratch (run from Lustre, or a fix if the cause
   is in the code), and for slip-to-traction a clean dedicated-node run;
   (ii) a few full earthquake cycles run with slip-to-traction and
   verified against exact at the comparator bar (1e-9 with the 1A floor);
   which cases and how many is proposed on this row first, with cost;
   (iii) the full test tier and CI green on the exact commit. Promotion
   then follows the same route as distributed direct: gated merges, a
   version bump, the code audit before any tag. NOT decided by this
   message: the option's name and parameter value (rule 1: proposed on
   this row and agreed before the merge-bound branch), a `vista` machine
   target, a compset for wide, the re-layout reader for the saved
   operator, any default change.
   Owner, 2026-11-18 14:00 CST: "then, let's promote slip to traction"
   (verbatim.) Read as (a reading): confirms that "it" above is
   slip-to-traction and that promotion comes after the conditions (i) to
   (iii) pass ("then"), not before them.
   Owner, 2026-11-18 14:00 CST: "and then optimize the slip to traction"
   (verbatim.) Read as (a reading): after the promotion lands, the
   slip-to-traction option is optimized; the measured costs on this row
   are the candidates (operator build time, 1.9 h at 32 ranks and 3 h at
   16 on wide; memory, about 1.3x summed and 1.5x peak rank against
   exact; the saved operator being readable only at its own rank count).
   Each optimization is proposed on this row with its before and after
   measurement and stays behaviour-preserving at the comparator bar; it
   comes after the promotion, not inside it.
   Owner, 2026-11-18 14:00 CST: "/autopilot go" (verbatim.) Read as (a
   reading): the unattended run continues on the scope above (Vista
   fixes, full cycles and verification, promotion of slip-to-traction if
   the conditions and CI pass, then its optimization, with 2A, 3A and 1A);
   "go" names no time budget, so it is read as until this row's items are
   done or blocked on the owner, replacing the 2026-11-19 16:20 CST end.
   The Vista ceiling of 150 SU and every "NOT decided" item above stand.
   Owner, 2026-11-18 16:18 CST, answering three numbered questions put by
   the invoking session: "1A. 2A. 3A" (verbatim.) The options as put were:
   question 1, refactor file name, "A) accept `src/onFaultHelpers.f90`
   (recommended), B) fold the helpers into an existing module"; question
   2, slip-to-traction option, "A) `par.solver = 4`, "slip to traction"
   (recommended), B) another value or name"; question 3, out-of-plane
   floor value, "A) 1e-12 m/s on `Slip_rate_3` and 1e-9 MPa on
   `Shear_stress_3`, with `Slip_3` measured first (recommended), B) other
   values". Read as (a reading): the new file `src/onFaultHelpers.f90`
   with module `onFaultHelpers` is agreed for refactor step 1 (rule 1);
   the slip-to-traction option is `par.solver = 4`, named "slip to
   traction", so its promotion stack may go on merge-bound branches; the
   comparator's absolute floors are 1e-12 m/s for `Slip_rate_3` and
   1e-9 MPa for `Shear_stress_3`, and the `Slip_3` floor is measured and
   put on this row before any value is used. NOT decided by these
   answers: the `Slip_3` floor value; how the floor applies to the log10
   columns; any new source file name the slip-to-traction stack may need;
   ICNTL(14) for native solver 1; a `vista` machine target; any default
   change.
   Owner, 2026-11-18 16:21 CST, after being told that the solver runs in
   double precision but the `disp.*.nc` snapshots are stored as 4-byte
   reals (`NF90_REAL` in `src/netcdf_io.f90`), which is why they fail the
   1e-9 comparator across machines: "weird. fix that. should be double
   precision" (verbatim.) Read as (a reading): `netcdf_write` stores the
   displacement snapshots as 8-byte reals, as a gated `src/` change with a
   version bump, performance and file size recorded before and after.
   NOT decided by this message: the same change for the other two
   single-precision writers in that file (`netcdf_write_on_fault`, whose
   `fault.*.nc` files are in the committed gold and are read back on
   restart, and `netcdf_write_roughness`); any regeneration of gold under
   `data/`; what readers outside this repo expect.
   Owner, 2026-11-18 16:22 to 16:24 CST, on the same subject: "there
   should be a setup of precision in globalvar.f90?" and "all f90 should
   follow that right?" (verbatim.) Told that `dp =
   selected_real_kind(15,307)` in `src/globalvar.f90` is followed
   everywhere except three `NF90_REAL` definitions in `src/netcdf_io.f90`
   and a few hard-coded but equivalent declarations, and asked two
   numbered questions, the owner answered: "for 1 and 2, yes, both should
   follow dp" (verbatim.) The options as put were: question 1, on-fault
   and roughness NetCDF writers, "A) switch them to double as well and
   regenerate the affected gold in a separate gated PR (recommended,
   since it matches your rule that everything follows `dp`), B)
   displacement only, leave those two in single"; question 2, the
   hard-coded `real(8)` and `double precision` lines, "A) change them to
   `kind=dp` in the same PR, behaviour-preserving (recommended), B) leave
   them". Read as (a reading): all three NetCDF writers store reals at
   the kind `dp` names, and the hard-coded `real(8)` and `double
   precision` declarations become `kind=dp`; the gold under `data/` that
   the on-fault change invalidates is regenerated in its own gated PR,
   with the old-against-new difference shown to be storage truncation
   only before any gold file is replaced. This supersedes the "NOT
   decided" items of the 16:21 paragraph for the two other writers and
   for that gold regeneration. NOT decided: the MPI calls that name
   `MPI_DOUBLE_PRECISION` (not asked); gold for anything other than the
   files this change invalidates; what readers outside this repo expect.
   Owner, 2026-11-18 16:25 CST, two messages in a row: "and rebless the
   affected disp" then "in ref" (verbatim.) Read as (a reading): the
   reference data that holds displacement snapshots affected by the
   change to 8-byte storage is regenerated ("reblessed") with the new
   binary, in the gated gold PR above, each replaced file shown to
   differ from the old one by storage truncation only. "ref" is read as
   the reference (gold) sets under `data/`; which sets hold `disp.*.nc`
   is listed on this row before any is replaced. NOT decided: adding
   displacement snapshots to a reference set that holds none today.
   Owner, 2026-11-18 16:36 CST, answering two numbered questions put
   after being told that no reference set holds `disp.*.nc` and that the
   bp1002 reference lacks the second-fault summary files: "1A. 2A."
   (verbatim.) The options as put were: question 1, displacement in the
   reference, "A) add the last-step snapshot per reference set
   (recommended), B) add every snapshot, C) keep none"; question 2,
   bp1002 second-fault summaries, "A) add them in the same gold PR
   (recommended), B) leave them out". Read as (a reading): when the gold
   is regenerated for the 8-byte storage, each reference set gains its
   last-step `disp.*.nc`, and the bp1002 set gains
   `cplot_ft2_EQquasi.txt`, `cplot_ft2_ruptarea_trac_slip.txt` and
   `cplot_ruptarea_trac_slip.txt`, all in that one gated gold PR, with
   the size the tree grows by recorded. This settles the "NOT decided"
   item just above. NOT decided: any other file absent from a reference
   set (`p1output.txt`, `p2output.txt`, `peak_sliprate_per_fault.dat` on
   bp1002); the tolerance the tier applies to the new files.
   Owner, 2026-11-18 16:38 CST: "/autopilot clear the board. a week
   budget. go" (verbatim.) Read as (a reading): the unattended run now
   covers every open row of this board, worked in `prio` order with
   state as the tiebreak, not only row 1; the budget is one week from
   this message, to 2026-11-25 16:38 CST, replacing the "until this
   row's items are done or blocked on the owner" end; the order already
   set on this row (refactor step 1, NetCDF output follows `dp`, gold
   regeneration, slip-to-traction promotion, then its optimization)
   keeps its place at the head of the queue. Unchanged: the Vista
   ceiling of 150 SU, every "NOT decided" item above, rule 1 for any new
   file, name, parameter or compset, no major version bump unattended,
   and that opening, closing and re-scoping rows stays with the rules
   owner.
   Owner, 2026-11-18 16:39 CST, answering one numbered question on the
   tolerance left open by the 16:36 paragraph: "1A" (verbatim.) The
   options as put were: "Bar for the double-precision NetCDF files in the
   test tier: A) 1e-10, same as the text columns (recommended), B) keep
   1e-6". Read as (a reading): once `disp.*.nc` and `fault.*.nc` are
   stored as 8-byte reals, the test tier compares them at a relative
   1e-10, the bar the double-precision text columns already use
   (`DD_RTOL`), in place of the 1e-6 kept for 4-byte storage
   (`DD_DISP_RTOL`). NOT decided: any change to the 1e-9 campaign
   comparator; what to do if a measured difference exceeds 1e-10 (it
   comes back to this row, the bar is not loosened unattended).
   Owner, 2026-11-18 16:45 CST, after being told that the 8-byte output
   fails the CI comparison against the 4-byte reference files: "though
   not bit identical, regenrated golds should match single precision with
   a relatively larger bar" (verbatim.) Read as (a reading): an
   acceptance condition on the gold regeneration, not a change of the
   1e-10 tier bar: every regenerated reference file must agree with the
   file it replaces at single-precision level (a looser bar, of the
   order of the 1e-6 the tier used for 4-byte storage), and the measured
   worst difference per file is recorded with the regeneration; a file
   that does not agree at that level is not replaced and comes back to
   this row. The stricter check already planned (new values cast to 4
   bytes equal the old values) satisfies it. NOT decided: the exact
   value of the looser bar (1e-6 is proposed here, as the bar already in
   the tier).
   Owner, 2026-11-18 16:45 CST, next message: "physics shouldn't change"
   (verbatim.) Read as (a reading): the stated reason for the condition
   above and a requirement on the whole change: storing the output at
   8 bytes may alter only the stored precision of each value, never the
   computed solution; the double-precision text output of the same case
   (station series, summary files) must stay within the 1e-10 tier bar
   between the binary before and the binary after, and any difference
   beyond storage truncation stops the gold regeneration and comes back
   to this row.
   Owner, 2026-11-19 17:07 CST, replying to the daily brief of 09:37 CST:
   "go" (verbatim.) The brief's three questions were: 1, slip-to-traction
   self-check floor (1e-14 of the operator's largest entry): "A) adopt,
   with a test option `-eq_s2t_plant_row` proving a wrong operator still
   aborts (recommended), B) no floor, which aborts bp7"; 2, round-off
   floors in the comparator: "A) `Slip_3` at 1e-15 m, the same floors
   for surface-station files, applied to un-logged values (recommended),
   B) no additions"; 3, already landed, confirm: "A) accept
   (recommended): replaced on-fault gold judged at the tier's 1e-6 floor;
   restart files unrounded, so cycle 2 differs by about 1e-6; B)
   revisit". Read as (a reading): the recommended option of each, that
   is owner-brief items 3 and 7 (floor adopted together with the
   planted-operator option `-eq_s2t_plant_row` and its run-level test),
   items 2 and 8 (`Slip_3` floor 1e-15 m; floors applied to the
   un-logged values, a log10 entry below the floor skipped), items 1
   and 6 (accepted). The surface-station part of question 2 is read as
   superseded: after the brief the invoking session told the owner that
   the bp7 8-rank failure did not reproduce in 7 repeats and no floor
   is needed for surface-station files, so item 10 is read as option D
   (no `srfst_` floor; re-run that test in the merge-head gate, and on a
   second failure keep the case and send it to the code auditor first).
   NOT decided by "go": item 11 (liu2020 reference, including any looser
   text bar for that set) and item 5 that depends on it; item 9
   (ICNTL(14) for native solver 1); step 2 of the operator-build plan (a
   new option name); any default change.
   Owner, 2026-11-19 17:20 CST, after being shown that the v1.24.7 and
   v1.24.8 liu2020 runs overlay at plot scale and differ only during the
   rupture (at most 3.6e-7 of a column's maximum), asked "so 1e-9 rel is
   not a good metrics?", and given the recommendation: "Keep 1e-9 / 1e-10
   for short regression cases. For full-cycle cases, compare event time,
   peak slip rate and event slip and moment at a stated tolerance. If you
   agree, the conductor would propose the metrics and tolerances on the
   board, measured on the existing full-cycle pairs, before any
   comparator change." The owner answered: "yes" (verbatim.) Read as (a
   reading): the 1e-9 / 1e-10 relative value bars stay for short
   regression cases; for full earthquake-cycle comparisons the metrics
   become physical ones (event times, peak slip rate, coseismic slip and
   moment per event); the conductor first proposes the exact metrics,
   their definitions (for example the slip-rate threshold that defines an
   event) and tolerances on this row, measured on the full-cycle pairs
   already run (liu2020 v1.24.7 against v1.24.8; the bp5, bp7 and bp1002
   full cycles of solver 4 against solver 3), and no comparator or test
   changes before the owner sees that proposal. NOT decided: the metric
   set, the definitions and tolerances, which cases they apply to, and
   item 11 (liu2020), which this proposal informs.
   Owner, 2026-11-20 15:58 CST, replying to the daily brief, whose
   question 1 was "Option name for a test that plants a NaN in slip to
   traction: A) `-eq_s2t_plant_nan` (recommended), B) another name":
   "1A." (verbatim; the same message asked why the liu2020 displacement
   file would be exempted and what ICNTL(14) is for, and answered
   nothing else.) Read as (a reading): the PETSc option
   `-eq_s2t_plant_nan <global row>` is agreed (rule 1), as conductor 37
   proposed it (m4, option A): it writes a NaN into that row of the
   operator output so a run-level test shows the non-finite abort
   fires; a gated `src/` change with a version bump and one e2e test.
   NOT decided: item 11 and item 5 (liu2020), item 9 (ICNTL(14) for
   native solver 1).
   Owner, 2026-11-20 16:00 CST, after being told that rule 22 caps every
   tracked file at 5 MB (enforced by `tests/contract/test_root_layout.py`,
   `MAX_TRACKED_FILE_BYTES`) and that liu2020's last-step displacement
   file is 6.3 MB, with options A) exempt that one file, B) commit
   liu2020 without a displacement file (recommended at that point), C)
   store displacement files compressed: "you can go above 5mb. use 10mb"
   (verbatim.) Read as (a reading): the per-file cap of rule 22 becomes
   10 MB, in the rule text and in the contract test, a rules-owner edit
   made on the owner's word; with it, item 5 needs no exemption, and the
   liu2020 regeneration of item 11 option A goes ahead from the v1.24.8
   run with its last-step displacement file, its README's "deterministically
   reproducible" claim corrected. The item 11 part is a reading of a
   message that answered the exemption question, not an explicit "2A";
   the invoking session states it to the owner when relaying. NOT
   decided: item 9; any further cap change; compression.
   Owner, 2026-11-20 16:01 CST, after being told what ICNTL(14) is (the
   extra MUMPS working memory, in percent, reserved beyond MUMPS's own
   estimate; already 100 on the PETSc paths, solvers 2 to 4; native
   solver 1 does not set it and no failure has been seen there; the
   recommendation was to leave it): "just make them 100?" (verbatim.)
   Read as (a reading): item 9 is decided as option B: native solver 1
   also sets ICNTL(14) = 100, so all solvers reserve the same MUMPS
   workspace; a gated `src/` change with a version bump that changes the
   default solver's memory reservation only, shown by gold results
   unchanged at the tier bar and peak memory per rank recorded before
   and after on the same cases. The message ends in a question mark; the
   invoking session tells the owner it is read as an instruction.
   Owner, 2026-11-23 11:39 CST, replying to a brief whose questions were:
   "1. Unmerged worktrees and branches: A) delete those already on
   master or marked superseded, copying cited runs out first
   (recommended); B) keep. 2. 27 uncited run folders and 39 scratch
   entries (about 7.7 GB): A) delete after a final citation check
   (recommended); B) keep. 3. Restart the unattended run in a fresh
   session (this one is past its context limit), with the week's budget
   ending 2026-11-25 16:38: A) yes, from `scratch/_logs/BRIEF_20261119.md`
   (recommended); B) wait." Answer: "1A. 2A, 3A" (verbatim.) Read as
   (a reading): (1) delete the worktrees and branches whose content is on
   master or that this board marks superseded, after copying any cited
   run output into `runs/`; anything neither on master nor marked
   superseded stays and is listed again; (2) delete the uncited `runs/`
   folders and `scratch/` entries listed by the 2026-11-20 cleanup report
   after a final grep of this board, `docs/dev/board_history.md` and
   `docs/` shows no citation; (3) the unattended run restarts from a fresh
   session seeded by that brief file, under the existing budget.
   Owner, 2026-11-23 11:42 CST, after "go" (11:41) on the restart: "no
   budget. keep working. /autopilot" (verbatim.) Read as (a reading):
   the run has no end date: the 2026-11-25 16:38 CST end is removed and
   the run continues until the board's rows are done or blocked on the
   owner ("until the board is green"); every other limit stands (Vista
   ceiling 150 SU, no major version bump unattended, rule 1, owner-held
   items).
   Owner, 2026-11-23 11:56 CST, on row 3 (BP8 discrepancy), after being
   offered A) compare EQquasi's slip-to-traction operator with the
   FEniCSx prototype's on BP8, B) audit the `src/qdct*.f90` force
   assembly, C) leave it parked: "bp8 hunt: it is OK. bc are different.
   pause it" (verbatim.) Read as (a reading): row 3 is paused, no work
   on it; the owner's view is that the difference from the reference
   comes from different boundary conditions, not a code bug. Row 4,
   which waits on row 3, stays waiting. NOT decided: whether row 3 is
   closed or rewritten (opening, closing and re-scoping rows stays with
   the rules owner, on the owner's word), and whether the uploaded BP8
   entries need a note.
   Owner, 2026-11-23 11:58 CST, on row 6 (time-step factor xi), after
   being offered "A) keep `par.xi_gamma` opt-in, add a warning when a
   constant `par.xi` exceeds the limit, and document it, with no default
   change (recommended); B) make per-node xi the default, which changes
   results for every case and needs new gold; C) leave as is": "1B"
   (verbatim.) Read as (a reading): the per-node time-step factor of
   Lapusta et al. (2000) becomes the default (`par.xi_gamma` default
   pi/4 in `scripts/defaultParameters.py`, today 0.0, with `par.xi` 0.2
   kept as the cap); this is an intended change of results, so every
   affected reference set under `data/` is regenerated and the tier
   passes against the new gold; each regenerated set records step
   count, run time and the event times and peak slip rates against the
   old gold, so the owner sees what moved. The owner's 2026-11-18
   "physics shouldn't change" was about 8-byte storage and does not
   apply. A changed default is a behaviour change: the version bump is
   at least minor (never major unattended), with a release-notes entry.
   NOT decided: whether any compset keeps a constant xi explicitly (for
   example a SEAS benchmark whose submitted results used xi 0.2); the
   value of `xi_gamma` other than pi/4.
   Owner, 2026-11-23 12:00 CST, on question 1 of the follow-up list,
   which read "Benchmark compsets under the new default (follows from
   1B): the SEAS cases (BP5, BP7, BP8, BP1002) were submitted with
   constant `xi`. A) Pin those compsets to their submitted constant
   `xi`, so the repo still reproduces the uploads; the new default
   applies to everything else (recommended). B) Switch them too and
   accept that the repo no longer reproduces the submissions.": "1:
   keep old xi for now" (verbatim.) Read as (a reading, the narrow
   one): the SEAS benchmark compsets keep their constant `xi`
   (`par.xi_gamma = 0` set explicitly in them) and their gold is not
   regenerated for this change. The words could also mean that the
   old constant `xi` stays the default everywhere for now, which would
   suspend the 11:58 "1B"; the invoking session asks the owner which,
   and until the answer NO default is changed and NO gold is
   regenerated for row 6. Either reading keeps the benchmark compsets
   on constant `xi`.
   Owner, 2026-11-23 12:01 CST, two messages clarifying the 12:00
   answer: "but also build the new xi gamma as refs as wel" and "then
   we can phase out" (verbatim.) Read as (a reading), replacing the two
   readings of the 12:00 paragraph: a transition in three steps. (1)
   Now: the constant `xi` stays the default and the existing reference
   sets stay as they are; nothing is regenerated in place. (2) Also
   now: a second set of references is built with the per-node factor
   (`par.xi_gamma = pi/4`, `par.xi` 0.2 as the cap) for the same cases,
   kept beside the existing ones, each recording step count, run time,
   event times and peak slip rates against its constant-`xi`
   counterpart. (3) Later, on the owner's word: the constant-`xi`
   references and default are phased out and the 11:58 "1B" default
   change takes effect. NOT decided: where the second set lives and
   what it is called (new names under `data/`, rule 1: proposed on this
   row first); which cases it covers (all nine sets, or the test-tier
   ones only: the liu2020 cycle costs about 5 h); whether the test tier
   runs the second set; when step 3 happens.
   Owner, 2026-11-23 12:02 CST, after being told that the FEniCSx BP8
   prototype differs from EQquasi mainly in its time integration
   (adaptive RK45 at tolerance 1e-6 against EQquasi's time loop with the
   `xi` step control) and in precomputing a slip-to-traction matrix:
   "eqquasi can have rk45" (verbatim.) Read as (a reading): the owner
   agrees to adding an adaptive Runge-Kutta 4(5) time integrator to
   EQquasi as an option beside the existing time loop, not as a default
   change. It is new work with no row yet: the rules owner opens a row
   on this word, and a design goes on the board before any code: how
   the integrator calls the elastic solve (several evaluations per step,
   cheap as operator products under `par.solver = 4`, a full solve each
   under the other solvers), the state it integrates, the error norm
   and tolerance, restart and output cadence, and the test (BP8 against
   the reference and against the prototype; a short case against the
   existing loop). NOT decided: the parameter name and values that
   select it (rule 1); whether it is limited to `par.solver = 4`; its
   priority against the rows already open; who implements it (a
   research-grade implementation needs the owner's word per dispatch).
   (h) end state, owner 2026-11-13: "in the end, mumps is a choice o petcs".
   One PETSc time loop; MUMPS, BDDC and AMG are solver choices inside it and
   `solveTimeLoopMUMPS.f90` retires. Before it retires: PETSc-with-MUMPS
   matches the gold on every gated compset and is no slower per step than
   native MUMPS; PETSc builds on every supported machine (today only
   `conda-linux` links it); one assembly path serves both `MATIS` and AIJ;
   what `par.solver` then means is agreed first (rule 1).
   Two design principles for all of (e)-(h), owner 2026-11-13: "things
   parallelize from beginning, and thing that can armotize once should be
   done that way".
   Parallel from the beginning: no rank-0 phase anywhere -- mesh, numbering,
   assembly, initial conditions, friction, pore pressure, station and fault
   output are all per-rank.
   Amortize once: whatever does not change is built once and reused -- the
   solver setup across steps (K is constant); mesh, partition and setup
   across the cycles of a run (each cycle rebuilds and refactorizes today,
   137 s per cycle at Case 3); and, as a candidate to measure, the fault
   slip-to-traction operator (row 4's note), which would amortize the volume
   solve itself across steps and across runs on one mesh and modulus.
   Design choices shown to the owner 2026-11-13, answer "go":
   a new distributed time loop beside the two existing ones, which retire
   only at (h), so the gold gate stays meaningful throughout; a regular
   x-y-z block decomposition from the analytic generator, not `DMPlex`;
   split fault nodes may cross subdomain boundaries, with a ghost exchange
   per step for friction; the fault slip-to-traction operator is a measured
   side experiment, outside (e)-(h); an iterative solve is held to parity
   within a stated tolerance, shown not to move event times. Still to
   agree at (h): what `par.solver` means, and PETSc as a hard dependency.
   Evidence: parity to the gold at each of (e)-(h), and the per-step table
   against native MUMPS at 8/16/32 ranks on bp1002 and 10 m BP8 wide.
   Unmeasured: whether BDDC beats MUMPS at these sizes (CG+GAMG did not).
   Host, owner 2026-11-13: "do it on cotopaxi" (64 cores, the
   `eqquasi-petsc` env is on the shared home). Timings there are not
   comparable with the theo4 numbers above: re-measure native MUMPS on
   cotopaxi as the baseline first.
   **Archived measurement history (row 1, 2026-11-14 to 2026-11-18).** The state tables, the per-step solver comparisons and the conductor 15-30 entries are in `docs/dev/board_history.md`, section "Row 1 measurement history, 2026-11-14 to 2026-11-18 (archived 2026-11-18)". What they established (figures copied from that text; MEASURED, none VERIFIED):
   - Exact 3-step substructuring passes the 1e-9 comparator (bp1002 25/25, bp5 120/120, gs10 33/33) and does not beat exact MUMPS: bp1002 at 8-32 ranks 3.4-4.1x the memory, 2.3-3.6x the setup, 1.7-2.3x the loop time; gs10 at 16 ranks 4.6x peak memory (20.7 against 4.47 GB), 3.3x setup (268 against 81 s), 1.6x loop (4.45 against 2.71 s per step).
   - BLR on the interface, deluxe scaling and two-level BDDC closed negative or gained 5-40 % of the gap only; CG + PCBDDC on wide at 16 ranks: 339.7 s/step against 7.18 s/step for MUMPS (~47x slower).
   - GAMG CG: loop ratio 2.4x slower than exact at 32 ranks on wide (8.74 against 3.63 s/step); on bp8.qdc.gs.10 at 16 ranks loop 49.3 s against 8.1 s (6.1x). BLR as a CG preconditioner does not win.
   - Out-of-plane roundoff at 400 steps is not s2t-specific: the three out-of-plane files fail against exact also between two exact runs differing only in rank count; every in-plane file passes at 1e-9.
   - Vista job 1152125 (spike-s2t d78d822, version 1.24.14, 4 nodes on gh-dev) COMPLETED 00:42:27, exit 0:0.
   Open items and every Owner paragraph of the archived range follow.

   #### Owner paragraphs and open items kept from the archived range

   Source: line 640, conductor 15 hand-off entry of the c83b771 board (kept verbatim, not edited).
    *Owner item (conductor 25, 2026-11-18): agent type names in older logs under `docs/dev/` (grep of the type names, line counts): `SESSION_LOG_2026-10-12_board-sweep.md` 15, `board_history.md` 3, `fenicsx_feasibility_2026-10-09.md` 1; `docs/user/gen_params.py` lists names only as a scrub regex (not a doc, kept). Rewriting past logs is the owner's call; not edited. This PR replaced the names in this board (hand-off text of conductors 15/16) by role names; git history is not rewritten.*

   Source: line 699, design block of the substructuring line of the c83b771 board (kept verbatim, not edited).
    **Owner-held, not touched by any step:** any new parameter name or value (including the spike `-eq_*` options and final names), retiring solver 1 or 2, renaming `solveTimeLoopPETScDD.f90`, what `par.solver` means after the merge, the default solver, relaxing the guard at :557-562 (this decides whether single precision, BLR and CG options are reachable at all), PETSc as a hard dependency (the makefile builds it only on `conda-linux`), promoting any option past spike.

   Source: lines 1041-1047, spike-run reporting defect and its owner questions of the c83b771 board (kept verbatim, not edited).
     Spike-run reporting defect, 2026-11-14 (audit): `runInfo.json` of the GAMG runs reads `ksp_type` preonly, `pc_type`
     cholesky (g64b checked), not the solver used; the real iteration counts (`ksp_iterations_*`) are right. Merge-bound
     solver 3 path: `solveTimeLoopPETSc.f90:496-509` sets the names from the whitelist booleans (`isPreonly`,
     `isCholesky`) and the guard at `:472` rejects other combinations, so it names what it ran by construction; read from
     source only, the merge-bound binary was not run on a GAMG option. Not changed (source change, owner). Owner
     questions, not decided here: a test pinning the INFO(2) meaning (new test); `KSPSetUp` reading INFOG(1)/(2) on
     failure and any default ICNTL(14) (changes merge-bound solver behaviour).

   Source: line 1148, conductor 24 entry of the c83b771 board (kept verbatim, not edited).

Owner question: a default ICNTL(14) above 35 for solver 3 at high rank counts changes a solver default, so it is the owner's; not changed.

   Source: lines 1182-1184, conductor 25 open questions of the c83b771 board (kept verbatim, not edited).
Open questions:
- cplot 9.8e-4 in the step-0 rerun: `runs/20261118_step0/rerun`, `test.bp8.qdc.gs.10_s3_np4_rerun` and `_np8_rerun` against their first runs, column 4 (index 3), entry 277770180767.1491 against 277770180767.1501, 2.4e-15 relative to the column maximum (about 16 ulp); no other file pair exceeds 5e-7 absolute, all at most 3.3e-15 of column scale; np1 reruns identical (bp1002, bp5, bp7, bp8 np1: 0 difference). Roundoff from reduction order at np greater than 1, not noise of another kind.
- local wide reference 5.45 s/step (218.1 s loop, 40 steps, 32 ranks, `runs/20261118_vista/local/wide40_np32`): host load 5.8 at start, 38.3 at end (`load_before.txt`, `load_after.txt`), so the second half ran loaded; one run, no repeat at stable load. the 3.3-3.5 s/step exact figure for wide at 32 ranks IS on this row (timing pairs at 32 ranks, "is met by exact itself in these runs (3.3-3.5 s/step)", `runs/20261117_c17-*`; my earlier sentence that it was absent was wrong, corrected here after the coordinator pointed to it). Restated: exact at 3.3-3.5 s/step in those batches against 5.45 s/step in the 2026-11-18 reference (same case, same ranks); host load (5.8 at start, 38.3 at end of the reference run) is a plausible cause, untested. UNEXPLAINED.

   Source: line 1190, conductor 25 entry of the c83b771 board (kept verbatim, not edited).

Owner items added: (1) wide cplot files lose rows on Vista NFS (the output path writes through rank 0; whether to write elsewhere or sync is the owner's); (2) cross-architecture bar for float32 `disp` files (1e-9 fails at 2.4e-8 of scale); (3) the comparator bar for out-of-plane columns, now with Vista evidence; (4) default ICNTL(14): default 35 failed at np64, 288, 576 on Vista and 100 passed, and passed at np32 and np144 at 35.

   Source: line 1457, conductor 29 entry of the c83b771 board (kept verbatim, not edited).
**OWNER ITEM, holds the merge:** the design on this row says step 1 moves the helpers into "a new or existing module" and names no file or module; `PROJECT_RULES.md` rule 1 requires agreement before a new file or name. 3A approved the design, not this name. Question for the owner: A) new `src/onFaultHelpers.f90`, module `onFaultHelpers`; B) an existing module (no new file; the makefile dependency order changes instead). No refactor PR is opened and v1.24.7 is not merged until answered; branch `worktree-agent-a42a732ee98ad8c43` is pushed, head 6b8693d, worktree `.claude/worktrees/agent-a42a732ee98ad8c43` kept.

   Source: line 1466, conductor 29 hand-off of the c83b771 board (kept verbatim, not edited).
OWNER ITEMS OPEN: (1) the s2t option's name and parameter value; (2) the file name for refactor step 1, A) new `src/onFaultHelpers.f90` or B) an existing module (section above; branch `worktree-agent-a42a732ee98ad8c43` head 6b8693d gated, held); (3) the 1A floor value; (4) ICNTL(14) for native solver 1 and `KSPSetUp` reading INFOG, unchanged.

   Source: line 1490, conductor 30 entry of the c83b771 board (kept verbatim, not edited).
Proposal 14 gaps still OPEN, unchanged: (a) no `disp.*.nc` in the committed gold, so the s2t displacement field is not compared on any case; (b) bp1002 second-fault `cplot_ft2_*` and `cplot_ruptarea_trac_slip.txt` have no gold file; (c) s2t-only extras have no gold file. Proposal 13 job 1152125 (started 16:02 CST on gh-dev) is not harvested; the invoking session relays its exit. Owner items open as in the conductor 29 hand-off.

3. [ ] **P1 -- EQquasi BP8 discrepancy (bug hunt).** Owner, 2026-10-09: "hunt
   the bug". An independent FEniCSx linear-hex model with EQquasi's box and BCs
   converges to Kim within ~2 % at 25 m (`docs/dev/fenicsx_feasibility_2026-10-09.md`);
   EQquasi does not (depth station -18 % at 50 m, -13 % at 10 m). Pressure
   matches; domain and dtmax ruled out. The uploaded BP8 entries are suspect
   along depth until this closes. Done when: the cause is found and fixed, and
   EQquasi at 25 m agrees with Kim within a few % at all 9 stations.
   **Escalation, 2026-10-12**: 6 code-audit findings on this hunt
   rejected after a source/gate check (4 before today: nodal area/dx^2, GS
   z-symmetry, and two more per the 2026-10-09 brief; 2 today, PRs #58 and
   #59, closed unmerged -- see `docs/dev/SESSION_LOG_2026-10-12_board-sweep.md`
   for each). Today's two both targeted `src/faulting.f90`'s split-node
   traction formula (lines 112-128): #58 claimed "ruled out" without ever
   computing the required EQquasi-vs-FEniCSx ratio; #59 computed only the
   FEniCSx half and flagged the mass-weighted cross term as a "structural
   bug," which does not survive the rule-16a sanity check -- that formula
   sits outside every `bp ==` branch (the only ones in the file are at
   lines 159/162/359, pore-pressure/BP7-specific), so a bias large enough to
   explain BP8's -13%/-18% would also break BP5/BP1002's passing regression
   gold, and it doesn't. The split-node coupling itself is very likely not
   the bug (3 independent rejected attempts have now landed on it from
   different angles). **Needs owner diagnosis, not a 3rd same-session
   retry**: live lead is still unconfirmed-but-unrefuted (tau/pressure match
   Kim within ~0.5% but V amplifies a small tau gap exponentially); next
   untried candidate is `src/qdct*.f90`'s element-to-fault force
   contribution (assembly path, not yet audited by anyone).
4. [ ] **P2 -- FEniCSx backend, EQdyna pattern.** Owner, 2026-10-09: "like
   eqdyna, we will have a fortran and fenicsx backend now ... follow eqdyna's
   success". One compset and one gold set, a test matrix that runs each case
   on each backend, a backend counts only at parity with the gold. Order:
   after row 3 (the Fortran gold must be right); BP5 and step-over parity
   first, BP8 against Kim. The prototype's precomputed slip-to-traction matrix
   (quasi-static, linear) is also a scaling idea for row 1.
6. [ ] **P2 -- Time-step factor xi above the stability limit.** Owner,
   2026-11-13: "that ksi lesson, put on board". Lesson: a constant `par.xi`
   above the Lapusta et al. (2000) eq. 15 limit makes the time step
   oscillate between events. Reported by the session that ran the
   comparison, unaudited here: Wang 2026 Case 3 at 1 km with xi 0.2
   (`runs/20261106_wang2026_case3/{c3a,c3b}`, `result/cycle*/global.dat`),
   90-97 % of interseismic steps reverse the sign of the step-to-step change
   in log10 Vmax; with per-node xi (PR #72, `par.xi_gamma = pi/4`) 0-1 %.
   The eq. 15 limit with k = gamma mu / h, a 0.005, b 0.01, L 0.05 m,
   sigma 25 MPa, 36 GPa is 0.097 at 1 km and 0.046 at 500 m, so xi 0.2 at
   1 km is about twice too large. Event times at 1 km moved by at most
   0.8 yr over five events; the defect is between events, not in the cycles.
   Open, owner's call: (a) `case.setup` computes the eq. 15 limit from
   material, h and friction and warns or stops when `par.xi` exceeds it;
   (b) per-node xi as the default with `par.xi` as the cap; (c) PR #72's
   open audit items -- `xi_gamma = NaN` silently off, the stop-13 check also
   runs at creeping nodes, stop 13 exits without `MPI_Abort`.
   Done when: PR #72 is merged or closed with (c) resolved, and (a) and (b)
   are each decided and, if yes, landed with a test.

Proposals for the owner (not decided): (1) `Slip_3` floor of 1e-15 m, `fltst` files only, measured by conductor 30 from the committed reference; the 1e-12 m/s and 1e-9 MPa floors were already decided. (2) The tolerance change to 1e-10 (disp and fault NetCDF at the repository tolerance after the 8-byte change) rides in PR 203 together with the regenerated reference, because a separate gold PR would leave master red.

Hand-off (conductor 30 to the next): (a) PR 203 (branch `c30-dp`, head b6c6e67, HOLD) waits on the gold-regeneration branch from the running implementer; run the fresh-clone fast tier and `-m e2e_ci` on the 1.24.8 build before folding it onto 203, merge only on `gh run watch --exit-status` at the exact head, remove the HOLD marker after; (b) then the slip-to-traction promotion stack off d78d822 (`par.solver = 4`), one src PR at a time, code audit before any tag, no new source file name until the owner names it; (c) remove the nested worktree `.claude/worktrees/c30-dp` after 203 lands and reap dead child worktrees only after copying any cited file out; (d) the remote branch `worktree-agent-a42a732ee98ad8c43` failed to delete, cause unread.

Hand-off addendum, conductor 30 (harness tool calls 88 of 100): the gold-regeneration implementer finished PARTIAL and is stopped. Branch `c30-gold`, head 05df35f, one commit on top of `c30-dp` b6c6e67, pushed; NOT gated by the conductor (its own report: fast tier 302 passed, `-m e2e_ci` 14 passed and 8 skipped for absent reference categories; unaudited).

- Replaced: 12 `fault.*.nc` in 4 sets (bp5 fast, bp5 dip90, stepover, stepover.con); `disp.00101.nc` added to dip90, stepover, stepover.con; data tree 94,260,435 to 100,994,278 bytes. Tolerance sites: `compare_netcdf` default 1e-10, `DD_DISP_RTOL` removed (disp and fault at `DD_RTOL`), a `disp` category and the bp1002 ft2/ruptarea cplot files added to the comparators, `"disp."` added to the contract read-glob.
- Not done, 5 sets: bp7 fast (`fault.00001.nc` fails the 1e-6 single-precision bar, max abs 3.8e-19 on a scale near 1e-18, relative 0.39: needs an owner decision, nothing replaced), bp8 (1.24.8 run done, 703 s; no 1.24.7 pair), bp7 full, bp1002 (including its 3 cplot files), bp5 kink, liu2020 (estimated over 4 h at 3 ranks, unmeasured). Estimates per set at 3-4 ranks: bp5 full about 50 min per binary, bp1002 about 45, bp5 kink about 40, each times 2 for the 1.24.7 pair.
- Open items for the next conductor: no 8-byte check exists on the reference files (a contract test is needed); `p1output.txt` fails a raw 1e-10 `compare_series` on stepover and passes the tier's own `_dd_profile_series` (a comparator choice to confirm); the bp5 and bp7 fast subdirectories hold no `disp` file (so `e2e_ci` skips that comparison there; one file per set was the instruction); worst difference at 1, 4 and 8 ranks against 1e-10, the restart run and write time are not yet recorded. Resume with `scratch/gold8/driver.py` and `scratch/gold8/check.py` in the implementer's worktree (uncommitted notes; copy out before reaping). PR 203 stays HOLD until the gold is complete or the owner accepts a partial set.

### Row 1, conductor 32: state recovered from the record (2026-11-18 22:20 CST; harness tool calls 10 of 100; budget read as wall clock to 2026-11-25 16:38 CST, row 1's order first, hand-off at about 80 calls)

Conductor 31 left no entry. Recovered from git, `ps` and `runs/` at 22:12 CST; every figure below is read this turn, nothing is gated or audited yet.

- Master `cc5137e`. Open PRs: #203 (`c30-dp`, head `4eca07d`, GitHub CLEAN) and the spike #87 (not for merge). Remote-only branches of conductor 31, none with a PR: `c31-p1` (`3b39a54` `par.solver = 4`, `00d6162` tests), `c31-p1fix` (head `09496fb`, version 1.24.10, stacked on `4eca07d`), `c31-goldA/B/D2/G`, `c30-gold8`, `worktree-agent-ad7e2b054bafedb8b`.
- #203 content (diff against master): the `dp` change `b6c6e67` (1.24.8); gold commits `05df35f` (4 sets), `f72f119` (new file `tests/contract/test_gold_netcdf_shape.py`, 82 lines, new file needing rule 1 agreement: asked below), `f03b8e2` (bp8 regenerated, bp7 full disp only), `2f822c0` (bp1002 with its 3 second-fault cplot files and last-step disp), `4eca07d` (bp7 fast, bp7 full, bp5 kink). 73 files changed, binary data only plus `src` (4 files), docs, tests. Not yet regenerated: liu2020.qdc.kink.600.
- liu2020 regeneration: the run in `runs/20261118_goldF/248/liu600` (3 ranks, nstep 10000, reference ends at `fault.10000.nc`) is NOT running: no `prterun` or solver process exists at 22:12, `cycle0.log` last written 22:10:35 at step 2821, `248.liu600.log` is 0 bytes, so no error text; cause unread (no `dmesg` access). Checkpoint files `fault/disp.{00001,01001,02001}.nc` are in its `scratch/`. A cycle cannot resume, so the set needs a fresh run (about 6 h at 2821 steps in 1 h 46 min).
- bp7 fast question: `4eca07d` DID replace the bp7 fast `fault.*.nc`, on the basis that each passes the tier comparator (`compare_netcdf`, rtol 1e-6 with an absolute floor of 1e-6) against the file it replaces. Conductor 30 had held it back because the variable of scale about 1e-18 differs relatively by 0.39 (absolute 3.8e-19). The owner's 16:45 condition says "agree at single-precision level, a looser bar of the order of the 1e-6 the tier used", worst difference recorded per file; whether the tier's absolute floor counts is the open reading. Being checked by me before #203 merges; if the file only passes through the floor it goes back to the old file and the question stays here.
- `par.solver = 5`: not a selectable value. `cb3c05e` adds only a test that `par.solver = 5` is refused with "par.solver must be 1, 2, 3 or 4". Nothing to remove from `src`.
- Operator-check floor (`570977e`: bound floored at 1e-14 of the operator's largest entry; motive: bp7 step 1 holds a traction of 8e-15, round-off, which a pure relative bound aborted) is a numerics choice for the owner, NOT agreed. Last gate at `cb3c05e`: fast 371 passed, `e2e_ci` 14 passed and 8 skipped, solver-4 tests 11 passed and 2 failed (bp7); `570977e` and `09496fb` never run. Not done: s/step of bp5 and bp1002 after the operator check, the planted-error check with the floor, the no-PETSc refusal test. No gate script holds `scratch/locks/heavy.lock` (a 0-byte file, no holder in `ps`); a process of another project (`mpirun -np 10 eqdyna`, PID 8034) loads the host.
- Worktrees with cited files, not to be reaped before copy-out: `agent-a2e8a419ac82b9331` (`runs/20261118_p1fix/` logs), `agent-a1a440daa6089e9de` (locked; goldF), `agent-af6f00fefea46947a` (`scratch/gold8/driver.py`, `check.py`).
- Queued for the owner's daily brief (09:37 CST): (1) the bp7 fast file reading above; (2) `Slip_3` floor of 1e-15 m (already on this row as proposal 1 of conductor 30); (3) the operator-check floor of 1e-14 of the operator's largest entry, with the alternatives (pure relative, which aborts bp7 step 1; a floor tied to the traction scale) once characterized; (4) the gold lands in the same PR as #203; (5) rule 1: `tests/contract/test_gold_netcdf_shape.py` is a new file, named by conductor 31 without agreement; it stays unless the owner objects, and the check can instead move into an existing contract file.

Order of work (invoking session): this entry; finish and gate the gold and land #203 with master never red; gate, audit and land the s2t stack as `par.solver = 4`; floors in the comparator; s2t optimization and the rest of the board in `prio` order. Owner-held and untouched: ICNTL(14) for native solver 1, a `vista` machine target, any default change, any new file or name beyond those agreed, other absent bp1002 reference files.

### Row 1, conductor 32: #203 audit, gold measured, liu2020 relaunched (2026-11-18 22:35 CST; harness tool calls 43 of 100)

- Correction to the entry above: `tests/contract/test_gold_netcdf_shape.py` is an existing file that #203 extends by 82 lines (`git diff --diff-filter=A` lists no new file outside `data/`); rule 1 is not engaged and the fifth item of the owner-brief list is struck.
- bp7 fast and the reading of the owner's 16:45 condition, MEASURED by me this turn (all 53 replaced `fault.*.nc` of #203, each old (master) against new (`4eca07d`) with the tier's own `compare_netcdf`): 0 fail at rtol 1e-6 with floor 1e-6; 31 fail with the floor at 0 (bp5 kink 22 of 22 after step 600, bp8, stepover, stepover.con). The only variables above 1e-6 of their own scale are zero by physics: `shear_dip` (own scale 9.7e-19, group scale 1.6e7), `slipd` (7.1e-38), `vzm` and `vzs` (3.5e-35), and `slipn` (own scale 5.6e-12, abs 2.3e-15, 4.0e-4 of own scale in the bp5 kink set); largest difference over the group scale of its field: 1.0e-13. Every signal variable agrees within 1e-6 of its own scale; the worst per-file figure (difference over the file's own maximum) is 5.7e-8 for bp1002, bp5 fast, dip90, stepover and bp7 `fault.00101`, 4.6e-5 for bp7 full from step 1001, 2.4e-7 for bp8, and 3.9e-1 (bp7 `fault.00001`, the zero-by-physics fields above) and 3.9e-4 (kink `slipn`). Reading kept: the replacement meets the condition as the tier defines its 1e-6 bar (group scale, absolute floor); under an own-scale reading the five zero-by-physics variables would not, and that is the question for the owner brief. The files stay replaced.
- Tree growth, #203 head: `data/` 94,260,007 to 112,767,711 bytes (+18.5 MB, `git ls-tree -r -l`). A liu2020 last-step `disp` file is about 6.3 MB, above rule 22's 5 MB per-file cap: owner question (exempt the file, or keep liu2020 without its `disp`); liu2020 therefore waits and does not block #203.
- Audit of #203 (read-only, role auditor, 12 tool calls): no BLOCKER, verdict "merge after fixes". MAJOR 1: `data/test.bp5.qdc.2000/cycle0` (full tier) still holds 4-byte gold compared at 1e-10, so the next `-m e2e` fails it; not on the width list. MAJOR 2: the CI case `test.bp5.qdc.2000/cycle0-step101-fast` has no `disp.*.nc`, so the new `disp` category skips there (rule 3a); owner decision fd153ba wants one per set. MAJOR 3: restart files now hold unrounded values, so cycle 2 and later start from a different state than under 1.24.7 (cycle 1 unaffected); the owner was told fault files are read back on restart when answering "both follow dp", and the doc line was missing. MINOR: bp8 absent from `REGENERATED_8BYTE`; `FAULT_SNAPSHOTS` globs `data/*/fault.*.nc` and finds only bp8 (older than #203); the string "1.24.8" is also the version of an uncommitted spike binary. Fixes in preparation on branch `c30-dp`: restart paragraph in `docs/user/outputs.md`, bp8 on the list, bp5 fast `disp.00101.nc`, bp5 full regeneration (50 min per binary at 4 ranks, queued behind `scratch/locks/heavy.lock`). Multi-cycle difference not measured; queued for the owner brief as a consequence of the dp decision.
- Runs (detached with `setsid`; conductor 31's run died with its shell at 22:10, likely cause, unproven): liu2020 relaunched 22:17:33 (`runs/20261118_goldF/248/liu600`, run.sh PID 719210, prterun 719237, 3 ranks, binary copy `scratch/c32/bin/eqquasi-1.24.8` md5 dc714b19, same as the one conductor 31 used; dead-run output kept in `died_at_step2821_22h10/`; about 6 h, finish about 04:30); rank study at 4, 1 and 8 ranks, `runs/20261118_c32_rank/` (driver `scratch/c32/rank_study.py`, PID in `study.pid`); text pair 1.24.7 against 1.24.8 at 4 ranks, `runs/20261118_c32_pair/` (`scratch/c32/pair.py`). Gate trees: `.claude/worktrees/c32-gate` (head `4eca07d`, with `bin/eqquasi-1.24.8`), `c32-m247` (master, `bin/eqquasi-1.24.7`), `c32-fix` (the PR fix edits).

### Row 1, conductor 33: board split landed, budget reading, owner-brief block (2026-11-18 23:10 CST; harness tool calls 26 of 100; unaudited)

- Budget read as: the owner's "/autopilot clear the board. a week budget. go" (16:38 CST) covers every open row in `prio` order until 2026-11-25 16:38 CST with this row's order first (#203 with its gold as one PR, the slip-to-traction stack as `par.solver = 4`, comparator floors, then optimization and the rest). Owner-held and untouched: ICNTL(14) for native solver 1, a `vista` machine target, any default change, any new file or name beyond those agreed, the operator-check floor until answered. Hand-off at about 80 calls.
- Landed: #210 (`65b3e5c`, CI run 44460871482 success on head `8359e09`): the row 1 measurement history (2026-11-14 to conductor 30) moved verbatim into `docs/dev/board_history.md`, section "Row 1 measurement history, 2026-11-14 to 2026-11-18 (archived 2026-11-18)"; board 345,997 to 52,349 bytes. Proof run by me on the pushed head: every Owner paragraph and open-item paragraph byte-identical, the archive's old text an untouched prefix, 0 non-blank old lines lost. Existing file, no new name (rule 1 not engaged). Pre-existing in the moved text, not edited: a repository URL in the archive that carries a personal handle (`board_history.md:979`); the owner decides whether past logs are rewritten (item kept on row 1).
- #203: head `9486225` (restart paragraph, bp8 on the 8-byte list, bp5 fast `disp.00101.nc`, pushed by me, fast-forward). Still to do before merge: bp5 full regeneration (1.24.8 run prterun 735604 started 22:25, about 50 min), then fast tier and `-m e2e_ci` at the merge head from a fresh clone, then `gh run watch --exit-status` at the exact head. liu2020: run.sh 719210 / prterun 719237, started 22:17, expected about 05:00; it does not block #203.

Owner brief, 09:37 CST (questions wait for it; every item has its options and my recommendation; nothing below is decided):
1. bp7 fast and the other replaced `fault.*.nc` (53 files, measured by conductor 32: 0 fail at rtol 1e-6 with the tier's absolute floor of 1e-6, 31 fail with the floor at 0; the failing variables are zero by physics, scale 1e-18 to 1e-38). A) keep the replacement, the tier's own bar being the single-precision bar (recommended: the signal variables agree within 1e-6 of their own scale and the physics is unchanged). B) require agreement at 1e-6 of each variable's own scale: five zero-by-physics variables fail, so bp7 fast and the bp5 kink would keep their old files and fail the 1e-10 tier bar once the output is 8 bytes.
2. `Slip_3` floor of 1e-15 m, `fltst` files only (measured by conductor 30 from the committed reference). A) adopt (recommended; the 1e-12 m/s and 1e-9 MPa floors were already decided). B) no floor: the series fails on round-off-level values.
3. Operator-check floor of the slip-to-traction check (`570977e`: bound floored at 1e-14 of the operator's largest entry; bp7 step 1 holds a traction of 8e-15, round-off). A) adopt the 1e-14 floor (recommended only together with a test that plants a wrong operator and shows the check still aborts it; that test is audit finding M2 and is not yet written). B) pure relative bound: aborts bp7 step 1. C) a floor tied to the traction scale: not yet characterized. Not merged before the answer.
4. The gold lands in the same PR as #203. A) one PR (recommended: a separate gold PR leaves master red for the tolerance change). B) two PRs, master red between them.
5. liu2020 last-step `disp` file, about 6.3 MB against rule 22's 5 MB per-file cap. A) exempt this one file, recorded by the rules owner (recommended: the owner's decision wants one `disp` per set). B) keep liu2020 without its `disp`.
6. Restart files now hold unrounded values (dp), so cycle 2 and later differ from 1.24.7 at the level of the 4-byte rounding; cycle 1 is unchanged; the docs say so. Not measured. A) accept, after I measure the cycle-2 difference on bp5 fast with both binaries and log it (recommended). B) keep 4-byte restart rounding: a second code change.

### Row 1, conductor 33: #203 landed, HAND-OFF (2026-11-18 23:42 CST; harness tool calls about 50 of 100, a context summary at 23:40 voids the count: no new work, no further merge; unaudited)

- Landed: #210 (`65b3e5c`, board split) and #211 (`231f8f1`, owner-brief block, six items with options and recommendations) by board PR; #203 (`6934a28`, NetCDF reals follow dp, v1.24.8, with the gold of every set but liu2020; its bp5 full regeneration added at head `4a170d4`: faults at most 2.6e-11 of scale from the file replaced, 18 text files at most 1.9e-14). Gate on `4a170d4` from a fresh bare clone under `env -i` (`runs/20261118_c33_gate203/`): fast tier 384 passed with the binary built (383 passed, 1 failed before the build, the binary-on-disk test), `-m e2e_ci` 15 passed, 7 skipped (data-absence skips: no offfault, perfault, profile files, no `summary.json`, no section 4.3 profiles; none new). GitHub CI run 44463991495 success on that head; master CI run 44465226691 success on `6934a28`. The gate and the reference-data record are the comment on #203; the PR body's HOLD text was replaced by that record.
- NOT done on #203's follow-through: liu2020 regeneration still running (run.sh 719210 / prterun 719237 started 22:17, 1 h 22 min at 23:40, expected about 05:00; its `disp` of about 6.3 MB exceeds rule 22's 5 MB cap: owner brief item 5); the multi-cycle restart difference is not measured (brief item 6). The 1.24.7 pair driver (`scratch/c32/pair.py`, PID 733568) and the rank study (PID 732912) still run; their outputs are unread by me.
- Slip-to-traction stack: branch `c33-s2t-fix`, head `9903acd` (pushed; on top of `c31-p1fix` `09496fb`; version 1.24.11), UNAUDITED, UNMERGED, NOT rebased onto `6934a28`. The implementer's report (not checked by me): B1 fixed (the test uses the solver's bound `max(1e-8*field_max, 1e-14*operator_max)`, marked "PROPOSAL, owner decision pending"), M1 fixed (NaN flag per value, reduced across ranks), the abort names both bounds, a planted-row test of the assert alone; M2 only partly done (no run-level planted-operator test: it needs a solver hook, a new name, rule 1); minors (b) allocation message, (c) non-MUMPS override refusal test, (d) no-PETSc refusal test NOT done; 1-, 3- and empty-rank layouts open. Its gate (own worktree, 4 ranks): fast tier 371 passed; solver-4 tests passed for bp1002 zone, bp5 and bp7.a.10, the unsupported-combination refusals (krylov3, krylov4, solver5) and the planted-row test; the 8- and 12-rank cases (bp7.8, bp1002.12, bp8 gs.8 and pw.8) were not run (over the 4-rank cap). Next: rebase onto master, re-run those cases at the allowed rank count on a quiet host, a code audit of the final head, no merge before the owner answers brief item 3 (operator-check floor).
- Next for a fresh conductor, in order: (1) read the owner's answers to the brief and apply them; (2) liu2020 gold when its run ends (copy out, compare against the old file with the tier comparator, apply per the owner's answer on item 5), one gated PR; (3) the stack above as `par.solver = 4`; (4) comparator floors (`Slip_3`, brief item 2); (5) s2t optimization and the rest of the board in `prio` order. Owner-held: ICNTL(14) for native solver 1, a `vista` machine target, any default change, any new file or name (the planted-operator hook included).
- Worktrees and branches I made: `.claude/worktrees/c33-board` (branch `c33-handoff`, this entry), `agent-a94de11514ecea170` (finished, reapable after a listing for gitignored files), `agent-a52c2bbfc0ba82b73` (finished, branch `c33-s2t-fix`; keep until the stack lands); `.claude/worktrees/c32-fix` holds gate copies under `scratch/gold8/bp5full` (ignored, copy-out not needed: the run dir `runs/20261118_c32_pair/v1248_np4/bp5full` stays). Not reaped: all others (56 at the start; none shown unreferenced this run).

### Row 1, conductor 34: #212 landed, dispatches, checkpoint (2026-11-19 00:06 CST; harness tool calls 17 of 100; unaudited)

- Budget read as: row 1's order first, then every open row in `prio` order, to 2026-11-25 16:38 CST; items 1-6 of the owner brief wait for the owner's 09:37 CST answer and nothing past them merges. Hand-off at about 80 calls.
- Landed: #212 (`3aae60b`; CI run 44466783010 conclusion success on head `a90f2a6`, state MERGED read after the merge). The branch-delete error after merges is the explicit `git push origin --delete` finding the branch already gone (`error: unable to delete 'c33-handoff': remote ref does not exist`): the repository deletes head branches on merge, so the explicit delete is redundant. `c33-brief` is not on origin (`git ls-remote --heads origin 'c33*'` lists `c33-s2t-fix` only).
- Dispatched, both cap-bound, worktree-isolated, heavy runs behind `scratch/locks/heavy.lock`: (1) the slip-to-traction stack, branch `c34-s2t` from `6934a28` (8 commits of `c33-s2t-fix` cherry-picked, one bump to 1.24.9, minors (b)(c)(d) without a new name, gates at 1/3/4/8/12 ranks, setup, s/step and memory per case, draft PR marked HOLD; cap 60 calls); the planted-operator hook (M2) is NOT built: its name goes in the report for the owner (rule 1). (2) the cycle-2 restart difference between 1.24.7 and 1.24.8 on bp5 fast (measurement only, brief item 6; cap 40 calls).
- Still running, `ps` this turn: prterun 719237 (liu2020 gold, 3 ranks, started 22:17), prterun 945976 (rank study bp8, 4 ranks, started 23:43), python 732912 (rank study driver), python 733568 (1.24.7 pair driver). Outputs unread; to be harvested onto this row when they end.
- Reaped (clean, merged, no ignored output beyond `__pycache__`): worktrees `c25-board2`, `c25-board3`, `c26-board`, `c32-board`, nested `owner-relays-2/.claude/worktrees/c30-dp`; local branches `board-c25-2`, `board-c25-3`, `board-c26-5`, `board-c32-2`, `c30-dp`. 56 worktrees remain; every other one is either a live agent, holds ignored output not yet listed, or is referenced.

### Row 1, conductor 34: audit of the slip-to-traction stack, fixes, gates queued (2026-11-19 00:52 CST; harness tool calls 48 of 100; unaudited)

- Stack PR #214 (draft, HOLD), branch `c34-s2t`: 8 commits of `c33-s2t-fix` cherry-picked onto `6934a28` with no conflict, one version 1.24.9 (head `3daa5ec`, then `71a09ed`). Diff against master: 7 files, no `data/`, no reverted line (implementer's report; I have not re-read the diff).
- Read-only code audit of `3daa5ec` (role auditor, 14 tool calls, unaudited by me beyond reading it): "merge after fixes". Rank layouts 1, 3, zero-row ranks: NONE (code read, every collective entered by all ranks, zero-row guards at `solveTimeLoopPETScDD.f90:378,552,570`); no integer overflow reachable; solvers 1-3 untouched outside solver-4 branches. MAJOR F1 (`:500,503`): a NaN in the pair displacement or in S*g on a non-check step is not flagged. MINOR: non-MUMPS factor never refused although the docs say so; `nonmumps4` asserts only a nonzero exit; troubleshooting wording of the floor; the planted-operator test checks the Python copy of the bound, not the Fortran abort (M2, needs the hook named below).
- Fixes pushed as `71a09ed` (implementer report, not yet gated by me): F1 finite-check of both gathered arrays on every step (rank 0 prints, all ranks abort), F2 `PCFactorGetMatSolverType` check with message naming the type and `nonmumps4` asserting that text, F3 wording, F4 docstrings. Fast tier 384 passed, no-PETSc link rc 0 (that copy does not compile the PETSc loop file, so it does not exercise the new Fortran). NOT MEASURED: s/step before and after F1.
- Owner-held, unchanged: the operator-check floor (brief item 3); the planted-operator hook, proposed name `-eq_s2t_plant_row <global row>` (a PETSc option read next to `-eq_s2t_build`, scales one row of S by 1+1e-6 after the build, before the check): a new name under rule 1, NOT built.
- NOT DONE and why: the solver-4 gate at 1, 3 and an empty-rank layout and the three full-cycle solver-4 cases: the tests hard-code 4, 8 and 12 ranks (`DD_LAYOUT_CASES`), 1/3/empty-rank/full-cycle cases need a driver that does not exist yet; they wait behind the lock.
- Queue on `scratch/locks/heavy.lock` (heavy runs one at a time; ps this turn): 1.24.7 pair `bp5full` (python 733568), rank study np 1 and 8 (python 732912), the cycle-2 restart measurement (`scratch/c34_restart/run.sh` PID 979071; runs `runs/20261119_c34_restart/{v1247,v1248}`, 3 cycles each, 4 ranks), light solver-4 check of `71a09ed` (`scratch/c34_s2t/fix_light.sh` PID 1108034), and the gate driver for `71a09ed` (`scratch/c34_s2t/gate_71a09ed.sh`, PID/PGID 1112415; status `runs/20261119_c34_s2t/status_71a09ed.txt`: fast tier, solver-4 e2e_fast tests, rest of e2e_fast, each step takes the lock for its own run, caps 5400 and 7200 s). The earlier driver for `3daa5ec` (PID 1024005) was stopped by me by PID after the audit, status kept as `status_3daa5ec_fastonly.txt`. Waiters are not FIFO (rule 15b).
- Stale board text found, NOT edited (re-scoping is the rules owner's): the "OWNER ITEM, holds the merge" paragraph and "OWNER ITEMS OPEN (2)" about the refactor file name describe step 1 as waiting; the move landed as `c86bafc` (#202, v1.24.7) and the owner agreed the name (answers quoted at lines 244-252). A dispatch of mine on that stale text cost one agent 9 tool calls and delivered nothing; the check `git log origin/master --grep onFaultHelpers` before dispatching would have caught it. Worktree `agent-a42a732ee98ad8c43` and its local branch are reaped (clean, `git cherry` empty).
- Comparator floors: `Slip_rate_3` 1e-12 m/s and `Shear_stress_3` 1e-9 MPa are agreed, but the owner's text leaves open how a floor applies to the log10 columns and the `Slip_3` value (brief item 2); no floor is coded in `tests/e2e/cases.py` `compare_series` (line 655) yet. Held for the 09:37 answer so the change lands once.
- Rank study (`runs/20261118_c32_rank/study.log`, read this turn): np 4 finished for bp5 fast, bp7 fast, dip90, stepover, con and bp8, "0 bad" in each (rows 13, 13, 13, 13, 13 and 24); np 1 and np 8 still queued. The 1.24.7 pair: bp5 fast, bp7 fast, dip90, stepover and con done, `bp5full` was running; no comparison read yet.
- Still owned: prterun 719237 (liu2020 gold; expected about 05:00-05:15), the drivers above. Specialists live (`ListAgents` not available to me; completion notices received as "waiting on its own background work"): the restart-measurement agent `ad5d317dcab37a765` (idle, TaskStop refused for me) and the fixer `a84adddbe0e32be5b` (idle). Worktrees held: `agent-aad3cc924cc4cdfe4` (old gate dir; the fix driver uses `agent-a84adddbe0e32be5b`), do not reap.

### Row 1, conductor 34: restart measurement (brief item 6), pair and rank evidence, fix round 2 (2026-11-19 01:25 CST; harness tool calls 70 of 100; unaudited)

- Brief item 6 MEASURED (bp5 fast, 4 ranks, 3 cycles each, one run per binary, both versions read from each cycle's `runInfo.json`: 1.24.7 and 1.24.8; case inputs and `user_defined_params.py` identical; repository comparators at rtol 1e-9, floor 1e-9, group scale = largest magnitude in the 1.24.7 file; `scratch/c34_restart/cmp.out` written 01:07:32; 1.24.7 `result/cycle1/runInfo.json` 00:46:42, 1.24.8 01:06:05; I re-read three lines of `cmp.out` myself and they match the report). Cycle 1: text series pass (worst 1.7e-14 on `global.dat`), NetCDF differ by 4.2e-8 of scale (float32 storage only). Cycle 2: `global.dat` 4.8e-7 (peak slip rate 0.046756336 vs 0.04675634), `fltst` files 1e-7 to 1.2e-5, fault and disp NetCDF 1.2e-6 of scale, slip rate 7.1e-7. Cycle 3: `global.dat` 2.7e-6, NetCDF 1.3e-6, slip rate 9.0e-7. Both finish every cycle at 101 steps; the fast case has no event, so event count and event times are NOT compared (a longer case is needed); no same-binary repeat, so no reduction-order noise floor. Wall per step 0.651 s (1.24.7) against 0.674 to 0.694 s (1.24.8), 3 to 7 %, host shared with the liu2020 run: suggestive only. Reading for the owner (option A of item 6): the drift is the 4-byte restart rounding, about 1e-6 of scale after two restarts, growing 1.0 to 1.9x from cycle 2 to cycle 3.
- 1.24.7 against 1.24.8, one cycle, same cases and ranks (4), tier comparators, MEASURED by me this turn with `scratch/c34_conductor/pair_cmp.py` on `runs/20261118_c32_pair` and `runs/20261118_c32_rank`: bp5 fast, bp7 fast, dip90: text series 0 of 26 fail at 1e-9, NetCDF 0 of 6 fail at 1e-6; bp5 full (4483 steps both): 0 of 26 text, 0 of 14 NetCDF; stepover and con: 2 of 16 text fail, both `p1output.txt` (216 of 1080 entries outside 1e-9, worst in column 5 at about 2e-10, the near-zero profile column: the known raw-comparator failure that the tier's own `_dd_profile_series` passes; not re-checked here). Physics unchanged between the binaries at the level the tier measures; this is the evidence for the owner's "physics shouldn't change" condition.
- Rank study against the regenerated reference (`runs/20261118_c32_rank/*_np*.json`, tier comparators at 1e-10 for disp and fault files): np 4: bp5 fast, bp7 fast, dip90, stepover, con and bp8, 0 bad each (rows 24, 25, 25, 13, 13, 24); np 1: bp5 fast and bp7 fast 0 bad (24, 25 rows), the other np 1 cases and np 8 still queued.
- Slip-to-traction stack, PR #214, head `65a8435` (pushed, fast-forward): `1184397` (the `nonmumps4` test asserts the refusal PETSc itself gives at `KSPSetUp`, because the `ddS2tBuild` factor check is unreachable with `-pc_factor_mat_solver_type petsc`) and `65a8435` (new case `superlu4`: `-pc_factor_mat_solver_type superlu_dist` with Cholesky reaches the new check and prints `EQ_S2T: the factor solver type is 'superlu_dist', not mumps`; `-pc_type lu` is stopped earlier by the existing PC-type guard). Fixer's report, not re-run by me: `superlu4` and `nonmumps4` 2 passed, fast tier rc 0, F1 s/step before and after NOT MEASURED (the baseline never got the lock). Gate driver `scratch/c34_s2t/gate_65a8435.sh` PID/PGID 1167975 (ps this turn), status `runs/20261119_c34_s2t/status_65a8435.txt` (fast rc=0 at 01:12:04; the s2t e2e_fast and rest of e2e_fast steps wait for the lock). A dispatch for the 1-rank, 3-rank, empty-rank and three full-cycle pairs (detached driver, `scratch/c34_layout/driver.pid`, results `runs/20261119_c34_layout/`) is running.
- Process notes: the fixer used `pkill -f` once and it killed only its own launching call (it reports no other job touched; my processes 719237 and 732912 and the new gate driver are alive by `ps` at 01:14); it also ran two small tests under a lock named `light.lock` instead of the heavy lock (25 s). 1.24.7 pair driver (python 733568) finished: its six cases are in `runs/20261118_c32_pair/v1247_np4/`. Agent costs (tool_uses from completion notices): first stack agent 42 calls (staged, gates NOT RUN: lock held), audit 14, refactor dispatch 9 (nothing delivered: the move was already #202, stale board text), restart measurement 32 over about 74 min of lock wait, fixer 74 against a cap of 40 plus 12 (over the cap).
- Unblock events: heavy.lock queue drain (drivers above); the owner's 09:37 answer (items 1-6, comparator floors, planted-operator hook name); prterun 719237 end (about 05:00 to 05:15, liu2020 gold).

### Row 1, conductor 34: HAND-OFF (2026-11-19 01:30 CST; harness tool calls 76 of 100; no live agent; unaudited)

- Landed this run: #212 (`3aae60b`, run 44466783010 success), #213, #215, #216 (board entries, CI success each, state MERGED read after). No code PR landed: #214 is HOLD (owner item 3).
- Still owned by me, no agent behind them (`ps` at 01:30): gate driver for `65a8435` PID/PGID 1167975 (`runs/20261119_c34_s2t/status_65a8435.txt`; steps: s2t e2e_fast then rest of e2e_fast, caps 5400 and 7200 s); layout and full-cycle driver PID/PGID 1184934 (`scratch/c34_layout/driver.py`, log `driver.log`, results `runs/20261119_c34_layout/<case>_np<k>.json`, `DONE` file at the end; solver 3 against solver 4 on the 1.24.9 binary copy; 1, 3 ranks, an empty-rank layout derived after A, bp5 and bp7 full cycles at 4 ranks, bp1002 zone is only 101 steps in its compset so it is not a full cycle; estimated 5 to 6 h of lock time; the baseline is solver 3 where the repo test compares against solver 2); rank study python 732912 (np 1 and 8 remaining, `runs/20261118_c32_rank/study.log`); liu2020 gold prterun 719237 (expected about 05:00 to 05:15; `runs/20261118_goldF/248/liu600`; its `disp` of about 6.3 MB is not committed before the owner's item 5). Kill by PID from the pidfiles after `ps -o pid,args -p`, never `pkill -f`.
- Next for a fresh conductor, in order: (1) the owner's 09:37 answers (items 1-6, `Slip_3` and log10 floors, the planted-operator hook name `-eq_s2t_plant_row`); (2) harvest the drivers above onto this row (the gate status files, `runs/20261119_c34_layout/*.json`, the rank study np 1 and 8, the liu2020 run), reading each result file with its mtime; (3) if the owner accepts the floor: update PR #214 with the evidence table and an audit of the exact head (the earlier read-only audit covered `3daa5ec`; `71a09ed`, `1184397`, `65a8435` are fixer-reported only), merge on `gh run watch --exit-status` at the exact head; (4) comparator floors in `tests/e2e/cases.py` `compare_series` once, after the answer; (5) liu2020 gold as one gated PR; (6) s2t optimization (operator build 2131 s at the last harvest) after #214; (7) the rest of the board in `prio` order (rows 3, 4, 6 are owner-decision rows; row 1 steps a to d need the host to ourselves).
- Worktrees: kept `agent-a84adddbe0e32be5b` (branch `c34-s2t-fix`, the gate drivers run from it; do not reap until both drivers end), `agent-aad3cc924cc4cdfe4` (stale head `3daa5ec`, reapable after the drivers), `agent-ad5d317dcab37a765` (restart measurement; `scratch/c34_restart/cmp.py` imports the comparators from it, keep until item 6 is answered), `agent-ae69cb9f40076717f` (layout dispatch, its driver may import from it; keep). Repository auto-deletes merged head branches; `c34-s2t` stays open for #214.

### Row 1, conductor 35: #217 landed, audit of `65a8435`, brief refresh, optimization plan (2026-11-19 01:40 CST; harness tool calls 12 of 100; unaudited)

- Budget read as: the owner's "/autopilot clear the board. a week budget. go" (16:38 CST) covers every open row in `prio` order until 2026-11-25 16:38 CST, row 1's order first. Items 1-6 of the owner brief, the `Slip_3` floor and its log10 handling, the planted-operator hook name, ICNTL(14) for solver 1, a `vista` machine target, any default change and any new file or name stay held; nothing past them merges. Hand-off at about 80 calls.
- Landed: #217 (`a010b22`, state MERGED read from `gh pr view` at head `6057f83`, merged by the previous conductor's watcher).
- Code audit of the exact head `65a8435` of #214 (read-only, role auditor, 7 tool calls, nothing built or run; unaudited by me beyond reading it): "merge after fixes", 0 BLOCKER, 0 MAJOR. Closed: the finite check of both gathered arrays on every step (`src/solveTimeLoopPETScDD.f90:1549-1566`, `ddS2tApply` runs every step), the non-MUMPS factor refusal (`:1171-1177`, message matches the `superlu4` test at `tests/e2e/test_benchmarks.py:987`), `nonmumps4` asserting PETSc's own text (`:980-982`; the text is PETSc's and was not confirmed against PETSc). No new collective (both aborts test arrays identical on every rank after the gather; two local scans per step); solvers 1-3 untouched outside solver-4 and shared-refusal hunks; version 1.24.9 appears once. Remaining MINOR, go to the follow-up list, not a fix round: (m1) the station displacement `s2tU` (`:1577-1579`) has no finite check on non-full steps; (m2) the first abort line still says only "beyond 1e-8 of the field maximum" and understates the floored bound; (m3) the planted-operator test checks the Python copy of the bound (`test_benchmarks.py:847-862`), not the Fortran abort (needs the hook, item 7 below); (m4) no test plants a NaN, so the new finite checks have no run-level test. Not checked by the auditor: compile, any run, 1-rank and 3-rank layouts, PETSc behaviour of `superlu_dist` under Cholesky. #214 stays draft and HOLD (item 3).
- Row 6 (read from source and `gh pr view` this turn): PR #72 state MERGED; its open audit items (c) are fixed on master: `xi_gamma` NaN is rejected (`src/read_input.f90:232`, `.not. (xi_gamma >= 0 .and. xi_gamma <= 1e6)`), the per-node xi applies at frictional nodes only (`src/faulting.f90:461-462`), stop 13 is `MPI_Abort(MPI_COMM_WORLD, 13, ...)` (`src/faulting.f90:470-473`). Not checked: a test pinning each. Request to the rules owner: row 6's "Done when" now waits only on the owner's decisions (a) and (b); re-scoping the row is theirs.
- Running, `ps` at 01:25 (no agent behind any): gate driver PID 1167975 (`runs/20261119_c34_s2t/status_65a8435.txt` reads only "fast rc=0 01:12:04"), layout and full-cycle driver PID 1184934 (`runs/20261119_c34_layout/` holds `test.bp5.qdc.2000_np1_s3` only), rank study PID 732912 (np 1: bp5 fast, bp7 fast, dip90 0 bad each; `runs/20261118_c32_rank/study.log`), liu2020 gold prterun 719237 (expected about 05:00-05:15).
- Owner brief, refreshed at 01:40 (options and recommendation per item; nothing decided; the 09:37 brief is the place to put them):
  - Items 1, 2, 4, 5: unchanged from the entry above (1 A keep; 2 A adopt 1e-15 m; 4 A one PR, already done for #203 as landed; 5 A exempt the one liu2020 `disp`).
  - Item 3, operator-check floor: the audit of `65a8435` says the code implements what its message says (`max(1e-8*field_max, 1e-14*operator_max)`); it does not bear on the value. Recommendation unchanged: A adopt 1e-14 only together with a run-level planted-operator test (item 7); the layout and full-cycle sweep (np 1, 3, an empty-rank layout, bp5 and bp7 full cycles) is not finished, so the floor has no new evidence yet.
  - Item 6, restart rounding: now measured (conductor 34 entry: cycle 2 about 1e-6 of scale, cycle 3 about 1.3e-6; event times not compared on the fast case). Recommendation A accept, with a longer case run for event times if the owner wants them; B (keep 4-byte restart rounding) is a second code change.
  - Item 7 (new), planted-operator hook: A) `-eq_s2t_plant_row <global row>` as a PETSc option next to `-eq_s2t_build` (recommended: it is the only way a run-level test shows the check aborts a wrong operator); B) another name; C) no hook, keep the Python-copy test only (then item 3 A loses its precondition and the floor is unprotected at run level).
  - Item 8 (new), the floor on log10 columns (`Slip_3`, `Slip_rate_3`, `Shear_stress_3` in the `fltst` files): A) compare the un-logged values and apply the absolute floor to them, skipping a log10 entry whose value is below the floor (recommended; no unit mixing); B) apply no floor to the log10 columns (round-off values fail); C) a floor in log10 units. Not measured on the reference; no code until answered, so the comparator changes once.
  - Item 9 (new), ICNTL(14) for native solver 1: A) leave at the default until a native-solver failure is shown (recommended: the Vista failures were on the PETSc paths and no native-solver failure is on this board); B) default 100 as on the PETSc paths.
- Measurement plan, s2t operator-build optimization (proposal; nothing run; after promotion, behaviour-preserving, each step before and after on the same case, ranks and host, with comparator against exact at 1e-9, setup, s/step, peak memory per rank, machine and library versions). Baseline from the Vista harvest `runs/20261118_vista/p13harvest/S2T_build/case/run.log`: build 2130.93 s, 61206 columns in 1913 batches of 32 at 144 ranks (about 1.1 s per batch), 1,231,276 nodes. The build (`ddS2tBuild`, `:1152`) calls `MatMatTransposeSolve` on the MUMPS factor per batch of 32 (`nbat0`) with the sparse right-hand sides gathered to rank 0 each batch. Step 0 (no source change): bp5 and bp1002 at 4-8 ranks with PETSc `-log_view`, splitting gather and assembly, solve, and column copy; this says which term to attack. Step 1 (no source change): `-mat_mumps_icntl_27` (blocking size of the multiple-RHS solve) and `-mat_mumps_icntl_20` through `PETSC_OPTIONS`, one run each. Step 2 (spike branch, new option name owner-held): batch size as a runtime value instead of 32. Step 3 (design question, not a plan): exploit the symmetry of the operator to build half the columns. Runs wait for `scratch/locks/heavy.lock`; step 0 and 1 are the first to run after the sweep drains.
- Unblock events: heavy.lock drain (gate, layout, rank-study drivers); the owner's 09:37 answers; prterun 719237 end (liu2020).

### Row 1, conductor 35: solver-4 e2e_fast failure on `65a8435` classified (2026-11-19 02:00 CST; harness tool calls 22 of 100; unaudited)

- Fact (`runs/20261119_c34_s2t/status_65a8435.txt`: `s2t_e2e_fast rc=1 01:37:46`; `s2t_e2e_fast_65a8435.log`: 1 failed, 16 passed, 844 s): the failing test is `test_slip_to_traction_matches_petsc_loop_z_and_corner_blocks[bp7.8]` (solver 4 at 8 ranks against solver 2 at 4, 31 steps); assert text "srfst_strk000st400dp300.txt: 25 of 217 entries outside rtol=2e-06; worst at row 30, column 4: 7.003613e-17 vs -5.839928e-17". The log of head `1184397` also reads 1 failed, 16 passed (`s2t_e2e_fast_1184397.log`), the head before `65a8435`; not read further. bp7.8 had not run before (over the 4-rank cap), so no earlier head passed it.
- Classification, a reading: test (comparator) defect, not a code defect and not a lock artefact. Evidence: (1) the solver-2 output of the failing run (`scratch/e2e.dd_s2.s4.bp7.8/result/cycle0/srfst_strk000st400dp300.txt` in the gate worktree) is a 31 x 7 file whose columns 3 and 5 peak at 5.8e-17 and 5.5e-17 and column 4 at 6.8e-22, round-off of an out-of-plane component, so each of those columns is judged against its own round-off maximum; (2) the compared values change sign (7.0e-17 against -5.8e-17); (3) `_dd_differences` (`tests/e2e/test_benchmarks.py:587`) sends `srfst_` through `compare_series` with only `rtol`, while `fltst_` goes through `_dd_station_series`, which has the dip-column floor against the vector's largest component; `srfst_` has none; (4) the run completed, the same failure appears at two heads, and the result is deterministic in kind, so it is not contention. NOT falsified: the control, solver 3 against solver 2 on bp7.8 with the same comparator, is in the `rest_e2e_fast` step (still running) and in the layout driver (solver 3 against 4); if solver 3 fails the same way the defect is older than #214.
- Consequence: #214's HOLD has a second reason. The fix is a comparator change (a floor for `srfst_` out-of-plane columns, like the `fltst_` one), which is a tolerance decision and is owner-held with items 2 and 8: the owner's 13:57 "1A" floor covers `fltst` columns, not `srfst`. Owner item 10 (new): A) extend the `fltst_` dip-column floor mechanism to `srfst_` (recommended: same cause, same bound, DD_RTOL of the vector's largest component); B) loosen only this test's `srfst_` bar (not recommended: a bump); C) exclude the out-of-plane `srfst` columns from the layout test. No tolerance changed; no new head; a fix is a new head and the audit of its delta is owed.

### Row 1, conductor 35: gate on `65a8435` complete, control read, classification corrected (2026-11-19 02:15 CST; harness tool calls 26 of 100; unaudited)

- Gate driver for `65a8435` finished (`runs/20261119_c34_s2t/status_65a8435.txt`: `fast rc=0 01:12:04`, `s2t_e2e_fast rc=1 01:37:46`, `rest_e2e_fast rc=1 02:08:31`, `done`). `rest_e2e_fast_65a8435.log`: 1 failed, 69 passed, 21 skipped, 1696 s; the one failure is `test_install_script_produces_a_working_binary` (`tests/e2e/test_clean_build.py:68`, the machine defaults to `MACHINE=ubuntu`): "Fatal Error: Cannot open included file 'dmumps_struc.h'". Cause read: `scratch/c34_s2t/gate_65a8435.sh` does not export `MACHINE=conda-linux` (CLAUDE.md: build with `-m conda-linux` on this host; `/usr/include/dmumps_struc.h` does not exist), so this is a provisioning artefact of the gate script, not a code regression; the test was not re-run (it rebuilds the binary in place; it runs in the merge-head gate from a fresh clone with `MACHINE=conda-linux`, not now). The gate worktree's `bin/eqquasi-1.24.9` is dated 00:39, so the failed build did not replace it. 21 skips, reasons (data absence, same categories as the #203 gate, none new): no restart files 1, no offfault files 3, no perfault files 4, no profile files 5, no `summary.json` 3, no section 4.3 profiles 5.
- Control for the bp7.8 failure, read from the same log: `test_distributed_loop_matches_petsc_loop_z_and_corner_blocks[bp7.8]` (solver 3 against solver 2, same layout, same `_dd_differences`) PASSED (40.58 s). Correction to my entry of 02:00: the defect is not older than #214, and "test defect" alone is too strong. Reading now: the comparator has no floor for `srfst_` out-of-plane columns (that part stands), and solver 4 produces different round-off than solver 3 in a column that is zero by symmetry (values near 1e-17, own-column maximum 5.8e-17); the in-plane files of the same run pass (25 of 217 entries failed, all in that file family). Not a physics difference at the level seen, but it is solver-4-specific round-off, which the owner's item 10 decides; no tolerance changed. Control still owed: the layout driver's solver 3 against solver 4 on bp5 and bp7 full cycles.
- Not yet read: rank study np 1 (rest) and np 8, layout driver results (`runs/20261119_c34_layout/` holds `test.bp5.qdc.2000_np1_s3` and `_s4` only), liu2020 gold. The gate is not green on `65a8435`; the merge gate is a fresh-clone run at the final head after the owner's answers.

### Row 1, conductor 35: harvest of the rank study and the layout sweep so far (2026-11-19 03:20 CST; harness tool calls 30 of 100; unaudited)

- Landed: #219 (`bc6b283`, run 44482992283 success) and #220 (run 44485244015 success), both board-only, state MERGED read after each merge.
- Rank study finished (`runs/20261118_c32_rank/study.log`, driver exited 03:12): 18 lines, each "0 bad" (6 cases at np 1, 4 and 8 against the regenerated reference, tier comparators, 1e-10 on disp and fault files; np 8 rows: dip90 25, stepover 13, con 13, bp8 24 compared entries). Tier statement: MEASURED, one run per case and rank count; no timing in this reading.
- Layout and full-cycle sweep of the 1.24.9 binary copy, solver 4 against solver 3, `_dd_differences` at DD_RTOL 1e-10 with the 1e-8 operator bar and 1e-14 operator floor (results `runs/20261119_c34_layout/*.json`, mtimes 02:14-03:02, host cotopaxi, bp5 `test.bp5.qdc.2000`, 204,228 equations):

  | layout | steps | differences | s4 s/step | s3 s/step | peak rank RSS s4 / s3 (MB) | sum RSS s4 / s3 (MB) |
  |---|---|---|---|---|---|---|
  | np 1 | 31 | 0 | 0.144 | 0.652 | 4007 / 2108 | 4007 / 2108 |
  | np 3 | 31 | 0 | 0.0669 | 0.301 | 1584 / 897 | 4645 / 2579 |
  | empty-rank (np 4, derived layout) | 31 | 0 | 0.0690 | 0.228 | 1406 / 762 | 4724 / 2788 |
  | full cycle, np 4 | 4483 (stoptag, both) | 0 | 0.0428 | 0.1906 | 1560 / 929 | 5361 / 3438 |

  Reads: the operator check ran and passed in every row (`s2t_ran_and_operator_check: ok`); the distributed-mesh assertion is not met at np 1 and 3 by construction (one rank holds the whole mesh at np 1; the message calls it expected). Setup (`factorization_seconds`) is the same for both solvers (np 4 full cycle 4.68 against 4.71 s); the operator build time (`s2t_build`) is null in the JSONs, the build is inside `wall_in_lock_s` (371 s s4 against 863 s s3 for the whole full-cycle run, wall including lock-free run time only). Time loop at the full cycle: 191.7 s against 854.5 s (4.5x); peak rank memory 1.68x and summed memory 1.56x of solver 3. Not covered by these rows: out-of-plane `srfst_` and `fltst_` families at bp5 pass at 1e-10 without a floor (0 differences), so the item-10 question is bp7-specific; the bp7 full cycle and the bp1002 zone cases are still queued (`test.bp7.qdc.a.10.full_np4_s3` waiting on the heavy lock at 03:02). MEASURED, one run each, host shared with the liu2020 run (4 of 4 ranks plus 3 busy); the s/step ratios are suggestive of speed and not a benchmark.

### Row 1, conductor 35: layout sweep complete, liu2020 regeneration read, 1.24.7 pair launched (2026-11-19 04:40 CST; harness tool calls 42 of 100; unaudited)

- Layout sweep finished (driver wrote `DONE` 03:22; `runs/20261119_c34_layout/*.json`, solver 4 against solver 3, 1.24.9 binary copy, host cotopaxi shared, one run each, `_dd_differences` at DD_RTOL 1e-10, operator check `ok` in both rows): bp7 `test.bp7.qdc.a.10` full cycle at 4 ranks, 1358 steps (stoptag both), 0 differences, 0.0348 against 0.2579 s/step, peak rank RSS 1442 against 844 MB, summed 4873 against 3325 MB; bp1002 zone at 4 ranks, 101 steps (nstep cap, not a full cycle), 0 differences, 0.0119 against 0.1750 s/step, peak rank 793 against 706 MB. With the bp5 rows of the 03:20 entry: 0 differences in 6 of 6 rows. The failure of item 10 (`bp7.8`, solver 4 against solver 2 at 8 ranks) is therefore not reproduced at bp7 np 4 against solver 3 at a full cycle; no solver 4 against solver 3 run exists at 8 ranks, so whether it is 8-rank-specific is NOT measured. bp8 cases and the 12-rank bp1002 layout are not in this sweep (they are in the e2e_fast tests, which ran in the gate: bp1002.12, bp8.gs.8 and bp8.pw.8 passed in the s2t step, only bp7.8 failed).
- liu2020 regeneration finished (prterun 719237 exited 04:29; `runs/20261118_goldF/248/liu600/result/cycle0`, 10000 steps, 1.24.8, 3 ranks, solver 1, time loop 22285.9 s on cotopaxi; fault and disp files written 04:29, `disp.10000.nc` 6,346,312 bytes). Compared with the committed reference by the tier comparators (script `scratch/c35`-side, scratchpad `liucmp.py`, read by me this turn; unaudited): the committed reference was made by `eqquasi-1.14.0` on knox (2026-08-15, `data/liu2020.qdc.kink.600/README.md`, `runInfo.json` version 1.14.0), ten minor versions before the other regenerated sets' reference binary (1.24.7). `fault.*.nc` at the tier's 1e-6 (floor 1e-6): 11 of 12 pass; `fault.05001.nc` FAILS the tier bar, and the worst difference over each variable's own scale rises from 0 (`fault.01001`) through 3.2e-10, 6.1e-8, 7.3e-7 to 2.32e-6 at step 5001, then 5.7e-7, 6.6e-7, 4.9e-7, 3.1e-7, 1.7e-7 (`fault.10000`); `fault.00001.nc` shows 2.9e-1 of own scale (a near-zero variable, as in the bp7 fast case). Text series at the 1e-9 comparator: `global.dat` 17171 of 70000 entries outside (worst 0.25604357 against 0.25604505, row 4482), `fltst_strk000dp009.txt` 8911 of 90000, `cplot_EQquasi.txt` 210 of 82416, `cplot_ruptarea_trac_slip.txt` 29 of 20604. Reading (a reading, not shown): the drift grows with step as in a chaotic interseismic cycle across ten versions, and the reference is 1.14.0, so the comparison cannot separate the 8-byte change from version drift. Under the owner's 16:45 condition (a file that does not agree at single-precision level is not replaced and comes back to this row) NOTHING of liu2020 is replaced, and the 6.3 MB `disp` question (item 5) is moot until the set agrees. NOT committed, nothing under `data/` touched.
- Control launched to separate the two: the same case, 3 ranks, binary `eqquasi-1.24.7` (md5 65b46115..., identical to the `c32-m247` build; copy `scratch/c35/bin/eqquasi-1.24.7`), `runs/20261119_c35_liu247/liu600` (case files copied from the 1.24.8 run, `run.sh` binary name the only edit), launcher PID 1572898, `runs/20261119_c35_liu247/status.txt` ("start 04:31:04"), `flock` on `scratch/locks/heavy.lock`; expected about 6 h (to about 10:30). Decision it enables: if 1.24.7 against 1.24.8 agrees at the tier bar, the 8-byte change is clean and the reference's age is the cause, so liu2020 can be regenerated against the 1.24.7 pair (owner items 5 stays for the `disp` size); if it does not, the 8-byte change moves a chaotic run and comes back to the owner.
- Heavy lock: free except for the run above; step 0 of the operator-build plan (PETSc `-log_view`, no source change) waits for a quiet host and is not started.

### Row 1, conductor 35: checkpoint and owner-brief refresh (2026-11-19 06:20 CST; harness tool calls 47 of 100; unaudited)

- Landed since the 04:40 entry: #222 (`dc948c8`, run 44506518308 success, state MERGED read after). No code PR: #214 is HOLD (items 3 and 10), the gold PR for liu2020 waits on the control below.
- Running, `ps` at 06:14: the 1.24.7 control for liu2020 (`eqquasi-1.24.7`, PIDs 1572946 and 1572947, 1 h 43 min elapsed, `runs/20261119_c35_liu247/liu600/log/cycle0.log` 1744 lines and growing, `status.txt` "start 04:31:04"); expected end about 10:30. No other process of mine.
- Owner brief, refreshed 06:20 (the 09:37 brief is the place; nothing decided; every item has options and a recommendation):
  - Items 1, 2, 4, 6, 7, 8, 9, 10: as on the board (6: MEASURED, accept; 10: the bp7 np 4 full cycle and bp1002 zone pass solver 4 against solver 3 with 0 differences, only bp7.8 at 8 ranks against solver 2 fails, so item 10 stays a comparator-floor question for the `srfst_` out-of-plane columns).
  - Item 3, operator-check floor: the layout sweep is now complete (np 1, 3, an empty-rank layout, bp5 and bp7 full cycles, bp1002 zone: 0 differences, the operator check passed in every row), so the floor 1e-14 of the operator's largest entry was exercised on bp5 and bp7 without a false abort; it was NOT tested against a planted wrong operator (item 7). Recommendation unchanged: A adopt, together with item 7.
  - Item 5, liu2020 `disp` against the 5 MB cap: now conditional on item 11; no `disp` is committed unless liu2020 is regenerated.
  - Item 11 (new), liu2020 reference: the committed files come from `eqquasi-1.14.0` on knox and the 1.24.8 regeneration fails the 1e-6 tier bar at `fault.05001.nc` (2.32e-6 of own scale) and the 1e-9 text bar on `global.dat` (17171 of 70000). A) wait for the 1.24.7 control (about 10:30) and decide on its result: if 1.24.7 against 1.24.8 passes the tier bar, regenerate against the 1.24.7 pair and then apply item 5; if it fails, the 8-byte change moves a chaotic run and the owner decides (recommended: it separates version drift from the dp change at about 6 h of host time already spent). B) leave the 1.14.0 reference and its 4-byte files as they are and record liu2020 as the one set not regenerated (no `disp`, item 5 moot). C) replace the reference with the 1.24.8 run and loosen the bar for this set (not recommended: a loosened bar, and the owner's condition says such a file is not replaced).

### Row 1, conductor 36: start, #223 landed, item 10 control queued (2026-11-19 09:05 CST; harness tool calls 19 of 100; unaudited)

- Budget read as: the owner's "/autopilot clear the board. a week budget. go" (2026-11-18 16:38 CST) covers every open row in `prio` order until 2026-11-25 16:38 CST, row 1's order first. Held until the owner's answers are on the board: items 1-11, ICNTL(14) for native solver 1, a `vista` machine target, any default change, any new file or name. Hand-off at about 80 calls.
- Landed: #223 (squash, head `de1398e`, `gh pr view` state MERGED, merge commit `6963793`; CI `build` SUCCESS after the invoking session's rerun of an apt-mirror failure in "Install dependencies").
- Running, `ps` at 08:53: the 1.24.7 liu2020 control (launcher PID 1572898, `flock` PID 1572922, `eqquasi-1.24.7` ranks PID 1572946-1572948, `runs/20261119_c35_liu247/status.txt` "start 04:31:04", expected about 10:30). Queued behind it on `scratch/locks/heavy.lock`: item 10 control, bp7 `test.bp7.qdc.a.10` at 8 ranks, 31 steps, 1.24.9 binary copy, solvers 4, 3 and 2, comparator `_dd_differences` for 3 against 4, 2 against 4 and 2 against 3 (script `scratch/c36/item10.py`, launcher PID 1960810, results `runs/20261119_c36_item10/result.json`, `DONE` on completion). It answers whether the 1e-17 `srfst_` difference is specific to solver 4, to the solver-2 comparison or to 8 ranks.
- Rows read: 3 (BP8 discrepancy) says "Needs owner diagnosis, not a 3rd same-session retry" after six rejected audit findings; 4 waits on 3; 6 waits on the owner's (a) and (b). No light decided row is workable; row 1 step 0 of the operator-build plan waits for a quiet host after the liu2020 control.

### Row 1, conductor 36: liu2020 1.24.7 control read (2026-11-19 10:50 CST; harness tool calls 34 of 100; unaudited)

- Control finished: `runs/20261119_c35_liu247/status.txt` "rc=0 end 10:32:12" (`eqquasi-1.24.7`, md5 65b46115..., 3 ranks, solver 1, 10000 steps, `runs/20261119_c35_liu247/liu600/result/cycle0`; `disp.10000.nc` 3,626,584 bytes, float32).
- Comparison, 1.24.7 (new slot) against the 1.24.8 run (`runs/20261118_goldF/248/liu600`), tier comparators as in the 04:40 entry (script in the conductor scratchpad, read by me; unaudited): all 25 NetCDF files pass the tier bar at rtol 1e-6 (worst of own scale 1.91e-7 at `fault.05001.nc`, 9.3e-2 at `fault.00001.nc`, a near-zero variable, as before; every other file 3e-8 to 8e-8, which is the float32 rounding of the 8-byte values). Text series at the 1e-9 comparator FAIL: `global.dat` 5326 of 70000 entries (worst 0.25604371 against 0.25604357, row 4482), `fltst_strk000dp009.txt` 1677 of 90000, `peak_sliprate_per_fault.dat` 3363 of 20000, `srfst_strk000st005dp000.txt` 39 of 70000; `cplot_EQquasi.txt` not in the failing list.
- Why, a reading with its evidence: the 1.24.8 change touched computation only by renaming `real(kind=8)` to `real(kind=dp)` (`git show 6934a28 -- src/mesh4num.f90 src/shl.f90`; `dp = selected_real_kind(15,307)` at `src/globalvar.f90:22`, which is kind 8), the two `run.sh` differ only in the binary name, and the case files are copies. The two text series nevertheless differ from row 5 of `global.dat` (9976 of 10000 rows differ exactly) and cross 1e-9 of the column scale at row 3936 (`global.dat`), 4341 (`fltst`), 3495 (`cplot`), growing toward the event at row 4482. So the pair is not bitwise reproducible at 3 ranks while the source arithmetic is the same: round-off seeds that a chaotic interseismic cycle amplifies. NOT measured: whether the same binary run twice also differs (needs another ~6 h of host time; the owner decides, item 11).
- Consequence for item 11: the 8-byte change did not move the NetCDF fields beyond the 1e-6 tier bar (1.24.7 against 1.24.8, 25 of 25); the 1.14.0 reference's 2.3e-6 at `fault.05001.nc` is the same size of drift as this pair's 1.9e-7 plus ten versions of change. The text bar of 1e-9 cannot be met by any two runs of this case on this evidence, so the liu2020 text series would need a different bar or a different comparison point (owner). Recommendation changed from the 06:20 entry: A) regenerate liu2020 from the 1.24.8 run and gate its NetCDF files at the tier's 1e-6 and its text series at the 5e-6-of-scale level the pair shows (a loosened bar for this one chaotic set; owner's call); B) as before, leave the 1.14.0 reference and record the set as not regenerated. Nothing under `data/` touched.

### Row 1, conductor 36: item 10 control, step 0 and 1 of the operator-build plan (2026-11-19 11:10 CST; harness tool calls 55 of 100; unaudited)

- Item 10 control (MEASURED, `runs/20261119_c36_item10/`, scripts `scratch/c36/item10*.py`, 1.24.9 binary copy md5 `b73aa4016e1d`, bp7 `test.bp7.qdc.a.10`, 31 steps, gate-tree `_dd_differences`, DD_RTOL): solver 4 at 8 ranks against solver 3 at 8 ranks 0 differences; against solver 2 at 8 ranks 0; solver 2 at 4 ranks (the test's baseline) against solver 4 at 8 ranks 0, against solver 3 at 8 ranks 0; operator check `ok`. Solver 4 at 8 ranks run 7 more times (3 of them with 60 busy loops added to the host, load about 14 to 60): every run 0 differences against solver 2 at 4 ranks, and row 30 of `srfst_strk000st400dp300.txt` is bit-identical in all of them (`0.3020896E-10 0.4206606E-16 -0.4168300E-16 ...`).
- The gate's failing run is NOT reproduced: its solver-4 output (`.claude/worktrees/agent-a84adddbe0e32be5b/scratch/e2e.dd_s4.s4.bp7.8`, written 01:36) reads `0.3020914E-10 0.4206825E-16 0.4998864E-16 ...` on the same row, while the gate's own solver-2 and solver-3 outputs equal mine. Same binary md5, same `user_defined_params.py` (only the line order of `par.nstep` differs), same `run.sh` apart from nothing, same launcher settings; `runInfo.json` differs only in timings and timestamp. Cause not found: not the binary, config, rank count or host load as far as these 7 runs go. The in-plane difference in that run (6e-6 on one early entry) is larger than round-off, so the item 10 reading "out-of-plane round-off, add a floor" is not established; a floor would be a tolerance change on an unexplained failure. Item 10 recommendation changed: D) no floor; re-run the `[bp7.8]` test alone in the merge-head gate on the final head and, if it fails again, keep the failed case directory and the process environment and hand it to the code auditor before any floor is considered.
- Step 0, bp5 `test.bp5.qdc.2000`, solver 4, 4 ranks, 31 steps, `PETSC_OPTIONS=-log_view`, PETSc 3.25.5 (`runs/20261119_c36_step0/`): operator build 169.8 s of 177.0 s wall, 355 batches of 32 for 11,346 columns; `MatMatTrSolve` 141.1 s (83 % of the build), `SFSetGraph` 25.2 s plus `SFSetUp` 12.8 s (714 calls, two per batch; the build's own SF creation, whether all of it is inside `MatMatTrSolve` not checked), factor 4.5 s, scatter 7 s. Host load 44 at the start, 12 at the end: indicative, not host-exclusive.
- Step 1, same case (`runs/20261119_c36_step1/rec.json`), `-mat_mumps_icntl_27`: baseline repeat build 178.4 s (step 0: 169.8 s, so run-to-run spread about 5 %), `icntl_27 16` 197.7 s (+11 %), `icntl_27 512` 179.6 s (no change); operator check `ok` and 0 differences against the step 0 run in all three; peak rank RSS 1399-1401 MB in all. Reading: ICNTL(27) is not a lever; per-batch SF creation (about 38 s, 21 % of the build, constant across the three runs) and the MUMPS multi-RHS solve itself (about 110 s) are the two terms. `-mat_mumps_icntl_20` not run. Step 2 (batch size as a runtime value, spike branch) needs a new option name: owner-held.

### Row 1, conductor 36: owner "go" read, #214 finish dispatched (2026-11-19 17:20 CST; harness tool calls 60 of 100; unaudited)

- Relay verified: branch `board-owner-relay-11` head `d8e585c` (PR #227, open when read) holds the paragraph "Owner, 2026-11-19 17:07 CST" ("go"). Acting on the paragraph's own reading: operator-check floor 1e-14 of the operator's largest entry with `-eq_s2t_plant_row` and a run-level planted-operator test; `Slip_3` floor 1e-15 m on un-logged values, a log10 entry below the floor on both sides skipped; no `srfst_` floor (item 10 option D); items 1 and 6 accepted. Not touched: item 11 (liu2020), item 5, item 9, step 2 option name, any default.
- Dispatched, one specialist with a 45-call cap, own worktree, branch `c36-s2t-hook` from `origin/c34-s2t` (`65a8435`): the hook, its run-level test, the comparator floors with a measured catch-rate, push as a fast-forward of `c34-s2t` (one new head), no PR action, no tag. Then, mine, in order: audit of the delta, fresh-clone merge-head gate with `MACHINE=conda-linux` (the `[bp7.8]` test re-run alone: a second failure keeps the case and goes to the code auditor first), CI on the exact head, merge #214, release handling (a code audit quoted before any tag).
- Hand-off at about 85 calls with a board checkpoint first; a fresh conductor resumes from this entry.

### Row 1, conductor 36: #214 new head 4f6743b, audit and merge-head gate running (2026-11-19 17:38 CST; harness tool calls 70 of 100; unaudited)

- Landed: #226 and #228 (board only, runs 44572573307 and 44635099783, state MERGED). Owner decisions read: 17:07 "go" and 17:20 "yes" (master `118726c`). Part C reading: the 17:20 paragraph ("no comparator or test changes before the owner sees that proposal") concerns the full-cycle physical metrics only; the 17:07 floors (items 2 and 8) are decided and land as their own small PR right after #214 unless a test needs them in #214 (the invoking session confirmed this reading; my first hold was a stricter reading of the 17:20 text, withdrawn).
- #214 head `4f6743b` (fast-forward of `65a8435`, v1.24.9, 3 files +89/-13, `src/solveTimeLoopPETScDD.f90`, `tests/e2e/test_benchmarks.py`, `docs/user/troubleshooting.md`): `-eq_s2t_plant_row` hook (owning rank doubles its row of the operator before step 1; out-of-range row aborts every rank), first abort line states the floored bound, finite check of `s2tU` on non-full steps (audit minors m1 and m2), run-level test `test_slip_to_traction_planted_operator_row_aborts[row0, out_of_range]` (bp5, 4 ranks, 31 steps, about 175 s each). Specialist's own gate (reported, not mine): build rc 0, fast tier 384 passed 0 skipped, s2t e2e_fast 19 passed 0 skipped 0 failed (1151.7 s, includes bp7.8, bp1002.12, bp8.gs.8, bp8.pw.8). Not done by it: NaN-plant test (m4, needs a new option name, owner-held), `docs/user` option entry (no option list exists). Specialist used 52 tool calls against a cap of 45 (its own report says about 43). Its commit carries an agent-name trailer (the brief asked for it; the standing constraint bars agent names in committed text); history is not rewritten, raised to the invoking session.
- Part C read-only catch-rate (specialist, 90 gold `fltst` files, un-logged, not committed): `Slip_3` below 1e-15 m: 32,606 of 186,872 entries (17.4 %); `Slip_rate_3` below 1e-12: 65,617 (35.1 %); `Shear_stress_3` below 1e-9: 28,774 (15.4 %). Not run: whether any gold-comparator test fails without the floors.
- Running: my merge-head gate on `4f6743b` (fresh clone `scratch/c36/clone_4f6743b`, `env -i`, MACHINE=conda-linux; driver PID 3234183; `runs/20261119_c36_gate_4f6743b/status.txt`; steps build, fast tier, whole e2e_fast under heavy.lock, expected about 1.5 h); the delta audit by one read-only specialist (cap 20 calls). Open after them, in order: CI on the exact head, merge #214 (squash), release handling (code audit quoted before any tag), part C as its own PR, then the full-cycle metrics proposal (owner order: #214 first).
- Hand-off: at about 85 calls; a fresh conductor resumes from this entry.

### Row 1, conductor 37: #214 landed, release cycle started (2026-11-19 18:40 CST; harness tool calls 10 of 100; unaudited)

- Budget read as: the owner's "/autopilot clear the board. a week budget. go" (2026-11-18 16:38 CST) covers every open row in `prio` order until 2026-11-25 16:38 CST, row 1's order first. Held: owner items 5, 9, 11, the operator-build step 2 option name, any default change, any new file or name not agreed; the 17:20 hold covers full-cycle comparator and test changes only. Hand-off at about 80 calls.
- Landed: #214 (`94c0a87`, `par.solver = 4`, v1.24.9; `gh pr view` state MERGED). Master CI run 44646078526 on `94c0a87`: in progress at 18:36 (a conclusion is read before the tag, not assumed). `git diff --stat 4f6743b 94c0a87` lists `PATHWAY_FORWARD.md` only (read by the code auditor), so the audited head is what merged.
- Release cycle for the slip-to-traction promotion (last tag v1.24.5; versions 1.24.6 to 1.24.9 untagged):
  - Code audit of v1.24.5..94c0a87 (read-only, 14 tool calls, unaudited by me beyond reading): 0 BLOCKER, 0 MAJOR. Covered: #194 (1.24.6), #202 (helper move identical to the original apart from whitespace), #203 (no `NF90_REAL`, `real(kind=4)`, `real(kind=8)` left in `src`; tolerances only tightened, 0 removed checks; restart reader reads old 4-byte files), version string 1.24.9; #214's own logic is the audit on that PR (65a8435..4f6743b). MINOR: restart values unrounded from cycle 2 (release notes; measured, see conductor 34 entry: cycle 2 about 1e-6 of scale); liu2020 gold still 4-byte and not in `REGENERATED_8BYTE` (`tests/contract/test_gold_netcdf_shape.py:83-94`, no case reads it); versions 1.24.6 to 1.24.8 never tagged; the ICNTL(14) default is set in both PETSc loops before `KSPSetFromOptions` regardless of the PC type (harmless, workspace only).
  - Rules audit at `94c0a87` (read-only, 8 tool calls): no mechanical blocker; findings: V1 rule 8, 60 gold files rewritten since v1.24.5 (owner-decided regeneration; the board holds a per-file measurement for 53 `fault.*.nc` and the bp5 full set, the 60 against 53 count is not reconciled: release notes quote it); V2 13 added `disp.*.nc` with a reader only through the `disp` category in `tests/e2e/cases.py`; V3 no `docs/user` entry for `-eq_s2t_plant_row` (fixed by this PR, `docs/user/troubleshooting.md`); V4 older machine-local paths in `scripts/bp8Campaign.sh:9`, `scripts/bp8Phase3b.sh:10`; V5 a repository URL carrying a personal handle in `docs/dev/board_history.md:305` (owner's call, as before); V6 a dangling "rule 26" in row 1 text (rule book has 25). Not in the audit's reach and unchecked by me: it reports `-eq_s2t_build` appears nowhere in `src`, `tests` or `docs/user` (this board's older text names it; the option name in source is the check).
  - Unenforceable as written, per the audit (for the rules owner, nothing edited): rule 1 (no command decides agreement), rule 3 per-PR audit and rule 21a quoted audit (no grep), rule 21b stranger-clone evidence (nothing checks it), rule 8 reader check (passes when a test merely references the directory), the weekly `e2e` tier of rule 3 has no schedule in the workflows; the board header "Updated 2026-11-13, at v1.22.0" is stale. Request to the rules owner: re-scope those and the header; the write is theirs.
- Dispatched (2 specialists, worktree-isolated): the `fltst_` out-of-plane floors, first round landed as branch `c37-floors` head `089c49c` (unaudited; `_dd_station_series` only, 11 unit tests, fail with the floor at 0; fast tier 383 passed, 1 failed = the no-binary build-script test in a worktree) and a second round on the gold-vs-run comparator `compare_series` (`tests/e2e/cases.py:655`, the path the owner's floors were meant for) on branch `c37-floors-2`; no `srfst_` floor, no full-cycle change, no tolerance elsewhere.
- Next: this PR; release by the release specialist after the docs PR merges (tag SHA = a master SHA with green CI, stranger clone under `env -i` first); floors PR; the full-cycle metrics proposal on this row; audit follow-ups F2 to F7 and a NaN-plant test.

### Row 1, conductor 37: full-cycle metrics PROPOSAL for the owner (2026-11-19 19:10 CST; harness tool calls 21 of 100; MEASURED once, unaudited; nothing changed in any comparator or test)

Owner 17:20 "yes": 1e-9 / 1e-10 stay for short regression cases; full-cycle comparisons move to physical metrics, proposed here first. Inputs (read by the measuring agent, one mtime per directory): liu2020 v1.24.7 `runs/20261119_c35_liu247/liu600/result/cycle0` (10:32:12) against v1.24.8 `runs/20261118_goldF/248/liu600/result/cycle0` (04:29:43), 3 ranks, solver 1, 10000 steps; bp5 `test.bp5.qdc.2000` (4483 steps) and bp7 `test.bp7.qdc.a.10` (1358 steps), solver 3 against solver 4, 4 ranks, v1.24.9, `runs/20261119_c34_layout/*full_np4_s{3,4}` (02:55 to 03:21). bp1002 has no full-cycle pair (101 steps) and is not measured. Tables and scripts: NOT KEPT. The measuring agent wrote its CSVs and three scripts under its own worktree, which was removed when it finished (gitignored files do not survive), and the copy I attempted found nothing (`cp: cannot stat`); only the numbers in its report survive, and I have not re-derived them. The rerun proposed below regenerates them under `runs/` with the scripts committed or kept.
- Proposed definitions: event = contiguous samples of `global.dat` column 2 (peak slip rate) at or above 1e-3 m/s (the 1 mm/s the solver's own rupture-time output uses); onset and end by log-linear threshold crossing; per event: onset and end time difference (seconds and local time steps), peak slip rate, moment (trapezoid integral of `global.dat` column 3 over the event window), coseismic slip at the station with the largest slip (`fltst` columns 2 and 3, S(end) minus S(onset)); fault-wide max or mean slip per event does not exist (snapshots only every `nt_out` steps), only the end-of-run rupture map `cplot_ruptarea_trac_slip.txt` column 2, compared nodewise.
- Measured worst |B minus A| over three thresholds (3e-4, 1e-3, 3e-3 m/s) and all events: liu2020 onset 3.6e-7 s, end 1.2e-7 s (one float64 spacing at t = 5.87e8 s, interpolation rounding), peak slip rate 4.75e-8 relative, moment 1.09e-9, station slip 1.96e-9 relative, rupture-map slip 2.6e-10 of max; bp5 peak slip rate 9.5e-16, moment 1.8e-15, station slip 1.3e-15; bp7 onset 2.2e-16 s, peak 6.0e-15, moment 1.5e-16, station slip 3.3e-16. `tdyna.txt` identical in all three pairs.
- Proposed tolerances = 10 x measured worst, one significant figure: liu2020 onset 4e-6 s, end 2e-6 s, peak slip rate 5e-7, moment 2e-8, station slip 2e-8 relative; solver 4 against solver 3 (bp5, bp7): peak slip rate 1e-14, moment 2e-14, station slip 2e-14 relative (these are round-off scale, not physical tolerances). The factor 10 is arbitrary and NO same-binary rerun noise floor was measured, so the tolerances have no noise-floor control: the rerun that would give it (same binary, case and rank layout, at least 3 repeats on an idle host, same script) is proposed and not run.
- Problems the owner should see: (1) liu2020's event count is not stable under a factor 3 on the threshold (2, 6, 8 events at 3e-4, 1e-3, 3e-3 m/s; identical in both runs): its events are slip-rate fluctuations inside one rupture lasting about 2.2e3 s at the end of the 5.87e8 s run (the "18.6-year rupture" wording is the whole simulated time), so a threshold event is not a stable unit there; the choice is to merge fluctuations into one event or to use the solver's own `tdyna.txt` window (587154730 to 587154974 s); (2) bp5 starts at t = 0.065 s already at 0.029 m/s, so its event is start-truncated and has no onset to compare; bp5 and bp7 are single-event runs, not multi-cycle; (3) the pairs differ in kind (two binaries on one solver against two solvers on one binary).
- Owner options (to the 09:37 brief): A) adopt these definitions with the liu2020 event defined by the `tdyna.txt` window, tolerances as above after the noise-floor rerun (recommended: the rerun is about 6 h of host time for liu2020, less for bp5 and bp7, and without it the numbers are untested thresholds); B) adopt now at 10 x worst with no rerun; C) different metrics. This informs item 11 (liu2020).
### Row 1, conductor 37: owner-session order on the metrics proposal, lesson for the rules owner (2026-11-19 19:40 CST; harness tool calls 29 of 100; unaudited)

- Order from the invoking session (no new board commit: the invoking session states it is its verification requirement on work the owner's 17:20 paragraph, master 118726c, already ordered, a proposal "measured on the full-cycle pairs already run"; evidence that cannot be reproduced does not meet "measured"; queueing after the release accepted): the metrics proposal cannot rest on numbers whose scripts and tables are gone. After the release, rerun the measurement with script and CSVs written under `runs/<date>_metrics/` (outside any agent worktree), cite those paths on this row, and mark the numbers of the 19:10 entry superseded or confirmed; only post-processing repeats, the runs exist. Status: queued behind the release (gate on `f463902` running, `runs/20261119_c37_gate_f463902/status.txt`); the 19:10 numbers stay marked unreproduced until then.
- Lesson, request to the rules owner (the write is theirs; nothing edited): an agent's worktree was removed when the agent finished and took its gitignored evidence with it (metrics CSVs and three scripts; `cp` found nothing 19:10). Failure chain: the brief asked for the files "in your worktree, not committed", the report pasted summaries, no one copied the files out before the removal; the conductor's copy came after the removal. Proposed rule: a brief that asks for evidence files names a destination under `runs/<date>_<slug>/` in the main checkout's ignored tree (or a committed path), the agent writes there directly, and the conductor lists the destination before accepting the report; an evidence file that exists only inside an agent worktree is treated as not kept.

**Conductor 37, 2026-11-19 22:50 CST [44/100]: v1.24.9 released.** Tag `v1.24.9` on f463902 (lightweight, created with its Release): https://github.com/dunyuliu/EQquasi/releases/tag/v1.24.9. CI run 44647749988 success on that exact SHA. Full tier from a fresh `env -i` clone: `runs/20261119_c37_gate_f463902/status.txt` (build, fast, e2e rc 0); 119 passed, 46 skipped (data absence: bp1002 zone reference 11, restart 1, offfault 4, perfault 8, summary.json 4, section 4.3 profiles 9, profile files 9; one more zone skip than the v1.24.5 note), 0 failed, 13015 s. Stranger clone under `env -i`: README install and `create.newcase` + `case.setup` rc 0 in 397 s, `bash run.sh` printed `Welcome to EQquasi 1.24.9` and reached step 201 before a 150 s timeout (full cycle not run). README gap: conda is not listed as a prerequisite, and `conda env create` stops on a host that already has `eqquasi-petsc`. Notes quote the code audit (0 BLOCKER, 0 MAJOR) and the rules audit (V1 to V6) with open items as open issues. Source read at f463902: a non-finite `s2tU` on a non-full step calls `MPI_ABORT` (`src/solveTimeLoopPETScDD.f90:1692-1699`), closing audit m1. Not found: the 65a8435 audit with m1 to m4 is not on PR 214 (only the 65a8435..4f6743b delta audit is), so m2 to m4 stay unverified. No `tests/release_gate.sh` exists at that SHA. Metrics rerun under `runs/20261120_metrics/` in progress; the proposal numbers stay unreproduced until it reports.

**Conductor 37, 2026-11-19 23:10 CST [46/100]: full-cycle metrics rerun (answers the unreproduced proposal).** Script `runs/20261120_metrics/metrics.py`; CSVs `events.csv`, `summary.csv`, `provenance.csv`, `rupture_map.csv`, `same_binary_groups.csv` and `output.txt` in the same directory (outside any agent worktree; CSV writer, read back). Pairs: liu2020 1.24.7 vs 1.24.8, bp5 and bp7 solver 3 vs 4 at 1.24.9 np 4 (run dirs and `global.dat` mtimes in `provenance.csv`). Worst difference over thresholds 3e-4, 1e-3, 3e-3 m/s (B minus A; moment, peak, slip relative): liu2020 onset 3.576e-7 s, end 1.192e-7 s, peak 4.752e-8, moment 1.577e-7, station slip 1.583e-7, rupture map 2.579e-10; bp5 end 9.55e-12 s, peak 9.534e-16, moment 1.455e-15, slip 1.293e-15 (onset undefined, event touches the first sample); bp7 onset 3.331e-16 s, end 9.104e-15 s, peak 5.996e-15, moment 6.004e-16, slip 4.938e-16. Event counts and sample times identical in every pair; `tdyna.txt` identical. Against the 19:10 table: CONFIRMED liu2020 onset, end, peak, rupture map, event counts 2/6/8, bp5 peak and slip, bp7 peak; SUPERSEDED liu2020 moment (row 1.09e-9, rerun 1.577e-7) and station slip (1.96e-9 vs 1.583e-7), bp5 moment (1.8e-15 vs 1.455e-15), bp7 onset (2.2e-16 vs 3.331e-16), moment (1.5e-16 vs 6.004e-16) and slip (3.3e-16 vs 4.938e-16). Cause of the liu2020 difference is unknown: the lost script's moment and slip definitions are not recoverable, and the rerun depends on the window (at 1e-3 the five complete events give moment 1.35e-9 and slip 1.10e-9; the sixth is truncated by the run end). So the moment and station-slip metrics need a stated window and truncation rule before the owner chooses A/B/C; the other metrics are reproduced. No same-binary noise floor exists for liu2020, bp5 or bp7. Only bp8 reruns exist: 1.24.9 np 2 solver 1 pair bitwise identical (md5 db3309af); 1.24.8 np 1 pair identical (md5 9cb0634e) with one same-version outlier (`20261118_c32_rank/bp8_np1`, md5 a5b50d20, cause unknown). Cheapest floor: one repeat of bp7 `test.bp7.qdc.a.10` solver 4 np 4 on an idle host (not launched; needs a quiet host). No comparator or test was changed.

**Conductor 37, 2026-11-20 00:05 CST [69/100]: audit follow-ups classified, handoff.** PR 236 (docs wrap, F6) merged, CI success at its exact head. Items from the slip-to-traction audit, each checked against the source at master: m1 closed (`src/solveTimeLoopPETScDD.f90:1692-1699` aborts on a non-finite `s2tU`); m2 closed (message at 1678-1685 names the 1e-14 floor); m3 closed as stated (Python-copy test; the Fortran abort has `test_slip_to_traction_planted_operator_row_aborts[row0]`, 184.57 s, e2e_fast, kept out of e2e_ci because the budget was 60 s; a shorter planted case could fit, unmeasured); F2 withdrawn (the abort line already prints once per offending rank with the rank id); F3 in the v1.24.9 notes; F4 deferred (a second planted row costs about 185 s and no case that row 0 misses is known); F7 closed (no 'must not merge' comment remains in `tests/`). OWNER QUESTION (m4, new option name, so owner-held): add `-eq_s2t_plant_nan <global row>` that writes a NaN into that row of the operator output, so a run-level test proves the non-finite abort fires on a non-full step. Options: A) add it as proposed (about 10 src lines, version bump, one e2e test, about 185 s per run); B) reuse `-eq_s2t_plant_row` with a negative row to mean NaN (no new name, less readable); C) leave m4 open (the abort is read in source but not exercised). Recommend A. Still open and not mine: s2t operator-build step 2 (option name owner-held), rows 3, 4, 6 (owner decisions), row 1 steps a to d (need a quiet host; the bp7 same-binary repeat for the noise floor belongs there). Worktrees: `c37-board` is mine; agent worktrees of this run are reaped by the next conductor after a liveness check (`git worktree list`, then `git status --ignored` in each: `scratch/c37_followups/` and `scratch/c37_release/` hold logs and the release notes; `runs/20261120_metrics/` and `runs/20261119_c37_gate_f463902/` are evidence and stay). Remote branches `c37-floors`, `c37-floors-2`, `c37-floors-3`, `docs-wrap-s2t-line` are merged or superseded and can be deleted. No agent is running. Tool calls used 69 of 100.

**Conductor 38, 2026-11-23 [15/100]: owner cleanup (1A, 2A) done; run has no end date.** Relay verified: board PR 239 head `ed0f8b5` holds "Owner, 2026-11-23 11:42 CST" ("no budget. keep working"); the 2026-11-25 end is removed, every other limit stands. Worktrees: classified mechanically (HEAD an ancestor of origin/master, or equal to a merged PR's head; locks held by the host session pid 495646 are live and kept). Removed 9 worktrees (`git worktree remove`, no force; their ignored contents were build products, caches and gate copies; the one board-cited item, `c32-fix` `scratch/gold8`, is already marked not needing copy-out, and the `bin/eqquasi-1.24.7` of `c32-m247` is md5-identical to `scratch/c35/bin/eqquasi-1.24.7`): `agent-a295f81153c12e0d9`, `agent-a3053f613f9486151`, `agent-a8cd03df08d8b00de`, `agent-a94de11514ecea170`, `c32-fix`, `c32-m247`, and the nested `build247`, `disp_before`, `dpbase`. Deleted 6 local branches whose PR head or tip is on master (board-owner-relay-12, c29-194, worktree-agent-a295f81153c12e0d9, fix-skip-gate-r2, c37-floors-3, c33-board-split). Remote: no branch has a merged PR head or master ancestry (24 remain, all unmerged; nothing deleted). Kept and listed: protected `owner-relays-2`, `c23-read`, `lapusta-xi`, `spike-blr`; 11 locked worktrees (host session pid 495646 alive) incl. `agent-a0e30f2d39f683bdd` (live cwd) and merged-branch ones `agent-a2c03a4c132d7d303`, `agent-a07d3a48124b0f9b0`, `agent-ac1731c9d8fd8d5c5`, `agent-ad5d317dcab37a765`; and 33 unlocked worktrees whose branch is neither on master nor a merged PR head (c31-*, c33-s2t-fix, c34-*, c37-floors*, c30-gold8, spike-s2t-*, dd-*, fix-skip-gate with 2 dirty files, gate-full with 1 dirty file): the board marks none superseded, so they stay. Files deleted (55, final grep of each name over PATHWAY_FORWARD.md, docs/, tests/, scripts/, PROJECT_RULES.md found no citation; `rm -rf` was not refused): runs/ 27 folders, 2.13 GB: 20261113_dd-2b 135M; 20261113_dd-2c 209M; 20261113_dd-ownership 213M; 20261113_gate-83 0M; 20261114_bddc-inexact2 69M; 20261114_bddc-inexact-cg 44M; 20261114_blr-wide2 381M; 20261114_dd-2c-fix 1M; 20261114_gate-full_nobuild 0M; 20261114_gate-profser.rc 0M; 20261114_krylov2 63M; 20261114_krylov3 44M; 20261114_krylov4 204M; 20261116_c13-blr-bp1002 182M; 20261116_c13-blr-gs10 51M; 20261116_c13-gamg-guess 52M; 20261116_c14-guess-tight 55M; 20261116_c14-single 27M; 20261117_c20-s2t12 11M; 20261118_gold8 128M; 20261118_goldA 0M; 20261118_goldB 98M; 20261118_goldC 138M; 20261118_goldD 4M; 20261118_goldD2 26M; 20261118_goldE 33M; 20261118_p1gate_aborted_prefix 0M; scratch/ 28 entries, 1.56 GB: bp8.prod10 108M; bp8.pw.10m.wide.v1.21.0 563M; bp8.submit.v1.21.0 19M; build_1.22.0.log 0M; build_r2.log 0M; c32gate 123M; c3.xi0.05 5M; c3.xi0.1 5M; camp.dom1000 9M; camp.dom2000 0M; camp.dom500 5M; camp.xi0.05 2M; deorphan.con 8M; deorphan.dip90 11M; deorphan.stepover 10M; dliu_eqquasi-1-7-0-50m.zip 4M; gate85 416M; gold.bp5.np4 28M; kink.dip90ctrl.dx2000 185M; Liu_coarse50m.zip 2M; NOTES_dd2c.md 0M; ref.bp5.step101 7M; ref.bp7.cycle0 11M; ref.bp7.step101 6M; stranger_v1.24.5.log 0M; sub50.gs 6M; verify.bp8 6M; vista_src_bee60f1 50M; runs/ now 99 GB, scratch/ 13 GB. The remaining 42 scratch/ entries are cited on the board and stay. Next: NaN-plant option PR (v bump, e2e test).
