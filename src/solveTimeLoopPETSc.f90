subroutine solveTimeLoopPETSc
! PETSc-KSP counterpart of solveTimeLoopMUMPS.f90 (PATHWAY_FORWARD row 8a).
!
! Row 8b (added on top): elasticity
! near-null-space (MatSetNearNullSpace via MatNullSpaceCreateRigidBody on
! nodal coordinates) so `-pc_type gamg` has the 6 rigid-body modes it needs,
! a warm start from the previous step's displacement
! (KSPSetInitialGuessNonzero), and per-solve KSP iteration counts logged to
! RUN SUMMARY.
!
! CORRECTED after a source audit caught two false claims in an earlier
! version of this comment, before either was ever compiled:
! (1) KSPSetInitialGuessNonzero is NOT harmless for the hardcoded KSPPREONLY
!     default -- PETSc explicitly rejects a nonzero initial guess on it
!     ("Running KSP of preonly doesn't make sense with nonzero initial
!     guess", confirmed against the installed libpetsc.so 3.25.5). It is
!     now called only when the resolved KSP type is NOT KSPPREONLY.
! (2) `-ksp_type cg -pc_type gamg` is NOT a drop-in runtime option with "no
!     code branch needed": Amat below stores only the upper triangle, which
!     MUMPS's Cholesky path reads correctly by SYM=1 convention but which a
!     real MatMult (what CG/GAMG do) reads as a wrong, non-symmetric
!     matrix. Not fixed in this pass (needs MATSBAIJ storage or both
!     triangles inserted) -- guarded instead: selecting any KSP/PC other
!     than the hardcoded default now aborts loudly rather than silently
!     solving the wrong system (rule 2). See the guard right after
!     KSPSetFromOptions below.
!
! CAUGHT BY THE ROW 8C SWEEP, FIXED IN A LATER COMMIT ON THIS SAME PR: this
! Mat was originally MatCreateSeqAIJ on PETSC_COMM_SELF (rank 0 only), which
! the bp5 sweep showed made solver=2 genuinely serial regardless of -np
! (flat ~1.015-1.018 s/step at ranks 1/2/4/8, while solver=1/MUMPS scaled
! 1.015 -> 0.524). Amat is now MATAIJ on PETSC_COMM_WORLD (see "Row 8c
! distribution" below), still built from the same rank-0-only CRS arrays.
! Verified on theo4/conda-linux: solver=2 now scales too (factorization
! 11.03 -> 3.86s, per-step 1.02 -> 0.59s at ranks 1/2/4/8) -- tracks
! solver=1 closely at 2 and 4 ranks (within ~2%), less closely at 8 ranks
! (+22% factorization, +13% per-step) -- fixed, not just guarded, though
! not yet as fast as the raw MUMPS path at higher rank counts.
!
! Same time-loop physics, same CRS assembly (createMatrixHolderInCRSFormat /
! elemAssembleInCRS, both reused UNMODIFIED and rank-0-only, exactly as in
! the MUMPS path -- confirmed by reading both files: neither is per-rank,
! both operate on the global neq/numel and are only ever called from inside
! an `if (me == 0)` guard). The only thing that differs from
! solveTimeLoopMUMPS.f90 is how the linear system is factorized and solved:
! KSP with -pc_type cholesky -pc_factor_mat_solver_type mumps, i.e. PETSc
! calling the SAME MUMPS library through its own wrapper, so parity with the
! direct DMUMPS path is a same-library, same-algorithm question, gated in
! Phase C.
!
! NOT "-pc_type lu" (PATHWAY_FORWARD row 8a's literal wording): PETSc 3.25.5's
! MatGetFactor_aij_mumps hardcodes mumps->sym = 0 for MAT_FACTOR_LU on a
! MATAIJ regardless of MAT_SPD/MAT_SYMMETRIC -- confirmed by reading
! src/mat/impls/aij/mpi/mumps/impl/imumps.c, and empirically: PCLU's
! -ksp_view reported MUMPS "structural symmetry ... 3%", i.e. it silently
! factorized our upper-triangular-only input as if it were the full
! unsymmetric matrix, which is a different, wrong linear system, not a MUMPS
! roundoff difference. Only MAT_FACTOR_CHOLESKY honors A->spd to select
! mumps->sym = 1, the actual match for mumps_par%SYM = 1 in the MUMPS path.
! PCCHOLESKY is therefore the PETSc call that reaches the same MUMPS code
! path; parity below is gated on that choice, not on the literal PCLU text.
!
! Mirrors the MUMPS JOB=4-once / JOB=3-per-step split as KSPSetUp-once (factor)
! / KSPSolve-per-step (solve), and reports timeUsedInFactorization /
! timeUsedInComputing through the SAME output_run_metadata call and the same
! RUN SUMMARY block, unchanged, so the two solvers' logs are diffable line
! for line.
!
! SUPERSEDED by the row 8c note near MatCreate below (rule 10: kept only for
! history) -- row 8a originally built Amat as a SeqAIJ Mat on
! PETSC_COMM_SELF, rank 0 only, since the CRS arrays feeding it
! (kstiff/ia/ja/num) are still rank-0-only globals. The row 8c sweep showed
! that made solver=2 genuinely serial regardless of -np, so Amat is now
! MATAIJ on PETSC_COMM_WORLD (auto Seq/MPI by comm size) -- see the "Row 8c
! distribution" comment before MatCreate for the current architecture.
! kstiff/ia/ja/num themselves are still rank-0-only, unchanged.
!
! kstiff/ia/ja hold only the upper-triangular half (i <= j) of the symmetric
! stiffness matrix -- createMatrixHolderInCRSFormat sorts every element pair
! (inv,jnv) so jnv >= inv before storing it. That is exactly what MUMPS's
! SYM=1 (symmetric positive definite) expects as input, and is exactly what
! is inserted into the PETSc Mat below (row i-1, col ja(:)-1 with col >= row
! throughout). MAT_SPD is set for the same reason mumps_par%SYM = 1 is set in
! the MUMPS path: it selects MUMPS's SYM=1 code path (not SYM=2's general
! symmetric factorization), which is the bit-exact match the parity gate
! needs. Confirm at runtime with `-ksp_view`, which prints the MUMPS
! "SYM (matrix type)" it actually received.

#include <petsc/finclude/petscksp.h>
    use petscksp
    use globalvar, gv_mat => mat
    implicit none

    Mat                        :: Amat
    Vec                        :: bvec, xvec, coordsVec, xvecSeq
    VecScatter                 :: xScatter
    KSP                        :: ksp
    PC                         :: pc
    MatNullSpace               :: nullsp
    PetscErrorCode             :: perr
    PetscInt                   :: prow, pncols
    PetscInt, allocatable      :: pnnz(:)
    PetscInt, allocatable      :: nnzArr(:)
    PetscInt, allocatable      :: allIdx(:)
    PetscInt                   :: pcols(100)
    PetscScalar                :: pvals(100)
    PetscScalar, pointer       :: parr(:)
    PetscInt                   :: kspIts, totalKSPIts, nKSPSolves
    PetscReal                  :: kspRtol, kspAtol, kspDtol
    PetscInt                   :: kspMaxIts
    PetscBool                  :: isPreonly, isCholesky, isCG, isGAMG
    KSPConvergedReason         :: kspReason
    PetscInt                   :: maxRowNnz
    PetscInt                   :: matLocalRows, matLocalCols

    integer (kind = 4) :: i,j,inv,jnv,ntag,node_num,var,l,k, iiTag
    integer (kind = 4) :: maxRowNnzI4
    real (kind = dp) :: startTime, endTime, timeUsedInFactorization,&
        timeUsedInComputing
    character (len = 50) :: netcdf_outfile, output_type

    call PetscInitialize(PETSC_NULL_CHARACTER, perr)
    CHKERRA(perr)

    call getScalarOnFaultQuant

    if (bp == 8) call pore_pressure_init
    if (bp == 8) call bp8_profile_init

    ! Row 8c distribution (verified on theo4/conda-linux: parity holds and
    ! solver=2 now scales with ranks -- see the top-of-file note). Mat/Vec/
    ! KSP creation, assembly and the KSPSolve/scatter in the time loop below
    ! are now COLLECTIVE PETSc calls on PETSC_COMM_WORLD, run by every rank
    ! -- row 8a's MatCreateSeqAIJ/PETSC_COMM_SELF Mat was confirmed genuinely
    ! serial regardless of -np by the row 8c timing sweep (solver=2 flat at
    ! ~1.015-1.018 s/step for ranks 1/2/4/8, while solver=1/MUMPS scaled
    ! 1.015 -> 0.524). createMatrixHolderInCRSFormat/elemAssembleInCRS and
    ! the MatSetValues/VecSetValues loops below stay rank-0-only, unchanged
    ! from row 8a: kstiff/ia/ja/num/right/resu are still rank-0-only globals,
    ! and distributing THOSE is a bigger rewrite than this pass attempts.
    ! A rank that inserts nothing into a collective Mat/Vec assembly is
    ! valid PETSc usage, not an error.
    !
    ! neq/ndof/numnp/x/id are valid, identical values on every rank at this
    ! point, not just rank 0: mesh4num (sets neq, mesh4num.f90:221) and
    ! meshgen (sets x/id) are both called unconditionally from eqquasi.f90
    ! (lines 44/46, no `if (me == 0)` guard), before solveTimeLoopPETSc is
    ! ever reached -- confirmed by reading eqquasi.f90, not assumed.
    if (me.eq.0) then
        write(*,*) '= Building Stiffness Matrix in CRS format      ='
        call createMatrixHolderInCRSFormat
        call elemAssembleInCRS
        write(*,*) '= Converting from CRS to PETSc Mat/Vec format  ='
    endif

    ! Preallocation, FIRST PASS: every row gets the SAME maxRowNnz slots on
    ! every rank (both the Seq and MPI preallocation calls below -- whichever
    ! does not match Amat's actual runtime type, decided by MatSetType/comm
    ! size, is a documented no-op, so both are safe to call unconditionally).
    ! This is deliberately not tight: exact per-rank diagonal/off-diagonal
    ! counts need PETSc's row-ownership split, which is only decided inside
    ! MatSetSizes/MatSetType below, so rank 0's pnnz(:) can't be sliced for
    ! it ahead of time without more bookkeeping than this pass attempts.
    ! Correct first (uniform over-allocation wastes memory, never produces a
    ! wrong answer), fast later.
    maxRowNnzI4 = 0
    if (me == 0) then
        allocate(pnnz(neq))
        do i = 1, neq
            pnnz(i) = num(i)
        enddo
        maxRowNnzI4 = maxval(pnnz)
        deallocate(pnnz)
    endif
    call MPI_BCAST(maxRowNnzI4, 1, MPI_INT, 0, MPI_COMM_WORLD, IERR)
    maxRowNnz = maxRowNnzI4

    call MatCreate(PETSC_COMM_WORLD, Amat, perr)
    CHKERRA(perr)
    call MatSetSizes(Amat, PETSC_DECIDE, PETSC_DECIDE, neq, neq, perr)
    CHKERRA(perr)
    call MatSetType(Amat, MATAIJ, perr)
    CHKERRA(perr)
    ! Row 8b, real bug found by actually running the 8-rank config (not
    ! assumed from the 1/2/4-rank passes): with no explicit block size,
    ! PETSc's PETSC_DECIDE row partitioning does not have to land on a
    ! multiple of ndof=3 per rank, and coordsVec (given block size 3
    ! explicitly, below) then gets ITS OWN independent PETSC_DECIDE split
    ! that can disagree with Amat's -- PCSetData_AGG aborted at 8 ranks
    ! with "30 != matrix size 25529", a local-size mismatch between the
    ! near-null-space vector and the matrix GAMG was actually handed.
    ! Fixed by making Amat's own partition block-respecting, then deriving
    ! every other PETSC_COMM_WORLD vector's local size FROM Amat's actual
    ! partition (below) instead of each independently guessing PETSC_DECIDE.
    call MatSetBlockSize(Amat, ndof, perr)
    CHKERRA(perr)
    ! Explicit arrays, not PETSC_NULL_INTEGER: a second audit pass flagged
    ! that this PETSc Fortran interface (3.25.5) may require the
    ! array-specific PETSC_NULL_INTEGER_ARRAY sentinel here rather than the
    ! scalar PETSC_NULL_INTEGER, and could not confirm which without
    ! compiling. Sidestepped entirely by passing a real, generously-sized
    ! (length neq, safe for any rank's actual local row count) uniform-fill
    ! array instead of a null placeholder -- unambiguous across Fortran
    ! interface versions. The leading scalar `0` in both calls is the
    ! "uniform nz" convenience argument, ignored once the array is
    ! non-null; using it too would have hit the exact same sentinel
    ! question, so both preallocation calls stay fully explicit.
    ! Row 8b: rows now also receive mirrored off-diagonal entries from OTHER
    ! rows' upper-triangle lists (see the MatSetValues loop below), so a
    ! row's true final count can exceed its own upper-triangle-only
    ! maxRowNnz. Doubled as a generous (not exact) safety margin --
    ! PETSc reallocates and warns rather than corrupting anything if this
    ! is still short, so "correct first" holds even if this bound is loose.
    allocate(nnzArr(neq))
    nnzArr = 2 * maxRowNnz
    call MatSeqAIJSetPreallocation(Amat, 0, nnzArr, perr)
    CHKERRA(perr)
    call MatMPIAIJSetPreallocation(Amat, 0, nnzArr, 0, nnzArr, perr)
    CHKERRA(perr)
    deallocate(nnzArr)

    if (me.eq.0) then
        do i = 1, neq
            prow   = i - 1
            pncols = num(i)
            ! Rule 2: pcols/pvals are fixed-size(100) buffers; an
            ! unchecked overflow here would silently corrupt memory rather
            ! than fail. Upper-triangle hex8 rows fit today (at most 81),
            ! but a future mesh/element change should not find out by
            ! crashing somewhere unrelated.
            if (pncols > 100) then
                write(*,*) 'PETSc solver: row', i, 'has', pncols, &
                    'entries, more than the fixed pcols/pvals(100) ', &
                    'buffers hold. Stopping rather than overflowing them.'
                call MPI_ABORT(MPI_COMM_WORLD, 1, IERR)
            endif
            do j = 1, pncols
                pcols(j) = ja(ia(i)+j-1) - 1
                pvals(j) = kstiff(ia(i)+j-1)
            enddo
            call MatSetValues(Amat, 1, [prow], pncols, pcols(1:pncols), &
                pvals(1:pncols), INSERT_VALUES, perr)
            CHKERRA(perr)
            ! Row 8b storage fix: kstiff/ia/ja hold only the upper triangle
            ! (col >= row), which is exactly what MUMPS's SYM=1 factorization
            ! needs and exactly what row 8a verified byte-for-byte -- PETSc's
            ! aij-to-mumps symmetric conversion extracts the col>=row half of
            ! whatever is stored regardless of what else is present, so
            ! mirroring the OFF-DIAGONAL entries below into their transposed
            ! (col,row) position adds information for a real MatMult
            ! (CG/GAMG) without changing what the Cholesky/MUMPS path reads.
            ! That claim is exactly what the immediately-following parity
            ! re-run checks -- this is not assumed to hold, it is verified
            ! every time this file is tested. The diagonal (j=1, since
            ! createMatrixHolderInCRSFormat sorts col>=row ascending) is
            ! never mirrored -- it would double-insert the same value at
            ! the same (row,row) position, which INSERT_VALUES tolerates
            ! (last write wins, same value) but which is pointless and
            ! confusing to read.
            do j = 1, pncols
                if (pcols(j) /= prow) then
                    call MatSetValues(Amat, 1, [pcols(j)], 1, [prow], &
                        [pvals(j)], INSERT_VALUES, perr)
                    CHKERRA(perr)
                endif
            enddo
        enddo
    endif
    call MatAssemblyBegin(Amat, MAT_FINAL_ASSEMBLY, perr)
    CHKERRA(perr)
    call MatAssemblyEnd(Amat, MAT_FINAL_ASSEMBLY, perr)
    CHKERRA(perr)
    ! Both MAT_SPD (unchanged -- still selects MUMPS's SYM=1 code path
    ! exactly as before) and MAT_SYMMETRIC (new -- now honest, since both
    ! triangles are genuinely stored, and this is what lets GAMG/CG treat
    ! the matrix as symmetric for their own algorithms).
    call MatSetOption(Amat, MAT_SPD, PETSC_TRUE, perr)
    CHKERRA(perr)
    call MatSetOption(Amat, MAT_SYMMETRIC, PETSC_TRUE, perr)
    CHKERRA(perr)

    ! The actual local row/col count PETSc gave THIS rank of Amat, now that
    ! it is block-respecting (MatSetBlockSize above) -- every other
    ! PETSC_COMM_WORLD Vec below is sized from THIS, not its own
    ! independent PETSC_DECIDE guess, so they are all guaranteed to agree
    ! with Amat's partition (the fix for the 8-rank PCSetData_AGG crash).
    call MatGetLocalSize(Amat, matLocalRows, matLocalCols, perr)
    CHKERRA(perr)

    ! Row 8b: elasticity near-null-space for -pc_type gamg (a runtime
    ! option; harmless/ignored by the hardcoded PCCHOLESKY default below).
    ! GAMG needs the 6 rigid-body modes (3 translations + 3 rotations,
    ! ndof=3) to build good coarse grids; MatNullSpaceCreateRigidBody
    ! wants one Vec of nodal coordinates, block size ndof, laid out in
    ! the SAME order as the matrix rows it will be attached to.
    !
    ! That ordering assumption holds here: meshgen.f90's equation
    ! numbering (search `neq0=neq0+1; id(i1,nnode)=neq0`) assigns all
    ! ndof=3 equations of a free node consecutively, in one inner loop,
    ! before moving to the next node -- a node either contributes a full
    ! contiguous 3-block or none at all (boundary/fault-boundary nodes
    ! get negative id() codes and never enter the equation count). So
    ! id(1,i), id(1,i)+1, id(1,i)+2 are exactly id(1,i), id(2,i), id(3,i)
    ! for every free node i, and a length-neq vector indexed by id(1,i)
    ! is a valid blocked coordinate vector -- checked defensively below
    ! rather than assumed silently (rule 2: no silent wrong behavior).
    if (mod(neq, ndof) /= 0) then
        ! Rule 2: no warn-and-proceed. An audit flagged the original
        ! "skip and continue with a weaker null space" here as exactly the
        ! disallowed pattern -- fail loudly instead, since a silently
        ! degraded GAMG hint is a correctness question a caller must decide,
        ! not one this code should quietly decide for them.
        if (me == 0) write(*,*) 'PETSc solver: neq is not a multiple of ', &
            'ndof; the rigid-body near-null-space assumes every free ', &
            'node contributes a full ndof-block of equations, which does ', &
            'not hold here. Stopping rather than silently building a ', &
            'near-null-space that GAMG would use incorrectly.'
        call MPI_ABORT(MPI_COMM_WORLD, 1, IERR)
    else
        call VecCreate(PETSC_COMM_WORLD, coordsVec, perr)
        CHKERRA(perr)
        call VecSetSizes(coordsVec, matLocalRows, neq, perr)
        CHKERRA(perr)
        call VecSetBlockSize(coordsVec, ndof, perr)
        CHKERRA(perr)
        call VecSetFromOptions(coordsVec, perr)
        CHKERRA(perr)
        if (me == 0) then
            do i = 1, numnp
                if (id(1,i) > 0) then
                    do j = 1, ndof
                        call VecSetValue(coordsVec, id(1,i)-2+j, x(j,i), &
                            INSERT_VALUES, perr)
                        CHKERRA(perr)
                    enddo
                endif
            enddo
        endif
        call VecAssemblyBegin(coordsVec, perr)
        CHKERRA(perr)
        call VecAssemblyEnd(coordsVec, perr)
        CHKERRA(perr)
        call MatNullSpaceCreateRigidBody(coordsVec, nullsp, perr)
        CHKERRA(perr)
        call MatSetNearNullSpace(Amat, nullsp, perr)
        CHKERRA(perr)
        call MatNullSpaceDestroy(nullsp, perr)
        CHKERRA(perr)
        call VecDestroy(coordsVec, perr)
        CHKERRA(perr)
    endif

    call VecCreate(PETSC_COMM_WORLD, bvec, perr)
    CHKERRA(perr)
    call VecSetSizes(bvec, matLocalRows, neq, perr)
    CHKERRA(perr)
    call VecSetFromOptions(bvec, perr)
    CHKERRA(perr)
    call VecDuplicate(bvec, xvec, perr)
    CHKERRA(perr)
    ! Rank 0 still does all the serial FEM bookkeeping (bound_load, faulting,
    ! output) against plain Fortran arrays (`right`, `resu`), unchanged from
    ! row 8a -- only the linear solve is distributed. allIdx is the fixed
    ! 0..neq-1 global index list rank 0 uses to push `right` through
    ! VecSetValues each step; built once here, not per step.
    if (me == 0) then
        allocate(allIdx(neq))
        do i = 1, neq
            allIdx(i) = i - 1
        enddo
    endif
    ! Gathers the distributed xvec back to a full copy on rank 0 after every
    ! solve (xvecSeq) -- the standard PETSc idiom for a solver embedded in an
    ! otherwise-serial legacy driver. Built once; reused via
    ! VecScatterBegin/End every step in the time loop below.
    call VecScatterCreateToZero(xvec, xScatter, xvecSeq, perr)
    CHKERRA(perr)
    ! Row 8b warm start: KSPSetInitialGuessNonzero makes KSPSolve read
    ! xvec's current contents as the initial guess instead of starting
    ! from zero every step. xvec is never zeroed or overwritten between
    ! solves in the time loop below (only read via the scatter above into
    ! `resu`), so it already carries the previous step's solution into
    ! the next KSPSolve call once this is on -- the one explicit
    ! VecZeroEntries below is only for the very first solve, so that
    ! "previous step's solution" is a defined zero rather than
    ! uninitialized memory the first time through.
    call VecZeroEntries(xvec, perr)
    CHKERRA(perr)

    call KSPCreate(PETSC_COMM_WORLD, ksp, perr)
    CHKERRA(perr)
    call KSPSetOperators(ksp, Amat, Amat, perr)
    CHKERRA(perr)
    call KSPSetType(ksp, KSPPREONLY, perr)
    CHKERRA(perr)
    call KSPGetPC(ksp, pc, perr)
    CHKERRA(perr)
    ! NOT PCLU: read against PETSc 3.25.5 src/mat/impls/aij/mpi/mumps
    ! (MatGetFactor_aij_mumps) shows MAT_FACTOR_LU on a MATAIJ hardcodes
    ! mumps->sym = 0 (general unsymmetric) UNCONDITIONALLY, ignoring
    ! MAT_SPD/MAT_SYMMETRIC entirely -- confirmed empirically too: with
    ! PCLU, MUMPS's own -ksp_view reported "structural symmetry ... 3%"
    ! (i.e. it silently treated our upper-triangular-only input as if it
    ! were the full unsymmetric matrix, factorizing the wrong system and
    ! producing a plausible-looking but physically wrong solution). Only
    ! MAT_FACTOR_CHOLESKY honors A->spd to select mumps->sym = 1, which is
    ! the SYM=1 the raw mumps_par%SYM = 1 call actually uses. PCCHOLESKY
    ! is therefore the PETSc equivalent that reaches the same MUMPS code
    ! path, not PCLU as PATHWAY_FORWARD's row 8a note literally says.
    call PCSetType(pc, PCCHOLESKY, perr)
    CHKERRA(perr)
    call PCFactorSetMatSolverType(pc, MATSOLVERMUMPS, perr)
    CHKERRA(perr)
    ! Runtime overrides (e.g. -ksp_view, or a different -pc_type for row
    ! 8b) are applied on top of the hardcoded defaults above, never
    ! instead of them.
    call KSPSetFromOptions(ksp, perr)
    CHKERRA(perr)

    call PetscObjectTypeCompare(ksp, KSPPREONLY, isPreonly, perr)
    CHKERRA(perr)
    call PetscObjectTypeCompare(pc, PCCHOLESKY, isCholesky, perr)
    CHKERRA(perr)
    call PetscObjectTypeCompare(ksp, KSPCG, isCG, perr)
    CHKERRA(perr)
    call PetscObjectTypeCompare(pc, PCGAMG, isGAMG, perr)
    CHKERRA(perr)
    ! CORRECTNESS GUARD, WHITELIST (widened again after the storage fix
    ! above): Amat now stores BOTH triangles, so a real MatMult (CG/GAMG)
    ! reads the true symmetric stiffness matrix, not half of it -- the
    ! silent-wrong-answer risk this guard existed to prevent is gone for
    ! that combination specifically. Still a whitelist, not a blanket
    ! allow: only the two combinations actually verified against a
    ! reference are permitted --
    !   (a) KSPPREONLY + PCCHOLESKY (+ MUMPS via PCFactorSetMatSolverType
    !       above): the original row 8a bit-exact parity path, unaffected
    !       by the storage fix (PETSc's aij-mumps symmetric conversion
    !       extracts col>=row regardless of what else is stored -- verified
    !       by re-running row 8a's exact parity gate after this change, not
    !       assumed).
    !   (b) KSPCG + PCGAMG (+ the near-null-space below): row 8b's own
    !       gate, tolerance-level parity against the fast references,
    !       iteration count and ksp_rtol reported in RUN SUMMARY.
    ! Anything else (-pc_type lu, -pc_type jacobi, -ksp_type gmres with the
    ! default PC, etc.) is unverified against this matrix's actual
    ! numerical behavior and stays refused (rule 2) rather than assumed
    ! safe by analogy.
    if (.not. ((isPreonly .and. isCholesky) .or. (isCG .and. isGAMG))) then
        if (me == 0) write(*,*) 'PETSc solver: this KSP/PC combination ', &
            'has not been verified against a reference. Only ', &
            'KSPPREONLY+PCCHOLESKY (MUMPS, bit-exact parity) and ', &
            'KSPCG+PCGAMG (tolerance-level parity) are. Refusing rather ', &
            'than computing an unverified answer.'
        call MPI_ABORT(MPI_COMM_WORLD, 1, IERR)
    endif
    ! Row 8b warm start: a no-op under KSPPREONLY (PETSc explicitly rejects
    ! a nonzero initial guess there, "Running KSP of preonly doesn't make
    ! sense with nonzero initial guess", confirmed against libpetsc.so
    ! 3.25.5) and the real point once CG+GAMG is selected.
    if (.not. isPreonly) then
        call KSPSetInitialGuessNonzero(ksp, PETSC_TRUE, perr)
        CHKERRA(perr)
    endif
    ! Reported in RUN SUMMARY (rank 0 only reads these local copies below;
    ! the call itself just reads ksp's already-collectively-set options,
    ! not a collective operation needing every rank's result).
    call KSPGetTolerances(ksp, kspRtol, kspAtol, kspDtol, kspMaxIts, perr)
    CHKERRA(perr)

    if (me == 0) call initOnFaultKinematics

    ! KSPSetUp (analysis + factorization) is collective on PETSC_COMM_WORLD
    ! now -- gating it behind `if (me.eq.0)` (row 8c's first pass) hangs
    ! every other rank waiting on a call that never comes, caught by a
    ! second audit pass before this was ever compiled. cpu_time is a
    ! per-process wall/CPU clock, still meaningful read on rank 0 alone.
    call cpu_time(startTime)
    call KSPSetUp(ksp, perr)   ! Combines analysis + factorization, once.
    CHKERRA(perr)
    call cpu_time(endTime)
    timeUsedInFactorization = endTime - startTime

    stoptag = 0 ! set stoptag to FALSE.
    totalKSPIts = 0
    nKSPSolves = 0

    call cpu_time(startTime)
    do it = 1, nstep
        if (stoptag == 1) exit ! exit EQquasi if stoptag is TRUE.

        if (me == 0) then
            resu_1 = resu

            if (bp == 7 .and. icstart == 1) then ! for the first cycle of bp7, if time<nuct, use dt for dtev1.
                if (it < (nuct/dt)) then
                    dtev1 = dt
                endif
            endif

            if (dtmax > 0.0d0) dtev1 = min(dtev1, dtmax)

            time = time + dtev1

            ! bp8: advance the along-fault pore pressure to t+dtev1 before
            ! faulting consumes it through fric(FR_PORE_DP,:,:).
            if (bp == 8) call pore_pressure_update(dtev1)

            if (mod(it,nhplt) == 1 .and. me ==0) then
                write(*,*) '=                                                                   ='
                write(*,*) '=     Current time =                                                ='
                write(*,'(X,A,40X,E15.7,4X,A)') '=',  time/60.0d0/60.0d0/24.0d0/365.0d0, 'year'
                write(*,*) '=     time step =                                                   ='
                write(*,'(X,A,40X,i7,4X,A)') '=',  it
                write(*,*) '=     dt =                                                          ='
                write(*,'(X,A,40X,E15.7,4X,A)') '=',  dtev1, 'seconds'
                write(*,*) '=     pma = MA/KU =                                                 ='
                write(*,'(X,A,40X,E15.7,4X,A)') '=',  pma
                write(*,*) '=     maximum sliprate =                                            ='
                write(*,'(X,A,40X,E15.7,4X,A)') '=',  maxSlipRate, 'm/s'
            endif

            call exitCriteria

            if(ndout>0) then
                dout(1,it)=time
                do i=1,ndout
                    j=idhist(1,i)
                    if(j<=0) j=1  !avoid zero that cannot be used below
                        k=idhist(2,i)
                        l=idhist(3,i)
                    if(l==1) then
                        dout(i+1,it)=cons(k,j)
                    elseif(l==2) then
                        dout(i+1,it)=consv(k,j)
                    elseif(l==3) then
                        dout(i+1,it)=consa(k,j)
                    endif
                enddo
            endif
        endif

        do iiTag = 0, 1
            if (me == 0) then
                ! When iiTag==1, !V*(t+1) OBTAINED.
                ! UPDATE BOUNDARY U**(t+1),THEN DECLARE [CONSTRAIN] BY U**(t+1).

                itag = iiTag
                ! COMPUTE U*(t+1) FOR THE WHOLE VOLUME'
                ! COMPUTE U**(t+1) FOR THE WHOLE VOLUME.
                call bound_load
                ! Row 8c distribution: bvec is now a PETSC_COMM_WORLD Vec, so
                ! VecGetArray only exposes rank 0's OWN local slice, not the
                ! whole length-neq array -- `parr = right` (row 8a) silently
                ! wrote into the wrong-sized buffer once bvec stopped being
                ! sequential. VecSetValues by explicit global index (allIdx,
                ! built once above) is the correct way for rank 0 to push the
                ! full RHS into a distributed vector; the collective
                ! Assembly call below routes each entry to its owning rank.
                call VecSetValues(bvec, neq, allIdx, right(1:neq), &
                    INSERT_VALUES, perr)
                CHKERRA(perr)
            endif
            call VecAssemblyBegin(bvec, perr)
            CHKERRA(perr)
            call VecAssemblyEnd(bvec, perr)
            CHKERRA(perr)

            call MPI_BCAST(stoptag, 1, MPI_INT, 0, MPI_COMM_WORLD, IERR)

            ! KSPSolve is collective -- every rank must call it now that Amat/
            ! bvec/xvec are PETSC_COMM_WORLD objects, not just rank 0 (row
            ! 8a's `if (me.eq.0)` guard here would hang every other rank
            ! waiting on a call that never comes).
            call KSPSolve(ksp, bvec, xvec, perr)
            CHKERRA(perr)
            ! Rule 2: a solve that did not converge is not a solution.
            call KSPGetConvergedReason(ksp, kspReason, perr)
            CHKERRA(perr)
            ! KSPConvergedReason is a Fortran derived type (type(eKSPConvergedReason),
            ! petsc/finclude/petscksp.h), not a plain integer -- petscksp.mod
            ! only defines ==/ /= for it, not <, so the raw integer code
            ! must be read via its %v component. Caught by an audit reading
            ! the actual .mod/finclude headers before this was compiled.
            if (kspReason%v < 0) then
                if (me == 0) write(*,*) 'PETSc solver: KSPSolve did NOT ', &
                    'converge (KSPConvergedReason =', kspReason%v, &
                    ') at step', it, 'iiTag', iiTag, &
                    '-- stopping rather than using an unconverged solve.'
                call MPI_ABORT(MPI_COMM_WORLD, 1, IERR)
            endif
            ! Row 8b: iterations/step, meaningful once -ksp_type cg is
            ! selected at runtime (KSPPREONLY's default below always
            ! reports 1). Logged even in the LU/Cholesky-default case,
            ! since a stray non-1 there would itself be a bug worth
            ! seeing.
            call KSPGetIterationNumber(ksp, kspIts, perr)
            CHKERRA(perr)

            ! Gather the distributed solution back to a full copy on rank 0
            ! (xvecSeq, via the VecScatterCreateToZero context built once
            ! above) -- VecGetArray directly on xvec (row 8a) only ever
            ! exposed rank 0's own local slice, which is a real bug once
            ! xvec is no longer sequential, not merely a style change.
            call VecScatterBegin(xScatter, xvec, xvecSeq, INSERT_VALUES, &
                SCATTER_FORWARD, perr)
            CHKERRA(perr)
            call VecScatterEnd(xScatter, xvec, xvecSeq, INSERT_VALUES, &
                SCATTER_FORWARD, perr)
            CHKERRA(perr)

            if (me.eq.0) then
                totalKSPIts = totalKSPIts + kspIts
                nKSPSolves = nKSPSolves + 1

                call VecGetArray(xvecSeq, parr, perr)
                CHKERRA(perr)
                resu = parr
                call VecRestoreArray(xvecSeq, parr, perr)
                CHKERRA(perr)

                if (iiTag==0) then
                    do i=1,numnp
                        do j=1,ndof
                            if (id(j,i)>0) then !Only nodes that have equation number have been updated here.
                                constmp(j,i)=resu(id(j,i))
                            endif
                        enddo
                    enddo
                elseif (iiTag==1) then
                    do i=1,numnp
                        do j=1,ndof
                            if (id(j,i)>0) then
                                cons(j,i) = resu(id(j,i))
                            else
                                cons(j,i) = constmp(j,i)
                            endif
                        enddo
                    enddo
                endif

                !COMPUTE F*=KU*(t+1)=[KSTIFF].DOT.[CONSTRAINTMP]'
                !COMPUTE F**=KU**(t+1)=[KSTIFF].DOT.[CONSTRAIN]'
                call bound_ft_ku
                !GET FIRST PREDICTIONS OF V*(t+1).
                !GET SECOND PREDICTIONS OF V**(t+1), AND DECLARE V(t+1)=V**(t+1).'
                call faulting
            endif
        enddo

        if (me == 0) then
            if (bp == 8) call bp8_profile_record

            status0 = status1
            dtev1 = max(dt,int8(dtev/dt)*dt)

            do i=1,numnp
                do j=1,ndof
                    if (id(j,i)>0) then
                        consv(j,i) = (resu(id(j,i)) - resu_1(id(j,i)))/dtev1
                    endif
                enddo
            enddo
            globaldat(1,it) = time
            globaldat(2,it) = maxSlipRate
            globaldat(3,it) = totMomRate
            globaldat(4,it) = totTaoRuptArea
            globaldat(5,it) = totSlipRuptArea
            globaldat(6,it) = totRuptArea
            globaldat(7,it) = totMomRateVW

            if (mod(it,nt_output_stress) == 1 .or. (stoptag == 1) .or. (it==nstep)) then
                write(proc_str,'(I5.5)') it
                netcdf_outfile = trim(outDir)//'disp.'//trim(proc_str)//'.nc'
                output_type = 'disp'
                call netcdf_write(netcdf_outfile, output_type)

                netcdf_outfile = trim(outDir)//'fault.'//trim(proc_str)//'.nc'
                call netcdf_write_on_fault(netcdf_outfile)
            endif
            ! If exiting, write again the restart files.
            if ((stoptag == 1) .or. (it==nstep)) then
                netcdf_outfile = trim(outDir)//'disp.r.nc'
                output_type = 'disp'
                call netcdf_write(netcdf_outfile, output_type)

                netcdf_outfile = trim(outDir)//'fault.r.nc'
                call netcdf_write_on_fault(netcdf_outfile)
            endif
        endif
    enddo
    call cpu_time(endTime)
    timeUsedInComputing = endTime-startTime

    if (me == 0) then
        write(*,*) '====================================================================='
        write(*,*) '=                        RUN SUMMARY                                ='
        write(*,'(X,A,40X,i7,4X,A)')     '= MPI ranks                = ', nprocs, '='
        write(*,'(X,A,40X,i7,4X,A)')     '= Nodes                    = ', numnp, '='
        write(*,'(X,A,40X,i7,4X,A)')     '= Elements                 = ', numel, '='
        write(*,'(X,A,40X,i7,4X,A)')     '= Equations                = ', neq, '='
        do i = 1, ntotft
            write(*,'(X,A,I3,A,E13.5,A,I7,4X,A)') '= Fault ', i, ' (y = ', fltxyz(1,2,i), &
                ' m) nodes = ', nftnd(i), '='
        enddo
        write(*,'(X,A,40X,i7,4X,A)')     '= Time steps completed     = ', it-1, '='
        write(*,'(X,A,40X,E15.7,4X,A)')  '= Simulated time           = ', time/86400.0d0, 'days'
        write(*,'(X,A,40X,E15.7,4X,A)')  '= Time loop                = ', timeUsedInComputing, 'seconds'
        write(*,'(X,A,40X,E15.7,4X,A)')  '= Factorization            = ', timeUsedInFactorization, 'seconds'
        write(*,'(X,A,40X,E15.7,4X,A)')  '= Seconds per step         = ', timeUsedInComputing/max(1,it-1), 'seconds'
        write(*,'(X,A,40X,E15.7,4X,A)')  '= Final max slip rate      = ', maxSlipRate, 'm/s'
        write(*,'(X,A,40X,E15.7,4X,A)')  '= Avg KSP iterations/solve = ', &
            dble(totalKSPIts)/dble(max(1,nKSPSolves)), '='
        write(*,'(X,A,40X,E15.7,4X,A)')  '= KSP rtol                 = ', &
            kspRtol, '='
        write(*,*) '====================================================================='
        call output_run_metadata(timeUsedInComputing, timeUsedInFactorization)
    endif

    ! Row 8c: Amat/bvec/xvec/ksp are PETSC_COMM_WORLD objects now, so their
    ! destroy calls are collective too -- gating them behind `if (me==0)`
    ! (row 8a) would hang every other rank.
    call VecScatterDestroy(xScatter, perr)
    CHKERRA(perr)
    call VecDestroy(xvecSeq, perr)
    CHKERRA(perr)
    call VecDestroy(bvec, perr)
    CHKERRA(perr)
    call VecDestroy(xvec, perr)
    CHKERRA(perr)
    call KSPDestroy(ksp, perr)
    CHKERRA(perr)
    call MatDestroy(Amat, perr)
    CHKERRA(perr)

    call PetscFinalize(perr)
    CHKERRA(perr)

end subroutine solveTimeLoopPETSc
