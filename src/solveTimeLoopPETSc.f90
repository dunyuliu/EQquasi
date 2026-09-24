subroutine solveTimeLoopPETSc
! PETSc-KSP counterpart of solveTimeLoopMUMPS.f90 (PATHWAY_FORWARD row 8a).
!
! Row 8b (added on top, NOT YET COMPILED OR TESTED as of this commit -- the
! host was mid-timing-sweep for row 8c and a build would have perturbed it;
! compile and verify before trusting anything below this note): elasticity
! near-null-space (MatSetNearNullSpace via MatNullSpaceCreateRigidBody on
! nodal coordinates) so `-pc_type gamg` has the 6 rigid-body modes it needs,
! a warm start from the previous step's displacement
! (KSPSetInitialGuessNonzero), and per-solve KSP iteration counts logged to
! RUN SUMMARY. `-ksp_type cg -pc_type gamg` itself is a RUNTIME PETSc option
! (e.g. via PETSC_OPTIONS or -ksp_type/-pc_type on the command line,
! consumed by the existing KSPSetFromOptions call below) -- no code branch
! needed to select it, and the hardcoded PCCHOLESKY/MUMPS default above is
! unaffected when it is not passed (MatSetNearNullSpace is ignored by a
! direct solver; KSPSetInitialGuessNonzero is harmless for KSPPREONLY, whose
! one "iteration" is the direct solve itself).
!
! CAUGHT BY THE ROW 8C SWEEP, NOT FIXED HERE: this Mat is still
! MatCreateSeqAIJ on PETSC_COMM_SELF (rank 0 only), so this path is
! genuinely serial regardless of -np -- the bp5 sweep shows solver=2 flat at
! ~1.015-1.018 s/step for ranks 1/2/4/8 while solver=1 (MUMPS, parallel-root)
! scales 1.015 -> 0.524 s/step. Distributing the matrix (MPIAIJ on
! PETSC_COMM_WORLD, row-owned per rank) is real work, not a follow-on
! polish item -- it is the FIRST thing row 8b must do, ahead of the
! GAMG/CG/near-null-space code below, which was written and left uncompiled
! pending that decision. Row 8a's own gate (parity) does not need it: a
! serial solve is still the same solve.
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
! Matrix is built as a SeqAIJ Mat on PETSC_COMM_SELF, rank 0 only: CRS here is
! rank-0-centralized (not distributed across ranks), so MatCreateMPIAIJ would
! be the wrong tool -- see the comment above, and PATHWAY_FORWARD row 8b/8c
! for what changes when that stops being true.
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
    Vec                        :: bvec, xvec, coordsVec
    KSP                        :: ksp
    PC                         :: pc
    MatNullSpace               :: nullsp
    PetscErrorCode             :: perr
    PetscInt                   :: prow, pncols
    PetscInt, allocatable      :: pnnz(:)
    PetscInt                   :: pcols(100)
    PetscScalar                :: pvals(100)
    PetscScalar, pointer       :: parr(:)
    PetscInt                   :: kspIts, totalKSPIts, nKSPSolves

    integer (kind = 4) :: i,j,inv,jnv,ntag,node_num,var,l,k, iiTag
    real (kind = dp) :: startTime, endTime, timeUsedInFactorization,&
        timeUsedInComputing
    character (len = 50) :: netcdf_outfile, output_type

    call PetscInitialize(PETSC_NULL_CHARACTER, perr)
    CHKERRA(perr)

    call getScalarOnFaultQuant

    if (bp == 8) call pore_pressure_init
    if (bp == 8) call bp8_profile_init

    if (me.eq.0) then
        write(*,*) '= Building Stiffness Matrix in CRS format      ='
        call createMatrixHolderInCRSFormat
        call elemAssembleInCRS

        write(*,*) '= Converting from CRS to PETSc Mat/Vec format  ='
        allocate(pnnz(neq))
        do i = 1, neq
            pnnz(i) = num(i)
        enddo

        call MatCreateSeqAIJ(PETSC_COMM_SELF, neq, neq, 0, pnnz, Amat, perr)
        CHKERRA(perr)
        deallocate(pnnz)

        do i = 1, neq
            prow   = i - 1
            pncols = num(i)
            do j = 1, pncols
                pcols(j) = ja(ia(i)+j-1) - 1
                pvals(j) = kstiff(ia(i)+j-1)
            enddo
            call MatSetValues(Amat, 1, [prow], pncols, pcols(1:pncols), &
                pvals(1:pncols), INSERT_VALUES, perr)
            CHKERRA(perr)
        enddo
        call MatAssemblyBegin(Amat, MAT_FINAL_ASSEMBLY, perr)
        CHKERRA(perr)
        call MatAssemblyEnd(Amat, MAT_FINAL_ASSEMBLY, perr)
        CHKERRA(perr)
        ! SYM=1 (SPD) in MUMPS-speak -- matches mumps_par%SYM = 1 in
        ! solveTimeLoopMUMPS.f90 exactly; only the upper triangle was
        ! inserted above, matching what SYM=1 expects.
        call MatSetOption(Amat, MAT_SPD, PETSC_TRUE, perr)
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
            write(*,*) 'PETSc solver: neq is not a multiple of ndof; the ', &
                'rigid-body near-null-space assumes every free node ', &
                'contributes a full ndof-block of equations. Skipping ', &
                'MatSetNearNullSpace -- GAMG will still run (with a ', &
                'weaker default null space), just without this hint.'
        else
            call VecCreateSeq(PETSC_COMM_SELF, neq, coordsVec, perr)
            CHKERRA(perr)
            call VecSetBlockSize(coordsVec, ndof, perr)
            CHKERRA(perr)
            do i = 1, numnp
                if (id(1,i) > 0) then
                    do j = 1, ndof
                        call VecSetValue(coordsVec, id(1,i)-2+j, x(j,i), &
                            INSERT_VALUES, perr)
                        CHKERRA(perr)
                    enddo
                endif
            enddo
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

        call VecCreateSeq(PETSC_COMM_SELF, neq, bvec, perr)
        CHKERRA(perr)
        call VecDuplicate(bvec, xvec, perr)
        CHKERRA(perr)
        ! Row 8b warm start: KSPSetInitialGuessNonzero makes KSPSolve read
        ! xvec's current contents as the initial guess instead of starting
        ! from zero every step. xvec is never zeroed or overwritten between
        ! solves in the time loop below (only read via VecGetArray into
        ! `resu`), so it already carries the previous step's solution into
        ! the next KSPSolve call once this is on -- the one explicit
        ! VecZeroEntries below is only for the very first solve, so that
        ! "previous step's solution" is a defined zero rather than
        ! uninitialized memory the first time through.
        call VecZeroEntries(xvec, perr)
        CHKERRA(perr)

        call KSPCreate(PETSC_COMM_SELF, ksp, perr)
        CHKERRA(perr)
        call KSPSetOperators(ksp, Amat, Amat, perr)
        CHKERRA(perr)
        call KSPSetInitialGuessNonzero(ksp, PETSC_TRUE, perr)
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

        call initOnFaultKinematics
    endif

    call cpu_time(startTime)
    if (me.eq.0) then
        call KSPSetUp(ksp, perr)   ! Combines analysis + factorization, once.
        CHKERRA(perr)
    endif
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
                call VecGetArray(bvec, parr, perr)
                CHKERRA(perr)
                parr = right
                call VecRestoreArray(bvec, parr, perr)
                CHKERRA(perr)
            endif

            call MPI_BCAST(stoptag, 1, MPI_INT, 0, MPI_COMM_WORLD, IERR)

            if (me.eq.0) then
                call KSPSolve(ksp, bvec, xvec, perr)
                CHKERRA(perr)
                ! Row 8b: iterations/step, meaningful once -ksp_type cg is
                ! selected at runtime (KSPPREONLY's default below always
                ! reports 1). Logged even in the LU/Cholesky-default case,
                ! since a stray non-1 there would itself be a bug worth
                ! seeing.
                call KSPGetIterationNumber(ksp, kspIts, perr)
                CHKERRA(perr)
                totalKSPIts = totalKSPIts + kspIts
                nKSPSolves = nKSPSolves + 1

                call VecGetArray(xvec, parr, perr)
                CHKERRA(perr)
                resu = parr
                call VecRestoreArray(xvec, parr, perr)
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
        write(*,*) '====================================================================='
        call output_run_metadata(timeUsedInComputing, timeUsedInFactorization)

        call VecDestroy(bvec, perr)
        CHKERRA(perr)
        call VecDestroy(xvec, perr)
        CHKERRA(perr)
        call KSPDestroy(ksp, perr)
        CHKERRA(perr)
        call MatDestroy(Amat, perr)
        CHKERRA(perr)
    endif

    call PetscFinalize(perr)
    CHKERRA(perr)

end subroutine solveTimeLoopPETSc
