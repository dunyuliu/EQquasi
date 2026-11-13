
subroutine solveTimeLoopPETScDD
! Distributed-ownership time loop, par.solver = 3 (PATHWAY_FORWARD row 1,
! step (e)). Beside solveTimeLoopMUMPS / solveTimeLoopPETSc, which it does
! not touch.
!
! Each rank owns one block of a regular x-y-z decomposition of the analytic
! grid (cells split evenly per direction; no partitioner). A rank
!   - assembles the stiffness of its own elements through the local-to-global
!     mapping (MatSetValuesLocal), so the same loop can later fill a MATIS;
!   - builds its own part of the right-hand side each step (VecSetValuesLocal);
!   - reads back only the displacements of its own elements' nodes (ghosted
!     Vec, no scatter-to-zero);
!   - runs friction on the fault nodes it owns. A split-node pair belongs to
!     the block of its grid point; elements of other ranks that touch it send
!     their part of the nodal force to the owner (reverse ghost update), and
!     the owner returns the new split-node velocities (forward ghost update).
! Scalars that steer the step (dtev, max slip rate, moment rate) are
! MPI_Allreduce'd, so every rank advances time itself; nothing is broadcast
! from rank 0.
!
! Built once and reused (K is constant): the factorization, the stiffness of
! every boundary/fault element (et = 2, 3) that the RHS and the fault force
! need each step, and the fault-node mass.
!
! The solve is unchanged from solveTimeLoopPETSc: KSPPREONLY + PCCHOLESKY
! through MUMPS, so the two agree to round-off (equations are renumbered so
! each rank's are contiguous, which changes the MUMPS ordering, not the
! system).
!
! What is NOT distributed yet (stage 2 of step (e)):
!   - mesh storage: mesh4num/meshgen/allocAndInit still run on every rank and
!     hold the full mesh, numbering and fault arrays (replicated, not rank-0);
!     ownership below is derived from that replicated mesh;
!   - initial conditions are read by every rank for the whole fault;
!   - file output: owners' values are gathered to rank 0, which writes
!     disp/fault netCDF and the station/global text files;
!   - pore pressure (bp == 8) is refused here, not distributed.

#include <petsc/finclude/petscksp.h>
    use petscksp
    use globalvar, gv_mat => mat
    implicit none

    Mat                        :: Amat
    Vec                        :: bvec, xvec, xloc, fvec, floc
    KSP                        :: ksp
    PC                         :: pc
    ISLocalToGlobalMapping     :: l2g
    PetscErrorCode             :: perr
    PetscInt                   :: pnOwn, pnGhost, pneq, pnLoc, pnOwnF, pnGhostF, p6, p24
    PetscInt, allocatable      :: ghostIdx(:), ghostF(:), nnzD(:), nnzO(:), allLoc(:)
    PetscInt                   :: lidx(24)
    PetscScalar                :: vals(576)
    PetscScalar, pointer       :: parr(:)
    PetscInt                   :: kspIts, totalKSPIts, nKSPSolves, minKSPIts, maxKSPIts
    PetscInt                   :: itsAtLastProgress, solvesAtLastProgress, kspMaxIts
    PetscBool                  :: isPreonly, isCholesky
    KSPConvergedReason         :: kspReason
    PetscReal                  :: kspRtol, kspAtol, kspDtol
    character (len = 16)       :: kspTypeStr, pcTypeStr

    integer (kind = 4) :: pdim(3), ncell(3), dimsS(3), ord(3)
    integer (kind = 4), allocatable :: nodeOwner(:), lnode(:), lbase(:), pidn(:), &
        myNodes(:), myElems(:), beEl(:), beIdx(:,:), locFlt(:,:), fpid(:,:), &
        cntR(:), offR(:), cntF(:), offF(:)
    logical, allocatable :: isMaster(:), fltLocal(:,:)
    real (kind = dp), allocatable :: beEs(:,:,:), uloc(:), uprev(:), rloc(:)
    integer (kind = 4) :: nMyNodes, nMyElems, nBE, nOwnEq, nGhostNodes, nLocDof, &
        nOwnF, nGhostF, nLocF

    integer (kind = 4) :: i, j, k, l, m, n, ne, r, ix, iy, iz, ift, ie, ib, nel, &
        node_num, ntag, iiTag, itmp
    real (kind = dp) :: elemmat(5), elemx(ndof,nen), es(nee,nee), em(nee), ef(nee), &
        elemu(nee), Vol, rmax(2+ntotft), rmaxG(2+ntotft), rsum(5), rsumG(5), dtevG
    real (kind = dp) :: startTime, endTime, timeUsedInFactorization, &
        timeUsedInComputing, progressStartTime
    character (len = 50) :: netcdf_outfile, output_type
    logical :: wantProf

    call PetscInitialize(PETSC_NULL_CHARACTER, perr)
    CHKERRA(perr)

    if (bp == 8) then
        if (me == 0) write(*,*) 'par.solver=3 (distributed PETSc) does not ', &
            'distribute pore pressure yet (bp == 8). Use par.solver = 1 or 2.'
        flush(6)
        call MPI_ABORT(MPI_COMM_WORLD, 1, IERR)
    endif

    call getScalarOnFaultQuant
    call initOnFaultKinematics

    ! ---- 1. Block decomposition of the analytic grid ----------------------
    ! Cells per direction are (nxt-1, nyt-1, nzt-1). MPI_Dims_create factors
    ! nprocs into three non-increasing factors; the largest goes to the
    ! direction with the most cells.
    ncell = (/ nxt-1, nyt-1, nzt-1 /)
    dimsS = 0
    call MPI_Dims_create(nprocs, 3, dimsS, IERR)
    ord = (/ 1, 2, 3 /)
    do i = 1, 2
        do j = i+1, 3
            if (ncell(ord(j)) > ncell(ord(i))) then
                itmp = ord(i); ord(i) = ord(j); ord(j) = itmp
            endif
        enddo
    enddo
    do i = 1, 3
        pdim(ord(i)) = dimsS(i)
    enddo
    if (any(pdim > ncell)) then
        if (me == 0) write(*,*) 'par.solver=3: more blocks than cells in a ', &
            'direction. blocks (x,y,z) =', pdim, ' cells (x,y,z) =', ncell, &
            '. Use fewer MPI ranks.'
        flush(6)
        call MPI_ABORT(MPI_COMM_WORLD, 1, IERR)
    endif

    allocate(nodeOwner(numnp), isMaster(numnp), myElems(numel))
    isMaster = .false.
    do ift = 1, ntotft
        do i = 1, nftnd(ift)
            isMaster(nsmp(2,i,ift)) = .true.
        enddo
    enddo
    ! Replay meshgen's node/element numbering (ix, iz, iy; a master node
    ! follows its slave immediately). A grid point belongs to the block of
    ! the cell it is the lower corner of; an element to the block of its cell.
    n = 0
    ne = 0
    nMyElems = 0
    do ix = 1, nxt
        do iz = 1, nzt
            do iy = 1, nyt
                n = n + 1
                if (n > numnp) exit
                r = blockRank(min(ix,ncell(1)), min(iy,ncell(2)), min(iz,ncell(3)))
                nodeOwner(n) = r
                if (n < numnp) then
                    if (isMaster(n+1)) then
                        n = n + 1
                        nodeOwner(n) = r
                    endif
                endif
                if (ix >= 2 .and. iy >= 2 .and. iz >= 2) then
                    ne = ne + 1
                    if (ne <= numel) then
                        if (blockRank(ix-1, iy-1, iz-1) == me) then
                            nMyElems = nMyElems + 1
                            myElems(nMyElems) = ne
                        endif
                    endif
                endif
            enddo
        enddo
    enddo
    if (n /= numnp .or. ne /= numel) then
        write(*,*) 'par.solver=3: ownership replay disagrees with meshgen. nodes', &
            n, numnp, ' elements', ne, numel
        flush(6)
        call MPI_ABORT(MPI_COMM_WORLD, 1, IERR)
    endif

    ! ---- 2. Equation and fault-pair numbering, contiguous per rank --------
    allocate(cntR(0:nprocs-1), offR(0:nprocs-1), cntF(0:nprocs-1), offF(0:nprocs-1), &
        pidn(numnp), fpid(nftmx,ntotft), fltLocal(nftmx,ntotft))
    if (.not. allocated(fltMine)) allocate(fltMine(nftmx,ntotft))
    cntR = 0
    do n = 1, numnp
        if (id(1,n) > 0) then
            if (id(2,n) /= id(1,n)+1 .or. id(3,n) /= id(1,n)+2) then
                write(*,*) 'par.solver=3: a free node without 3 consecutive equations, node', n
                flush(6)
                call MPI_ABORT(MPI_COMM_WORLD, 1, IERR)
            endif
            cntR(nodeOwner(n)) = cntR(nodeOwner(n)) + 1
        elseif (id(2,n) > 0 .or. id(3,n) > 0) then
            write(*,*) 'par.solver=3: a partly constrained node, node', n
            flush(6)
            call MPI_ABORT(MPI_COMM_WORLD, 1, IERR)
        endif
    enddo
    if (3*sum(cntR) /= neq) then
        write(*,*) 'par.solver=3: equation count mismatch', 3*sum(cntR), neq
        flush(6)
        call MPI_ABORT(MPI_COMM_WORLD, 1, IERR)
    endif
    offR(0) = 0
    do r = 1, nprocs-1
        offR(r) = offR(r-1) + cntR(r-1)
    enddo
    nOwnEq = 3*cntR(me)
    cntR = 0
    pidn = -1
    do n = 1, numnp
        if (id(1,n) > 0) then
            r = nodeOwner(n)
            pidn(n) = 3*(offR(r) + cntR(r))
            cntR(r) = cntR(r) + 1
        endif
    enddo

    cntF = 0
    do ift = 1, ntotft
        do i = 1, nftnd(ift)
            r = nodeOwner(nsmp(1,i,ift))
            cntF(r) = cntF(r) + 1
        enddo
    enddo
    offF(0) = 0
    do r = 1, nprocs-1
        offF(r) = offF(r-1) + cntF(r-1)
    enddo
    nOwnF = cntF(me)
    cntF = 0
    fltMine = .false.
    fpid = -1
    do ift = 1, ntotft
        do i = 1, nftnd(ift)
            r = nodeOwner(nsmp(1,i,ift))
            fpid(i,ift) = offF(r) + cntF(r)
            cntF(r) = cntF(r) + 1
            fltMine(i,ift) = (r == me)
        enddo
    enddo

    ! ---- 3. Local sets: nodes of own elements, plus owned and touched fault pairs
    allocate(lnode(numnp), lbase(numnp))
    lnode = 0
    nMyNodes = 0
    do ie = 1, nMyElems
        do k = 1, nen
            n = ien(k, myElems(ie))
            if (lnode(n) == 0) then
                nMyNodes = nMyNodes + 1
                lnode(n) = nMyNodes
            endif
        enddo
    enddo
    fltLocal = .false.
    nGhostF = 0
    do ift = 1, ntotft
        do i = 1, nftnd(ift)
            if (fltMine(i,ift) .or. lnode(nsmp(1,i,ift)) > 0 .or. lnode(nsmp(2,i,ift)) > 0) then
                fltLocal(i,ift) = .true.
                if (.not. fltMine(i,ift)) nGhostF = nGhostF + 1
            endif
        enddo
    enddo
    do ift = 1, ntotft
        do i = 1, nftnd(ift)
            if (fltLocal(i,ift)) then
                do k = 1, 2
                    n = nsmp(k,i,ift)
                    if (lnode(n) == 0) then
                        nMyNodes = nMyNodes + 1
                        lnode(n) = nMyNodes
                    endif
                enddo
            endif
        enddo
    enddo
    allocate(myNodes(max(1,nMyNodes)))
    do n = 1, numnp
        if (lnode(n) > 0) myNodes(lnode(n)) = n
    enddo

    nLocF = nOwnF + nGhostF
    allocate(locFlt(2,max(1,nLocF)), ghostF(max(1,nGhostF)))
    l = 0
    do ift = 1, ntotft
        do i = 1, nftnd(ift)
            if (fltMine(i,ift)) then
                k = fpid(i,ift) - offF(me) + 1
                locFlt(1,k) = i
                locFlt(2,k) = ift
            elseif (fltLocal(i,ift)) then
                l = l + 1
                locFlt(1,nOwnF+l) = i
                locFlt(2,nOwnF+l) = ift
                ghostF(l) = fpid(i,ift)
            endif
        enddo
    enddo

    lbase = -1
    nGhostNodes = 0
    k = 0
    do i = 1, nMyNodes
        n = myNodes(i)
        if (id(1,n) > 0) then
            if (nodeOwner(n) == me) then
                k = k + 1
            else
                nGhostNodes = nGhostNodes + 1
            endif
        endif
    enddo
    if (3*k /= nOwnEq) then
        write(*,*) 'par.solver=3: rank', me, 'owns free nodes outside its own elements', 3*k, nOwnEq
        flush(6)
        call MPI_ABORT(MPI_COMM_WORLD, 1, IERR)
    endif
    allocate(ghostIdx(max(1,3*nGhostNodes)))
    l = 0
    do i = 1, nMyNodes
        n = myNodes(i)
        if (id(1,n) > 0) then
            if (nodeOwner(n) == me) then
                lbase(n) = pidn(n) - 3*offR(me)
            else
                lbase(n) = nOwnEq + 3*l
                do j = 1, ndof
                    ghostIdx(3*l+j) = pidn(n) + j - 1
                enddo
                l = l + 1
            endif
        endif
    enddo
    nLocDof = nOwnEq + 3*nGhostNodes
    allocate(uloc(max(1,nLocDof)), uprev(max(1,nLocDof)), rloc(max(1,nLocDof)), &
        allLoc(max(1,nLocDof)))
    uloc = 0.0d0
    uprev = 0.0d0
    do i = 1, nLocDof
        allLoc(i) = i - 1
    enddo

    nBE = 0
    do ie = 1, nMyElems
        if (et(myElems(ie)) == 2 .or. et(myElems(ie)) == 3) nBE = nBE + 1
    enddo
    allocate(beEl(max(1,nBE)), beIdx(nee,max(1,nBE)), beEs(nee,nee,max(1,nBE)))

    ! ---- 4. PETSc objects -------------------------------------------------
    pnOwn = nOwnEq
    pnGhost = 3*nGhostNodes
    pneq = neq
    pnLoc = nLocDof
    pnOwnF = 6*nOwnF
    pnGhostF = nGhostF
    p6 = 6
    p24 = nee
    call VecCreateGhost(PETSC_COMM_WORLD, pnOwn, pneq, pnGhost, ghostIdx, xvec, perr)
    CHKERRA(perr)
    call VecDuplicate(xvec, bvec, perr)
    CHKERRA(perr)
    call VecGetLocalToGlobalMapping(xvec, l2g, perr)
    CHKERRA(perr)
    call VecSetLocalToGlobalMapping(bvec, l2g, perr)
    CHKERRA(perr)
    call VecZeroEntries(xvec, perr)
    CHKERRA(perr)
    ! Fault-pair exchange vector: 6 values per pair (slave xyz, master xyz).
    call VecCreateGhostBlock(PETSC_COMM_WORLD, p6, pnOwnF, PETSC_DECIDE, pnGhostF, &
        ghostF, fvec, perr)
    CHKERRA(perr)

    call MatCreate(PETSC_COMM_WORLD, Amat, perr)
    CHKERRA(perr)
    call MatSetSizes(Amat, pnOwn, pnOwn, pneq, pneq, perr)
    CHKERRA(perr)
    call MatSetType(Amat, MATAIJ, perr)
    CHKERRA(perr)
    ! A hex8 row couples to at most 27 nodes = 81 dofs. Uniform bound, not
    ! exact; an insufficient bound fails loudly in MatSetValuesLocal.
    allocate(nnzD(max(1,nOwnEq)), nnzO(max(1,nOwnEq)))
    nnzD = min(81, nOwnEq)
    nnzO = min(81, neq - nOwnEq)
    call MatSeqAIJSetPreallocation(Amat, 0, nnzD, perr)
    CHKERRA(perr)
    call MatMPIAIJSetPreallocation(Amat, 0, nnzD, 0, nnzO, perr)
    CHKERRA(perr)
    deallocate(nnzD, nnzO)
    call MatSetLocalToGlobalMapping(Amat, l2g, l2g, perr)
    CHKERRA(perr)

    ! ---- 5. Assembly of own elements; cache what each step reuses ----------
    if (me == 0) write(*,*) '= Assembling own elements on every rank        ='
    do i = 1, nMyNodes
        consm(:,myNodes(i)) = 0.0d0
    enddo
    ib = 0
    do ie = 1, nMyElems
        nel = myElems(ie)
        do j = 1, nen
            node_num = ien(j,nel)
            do i = 1, ndof
                elemx(i,j) = x(i,node_num)
            enddo
        enddo
        do j = 1, 5
            elemmat(j) = gv_mat(nel,j)
        enddo
        call c8g2(elemx, elemmat, es, em, ef, Vol)
        ntag = 0
        do i = 1, nen
            node_num = ien(i,nel)
            do j = 1, ndof
                ntag = ntag + 1
                if (id(j,node_num) > 0) then
                    lidx(ntag) = lbase(node_num) + j - 1
                else
                    lidx(ntag) = -1
                endif
            enddo
        enddo
        ! Both triangles from the upper one, as the CRS path stores it.
        do i = 1, nee
            do j = 1, nee
                vals((i-1)*nee+j) = es(min(i,j), max(i,j))
            enddo
        enddo
        call MatSetValuesLocal(Amat, p24, lidx, p24, lidx, vals, ADD_VALUES, perr)
        CHKERRA(perr)
        if (et(nel) == 2 .or. et(nel) == 3) then
            ib = ib + 1
            beEl(ib) = nel
            beEs(:,:,ib) = es
            do i = 1, nee
                beIdx(i,ib) = lidx(i) + 1   ! 1-based local dof, 0 if constrained
            enddo
            if (et(nel) == 2) then
                do i = 1, nen
                    node_num = ien(i,nel)
                    do j = 1, ndof
                        if (id(j,node_num) < 0) consm(j,node_num) = consm(j,node_num) + em((i-1)*3+j)
                    enddo
                enddo
            endif
        endif
    enddo
    call MatAssemblyBegin(Amat, MAT_FINAL_ASSEMBLY, perr)
    CHKERRA(perr)
    call MatAssemblyEnd(Amat, MAT_FINAL_ASSEMBLY, perr)
    CHKERRA(perr)
    call MatSetOption(Amat, MAT_SPD, PETSC_TRUE, perr)
    CHKERRA(perr)
    call MatSetOption(Amat, MAT_SYMMETRIC, PETSC_TRUE, perr)
    CHKERRA(perr)
    ! Fault-node mass is constant: sum the per-rank parts at the owner, once.
    call ddFltReverse(consm)

    call KSPCreate(PETSC_COMM_WORLD, ksp, perr)
    CHKERRA(perr)
    call KSPSetOperators(ksp, Amat, Amat, perr)
    CHKERRA(perr)
    call KSPSetType(ksp, KSPPREONLY, perr)
    CHKERRA(perr)
    call KSPGetPC(ksp, pc, perr)
    CHKERRA(perr)
    call PCSetType(pc, PCCHOLESKY, perr)
    CHKERRA(perr)
    call PCFactorSetMatSolverType(pc, MATSOLVERMUMPS, perr)
    CHKERRA(perr)
    call KSPSetFromOptions(ksp, perr)
    CHKERRA(perr)
    call PetscObjectTypeCompare(ksp, KSPPREONLY, isPreonly, perr)
    CHKERRA(perr)
    call PetscObjectTypeCompare(pc, PCCHOLESKY, isCholesky, perr)
    CHKERRA(perr)
    call KSPGetTolerances(ksp, kspRtol, kspAtol, kspDtol, kspMaxIts, perr)
    CHKERRA(perr)
    if (.not. (isPreonly .and. isCholesky)) then
        if (me == 0) write(*,*) 'par.solver=3: only KSPPREONLY + PCCHOLESKY ', &
            '(MUMPS) is verified for the distributed loop. Refusing.'
        flush(6)
        call MPI_ABORT(MPI_COMM_WORLD, 1, IERR)
    endif
    kspTypeStr = 'preonly'
    pcTypeStr = 'cholesky'

    call cpu_time(startTime)
    call KSPSetUp(ksp, perr)
    CHKERRA(perr)
    call cpu_time(endTime)
    timeUsedInFactorization = endTime - startTime

    stoptag = 0
    totalKSPIts = 0
    nKSPSolves = 0
    minKSPIts = huge(minKSPIts)
    maxKSPIts = 0
    itsAtLastProgress = 0
    solvesAtLastProgress = 0
    ddActive = .true.

    call cpu_time(startTime)
    progressStartTime = startTime
    do it = 1, nstep
        if (stoptag == 1) exit

        uprev = uloc

        if (bp == 7 .and. icstart == 1) then
            if (it < (nuct/dt)) then
                dtev1 = dt
            endif
        endif
        if (dtmax > 0.0d0) dtev1 = min(dtev1, dtmax)
        time = time + dtev1

        if (mod(it,nhplt) == 1 .and. me == 0) then
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

        ! Off-fault stations: recorded by the rank that owns the node; the
        ! rest stay zero and are summed onto rank 0 after the loop.
        if (allocated(dout)) then
            if (me == 0) dout(1,it) = time
            do i = 1, ndout
                j = idhist(1,i)
                if (j <= 0) j = 1
                if (nodeOwner(j) == me) then
                    k = idhist(2,i)
                    l = idhist(3,i)
                    if (l == 1) then
                        dout(i+1,it) = cons(k,j)
                    elseif (l == 2) then
                        dout(i+1,it) = consv(k,j)
                    elseif (l == 3) then
                        dout(i+1,it) = consa(k,j)
                    endif
                endif
            enddo
        endif

        do iiTag = 0, 1
            itag = iiTag
            call ddBoundLoad
            call VecZeroEntries(bvec, perr)
            CHKERRA(perr)
            call VecSetValuesLocal(bvec, pnLoc, allLoc, rloc, ADD_VALUES, perr)
            CHKERRA(perr)
            call VecAssemblyBegin(bvec, perr)
            CHKERRA(perr)
            call VecAssemblyEnd(bvec, perr)
            CHKERRA(perr)

            call KSPSolve(ksp, bvec, xvec, perr)
            CHKERRA(perr)
            call KSPGetConvergedReason(ksp, kspReason, perr)
            CHKERRA(perr)
            if (kspReason%v < 0) then
                if (me == 0) write(*,*) 'par.solver=3: KSPSolve did NOT converge ', &
                    '(KSPConvergedReason =', kspReason%v, ') at step', it, 'iiTag', iiTag
                flush(6)
                call MPI_ABORT(MPI_COMM_WORLD, 1, IERR)
            endif
            call KSPGetIterationNumber(ksp, kspIts, perr)
            CHKERRA(perr)
            totalKSPIts = totalKSPIts + kspIts
            nKSPSolves = nKSPSolves + 1
            minKSPIts = min(minKSPIts, kspIts)
            maxKSPIts = max(maxKSPIts, kspIts)

            ! Displacements of own elements' nodes: owned + ghost values.
            call VecGhostUpdateBegin(xvec, INSERT_VALUES, SCATTER_FORWARD, perr)
            CHKERRA(perr)
            call VecGhostUpdateEnd(xvec, INSERT_VALUES, SCATTER_FORWARD, perr)
            CHKERRA(perr)
            call VecGhostGetLocalForm(xvec, xloc, perr)
            CHKERRA(perr)
            call VecGetArray(xloc, parr, perr)
            CHKERRA(perr)
            do i = 1, nLocDof
                uloc(i) = parr(i)
            enddo
            call VecRestoreArray(xloc, parr, perr)
            CHKERRA(perr)
            call VecGhostRestoreLocalForm(xvec, xloc, perr)
            CHKERRA(perr)

            if (iiTag == 0) then
                do i = 1, nMyNodes
                    n = myNodes(i)
                    if (id(1,n) > 0) then
                        do j = 1, ndof
                            constmp(j,n) = uloc(lbase(n)+j)
                        enddo
                    endif
                enddo
            else
                do i = 1, nMyNodes
                    n = myNodes(i)
                    if (id(1,n) > 0) then
                        do j = 1, ndof
                            cons(j,n) = uloc(lbase(n)+j)
                        enddo
                    else
                        do j = 1, ndof
                            cons(j,n) = constmp(j,n)
                        enddo
                    endif
                enddo
            endif

            call ddBoundFtKu
            call faulting
            ! Owner -> ghost copies of the split-node velocities the next
            ! RHS needs.
            call ddFltForward(consv)
            call ddFltForward(consvtmp)

            call MPI_Allreduce(dtev, dtevG, 1, MPI_DOUBLE_PRECISION, MPI_MIN, &
                MPI_COMM_WORLD, IERR)
            dtev = dtevG
            if (itag == 1) then
                rmax(1) = pma
                rmax(2) = maxSlipRate
                rmax(3:2+ntotft) = peakSlipRatePerFault(1:ntotft,it)
                call MPI_Allreduce(rmax, rmaxG, 2+ntotft, MPI_DOUBLE_PRECISION, MPI_MAX, &
                    MPI_COMM_WORLD, IERR)
                pma = rmaxG(1)
                maxSlipRate = rmaxG(2)
                peakSlipRatePerFault(1:ntotft,it) = rmaxG(3:2+ntotft)
                rsum = (/ totMomRate, totMomRateVW, totRuptArea, totTaoRuptArea, totSlipRuptArea /)
                call MPI_Allreduce(rsum, rsumG, 5, MPI_DOUBLE_PRECISION, MPI_SUM, &
                    MPI_COMM_WORLD, IERR)
                totMomRate = rsumG(1)
                totMomRateVW = rsumG(2)
                totRuptArea = rsumG(3)
                totTaoRuptArea = rsumG(4)
                totSlipRuptArea = rsumG(5)
            endif
        enddo

        ! bp5 profile file (output_prof, same schedule as inside faulting).
        wantProf = (bp == 5) .and. ((status1 == 0 .and. mod(it,20) == 1) .or. &
            (status1 == 1 .and. mod(it,10) == 0))
        if (wantProf) then
            call ddGatherFault
            if (me == 0) call output_prof
        endif

        status0 = status1
        dtev1 = max(dt,int8(dtev/dt)*dt)

        do i = 1, nMyNodes
            n = myNodes(i)
            if (id(1,n) > 0) then
                do j = 1, ndof
                    consv(j,n) = (uloc(lbase(n)+j) - uprev(lbase(n)+j))/dtev1
                enddo
            endif
        enddo
        globaldat(1,it) = time
        globaldat(2,it) = maxSlipRate
        globaldat(3,it) = totMomRate
        globaldat(4,it) = totTaoRuptArea
        globaldat(5,it) = totSlipRuptArea
        globaldat(6,it) = totRuptArea
        globaldat(7,it) = totMomRateVW

        if (mod(it,nt_output_stress) == 1 .or. (stoptag == 1) .or. (it==nstep)) then
            call ddGatherFault
            call ddGatherNodal(cons)
            if (me == 0) then
                write(proc_str,'(I5.5)') it
                netcdf_outfile = trim(outDir)//'disp.'//trim(proc_str)//'.nc'
                output_type = 'disp'
                call netcdf_write(netcdf_outfile, output_type)
                netcdf_outfile = trim(outDir)//'fault.'//trim(proc_str)//'.nc'
                call netcdf_write_on_fault(netcdf_outfile)
                if ((stoptag == 1) .or. (it==nstep)) then
                    netcdf_outfile = trim(outDir)//'disp.r.nc'
                    output_type = 'disp'
                    call netcdf_write(netcdf_outfile, output_type)
                    netcdf_outfile = trim(outDir)//'fault.r.nc'
                    call netcdf_write_on_fault(netcdf_outfile)
                endif
            endif
        endif

        if (mod(it, 100) == 0 .and. me == 0) then
            call cpu_time(endTime)
            write(*,'(A,I8,A,E12.5,A,E12.5,A,E12.5,A,F10.2,A,F10.4,A,F8.2)') &
                'PROGRESS step=', it, ' t=', time, ' dt=', dtev1, &
                ' Vmax=', maxSlipRate, ' wall=', endTime - progressStartTime, &
                ' s/step=', (endTime - progressStartTime) / 100.0d0, &
                ' ksp_its=', dble(totalKSPIts - itsAtLastProgress) / &
                    dble(max(1, nKSPSolves - solvesAtLastProgress))
            progressStartTime = endTime
            itsAtLastProgress = totalKSPIts
            solvesAtLastProgress = nKSPSolves
        endif
    enddo
    call cpu_time(endTime)
    timeUsedInComputing = endTime-startTime

    ! Everything writeResults (rank 0) reads: owners' values onto rank 0.
    call ddGatherFault
    if (allocated(fltsta)) call ddReduceFlat(fltsta, size(fltsta))
    if (allocated(dout)) call ddReduceFlat(dout, size(dout))
    ddActive = .false.

    if (me == 0) then
        write(*,*) '====================================================================='
        write(*,*) '=                        RUN SUMMARY                                ='
        write(*,'(X,A,40X,i7,4X,A)')     '= MPI ranks                = ', nprocs, '='
        write(*,'(X,A,28X,3i6,5X,A)')    '= Blocks (x, y, z)         = ', pdim, '='
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
        if (it-1 < nstep) then
            call output_run_metadata(timeUsedInComputing, timeUsedInFactorization, &
                kspTypeStr, pcTypeStr, 'stoptag', minKSPIts, &
                dble(totalKSPIts)/dble(max(1,nKSPSolves)), maxKSPIts)
        else
            call output_run_metadata(timeUsedInComputing, timeUsedInFactorization, &
                kspTypeStr, pcTypeStr, 'nstep_reached', minKSPIts, &
                dble(totalKSPIts)/dble(max(1,nKSPSolves)), maxKSPIts)
            if (eqquasi_mode == 2) then
                write(*,*) 'WARNING: eqquasi_mode 2 reached nstep without a rupture'
                write(*,*) '         (no slip rate above exit_slip_rate); raise par.nstep.'
            endif
        endif
    endif

    call VecDestroy(fvec, perr)
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

contains

    integer (kind = 4) function blockRank(cx, cy, cz)
        ! Rank of the block that holds cell (cx, cy, cz), 1-based cell indices.
        integer (kind = 4), intent(in) :: cx, cy, cz
        blockRank = ((cx-1)*pdim(1))/ncell(1) + pdim(1)*(((cy-1)*pdim(2))/ncell(2) &
            + pdim(2)*(((cz-1)*pdim(3))/ncell(3)))
    end function blockRank

    subroutine ddBoundLoad
        ! bound_load on this rank: prescribed displacements of the local
        ! constrained nodes, then this rank's part of the RHS into rloc.
        integer (kind = 4) :: ii, jj, nn, code, b
        rloc = 0.0d0
        do ii = 1, nMyNodes
            nn = myNodes(ii)
            if (id(1,nn) > 0) then
                constmp(:,nn) = 0.0d0
                cycle
            endif
            do jj = 1, ndof
                code = id(jj,nn)
                if (code == -2) then
                    consv(jj,nn) = -far_load_rate
                elseif (code == -3) then
                    consv(jj,nn) = far_load_rate
                elseif (code == -21) then
                    consv(jj,nn) = -far_norm_load_vel
                elseif (code == -31) then
                    consv(jj,nn) = far_norm_load_vel
                elseif (code == -1) then
                    consv(jj,nn) = 0.0d0
                endif
                consa(jj,nn) = 0.0d0
                if (itag == 0) then
                    consvtmp(jj,nn) = consv(jj,nn)
                    constmp(jj,nn) = cons(jj,nn) + consv(jj,nn)*dtev1
                else
                    constmp(jj,nn) = cons(jj,nn) + 0.5d0*(consv(jj,nn)+consvtmp(jj,nn))*dtev1
                endif
            enddo
        enddo
        do b = 1, nBE
            ntag = 0
            do ii = 1, nen
                nn = ien(ii,beEl(b))
                do jj = 1, ndof
                    ntag = ntag + 1
                    elemu(ntag) = constmp(jj,nn)
                enddo
            enddo
            do jj = 1, nee
                if (beIdx(jj,b) == 0) then
                    do ii = 1, nee
                        if (beIdx(ii,b) > 0) then
                            rloc(beIdx(ii,b)) = rloc(beIdx(ii,b)) - beEs(ii,jj,b)*elemu(jj)
                        endif
                    enddo
                endif
            enddo
        enddo
    end subroutine ddBoundLoad

    subroutine ddBoundFtKu
        ! bound_ft_ku on this rank: -K u at the constrained nodes of own
        ! fault elements; the parts are then summed at each pair's owner.
        integer (kind = 4) :: ii, jj, mm, nn, b
        do ii = 1, nMyNodes
            consf(:,myNodes(ii)) = 0.0d0
        enddo
        do b = 1, nBE
            if (et(beEl(b)) /= 2) cycle
            ntag = 0
            do ii = 1, nen
                nn = ien(ii,beEl(b))
                do jj = 1, ndof
                    ntag = ntag + 1
                    if (itag == 0) then
                        elemu(ntag) = constmp(jj,nn)
                    else
                        elemu(ntag) = cons(jj,nn)
                    endif
                enddo
            enddo
            do ii = 1, nen
                nn = ien(ii,beEl(b))
                do jj = 1, ndof
                    if (id(jj,nn) < 0) then
                        do mm = 1, nee
                            consf(jj,nn) = consf(jj,nn) - beEs((ii-1)*3+jj,mm,b)*elemu(mm)
                        enddo
                    endif
                enddo
            enddo
        enddo
        call ddFltReverse(consf)
    end subroutine ddBoundFtKu

    subroutine ddFltPack(a)
        real (kind = dp), intent(in) :: a(ndof,numnp)
        integer (kind = 4) :: ll, kk
        call VecGhostGetLocalForm(fvec, floc, perr)
        CHKERRA(perr)
        call VecGetArray(floc, parr, perr)
        CHKERRA(perr)
        do ll = 1, nLocF
            do kk = 1, 3
                parr(6*(ll-1)+kk)   = a(kk, nsmp(1,locFlt(1,ll),locFlt(2,ll)))
                parr(6*(ll-1)+3+kk) = a(kk, nsmp(2,locFlt(1,ll),locFlt(2,ll)))
            enddo
        enddo
        call VecRestoreArray(floc, parr, perr)
        CHKERRA(perr)
        call VecGhostRestoreLocalForm(fvec, floc, perr)
        CHKERRA(perr)
    end subroutine ddFltPack

    subroutine ddFltUnpack(a, l1, l2)
        real (kind = dp), intent(inout) :: a(ndof,numnp)
        integer (kind = 4), intent(in) :: l1, l2
        integer (kind = 4) :: ll, kk
        call VecGhostGetLocalForm(fvec, floc, perr)
        CHKERRA(perr)
        call VecGetArray(floc, parr, perr)
        CHKERRA(perr)
        do ll = l1, l2
            do kk = 1, 3
                a(kk, nsmp(1,locFlt(1,ll),locFlt(2,ll))) = parr(6*(ll-1)+kk)
                a(kk, nsmp(2,locFlt(1,ll),locFlt(2,ll))) = parr(6*(ll-1)+3+kk)
            enddo
        enddo
        call VecRestoreArray(floc, parr, perr)
        CHKERRA(perr)
        call VecGhostRestoreLocalForm(fvec, floc, perr)
        CHKERRA(perr)
    end subroutine ddFltUnpack

    subroutine ddFltReverse(a)
        ! Sum every rank's part of a fault-pair quantity at the pair's owner.
        real (kind = dp), intent(inout) :: a(ndof,numnp)
        call ddFltPack(a)
        call VecGhostUpdateBegin(fvec, ADD_VALUES, SCATTER_REVERSE, perr)
        CHKERRA(perr)
        call VecGhostUpdateEnd(fvec, ADD_VALUES, SCATTER_REVERSE, perr)
        CHKERRA(perr)
        call ddFltUnpack(a, 1, nOwnF)
    end subroutine ddFltReverse

    subroutine ddFltForward(a)
        ! Copy a fault-pair quantity from its owner to the ranks that touch it.
        real (kind = dp), intent(inout) :: a(ndof,numnp)
        call ddFltPack(a)
        call VecGhostUpdateBegin(fvec, INSERT_VALUES, SCATTER_FORWARD, perr)
        CHKERRA(perr)
        call VecGhostUpdateEnd(fvec, INSERT_VALUES, SCATTER_FORWARD, perr)
        CHKERRA(perr)
        call ddFltUnpack(a, nOwnF+1, nLocF)
    end subroutine ddFltForward

    subroutine ddReduceFlat(a, nn)
        ! Sum a over ranks onto rank 0 (entries are zero off their owner).
        integer (kind = 4), intent(in) :: nn
        real (kind = dp), intent(inout) :: a(nn)
        real (kind = dp), allocatable :: g(:)
        allocate(g(nn))
        call MPI_Reduce(a, g, nn, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
        if (me == 0) a = g
        deallocate(g)
    end subroutine ddReduceFlat

    subroutine ddGatherFault
        ! Owners' on-fault state onto rank 0, for output. A rank never reads
        ! the pairs it does not own, so its copies of them are zeroed here.
        integer (kind = 4) :: ii, jft
        do jft = 1, ntotft
            do ii = 1, nftmx
                if (ii > nftnd(jft) .or. .not. fltMine(ii,jft)) then
                    fric(:,ii,jft) = 0.0d0
                    fnft(ii,jft) = 0.0d0
                    slp4fri(ii,jft) = 0.0d0
                endif
            enddo
        enddo
        call ddReduceFlat(fric, size(fric))
        call ddReduceFlat(fnft, size(fnft))
        call ddReduceFlat(slp4fri, size(slp4fri))
    end subroutine ddGatherFault

    subroutine ddGatherNodal(a)
        ! Owners' nodal values onto rank 0, for output.
        real (kind = dp), intent(inout) :: a(ndof,numnp)
        real (kind = dp), allocatable :: s(:,:)
        integer (kind = 4) :: nn
        allocate(s(ndof,numnp))
        s = 0.0d0
        do nn = 1, numnp
            if (nodeOwner(nn) == me) s(:,nn) = a(:,nn)
        enddo
        call ddReduceFlat(s, size(s))
        if (me == 0) a = s
        deallocate(s)
    end subroutine ddGatherNodal

end subroutine solveTimeLoopPETScDD
