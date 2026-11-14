
subroutine solveTimeLoopPETScDD
! Distributed-ownership time loop, par.solver = 3 (PATHWAY_FORWARD row 1,
! step (e)). Beside solveTimeLoopMUMPS / solveTimeLoopPETSc, which it does
! not touch.
!
! Each rank owns one block of a regular x-y-z decomposition of the analytic
! grid (cells split evenly per direction; no partitioner) and holds only that
! block's mesh: ddMeshBox (below) sets the block before allocAndInit, meshgen
! then generates the block's own nodes and elements, and the initial
! conditions are read for the fault nodes the rank owns. A rank
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
! Output, the same files as solvers 1 and 2:
!   - a station file (fltst_, srfst_) is written by the rank that owns the
!     station (writeResults);
!   - the whole-fault text files (cplot_*, p1output, p2output) have lines of
!     one fixed length, so every rank writes its own nodes' lines at their
!     byte offsets (MPI-IO, ddTextRecords);
!   - the netCDF files (mesh, disp, fault) are netCDF classic and written by
!     rank 0, serially: every rank sends only its own nodes', elements' or
!     fault nodes' values, and rank 0 writes one rank's part at a time
!     (ddColsToFile, netcdf_write_on_fault), so it never holds a whole field;
!   - global.dat and the other step series are scalars every rank already
!     holds (MPI_Allreduce above); rank 0 writes them.
!
! Numbering. A grid point is owned by the block of the cell it is the lower
! corner of. Equations are contiguous per rank, in the rank's own node order;
! a rank asks the owners for the equation numbers of its ghost nodes, once.
! Whole-mesh node and element numbers (meshgen's, for the output files) are
! closed-form in the grid indices.
!
! Built once and reused (K is constant): the factorization, the stiffness of
! every boundary/fault element (et = 2, 3) that the RHS and the fault force
! need each step, the fault-node mass and the fault-node area.
!
! The solve is unchanged from solveTimeLoopPETSc: KSPPREONLY + PCCHOLESKY
! through MUMPS, so the two agree to round-off (equations are renumbered so
! each rank's are contiguous, which changes the MUMPS ordering, not the
! system).
!
! What is NOT distributed yet:
!   - the fault arrays (fric, nsmp, un, ...) keep the whole fault's index
!     range on every rank; a rank fills only its own block's entries;
!   - the netCDF files pass through rank 0 (above): no whole field in
!     memory, but one writer;
!   - mesh4num still counts the whole grid on every rank (a loop, no arrays);
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

    integer (kind = 4) :: pdim(3), ncell(3), icnt(6), icntG(6), nloc(3), profN(2)
    integer (kind = 4), allocatable :: nodeAt(:,:,:), nodeOwner(:), gnode(:), lbase(:), &
        pidn(:), beEl(:), beIdx(:,:), locFlt(:,:), fpid(:,:), cntR(:), offR(:), &
        cntF(:), offF(:), scnt(:), sdsp(:), rcnt(:), rdsp(:), sbuf(:), rbuf(:), &
        reqNode(:), stnNode(:), stnOwn(:), stnOwnG(:), profI(:,:), profP(:,:)
    logical, allocatable :: fltLocal(:,:)
    real (kind = dp), allocatable :: beEs(:,:,:), uloc(:), uprev(:), rloc(:)
    integer (kind = 4) :: nBE, nOwnEq, nGhostNodes, nLocDof, nOwnF, nGhostF, nLocF, &
        nOwnNodes, nOwnFree

    integer (kind = 4) :: i, j, k, l, n, r, ix, iy, iz, ift, ib, nel, &
        node_num, ntag, iiTag
    real (kind = dp) :: elemmat(5), elemx(ndof,nen), es(nee,nee), em(nee), ef(nee), &
        elemu(nee), Vol, rmax(2+ntotft), rmaxG(2+ntotft), rsum(5), rsumG(5), dtevG, &
        hwm(2), hwmG(2)
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
    if (.not. ddActive) then
        if (me == 0) write(*,*) 'par.solver=3: the per-rank mesh was not set up ', &
            '(ddMeshBox did not run).'
        flush(6)
        call MPI_ABORT(MPI_COMM_WORLD, 1, IERR)
    endif

    ! Peak resident memory of this rank with its mesh and initial conditions
    ! in place, before anything this routine allocates.
    hwm(1) = dble(ddPeakRSSkB())/1024.0d0
    hwm(2) = hwm(1)

    ! ---- 1. This block's nodes: owner and whole-mesh number ----------------
    ! Replays meshgen's numbering over the block (ix, iz, iy; a master node
    ! follows its slave immediately).
    call ddBlocks(pdim, ncell)
    if (dble(nxt)*dble(nyt)*dble(nzt) + 2.0d0*sum(nftnd(1:ntotft)) > dble(huge(1))) then
        if (me == 0) write(*,*) 'par.solver=3: more grid points than a 4-byte node number holds.'
        flush(6)
        call MPI_ABORT(MPI_COMM_WORLD, 1, IERR)
    endif
    allocate(nodeAt(ddBox(1,1):ddBox(2,1), ddBox(1,2):ddBox(2,2), ddBox(1,3):ddBox(2,3)), &
        nodeOwner(numnp), gnode(numnp))
    nodeOwner = -1
    gnode = 0
    n = 0
    nOwnNodes = 0
    do ix = ddBox(1,1), ddBox(2,1)
        do iz = ddBox(1,3), ddBox(2,3)
            do iy = ddBox(1,2), ddBox(2,2)
                r = blockRank(min(ix,ncell(1)), min(iy,ncell(2)), min(iz,ncell(3)))
                l = 1
                if (faultAt(ix, iy, iz) > 0) l = 2
                do k = 1, l
                    n = n + 1
                    nodeOwner(n) = r
                    gnode(n) = ((ix-1)*nzt + (iz-1))*nyt + iy + splitsBefore(ix, iy, iz) + k - 1
                    if (r == me) nOwnNodes = nOwnNodes + 1
                enddo
                nodeAt(ix,iy,iz) = n - l + 1
            enddo
        enddo
    enddo
    if (n /= numnp) then
        write(*,*) 'par.solver=3: rank', me, 'numbering replay disagrees with meshgen', n, numnp
        flush(6)
        call MPI_ABORT(MPI_COMM_WORLD, 1, IERR)
    endif

    ! Every node and element of the whole mesh on exactly one rank; element
    ! types summed over ranks (report_element_types leaves this to us).
    icnt = (/ nOwnNodes, numel, count(et(1:numel) == 1), count(et(1:numel) == 2), &
        count(et(1:numel) == 3), 0 /)
    call MPI_Allreduce(icnt, icntG, 6, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, IERR)
    if (icntG(1) /= ddGlobal(1) .or. icntG(2) /= ddGlobal(2)) then
        if (me == 0) write(*,*) 'par.solver=3: the blocks do not add up to the mesh. nodes', &
            icntG(1), ddGlobal(1), ' elements', icntG(2), ddGlobal(2)
        flush(6)
        call MPI_ABORT(MPI_COMM_WORLD, 1, IERR)
    endif
    nloc = (/ numnp, numel, neq /)
    call MPI_Allreduce(nloc, icnt, 3, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, IERR)
    call MPI_Allreduce(hwm(1), hwmG(1), 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, IERR)
    call MPI_Allreduce(hwm(2), hwmG(2), 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, IERR)
    if (me == 0) then
        write(*,'(X,A,40X,i7,4X,A)') '= Interior elements (et=1)  = ', icntG(3), '='
        write(*,'(X,A,40X,i7,4X,A)') '= Fault elements    (et=2)  = ', icntG(4), '='
        write(*,'(X,A,40X,i7,4X,A)') '= Far-field elements(et=3)  = ', icntG(5), '='
        write(*,'(X,A,2i12)')   '= Nodes on the largest rank, in the mesh    = ', icnt(1), ddGlobal(1)
        write(*,'(X,A,2i12)')   '= Elements on the largest rank, in the mesh = ', icnt(2), ddGlobal(2)
        ! No line where /proc/self/status is not there (ddPeakRSSkB = -1).
        if (hwmG(1) > 0.0d0) write(*,'(X,A,2f12.1)') &
            '= Peak RSS after mesh setup, MB: largest rank, sum = ', hwmG(1), hwmG(2)
    endif
    if (icntG(4) == 0 .and. any(nftnd(1:ntotft) > 0)) then
        if (me == 0) write(*,*) 'par.solver=3: no fault-adjacent elements (et=2) were tagged, ', &
            'so every fault node would carry zero mass. Reduce dx or widen the model along y.'
        flush(6)
        call MPI_ABORT(MPI_COMM_WORLD, 3, IERR)
    endif

    ! ---- 2. Equation and fault-pair numbering, contiguous per rank --------
    allocate(cntR(0:nprocs-1), offR(0:nprocs-1), cntF(0:nprocs-1), offF(0:nprocs-1), &
        pidn(numnp), lbase(numnp), fpid(nftmx,ntotft), fltLocal(nftmx,ntotft), &
        scnt(0:nprocs-1), sdsp(0:nprocs-1), rcnt(0:nprocs-1), rdsp(0:nprocs-1))
    nOwnFree = 0
    scnt = 0
    do n = 1, numnp
        if (id(1,n) > 0) then
            if (id(2,n) /= id(1,n)+1 .or. id(3,n) /= id(1,n)+2) then
                write(*,*) 'par.solver=3: a free node without 3 consecutive equations, node', n
                flush(6)
                call MPI_ABORT(MPI_COMM_WORLD, 1, IERR)
            endif
            if (nodeOwner(n) == me) then
                nOwnFree = nOwnFree + 1
            else
                scnt(nodeOwner(n)) = scnt(nodeOwner(n)) + 1
            endif
        elseif (id(2,n) > 0 .or. id(3,n) > 0) then
            write(*,*) 'par.solver=3: a partly constrained node, node', n
            flush(6)
            call MPI_ABORT(MPI_COMM_WORLD, 1, IERR)
        endif
    enddo
    call MPI_Allgather(nOwnFree, 1, MPI_INTEGER, cntR, 1, MPI_INTEGER, MPI_COMM_WORLD, IERR)
    if (3*sum(cntR) /= ddGlobal(3)) then
        if (me == 0) write(*,*) 'par.solver=3: equation count mismatch', 3*sum(cntR), ddGlobal(3)
        flush(6)
        call MPI_ABORT(MPI_COMM_WORLD, 1, IERR)
    endif
    offR(0) = 0
    do r = 1, nprocs-1
        offR(r) = offR(r-1) + cntR(r-1)
    enddo
    nOwnEq = 3*nOwnFree
    pidn = -1
    k = 0
    do n = 1, numnp
        if (id(1,n) > 0 .and. nodeOwner(n) == me) then
            pidn(n) = 3*(offR(me) + k)
            k = k + 1
        endif
    enddo

    ! Ghost nodes: ask each owner, by grid point, for the equation number.
    nGhostNodes = sum(scnt)
    call MPI_Alltoall(scnt, 1, MPI_INTEGER, rcnt, 1, MPI_INTEGER, MPI_COMM_WORLD, IERR)
    sdsp(0) = 0
    rdsp(0) = 0
    do r = 1, nprocs-1
        sdsp(r) = sdsp(r-1) + scnt(r-1)
        rdsp(r) = rdsp(r-1) + rcnt(r-1)
    enddo
    allocate(sbuf(max(1,nGhostNodes)), reqNode(max(1,nGhostNodes)), rbuf(max(1,sum(rcnt))))
    scnt = 0
    do ix = ddBox(1,1), ddBox(2,1)
        do iz = ddBox(1,3), ddBox(2,3)
            do iy = ddBox(1,2), ddBox(2,2)
                n = nodeAt(ix,iy,iz)
                if (id(1,n) > 0 .and. nodeOwner(n) /= me) then
                    r = nodeOwner(n)
                    scnt(r) = scnt(r) + 1
                    sbuf(sdsp(r)+scnt(r)) = ((ix-1)*nzt + (iz-1))*nyt + iy
                    reqNode(sdsp(r)+scnt(r)) = n
                endif
            enddo
        enddo
    enddo
    call MPI_Alltoallv(sbuf, scnt, sdsp, MPI_INTEGER, rbuf, rcnt, rdsp, MPI_INTEGER, &
        MPI_COMM_WORLD, IERR)
    do i = 1, sum(rcnt)
        k = rbuf(i) - 1
        iy = mod(k, nyt) + 1
        k = k/nyt
        iz = mod(k, nzt) + 1
        ix = k/nzt + 1
        n = 0
        if (ix >= ddBox(1,1) .and. ix <= ddBox(2,1) .and. iy >= ddBox(1,2) .and. &
            iy <= ddBox(2,2) .and. iz >= ddBox(1,3) .and. iz <= ddBox(2,3)) n = nodeAt(ix,iy,iz)
        rbuf(i) = -1
        if (n > 0) rbuf(i) = pidn(n)
    enddo
    call MPI_Alltoallv(rbuf, rcnt, rdsp, MPI_INTEGER, sbuf, scnt, sdsp, MPI_INTEGER, &
        MPI_COMM_WORLD, IERR)
    do i = 1, nGhostNodes
        if (sbuf(i) < 0) then
            write(*,*) 'par.solver=3: rank', me, 'has a ghost node its owner does not solve for'
            flush(6)
            call MPI_ABORT(MPI_COMM_WORLD, 1, IERR)
        endif
        pidn(reqNode(i)) = sbuf(i)
    enddo
    deallocate(sbuf, rbuf, reqNode, scnt, sdsp, rcnt, rdsp)

    ! Fault pairs. fltMine was set with the mesh (ddMeshFinish); a pair is
    ! local when its grid point is in this block.
    fltLocal = .false.
    do ift = 1, ntotft
        iy = ddFlt(1,2,ift)
        if (iy < ddBox(1,2) .or. iy > ddBox(2,2)) cycle
        do ix = max(ddFlt(1,1,ift), ddBox(1,1)), min(ddFlt(2,1,ift), ddBox(2,1))
            do iz = max(ddFlt(1,3,ift), ddBox(1,3)), min(ddFlt(2,3,ift), ddBox(2,3))
                fltLocal(fltIndex(ix, iz, ift), ift) = .true.
            enddo
        enddo
    enddo
    nOwnF = 0
    nGhostF = 0
    do ift = 1, ntotft
        do i = 1, nftnd(ift)
            if (fltMine(i,ift)) then
                nOwnF = nOwnF + 1
                if (.not. fltLocal(i,ift)) then
                    write(*,*) 'par.solver=3: rank', me, 'owns a fault node outside its block'
                    flush(6)
                    call MPI_ABORT(MPI_COMM_WORLD, 1, IERR)
                endif
            elseif (fltLocal(i,ift)) then
                nGhostF = nGhostF + 1
            endif
        enddo
    enddo
    call MPI_Allgather(nOwnF, 1, MPI_INTEGER, cntF, 1, MPI_INTEGER, MPI_COMM_WORLD, IERR)
    offF(0) = 0
    do r = 1, nprocs-1
        offF(r) = offF(r-1) + cntF(r-1)
    enddo
    fpid = -1
    k = 0
    do ift = 1, ntotft
        do i = 1, nftnd(ift)
            if (fltMine(i,ift)) then
                fpid(i,ift) = offF(me) + k
                k = k + 1
            endif
        enddo
    enddo
    call MPI_Allreduce(MPI_IN_PLACE, fpid, nftmx*ntotft, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, IERR)
    do ift = 1, ntotft
        if (sum(cntF) /= sum(nftnd(1:ntotft)) .or. any(fpid(1:nftnd(ift),ift) < 0)) then
            if (me == 0) write(*,*) 'par.solver=3: a fault node without exactly one owner, fault', ift
            flush(6)
            call MPI_ABORT(MPI_COMM_WORLD, 1, IERR)
        endif
    enddo

    ! ---- 3. Local sets ----------------------------------------------------
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
    allocate(ghostIdx(max(1,3*nGhostNodes)))
    l = 0
    do n = 1, numnp
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

    nBE = count(et(1:numel) == 2 .or. et(1:numel) == 3)
    allocate(beEl(max(1,nBE)), beIdx(nee,max(1,nBE)), beEs(nee,nee,max(1,nBE)))

    ! Off-fault stations: each is recorded by the rank that owns its node, and
    ! a station no rank owns stops the run here.
    if (allocated(dout)) then
        allocate(stnNode(n4out), stnOwn(n4out), stnOwnG(n4out))
        stnNode = 0
        do i = 1, n4out
            j = an4nds(2,i)
            if (j > 0) then
                if (nodeOwner(j) == me) stnNode(i) = j
            endif
        enddo
        stnOwn = merge(1, 0, stnNode > 0)
        call MPI_Allreduce(stnOwn, stnOwnG, n4out, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, IERR)
        if (any(stnOwnG /= 1)) then
            if (me == 0) write(*,*) 'par.solver=3: an off-fault station without exactly one ', &
                'owning rank. owners per station =', stnOwnG
            flush(6)
            call MPI_ABORT(MPI_COMM_WORLD, 1, IERR)
        endif
        deallocate(stnOwn, stnOwnG)
        ! output_offfault_st writes a station on the rank that owns it.
        an4nds(2,1:n4out) = stnNode
    endif

    ! ---- 4. PETSc objects -------------------------------------------------
    pnOwn = nOwnEq
    pnGhost = 3*nGhostNodes
    pneq = ddGlobal(3)
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

    ! Initial state. The largest initial slip rate sets the first time step on
    ! every rank; each owner sets its pairs' velocities and sends them to the
    ! ranks whose elements touch them.
    call getScalarOnFaultQuant
    call MPI_Allreduce(maxSlipRate, dtevG, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, IERR)
    maxSlipRate = dtevG
    call initOnFaultKinematics
    call ddFltForward(consv)
    ! Fault-node area: meshgen left each rank's own quads' parts; sum them at
    ! the owner.
    call ddFltArea

    ! Mesh files, whole-mesh numbering; which lines of p1/p2output are this
    ! rank's.
    call ddWriteMeshFiles
    if (bp == 5 .or. bp == 7) call ddProfSetup
    deallocate(nodeAt)

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
    nnzO = min(81, ddGlobal(3) - nOwnEq)
    call MatSeqAIJSetPreallocation(Amat, 0, nnzD, perr)
    CHKERRA(perr)
    call MatMPIAIJSetPreallocation(Amat, 0, nnzD, 0, nnzO, perr)
    CHKERRA(perr)
    deallocate(nnzD, nnzO)
    call MatSetLocalToGlobalMapping(Amat, l2g, l2g, perr)
    CHKERRA(perr)

    ! ---- 5. Assembly of own elements; cache what each step reuses ----------
    if (me == 0) write(*,*) '= Assembling own elements on every rank        ='
    consm = 0.0d0
    ib = 0
    do nel = 1, numel
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

        ! Off-fault stations: recorded, and written after the loop, by the
        ! rank that owns the node.
        if (allocated(dout)) then
            dout(1,it) = time
            do i = 1, ndout
                j = stnNode((i-1)/(ndof*noid) + 1)
                if (j > 0) then
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
                do n = 1, numnp
                    if (id(1,n) > 0) then
                        do j = 1, ndof
                            constmp(j,n) = uloc(lbase(n)+j)
                        enddo
                    endif
                enddo
            else
                do n = 1, numnp
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

        ! bp5/bp7 profile files (output_prof, same schedule as inside faulting).
        wantProf = (bp == 5 .or. bp == 7) .and. ((status1 == 0 .and. mod(it,20) == 1) .or. &
            (status1 == 1 .and. mod(it,10) == 0))
        if (wantProf) call ddWriteProf

        status0 = status1
        dtev1 = max(dt,int8(dtev/dt)*dt)

        do n = 1, numnp
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
            ! Every rank sends its own nodes' and fault nodes' values; rank 0
            ! writes them as they arrive.
            write(proc_str,'(I5.5)') it
            output_type = 'disp'
            netcdf_outfile = trim(outDir)//'disp.'//trim(proc_str)//'.nc'
            call ddWriteNodal(netcdf_outfile, output_type, 'disp', cons)
            netcdf_outfile = trim(outDir)//'fault.'//trim(proc_str)//'.nc'
            call netcdf_write_on_fault(netcdf_outfile)
            if ((stoptag == 1) .or. (it==nstep)) then
                netcdf_outfile = trim(outDir)//'disp.r.nc'
                call ddWriteNodal(netcdf_outfile, output_type, 'disp', cons)
                netcdf_outfile = trim(outDir)//'fault.r.nc'
                call netcdf_write_on_fault(netcdf_outfile)
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

    ! The whole-fault text files, each rank its own nodes' lines. ddActive
    ! stays set: writeResults then has every rank write the stations it owns.
    call ddWriteCplot

    hwm(1) = dble(ddPeakRSSkB())/1024.0d0
    hwm(2) = hwm(1)
    call MPI_Allreduce(hwm(1), hwmG(1), 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, IERR)
    call MPI_Allreduce(hwm(2), hwmG(2), 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, IERR)

    if (me == 0) then
        write(*,*) '====================================================================='
        write(*,*) '=                        RUN SUMMARY                                ='
        write(*,'(X,A,40X,i7,4X,A)')     '= MPI ranks                = ', nprocs, '='
        write(*,'(X,A,28X,3i6,5X,A)')    '= Blocks (x, y, z)         = ', pdim, '='
        write(*,'(X,A,40X,i7,4X,A)')     '= Nodes                    = ', ddGlobal(1), '='
        write(*,'(X,A,40X,i7,4X,A)')     '= Elements                 = ', ddGlobal(2), '='
        write(*,'(X,A,40X,i7,4X,A)')     '= Equations                = ', ddGlobal(3), '='
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
        if (hwmG(1) > 0.0d0) write(*,'(X,A,2f12.1)') &
            '= Peak RSS at the end, MB: largest rank, sum       = ', hwmG(1), hwmG(2)
        write(*,*) '====================================================================='
        ! The run metadata reports the whole mesh, not rank 0's block.
        nloc = (/ numnp, numel, neq /)
        numnp = ddGlobal(1)
        numel = ddGlobal(2)
        neq = ddGlobal(3)
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
        numnp = nloc(1)
        numel = nloc(2)
        neq = nloc(3)
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

    integer (kind = 4) function faultAt(jx, jy, jz)
        ! The fault whose node rectangle holds grid point (jx, jy, jz), or 0.
        integer (kind = 4), intent(in) :: jx, jy, jz
        integer (kind = 4) :: jft
        faultAt = 0
        do jft = 1, ntotft
            if (jy == ddFlt(1,2,jft) .and. jx >= ddFlt(1,1,jft) .and. jx <= ddFlt(2,1,jft) &
                .and. jz >= ddFlt(1,3,jft) .and. jz <= ddFlt(2,3,jft)) then
                faultAt = jft
                return
            endif
        enddo
    end function faultAt

    integer (kind = 4) function fltIndex(jx, jz, jft)
        ! Index of fault jft's node at grid lines (jx, jz): strike outer, dip
        ! inner, as meshgen counts them.
        integer (kind = 4), intent(in) :: jx, jz, jft
        fltIndex = (jx-ddFlt(1,1,jft))*(ddFlt(2,3,jft)-ddFlt(1,3,jft)+1) + jz-ddFlt(1,3,jft)+1
    end function fltIndex

    integer (kind = 4) function splitsBefore(jx, jy, jz)
        ! Number of fault nodes meshgen meets before grid point (jx, jy, jz) in
        ! its ix, iz, iy order: each adds one master node to the numbering.
        integer (kind = 4), intent(in) :: jx, jy, jz
        integer (kind = 4) :: jft, nzf
        splitsBefore = 0
        do jft = 1, ntotft
            nzf = ddFlt(2,3,jft) - ddFlt(1,3,jft) + 1
            splitsBefore = splitsBefore + nzf*max(0, min(jx-1, ddFlt(2,1,jft)) - ddFlt(1,1,jft) + 1)
            if (jx >= ddFlt(1,1,jft) .and. jx <= ddFlt(2,1,jft)) then
                splitsBefore = splitsBefore + max(0, min(jz-1, ddFlt(2,3,jft)) - ddFlt(1,3,jft) + 1)
                if (jz >= ddFlt(1,3,jft) .and. jz <= ddFlt(2,3,jft) .and. ddFlt(1,2,jft) < jy) &
                    splitsBefore = splitsBefore + 1
            endif
        enddo
    end function splitsBefore

    integer (kind = 4) function ddPeakRSSkB()
        ! This process's peak resident set size (VmHWM of /proc/self/status),
        ! in kB; -1 where the file is not there.
        character (len = 128) :: line
        integer (kind = 4) :: u, ios
        ddPeakRSSkB = -1
        open(newunit = u, file = '/proc/self/status', status = 'old', action = 'read', iostat = ios)
        if (ios /= 0) return
        do
            read(u, '(A)', iostat = ios) line
            if (ios /= 0) exit
            if (line(1:6) == 'VmHWM:') then
                read(line(7:), *, iostat = ios) ddPeakRSSkB
                exit
            endif
        enddo
        close(u)
    end function ddPeakRSSkB

    subroutine ddBoundLoad
        ! bound_load on this rank: prescribed displacements of the local
        ! constrained nodes, then this rank's part of the RHS into rloc.
        integer (kind = 4) :: ii, jj, nn, code, b
        rloc = 0.0d0
        do nn = 1, numnp
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
        consf = 0.0d0
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

    subroutine ddFltArea
        ! Sum every rank's part of the fault-node area arn at the node's owner
        ! (first of the six exchange slots; the rest travel as zero).
        integer (kind = 4) :: ll
        call VecGhostGetLocalForm(fvec, floc, perr)
        CHKERRA(perr)
        call VecGetArray(floc, parr, perr)
        CHKERRA(perr)
        do ll = 1, nLocF
            parr(6*(ll-1)+1:6*ll) = 0.0d0
            parr(6*(ll-1)+1) = arn(locFlt(1,ll),locFlt(2,ll))
        enddo
        call VecRestoreArray(floc, parr, perr)
        CHKERRA(perr)
        call VecGhostRestoreLocalForm(fvec, floc, perr)
        CHKERRA(perr)
        call VecGhostUpdateBegin(fvec, ADD_VALUES, SCATTER_REVERSE, perr)
        CHKERRA(perr)
        call VecGhostUpdateEnd(fvec, ADD_VALUES, SCATTER_REVERSE, perr)
        CHKERRA(perr)
        call VecGhostGetLocalForm(fvec, floc, perr)
        CHKERRA(perr)
        call VecGetArray(floc, parr, perr)
        CHKERRA(perr)
        do ll = 1, nOwnF
            arn(locFlt(1,ll),locFlt(2,ll)) = parr(6*(ll-1)+1)
        enddo
        call VecRestoreArray(floc, parr, perr)
        CHKERRA(perr)
        call VecGhostRestoreLocalForm(fvec, floc, perr)
        CHKERRA(perr)
    end subroutine ddFltArea

    subroutine ddColsToFile(fname, vname, idname, m, n, gi, gv)
        ! Columns gi(1:n) of the m-row variable vname of the netCDF file fname,
        ! and the same entries of its index variable idname, from every rank.
        ! gi are whole-mesh numbers, increasing. Rank 0 has written the file's
        ! header (netcdf_write); it takes one rank's columns at a time and
        ! writes each run of consecutive numbers, so it never holds more than
        ! one rank's part. Integers travel as doubles, which is exact.
        use netcdf
        character (len = *), intent(in) :: fname, vname, idname
        integer (kind = 4), intent(in) :: m, n, gi(n)
        real (kind = dp), intent(in) :: gv(m,n)
        integer (kind = 4), allocatable :: ri(:)
        real (kind = dp), allocatable :: rv(:,:)
        integer (kind = 4) :: ncid, vid, iid, rr, q1, q2, cnt(1)
        cnt(1) = n
        if (me /= 0) then
            call MPI_Send(cnt, 1, MPI_INTEGER, 0, 71, MPI_COMM_WORLD, IERR)
            call MPI_Send(gi, n, MPI_INTEGER, 0, 72, MPI_COMM_WORLD, IERR)
            call MPI_Send(gv, m*n, MPI_DOUBLE_PRECISION, 0, 73, MPI_COMM_WORLD, IERR)
            return
        endif
        lastNcPath = fname
        call check(nf90_open(trim(fname), NF90_WRITE, ncid))
        call check(nf90_inq_varid(ncid, vname, vid))
        call check(nf90_inq_varid(ncid, idname, iid))
        do rr = 0, nprocs-1
            if (rr == 0) then
                allocate(ri(n), rv(m,n))
                ri = gi
                rv = gv
            else
                call MPI_Recv(cnt, 1, MPI_INTEGER, rr, 71, MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
                allocate(ri(cnt(1)), rv(m,cnt(1)))
                call MPI_Recv(ri, cnt(1), MPI_INTEGER, rr, 72, MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
                call MPI_Recv(rv, m*cnt(1), MPI_DOUBLE_PRECISION, rr, 73, MPI_COMM_WORLD, &
                    MPI_STATUS_IGNORE, IERR)
            endif
            q1 = 1
            do while (q1 <= cnt(1))
                q2 = q1
                do while (q2 < cnt(1))
                    if (ri(q2+1) /= ri(q2) + 1) exit
                    q2 = q2 + 1
                enddo
                call check(nf90_put_var(ncid, vid, rv(:,q1:q2), start = (/ 1, ri(q1) /), &
                    count = (/ m, q2-q1+1 /)))
                call check(nf90_put_var(ncid, iid, ri(q1:q2), start = (/ ri(q1) /), &
                    count = (/ q2-q1+1 /)))
                q1 = q2 + 1
            enddo
            deallocate(ri, rv)
        enddo
        call check(nf90_close(ncid))
    end subroutine ddColsToFile

    subroutine ddWriteNodal(fname, ftype, vname, a)
        ! A whole-mesh nodal file, netcdf_write's 'coor' or 'disp': the header
        ! from rank 0, then every rank's own nodes.
        character (len = 50) :: fname, ftype
        character (len = *), intent(in) :: vname
        real (kind = dp), intent(in) :: a(:,:)
        integer (kind = 4), allocatable :: gi(:)
        real (kind = dp), allocatable :: gv(:,:)
        integer (kind = 4) :: nn, q, keep
        if (me == 0) then
            keep = numnp
            numnp = ddGlobal(1)
            call netcdf_write(fname, ftype)
            numnp = keep
        endif
        allocate(gi(nOwnNodes), gv(ndof,nOwnNodes))
        q = 0
        do nn = 1, numnp
            if (nodeOwner(nn) == me) then
                q = q + 1
                gi(q) = gnode(nn)
                gv(:,q) = a(:,nn)
            endif
        enddo
        call ddColsToFile(fname, vname, 'node_id', ndof, nOwnNodes, gi, gv)
        deallocate(gi, gv)
    end subroutine ddWriteNodal

    subroutine ddWriteMeshFiles
        ! eqquasi.mesh.{coor,ien,nsmp}.nc in whole-mesh numbering, as solvers
        ! 1 and 2 write them from checkAndReport. Each rank sends its own
        ! nodes, elements and fault-1 pairs.
        integer (kind = 4), allocatable :: gi(:)
        real (kind = dp), allocatable :: gv(:,:)
        integer (kind = 4) :: jx, jy, jz, ee, kk, ii, q, keep
        character (len = 50) :: fname, ftype

        fname = trim(outDir)//'eqquasi.mesh.coor.nc'
        ftype = 'coor'
        call ddWriteNodal(fname, ftype, 'coordinates', x)

        ! Elements: whole-mesh element number, and its eight whole-mesh nodes.
        fname = trim(outDir)//'eqquasi.mesh.ien.nc'
        ftype = 'ien'
        if (me == 0) then
            keep = numel
            numel = ddGlobal(2)
            call netcdf_write(fname, ftype)
            numel = keep
        endif
        allocate(gi(numel), gv(nen,numel))
        ee = 0
        do jx = ddBox(1,1)+1, ddBox(2,1)
            do jz = ddBox(1,3)+1, ddBox(2,3)
                do jy = ddBox(1,2)+1, ddBox(2,2)
                    ee = ee + 1
                    gi(ee) = ((jx-2)*(nzt-1) + (jz-2))*(nyt-1) + jy - 1
                    do kk = 1, nen
                        gv(kk,ee) = gnode(ien(kk,ee))
                    enddo
                enddo
            enddo
        enddo
        call ddColsToFile(fname, 'node_id', 'element_id', nen, numel, gi, gv)
        deallocate(gi, gv)

        ! Split-node pairs of fault 1 (all netcdf_write stores): whole-mesh
        ! node numbers from each pair's owner. The file has nftmx columns;
        ! rank 0 writes zero past the fault's last node.
        fname = trim(outDir)//'eqquasi.mesh.nsmp.nc'
        ftype = 'nsmp'
        if (me == 0) call netcdf_write(fname, ftype)
        q = count(fltMine(1:nftnd(1),1))
        if (me == 0) q = q + nftmx - nftnd(1)
        allocate(gi(q), gv(2,q))
        gv = 0.0d0
        q = 0
        do ii = 1, nftnd(1)
            if (fltMine(ii,1)) then
                q = q + 1
                gi(q) = ii
                gv(1,q) = gnode(nsmp(1,ii,1))
                gv(2,q) = gnode(nsmp(2,ii,1))
            endif
        enddo
        if (me == 0) then
            do ii = nftnd(1)+1, nftmx
                q = q + 1
                gi(q) = ii
            enddo
        endif
        call ddColsToFile(fname, 'node_id', 'fault_node_id', 2, q, gi, gv)
        deallocate(gi, gv)
    end subroutine ddWriteMeshFiles

    subroutine ddIOCheck(fname)
        ! MPI-IO returns its errors; none is survivable here.
        character (len = *), intent(in) :: fname
        if (IERR == MPI_SUCCESS) return
        write(*,*) 'par.solver=3: rank', me, 'MPI-IO failed on ', trim(fname)
        flush(6)
        call MPI_ABORT(MPI_COMM_WORLD, 1, IERR)
    end subroutine ddIOCheck

    subroutine ddTextRecords(fname, reclen, n, pos, buf, append)
        ! This rank's n lines of a text file whose lines all have reclen bytes
        ! (newline included): line q goes to record pos(q), 0-based and
        ! increasing, counted from the start of the file, or from its present
        ! end when append is set. Called by every rank; no rank sees another
        ! rank's lines.
        character (len = *), intent(in) :: fname
        integer (kind = 4), intent(in) :: reclen, n, pos(n)
        character (len = reclen), intent(in) :: buf(n)
        logical, intent(in) :: append
        integer (kind = 4) :: fh, q1, q2
        integer (kind = MPI_OFFSET_KIND) :: base, off
        call MPI_File_open(MPI_COMM_WORLD, trim(fname), ior(MPI_MODE_WRONLY, MPI_MODE_CREATE), &
            MPI_INFO_NULL, fh, IERR)
        call ddIOCheck(fname)
        base = 0
        if (append) then
            if (me == 0) then
                call MPI_File_get_size(fh, base, IERR)
                call ddIOCheck(fname)
            endif
            call MPI_Bcast(base, 1, MPI_OFFSET, 0, MPI_COMM_WORLD, IERR)
        else
            call MPI_File_set_size(fh, base, IERR)
            call ddIOCheck(fname)
        endif
        q1 = 1
        do while (q1 <= n)
            q2 = q1
            do while (q2 < n)
                if (pos(q2+1) /= pos(q2) + 1) exit
                q2 = q2 + 1
            enddo
            off = base + int(pos(q1), MPI_OFFSET_KIND)*reclen
            call MPI_File_write_at(fh, off, buf(q1:q2), reclen*(q2-q1+1), MPI_CHARACTER, &
                MPI_STATUS_IGNORE, IERR)
            call ddIOCheck(fname)
            q1 = q2 + 1
        enddo
        call MPI_File_close(fh, IERR)
        call ddIOCheck(fname)
    end subroutine ddTextRecords

    subroutine ddProfSetup
        ! Which fault-1 nodes go to p1output.txt (k = 1) and p2output.txt
        ! (k = 2) -- output_prof's own coordinate tests, bp 5 and 7 -- and for
        ! the ones this rank owns, profI(q,k), the line profP(q,k) they take
        ! in each block of lines: the number of profile nodes before them on
        ! the fault. One integer per fault node and file, summed over ranks,
        ! once.
        integer (kind = 4), allocatable :: flag(:)
        integer (kind = 4) :: ii, kk, nf, cnt
        real (kind = dp) :: xx, zz
        logical :: s1, s2
        nf = nftnd(1)
        cnt = count(fltMine(1:nf,1))
        allocate(flag(2*max(1,nf)), profI(max(1,cnt),2), profP(max(1,cnt),2))
        flag = 0
        do ii = 1, nf
            if (.not. fltMine(ii,1)) cycle
            xx = x(1,nsmp(1,ii,1))
            zz = x(3,nsmp(1,ii,1))
            if (bp == 5) then
                s1 = abs(xx)<=38.0d3.and.abs(zz-(-10.0d3))<0.01d0
                s2 = abs(xx)<=0.01d0.and.abs(zz)<=40.0d3
            else
                s1 = abs(zz-0.0d3)<0.01d0
                s2 = abs(xx)<=0.01d0
            endif
            if (s1) flag(ii) = 1
            if (s2) flag(nf+ii) = 1
        enddo
        call MPI_Allreduce(MPI_IN_PLACE, flag, 2*nf, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, IERR)
        profN = 0
        do kk = 1, 2
            cnt = 0
            do ii = 1, nf
                if (flag((kk-1)*nf+ii) == 0) cycle
                if (fltMine(ii,1)) then
                    profN(kk) = profN(kk) + 1
                    profI(profN(kk),kk) = ii
                    profP(profN(kk),kk) = cnt
                endif
                cnt = cnt + 1
            enddo
        enddo
        deallocate(flag)
    end subroutine ddProfSetup

    subroutine ddWriteProf
        ! One block of lines appended to p1output.txt and p2output.txt, as
        ! output_prof writes it: each rank its own profile nodes' lines.
        character (len = 162), allocatable :: buf(:)
        character (len = 12) :: pname(2)
        integer (kind = 4) :: kk, q, ii
        pname = (/ 'p1output.txt', 'p2output.txt' /)
        do kk = 1, 2
            allocate(buf(profN(kk)))
            do q = 1, profN(kk)
                ii = profI(q,kk)
                write(buf(q)(1:161),'(1x,5e32.21e4)') time,fric(FR_SLIP_S,ii,1),fric(FR_SLIP_D,ii,1),fric(FR_TSTK,ii,1),fric(FR_TDIP,ii,1)
                buf(q)(162:162) = achar(10)
            enddo
            call ddTextRecords(trim(outDir)//pname(kk), 162, profN(kk), profP(1:profN(kk),kk), buf, .true.)
            deallocate(buf)
        enddo
    end subroutine ddWriteProf

    subroutine ddWriteCplot
        ! cplot_<fault>EQquasi.txt and cplot_<fault>ruptarea_trac_slip.txt, one
        ! line per fault node in the fault's node order, as
        ! output_onfault_transfer and output_ruptarea_trac_slip write them:
        ! each rank its own nodes' lines.
        character (len = 514), allocatable :: b1(:)
        character (len = 130), allocatable :: b2(:)
        integer (kind = 4), allocatable :: pos(:)
        integer (kind = 4) :: ii, jft, q
        character (len = 8) :: fttmp, faultTag
        do jft = 1, ntotft
            if (nftnd(jft) <= 0) cycle
            fttmp = faultTag(jft)
            q = count(fltMine(1:nftnd(jft),jft))
            allocate(pos(q), b1(q), b2(q))
            q = 0
            do ii = 1, nftnd(jft)
                if (.not. fltMine(ii,jft)) cycle
                q = q + 1
                pos(q) = ii - 1
                write(b1(q)(1:513),'(1x,16e32.21e4)') x(1,nsmp(1,ii,jft)), -x(3,nsmp(1,ii,jft)), &
                    fric(FR_V_TRIAL,ii,jft), fric(FR_STATE,ii,jft), fric(FR_TSTK,ii,jft), &
                    fric(FR_TDIP,ii,jft), fric(FR_TNRM,ii,jft), &
                    fric(FR_VXM,ii,jft), fric(FR_VYM,ii,jft), fric(FR_VZM,ii,jft), &
                    fric(FR_VXS,ii,jft), fric(FR_VYS,ii,jft), fric(FR_VZS,ii,jft), &
                    fric(FR_SPARE44,ii,jft), fric(FR_SPARE45,ii,jft), fnft(ii,jft)
                b1(q)(514:514) = achar(10)
                write(b2(q)(1:129),'(1x,4e32.21e4)') fric(FR_RUPT_AREA,ii,jft), fric(FR_RUPT_SLIP,ii,jft), &
                    fric(FR_TRACT_START,ii,jft), fric(FR_TRACT_END,ii,jft)
                b2(q)(130:130) = achar(10)
            enddo
            call ddTextRecords(trim(outDir)//'cplot_'//trim(fttmp)//'EQquasi.txt', 514, q, pos, b1, .false.)
            call ddTextRecords(trim(outDir)//'cplot_'//trim(fttmp)//'ruptarea_trac_slip.txt', 130, q, pos, &
                b2, .false.)
            deallocate(pos, b1, b2)
        enddo
    end subroutine ddWriteCplot

end subroutine solveTimeLoopPETScDD

subroutine ddBlocks(pdim, ncell)
! par.solver = 3: the regular block decomposition of the analytic grid. Cells
! per direction are (nxt-1, nyt-1, nzt-1). MPI_Dims_create factors nprocs into
! three non-increasing factors; the largest goes to the direction with the
! most cells. Rank r holds block (bx, by, bz), r = bx + pdim(1)*(by + pdim(2)*bz).
    use globalvar
    implicit none
    include 'mpif.h'
    integer (kind = 4), intent(out) :: pdim(3), ncell(3)
    integer (kind = 4) :: dimsS(3), ord(3), i, j, itmp

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
end subroutine ddBlocks

subroutine ddMeshBox
! par.solver = 3, after mesh4num and before allocAndInit: sets this rank's
! block (ddBox), every fault's node rectangle (ddFlt), and replaces numnp,
! numel, neq by the block's own counts (the whole mesh's go to ddGlobal), so
! allocAndInit and meshgen size and fill this rank's part only. meshgen's
! consistency check then holds the counts made here against what it generates.
!
! The block's grid lines along a direction are its cells' lower lines plus
! the one line above (the ghost layer); a grid point is owned by the block of
! the cell it is the lower corner of.
    use globalvar
    implicit none
    include 'mpif.h'

    interface
        subroutine build_yline_belt(ylinet, nyt)
            use globalvar, only: dp
            real (kind = dp), allocatable, intent(out) :: ylinet(:)
            integer (kind = 4), intent(out) :: nyt
        end subroutine build_yline_belt
        subroutine build_xzlines(fltx1, fltx2, fltz1, fltz2, xlinet, nxt, zlinet, nzt)
            use globalvar, only: dp
            real (kind = dp), intent(in) :: fltx1, fltx2, fltz1, fltz2
            real (kind = dp), allocatable, intent(out) :: xlinet(:), zlinet(:)
            integer (kind = 4), intent(out) :: nxt, nzt
        end subroutine build_xzlines
    end interface

    integer (kind = 4) :: pdim(3), ncell(3), b(3), npt(3), idir, ift, jft, k, n1, n2, n3, &
        nb, nfb, nbnd
    real (kind = dp), allocatable :: xl(:), yl(:), zl(:)
    real (kind = dp) :: tol, fltx1, fltx2, fltz1, fltz2
    logical :: bad

    call ddBlocks(pdim, ncell)
    b = (/ mod(me, pdim(1)), mod(me/pdim(1), pdim(2)), me/(pdim(1)*pdim(2)) /)
    do idir = 1, 3
        ddBox(1,idir) = (b(idir)*ncell(idir) + pdim(idir) - 1)/pdim(idir) + 1
        ddBox(2,idir) = ((b(idir)+1)*ncell(idir) + pdim(idir) - 1)/pdim(idir) + 1
    enddo

    ! The node lines, as mesh4num and meshgen build them.
    fltx1 = minval(fltxyz(1,1,:))
    fltx2 = maxval(fltxyz(2,1,:))
    fltz1 = minval(fltxyz(1,3,:))
    fltz2 = maxval(fltxyz(2,3,:))
    call build_xzlines(fltx1, fltx2, fltz1, fltz2, xl, n1, zl, n3)
    call build_yline_belt(yl, n2)
    npt = (/ nxt, nyt, nzt /)
    if (n1 /= nxt .or. n2 /= nyt .or. n3 /= nzt) then
        if (me == 0) write(*,*) 'par.solver=3: node lines disagree with mesh4num', n1, n2, n3, npt
        flush(6)
        call MPI_ABORT(MPI_COMM_WORLD, 1, IERR)
    endif

    ! Each fault's node rectangle, by meshgen's own coordinate test. The
    ! closed-form numbering needs a fault to lie on one interior y line and
    ! no two faults to share a node; anything else is refused.
    tol = dx/100.d0
    allocate(ddFlt(2,3,ntotft))
    bad = .false.
    do ift = 1, ntotft
        ddFlt(1,:,ift) = huge(1)
        ddFlt(2,:,ift) = 0
        do k = 1, nxt
            if (xl(k) >= fltxyz(1,1,ift)-tol .and. xl(k) <= fltxyz(2,1,ift)+tol) then
                ddFlt(1,1,ift) = min(ddFlt(1,1,ift), k)
                ddFlt(2,1,ift) = max(ddFlt(2,1,ift), k)
            endif
        enddo
        do k = 1, nyt
            if (yl(k) >= fltxyz(1,2,ift)-tol .and. yl(k) <= fltxyz(2,2,ift)+tol) then
                ddFlt(1,2,ift) = min(ddFlt(1,2,ift), k)
                ddFlt(2,2,ift) = max(ddFlt(2,2,ift), k)
            endif
        enddo
        do k = 1, nzt
            if (zl(k) >= fltxyz(1,3,ift)-tol .and. zl(k) <= fltxyz(2,3,ift)+tol) then
                ddFlt(1,3,ift) = min(ddFlt(1,3,ift), k)
                ddFlt(2,3,ift) = max(ddFlt(2,3,ift), k)
            endif
        enddo
        if (ddFlt(1,2,ift) /= ddFlt(2,2,ift) .or. ddFlt(1,2,ift) <= 1 .or. ddFlt(1,2,ift) >= nyt) then
            bad = .true.
        elseif (abs(yl(ddFlt(1,2,ift)) - fltxyz(1,2,ift)) >= tol) then
            bad = .true.
        elseif ((ddFlt(2,1,ift)-ddFlt(1,1,ift)+1)*(ddFlt(2,3,ift)-ddFlt(1,3,ift)+1) /= nftnd(ift)) then
            bad = .true.
        endif
        do jft = 1, ift-1
            if (ddFlt(1,2,jft) == ddFlt(1,2,ift) .and. &
                ddFlt(1,1,jft) <= ddFlt(2,1,ift) .and. ddFlt(1,1,ift) <= ddFlt(2,1,jft) .and. &
                ddFlt(1,3,jft) <= ddFlt(2,3,ift) .and. ddFlt(1,3,ift) <= ddFlt(2,3,jft)) bad = .true.
        enddo
        if (bad) then
            if (me == 0) write(*,*) 'par.solver=3: fault', ift, 'is not one rectangle of nodes on ', &
                'a single interior y line, or shares nodes with another fault. Node lines (x, y, z):', &
                ddFlt(:,:,ift), ' nodes counted by mesh4num:', nftnd(ift), '. Use par.solver = 1 or 2.'
            flush(6)
            call MPI_ABORT(MPI_COMM_WORLD, 1, IERR)
        endif
    enddo

    ! This block's counts.
    nb = product(ddBox(2,:) - ddBox(1,:) + 1)
    nbnd = 0   ! grid points on the y edges of the model: no equations
    if (ddBox(1,2) == 1) nbnd = nbnd + 1
    if (ddBox(2,2) == nyt) nbnd = nbnd + 1
    nbnd = nbnd*(ddBox(2,1)-ddBox(1,1)+1)*(ddBox(2,3)-ddBox(1,3)+1)
    nfb = 0    ! fault nodes in the block
    do ift = 1, ntotft
        k = 1
        do idir = 1, 3
            k = k*max(0, min(ddFlt(2,idir,ift), ddBox(2,idir)) - max(ddFlt(1,idir,ift), ddBox(1,idir)) + 1)
        enddo
        nfb = nfb + k
    enddo
    ddGlobal = (/ numnp, numel, neq /)
    numnp = nb + nfb
    numel = product(ddBox(2,:) - ddBox(1,:))
    neq = 3*(nb - nbnd - nfb)
    ddActive = .true.
end subroutine ddMeshBox

subroutine ddMeshFinish
! par.solver = 3, after meshgen: what the block-wise mesh pass cannot settle
! on its own. Which fault nodes this rank owns (fltMine), and one station
! numbering for all ranks -- meshgen matched stations within each block, so
! the lists are merged here and every rank gets the same slots.
    use globalvar
    implicit none
    include 'mpif.h'
    integer (kind = 4) :: own(2,3), npt(3), idir, ift, i, k, l, ix, iz, nzf, mx, ns
    integer (kind = 4), allocatable :: hit(:), hitG(:), node4(:)
    real (kind = dp), allocatable :: key(:)
    real (kind = dp) :: kk

    ! Grid lines whose points this rank owns.
    npt = (/ nxt, nyt, nzt /)
    do idir = 1, 3
        own(1,idir) = ddBox(1,idir)
        own(2,idir) = ddBox(2,idir) - 1
        if (ddBox(2,idir) == npt(idir)) own(2,idir) = npt(idir)
    enddo
    allocate(fltMine(nftmx,ntotft))
    fltMine = .false.
    do ift = 1, ntotft
        if (ddFlt(1,2,ift) < own(1,2) .or. ddFlt(1,2,ift) > own(2,2)) cycle
        nzf = ddFlt(2,3,ift) - ddFlt(1,3,ift) + 1
        do ix = max(ddFlt(1,1,ift), own(1,1)), min(ddFlt(2,1,ift), own(2,1))
            do iz = max(ddFlt(1,3,ift), own(1,3)), min(ddFlt(2,3,ift), own(2,3))
                fltMine((ix-ddFlt(1,1,ift))*nzf + iz-ddFlt(1,3,ift)+1, ift) = .true.
            enddo
        enddo
    enddo

    ! On-fault stations: anonfs(1,k) is already the node's index on the whole
    ! fault. Merge the ranks' matches, then order the slots as the whole-mesh
    ! pass meets the nodes (ix, iz, iy), which is the order solvers 1 and 2 use.
    mx = maxval(nonfs)
    ns = max(1, mx*ntotft)
    allocate(hit(ns), hitG(ns), key(ns))
    hit = 0
    do k = 1, n4onf
        hit((anonfs(3,k)-1)*mx + anonfs(2,k)) = anonfs(1,k)
    enddo
    call MPI_Allreduce(hit, hitG, ns, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, IERR)
    n4onf = 0
    do ift = 1, ntotft
        nzf = ddFlt(2,3,ift) - ddFlt(1,3,ift) + 1
        do i = 1, nonfs(ift)
            l = hitG((ift-1)*mx + i)
            if (l == 0) cycle
            ix = ddFlt(1,1,ift) + (l-1)/nzf
            iz = ddFlt(1,3,ift) + mod(l-1, nzf)
            kk = (dble(ix-1)*nzt + dble(iz-1))*nyt + dble(ddFlt(1,2,ift))
            k = n4onf
            do while (k >= 1)
                if (key(k) <= kk) exit
                key(k+1) = key(k)
                anonfs(:,k+1) = anonfs(:,k)
                k = k - 1
            enddo
            key(k+1) = kk
            anonfs(1,k+1) = l
            anonfs(2,k+1) = i
            anonfs(3,k+1) = ift
            n4onf = n4onf + 1
        enddo
    enddo
    deallocate(hit, hitG, key)

    ! Off-fault stations: one slot per matched station, in station order, on
    ! every rank. an4nds(2,k) is this rank's node for it, 0 when the node is
    ! not in this block; solveTimeLoopPETScDD checks that exactly one rank
    ! owns each.
    if (n4nds > 0) then
        allocate(node4(n4nds), hitG(n4nds))
        node4 = 0
        do k = 1, n4out
            node4(an4nds(1,k)) = an4nds(2,k)
        enddo
        call MPI_Allreduce(n4yn, hitG, n4nds, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, IERR)
        n4out = 0
        do i = 1, n4nds
            if (hitG(i) > 0) then
                n4out = n4out + 1
                an4nds(1,n4out) = i
                an4nds(2,n4out) = node4(i)
            endif
        enddo
        deallocate(node4, hitG)
    endif
end subroutine ddMeshFinish
