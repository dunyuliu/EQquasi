! Copyright (C) 2026 Dunyu Liu <dliu@ig.utexas.edu>.
! MIT
!
! Along-fault pore fluid pressure diffusion for SCEC SEAS benchmark BP8
! (quasi-dynamic fluid injection in 3D). Only active when bp == 8.
!
! On the fluid/frictional domain Omega_f, defined by xminc, xmaxc and zminc,
! this solves
!     dp/dt = alpha * lap(p) + Q_inj(t) * w(x,z) / (fluid_beta*fluid_phi*A*Lfwid)
! with zero fluid flux across the boundary of Omega_f, where the hydraulic
! diffusivity alpha = k/(phi*beta*eta). See SEAS_BP8 eqs (18)-(25).
!
! The point source is regularized either as a normalized Gaussian of
! characteristic size Lgauss (fluid_src == 1, BP8-QD-GS) or with a Peaceman
! well model (fluid_src == 2, BP8-QD-PW, operator-split explicit well update).
!
! Fault nodes are ordered n = (ix-1)*nz1 + iz (nz1: fault 1's z lines), the same convention used by
! netcdf_read_on_fault in netcdf_io.f90.
!
! Everything here runs on the master rank only, like faulting. Under
! par.solver = 3 (ddActive) every rank runs it on the fault nodes it owns:
! the loops skip the other ranks' nodes, the sums and the well cell are
! reduced over ranks, and the one line of nodes the 5-point stencil reaches
! across a block boundary comes from the rank that owns it (ddPfHalo, below).
! pf, pfActive and the rest keep the whole fault's index range on every rank,
! like fric; a rank fills its own nodes and that line.

! Subroutine #1.
! pore_pressure_init allocates the pore pressure field, builds the Omega_f
! mask, and precomputes the normalized source weights.
subroutine pore_pressure_init

    use globalvar
    use mpi
    implicit none

    integer (kind = 4) :: i, isn, jm, jp, nact
    real (kind = dp)   :: xtmp, ztmp, wsum, re
    real (kind = dp), allocatable :: halo(:)

    ! pf_neighbor walks n = (ix-1)*nz1 + iz; that needs whole columns.
    if (mod(nftnd(1), int((fltxyz(2,3,1) - fltxyz(1,3,1))/dx + 0.5d0) + 1) /= 0) then
        write(*,*) 'BP8: fault 1 node count ', nftnd(1), ' is not a whole number of', &
            ' z columns; the pore-pressure stencil would pick wrong neighbours.'
        stop 12
    endif
    allocate(pf(nftnd(1)), pfActive(nftnd(1)), pfWt(nftnd(1)))
    allocate(pfFx(nftnd(1)), pfFz(nftnd(1)))

    pf       = 0.0d0
    pfWt     = 0.0d0
    pfFx     = 1.0d0
    pfFz     = 1.0d0
    pfActive = 0
    pfWell   = 0
    fluid_pwell = 0.0d0

    ! hydraulic diffusivity, m^2/s.
    fluid_alpha = fluid_perm / (fluid_phi * fluid_beta * fluid_eta)

    if (fluid_src == 2) then
        bp8tag = 'PW'
    else
        bp8tag = 'GS'
    endif

    do i = 1, nftnd(1)
        if (ddActive) then
            if (.not. fltMine(i,1)) cycle ! par.solver = 3: another rank owns it.
        endif
        isn  = nsmp(1,i,1)
        xtmp = x(1,isn)
        ztmp = x(3,isn)
        ! Omega_f: the same region that is governed by RSF in faulting.
        if (xtmp <= xmaxc .and. xtmp >= xminc .and. abs(ztmp) <= abs(zminc)) then
            pfActive(i) = 1
            pfWt(i)     = dexp(-(xtmp*xtmp + ztmp*ztmp) / (2.0d0*fluid_Lgauss**2))
            ! well cell: the fault cell containing the injection point (0,0).
            if (abs(xtmp) < 0.5d0*dx .and. abs(ztmp) < 0.5d0*dx) pfWell = i
        endif
    enddo

    ! par.solver = 3: pfActive of the neighbours' nodes this rank's stencil
    ! reaches, and the well cell from the rank that found it.
    if (ddActive) then
        allocate(halo(nftnd(1)))
        halo = dble(pfActive)
        call ddPfHalo(halo)
        pfActive = int(halo + 0.5d0)
        deallocate(halo)
        call MPI_Allreduce(MPI_IN_PLACE, pfWell, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, IERR)
    endif

    ! Mark the boundary of Omega_f and give those nodes half a cell in the
    ! direction normal to it, a quarter at a corner. Omega_f is the closed
    ! square [-l_f, l_f]^2, so with N nodes per side the storage area must be
    ! ((N-1)*dx)^2 = (2*l_f)^2, not (N*dx)^2. Treating every node as a full cell
    ! stretches Omega_f by dx/2 on each side: 850 m instead of 800 m at
    ! dx = 50 m, which traps the injected fluid in 12.9 % too much rock and
    ! leaves the late-time uniform pressure 11.4 % low.
    do i = 1, nftnd(1)
        if (ddActive) then
            if (.not. fltMine(i,1)) cycle ! par.solver = 3: another rank owns it.
        endif
        if (pfActive(i) == 0) cycle
        ! neighbour along dip (z) missing or outside Omega_f -> z boundary
        call pf_neighbor(i, 1, jm); call pf_neighbor(i, 2, jp)
        if (jm == 0 .or. jp == 0) pfFz(i) = 0.5d0
        ! neighbour along strike (x) missing or outside Omega_f -> x boundary
        call pf_neighbor(i, 3, jm); call pf_neighbor(i, 4, jp)
        if (jm == 0 .or. jp == 0) pfFx(i) = 0.5d0
    enddo

    ! Normalize the Gaussian weights discretely so that sum(w) == 1 exactly.
    ! This removes both the analytic normalization ambiguity in SEAS_BP8 eq (20)
    ! and the discretization error of the smeared source.
    wsum = sum(pfWt)
    if (ddActive) call MPI_Allreduce(MPI_IN_PLACE, wsum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
        MPI_COMM_WORLD, IERR)
    if (wsum > 0.0d0) pfWt = pfWt / wsum

    ! Peaceman well index, m^3/(Pa s). The equivalent radius is
    ! discretization-dependent: the 2026-08-12 benchmark description gives
    ! re = 0.198*dx for a five-point finite-difference stencil and
    ! re = 0.208*sqrt(A) for finite volume on the same mesh. This solver is
    ! finite volume, so the latter applies. A is the cell area, which is dx*dx
    ! in the interior but half or a quarter of that on the boundary of Omega_f.
    re = 0.208d0 * dsqrt(dx * dx)
    fluid_WI = 2.0d0 * pi * fluid_perm * fluid_Lfwid / (fluid_eta * dlog(re/fluid_rwell))

    if (fluid_src == 2) then
        if (pfWell == 0) then
            write(*,*) 'BP8-PW: no fault node found within 0.5*dx of the ', &
                'injection point (0,0); the well cell was not located. ', &
                'Check that the fault mesh is centred on the origin and ', &
                'that dx evenly divides the domain.'
            stop 9
        endif
        if (re <= fluid_rwell) then
            write(*,*) 'BP8-PW: re =', re, &
                ' is <= fluid_rwell =', fluid_rwell, &
                '; the Peaceman well index is undefined. Increase dx or ', &
                'decrease fluid_rwell.'
            stop 10
        endif
    endif

    if (ddActive) then
        ! From here on the well cell is known only to the rank that owns it,
        ! which alone advances the well pressure.
        if (pfWell > 0) then
            if (.not. fltMine(pfWell,1)) pfWell = 0
        endif
        nact = count(pfActive(1:nftnd(1)) == 1 .and. fltMine(1:nftnd(1),1))
        call MPI_Allreduce(MPI_IN_PLACE, nact, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, IERR)
    endif

    if (me == 0) then
        write(*,*) '=     BP8 pore fluid diffusion enabled                              ='
        write(*,'(X,A,40X,E15.7,4X,A)') '=', fluid_alpha, 'm^2/s hydraulic diffusivity'
        if (ddActive) write(*,'(X,A,40X,i7,4X,A)') '=', nact, 'fault nodes in Omega_f'
        if (.not. ddActive) &
        write(*,'(X,A,40X,i7,4X,A)')    '=', sum(pfActive), 'fault nodes in Omega_f'
        if (fluid_src == 2) write(*,'(X,A,40X,E15.7,4X,A)') '=', fluid_WI, 'm^3/(Pa s) well index'
    endif

end subroutine pore_pressure_init

! Subroutine #2.
! pore_pressure_update advances the pore pressure field by dtstep with an
! explicit FTCS scheme, sub-stepping to respect dt <= dx^2/(4*alpha), then
! stores p in fric(FR_PORE_DP,:,1) and the Darcy velocity in fric(51:52,:,1).
subroutine pore_pressure_update(dtstep)

    use globalvar
    implicit none

    real (kind = dp), intent(in) :: dtstep
    real (kind = dp) :: dsub, lap, lapx, lapz, qinj, src, area, acell, dpw, pe
    real (kind = dp), allocatable :: pnew(:)
    integer (kind = 4) :: i, isub, nsub, jm, jp

    if (bp /= 8 .or. fluid_src == 0) return

    area = dx * dx

    ! injection rate, m^3/s. Constant until t_off, then off. SEAS_BP8 eq (21).
    qinj = 0.0d0
    if (time < fluid_toff) qinj = fluid_q0

    ! sub-stepping for explicit stability, with a 0.8 safety factor.
    nsub = max(1, ceiling(dtstep / (0.2d0*dx*dx/fluid_alpha)))
    dsub = dtstep / dble(nsub)

    allocate(pnew(nftnd(1)))

    do isub = 1, nsub

        pnew = pf

        do i = 1, nftnd(1)
            if (pfActive(i) == 0) cycle
            if (ddActive) then
                if (.not. fltMine(i,1)) cycle ! par.solver = 3: another rank owns it.
            endif

            ! Finite-volume 5-point Laplacian on cells of area
            ! dx^2 * pfFx(i) * pfFz(i). A face on the boundary of Omega_f
            ! carries no flux and is simply absent; the remaining faces are
            ! divided by the cell fraction in their own direction. On an edge
            ! that doubles the single normal-direction term, which is exactly
            ! the mirrored-ghost-node form of a zero-flux condition. Dropping
            ! the absent neighbour without rescaling -- what this did before --
            ! instead models a full cell hanging outside Omega_f.
            call pf_neighbor(i, 1, jm); call pf_neighbor(i, 2, jp)
            lapz = 0.0d0
            if (jm > 0) lapz = lapz + (pf(jm) - pf(i))
            if (jp > 0) lapz = lapz + (pf(jp) - pf(i))
            call pf_neighbor(i, 3, jm); call pf_neighbor(i, 4, jp)
            lapx = 0.0d0
            if (jm > 0) lapx = lapx + (pf(jm) - pf(i))
            if (jp > 0) lapx = lapx + (pf(jp) - pf(i))
            lap = (lapx/pfFx(i) + lapz/pfFz(i)) / (dx*dx)

            ! source, Pa/s, over this node's own cell area so that
            ! sum_i src_i * A_i is exactly q_inj/(beta*phi*L_fwid).
            acell = area * pfFx(i) * pfFz(i)
            src = 0.0d0
            if (fluid_src == 1) then
                ! Gaussian source (GS).
                src = qinj * pfWt(i) / (fluid_beta * fluid_phi * acell * fluid_Lfwid)
            elseif (fluid_src == 2 .and. i == pfWell) then
                ! Peaceman well (PW): flow from the well into the well cell.
                src = fluid_WI * (fluid_pwell - pf(i)) &
                    / (fluid_beta * fluid_phi * acell * fluid_Lfwid)
            endif

            pnew(i) = pf(i) + dsub * (fluid_alpha * lap + src)
        enddo

        ! operator-split explicit well pressure update. SEAS_BP8 eq (23).
        if (fluid_src == 2 .and. pfWell > 0) then
            pe  = pf(pfWell)
            dpw = (fluid_WI * (pe - fluid_pwell) + qinj) / fluid_Swell
            fluid_pwell = fluid_pwell + dsub * dpw
        endif

        pf = pnew
        ! par.solver = 3: the neighbours' new values for the next sub-step and
        ! for the Darcy velocity below.
        if (ddActive) call ddPfHalo(pf)
    enddo

    deallocate(pnew)

    ! Hand the pore pressure to faulting through the existing fric(FR_PORE_DP,:,:) slot
    ! and record the Darcy velocity q = -(k/eta) grad(p) for benchmark output.
    do i = 1, nftnd(1)
        if (ddActive) then
            if (.not. fltMine(i,1)) cycle ! par.solver = 3: another rank owns it.
        endif
        fric(FR_PORE_DP,i,1) = pf(i)

        fric(FR_DARCY_S,i,1) = 0.0d0
        fric(FR_DARCY_D,i,1) = 0.0d0
        if (pfActive(i) == 1) then
            call pf_neighbor(i, 3, jm); call pf_neighbor(i, 4, jp)
            if (jm > 0 .and. jp > 0) &
                fric(FR_DARCY_S,i,1) = -fluid_perm/fluid_eta * (pf(jp)-pf(jm))/(2.0d0*dx)
            call pf_neighbor(i, 1, jm); call pf_neighbor(i, 2, jp)
            if (jm > 0 .and. jp > 0) &
                fric(FR_DARCY_D,i,1) = -fluid_perm/fluid_eta * (pf(jp)-pf(jm))/(2.0d0*dx)
        endif
    enddo

end subroutine pore_pressure_update

! Subroutine #3.
! bp8_profile_init locates the fault nodes lying on the two cross-section lines
! through the injection point -- along strike (x3 = 0) and along dip (x2 = 0) --
! and allocates storage for the section 4.3 output of SEAS BP8.
!
! The benchmark asks for nodes at exactly 10 m spacing over [-400, 400], i.e. 81
! columns. That is what the production compset at dx = 10 m provides; a coarser
! dx yields correspondingly fewer columns, and the header records the actual
! spacing so the file stays self-describing rather than silently wrong.
subroutine bp8_profile_init

    use globalvar
    use mpi
    implicit none

    integer (kind = 4) :: i, isn, is, idp, loc(2), tot(2), off(2)
    real (kind = dp)   :: xtmp, ztmp, tol

    tol = 0.01d0 * dx

    nProfStrike = 0
    nProfDepth  = 0
    do i = 1, nftnd(1)
        if (ddActive) then
            if (.not. fltMine(i,1)) cycle ! par.solver = 3: another rank owns it.
        endif
        if (pfActive(i) == 0) cycle
        isn  = nsmp(1,i,1)
        xtmp = x(1,isn)
        ztmp = x(3,isn)
        if (abs(ztmp) < tol) nProfStrike = nProfStrike + 1
        if (abs(xtmp) < tol) nProfDepth  = nProfDepth  + 1
    enddo

    ! par.solver = 3: the counts above are this rank's own nodes. A line runs
    ! along x or along z through blocks whose rank number grows with it
    ! (rank = bx + pdim(1)*(by + pdim(2)*bz)), so the nodes of the lower ranks
    ! are the ones before this rank's on the line.
    if (ddActive) then
        loc = (/ nProfStrike, nProfDepth /)
        off = 0
        call MPI_Exscan(loc, off, 2, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, IERR)
        if (me == 0) off = 0
        call MPI_Allreduce(loc, tot, 2, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, IERR)
        nProfStrike = tot(1)
        nProfDepth  = tot(2)
    endif

    if (nProfStrike == 0 .or. nProfDepth == 0) return

    ! Subsample in time so the files stay near 1000 rows, as the description asks.
    nProfSkip = max(1, nstep/1000)
    nProfMax  = nstep/nProfSkip + 2

    allocate(idProfStrike(nProfStrike), xProfStrike(nProfStrike))
    allocate(idProfDepth(nProfDepth),   xProfDepth(nProfDepth))
    allocate(profTime(nProfMax), profVmax(nProfMax))
    allocate(profS2s(nProfMax,nProfStrike), profS3s(nProfMax,nProfStrike), &
             profT2s(nProfMax,nProfStrike), profT3s(nProfMax,nProfStrike), &
             profPs (nProfMax,nProfStrike))
    allocate(profS2d(nProfMax,nProfDepth),  profS3d(nProfMax,nProfDepth),  &
             profT2d(nProfMax,nProfDepth),  profT3d(nProfMax,nProfDepth),  &
             profPd (nProfMax,nProfDepth))

    profTime = 0.0d0; profVmax = 0.0d0
    profS2s = 0.0d0; profS3s = 0.0d0; profT2s = 0.0d0; profT3s = 0.0d0; profPs = 0.0d0
    profS2d = 0.0d0; profS3d = 0.0d0; profT2d = 0.0d0; profT3d = 0.0d0; profPd = 0.0d0

    ! Fault nodes are ordered n = (ix-1)*nz1 + iz, so walking i in order gives
    ! increasing x3 along the dip line and increasing x2 along the strike line.
    is = 0; idp = 0
    ! par.solver = 3: a rank fills its own nodes' slots of the two lines; the
    ! other slots stay 0 (bp8_profile_record skips them).
    if (ddActive) then
        idProfStrike = 0; idProfDepth = 0
        xProfStrike = 0.0d0; xProfDepth = 0.0d0
        is = off(1); idp = off(2)
    endif
    do i = 1, nftnd(1)
        if (ddActive) then
            if (.not. fltMine(i,1)) cycle ! par.solver = 3: another rank owns it.
        endif
        if (pfActive(i) == 0) cycle
        isn  = nsmp(1,i,1)
        xtmp = x(1,isn)
        ztmp = x(3,isn)
        if (abs(ztmp) < tol) then
            is = is + 1
            idProfStrike(is) = i
            xProfStrike(is)  = xtmp
        endif
        if (abs(xtmp) < tol) then
            idp = idp + 1
            idProfDepth(idp) = i
            ! BP8 measures x3 positive downward; EQquasi's z is positive up.
            xProfDepth(idp)  = -ztmp
        endif
    enddo

    ! x3 = -z, so walking nodes in increasing z built the depth line in
    ! DECREASING x3. The benchmark expects the coordinate row to ascend, like
    ! the strike line. Reverse it.
    call reverse_profile_line(idProfDepth, xProfDepth, nProfDepth)

    ! par.solver = 3: the coordinate row of each line, from the owners.
    if (ddActive) then
        call MPI_Allreduce(MPI_IN_PLACE, xProfStrike, nProfStrike, MPI_DOUBLE_PRECISION, MPI_SUM, &
            MPI_COMM_WORLD, IERR)
        call MPI_Allreduce(MPI_IN_PLACE, xProfDepth, nProfDepth, MPI_DOUBLE_PRECISION, MPI_SUM, &
            MPI_COMM_WORLD, IERR)
    endif

    if (me == 0) then
        write(*,'(X,A,40X,i7,4X,A)') '= BP8 profile nodes along strike = ', nProfStrike, '='
        write(*,'(X,A,40X,i7,4X,A)') '= BP8 profile nodes along dip    = ', nProfDepth, '='
    endif

end subroutine bp8_profile_init

! Subroutine #4.
! bp8_profile_record stores one time slice along both cross-section lines.
subroutine bp8_profile_record

    use globalvar
    implicit none

    integer (kind = 4) :: k, i

    if (nProfStrike == 0 .or. nProfDepth == 0) return
    if (nProfRec >= nProfMax) return
    if (mod(it, nProfSkip) /= 0 .and. it /= 1) return

    nProfRec = nProfRec + 1
    k = nProfRec
    profTime(k) = time
    profVmax(k) = maxSlipRate

    do i = 1, nProfStrike
        if (ddActive) then
            if (idProfStrike(i) == 0) cycle ! par.solver = 3: another rank owns it.
        endif
        profS2s(k,i) =  fric(FR_SLIP_S, idProfStrike(i), 1)          ! slip along strike, m
        profS3s(k,i) = -fric(FR_SLIP_D, idProfStrike(i), 1)          ! slip, positive down
        profT2s(k,i) =  fric(FR_TSTK, idProfStrike(i), 1)/1.0d6    ! shear stress, MPa
        profT3s(k,i) = -fric(FR_TDIP, idProfStrike(i), 1)/1.0d6
        profPs (k,i) =  fric(FR_PORE_DP,  idProfStrike(i), 1)/1.0d6    ! pore pressure, MPa
    enddo
    do i = 1, nProfDepth
        if (ddActive) then
            if (idProfDepth(i) == 0) cycle ! par.solver = 3: another rank owns it.
        endif
        profS2d(k,i) =  fric(FR_SLIP_S, idProfDepth(i), 1)
        profS3d(k,i) = -fric(FR_SLIP_D, idProfDepth(i), 1)
        profT2d(k,i) =  fric(FR_TSTK, idProfDepth(i), 1)/1.0d6
        profT3d(k,i) = -fric(FR_TDIP, idProfDepth(i), 1)/1.0d6
        profPd (k,i) =  fric(FR_PORE_DP,  idProfDepth(i), 1)/1.0d6
    enddo

end subroutine bp8_profile_record

! reverse_profile_line flips a profile line so its coordinates ascend.
subroutine reverse_profile_line(ids, coords, n)

    use globalvar, only : dp
    implicit none

    integer (kind = 4), intent(in)    :: n
    integer (kind = 4), intent(inout) :: ids(n)
    real (kind = dp),   intent(inout) :: coords(n)

    integer (kind = 4) :: i, itmp
    real (kind = dp)   :: rtmp

    do i = 1, n/2
        itmp        = ids(i)
        ids(i)      = ids(n+1-i)
        ids(n+1-i)  = itmp
        rtmp          = coords(i)
        coords(i)     = coords(n+1-i)
        coords(n+1-i) = rtmp
    enddo

end subroutine reverse_profile_line

! pf_neighbor returns the fault-node index of the Omega_f neighbour of node i
! in the requested stencil direction, or 0 if there is none -- either the grid
! edge or a neighbour that lies outside Omega_f (pfActive == 0). This is the
! single place that turns the n = (ix-1)*nz1 + iz node ordering into a 4-point
! stencil; pore_pressure_init (boundary cell fractions) and pore_pressure_update
! (the Laplacian and the Darcy velocity) all call it instead of repeating the
! bounds-and-active checks inline.
! dir: 1 = -dip (iz-1), 2 = +dip (iz+1), 3 = -strike (i-nz1), 4 = +strike (i+nz1)
subroutine pf_neighbor(i, dir, j)

    use globalvar
    implicit none

    integer (kind = 4), intent(in)  :: i, dir
    integer (kind = 4), intent(out) :: j

    integer (kind = 4) :: iz, nz1

    ! Fault 1's own count of z lines. Equal to the mesh's nzt while the fault
    ! spans the whole domain in z; a buried fault (1.21, enlarging_ratio_xz)
    ! has fewer, and walking it with nzt would pick wrong neighbours silently.
    nz1 = int((fltxyz(2,3,1) - fltxyz(1,3,1))/dx + 0.5d0) + 1
    iz = mod(i-1, nz1) + 1
    j = 0
    select case (dir)
    case (1)
        if (iz > 1) j = i - 1
    case (2)
        if (iz < nz1) j = i + 1
    case (3)
        if (i > nz1) j = i - nz1
    case (4)
        if (i + nz1 <= nftnd(1)) j = i + nz1
    end select
    if (j > 0) then
        if (pfActive(j) == 0) j = 0
    endif

end subroutine pf_neighbor

! ddPfHalo, par.solver = 3 only. a is a fault-1 nodal field of which this rank
! holds its own nodes' values; on return it also holds the values of the other
! ranks' nodes that pf_neighbor reaches from this rank's nodes -- one line of
! nodes across each block boundary on the fault. Rank to rank (MPI_Isend /
! MPI_Irecv between blocks that share a boundary); nothing passes through
! rank 0. The lists are built on the first call: rank A sends rank B, in
! increasing node number, its nodes that have a stencil neighbour owned by B,
! which is the list B builds of what it needs from A.
subroutine ddPfHalo(a)

    use globalvar
    use mpi
    implicit none

    real (kind = dp), intent(inout) :: a(nftnd(1))

    integer (kind = 4), save :: nnb = -1
    integer (kind = 4), allocatable, save :: nbr(:), sOff(:), rOff(:), sIdx(:), rIdx(:)
    integer (kind = 4), allocatable :: cnt(:,:), last(:,:), keep(:), req(:)
    real (kind = dp), allocatable :: sbuf(:), rbuf(:)
    integer (kind = 4) :: pdim(3), ncell(3), nzf, j, k, dir, r, q, pass

    if (nnb < 0) then
        call ddBlocks(pdim, ncell)
        nzf = ddFlt(2,3,1) - ddFlt(1,3,1) + 1
        ! The owner of a node in closed form must be the owner fltMine says,
        ! and pf_neighbor's z columns the fault rectangle's.
        k = 0
        if (nzf /= int((fltxyz(2,3,1) - fltxyz(1,3,1))/dx + 0.5d0) + 1) k = 1
        do j = 1, nftnd(1)
            if ((pfOwner(j) == me) .neqv. fltMine(j,1)) k = 1
        enddo
        if (k /= 0) then
            write(*,*) 'par.solver=3: rank', me, 'cannot tell which rank owns each fault node ', &
                'for the pore-pressure exchange (bp == 8).'
            flush(6)
            call MPI_ABORT(MPI_COMM_WORLD, 1, IERR)
        endif
        ! pf_neighbor drops a neighbour outside Omega_f; the lists need every
        ! grid neighbour (pfActive itself is exchanged through them).
        allocate(cnt(2,0:nprocs-1), last(2,0:nprocs-1), keep(nftnd(1)))
        keep = pfActive
        pfActive = 1
        cnt = 0
        do pass = 1, 2
            if (pass == 2) then
                nnb = count(cnt(1,:) + cnt(2,:) > 0)
                allocate(nbr(nnb), sOff(nnb+1), rOff(nnb+1), sIdx(sum(cnt(1,:))), rIdx(sum(cnt(2,:))))
                q = 0
                sOff(1) = 0
                rOff(1) = 0
                do r = 0, nprocs-1
                    if (cnt(1,r) + cnt(2,r) == 0) cycle
                    q = q + 1
                    nbr(q) = r
                    sOff(q+1) = sOff(q) + cnt(1,r)
                    rOff(q+1) = rOff(q) + cnt(2,r)
                    cnt(1,r) = sOff(q)
                    cnt(2,r) = rOff(q)
                enddo
            endif
            last = 0
            do j = 1, nftnd(1)
                do dir = 1, 4
                    call pf_neighbor(j, dir, k)
                    if (k == 0) cycle
                    if (fltMine(j,1) .eqv. fltMine(k,1)) cycle
                    if (fltMine(j,1)) then
                        q = 1          ! mine, next to a node of rank r: send
                        r = pfOwner(k)
                    else
                        q = 2          ! rank r's, next to a node of mine: receive
                        r = pfOwner(j)
                    endif
                    if (last(q,r) == j) cycle
                    last(q,r) = j
                    cnt(q,r) = cnt(q,r) + 1
                    if (pass == 2 .and. q == 1) sIdx(cnt(1,r)) = j
                    if (pass == 2 .and. q == 2) rIdx(cnt(2,r)) = j
                enddo
            enddo
        enddo
        pfActive = keep
        deallocate(cnt, last, keep)
    endif

    allocate(sbuf(size(sIdx)), rbuf(size(rIdx)), req(2*nnb))
    req = MPI_REQUEST_NULL
    do q = 1, nnb
        if (rOff(q+1) > rOff(q)) call MPI_Irecv(rbuf(rOff(q)+1), rOff(q+1)-rOff(q), &
            MPI_DOUBLE_PRECISION, nbr(q), 81, MPI_COMM_WORLD, req(q), IERR)
    enddo
    sbuf = a(sIdx)
    do q = 1, nnb
        if (sOff(q+1) > sOff(q)) call MPI_Isend(sbuf(sOff(q)+1), sOff(q+1)-sOff(q), &
            MPI_DOUBLE_PRECISION, nbr(q), 81, MPI_COMM_WORLD, req(nnb+q), IERR)
    enddo
    call MPI_Waitall(2*nnb, req, MPI_STATUSES_IGNORE, IERR)
    a(rIdx) = rbuf
    deallocate(sbuf, rbuf, req)

contains

    integer (kind = 4) function pfOwner(n)
        ! Rank that owns fault-1 node n: the block of the cell its grid point
        ! is the lower corner of (solveTimeLoopPETScDD's blockRank).
        integer (kind = 4), intent(in) :: n
        integer (kind = 4) :: cx, cy, cz
        cx = min(ddFlt(1,1,1) + (n-1)/nzf, ncell(1))
        cy = min(ddFlt(1,2,1), ncell(2))
        cz = min(ddFlt(1,3,1) + mod(n-1, nzf), ncell(3))
        pfOwner = ((cx-1)*pdim(1))/ncell(1) + pdim(1)*(((cy-1)*pdim(2))/ncell(2) &
            + pdim(2)*(((cz-1)*pdim(3))/ncell(3)))
    end function pfOwner

end subroutine ddPfHalo

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
