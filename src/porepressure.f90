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
! Fault nodes are ordered n = (ix-1)*nzt + iz, the same convention used by
! netcdf_read_on_fault in netcdf_io.f90.
!
! Everything here runs on the master rank only, like faulting.

! Subroutine #1.
! pore_pressure_init allocates the pore pressure field, builds the Omega_f
! mask, and precomputes the normalized source weights.
subroutine pore_pressure_init

    use globalvar
    implicit none

    integer (kind = 4) :: i, isn, jm, jp
    real (kind = dp)   :: xtmp, ztmp, wsum, re

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

    ! Mark the boundary of Omega_f and give those nodes half a cell in the
    ! direction normal to it, a quarter at a corner. Omega_f is the closed
    ! square [-l_f, l_f]^2, so with N nodes per side the storage area must be
    ! ((N-1)*dx)^2 = (2*l_f)^2, not (N*dx)^2. Treating every node as a full cell
    ! stretches Omega_f by dx/2 on each side: 850 m instead of 800 m at
    ! dx = 50 m, which traps the injected fluid in 12.9 % too much rock and
    ! leaves the late-time uniform pressure 11.4 % low.
    do i = 1, nftnd(1)
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
    if (wsum > 0.0d0) pfWt = pfWt / wsum

    ! Peaceman well index, m^3/(Pa s). re = 0.198*dx for a 5-point stencil.
    re = 0.198d0 * dx
    fluid_WI = 2.0d0 * pi * fluid_perm * fluid_Lfwid / (fluid_eta * dlog(re/fluid_rwell))

    if (me == 0) then
        write(*,*) '=     BP8 pore fluid diffusion enabled                              ='
        write(*,'(X,A,40X,E15.7,4X,A)') '=', fluid_alpha, 'm^2/s hydraulic diffusivity'
        write(*,'(X,A,40X,i7,4X,A)')    '=', sum(pfActive), 'fault nodes in Omega_f'
        if (fluid_src == 2) write(*,'(X,A,40X,E15.7,4X,A)') '=', fluid_WI, 'm^3/(Pa s) well index'
    endif

end subroutine pore_pressure_init

! Subroutine #2.
! pore_pressure_update advances the pore pressure field by dtstep with an
! explicit FTCS scheme, sub-stepping to respect dt <= dx^2/(4*alpha), then
! stores p in fric(6,:,1) and the Darcy velocity in fric(51:52,:,1).
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
    enddo

    deallocate(pnew)

    ! Hand the pore pressure to faulting through the existing fric(6,:,:) slot
    ! and record the Darcy velocity q = -(k/eta) grad(p) for benchmark output.
    do i = 1, nftnd(1)
        fric(6,i,1) = pf(i)

        fric(51,i,1) = 0.0d0
        fric(52,i,1) = 0.0d0
        if (pfActive(i) == 1) then
            call pf_neighbor(i, 3, jm); call pf_neighbor(i, 4, jp)
            if (jm > 0 .and. jp > 0) &
                fric(51,i,1) = -fluid_perm/fluid_eta * (pf(jp)-pf(jm))/(2.0d0*dx)
            call pf_neighbor(i, 1, jm); call pf_neighbor(i, 2, jp)
            if (jm > 0 .and. jp > 0) &
                fric(52,i,1) = -fluid_perm/fluid_eta * (pf(jp)-pf(jm))/(2.0d0*dx)
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
    implicit none

    integer (kind = 4) :: i, isn, is, idp
    real (kind = dp)   :: xtmp, ztmp, tol

    tol = 0.01d0 * dx

    nProfStrike = 0
    nProfDepth  = 0
    do i = 1, nftnd(1)
        if (pfActive(i) == 0) cycle
        isn  = nsmp(1,i,1)
        xtmp = x(1,isn)
        ztmp = x(3,isn)
        if (abs(ztmp) < tol) nProfStrike = nProfStrike + 1
        if (abs(xtmp) < tol) nProfDepth  = nProfDepth  + 1
    enddo

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

    ! Fault nodes are ordered n = (ix-1)*nzt + iz, so walking i in order gives
    ! increasing x3 along the dip line and increasing x2 along the strike line.
    is = 0; idp = 0
    do i = 1, nftnd(1)
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
        profS2s(k,i) =  fric(71, idProfStrike(i), 1)          ! slip along strike, m
        profS3s(k,i) = -fric(72, idProfStrike(i), 1)          ! slip, positive down
        profT2s(k,i) =  fric(28, idProfStrike(i), 1)/1.0d6    ! shear stress, MPa
        profT3s(k,i) = -fric(29, idProfStrike(i), 1)/1.0d6
        profPs (k,i) =  fric(6,  idProfStrike(i), 1)/1.0d6    ! pore pressure, MPa
    enddo
    do i = 1, nProfDepth
        profS2d(k,i) =  fric(71, idProfDepth(i), 1)
        profS3d(k,i) = -fric(72, idProfDepth(i), 1)
        profT2d(k,i) =  fric(28, idProfDepth(i), 1)/1.0d6
        profT3d(k,i) = -fric(29, idProfDepth(i), 1)/1.0d6
        profPd (k,i) =  fric(6,  idProfDepth(i), 1)/1.0d6
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
! single place that turns the n = (ix-1)*nzt + iz node ordering into a 4-point
! stencil; pore_pressure_init (boundary cell fractions) and pore_pressure_update
! (the Laplacian and the Darcy velocity) all call it instead of repeating the
! bounds-and-active checks inline.
! dir: 1 = -dip (iz-1), 2 = +dip (iz+1), 3 = -strike (i-nzt), 4 = +strike (i+nzt)
subroutine pf_neighbor(i, dir, j)

    use globalvar
    implicit none

    integer (kind = 4), intent(in)  :: i, dir
    integer (kind = 4), intent(out) :: j

    integer (kind = 4) :: iz

    iz = mod(i-1, nzt) + 1
    j = 0
    select case (dir)
    case (1)
        if (iz > 1) j = i - 1
    case (2)
        if (iz < nzt) j = i + 1
    case (3)
        if (i > nzt) j = i - nzt
    case (4)
        if (i + nzt <= nftnd(1)) j = i + nzt
    end select
    if (j > 0) then
        if (pfActive(j) == 0) j = 0
    endif

end subroutine pf_neighbor
