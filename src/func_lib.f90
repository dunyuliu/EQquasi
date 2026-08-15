! Table of Contents of functions and subroutines.
! #1 insert_rough_fault
! #2 rsf_rd
! #3 build_yline_belt

! #1 insert_rough_fault
subroutine insert_rough_fault(xcoor, ycoor, zcoor, ycoort, pfx, pfz, ymax1, ymin1)
    ! This subroutine is to modify ycoor if a rough_fault interface is inserted.
    use globalvar
    implicit none
    real (kind = dp) :: xcoor, ycoor, zcoor, peak, ycoort, pfx, pfz, ymax1, ymin1
    real (kind = dp) :: fx1, fx2, fz1, tol
    integer (kind = 4) :: ixx, izz
    tol = dx/1.0d3
    
    fx1 = rough_fx_min
    fx2 = rough_fx_max
    fz1 = rough_fz_min
    if ((xcoor < fx2 + tol) .and. (xcoor > fx1 - tol) .and. (zcoor > fz1 - tol)) then 
        ixx = (xcoor - fx1)/dx + 1
        izz = (zcoor - fz1)/dx + 1
    elseif ((xcoor < fx1 - tol) .and. (zcoor > fz1 - tol) ) then
        ixx = 1
        izz = (zcoor - fz1)/dx + 1
    elseif ((xcoor > fx2 + tol) .and. (zcoor > fz1 - tol)) then 
        ixx = nnx
        izz = (zcoor - fz1)/dx + 1
    elseif ((xcoor < fx2 + tol) .and. (xcoor > fx1 - tol) .and. (zcoor < fz1 - tol)) then 
        ixx = (xcoor - fx1)/dx + 1
        izz = 1
    elseif ((xcoor < fx1 - tol) .and. (zcoor < fz1 - tol)) then 
        ixx = 1
        izz = 1 
    elseif ((xcoor > fx2 + tol) .and. (zcoor < fz1 - tol)) then 
        ixx = nnx
        izz = 1
    endif 
    
    peak = rough_geo(1,nnz*(ixx-1)+izz)
    pfx = rough_geo(2,nnz*(ixx-1)+izz)
    pfz = rough_geo(3,nnz*(ixx-1)+izz)    
    
    if (ycoor > -tol) then
        ycoort = ycoor*(ymax1 - peak)/ymax1 + peak
    elseif (ycoor < -tol) then 
        ycoort = ycoor*(peak - ymin1)/(-ymin1) + peak 
    endif 
    
end subroutine insert_rough_fault

! #2 rsf
subroutine rsf_rd(t_shear, t_norm, t_a, t_b, t_f0, t_v0, t_vs, t_rou, t_vload)
! rsf_rd stands for rate- and state- friction with radiative damping.
! It calculates the shear stress given necessary parameters and returns to t_shear.
    use globalvar
    implicit none
    real (kind = dp) :: t_shear, t_norm, t_a, t_b, t_f0, t_v0, t_vs, t_rou, t_vload, t_fric
    t_fric = t_a * dasinh(t_vload/2.0d0/t_v0 * dexp((t_f0 + t_b*dlog(t_v0/t_vload))/t_a)) 
    t_shear = - t_norm * t_fric + t_rou*t_vs/2.0d0*t_vload
end subroutine rsf_rd

! #3 build_yline_belt
! Builds the y node-line array (ylinet, length nyt) shared by mesh4num and
! meshgen: a stretched far-field region on each side of a single uniform-dy
! belt spanning [min fault y, max fault y] (unchanged from v1.5.0).
!
! v1.5.0 widened that belt to span every fault's y, not just fault 1's, but
! still stepped by a single dy from the minimum -- so a fault whose offset
! from the minimum was not an integer multiple of dy fell between node
! lines and meshed with zero fault nodes. Silently: the run started,
! completed thousands of steps, wrote both faults' slabs to fault.*.nc (one
! all zeros), and the only hint was "Fault nodes = 0" for fault 1 alone in
! the run summary -- a line that read as cosmetic, not fatal.
!
! Fix, by owner's direction: REFUSE, do not reshape the mesh. The belt
! geometry stays exactly what it always was; every declared fault's y must
! already be an integer number of dy steps from the belt origin (the
! minimum fault y), checked here before any of that geometry is built. A
! non-integer offset stops the run at setup with a boxed, named error
! (fault index, its y, the belt origin, dy, the resulting ratio) rather
! than silently producing an unmeshed fault -- pick dy so the offset
! divides evenly (compsets/bp1002.qdc.2500's 5000 m offset needs dy in
! {1000, 1250, 2500, ...}, not 2000).
!
! nftnd(ift) == 0 is guarded again, independently, in mesh4num.f90 once
! actual node counts are known: belt-and-braces on purpose, so a future
! change to either check cannot alone let a zero-node fault through.
subroutine build_yline_belt(ylinet, nyt)

    ! ONLY-import: the dummy nyt below is this routine's own output and must
    ! not collide with globalvar's module-level nyt (the caller assigns the
    ! result into that module variable itself, by argument association).
    use globalvar, only: dp, dx, dis4uniF, dis4uniB, rat, dymax, np, ymin, &
        ymax, fltxyz, ntotft
    implicit none

    real (kind = dp), allocatable, intent(out) :: ylinet(:)
    integer (kind = 4), intent(out) :: nyt

    real (kind = dp) :: dy, tol, ystep, ycoor, yflt_lo, yflt_hi
    real (kind = dp) :: offset, ratio
    integer (kind = 4) :: iy, edgey1, nyuni, ift, nearest_n

    dy  = dx
    tol = dx/100.d0

    yflt_lo = minval(fltxyz(1,2,1:ntotft))
    yflt_hi = maxval(fltxyz(2,2,1:ntotft))

    ! ---- refuse up front: every fault's y must be an integer number of dy
    ! steps from the belt origin (yflt_lo), or it will not land on a node
    ! line below. See faultgeom in the compset for how to change ift's y;
    ! see this message for how to change dy instead. ----
    do ift = 1, ntotft
        offset = fltxyz(1,2,ift) - yflt_lo
        ratio  = offset / dy
        nearest_n = int(ratio + 0.5d0)
        if (abs(offset - dble(nearest_n)*dy) > tol) then
            write(*,*) '====================================================================='
            write(*,*) '=                          MESH ERROR                               ='
            write(*,*) '= A declared fault y-offset is not an integer multiple of dy.        ='
            write(*,*) '= EQquasi refuses to mesh this rather than silently placing the      ='
            write(*,*) '= fault off a node line, which meshes it with ZERO fault nodes.      ='
            write(*,*) '=                                                                    ='
            write(*,'(X,A,I4,4X,A)')          '=   fault index          = ', ift, '='
            write(*,'(X,A,E16.8,A,4X,A)')     '=   fault y              = ', fltxyz(1,2,ift), ' m', '='
            write(*,'(X,A,E16.8,A,4X,A)')     '=   belt origin (min fault y) = ', yflt_lo, ' m', '='
            write(*,'(X,A,E16.8,A,4X,A)')     '=   dy                   = ', dy, ' m', '='
            write(*,'(X,A,E16.8,A,4X,A)')     '=   offset (fault y - belt origin) = ', offset, ' m', '='
            write(*,'(X,A,E16.8,4X,A)')       '=   offset / dy          = ', ratio, '='
            write(*,*) '=                                                                    ='
            write(*,*) '= Fix: choose dy so that offset/dy is an integer. Example valid dy   ='
            write(*,*) '= values for THIS offset:                                            ='
            write(*,'(X,A,F10.2,A,F10.2,A,F10.2,A,4X,A)') '=   dy = ', offset/2.0d0, ' m,  dy = ', &
                offset/4.0d0, ' m,  dy = ', offset/5.0d0, ' m', '='
            write(*,*) '====================================================================='
            stop 4
        endif
    enddo

    nyuni=int((yflt_hi-yflt_lo)/dy+0.5d0)+dis4uniF+dis4uniB+1
    ystep=dy
    ycoor=yflt_lo-dy*(dis4uniF)
    do iy=1,np
        ystep=ystep*rat
        if (ystep>=dymax) ystep = dymax
        ycoor=ycoor-ystep
        if(ycoor<=ymin) exit
    enddo
    edgey1=iy
    ystep=dy
    ycoor=yflt_hi+dy*(dis4uniB)
    do iy=1,np
        ystep=ystep*rat
        if (ystep>=dymax) ystep = dymax
        ycoor=ycoor+ystep
        if(ycoor>=ymax) exit
    enddo
    nyt=nyuni+edgey1+iy
    allocate(ylinet(nyt))
    ylinet(edgey1+1)=yflt_lo-dy*(dis4uniF)
    ystep=dy
    do iy=edgey1,1,-1
        ystep=ystep*rat
        if (ystep>=dymax) ystep = dymax
        ylinet(iy)=ylinet(iy+1)-ystep
    enddo
    do iy=edgey1+2,edgey1+nyuni
        ylinet(iy)=ylinet(iy-1)+dy
    enddo
    ystep=dy
    do iy=edgey1+nyuni+1,nyt
        ystep=ystep*rat
        if (ystep>=dymax) ystep = dymax
        ylinet(iy)=ylinet(iy-1)+ystep
    enddo

end subroutine build_yline_belt
