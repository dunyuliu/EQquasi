! Table of Contents of functions and subroutines. 
! #1 insert_rough_fault

! #1 insert_rough_fault
subroutine insert_rough_fault(xcoor, ycoor, zcoor, ycoort, pfx, pfz)
	! This subroutine is to modify ycoor if a rough_fault interface is inserted.
	
	use globalvar
	implicit none
	real (kind = dp) :: xcoor, ycoor, zcoor, peak, ycoort, pfx, pfz
	real (kind = dp) :: fx1, fx2, fz1, tol
	integer (kind = 4) :: ixx, izz
	tol = dx/1.0d3
	
	fx1 = fltxyz(1,1,1)
	fx2 = fltxyz(2,1,1)
	fz1 = fltxyz(1,3,1)
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
		ycoort = ycoor*(ymax - peak)/ymax + peak
	elseif (ycoor < -tol) then 
		ycoort = ycoor*(peak - ymin)/(-ymin) + peak 
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