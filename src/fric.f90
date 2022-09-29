SUBROUTINE slip_weak(slip,fricsgl,xmu)
  use globalvar
  implicit none
  !
  !### subroutine to implement linear slip-weakening
  ! friction law for fault dynamics. B.D. 8/19/06
  !...revised for ExGM 100runs. B.D. 8/10/10
  !...revised for SCEC TPV19. B.D. 1/8/12
  !  fricsgl(i,*),i=1 mus, 2 mud, 3 do, 4 cohesion, 
  !  5 time for fixed rutpure, 6 for pore pressure
  !
  real (kind = dp) :: xmu,slip
  real (kind = dp),dimension(6) :: fricsgl
  !
  if(abs(slip).lt.1.0d-10) then
    xmu = fricsgl(1)	!xmu is frictional coefficient, node by node on fault
  elseif(slip < fricsgl(3)) then
    xmu = fricsgl(1) - (fricsgl(1) - fricsgl(2))*slip/fricsgl(3)
  endif
  if(slip >= fricsgl(3)) then
    xmu = fricsgl(2)
  endif
  !
end SUBROUTINE slip_weak

!================================================

SUBROUTINE time_weak(trupt,fricsgl,xmu)
  use globalvar
  implicit none
  !
  !### subroutine to implement linear time-weakening
  ! friction law for fault dynamics. B.D. 8/19/06
  !
  real (kind = dp) :: xmu,trupt
  real (kind = dp),dimension(2) :: fricsgl
  !
  if(trupt <= 0.0d0) then
    xmu = fricsgl(1)
  elseif(trupt < critt0) then
    xmu = fricsgl(1) - (fricsgl(1) - fricsgl(2))*trupt/critt0
  else
    xmu = fricsgl(2)
  endif
  !
end SUBROUTINE time_weak

!2017.2.1 Incorporated RSF into EQQuasi. D.Liu
SUBROUTINE rate_state_ageing_law(V2,theta,fricsgl,xmu,dxmudv)
	use globalvar
	implicit none
	!
	!### subroutine to implement rate- and state- 
	! friction law for fault dynamics. Bin Luo 4/9/2014
	!
	real (kind = dp) :: xmu, dxmudv, V2,theta, A,B,L,f0,V0,tmp, tmpc, phi, sta
	real (kind = dp),dimension(100) :: fricsgl
	!
	A  = fricsgl(9)
	B  = fricsgl(10)
	L  = fricsgl(11)
	f0 = fricsgl(13)
	V0 = fricsgl(12)
	sta = fricsgl(20) ! state variable
	
	tmpc = 1.0d0 / (2.0d0 * V0) * dexp((f0 + B * dlog(V0*theta/L)) / A)
	tmp = (V2) * tmpc
	xmu = A * dlog(tmp + sqrt(tmp**2 + 1.0d0)) !arcsinh(z)= ln(z+sqrt(z^2+1))
	dxmudv = A * tmpc / sqrt(1.0d0 + tmp**2.0d0) ! d(arcsinh(z))/dz = 1/sqrt(1+z^2)
	
	phi = dlog(V0*sta/L)
	if (V2*dtev1/L <= 1.0d-6) then 
		theta = dlog(dexp(phi)*(1-V2*dtev1/L) + V0*dtev1/L)
		theta = L/V0*dexp(theta)
	elseif (V2*dtev1/L > 1.0d-6) then 
		theta = dlog(V0/V2 + (dexp(phi)-V0/V2)*dexp(-V2*dtev1/L)) 
		theta = L/V0*dexp(theta)
	endif

end SUBROUTINE rate_state_ageing_law
!================================================
SUBROUTINE state_evolution_ageing(V2,theta,fricsgl)
  use globalvar
  implicit none
  real (kind = dp) :: V2,theta
  real (kind = dp) :: L
  real (kind = dp),dimension(10) :: fricsgl

  L  = fricsgl(3)

  theta = L/V2 + (theta - L/V2)*dexp(-V2*dt/L)
end SUBROUTINE state_evolution_ageing
!================================================

SUBROUTINE rate_state_slip_law(V2,psi,fricsgl,xmu,dxmudv)
  use globalvar
  implicit none
  !
  !### subroutine to implement rate- and state- 
  ! friction law for fault dynamics. Bin Luo 4/9/2014
  !
  real (kind = dp) :: xmu, dxmudv
  real (kind = dp) :: V2,psi,psiss,fLV,fss
  real (kind = dp) :: A,B,L,f0,V0,fw,Vw
  real (kind = dp),dimension(30) :: fricsgl
  real (kind = dp) :: tmp, tmpc,vold,V1,psi0
  !
  A  = fricsgl(9)
  B  = fricsgl(10)
  L  = fricsgl(11)
  f0 = fricsgl(13)
  V0 = fricsgl(12)
  fw = fricsgl(14)
  Vw = fricsgl(15)
  tmpc = 1.0d0 / (2.0d0 * V0) * dexp(psi/A)
  tmp = (V2+1.d-30) * tmpc
  xmu = A * dlog(tmp + sqrt(tmp**2 + 1.0d0)) !arcsinh(z)= ln(z+sqrt(z^2+1))
  dxmudv = A * tmpc / sqrt(1.0d0 + tmp**2)  ! d(arcsinh(z))/dz = 1/sqrt(1+z^2)
  fLV = f0 - (B - A) * dlog(V2/V0)
  fss = fw + (fLV - fw) / ((1.0d0 + (V2/Vw)**8)**0.125d0)
  psiss = A * dlog(2.0d0 * V0 / V2 * dsinh(fss/A))
  psi = psiss + (psi - psiss) * dexp(-V2*dt/L)
  


end SUBROUTINE rate_state_slip_law

