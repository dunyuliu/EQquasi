SUBROUTINE faulting(ift,nftnd,numnp,neq,constrainm,right,constrain,constrainv,constraina,constrainf,id,&
					n4onf,fltsta,&
					nsmp,fnft,fltslp,un,us,ud,fric,arn,r4nuc,slp4fri,anonfs,nonmx,it,x,&
					kstiff,ia,maxa,right2,itag,dtev1,dtev,status1,constrainvtmp, &
                                        maxsliprate,totmomrate,totruptarea,tottaoruptarea,totslipruptarea)
use globalvar
implicit none
!
!### program to implement faulting boundary conditions.
!	Call this routine after "form_rhs" to correct right-hand-side force
!	vector by taking into fault boundary. Basic coding was done in Feb,
!	2005. Now rewrite. B.D. 7/3/05
! ...totally rewrite this routine in November of 2006 by implementing
!     Day et al. (2005) formulation which does not require to treat
!     slipping, healing seperately. Much nicer! B.D. 11/23/06
! ...extend to 3D case. B.D. 1/26/07
!  
integer (kind=4)::ift,nftnd,numnp,neq,i,i1,j,k,n,isn,imn,n4onf,nonmx,it
real (kind=8) ::slipn,slips,slipd,slip,slipraten,sliprates,sliprated,&
sliprate,xmu,mmast,mslav,mtotl,fnfault,fsfault,fdfault,tnrm,tstk, &
tdip,taox,taoy,taoz,ttao,taoc,ftix,ftiy,ftiz,trupt,tr,&
tmp1,tmp2,tmp3,tmp4,tnrm0,xmu1,xmu2,dtev1,dtev,kapa
integer (kind=4),dimension(3,nonmx) :: anonfs
integer (kind=4),dimension(2,nftnd) :: nsmp
real (kind=8),dimension(nftnd) :: fnft,arn,r4nuc,slp4fri,dtev1D
real (kind=8),dimension(3,nftnd) :: un,us,ud,fltslp,fltslr
real (kind=8),dimension(50,nftnd) :: fric
real (kind=8),dimension(10,nstep,n4onf) :: fltsta
real (kind=8),dimension(6,2,4)::fvd=0.0d0
real (kind=8),dimension(neq) :: right,right2,mass,d,v,acc,cRf
integer (kind=4),dimension(ndof,numnp) :: id
real (kind=8),dimension(ndof,numnp) :: x, constrain, constrainv, constrainvtmp, constraina, constrainm, constrainf
real(kind=8)::alfa,delt,aa0,aa1,aa2,aa3,aa4,aa5,aa6,aa7
real(kind=8)::fa,fb
!RSF from Bin. 2016.08.28
real (kind=8),dimension(nftnd) :: state, A, Vw,aaa,bbb  !RSF
real (kind=8) :: rr,R0,T,F,G,dtao0,dtao,mr!RSF
real (kind=8) :: statetmp, T_coeff!RSF
integer (kind=4) :: iv,ivmax=1000,signal(2)  !RSF
real (kind=8) :: tstk0, tdip0,ttao0, tstk1, tdip1, ttao1, taoc_old, taoc_new !RSF
real (kind=8) :: dxmudv, rsfeq, drsfeqdv, vtmp,caph,caphprime,newtonstep !RSF

real (kind=8) :: frica,fricb,fricl,fricf0,fricv0,fricfw,fricvw,&
	vup,vdo,v_trial,v_trial_new,sliprs_trial,sliprd_trial, &
	tstk_trial,tdip_trial,tau_fem_up,tau_fem_do,tau_fem_trial, &
	tau_fric_up,tau_fric_do,tau_fric_trial,&
	gtrial,gup,gdo,&
	fLV,fss,psiss,psitmp,psi,tmp,xmutmp,&
	v_s_new_mast,v_s_new_slav,v_d_new_mast,v_d_new_slav,&
	uxm,uym,uzm,uxs,uys,uzs,&
	vxm,vym,vzm,vxs,vys,vzs,&
	axm,aym,azm,axs,ays,azs,&
	anm,asm,adm,ans,ass,ads,&
	slipratemast,sliprateslav,&
	slipaccn,slipaccs,slipaccd,a_trial,time0
	
integer(kind=4)::maxa,ia(neq+1),itag,status1
real(kind=8)::kstiff(maxa),acct(neq),phi,signres
real(kind = 8) :: ma_bar_ku_arr(nftnd), sliprate_arr(nftnd),momrate_arr(nftnd),maxsliprate,totmomrate,ruptarea_arr(nftnd),totruptarea,taoruptarea_arr(nftnd), &
                tottaoruptarea,slipruptarea_arr(nftnd),totslipruptarea! calculate abs(ma/ku) for all fault nodes.
!Parameters for Newmark Integration   
	delt=0.5d0  
	alfa=0.25d0*(0.5d0+delt)**2
 	! aa0=1.d0/alfa/dt/dt
	! aa1=delt/alfa/dt
	! aa2=1.d0/alfa/dt
	! aa3=1.d0/alfa*0.5d0-1.d0
	! aa4=delt/alfa-1.d0
	! aa5=dt*(delt/alfa-2.d0)*0.5d0
	! aa6=dt*(1.d0-delt)
	! aa7=delt*dt	       
	aa0=1.d0/alfa/dtev1/dtev1
	aa1=delt/alfa/dtev1
	aa2=1.d0/alfa/dtev1
	aa3=1.d0/alfa*0.5d0-1.d0
	aa4=delt/alfa-1.
	aa5=dtev1*(delt/alfa-2.d0)*0.5d0
	aa6=dtev1*(1.d0-delt)
	aa7=delt*dtev1
!...do not use OpenMP for better output for ExGM 100runs. B.D. 8/12/10
!*** loop over slave nodes ***
!!$omp parallel do default(shared) private(i,j,k,fnfault,fsfault,fdfault,isn,imn,fvd, &
!!$omp	tmp1,tmp2,tmp3,tmp4,slipn,slips,slipd,slip,slipraten,sliprates,sliprated,sliprate, & 
!!$omp	mslav,mmast,mtotl,tnrm,tstk,tdip,ttao,taoc,taox,taoy,taoz,ftix,ftiy,ftiz,&
!!$omp	xmu,trupt)
!...the above definition of private is very important. xmu was not defined as private
!	ealier and resulted in problems: it can be imaged that it should if different
!	OpenMP threads mess up xmu! B.D. 10/31/09
ruptarea_arr = 0.0d0
taoruptarea_arr = 0.0d0 
slipruptarea_arr = 0.0d0
do i=1,nftnd	!just fault nodes
    fnfault = fric(7,i) !initial forces on the fault node
    fsfault = fric(8,i) !norm, strike, dip components directly
    fdfault = 0.0d0

	isn = nsmp(1,i)
	imn = nsmp(2,i)	
	!if (x(1,isn)==-18.0d3.and.x(3,isn)==-21.0d3)then
	!	write(*,*) isn,imn
	!endif
	!...get nodal force,velocity, and diplacement in x,y,z.
	!   B.D. 1/26/07
	!...aslo add Rayleigh stiffness damping before using d.
	!   assume stifness coefficient is first material: rdampk(1).
	!   B.D. 11/26/06
	!...it seems the damping should not be used here.
	!   B.D. 1/28/07
	do j=1,2  !1-slave, 2-master
		do k=1,3  !1-x comp, 2-y comp, 3-z comp
			fvd(k,j,1) = constrainf(k,nsmp(j,i))  !1-force !DL 
			fvd(k,j,2) = constrainv(k,nsmp(j,i)) !2-vel
			fvd(k,j,3) = constrain(k,nsmp(j,i)) !3-disp
			fvd(k,j,4) = constraina(k,nsmp(j,i)) !4-acc
		enddo
	enddo
	!...resolve x,y,z components onto normal, strike and dip components.
	!   B.D. 1/26/07
	do j=1,4    !1-force,2-vel,3-disp,4-acc
		do k=1,2  !1-slave,2-master
			fvd(4,k,j) = fvd(1,k,j)*un(1,i) + fvd(2,k,j)*un(2,i) + fvd(3,k,j)*un(3,i)  !4-norm
			fvd(5,k,j) = fvd(1,k,j)*us(1,i) + fvd(2,k,j)*us(2,i) + fvd(3,k,j)*us(3,i)  !5-strike
			fvd(6,k,j) = fvd(1,k,j)*ud(1,i) + fvd(2,k,j)*ud(2,i) + fvd(3,k,j)*ud(3,i)  !6-dip
		enddo
	enddo	
!
	!...no opening and no penetrating constraints for SCEC TPV12/13. B.D. 10/22/09
	!   should apply before trial traction.
	!   use the average for normal v and d. B.D. 10/22/09
	!Actually, my implementation of Day et al. (2005) here may already take care
	!of this and an extra constraint caused incorrect slip (rate) at surface fault
	!station in TPV13 (dip-slip). So should not apply again here. B.D. 11/2/09
	!tmp1 = fvd(4,1,2)
	!tmp2 = fvd(4,2,2)
	!tmp3 = fvd(4,1,3)
	!tmp4 = fvd(4,2,3)
	!fvd(4,1,2) = 0.5*(tmp1+tmp2)
	!fvd(4,2,2) = fvd(4,1,2)
	!fvd(4,1,3) = 0.5*(tmp3+tmp4)
	!fvd(4,2,3) = fvd(4,1,3)
	!...slip and slip rate should be calculated from norm, strike, and dip so that
	!    the above no opening or penetrating constraints included. B.D. 10/22/09
	slipn = fvd(4,2,3) - fvd(4,1,3)
	slips = fvd(5,2,3) - fvd(5,1,3)
	slipd = fvd(6,2,3) - fvd(6,1,3)
	slip = sqrt(slipn**2 + slips**2 + slipd**2) !slip mag
	fltslp(1,i) = slips  !save for final slip output
	fltslp(2,i) = slipd
	!fltslp(3,i) = slipn  !normal should be zero, but still keep to ensure
	slipraten = fvd(4,2,2) - fvd(4,1,2)
	sliprates = fvd(5,2,2) - fvd(5,1,2)
	sliprated = fvd(6,2,2) - fvd(6,1,2)
	fltslr(1,i) = sliprates  !save for final slip output
	fltslr(2,i) = sliprated
	sliprate = sqrt(slipraten**2+sliprates**2+sliprated**2)
        sliprate_arr(i) = sliprate
	if (sliprate>fltslp(3,i)) then 
		fltslp(3,i)=sliprate
	endif
	slipaccn=fvd(4,2,4)-fvd(4,1,4)
	slipaccs=fvd(5,2,4)-fvd(5,1,4)
	slipaccd=fvd(6,2,4)-fvd(6,1,4)
	!...path-itegrated slip for slip-weakening. B.D. 8/12/10
	slp4fri(i) = slp4fri(i) + sliprate * dtev1
	!...calculate moment rate and moment if needed. B.D. 8/11/10
	!  also, max slip rate for early termination.
	!...for homogeneous material, i.e., only myshr(1). B.D. 1/3/12
	! or for heterogeneous case, but mushr(1) for rupture fault.
	! if(time > term-dt) then  !only at the end, do this for LVFZ3D Plastic.
		! momnt = momnt + miuonf(i) * arn4m(i) * slip
		! momntrat = momntrat + miuonf(i) * arn4m(i) *sliprate
		! if(sliprate>maxslprat) maxslprat=sliprate
	! endif
	!
	!...nodal mass. Mass of each element may not be distributed among its 
	! nodes evenly. Instead, distribution is related to element shape. 
	!   Note: nodal mass should not be directly obtained from left-hand-side
	! diagnoal mass matrix, because that's the effective mass, which takes 
	! damping coefficient into accout. Instead, I computed nodal mass from 
	! element mass and assembled in "qdct2.f90".B.D.7/3/05
	mslav = constrainm(1,isn)!mass(id(1,isn))		
	mmast = constrainm(1,imn)!mass(id(1,imn))
	mtotl = mslav + mmast
	if (mmast==0.0d0.or.mslav==0.0d0) then 
		write(*,*) 'ZERO MASS IN FAULTING, MMAST & MSLAV',mmast,mslav
		write(*,*) 'PROBLEMATIC COORDS ARE, X,Z =',x(1,isn)/1.0d3,x(3,isn)/1.0d3
		stop 501
	endif 
	!
	!...trial traction to enforce continuity. B.D. 11/23/06
	!...divided by the associated area to get traction from force for EQdyna3d v2.1.2.
	!   initial stress, not initial force (f*fault) used here. B.D. 2/28/08
	!...no fault initial stress in elastoplastic rheology. B.D. 1/8/12
	mtotl = mtotl * arn(i)
	if (IS==1) then
		tnrm = (mslav*mmast*((fvd(4,2,2)-fvd(4,1,2))+(fvd(4,2,3)-fvd(4,1,3))/dt)/dt &
			+ mslav*fvd(4,2,1) - mmast*fvd(4,1,1)) / mtotl + fnfault         
		tstk = (mslav*mmast*(fvd(5,2,2)-fvd(5,1,2))/dt + mslav*fvd(5,2,1) &
			- mmast*fvd(5,1,1)) / mtotl + fsfault
		tdip = (mslav*mmast*(fvd(6,2,2)-fvd(6,1,2))/dt + mslav*fvd(6,2,1) &
			- mmast*fvd(6,1,1)) / mtotl + fdfault
	elseif (IS==2) then
		tnrm=(mslav*mmast*(aa0*(fvd(4,2,3)-fvd(4,1,3))+(fvd(4,2,2)-fvd(4,1,2))*aa2+aa3*(fvd(4,2,4)-fvd(4,1,4)))+mslav*fvd(4,2,1)-mmast*fvd(4,1,1))/mtotl+fnfault
		tstk=(mslav*mmast/aa7*((fvd(5,2,2)-fvd(5,1,2))+(fvd(5,2,4)-fvd(5,1,4))*aa6)+mslav*fvd(5,2,1)-mmast*fvd(5,1,1))/mtotl+fsfault
		tdip=(mslav*mmast/aa7*((fvd(6,2,2)-fvd(6,1,2))+(fvd(6,2,4)-fvd(6,1,4))*aa6)+mslav*fvd(6,2,1)-mmast*fvd(6,1,1))/mtotl+fdfault
	endif	
	ttao = sqrt(tstk*tstk + tdip*tdip) !total shear magnitude	    
	!
	!...friction law to determine friction coefficient
	!   slip-weakening only so far. B.D. 1/26/07
	!... based on choices, call corresponding friction laws.
	! B.D. 10/8/08
if (friclaw==1.or.friclaw==2)then!Differ 1&2 and 3&4	
!PART2:SLIP-WEAKENING
	if(friclaw == 1) then
		call slip_weak(slp4fri(i),fric(1,i),xmu)
	elseif(friclaw == 2) then
		trupt =  time - fnft(i)
		call time_weak(trupt,fric(1,i),xmu)
	endif
	! !......for nucleation zone of the nucleation fault,which initiates rupture,
	! !	rupture propagates at a fixed speed to drop "xmu". B.D. 8/31/06
	! if(ift == nucfault .and. xmu > fric(2,i)) then	
		! !only nucleation fault and before finishing dropping, do...
		! if(r4nuc(i) <= srcrad0) then !only within nucleation zone, do...
			! tr = r4nuc(i) / vrupt0
			! if(tr <= time) then !only ready or already fail, do...
				! trupt = time - tr
			! call time_weak(trupt,fric(1,i),xmu)
			! endif
		! endif
      ! !only nucleation fault and before finishing dropping, do...
      ! ! if(r4nuc(i) <= srcrad0) then !only within nucleation zone, do...
        ! ! tr=(r4nuc(i)+0.081*srcrad0*(1./(1-(r4nuc(i)/srcrad0)* &
         ! ! (r4nuc(i)/srcrad0))-1))/(0.7*3464.)
        ! ! tr = r4nuc(i) / vrupt0
        ! !if(tr <= time) then !only ready or already fail, do...
          ! ! trupt = time - tr
          ! ! call time_weak(trupt,fric(1,i),xmu1)
        ! ! call slip_weak(slp4fri(i),fric(1,i),xmu2)
        ! ! xmu=min(xmu1,xmu2)  !minimum friction used. B.D. 2/16/13
        ! ! endif
      ! !endif		
	! endif
	if (C_Nuclea==1) then	
		if(r4nuc(i)<=srcrad0) then !only within nucleation zone, do...
			tr=(r4nuc(i)+0.081d0*srcrad0*(1.d0/(1.d0-(r4nuc(i)/srcrad0)*(r4nuc(i)/srcrad0))-1.d0))/(0.7d0*3464.d0)
		else
			tr=1.0d9 
		endif
		if(time<tr) then 
			fb=0.0d0
		elseif ((time<(tr+critt0)).and.(time>=tr)) then 
			fb=(time-tr)/critt0
		else 
			fb=1.0d0
		endif
		tmp1=fric(1,i)+(fric(2,i)-fric(1,i))*fb
		tmp2=xmu
		xmu=min(tmp1,tmp2)  !minimum friction used. B.D. 2/16/13	
	endif
	!...adjust tstk,tdip and tnrm based on jump conditions on fault.
	!   before calculate taoc, first adjust tnrm if needed. 
	!   after this, they are true (corrected) values. B.D. 11/23/06
	!   cohesion is added here. B.D. 2/28/08
	!...for SCEC TPV10/11, no opening is allowed. B.D. 11/24/08
	!if(tnrm > 0) tnrm = 0   !norm must be <= 0, otherwise no adjust
	!taoc = cohes - xmu * tnrm
	!...for ExGM 100 runs, no opening allowed means following.
	!  B.D. 8/12/10
	if((tnrm+fric(6,i))>0.d0) then
		tnrm0 = 0.0d0
	else
		tnrm0 = tnrm+fric(6,i)
	endif
	taoc = fric(4,i) - xmu *tnrm0
	!taoc = cohes - xmu * tnrm0
	!if(tnrm > 0) tnrm = 0   !norm must be <= 0, otherwise no adjust
	!taoc = fistr(5,i) - xmu * tnrm
	if(ttao > taoc) then
		tstk = tstk * taoc / ttao
		tdip = tdip * taoc / ttao
		if(fnft(i)>600.d0) then	!fnft should be initialized by >10000
			if(sliprate >= 0.001d0) then	!first time to reach 1mm/s
				fnft(i) = time	!rupture time for the node
			endif
		endif
	endif
	!
	!...add the above fault boundary force and initial force to elastic
	!	force of the split nodes. 
	!   first resolve normal, strike and dip back to x-,y-,z-. 
	!   then subtract them from slave, add to master as the above calculation
	!   based on this convention. see Day et al. (2005). B.D. 11/23/06
	!...due to traction, not force used in friction law above, need area to 
	!   convert traction to force for v2.1.2. B.D. 2/28/08
	taox = (tnrm0*un(1,i) + tstk*us(1,i) + tdip*ud(1,i))*arn(i)
	taoy = (tnrm0*un(2,i) + tstk*us(2,i) + tdip*ud(2,i))*arn(i)
	taoz = (tnrm0*un(3,i) + tstk*us(3,i) + tdip*ud(3,i))*arn(i)
!*.*Sep.13.2015/D.L.	
    ftix = (fnfault*un(1,i) + fsfault*us(1,i) + fdfault*ud(1,i))*arn(i)
    ftiy = (fnfault*un(2,i) + fsfault*us(2,i) + fdfault*ud(2,i))*arn(i)
    ftiz = (fnfault*un(3,i) + fsfault*us(3,i) + fdfault*ud(3,i))*arn(i)  
	if (IS==1)then 
	right(id(1,isn)) = right(id(1,isn)) + taox - ftix
	right(id(2,isn)) = right(id(2,isn)) + taoy - ftiy
	right(id(3,isn)) = right(id(3,isn)) + taoz - ftiz
	right(id(1,imn)) = right(id(1,imn)) - taox + ftix
	right(id(2,imn)) = right(id(2,imn)) - taoy + ftiy
	right(id(3,imn)) = right(id(3,imn)) - taoz + ftiz
	elseif (IS==2)then
	right(id(1,isn)) = right(id(1,isn)) + taox - ftix
	right(id(2,isn)) = right(id(2,isn)) + taoy - ftiy
	right(id(3,isn)) = right(id(3,isn)) + taoz - ftiz
	right(id(1,imn)) = right(id(1,imn)) - taox + ftix
	right(id(2,imn)) = right(id(2,imn)) - taoy + ftiy
	right(id(3,imn)) = right(id(3,imn)) - taoz + ftiz
		! axs=(right(id(1,isn))+taox-ftix)/mass(id(1,isn))
		! ays=(right(id(2,isn))+taoy-ftiy)/mass(id(2,isn))
		! azs=(right(id(3,isn))+taoz-ftiz)/mass(id(3,isn))
		! d(id(1,isn))=d(id(1,isn))+v(id(1,isn))*dt+(axs*alfa+(0.5-alfa)*acc(id(1,isn)))*dt*dt
		! d(id(2,isn))=d(id(2,isn))+v(id(2,isn))*dt+(ays*alfa+(0.5-alfa)*acc(id(2,isn)))*dt*dt
		! d(id(3,isn))=d(id(3,isn))+v(id(3,isn))*dt+(azs*alfa+(0.5-alfa)*acc(id(3,isn)))*dt*dt
		
		! axm=(right(id(1,imn))-taox+ftix)/mass(id(1,imn))
		! aym=(right(id(2,imn))-taoy+ftiy)/mass(id(2,imn))
		! azm=(right(id(3,imn))-taoz+ftiz)/mass(id(3,imn))			
		! d(id(1,imn))=d(id(1,imn))+v(id(1,imn))*dt+(axm*alfa+(0.5-alfa)*acc(id(1,imn)))*dt*dt
		! d(id(2,imn))=d(id(2,imn))+v(id(2,imn))*dt+(aym*alfa+(0.5-alfa)*acc(id(2,imn)))*dt*dt
		! d(id(3,imn))=d(id(3,imn))+v(id(3,imn))*dt+(azm*alfa+(0.5-alfa)*acc(id(3,imn)))*dt*dt		
	endif
elseif (friclaw==3.or.friclaw==4)then 	
!PART3: RSF
!---3.1: ARTIFICIAL NUCLEATION FOR DYNAMIC PROCESS
	!dtao=0.0d0
	! time0=time
    ! if(C_Nuclea==1)then
       ! R0 = 3000.0d0
       ! dtao0 = 20.0d6
       ! T = 1.0d0
       ! F = 0.0d0
       ! rr=sqrt((x(1,nsmp(1,i))-xsource)**2+(x(3,nsmp(1,i))-zsource)**2)    
       ! if (rr<R0)    F=dexp(rr**2/(rr**2-R0**2))
       ! G = 1.0d0
       ! if (time0<=T)  G=dexp((time0-T)**2/(time0*(time0-2*T)))
       ! dtao=dtao0*F*G
    ! endif
	!fsfault = fric(8,i)+dtao
!---3.2: DECLARE THE TIME FOR RUPTURING.	
	if(fnft(i)<0.0d0) then	!fnft should be initialized by >10000
		if(sliprate >= 0.001d0) then	!first time to reach 1mm/s
			fnft(i) = time	!rupture time for the node
		endif
	endif	
!---3.3: INITIATE [V_TRIAL],[SLIPRS_TRIAL],,[SLIPRD_TRIAL]
	v_trial=sliprate
	sliprs_trial=sliprates
	sliprd_trial=sliprated
	statetmp=fric(20,i)
!---3.4: FRICTIONAL/NON FRIVTIONAL REGION	
	if 	(abs(x(3,isn)--60.0d3) <= 40.0d3) then 
!---3.4.1: FRICTIONAL REGION CONTROLLED BY RSF. 
		if (status1 == 0.or.status1 == 1) then 
!---3.4.1.2: Dynamic + NON-DYNAMIC PROCESS: [STATUS1]==0		
			tstk0 = (mslav * fvd(5,2,1) - mmast * fvd(5,1,1)) / mtotl + fsfault
			tdip0 = (mslav * fvd(6,2,1) - mmast * fvd(6,1,1)) / mtotl + fdfault	
			tnrm = (mslav * fvd(4,2,1) - mmast * fvd(4,1,1)) / mtotl + fnfault
            ! fric(41,i) is abs(KU)
			fric(41,i) = sqrt((mslav * fvd(5,2,1) - mmast * fvd(5,1,1))**2 + (mslav * fvd(6,2,1) - mmast * fvd(6,1,1))**2) / (mmast + mslav) 
            ttao0 = sqrt(tstk0 * tstk0 + tdip0 * tdip0)		
			if((tnrm + fric(6,i))>-10.0d6) then
				tnrm0 = -10.0d6
			else
				tnrm0 = tnrm + fric(6,i)
			endif	
			!theta = L/V + (theta - L/V)*dexp(-V*dtev1/L)
			!- fric(22): theta*
			!- fric(21): theta_dot
			!- [fric(11): L],[fric(12,i): V0]
!---3.4.1.3: TO COMPUTE THETA*(t+1) [FRIC(22,:)] FOR ITAG==0, USING [V_TRIAL]=[CONSTRAINV] 
!----------: THETA**(t+1) [FRIC(22,:)] FOR ITAG==1, , USING [V_TRIAL]={[CONSTRAINV]+[CONSTRAINVTMP]}*0.5	
			if (itag == 0) then 
				v_trial = sliprate
                fric(42,i) = sliprate !Only record the final sliprate in fric(42,i) in the last time step.
                phi = dlog(fric(12,i) * fric(20,i) / fric(11,i))
				if (v_trial * dtev1 / fric(11,i) <= 1.0d-6) then 
					fric(22,i) = dlog(dexp(phi)*(1-v_trial*dtev1/fric(11,i)) + fric(12,i)*dtev1/fric(11,i))
					fric(22,i) = fric(11,i)/fric(12,i)*dexp(fric(22,i))
				elseif (v_trial * dtev1 / fric(11,i) > 1.0d-6) then 
					fric(22,i) = dlog(fric(12,i)/v_trial + (dexp(phi)-fric(12,i)/v_trial)*dexp(-v_trial*dtev1/fric(11,i))) 
					fric(22,i) = fric(11,i)/fric(12,i)*dexp(fric(22,i))
				endif
			elseif (itag == 1) then 
				do j=1,2  !1-slave, 2-master
					do k=1,3  !1-x comp, 2-y comp, 3-z comp 
						fvd(k,j,2) = constrainvtmp(k,nsmp(j,i)) !2-vel
					enddo
				enddo	
				do k=1,2  !1-slave,2-master
					fvd(4,k,2) = fvd(1,k,2)*un(1,i) + fvd(2,k,2)*un(2,i) + fvd(3,k,2)*un(3,i)  !4-norm
					fvd(5,k,2) = fvd(1,k,2)*us(1,i) + fvd(2,k,2)*us(2,i) + fvd(3,k,2)*us(3,i)  !5-strike
					fvd(6,k,2) = fvd(1,k,2)*ud(1,i) + fvd(2,k,2)*ud(2,i) + fvd(3,k,2)*ud(3,i)  !6-dip
				enddo	
				slipraten = fvd(4,2,2) - fvd(4,1,2)
				sliprates = fvd(5,2,2) - fvd(5,1,2)
				sliprated = fvd(6,2,2) - fvd(6,1,2)
				v_trial = sqrt(slipraten**2+sliprates**2+sliprated**2)		
				phi = dlog(fric(12,i) * fric(20,i) / fric(11,i))
				if (v_trial * dtev1 / fric(11,i) <= 1.0d-6) then 
					fric(22,i) = dlog(dexp(phi)*(1-v_trial*dtev1/fric(11,i)) + fric(12,i)*dtev1/fric(11,i))
					fric(22,i) = fric(11,i)/fric(12,i)*dexp(fric(22,i))
				elseif (v_trial * dtev1 / fric(11,i) > 1.0d-6) then 
					fric(22,i) = dlog(fric(12,i)/v_trial + (dexp(phi)-fric(12,i)/v_trial)*dexp(-v_trial*dtev1/fric(11,i))) 
					fric(22,i) = fric(11,i)/fric(12,i)*dexp(fric(22,i))
				endif		
			endif
!---3.4.1.4: NEWTON METHOD TO GET NEW V_TRIAL BASED ON TRACTION=ELASTIC FORCE.			
			if (x(1,isn)==-20.0d3.and.x(3,isn)==-70.0d3) then 
				write(*,*) 'Before Newton, ITAG',itag,'V_trial',v_trial
                                write(*,*) 'State Variable',fric(22,i)
				write(*,*) 'tstk0,tdip0,tnrm0',tstk0/1d6,tdip0/1d6,tnrm0/1d6
				write(*,*) 'fvd511,521',fvd(5,1,1),fvd(5,2,1)
				write(*,*) 'mmast,mslav',mmast,mslav
			endif 
			!Modify ttao0 for quasi-dynamic. 11252019.
			!ttao0 = ttao0 - v_trial*2670.0d0*3464.0d0/2.0d0
                        !newtonstep = 1.0d0
			do iv = 1,ivmax	
				if(friclaw == 3) then
					call rate_state_ageing_law(v_trial,fric(22,i),fric(1,i),xmu,dxmudv,dtev1) !RSF
				elseif(friclaw == 4) then
					call rate_state_slip_law(v_trial,fric(22,i),fric(1,i),xmu,dxmudv) !RSF
				endif 	
				tau_fric_trial = fric(4,i) - xmu * tnrm0
				
				rsfeq = (tau_fric_trial-ttao0 + 3464.d0*2670.d0/2.0d0*v_trial)
				drsfeqdv = (-dxmudv * tnrm0 +  3464.d0*2670.d0/2.0d0)
                                
			!if(abs(rsfeq/drsfeqdv) < 1.d-14 * abs(v_trial).and. abs(rsfeq) < 1.d-6 * abs(v_trial)) exit 
				if (abs(rsfeq) < 1.0d-7 * ttao0) exit
				vtmp = v_trial -  rsfeq / drsfeqdv
				if(vtmp <= 0.0d0) then
				  v_trial = v_trial/1.1d0
                                  !newtonstep = newtonstep/2.0d0
			        else
				  v_trial = vtmp
                                  !newtonstep = 1.0d0
				endif
                                signres = rsfeq
                                
				if (x(1,isn)==-20d3.and.x(3,isn)==-70d3)then
					write(*,*) 'FAULTING=======AAA======ITAG=',itag
					write(*,*) 'iv,vtrial,phi*',iv,v_trial,fric(22,i)
                                        write(*,*) 'rsfeq',rsfeq,'drs',drsfeqdv
					write(*,*) 'fvd,mast,slav',fvd(5,2,1),fvd(5,1,1)
					write(*,*) 'tau_fric_trial,ttao0',tau_fric_trial/1.0d6,ttao0/1.0e6
					write(*,*) 'tstk0,tdip0',tstk0/1.0d6,tdip0/1.0e6
				endif	
				!if (x(1,isn)==-12d3.and.x(3,isn)==-21d3)then
				!	write(*,*) 'FAULTING=======BBB======ITAG=',itag
				!	write(*,*) 'iv,vtrial,phi*',iv,v_trial,fric(22,i)
				!	write(*,*) 'fvd,mast,slav',fvd(5,2,1),fvd(5,1,1)
				!	write(*,*) 'tau_fric_trial,ttao0',tau_fric_trial/1.0d6,ttao0/1.0e6
				!	write(*,*) 'tstk0,tdip0',tstk0/1.0d6,tdip0/1.0e6
				!endif					
			enddo !iv	
			if(v_trial < fric(19,i)) then
				v_trial = fric(19,i)
				tau_fric_trial = ttao0
			endif	
			
			tstk=tau_fric_trial*tstk0/ttao0
			tdip=tau_fric_trial*tdip0/ttao0
			ttao=sqrt(tstk**2+tdip**2)
			fric(26,i)=v_trial
			fric(28,i)=tstk
			fric(29,i)=tdip
			fric(30,i)=tnrm0
			slipratemast=(v_trial)*mslav/(mmast+mslav)
			sliprateslav=-(v_trial)*mmast/(mmast+mslav)

			v_s_new_mast=slipratemast*tstk/ttao
			v_d_new_mast=slipratemast*tdip/ttao
			v_s_new_slav=sliprateslav*tstk/ttao
			v_d_new_slav=sliprateslav*tdip/ttao
			vxm=v_s_new_mast*us(1,i)+v_d_new_mast*ud(1,i)
			vym=v_s_new_mast*us(2,i)+v_d_new_mast*ud(2,i)
			vzm=v_s_new_mast*us(3,i)+v_d_new_mast*ud(3,i)
			vxs=v_s_new_slav*us(1,i)+v_d_new_slav*ud(1,i)
			vys=v_s_new_slav*us(2,i)+v_d_new_slav*ud(2,i)
			vzs=v_s_new_slav*us(3,i)+v_d_new_slav*ud(3,i)
!---3.4.1.5: V* FOR ITAG==0 AND V** FOR ITAG==1 OBTAINED. 			
			if (itag == 0) then 		
				constrainvtmp(1,imn)=vxm
				constrainvtmp(2,imn)=vym 
				constrainvtmp(3,imn)=vzm 
				constrainvtmp(1,isn)=vxs
				constrainvtmp(2,isn)=vys 
				constrainvtmp(3,isn)=vzs 	
			elseif (itag == 1) then 
				fric(20,i) = fric(22,i)
                                ma_bar_ku_arr(i) = (v_trial - fric(42,i)) / dtev1 * mmast * mslav / (mmast + mslav) / fric(41,i)
                                ma_bar_ku_arr(i) = abs(ma_bar_ku_arr(i))
                                momrate_arr(i) =0.0d0 
                                if ((abs(x(1,isn))<=30.0d3+6.0d3).and.(abs(x(3,isn)--60.0d3)<=15.0d3+6.0d3)) then
                                        momrate_arr(i) = 3464.0d0**2*2670.0d0*v_trial*dx*dx
                                endif
                                ruptarea_arr(i) = 0.0d0
                                taoruptarea_arr(i) = 0.0d0
                                slipruptarea_arr(i) = 0.0d0
                                if (v_trial>=1.0d-3) then
                                        ruptarea_arr(i) = dx*dx
                                        taoruptarea_arr(i) = ttao*dx*dx
                                        slipruptarea_arr(i) = slip*dx*dx
                                endif
                                
				constrainv(1,imn)=vxm
				constrainv(2,imn)=vym 
				constrainv(3,imn)=vzm 
				constrainv(1,isn)=vxs
				constrainv(2,isn)=vys 
				constrainv(3,isn)=vzs 		
				fric(31,i) = vxm 
				fric(32,i) = vym 
				fric(33,i) = vzm 
				fric(34,i) = vxs
				fric(35,i) = vys 
				fric(36,i) = vzs 
			endif 
!---3.4.1.6: WHEN ITAG==0, SIMPLY STORE V* INTO [CONSTRAINVTMP]
!----------: WHEN ITAG==1, DECLARE [FRIC(22,:)] AND FINAL V** INTO [CONSTRAINV]			
		endif!status1 == 0 / 1
!---3.4.1.7: CHECKING BLOCK		
		if (isnan(sliprate).or.isnan(v_trial)) then 
			write(*,*) 'LOCATION:x,z',x(1,isn)/1.0d3,x(3,isn)/1.0d3
			write(*,*) 'sliprate =',sliprate,'v_trial',v_trial
			stop '"sliprate/v_trial" is a NaN'
		endif
		if (isnan(ttao/1.0d6)) then 
			write(*,*) 'iv,T_coeff',iv,T_coeff
			write(*,*) 'sliprates,sliprated,v_trial',sliprates,sliprated,v_trial
			write(*,*) 'Tstk0,tdip0',tstk,tdip
			write(*,*) 'x,z',x(1,isn)/1.0d3,x(3,isn)/1.0d3
			stop '"ttao" is a NaN'
		endif		
		if (isnan(fric(22,i))) then 
			
			write(*,*) 'x,z',x(1,isn)/1.0d3,x(3,isn)/1.0d3
			stop '"psi" is a NaN'
		endif		
	 elseif (abs(x(3,isn)--60.0d3)>40.0d3) then
! !---3.4.2: LOADING BOTTOM AT A FIXED SLIDING RATE	
		 tstk0=2.585534683723515d7
		 tdip0=0.0d6
		 tnrm=-50.0d6	 
		 if((tnrm+fric(6,i))>0) then
			 tnrm0 = 0.0d0
		 else
			 tnrm0 = tnrm+fric(6,i)
		 endif
		 fric(28,i)=tstk0
		 fric(29,i)=tdip0
		 v_trial = 1.0d-9
		 slipratemast=(v_trial)*mslav/(mmast+mslav)
		 sliprateslav=-(v_trial)*mslav/(mmast+mslav)
		 v_s_new_mast=slipratemast
		 v_d_new_mast=0.0d0
		 v_s_new_slav=sliprateslav
		 v_d_new_slav=0.0d0
		 vxm=v_s_new_mast*us(1,i)+v_d_new_mast*ud(1,i)
		 vym=v_s_new_mast*us(2,i)+v_d_new_mast*ud(2,i)
		 vzm=v_s_new_mast*us(3,i)+v_d_new_mast*ud(3,i)
		 vxs=v_s_new_slav*us(1,i)+v_d_new_slav*ud(1,i)
		 vys=v_s_new_slav*us(2,i)+v_d_new_slav*ud(2,i)
		 vzs=v_s_new_slav*us(3,i)+v_d_new_slav*ud(3,i)
		 
		 constrainvtmp(1,imn)=vxm
		 constrainvtmp(2,imn)=vym 
		 constrainvtmp(3,imn)=vzm 
		 constrainvtmp(1,isn)=vxs
		 constrainvtmp(2,isn)=vys 
		 constrainvtmp(3,isn)=vzs 
		
		 constrainv(1,imn)=vxm
		 constrainv(2,imn)=vym 
		 constrainv(3,imn)=vzm 
		 constrainv(1,isn)=vxs
		 constrainv(2,isn)=vys 
		 constrainv(3,isn)=vzs 	
	
	     !v_trial = 1.0d-20! reset v_trial to deactivate its effect in variable time step.
	 endif	
	!v(t+1)=v(t)+a6*a(t)+a7*a(t+1)====>>a(t+1)=(v(t+1)-v(t)-a6*a(t))/a7 
	!a(t+1)=a0*(u(t+1)-u(t))-a2*v(t)-a3*a(t)====>>u(t+1)=(a(t+1)+a3*a(t)+a2*v(t))/a0+u(t)	

	!Varaible time step 
	kapa=0.25d0*(3.0d10*fric(11,i)/fric(9,i)/abs(fnfault)-(fric(10,i)-fric(9,i))/fric(9,i))**2.0d0-3.0d10*fric(11,i)/fric(9,i)/abs(fnfault)

	dtev1D(i)=ksi*fric(11,i)/v_trial
	if (dtev1D(i) < 0) then
		write(*,*) 'NEGATIVE SLIPRATE, ITS LOC = ', x(1,isn),x(3,isn)
		write(*,*) 'PROBLEMATIC V_TRIAL = ', v_trial
		stop 502
	endif
		! if (x(1,isn)==-12.15d3.and.x(3,isn)==-11.70d3)then
			! write(*,*) 'STATION BBB'
			! write(*,*) 'sliprate',fvd(5,2,2),fvd(5,1,2)
			! write(*,*) 'slipacc',fvd(5,2,4),fvd(5,1,4)
			! write(*,*) 'Force',fvd(5,2,1),fvd(5,1,1)			
			! write(*,*) 'v_trial,psi,tstk',v_trial,psi,tstk/1e6
		! endif		
		! if (x(1,isn)==-9.15d3.and.x(3,isn)==-11.70d3)then
			! write(*,*) 'STATION AAA'
			! write(*,*) 'sliprate',fvd(5,2,2),fvd(5,1,2)
			! write(*,*) 'slipacc',fvd(5,2,4),fvd(5,1,4)
			! write(*,*) 'Force',fvd(5,2,1),fvd(5,1,1)			
			! write(*,*) 'v_trial,psi,tstk',v_trial,psi,tstk/1e6
		! endif			
	!				
	if (IS==1)then 
	right(id(1,isn)) = right(id(1,isn)) + taox - ftix
	right(id(2,isn)) = right(id(2,isn)) + taoy - ftiy
	right(id(3,isn)) = right(id(3,isn)) + taoz - ftiz
	right(id(1,imn)) = right(id(1,imn)) - taox + ftix
	right(id(2,imn)) = right(id(2,imn)) - taoy + ftiy
	right(id(3,imn)) = right(id(3,imn)) - taoz + ftiz
	elseif (IS==2)then
	!fric(20,i)=statetmp
	! if (dtev1>dt) then
		! right(id(1,isn)) = right(id(1,isn)) - mmast * axs
		! right(id(2,isn)) = right(id(2,isn)) - mmast * ays
		! right(id(3,isn)) = right(id(3,isn)) - mmast * azs
		! right(id(1,imn)) = right(id(1,imn)) - mmast * axm
		! right(id(2,imn)) = right(id(2,imn)) - mmast * aym
		! right(id(3,imn)) = right(id(3,imn)) - mmast * azm
	! endif		
	! right(id(1,isn)) = uxs*kstiff(ia(id(1,isn)))
	! right(id(2,isn)) = uys*kstiff(ia(id(2,isn)))
	! right(id(3,isn)) = uzs*kstiff(ia(id(3,isn)))
	! right(id(1,imn)) = uxm*kstiff(ia(id(1,imn)))
	! right(id(2,imn)) = uym*kstiff(ia(id(2,imn)))
	! right(id(3,imn)) = uzm*kstiff(ia(id(3,imn)))
	endif	
endif
	!...Store fault forces and slip/slipvel for fault nodes 
	!		at set time interval.
	! note: forces will be transferred to stress later
	! B.D. 8/21/05
	!...now, they are directly traction (stress) in version 2.1.2.
	!   and can be written out here. B.D. 2/28/08
	!...Store only, no write out. B.D. 10/25/09
	if ((status1==0.and.itag==1).or.(status1==1)) then
		if(n4onf>0) then
			do j=1,n4onf
				if(anonfs(1,j)==i.and.anonfs(3,j)==ift) then !only selected stations. B.D. 10/25/09    
					fltsta(1,it,j) = time
					fltsta(2,it,j) = v_trial
					fltsta(3,it,j) = sliprated
					fltsta(4,it,j) = fric(20,i)
					fltsta(5,it,j) = slips
					fltsta(6,it,j) = slipd
					fltsta(7,it,j) = slipn
					fltsta(8,it,j) = fric(28,i)!tstk
					fltsta(9,it,j) = fric(29,i)!tdip
					fltsta(10,it,j) = tnrm0
				endif
			enddo 
		endif   
	endif
	!if (x(1,isn)==xsource.and.x(3,isn)==zsource)then
		!write(*,*) '************Sliprates',sliprates,sliprate
		!write(*,*) '&&&&&&&&&&&&Traction',tnrm,tstk
		!write(*,*) 'Force',taox,ftix,xmu
		!write(*,*) 'dacc',(fvd(4,2,4)-fvd(4,1,4)),(fvd(4,2,2)-fvd(4,1,2)),(fvd(4,2,3)-fvd(4,1,3))
		!write(*,*) 'right',right(id(1,isn)),right(id(2,isn)),right(id(3,isn))
		!write(*,*) 'Acc',axs,ays,azs
		!write(*,*) 'U',d(id(1,isn)),d(id(2,isn)),d(id(3,isn))
		!write(*,*) 'neq',id(1,isn),id(2,isn),id(3,isn)
	!endif
enddo	!ending i
	if (((status1==0.and.itag==1).and.(mod(it,20)==1)).or.((status1==1.and.itag==1).and.(mod(it,10)==0))) then 
	        !open(9002,file='slip.txt',form='formatted',status='unknown',position='append')
		     !   write(9002,'(1x,6e32.21e4)') (fltslp(1,i),fltslr(1,i),fric(26,i),fric(20,i),fric(28,i),fric(30,i),i=1,nftnd)
		open(9002,file='p1output.txt',form='formatted',status='unknown',position='append')
			write(9002,'(1x,5e32.21e4)') (time, fltslp(1,p1nodearr(i)), fltslp(2,p1nodearr(i)), fric(28,p1nodearr(i)),fric(29,p1nodearr(i)),i=1,p1nnode)	
		open(9003,file='p2output.txt',form='formatted',status='unknown',position='append')
			write(9003,'(1x,5e32.21e4)') (time, fltslp(1,p2nodearr(i)), fltslp(2,p2nodearr(i)), fric(28,p2nodearr(i)),fric(29,p2nodearr(i)),i=1,p2nnode)			
    endif
	if (itag==1) then
                pma = maxval(ma_bar_ku_arr)
                maxsliprate = maxval(sliprate_arr) 
                totmomrate = sum(momrate_arr)
                totruptarea = sum(ruptarea_arr)
                tottaoruptarea = sum(taoruptarea_arr)
                totslipruptarea = sum(slipruptarea_arr)
	        !open(9003,file='timehis.txt',form='formatted',status='unknown',position='append')
		!        write(9003,'(1x,i10,3e32.21e4)') it,time, pma, sliprmax	
	endif 		
	dtev=minval(dtev1D)
!-------------------------------------------------------------------!	
end SUBROUTINE faulting	 
