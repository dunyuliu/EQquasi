subroutine faulting

use globalvar
implicit none

!  
integer (kind=4)::ift, i,i1,j,k,n,isn,imn,loc(2)
real (kind = dp) ::slipn,slips,slipd,slip,slipraten,sliprates,sliprated,&
sliprate,xmu,mmast,mslav,mtotl,fnfault,fsfault,fdfault,tnrm,tstk, &
tdip,taox,taoy,taoz,ttao,taoc,ftix,ftiy,ftiz,trupt,tr,&
tmp1,tmp2,tmp3,tmp4,tnrm0,xmu1,xmu2

real (kind = dp) :: dtev1D(maxval(nftnd),ntotft)
real (kind = dp),dimension(6,2,4)::fvd=0.0d0
real (kind = dp)::fa, fb, phi
real (kind = dp) :: rr,R0,T,dtao0, dtao, mr!RSF
real (kind = dp) :: statetmp, eta!RSF
integer (kind=4) :: iv
real (kind = dp) :: tstk0, tdip0,ttao0 !RSF
real (kind = dp) :: dxmudv, rsfeq, drsfeqdv, vtmp, theta_pc_dot, theta_pc_tmp
real (kind = dp) :: v_trial,v_trial_new,sliprs_trial,sliprd_trial, tau_fric_trial,&
    v_s_new_mast,v_s_new_slav,v_d_new_mast,v_d_new_slav,&
    uxm,uym,uzm,uxs,uys,uzs,&
    vxm,vym,vzm,vxs,vys,vzs,&
    axm,aym,azm,axs,ays,azs,&
    anm,asm,adm,ans,ass,ads,&
    slipratemast,sliprateslav,&
    slipaccn,slipaccs,slipaccd
! Sized (maxval(nftnd), ntotft) and indexed (i,ift), not (nftnd(1)) indexed by
! local node number alone: with the latter, fault 2's local node i aliased
! fault 1's slot i, and since these feed pma and maxSlipRate -- which set the
! adaptive step dtev -- the corruption reached the time integrator.
real (kind = dp) :: ma_bar_ku_arr(maxval(nftnd),ntotft), sliprate_arr(maxval(nftnd),ntotft),    momrate_arr(maxval(nftnd),ntotft),      &
                    momRateVW(maxval(nftnd),ntotft),     ruptarea_arr(maxval(nftnd),ntotft),    taoruptarea_arr(maxval(nftnd),ntotft),  &
                    slipruptarea_arr(maxval(nftnd),ntotft)

! All seven accumulators must be zeroed, not just three. momrate_arr,
! momRateVW and ma_bar_ku_arr are assigned only inside the frictional branch
! below, so every fault node in the creeping / no-slip region would otherwise
! contribute uninitialised stack memory to sum(momrate_arr) and
! maxval(ma_bar_ku_arr) -- i.e. to global.dat and the pma printout. That is not
! caught by the regression oracle, which only compares fault.*.nc. The same is
! true of any node beyond nftnd(ift) up to maxval(nftnd), now that every fault
! shares the maxval(nftnd) capacity.
ma_bar_ku_arr = 0.0d0
sliprate_arr = 0.0d0
momrate_arr = 0.0d0
momRateVW = 0.0d0
ruptarea_arr = 0.0d0
taoruptarea_arr = 0.0d0
slipruptarea_arr = 0.0d0
! dtev1D likewise must not carry stack garbage into minval() below for any
! node beyond nftnd(ift); huge() keeps padding out of the minimum.
dtev1D = huge(1.0d0)
do ift = 1, ntotft
    do i=1,nftnd(ift)    !just fault nodes
        fnfault = fric(FR_TNRM0,i,ift) !initial forces on the fault node
        fsfault = fric(FR_TSTK0,i,ift) !norm, strike, dip components directly
        if (icstart ==1) then
          fdfault = 0.0d0
        else
          fdfault = fric(FR_TDIP0,i,ift)
        endif

        isn = nsmp(1,i,ift) 
        imn = nsmp(2,i,ift)     

        do j=1,2  !1-slave, 2-master
            do k=1,3  !1-x comp, 2-y comp, 3-z comp
                fvd(k,j,1) = consf(k,nsmp(j,i,ift))  !1-force !DL 
                fvd(k,j,2) = consv(k,nsmp(j,i,ift)) !2-vel
                if (itag == 1) then
                    fvd(k,j,2) = consvtmp(k,nsmp(j,i,ift)) !2-vel
                endif
                fvd(k,j,3) = cons(k,nsmp(j,i,ift)) !3-disp
                fvd(k,j,4) = consa(k,nsmp(j,i,ift)) !4-acc
            enddo
        enddo

        do j=1,4    !1-force,2-vel,3-disp,4-acc
            do k=1,2  !1-slave,2-master
                fvd(4,k,j) = fvd(1,k,j)*un(1,i,ift) + fvd(2,k,j)*un(2,i,ift) + fvd(3,k,j)*un(3,i,ift)  !4-norm
                fvd(5,k,j) = fvd(1,k,j)*us(1,i,ift) + fvd(2,k,j)*us(2,i,ift) + fvd(3,k,j)*us(3,i,ift)  !5-strike
                fvd(6,k,j) = fvd(1,k,j)*ud(1,i,ift) + fvd(2,k,j)*ud(2,i,ift) + fvd(3,k,j)*ud(3,i,ift)  !6-dip
            enddo
        enddo    

        slipn = fvd(4,2,3) - fvd(4,1,3)
        slips = fvd(5,2,3) - fvd(5,1,3)
        slipd = fvd(6,2,3) - fvd(6,1,3)
        slip = sqrt(slipn**2 + slips**2 + slipd**2) !slip mag
        fric(FR_SLIP_S,i,ift) = slips  !save for final slip output
        fric(FR_SLIP_D,i,ift) = slipd
        fric(FR_SLIP_N,i,ift) = slipn  !normal should be zero, but still keep to ensure
        
        slipraten = fvd(4,2,2) - fvd(4,1,2)
        sliprates = fvd(5,2,2) - fvd(5,1,2)
        sliprated = fvd(6,2,2) - fvd(6,1,2)
        fric(FR_VS_FINAL,i,ift) = sliprates  !save for final slip output
        fric(FR_VD_FINAL,i,ift) = sliprated
        sliprate = sqrt(slipraten**2+sliprates**2+sliprated**2)
            sliprate_arr(i,ift) = sliprate
        if (sliprate>fric(FR_V_PEAK,i,ift)) then 
            fric(FR_V_PEAK,i,ift)=sliprate
        endif
        slipaccn=fvd(4,2,4)-fvd(4,1,4)
        slipaccs=fvd(5,2,4)-fvd(5,1,4)
        slipaccd=fvd(6,2,4)-fvd(6,1,4)
        !...path-itegrated slip for slip-weakening. B.D. 8/12/10
        slp4fri(i,ift) = slp4fri(i,ift) + sliprate * dtev1

        mslav = consm(1,isn)!mass(id(1,isn))        
        mmast = consm(1,imn)!mass(id(1,imn))
        mtotl = mslav + mmast
        if (mmast==0.0d0.or.mslav==0.0d0) then 
            write(*,*) 'ZERO MASS IN FAULTING, MMAST & MSLAV',mmast,mslav
            write(*,*) 'PROBLEMATIC COORDS ARE, X,Z =',x(1,isn)/1.0d3,x(3,isn)/1.0d3
            stop 501
        endif 
        !
        mtotl = mtotl * arn(i,ift)

        tnrm = (mslav*mmast*((fvd(4,2,2)-fvd(4,1,2))+(fvd(4,2,3)-fvd(4,1,3))/dt)/dt &
            + mslav*fvd(4,2,1) - mmast*fvd(4,1,1)) / mtotl + fnfault         
        tstk = (mslav*mmast*(fvd(5,2,2)-fvd(5,1,2))/dt + mslav*fvd(5,2,1) &
            - mmast*fvd(5,1,1)) / mtotl + fsfault
        tdip = (mslav*mmast*(fvd(6,2,2)-fvd(6,1,2))/dt + mslav*fvd(6,2,1) &
            - mmast*fvd(6,1,1)) / mtotl + fdfault
    
        ttao = sqrt(tstk*tstk + tdip*tdip) !total shear magnitude        
        !
        !...friction law to determine friction coefficient
        !   slip-weakening only so far. B.D. 1/26/07
        !... based on choices, call corresponding friction laws.
        ! B.D. 10/8/08
        if (friclaw==1.or.friclaw==2) then!Differ 1&2 and 3&4    
        elseif (friclaw==3.or.friclaw==4) then     
        !---3.2: DECLARE THE TIME FOR RUPTURING.    
            if(fnft(i,ift)<0.0d0) then    !fnft should be initialized by >10000
                if(sliprate >= 0.001d0) then    !first time to reach 1mm/s
                    fnft(i,ift) = time    !rupture time for the node
                endif
            endif    
        !---3.3: INITIATE [V_TRIAL],[SLIPRS_TRIAL],,[SLIPRD_TRIAL]
            v_trial=sliprate
            sliprs_trial=sliprates
            sliprd_trial=sliprated
            
        !---3.4: FRICTIONAL/NON FRIVTIONAL REGION    
            if (abs(x(3,isn)) <= abs(zminc) .and. x(1,isn) <= xmaxc .and. x(1,isn) >= xminc) then 
        !---3.4.1: FRICTIONAL REGION CONTROLLED BY RSF. 
    !---3.4.1.2: Dynamic + NON-DYNAMIC PROCESS: [STATUS1]==0        
                tstk0         = (mslav * fvd(5,2,1) - mmast * fvd(5,1,1)) / mtotl + fsfault
                tdip0         = (mslav * fvd(6,2,1) - mmast * fvd(6,1,1)) / mtotl + fdfault    
                tnrm0         = (mslav * fvd(4,2,1) - mmast * fvd(4,1,1)) / mtotl + fnfault

                ! sigma_bar = sigma_bar_0 - p. tnrm0 is negative in compression
                ! and fric(FR_PORE_DP,:,:) holds p >= 0.
                if (bp == 8) tnrm0 = tnrm0 + fric(FR_PORE_DP,i,ift)

                ! Nucleation
                if (bp == 7 .and. icstart == 1) then
                    call nucleation1(x(1,isn), x(3,isn), dtao)
                    tstk0     = tstk0 + dtao
                endif
                
                ! If non-planar fault geometry and elastic material, enforce normal stress caps.
                ! min_norm/max_norm now come from model.txt (par.min_norm,
                ! par.max_norm), defaulting to the -10/-40 MPa that used to be
                ! hardcoded here. They were literals no compset could see, and
                ! they silently overrode the initial condition of every case
                ! that reached this branch -- see the guard in eqquasi.f90.
                if (rough_fault == 1 .and. C_elastic == 1) then
                    if (tnrm0>=min_norm) then 
                        tnrm0 = min_norm
                    elseif (tnrm0<=max_norm) then
                        tnrm0 = max_norm
                    endif
                endif 
                ! Effective normal stress must stay compressive. Applies to every
                ! benchmark, not just bp == 8: whatever drives it -- pore pressure,
                ! normal-stress changes on a non-planar fault, anything added later
                ! -- a non-compressive sigma_bar means the fault is fully unclamped
                ! and the rate-and-state balance tau = f*sigma_bar has no solution.
                ! Clamping it to zero instead, as this used to do, hands the Newton
                ! solve a vanishing normal stress; it returns NaN and then silently
                ! corrupts every subsequent station, profile and global.dat row.
                call check_effective_normal(tnrm0, i, ift, isn)

                ! fric(FR_ABS_KU,i,ift) is abs(KU)
                fric(FR_ABS_KU,i,ift)= sqrt((mslav * fvd(5,2,1) - mmast * fvd(5,1,1))**2 + (mslav * fvd(6,2,1) - mmast * fvd(6,1,1))**2) / (mmast + mslav) 
                ttao0         = sqrt(tstk0 * tstk0 + tdip0 * tdip0)        
                

                !theta = L/V + (theta - L/V)*dexp(-V*dtev1/L)
                !- fric(22): theta*
                !- fric(21): theta_dot
                !- [fric(11): L],[fric(FR_RSF_V0,i): V0]
    !---3.4.1.3: TO COMPUTE THETA*(t+1) [FRIC(22,:)] FOR ITAG==0, USING [V_TRIAL]=[CONSTRAINV] 
    !----------: THETA**(t+1) [FRIC(22,:)] FOR ITAG==1, , USING [V_TRIAL]={[CONSTRAINV]+[CONSTRAINVTMP]}*0.5    
                !if (itag == 0) then 
                    v_trial   = sliprate
                    fric(FR_V_FINAL,i,ift) = sliprate !Only record the final sliprate in fric(FR_V_FINAL,i) in the last time step.
                    ! phi = dlog(fric(FR_RSF_V0,i,ift) * fric(FR_STATE,i,ift) / fric(FR_RSF_DC,i,ift))
                    ! if (v_trial * dtev1 / fric(FR_RSF_DC,i,ift) <= 1.0d-6) then 
                        ! fric(FR_STATE_TRIAL,i,ift) = dlog(dexp(phi)*(1-v_trial*dtev1/fric(FR_RSF_DC,i,ift)) + fric(FR_RSF_V0,i,ift)*dtev1/fric(FR_RSF_DC,i,ift))
                        ! fric(FR_STATE_TRIAL,i,ift) = fric(FR_RSF_DC,i,ift)/fric(FR_RSF_V0,i,ift)*dexp(fric(FR_STATE_TRIAL,i,ift))
                    ! elseif (v_trial * dtev1 / fric(FR_RSF_DC,i,ift) > 1.0d-6) then 
                        ! fric(FR_STATE_TRIAL,i,ift) = dlog(fric(FR_RSF_V0,i,ift)/v_trial + (dexp(phi)-fric(FR_RSF_V0,i,ift)/v_trial)*dexp(-v_trial*dtev1/fric(FR_RSF_DC,i,ift))) 
                        ! fric(FR_STATE_TRIAL,i,ift) = fric(FR_RSF_DC,i,ift)/fric(FR_RSF_V0,i,ift)*dexp(fric(FR_STATE_TRIAL,i,ift))
                    ! endif
                ! elseif (itag == 1) then 
                    ! do j=1,2  !1-slave, 2-master
                        ! do k=1,3  !1-x comp, 2-y comp, 3-z comp 
                            ! fvd(k,j,2) = consvtmp(k,nsmp(j,i,ift)) !2-vel
                        ! enddo
                    ! enddo    
                    ! do k=1,2  !1-slave,2-master
                        ! fvd(4,k,2) = fvd(1,k,2)*un(1,i,ift) + fvd(2,k,2)*un(2,i,ift) + fvd(3,k,2)*un(3,i,ift)  !4-norm
                        ! fvd(5,k,2) = fvd(1,k,2)*us(1,i,ift) + fvd(2,k,2)*us(2,i,ift) + fvd(3,k,2)*us(3,i,ift)  !5-strike
                        ! fvd(6,k,2) = fvd(1,k,2)*ud(1,i,ift) + fvd(2,k,2)*ud(2,i,ift) + fvd(3,k,2)*ud(3,i,ift)  !6-dip
                    ! enddo    
                    ! slipraten = fvd(4,2,2) - fvd(4,1,2)
                    ! sliprates = fvd(5,2,2) - fvd(5,1,2)
                    ! sliprated = fvd(6,2,2) - fvd(6,1,2)
                    ! v_trial = sqrt(slipraten**2+sliprates**2+sliprated**2)        
                    ! phi = dlog(fric(FR_RSF_V0,i,ift) * fric(FR_STATE,i,ift) / fric(FR_RSF_DC,i,ift))
                    ! if (v_trial * dtev1 / fric(FR_RSF_DC,i,ift) <= 1.0d-6) then 
                        ! fric(FR_STATE_TRIAL,i,ift) = dlog(dexp(phi)*(1-v_trial*dtev1/fric(FR_RSF_DC,i,ift)) + fric(FR_RSF_V0,i,ift)*dtev1/fric(FR_RSF_DC,i,ift))
                        ! fric(FR_STATE_TRIAL,i,ift) = fric(FR_RSF_DC,i,ift)/fric(FR_RSF_V0,i,ift)*dexp(fric(FR_STATE_TRIAL,i,ift))
                    ! elseif (v_trial * dtev1 / fric(FR_RSF_DC,i,ift) > 1.0d-6) then 
                        ! fric(FR_STATE_TRIAL,i,ift) = dlog(fric(FR_RSF_V0,i,ift)/v_trial + (dexp(phi)-fric(FR_RSF_V0,i,ift)/v_trial)*dexp(-v_trial*dtev1/fric(FR_RSF_DC,i,ift))) 
                        ! fric(FR_STATE_TRIAL,i,ift) = fric(FR_RSF_DC,i,ift)/fric(FR_RSF_V0,i,ift)*dexp(fric(FR_STATE_TRIAL,i,ift))
                    ! endif        
                ! endif
                
                ! retrieve state_t+1 calculated from the last time step.
                statetmp = fric(FR_STATE,i,ift)
                fric(FR_STATE_TRIAL,i,ift) = statetmp
                ! initialize statetmp for the first Newton-Raphson update.
                if      (friclaw == 3) then
                    call rate_state_ageing_law(v_trial,fric(FR_STATE_TRIAL,i,ift),fric(1,i,ift),xmu,dxmudv) !RSF
                elseif  (friclaw == 4) then
                    call rate_state_slip_law(v_trial,fric(FR_STATE_TRIAL,i,ift),fric(1,i,ift),xmu,dxmudv) !RSF
                endif
                
                ! retrieve theta_pc_t+1 calculated from the last time step for the state var of effective normal stress.
                theta_pc_tmp   = fric(FR_THETA_PC,i,ift)
                fric(FR_THETA_PC_TRIAL,i,ift) = theta_pc_tmp
                call rate_state_normal_stress(v_trial, fric(FR_THETA_PC_TRIAL,i,ift), theta_pc_dot, tnrm0, fric(1,i,ift))
                
                eta            = mat0(1,2)*mat0(1,3)/2.0d0 
                
                ! Newton-Raphson method to solve the slip_rate from RSF given the tractions and other params for every on-fault grids. Tractions are calculated from the finite element volume and projected to the fault surface. 
                ! x_n+1 = x_n - f(x_n)/f'(x_n)
                ! Here, residual f(v) = shear_traction - frictional_strength. eq(8), Jiang et al. (2022).
                ! shear_traction = static_shear_traction - eta*v. 
                ! Use theta_pc, the state variable that accounts for normal stress change for the effective normal stress.
                do iv = 1,ivmax
                    fric(FR_STATE_TRIAL,i,ift) = statetmp 
                    ! update fric(22) and v_trial for the next time step.
                    if(friclaw == 3) then
                        call rate_state_ageing_law(v_trial,fric(FR_STATE_TRIAL,i,ift),fric(1,i,ift),xmu,dxmudv) !RSF
                    elseif(friclaw == 4) then
                        call rate_state_slip_law(v_trial,fric(FR_STATE_TRIAL,i,ift),fric(1,i,ift),xmu,dxmudv) !RSF
                    endif
                    
                    fric(FR_THETA_PC_TRIAL,i,ift) = theta_pc_tmp 
                    ! update fric(24) and v_trial for the next time step.
                    call rate_state_normal_stress(v_trial, fric(FR_THETA_PC_TRIAL,i,ift), theta_pc_dot, tnrm0, fric(1,i,ift))
                    
                    ! trial shear traction = xmu * theta_pc_tmp 
                    tau_fric_trial = xmu * theta_pc_tmp ! new frictional traction.
                    ! [NOTE] what we input here - tstk0, tdip0, tnrm0, ttao0 are already new static tractions we want to match. We try to seek solution of v_trial that allows tau_fric_trial = (ttao0-eta*v_trial) to match ttao0.
                    ! calculate the residual of RSF friction equation. See eq (8) in Jiang et al. (2022).
                    ! default case: residual = trial shear traction - (friction - eta*sliprate), where eta*sliprate is the radiation damping approximation of inertia (Rice, 1993).
                    ! static case: residual = trial shear traction - old shear traction. When slip rate is low, it should be equivalent to the default case.
                    
                    !rsfeq = (tau_fric_trial-ttao0) + v_trial*eta
                    !drsfeqdv = -dxmudv * tnrm0 +  1.0d0*eta
                    
                    ! for mode 1, quasi-dynamic with radiating damping term 'eta'.
                    if (eqquasi_mode == 1) then 
                        rsfeq = tau_fric_trial - (ttao0 - v_trial*eta) ! f(v_n)
                        drsfeqdv = dxmudv * theta_pc_tmp + eta ! f'(v_n)
                    ! for mode 2, quasi-static without radiating damping term 'eta'.
                    elseif (eqquasi_mode == 2) then 
                        rsfeq = tau_fric_trial - ttao0 ! f(v_n)
                        drsfeqdv = dxmudv * theta_pc_tmp! f'(v_n)
                    endif 
                    
                    ! exiting critera for the Newton-Raphson solver.
                    ! residul stress is 1e-10 time of the new static traction.
                    if (abs(rsfeq) < 1.0d-10 * ttao0) exit 
                    
                    vtmp = v_trial -  rsfeq / drsfeqdv ! v_n+1 = v_n - f(v_n)/f'(v_n+1)
                    
                    if(vtmp <= 0.0d0) then ! avoid negative slip_rate
                        v_trial = v_trial/2.0d0
                    else
                        v_trial = vtmp ! renew v_n with v_n+1
                    endif                
                enddo ! end looping for the Newton-Raphson solver.    
                
                if (iv == ivmax) then 
                    write(*,*) "WARNING, Newton-Raphson reaches the maximum iteration ... ..."
                    write(*,*) "Node (x,y) location", x(1,isn), x(3,isn)
                    write(*,*) "Slip rate solution is, ", v_trial 
                    write(*,*) "Error of stress (in MPa) is ", rsfeq/1.0d6
                    write(*,*) 'sliprate, norm(MPa), tstk, tdip, state2(MPa)'
                    write(*,*) v_trial, tnrm0/1e6, tstk0/1e6, tdip0/1e6, theta_pc_tmp/1e6
                endif 
                
                if (v_trial > 20.0d0) then 
                    write(*,*) 'WARNING: node (x,y) yields faulty slip rate', x(1,isn), x(3,isn)
                    write(*,*) 'sliprate, norm(MPa), tstk, tdip, state2(MPa)'
                    write(*,*) v_trial, tnrm0/1e6, tstk0/1e6, tdip0/1e6, theta_pc_tmp/1e6
                endif
                tstk=tau_fric_trial*tstk0/ttao0 ! new shear traction along strike.
                tdip=tau_fric_trial*tdip0/ttao0 ! new shear traction along dip.
                ttao=sqrt(tstk**2+tdip**2) ! new total shear traction.
                fric(FR_V_TRIAL,i,ift)=v_trial 
                fric(FR_TSTK,i,ift)=tstk
                fric(FR_TDIP,i,ift)=tdip
                fric(FR_TNRM,i,ift)=tnrm0
                
                ! distribute slip rate to master-slave nodes according to masses to preserve moments. 
                ! [NOTE]: this is for right-lateral strike-slip fault. How could we make it general?
                slipratemast=(v_trial)*mslav/(mmast+mslav) 
                sliprateslav=-(v_trial)*mmast/(mmast+mslav)
                
                ! convert slip rates to velocities in the FEM coordinate system.
                v_s_new_mast=slipratemast*tstk/ttao
                v_d_new_mast=slipratemast*tdip/ttao
                v_s_new_slav=sliprateslav*tstk/ttao
                v_d_new_slav=sliprateslav*tdip/ttao
                
                vxm=v_s_new_mast*us(1,i,ift)+v_d_new_mast*ud(1,i,ift)
                vym=v_s_new_mast*us(2,i,ift)+v_d_new_mast*ud(2,i,ift)
                vzm=v_s_new_mast*us(3,i,ift)+v_d_new_mast*ud(3,i,ift)
                vxs=v_s_new_slav*us(1,i,ift)+v_d_new_slav*ud(1,i,ift)
                vys=v_s_new_slav*us(2,i,ift)+v_d_new_slav*ud(2,i,ift)
                vzs=v_s_new_slav*us(3,i,ift)+v_d_new_slav*ud(3,i,ift)

                if (itag == 0) then         
                    consvtmp(1,imn)=vxm
                    consvtmp(2,imn)=vym 
                    consvtmp(3,imn)=vzm 
                    consvtmp(1,isn)=vxs
                    consvtmp(2,isn)=vys 
                    consvtmp(3,isn)=vzs     
                elseif (itag == 1) then 
                    fric(FR_STATE,i,ift)      = fric(FR_STATE_TRIAL,i,ift)  ! update state at itag == 1
                    fric(FR_THETA_PC,i,ift)      = fric(FR_THETA_PC_TRIAL,i,ift)  ! update state2 at itag == 1
                    ma_bar_ku_arr(i,ift)    = (v_trial - fric(FR_V_FINAL,i,ift)) / dtev1 * mmast * mslav / (mmast + mslav) / fric(FR_ABS_KU,i,ift)
                    ma_bar_ku_arr(i,ift)    = abs(ma_bar_ku_arr(i,ift))
                    momrate_arr(i,ift)      = mat0(1,2)**2*mat0(1,3)*v_trial*dx*dx
                                        momRateVW(i,ift)        = 0.0d0
                                        if (bp == 7) then 
                                          if (sqrt(x(1,isn)**2+x(3,isn)**2)<200.0d0) then 
                                            momRateVW(i,ift)    = mat0(1,2)**2*mat0(1,3)*v_trial*dx*dx
                                          endif
                                        endif
                    ruptarea_arr(i,ift)     = 0.0d0
                    taoruptarea_arr(i,ift)  = 0.0d0
                    slipruptarea_arr(i,ift) = 0.0d0
                    
                    if (v_trial>=slipr_thres) then
                        ruptarea_arr(i,ift)     = dx*dx
                        taoruptarea_arr(i,ift)  = ttao*dx*dx
                        slipruptarea_arr(i,ift) = slip*dx*dx
                    endif
                                    
                    consv(1,imn)=vxm
                    consv(2,imn)=vym 
                    consv(3,imn)=vzm 
                    consv(1,isn)=vxs
                    consv(2,isn)=vys 
                    consv(3,isn)=vzs     
                    
                    fric(FR_VXM,i,ift) = vxm 
                    fric(FR_VYM,i,ift) = vym 
                    fric(FR_VZM,i,ift) = vzm 
                    fric(FR_VXS,i,ift) = vxs
                    fric(FR_VYS,i,ift) = vys 
                    fric(FR_VZS,i,ift) = vzs 
                endif 
        !---3.4.1.6: WHEN ITAG==0, SIMPLY STORE V* INTO [CONSTRAINVTMP]
        !----------: WHEN ITAG==1, DECLARE [FRIC(22,:)] AND FINAL V** INTO [CONSTRAINV]            

            else ! The if for fault regions.
        ! !---3.4.2: LOADING BOTTOM & SIDES AT A FIXED SLIDING RATE    
                !tstk0 = 2.585534683723515d7/2.0d0
                tdip0 = 0.0d6
                tnrm = init_norm ! -25.0d6 for bp5. No change for creeping region.
                tnrm0 = tnrm + fric(FR_PORE_DP,i,ift)
                call check_effective_normal(tnrm0, i, ift, isn)
                ! shear stress tstk0 at steady state.
                ! NOTE: outside the frictional region this is the steady-state
                ! stress for load_slip_rate, NOT the traction the elastic solve
                ! produced. For bp8 load_slip_rate is V_zero = 1e-20, so it
                ! evaluates to a constant far below tau^0 and shows up in
                ! diagnostic output as an apparent uniform stress drop at
                ! |x| > l_f. It is a placeholder for a locked region, not a
                ! boundary artefact -- do not read it as one. Section 4.3
                ! profiles span exactly +-400 m, and the branch above takes
                ! x = +-l_f inclusive, so submitted output never contains it.
                call rsf_rd(tstk0, tnrm0, fric(FR_RSF_A,i,ift), fric(FR_RSF_B,i,ift), fric(FR_RSF_F0,i,ift), fric(FR_RSF_V0,i,ift), mat0(1,2), mat0(1,3), load_slip_rate)
                
                v_trial = load_slip_rate
                fric(FR_V_TRIAL,i,ift) = v_trial
                fric(FR_TSTK,i,ift) = tstk0
                fric(FR_TDIP,i,ift) = tdip0
                fric(FR_TNRM,i,ift) = tnrm0
                slipratemast=(v_trial)*mslav/(mmast+mslav)
                sliprateslav=-(v_trial)*mslav/(mmast+mslav)
                v_s_new_mast=slipratemast
                v_d_new_mast=0.0d0
                v_s_new_slav=sliprateslav
                v_d_new_slav=0.0d0
                vxm=v_s_new_mast*us(1,i,ift)+v_d_new_mast*ud(1,i,ift)
                vym=v_s_new_mast*us(2,i,ift)+v_d_new_mast*ud(2,i,ift)
                vzm=v_s_new_mast*us(3,i,ift)+v_d_new_mast*ud(3,i,ift)
                vxs=v_s_new_slav*us(1,i,ift)+v_d_new_slav*ud(1,i,ift)
                vys=v_s_new_slav*us(2,i,ift)+v_d_new_slav*ud(2,i,ift)
                vzs=v_s_new_slav*us(3,i,ift)+v_d_new_slav*ud(3,i,ift)
                 
                 consvtmp(1,imn)=vxm
                 consvtmp(2,imn)=vym 
                 consvtmp(3,imn)=vzm 
                 consvtmp(1,isn)=vxs
                 consvtmp(2,isn)=vys 
                 consvtmp(3,isn)=vzs 
                
                 consv(1,imn)=vxm
                 consv(2,imn)=vym 
                 consv(3,imn)=vzm 
                 consv(1,isn)=vxs
                 consv(2,isn)=vys 
                 consv(3,isn)=vzs     
                fric(FR_VXM,i,ift) = vxm 
                fric(FR_VYM,i,ift) = vym 
                fric(FR_VZM,i,ift) = vzm 
                fric(FR_VXS,i,ift) = vxs
                fric(FR_VYS,i,ift) = vys 
                fric(FR_VZS,i,ift) = vzs
            endif
            
            ! update new time step size according to new slip rate v_trial.
            dtev1D(i,ift) = ksi * fric(FR_RSF_DC,i,ift)/v_trial
            
            ! [WARNING]: exit the code if time step size is negative.
            if (dtev1D(i,ift) < 0) then
                write(*,*) 'NEGATIVE SLIPRATE, ITS LOC = ', x(1,isn),x(3,isn)
                write(*,*) 'PROBLEMATIC V_TRIAL = ', v_trial
                stop 502
            endif
            ! Record rupture area, average stress vector, slip
            if (itag==1) then
                fric(FR_RUPT_SLIP,i,ift) = sqrt(slips**2 + slipn**2) ! Total slip
                if (fric(FR_V_TRIAL,i,ift) > 1.0d-3) then 
                    fric(FR_RUPT_AREA,i,ift) = arn(i,ift) ! Record ruptured area
                endif 
                if (abs(tdynastart-time)<dt/100.0) then 
                    fric(FR_TRACT_START,i,ift) = sqrt(fric(FR_TSTK,i,ift)**2 + fric(FR_TDIP,i,ift)**2) !Total shear traction when dyna starts.
                endif
                if (abs(tdynaend-time)<dt/100.0) then 
                    fric(FR_TRACT_END,i,ift) = sqrt(fric(FR_TSTK,i,ift)**2 + fric(FR_TDIP,i,ift)**2) !Total shear traction when dyna ends.    
                endif                    
            endif 
        endif

        if (itag==1) then
            if(n4onf>0) then
                do j=1,n4onf
                    if(anonfs(1,j)==i.and.anonfs(3,j)==ift) then !only selected stations. B.D. 10/25/09    
                        fltsta(1,it,j)  = time
                        fltsta(2,it,j)  = v_trial
                        fltsta(3,it,j)  = sliprated
                        fltsta(4,it,j)  = fric(FR_STATE,i,ift)
                        fltsta(5,it,j)  = slips
                        fltsta(6,it,j)  = slipd
                        fltsta(7,it,j)  = slipn
                        fltsta(8,it,j)  = fric(FR_TSTK,i,ift)!tstk
                        fltsta(9,it,j)  = fric(FR_TDIP,i,ift)!tdip
                        fltsta(10,it,j) = tnrm0
                        fltsta(11,it,j) = fric(FR_PORE_DP,i,ift)  ! pore pressure change, Pa
                        fltsta(12,it,j) = fric(FR_DARCY_S,i,ift) ! Darcy velocity along strike, m/s
                        fltsta(13,it,j) = fric(FR_DARCY_D,i,ift) ! Darcy velocity along dip, m/s
                        fltsta(14,it,j) = sliprates      ! slip rate along strike, m/s
                    endif
                enddo 
            endif   
        endif
    enddo    !ending i
enddo ! ending ift

    if (me == 0) call output_prof

    if (itag==1) then
        pma             = maxval(ma_bar_ku_arr)
        maxSlipRate     = maxval(sliprate_arr) 
        loc             = maxloc(sliprate_arr)
        !write(*,*) 'maxSlipRate at ', loc(1)
        totMomRate      = sum(momrate_arr)
                totMomRateVW    = sum(momRateVW)
        totRuptArea     = sum(ruptarea_arr)
        totTaoRuptArea  = sum(taoruptarea_arr)
        totSlipRuptArea = sum(slipruptarea_arr)
    endif         
    
    dtev=minval(dtev1D)
!-------------------------------------------------------------------!    
end subroutine faulting     

! Subroutine rate_state_normal_stress calculates the effect of normal stress change
! from RSF. The formulation follows Shi and Day (2013), eq B8. {"Frictional sliding experiments with variable normal stress show that the shear strength responds gradually to abrupt changes of normal stress (e.g., Prakash and Clifton, 1993; Prakash, 1998)."}

! theta_pc_dot = - V/L_pc*[theta_pc - abs(tnrm)]

! Input: slip_rate, L_pc, theta_pc, tnrm. 
! Output: theta_pc, the state variable which is used to calculate shear stress in eq B2
! B2: abs(traction) = friction * theta_pc.
subroutine rate_state_normal_stress(V2, theta_pc, theta_pc_dot, tnrm, fricsgl)
    use globalvar
    implicit none
    real (kind = dp) :: V2, theta_pc, theta_pc_dot, tnrm, L_pc
    real (kind = dp),dimension(100) :: fricsgl
    
    L_pc = fric_pc_L

    if (L_pc <= 0.0d0) then
        ! Default. No normal-stress state variable: strength is built from the
        ! instantaneous effective normal stress. Required by SEAS BP8 eq. (10),
        ! where pore pressure changes sigma_bar with essentially no slip -- the
        ! relaxation below is driven by slip, so it would otherwise stay pinned
        ! at its initial value and the fault would never feel the injection.
        theta_pc = abs(tnrm)
        return
    endif

    ! Exponential form: exact for constant V2 and tnrm over the step, and
    ! bounded for any dtev1. The linearised update it replaces overshoots once
    ! V2*dtev1/L_pc > 1, which became reachable once dtev1 could be capped and
    ! so no longer satisfied V2*dtev1/L_pc == xi identically.
    theta_pc = abs(tnrm) + (theta_pc - abs(tnrm))*dexp(-V2*dtev1/L_pc)

end subroutine rate_state_normal_stress

subroutine newton_solver(a_tmp, b_tmp)
    use globalvar
    implicit none
    real (kind = dp):: a_tmp, b_tmp

end subroutine newton_solver

subroutine nucleation1(xtmp, ztmp, dtao)
! Subroutine nucleation1 takes in on-fault node coordinates and return the stress perturbation
!   for the nucleation of first rupture in SCEC SEAS benchmark bp7-qd.
    use              globalvar
    implicit         none
    real (kind = dp) :: xtmp, ztmp, dtao, G1, G2, rtmp
    
    dtao       = 0.0d0 ! initialize 
    rtmp       = sqrt((xtmp - nucx)**2 + (ztmp - nucz)**2)
    
    if (rtmp < nucr) then 
        G1     = exp(rtmp**2/(rtmp**2-nucr**2))
    else
        G1     = 0.0d0
    endif 
    
    if (time < nuct) then 
        G2     = exp((time-nuct)**2/time/(time-2.0d0*nuct))
    else
        G2     = 1.0d0
    endif 
    
    dtao       = nucdtao0*G1*G2

end subroutine

! check_effective_normal stops the run if the effective normal stress is no
! longer compressive.
!
! Zero effective normal stress means the pore pressure has met the ambient
! effective normal stress and the fault is fully unclamped. The quasi-dynamic
! rate-and-state balance tau = f(V,theta)*sigma_bar has no solution there, so
! there is nothing meaningful to continue with. This is reachable with the
! literal SEAS BP8 parameters under the Peaceman well source: eq. (25) of the
! benchmark description gives 25.3 MPa at the equivalent radius by day 2 at
! dx = 50 m, and 44 MPa by t_off at dx = 10 m, against sigma_bar_0 = 25 MPa.
! A point source in 2D is logarithmically singular, so this gets worse, not
! better, as the mesh is refined.
subroutine check_effective_normal(tnrm0, i, ift, isn)

    use globalvar
    implicit none

    real (kind = dp), intent(in) :: tnrm0
    integer (kind = 4), intent(in) :: i, ift, isn

    if (tnrm0 < 0.0d0) return

    write(*,*) '===================================================================='
    write(*,*) '= EFFECTIVE NORMAL STRESS IS NO LONGER COMPRESSIVE                 ='
    write(*,*) '= The fault is fully unclamped and rate-and-state has no solution. ='
    write(*,'(X,A,E15.7,A)')  '=   time              = ', time, ' s'
    write(*,'(X,A,2E15.7,A)') '=   node x, z         = ', x(1,isn), x(3,isn), ' m'
    write(*,'(X,A,E15.7,A)')  '=   sigma_bar         = ', tnrm0/1.0d6, ' MPa'
    write(*,'(X,A,E15.7,A)')  '=   pore pressure p   = ', fric(FR_PORE_DP,i,ift)/1.0d6, ' MPa'
    write(*,'(X,A,E15.7,A)')  '=   sigma_bar_0       = ', fric(FR_TNRM0,i,ift)/1.0d6, ' MPa'
    write(*,*) '= Clamping this to zero would return NaN from the Newton solve and  ='
    write(*,*) '= silently corrupt every later output row, so the run stops here.   ='
    write(*,*) '===================================================================='
    stop 508

end subroutine check_effective_normal
