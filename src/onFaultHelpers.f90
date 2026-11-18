! Subroutines shared by the three time loops (solveTimeLoopMUMPS, solveTimeLoopPETSc,
! solveTimeLoopPETScDD). Moved unchanged from solveTimeLoopMUMPS.f90 so that solvers 2
! and 3 do not depend on the solver-1 loop file for them.
module onFaultHelpers
contains

subroutine exitCriteria
! Subroutine to compute if to exit EQquasi. 
! If mode == 1 (quasi-dynamic/quasi-static), the exit happens AFTER 
!    a rupture when maxSlipRate falls below slipr_thres.
! If mode == 2 (fully dynamic), the exit happens BEFORE a rupture when 
!    the maxSlipRate rises above the slipr_thres.

! The system is determined by the combination of two status variables - status0/1 and t_end_status.
! status0/status1: code status in the last/current time step.
! 0: interseismic; 1: rupture. 
! 
    use globalvar
    implicit none
    real (kind = dp) :: hang_time = 1e5

    ! bp8 is aseismic by construction (velocity-strengthening everywhere), so it
    ! runs to a prescribed final time instead of exiting on a slip rate criterion.
    if (bp == 8) then
        if (fluid_tend > 0.0d0 .and. time >= fluid_tend) stoptag = 1
        return
    endif

    if (eqquasi_mode == 1) then ! quasi-dynamic/quasi-static
        if (maxSlipRate>slipr_thres)then
        ! if maximum slip rate rises above the slipr_thres, entering/in the co-seismic phase.
            if (status0 == 0) then ! convert from inter-seismic to co-seismic.
                status1 = 1 ! current status is changed to co-seismic.
            elseif (status0 == 1) then ! convert from co-seismic to inter-seismic.
                if (t_end_status>0.0d0) then ! if t_end_status is TRUE, set it to FALSE and continue rupturing.
                    t_end_status = -1.0d0 ! t_end_status is FALSE.
                    tdynaend = -1000.0d0 ! don't need to record end time in the rupture phase.
                endif
            endif
        elseif (maxSlipRate<=slipr_thres) then
        ! if maximum slip rate falls below the slipr_thres, consider exiting the code. 
            if (status0 == 0) then 
                status1 = 0 ! if last step was in inter-seismic, no need to change.
            elseif (status0 == 1) then !if last step was in co-seismic, consider to exit.
                if (t_end_status<0.0d0) then 
                    t_end_status = 1.0d0 ! change exit status to be TRUE.
                    tdynaend = time ! record the end time of the rupture.
                else
                    if ((time-tdynaend)>=hang_time .and. tdynaend>0.0d0) then ! exit if the low slip rate status is kept for over 100 seconds.
                    ! NOTE: the 100 seconds threshold is subjective to change. 
                        stoptag = 1 ! stoptag is TRUE now. Direct the code to exit.
                    endif
                    !if ((time-tdynaend)>=1.0d0*7*24*60*60.and.tdynaend>0.0d0) exit
                endif
            endif
        endif
        
        if ((status1-status0)==1) then ! record time if changing from inter-seismic to co-seismic.
            tdynastart = time
            t_start_status = 1
        endif
    elseif (eqquasi_mode == 2) then ! fully dynamic
        if (maxSlipRate>slipr_thres)then
        ! In this mode, if maximum slip rate rises above the slipr_thres, exit. 
            if (status0 == 0) then ! convert from inter- to co-seismic.
                status1 = 1
                stoptag = 1 ! stoptag is TRUE now. Direct the code to exit.
            endif
        endif
    endif 
end subroutine exitCriteria

subroutine getScalarOnFaultQuant
    use globalvar
    implicit none
    integer (kind = 4) :: ift
    maxSlipRate = 0.0d0
    do ift = 1, ntotft
        if (nftnd(ift) > 0) maxSlipRate = max(maxSlipRate, maxval(fric(FR_V_CURRENT, 1:nftnd(ift), ift)))
    enddo
end subroutine getScalarOnFaultQuant

subroutine bound_load

    use globalvar
    implicit none
    
    integer (kind = 4) :: nel, node_num, i, j, ntag, inv, jnv, temp, elemvar(nee)
    real (kind = dp) :: elemmat(5), elemx(ndof,nen), es(nee,nee), em(nee), ef(nee), elemu(nee), Vol
    
    right=0.0d0
    constmp=0.0d0
!write(*,*) '!---3.5.1===UPDATE BOUNDARY U*(t+1) INTO [CONSTRAINTMP]'    
!--------===STORE EQUVILENT FORCES INTO [RIGHT]
    do nel=1,numel
        if (et(nel)==2.or.et(nel)==3) then 
            do j=1,nen
                node_num=ien(j,nel)
                do i=1,ndof
                    elemx(i,j)=x(i,node_num)
                enddo
            enddo
            do j=1,5
                elemmat(j)=mat(nel,j)
            enddo

            call c8g2(elemx,elemmat,es,em,ef,Vol)
            ntag=0
            do i=1,nen
                node_num=ien(i,nel)
                do j=1,ndof
                ntag=ntag+1
                elemvar(ntag)=id(j,node_num)
                !-2->-x; -3->+x; -5->fault
                if (elemvar(ntag)<0)then
                    if (itag == 0) then 
                        if (elemvar(ntag)==-2) then ! ymin, x- vel
                            consv(j,node_num) = -far_load_rate 
                        elseif (elemvar(ntag)==-3) then ! ymax, x- vel
                            consv(j,node_num) = far_load_rate 
                        elseif (elemvar(ntag)==-21) then ! ymin, y- vel
                            consv(j,node_num) = -far_norm_load_vel      
                        elseif (elemvar(ntag)==-31) then ! ymax, y+ vel
                            consv(j,node_num) = far_norm_load_vel         
                        elseif (elemvar(ntag)==-1) then
                            consv(j,node_num) = 0.0d0
                        endif 
                        consvtmp(j,node_num) = consv(j,node_num)
                        consa(j,node_num) = 0.0d0    
                        !First predictions of U*(t+1)[CONSTRAINTMP]
                        constmp(j,node_num)=cons(j,node_num)+consv(j,node_num)*dtev1
                    elseif (itag == 1) then 
                        !-2->-x; -3->+x; -5->fault
                        if (elemvar(ntag)==-2) then
                            consv(j,node_num) = -far_load_rate
                        elseif (elemvar(ntag)==-3) then
                            consv(j,node_num) = far_load_rate   
                        elseif (elemvar(ntag)==-21) then
                            consv(j,node_num) = -far_norm_load_vel
                        elseif (elemvar(ntag)==-31) then
                            consv(j,node_num) = far_norm_load_vel         
                        elseif (elemvar(ntag)==-1) then
                            consv(j,node_num) = 0.0d0
                        endif 
                        consa(j,node_num) = 0.0d0
!write(*,*) '!---3.5.5.1-FINALL UPDATION OF [CONSTRAIN].'
!-----------NO NEED TO WORRY NON-BOUNDARY NODES, SINCE THEIR CONSTRAINV&VTMP ARE ZERO.
!-----------ERROR, CANNOT UPDATE CONTRAIN ITSELF. EACH D.O.F HAS BEEN VISITIED SEVERAL TIMES FOR IN THE LOOP.
                        constmp(j,node_num) = cons(j,node_num) + 0.5d0*(consv(j,node_num)+consvtmp(j,node_num))*dtev1                    
                    endif 
                    elemu(ntag) = constmp(j,node_num)
                endif
                
                enddo
            enddo            
            do j=1,nee
                jnv=elemvar(j) 
                if (jnv.lt.0) then
                    do i=1,nee                   
                        inv=elemvar(i)
                        if (inv.gt.0) then
                            !ONLY ELEMU(J) WITH NEGATIVE J WILL CONTRIBUTE 
                            right(inv) = right(inv) - es(i,j)*elemu(j)
                        endif
                    enddo
                endif
            enddo                
        endif
        
    enddo

end subroutine bound_load


subroutine bound_ft_ku

    use globalvar
    implicit none
    
    integer (kind = 4) :: nel, node_num, i, j, m, ntag, inv, jnv, temp, elemvar(nee)
    real (kind = dp) :: elemmat(5), elemx(ndof,nen), es(nee,nee), em(nee), ef(nee), elemu(nee), Vol 
            
    consf=0.0d0
    if (itag == 0) then 
        consm=0.0d0!Updated only when itag == 0
    endif 
    
    do nel=1,numel!Compute traction. 
        if (et(nel)==2) then 
            do j=1,nen
                node_num=ien(j,nel)
                do i=1,ndof
                    elemx(i,j)=x(i,node_num)
                enddo
            enddo
            do j=1,5
                elemmat(j)=mat(nel,j)
            enddo
            ntag=0
            call c8g2(elemx,elemmat,es,em,ef,Vol)        
            do i=1,nen!F*=KU*(t+1)
                node_num=ien(i,nel)
                do j=1,ndof 
                    ntag=ntag+1
                    if (itag == 0) then 
                        elemu(ntag)=constmp(j,node_num)
                    elseif (itag == 1) then 
                        elemu(ntag)=cons(j,node_num)
                    endif 
                enddo 
            enddo        
            do i=1,nen!-KU(t+1)
                node_num=ien(i,nel)
                do j=1,ndof 
                    if (id(j,node_num)<0) then 
                        do m=1,nee
                            consf(j,node_num)=consf(j,node_num)-es((i-1)*3+j,m)*elemu(m)
                        enddo
                        if (itag == 0) then 
                            consm(j,node_num)=consm(j,node_num)+em((i-1)*3+j)
                        endif 
                    endif
                enddo 
            enddo                    
        endif
    enddo!nel
end subroutine bound_ft_ku

subroutine initOnFaultKinematics
    use globalvar
    implicit none
    integer (kind = 4) :: i, ift
    dtev       = ksi*minDc/maxSlipRate
    dtev1      = ksi*minDc/maxSlipRate!Initial dtev1>dt to enter static state.
    cons       = 0.0d0
    constmp    = 0.0d0
    consvtmp   = 0.0d0
    consv      = 0.0d0
    consa      = 0.0d0
    do ift = 1, ntotft
        if (nftnd(ift) > 0) then !RSF
            do i=1,nftnd(ift)
                if (ddActive) then
                    if (.not. fltMine(i,ift)) cycle ! par.solver = 3: the owner sets it and sends it on.
                endif
                if (icstart == 1) then
                    consv(1,nsmp(1,i,ift)) =  -fric(FR_VINIT,i,ift)/2.0d0
                    consv(1,nsmp(2,i,ift)) =  fric(FR_VINIT,i,ift)/2.0d0
                else
                consv(1:3,nsmp(1,i,ift)) = fric(34:36,i,ift)
                consv(1:3,nsmp(2,i,ift)) = fric(31:33,i,ift)
                endif
            enddo
        endif
    enddo
end subroutine initOnFaultKinematics
end module onFaultHelpers
