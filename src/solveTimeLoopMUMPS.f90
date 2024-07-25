subroutine solveTimeLoopMUMPS

    use globalvar
    implicit none
    include 'mpif.h'
    include 'dmumps_struc.h'
    
    TYPE (DMUMPS_STRUC) mumps_par
    integer (kind = 4) :: i,j,inv,jnv,ntag,node_num,var,l,k, iiTag
    real (kind = dp) :: startTime, endTime, timeUsedInFactorization,&
        timeUsedInComputing
    character (len = 50) :: netcdf_outfile, output_type
    
    mumps_par%COMM = MPI_COMM_WORLD
    mumps_par%JOB = -1
    mumps_par%SYM = 1! 0: unsymmetric; 1: Positive definite symmetric ; 2: general symmetric
    mumps_par%PAR = 1
    call DMUMPS(mumps_par)
 
    call getScalarOnFaultQuant
    
    if (me.eq.0) then
        write(*,*) '= Building Stiffness Matrix in CRS format      ='
        call createMatrixHolderInCRSFormat
        call elemAssembleInCRS
        
        write(*,*) '= Converting from CRS to MUMPS format          ='
        mumps_par%N   = neq 
        mumps_par%NNZ = maxa
        allocate( mumps_par%IRN (mumps_par%NNZ))
        allocate( mumps_par%JCN (mumps_par%NNZ))
        allocate( mumps_par%A   (mumps_par%NNZ))
        allocate( mumps_par%RHS (mumps_par%N))
        do i = 1, neq
            do j = ia(i),ia(i+1)-1
                mumps_par%IRN (j) = i
            enddo
        enddo
        mumps_par%JCN = ja
        mumps_par%A   = kstiff
        
        call initOnFaultKinematics
    endif

    call cpu_time(startTime)
    mumps_par%JOB = 4! Combines the actions of JOB=1, the analysis phase, 
    ! and JOB=2, the factorization phase.
    call DMUMPS(mumps_par)
    call cpu_time(endTime)
    timeUsedInFactorization = endTime - startTime
    
    stoptag = 0 ! set stoptag to FALSE.
    
    call cpu_time(startTime)
    do it = 1, nstep 
        if (stoptag == 1) exit ! exit EQquasi if stoptag is TRUE.
        
        if (me == 0) then
            resu_1 = resu
            
            if (bp == 7 .and. icstart == 1) then ! for the first cycle of bp7, if time<nuct, use dt for dtev1.
                if (it < (nuct/dt)) then 
                    dtev1 = dt 
                endif
            endif 
            
            time = time + dtev1
            
            if (mod(it,nhplt) == 1 .and. me ==0) then
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

            if(ndout>0) then
                dout(1,it)=time
                do i=1,ndout
                    j=idhist(1,i)
                    if(j<=0) j=1  !avoid zero that cannot be used below
                        k=idhist(2,i)
                        l=idhist(3,i)
                    if(l==1) then
                        dout(i+1,it)=cons(k,j)
                    elseif(l==2) then
                        dout(i+1,it)=consv(k,j)
                    elseif(l==3) then
                        dout(i+1,it)=consa(k,j)
                    endif
                enddo
            endif  
        endif 
        
        do iiTag = 0, 1
            if (me == 0) then 
                ! When iiTag==1, !V*(t+1) OBTAINED. 
                ! UPDATE BOUNDARY U**(t+1),THEN DECLARE [CONSTRAIN] BY U**(t+1). 
            
                itag = iiTag 
                ! COMPUTE U*(t+1) FOR THE WHOLE VOLUME'
                ! COMPUTE U**(t+1) FOR THE WHOLE VOLUME.
                call bound_load
                mumps_par%RHS = right
            endif
            
            call MPI_BCAST(stoptag, 1, MPI_INT, 0, MPI_COMM_WORLD, IERR)
            mumps_par%JOB = 3 ! Solve with JOB=3.
            mumps_par%ICNTL(3) = -1 
            ! ICNTL(3) Default value is 6.
            ! It controls the output stream of global information.
            ! Disabled here with -1. 
            call DMUMPS(mumps_par)
            
            if (me.eq.0) then
                resu = mumps_par%RHS
            
                if (iiTag==0) then 
                    do i=1,numnp
                        do j=1,ndof
                            if (id(j,i)>0) then !Only nodes that have equation number have been updated here.
                                constmp(j,i)=resu(id(j,i))                    
                            endif
                        enddo 
                    enddo   
                elseif (iiTag==1) then 
                    do i=1,numnp
                        do j=1,ndof
                            if (id(j,i)>0) then 
                                cons(j,i) = resu(id(j,i))
                            else
                                cons(j,i) = constmp(j,i)
                            endif
                        enddo 
                    enddo
                endif 
                
                !COMPUTE F*=KU*(t+1)=[KSTIFF].DOT.[CONSTRAINTMP]'  
                !COMPUTE F**=KU**(t+1)=[KSTIFF].DOT.[CONSTRAIN]' 
                call bound_ft_ku
                !GET FIRST PREDICTIONS OF V*(t+1).
                !GET SECOND PREDICTIONS OF V**(t+1), AND DECLARE V(t+1)=V**(t+1).'  
                call faulting
            endif 
        enddo 
        
        if (me == 0) then 
            status0 = status1
            dtev1 = max(dt,int8(dtev/dt)*dt)    

            do i=1,numnp
                do j=1,ndof
                    if (id(j,i)>0) then 
                        consv(j,i) = (resu(id(j,i)) - resu_1(id(j,i)))/dtev1
                    endif
                enddo 
            enddo                        
            globaldat(1,it) = time
            globaldat(2,it) = maxSlipRate
            globaldat(3,it) = totMomRate    
            globaldat(4,it) = totTaoRuptArea
            globaldat(5,it) = totSlipRuptArea
            globaldat(6,it) = totRuptArea
            globaldat(7,it) = totMomRateVW
            
            if (mod(it,nt_output_stress) == 1 .or. (stoptag == 1) .or. (it==nstep)) then 
                write(proc_str,'(I5.5)') it
                netcdf_outfile = 'disp.'//trim(proc_str)//'.nc'
                output_type = 'disp'
                call netcdf_write(netcdf_outfile, output_type)
                
                netcdf_outfile = 'fault.'//trim(proc_str)//'.nc'
                call netcdf_write_on_fault(netcdf_outfile)
            endif
            ! If exiting, write again the restart files.
            if ((stoptag == 1) .or. (it==nstep)) then
                netcdf_outfile = 'disp.r.nc'
                output_type = 'disp'
                call netcdf_write(netcdf_outfile, output_type)
                
                netcdf_outfile = 'fault.r.nc'
                call netcdf_write_on_fault(netcdf_outfile)
            endif
        endif
    enddo
    call cpu_time(endTime)
    timeUsedInComputing = endTime-startTime
    
    if (me == 0) then 
        write(*,*) it, ' steps use ', timeUsedInComputing, ' seconds ...'
        write(*,*) 'Factorization uses ', timeUsedInFactorization, ' seconds ...'
    endif
    
end subroutine solveTimeLoopMUMPS

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
    maxSlipRate = maxval(fric(47, 1:nftnd(1), 1)) 
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
                        if (elemvar(ntag)==-2) then
                            consv(j,node_num)=-far_load_rate!-5.0d-10
                            consvtmp(j,node_num) = consv(j,node_num)
                            consa(j,node_num)=0.0d0
                        elseif (elemvar(ntag)==-3) then
                            consv(j,node_num)=far_load_rate!5.0d-10
                            consvtmp(j,node_num) = consv(j,node_num)
                            consa(j,node_num)=0.0d0                    
                        elseif (elemvar(ntag)==-1) then
                            consv(j,node_num)=0.0d0
                            consvtmp(j,node_num) = consv(j,node_num)
                            consa(j,node_num)=0.0d0
                        endif 

                        !First predictions of U*(t+1)[CONSTRAINTMP]
                        constmp(j,node_num)=cons(j,node_num)+consv(j,node_num)*dtev1
                    elseif (itag == 1) then 
                        !-2->-x; -3->+x; -5->fault
                        if (elemvar(ntag)==-2) then
                            consv(j,node_num)=-far_load_rate!-5.0d-10
                            consa(j,node_num)=0.0d0
                        elseif (elemvar(ntag)==-3) then
                            consv(j,node_num)=far_load_rate!5.0d-10
                            consa(j,node_num)=0.0d0                    
                        elseif (elemvar(ntag)==-1) then
                            consv(j,node_num)=0.0d0
                            consa(j,node_num)=0.0d0
                        endif 
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
    integer (kind = 4) :: i
    dtev       = ksi*minDc/maxSlipRate
    dtev1      = ksi*minDc/maxSlipRate!Initial dtev1>dt to enter static state.
    cons       = 0.0d0
    constmp    = 0.0d0
    consvtmp   = 0.0d0
    consv      = 0.0d0
    consa      = 0.0d0
    if (nftnd(1) > 0) then !RSF
        do i=1,nftnd(1)
            if (icstart == 1) then 
                consv(1,nsmp(1,i,1)) =  -fric(46,i,1)/2.0d0
                consv(1,nsmp(2,i,1)) =  fric(46,i,1)/2.0d0
            else
            consv(1:3,nsmp(1,i,1)) = fric(34:36,i,1)
            consv(1:3,nsmp(2,i,1)) = fric(31:33,i,1)
            endif     
        enddo    
    endif
end subroutine initOnFaultKinematics