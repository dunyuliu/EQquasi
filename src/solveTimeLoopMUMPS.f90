subroutine solveTimeLoopMUMPS

    use globalvar
    use onFaultHelpers
    implicit none
    include 'mpif.h'
    include 'dmumps_struc.h'
    
    TYPE (DMUMPS_STRUC) mumps_par
    integer (kind = 4) :: i,j,inv,jnv,ntag,node_num,var,l,k, iiTag
    real (kind = dp) :: startTime, endTime, timeUsedInFactorization,&
        timeUsedInComputing
    real (kind = dp) :: progressStartTime
    character (len = 50) :: netcdf_outfile, output_type
    
    mumps_par%COMM = MPI_COMM_WORLD
    mumps_par%JOB = -1
    mumps_par%SYM = 1! 0: unsymmetric; 1: Positive definite symmetric ; 2: general symmetric
    mumps_par%PAR = 1
    call DMUMPS(mumps_par)
 
    call getScalarOnFaultQuant

    if (bp == 8) call pore_pressure_init
    if (bp == 8) call bp8_profile_init

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
    ! Owner request (enhanced log reporting): MUMPS's own error code
    ! (negative INFOG(1) on failure) was never checked here -- a failed
    ! factorization (singular matrix, out of memory inside MUMPS, etc.)
    ! would have silently continued into JOB=3 solves against a
    ! factorization that never happened. This is additive: on every run
    ! that has ever passed the parity gate, INFOG(1) is 0 and this branch
    ! never executes.
    if (mumps_par%INFOG(1) < 0) then
        if (me == 0) write(*,*) 'MUMPS factorization (JOB=4) failed, ', &
            'INFOG(1)=', mumps_par%INFOG(1), ' INFOG(2)=', mumps_par%INFOG(2), &
            ' -- stopping rather than solving against a factorization ', &
            'that never happened.'
        flush(6)  ! ensure the diagnostic just written reaches the log before abort tears the process down
        call MPI_ABORT(MPI_COMM_WORLD, 1, IERR)
    endif
    call cpu_time(endTime)
    timeUsedInFactorization = endTime - startTime
    
    stoptag = 0 ! set stoptag to FALSE.

    call cpu_time(startTime)
    progressStartTime = startTime
    do it = 1, nstep
        if (stoptag == 1) exit ! exit EQquasi if stoptag is TRUE.
        
        if (me == 0) then
            resu_1 = resu
            
            if (bp == 7 .and. icstart == 1) then ! for the first cycle of bp7, if time<nuct, use dt for dtev1.
                if (it < (nuct/dt)) then 
                    dtev1 = dt 
                endif
            endif 
            
            if (dtmax > 0.0d0) dtev1 = min(dtev1, dtmax)

            time = time + dtev1

            ! bp8: advance the along-fault pore pressure to t+dtev1 before
            ! faulting consumes it through fric(FR_PORE_DP,:,:).
            if (bp == 8) call pore_pressure_update(dtev1)

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
            ! Owner request: same additive check as the JOB=4 factorization
            ! above, for the per-step solve. INFOG(1) is MUMPS's own global
            ! (replicated-on-every-rank) error code -- checked on every
            ! rank so the MPI_ABORT below is genuinely collective, not
            ! rank-0-only.
            if (mumps_par%INFOG(1) < 0) then
                if (me == 0) write(*,*) 'MUMPS solve (JOB=3) failed at ', &
                    'step', it, 'iiTag', iiTag, ', INFOG(1)=', &
                    mumps_par%INFOG(1), ' INFOG(2)=', mumps_par%INFOG(2), &
                    ' -- stopping rather than using an unsolved RHS.'
                flush(6)  ! ensure the diagnostic just written reaches the log before abort tears the process down
                call MPI_ABORT(MPI_COMM_WORLD, 1, IERR)
            endif

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
            if (bp == 8) call bp8_profile_record

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
                netcdf_outfile = trim(outDir)//'disp.'//trim(proc_str)//'.nc'
                output_type = 'disp'
                call netcdf_write(netcdf_outfile, output_type)
                
                netcdf_outfile = trim(outDir)//'fault.'//trim(proc_str)//'.nc'
                call netcdf_write_on_fault(netcdf_outfile)
            endif
            ! If exiting, write again the restart files.
            if ((stoptag == 1) .or. (it==nstep)) then
                netcdf_outfile = trim(outDir)//'disp.r.nc'
                output_type = 'disp'
                call netcdf_write(netcdf_outfile, output_type)
                
                netcdf_outfile = trim(outDir)//'fault.r.nc'
                call netcdf_write_on_fault(netcdf_outfile)
            endif

            ! Owner request (enhanced log reporting): one line every 100
            ! steps, mean-since-last-line, matching solveTimeLoopPETSc.f90's
            ! PROGRESS line (no ksp_its here -- MUMPS is a direct solver).
            if (mod(it, 100) == 0) then
                call cpu_time(endTime)
                write(*,'(A,I8,A,E12.5,A,E12.5,A,E12.5,A,F10.2,A,F10.4)') &
                    'PROGRESS step=', it, ' t=', time, ' dt=', dtev1, &
                    ' Vmax=', maxSlipRate, ' wall=', endTime - progressStartTime, &
                    ' s/step=', (endTime - progressStartTime) / 100.0d0
                progressStartTime = endTime
            endif
        endif
    enddo
    call cpu_time(endTime)
    timeUsedInComputing = endTime-startTime
    
    if (me == 0) then
        write(*,*) '====================================================================='
        write(*,*) '=                        RUN SUMMARY                                ='
        write(*,'(X,A,40X,i7,4X,A)')     '= MPI ranks                = ', nprocs, '='
        write(*,'(X,A,40X,i7,4X,A)')     '= Nodes                    = ', numnp, '='
        write(*,'(X,A,40X,i7,4X,A)')     '= Elements                 = ', numel, '='
        write(*,'(X,A,40X,i7,4X,A)')     '= Equations                = ', neq, '='
        ! Per fault, with its y, not just "nftnd(1)" -- that alone hid the
        ! belt bug: a second fault with zero nodes went unreported here
        ! while fault 1's count printed a perfectly plausible number and
        ! the run continued for 4240 steps regardless.
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
        write(*,*) '====================================================================='
        ! Owner request (enhanced log reporting): 'N/A'/-1 for MUMPS, which
        ! has no KSP concept -- the direct-solve counterpart to
        ! solveTimeLoopPETSc.f90's real kspType/pcType/iteration reporting.
        ! exit_reason distinguishes the loop ending because nstep was
        ! reached from ending early via exitCriteria's stoptag (it-1, not
        ! it, since a stoptag exit is caught at the TOP of the following
        ! iteration, before that step's body ever runs).
        if (it-1 < nstep) then
            call output_run_metadata(timeUsedInComputing, timeUsedInFactorization, &
                'N/A', 'N/A', 'stoptag', -1, 0.0d0, -1)
        else
            call output_run_metadata(timeUsedInComputing, timeUsedInFactorization, &
                'N/A', 'N/A', 'nstep_reached', -1, 0.0d0, -1)
            ! fdc hands this cycle's state to the dynamic code as a rupture;
            ! say plainly that it is not one (exit_reason is nstep_reached).
            if (eqquasi_mode == 2) then
                write(*,*) 'WARNING: eqquasi_mode 2 reached nstep without a rupture'
                write(*,*) '         (no slip rate above exit_slip_rate); raise par.nstep.'
            endif
        endif
    endif
    
end subroutine solveTimeLoopMUMPS
