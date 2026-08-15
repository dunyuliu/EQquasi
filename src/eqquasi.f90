!/* Copyright (C) 2018-2023, Earthquake Modeling Lab @ Texas A&M University. 
! * and Institute for Geophysics, the University of Texas at Austin.
! * All Rights Reserved.
! * This code is part of software EQquasi, please see EQquasi License Agreement
! * attached before you copy, download, install or use EQquasi./

program eqquasi3d
    use globalvar
    implicit none
    include 'mpif.h'
    integer (kind = 4) :: i, l, k, j
    
    call MPI_INIT(IERR)
    call mpi_comm_rank(MPI_COMM_WORLD,me,IERR)
    call mpi_comm_size(MPI_COMM_WORLD,nprocs,IERR)
    
    if (me == 0) then 
    write(*,*) '====================================================================='
    write(*,*) '==================  Welcome to EQquasi '//EQQUASI_VERSION//'  ======================='
    write(*,*) '=====   Product of UTIG & Earthquake Modeling Lab@TAMU          ====='
    write(*,*) '========== GitHub: https://github.com/dunyuliu/EQquasi  ============='
    write(*,*) '=========== Contacts: dliu@ig.utexas.edu                ============='
    write(*,*) '=                                                                   ='
    write(*,*) '=   EQquasi is a parallel finite-element software to simulate       ='
    write(*,*) '=   earthquake cycles on geometrically realistic fault systems.     ='
    write(*,*) '=   It is part of the fully dynamic earthquake simulator            ='
    write(*,*) '=   EQsimu https://github.com/dunyuliu/EQsimu.                      ='
    write(*,*) '=                                                                   ='
    write(*,*) '=   All adjustable parameters are in user_defined_params.py.        ='
    write(*,*) '=   Simply run ./case.setup after each modification of              ='
    write(*,*) '=       user_defined_params.py                                      ='
    write(*,*) '====================================================================='
    endif 
    
    call setRunDate
    call set_io_dirs
    call readcurrentcycle 
    call readmodel
    call readstations1

    allocate(an4nds(2,n4nds), xonfs(2,maxval(nonfs),ntotft), x4nds(3,n4nds))
    call readstations2

    call mesh4num
    call allocAndInit
    call meshgen

    ! A requested station that matches no fault node used to vanish without a
    ! word. BP1002 asked for three and got one: two sat at x = +-2.0 km and
    ! z = -2.0/-18.0 km, none of which is a multiple of dx = 2500 m, so they
    ! matched nothing and the run reported three requested, one written, and
    ! no discrepancy. The reference was blessed from it.
    if (me == 0 .and. sum(nonfs) > 0 .and. n4onf < sum(nonfs)) then
        write(*,*) '====================================================='
        write(*,*) '= ON-FAULT STATIONS REQUESTED BUT NOT FOUND         ='
        write(*,*) '=   requested (stations.txt) =', sum(nonfs)
        write(*,*) '=   matched to a fault node  =', n4onf
        write(*,*) '= A station only exists if its (x, z) lands on a    ='
        write(*,*) '= node: both must be multiples of dx, measured from ='
        write(*,*) '= the fault corner. Fix par.st_coor_on_fault in     ='
        write(*,*) '= user_defined_params.py and rerun case.setup.      ='
        write(*,*) '====================================================='
    endif

    call checkAndReport(me)

    if(n4out>0) then 
        ndout=n4out*ndof*noid !3 components of 2 quantities: v and d
        allocate(idhist(3,ndout),dout(ndout+1,nstep+1))

        idhist=0
        dout=0.0d0
        l=0
        do i=1,n4out
            do j=1,ndof
                do k=1,noid
                    l = l + 1
                    idhist(1,l) = an4nds(2,i) !node number (>1, <NUMNP)
                    if(idhist(1,l)<=0) idhist(1,l)=1  !avoid zero that cannot be in array below
                    idhist(2,l) = j    !degree of freedom number (<=NDOF)
                    idhist(3,l) = k    !kinematic quantity specifier 
                    !(disp, vel, or acc)
                enddo
            enddo
        enddo            
    endif
    
    if (sol_op == 1) then 
        call solveTimeLoopMUMPS
    elseif (sol_op == 2) then 
        write(*,*) 'aztec is temporarily disabled.'
        !call main_aztec
    endif
    
    call writeResults(me)
    call MPI_FINALIZE(IERR)
    
end program eqquasi3d

subroutine setRunDate
! Stamp the run date into runDate, for the benchmark output headers.
    use globalvar
    implicit none
    integer (kind = 4) :: dtvals(8)
    call date_and_time(values = dtvals)
    write(runDate,'(i4.4,A,i2.2,A,i2.2)') dtvals(1), '/', dtvals(2), '/', dtvals(3)
end subroutine setRunDate

subroutine checkNormalStressCaps
    ! Stop if the initial normal stress lies outside the caps that will be
    ! applied to it.
    !
    ! The caps only act when rough_fault == 1 .and. C_elastic == 1. When they
    ! do, they clamp fric(FR_TNRM0,...) on the first evaluation -- but the initial
    ! SHEAR stress was computed by the compset from the uncapped normal
    ! stress, and nothing recomputes it. The initial condition is then no
    ! longer at steady state: tau/sigma_bar is too high by exactly the ratio
    ! the clamp introduced, and rate-and-state turns that into slip rate
    ! exponentially. liu2020.kink.qdc sets -50 MPa against a -40 MPa cap, a
    ! 25 per cent overshoot, and reached 2e-1 m/s on step 2 from a 1e-9 m/s
    ! start -- exp(0.15/0.007) ~ 2e9, which is what a = 0.007 does with it.
    !
    ! Silently clamping is the wrong answer here: the run continues and
    ! produces numbers. Refuse instead, and say which cap and by how much.
    use globalvar
    implicit none
    integer (kind = 4) :: ift, i
    real (kind = dp)   :: lo, hi

    if (.not. (rough_fault == 1 .and. C_elastic == 1)) return

    lo = min(min_norm, max_norm)      ! most compressive allowed, e.g. -40 MPa
    hi = max(min_norm, max_norm)      ! least compressive allowed, e.g. -10 MPa

    do ift = 1, ntotft
        do i = 1, nftnd(ift)
            if (fric(FR_TNRM0,i,ift) < lo .or. fric(FR_TNRM0,i,ift) > hi) then
                write(*,*) '====================================================='
                write(*,*) '= Initial normal stress lies outside the caps       ='
                write(*,*) '=                                                   ='
                write(*,'(A,E12.4,A)') ' =   initial sigma_n : ', fric(FR_TNRM0,i,ift), ' Pa'
                write(*,'(A,E12.4,A)') ' =   cap range       : ', lo, ' Pa (most compressive)'
                write(*,'(A,E12.4,A)') ' =                     ', hi, ' Pa (least compressive)'
                write(*,'(A,I0,A,I0)') ' =   first at fault ', ift, ', node ', i
                write(*,*) '=                                                   ='
                write(*,*) '= The caps would clamp this on the first step while ='
                write(*,*) '= the initial SHEAR stress stays as the compset      ='
                write(*,*) '= computed it, so the initial condition would no     ='
                write(*,*) '= longer be steady state and slip rate would run     ='
                write(*,*) '= away exponentially. Refusing rather than running.  ='
                write(*,*) '=                                                   ='
                write(*,*) '= Fix: widen par.min_norm/par.max_norm to cover the ='
                write(*,*) '= initial stress, or lower par.init_norm into range. ='
                write(*,*) '= The caps act only when rough_fault=1 and           ='
                write(*,*) '= C_elastic=1; most cases never reach them.          ='
                write(*,*) '====================================================='
                stop 7
            endif
        enddo
    enddo

end subroutine checkNormalStressCaps

subroutine checkAndReport(currentProcID)
    use globalvar
    implicit none 
    logical :: file_exists
    character (len = 50) :: filenametmp, filenametmp2, output_type
    
    integer (kind = 4) :: currentProcID
    if (icstart == 1) then
    INQUIRE(FILE=trim(inputDir)//"on_fault_vars_input.nc", EXIST=file_exists)
        if (.not. file_exists) then
            write(*,*) 'on_fault_vars_input.nc is required but missing ...'
        endif 
    
        ! Assign into a fixed-length local before the call. These are
        ! external subroutines with implicit interfaces, so an expression
        ! argument passes a temporary whose length the callee does not know --
        ! it declares character(len=50) and reads past the end, producing a
        ! path with 22 bytes of adjacent memory appended.
        filenametmp = trim(inputDir)//"on_fault_vars_input.nc"
        call netcdf_read_on_fault(filenametmp)
    else
        INQUIRE(FILE=trim(inputDir)//"on_fault_vars_input.nc", EXIST=file_exists)
            if (.not. file_exists) then
                write(*,*) 'on_fault_vars_input.nc is required but missing ...'
                stop 6
            endif 
        
        INQUIRE(FILE=trim(outDir)//"fault.r.nc", EXIST=file_exists)
            if (.not. file_exists) then
                write(*,*) 'icstart>1, fault.r.nc is required but missing ...'
                stop 6
            endif
        
        filenametmp = trim(inputDir)//"on_fault_vars_input.nc"
        filenametmp2 = trim(outDir)//"fault.r.nc"
        call netcdf_read_on_fault_restart(filenametmp, filenametmp2)
    endif
    ! Both branches above load fric(FR_TNRM0,...), so the caps check runs for every
    ! cycle. It used to sit inside the icstart==1 branch only, which left
    ! cycles 2..N of a multi-cycle run unchecked -- "refuses instead of
    ! silently clamping" held at cold start alone.
    call checkNormalStressCaps 
    
    if (currentProcID == 0) then 
        write(*,'(X,A,40X,i7,4X,A)') '= Total nodes = ',numnp,'='            
        write(*,'(X,A,40X,i7,4X,A)') '= Total elems = ',numel,'='            
        
        filenametmp = trim(outDir)//'eqquasi.mesh.coor.nc'
        output_type = 'coor'
        call netcdf_write(filenametmp, output_type)
        
        filenametmp = trim(outDir)//'eqquasi.mesh.ien.nc'
        output_type = 'ien'
        call netcdf_write(filenametmp, output_type)
        
        filenametmp = trim(outDir)//'eqquasi.mesh.nsmp.nc'
        output_type = 'nsmp'
        call netcdf_write(filenametmp, output_type)
     
    endif 
    
end subroutine checkAndReport
    
subroutine allocAndInit
    use globalvar 
    implicit none 
        
    allocate(x(ndof,numnp),id(ndof,numnp),ien(nen,numel),mat(numel,5), &
        et(numel), eledet(numel), ia(neq+1), num(neq), &
        mass(neq), f(neq), right(neq), resu(neq), resu_1(neq), &
        cons(ndof,numnp),constmp(ndof,numnp),consv(ndof,numnp), &
        consvtmp(ndof,numnp),consa(ndof,numnp),consf(ndof,numnp),consm(ndof,numnp))
    x        = 0.0d0 
    id       = 0
    ien      = 0
    mat      = 0
    et       = 0
    eledet   = 0.0d0
    ia       = 0
    num      = 0
    mass     = 0.0d0 
    f        = 0.0d0 
    right    = 0.0d0 
    resu     = 0.0d0 
    resu_1   = 0.0d0 
    cons     = 0.0d0 
    constmp  = 0.0d0
    consv    = 0.0d0 
    consvtmp = 0.0d0 
    consa    = 0.0d0 
    consf    = 0.0d0 
    consm    = 0.0d0 
    
    dimestimate = floor(real(neq)/nprocs)*2    
    nftmx=maxval(nftnd) !max fault nodel num for all faults, used for arrays.
    if(nftmx<=0) nftmx=1  !fortran arrays cannot be zero size,use 1 for 0
    nonmx=sum(nonfs)    !max possible on-fault stations number
    allocate(nsmp(2,nftmx,ntotft),fnft(nftmx,ntotft),un(3,nftmx,ntotft),&
                us(3,nftmx,ntotft),ud(3,nftmx,ntotft),fric(100,nftmx,ntotft),&
                arn(nftmx,ntotft),r4nuc(nftmx,ntotft),anonfs(3,nonmx),&
                slp4fri(nftmx,ntotft),fltslp(3,nftmx,ntotft), globaldat(10,nstep), fltsta(14,nstep,nonmx))

    fnft    = -1000.d0!Should be initialized over 600.
    fric    = 0.0d0
    un      = 0.0d0
    us      = 1000.0d0
    ud      = 0.0d0
    arn     = 0.0d0
    r4nuc   = 0.0d0
    anonfs  = 0
    slp4fri = 0.0d0
    fltslp  = 0.0d0
    fltsta  = 0.0d0
    status0 = 0
    status1 = 0
    fltsta  = 0.0d0
    globaldat     = 0.0d0
    totMomRate    = 0.0d0    
    tdynastart    = -1000.0d0
    tdynaend      = -1000.0d0
    tcheck        = 0.0d0
    t_end_status  = -1.0d0
    t_start_status  = -1.0d0
    
end subroutine allocAndInit

subroutine writeResults(currentProcID)
    integer (kind = 4) :: currentProcID
    if (currentProcID == 0) then 
        call output_onfault_st
        call output_offfault_st
        call output_onfault_transfer 
        call output_timedy
        call output_globaldat
        call output_bp8_profiles
        call output_ruptarea_trac_slip
        write(*,*) 'Done.'
    endif 
end subroutine writeResults