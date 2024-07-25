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
    write(*,*) '==================  Welcome to EQquasi 1.3.3  ======================='
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
    
    call readcurrentcycle 
    call readmodel
    call readstations1

    allocate(an4nds(2,n4nds), xonfs(2,maxval(nonfs),ntotft), x4nds(3,n4nds))
    call readstations2

    call mesh4num
    call allocAndInit
    call meshgen
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

subroutine checkAndReport(currentProcID)
    use globalvar
    implicit none 
    logical :: file_exists
    character (len = 50) :: filenametmp, output_type
    
    integer (kind = 4) :: currentProcID
    if (icstart == 1) then
    INQUIRE(FILE="on_fault_vars_input.nc", EXIST=file_exists)
        if (file_exists .eqv. .FALSE.) then
            write(*,*) 'on_fault_vars_input.nc is required but missing ...'
        endif 
    
        call netcdf_read_on_fault("on_fault_vars_input.nc")
    else
        INQUIRE(FILE="on_fault_vars_input.nc", EXIST=file_exists)
            if (file_exists .eqv. .FALSE.) then
                write(*,*) 'on_fault_vars_input.nc is required but missing ...'
            endif 
        
        INQUIRE(FILE="fault.r.nc", EXIST=file_exists)
            if (file_exists .eqv. .FALSE.) then
                write(*,*) 'icstart>1, fault.r.nc is required but missing ...'
            endif
        
        call netcdf_read_on_fault_restart("on_fault_vars_input.nc", "fault.r.nc")
    endif 
    
    if (currentProcID == 0) then 
        write(*,'(X,A,40X,i7,4X,A)') '= Total nodes = ',numnp,'='            
        write(*,'(X,A,40X,i7,4X,A)') '= Total elems = ',numel,'='            
        
        filenametmp = 'eqquasi.mesh.coor.nc'
        output_type = 'coor'
        call netcdf_write(filenametmp, output_type)
        
        filenametmp = 'eqquasi.mesh.ien.nc'
        output_type = 'ien'
        call netcdf_write(filenametmp, output_type)
        
        filenametmp = 'eqquasi.mesh.nsmp.nc'
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
                slp4fri(nftmx,ntotft),fltslp(3,nftmx,ntotft), globaldat(10,nstep), fltsta(10,nstep,nonmx))

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
        call output_ruptarea_trac_slip
        write(*,*) 'Done.'
    endif 
end subroutine writeResults