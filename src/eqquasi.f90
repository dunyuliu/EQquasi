!/* Copyright (C) 2018-2021, Earthquake Modeling Lab @ Texas A&M University. 
! * All Rights Reserved.
! * This code is part of software EQquasi, please see EQquasi License Agreement
! * attached before you copy, download, install or use EQquasi./

program eqquasi3d

	use globalvar
	implicit none
	include 'mpif.h'

	integer (kind = 4) :: itmp, i, l, k, j
	real (kind = dp) :: tmp
	character (len = 50) :: filenametmp, output_type

	CALL MPI_INIT(IERR)
	call mpi_comm_rank(MPI_COMM_WORLD,me,IERR)
	call mpi_comm_size(MPI_COMM_WORLD,nprocs,IERR)
	
	if (me == 0) then 
	write(*,*) '====================================================================='
	write(*,*) '==================  Welcome to EQquasi 1.3.2  ======================='
	write(*,*) '===== Product of UTIG & Earthquake Modeling Lab @ TAMU          ====='
	write(*,*) '========== Website https://seismotamu.wixsite.com/emlam ============='
	write(*,*) '=========== Contacts: dliu@ig.utexas.edu                ============='
	write(*,*) '=                                                                   ='
	write(*,*) '=   EQquasi uses FEM to simulate earthquake dynamic ruptures        ='
	write(*,*) '=   on geometrically realistic fault systems for quasi-static/      ='
	write(*,*) '=   quasi-dynamic deformation. It is part of the fully dynamic      ='
	write(*,*) '=   earthquake simulator EQsimu.                                    ='
	write(*,*) '=                                                                   ='
	write(*,*) '=   Model and system related parameters can be adjusted in          ='
	write(*,*) '=       model.txt,                                                  ='
	write(*,*) '=       fric.txt,                                                   ='
	write(*,*) '=       stations.txt,                                               ='
	write(*,*) '=                                                                   ='
	write(*,*) '====================================================================='
	endif 
	
	call readcurrentcycle 
	call readmodel
	!call readfric
	if (rough_fault == 1) then 
		call read_fault_rough_geometry 
		filenametmp = 'roughness.nc'
		call netcdf_write_roughness(filenametmp)
	endif 
	allocate(nonfs(ntotft))
	
	call readstations1
	itmp = maxval(nonfs)

	allocate(an4nds(2,n4nds), xonfs(2,itmp,ntotft), x4nds(3,n4nds))
	call readstations2

	fltxyz(1,1,1)=xmin
	fltxyz(2,1,1)=xmax
	fltxyz(1,2,1)=0.0d0
	fltxyz(2,2,1)=0.0d0
	fltxyz(1,3,1)=zmin
	fltxyz(2,3,1)=zmax
	fltxyz(1,4,1)=270.0d0/180.0d0*pi
	fltxyz(2,4,1)=90.0d0/180.0d0*pi

	dt = min(0.5d0*dx/mat0(1,1), 0.06d0) ! minimum time step size based on CFL criteria with alpha = 0.5.
	
	!...time hostories output steps
	nplpts = 0	!initialize number of time history plot
	if (nhplt > 0) then
		nplpts = int(nstep/nhplt) + 2
	endif

	allocate(nftnd(ntotft))
	
	call mesh4num

	allocate(x(ndof,numnp),id(ndof,numnp),ien(nen,numel),mat(numel,5),et(numel), eledet(numel), ia(neq+1), num(neq), &
			mass(neq), f(neq), right(neq), resu(neq), resu_1(neq), cons(ndof,numnp),constmp(ndof,numnp),consv(ndof,numnp),consvtmp(ndof,numnp),consa(ndof,numnp),consf(ndof,numnp),consm(ndof,numnp))
	x        = 0.0d0 
	id       = 0
	ien      = 0
	mat      = 0
	et       = 0
	eledet   = 0.0d0
	ia       = 0
	num      = 0
	mass     = 0.0d0 
	f        =0.0d0 
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
	
	tmp = neq/nprocs
	dimestimate = floor(tmp)*2	
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
    totmomrate    = 0.0d0	
	tdynastart    = -1000.0d0
	tdynaend      = -1000.0d0
	tcheck        = 0.0d0
	t_end_status  = -1.0d0
	t_start_status  = -1.0d0
	
	call meshgen
	
	if (icstart == 1) then 
		call netcdf_read_on_fault("on_fault_vars_input.nc")
	else
		call netcdf_read_on_fault_restart("on_fault_vars_input.nc", "fault.r.nc")
	endif 
	
	! Write out mesh information
	! coor, ien, and nsmp.
	if (me == 0) then 
		write(*,*) '=     Total nodes =                                                 ='
		write(*,'(X,A,40X,i7,4X,A)') '=',  numnp			
		write(*,*) '=     Total elements =                                              ='
		write(*,'(X,A,40X,i7,4X,A)') '=',  numel			
		
		filenametmp = 'eqquasi.mesh.coor.nc'
		output_type = 'coor'
		call netcdf_write(filenametmp, output_type)
		write(*,*) '=                                                                   ='
		write(*,*) '=       Writing out mesh.coor.nc                                    ='
		
		filenametmp = 'eqquasi.mesh.ien.nc'
		output_type = 'ien'
		call netcdf_write(filenametmp, output_type)
		write(*,*) '=                                                                   ='
		write(*,*) '=       Writing out mesh.ien.nc                                     ='
		
		filenametmp = 'eqquasi.mesh.nsmp.nc'
		output_type = 'nsmp'
		call netcdf_write(filenametmp, output_type)
		write(*,*) '=                                                                   ='
		write(*,*) '=       Writing out mesh.nsmp.nc                                    ='
	endif 
	
	if(n4out>0) then 
		ndout=n4out*ndof*noid !3 components of 2 quantities: v and d
		!   write(*,*) 'ndout= ',ndout    
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
					idhist(2,l) = j	!degree of freedom number (<=NDOF)
					idhist(3,l) = k	!kinematic quantity specifier 
					!(disp, vel, or acc)
				enddo
			enddo
		enddo			
	endif
	
	if (sol_op == 1) then 
		call main
	elseif (sol_op == 2) then 
		call main_aztec
	endif
	
	if (me == 0) then
		call output_onfault_st
		
		call output_offfault_st
		
		! phasing out
		!call output_onfault_transfer 
		
		call output_timedy
		
		call output_globaldat
		
		call output_ruptarea_trac_slip
	endif
	
	CALL MPI_FINALIZE ( IERR )
	
	write(*,*) '====================================================================='		
	write(*,*) '=       Job done, exiting...                                        ='
	write(*,*) '====================================================================='	

end
