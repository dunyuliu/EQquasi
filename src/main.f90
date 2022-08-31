subroutine main

	use globalvar
	implicit none
	include 'mpif.h'
	include 'dmumps_struc.h'
	
	TYPE (DMUMPS_STRUC) mumps_par
	integer (kind = 8) :: nev
	integer (kind = 4) :: i,j,inv,jnv,ntag,node_num,var,l,k
	real (kind = dp) :: tmparr(nftmx)
	character (len = 50) :: netcdf_outfile, output_type
	
	![MUMPS]-2018/03
	!Define a communicator for the package
	mumps_par%COMM = MPI_COMM_WORLD
	!Initialize an instance of the package 
	!for LU factorization (sym=0, with working host)
	mumps_par%JOB = -1
	mumps_par%SYM = 1! 0: unsymmetric; 1: Positive definite symmetric ; 2: general symmetric
	mumps_par%PAR = 1
	CALL DMUMPS(mumps_par)
 
	tmparr = 0.0d0
	do i = 1,nftnd(1)
		tmparr(i) = fric(47,i,1)
	enddo         
	maxsliprate = maxval(tmparr) 
	if (me.eq.0) then
		!LOOP 2000 ![MUMPS]
		write(*,*) '=                                                                   ='
		write(*,*) '=       Building KSTIFF in CRS format                               ='
		write(*,*) '=                                                                   ='
		write(*,*) '=       Converting CRS to MUMPS format                              ='
		!Define problem on the host (processor 0)
		dtev = ksi*minDc/maxsliprate
		dtev1 = ksi*minDc/maxsliprate!Initial dtev1>dt to enter static state.
		cons = 0.0d0!Initialize cons, displacements.
		constmp = 0.0d0
		consvtmp = 0.0d0
		consv = 0.0d0
		consa = 0.0d0

		call crs

		!LOOP 2001 : Calculating element stiff, mass, damp and constructing the overall stiffness
		call elemcal
		
		mumps_par%N = neq ![MUMPS]
		mumps_par%NNZ = maxa
		ALLOCATE( mumps_par%IRN (mumps_par%NNZ))
		ALLOCATE( mumps_par%JCN (mumps_par%NNZ))
		ALLOCATE( mumps_par%A (mumps_par%NNZ))
		ALLOCATE( mumps_par%RHS (mumps_par%N))
		do i = 1, neq
			if (i >= 1.and. i<= neq) then
				do j = ia(i),ia(i+1)-1
					mumps_par%IRN (j) = i
				enddo
			endif
		enddo
		mumps_par%JCN = ja
		mumps_par%A = kstiff

		if (nftnd(1) > 0) then !RSF
			do l=1,nftnd(1)
			  if (icstart == 1) then 
				consv(1,nsmp(1,l,1)) =  -fric(46,l,1)/2.0d0!-mass(id(1,nsmp(2,l,1)))/(mass(id(1,nsmp(1,l,1)))+mass(id(1,nsmp(2,l,1)))) * fric(46,l,1) *us(1,l,1)!Slave
				consv(2,nsmp(1,l,1)) =  0.0d0!-mass(id(1,nsmp(2,l,1)))/(mass(id(1,nsmp(1,l,1)))+mass(id(1,nsmp(2,l,1)))) * fric(46,l,1) *us(2,l,1)
				consv(3,nsmp(1,l,1)) =  0.0d0!-mass(id(1,nsmp(2,l,1)))/(mass(id(1,nsmp(1,l,1)))+mass(id(1,nsmp(2,l,1)))) * fric(46,l,1) *us(3,l,1)
				consv(1,nsmp(2,l,1)) =  fric(46,l,1)/2.0d0!mass(id(1,nsmp(1,l,1)))/(mass(id(1,nsmp(1,l,1)))+mass(id(1,nsmp(2,l,1)))) * fric(46,l,1) *us(1,l,1)!Master 
				consv(2,nsmp(2,l,1)) =  0.0d0!mass(id(1,nsmp(1,l,1)))/(mass(id(1,nsmp(1,l,1)))+mass(id(1,nsmp(2,l,1)))) * fric(46,l,1) *us(2,l,1) 
				consv(3,nsmp(2,l,1)) =  0.0d0!mass(id(1,nsmp(1,l,1)))/(mass(id(1,nsmp(1,l,1)))+mass(id(1,nsmp(2,l,1)))) * fric(46,l,1) *us(3,l,1) 
				! consv(1,nsmp(1,l,1)) =  -mass(id(1,nsmp(2,l,1)))/(mass(id(1,nsmp(1,l,1)))+mass(id(1,nsmp(2,l,1)))) * fric(46,l,1) *us(1,l,1)!Slave
				! consv(2,nsmp(1,l,1)) =  -mass(id(1,nsmp(2,l,1)))/(mass(id(1,nsmp(1,l,1)))+mass(id(1,nsmp(2,l,1)))) * fric(46,l,1) *us(2,l,1)
				! consv(3,nsmp(1,l,1)) =  -mass(id(1,nsmp(2,l,1)))/(mass(id(1,nsmp(1,l,1)))+mass(id(1,nsmp(2,l,1)))) * fric(46,l,1) *us(3,l,1)
				! consv(1,nsmp(2,l,1)) =  mass(id(1,nsmp(1,l,1)))/(mass(id(1,nsmp(1,l,1)))+mass(id(1,nsmp(2,l,1)))) * fric(46,l,1) *us(1,l,1)!Master 
				! consv(2,nsmp(2,l,1)) =  mass(id(1,nsmp(1,l,1)))/(mass(id(1,nsmp(1,l,1)))+mass(id(1,nsmp(2,l,1)))) * fric(46,l,1) *us(2,l,1) 
				! consv(3,nsmp(2,l,1)) =  mass(id(1,nsmp(1,l,1)))/(mass(id(1,nsmp(1,l,1)))+mass(id(1,nsmp(2,l,1)))) * fric(46,l,1) *us(3,l,1) 
				! write(*,*) 'l,',l,consv(1,nsmp(1,l,1)),consv(2,nsmp(1,l,1)),consv(3,nsmp(1,l,1))
				! write(*,*) 'mass,',mass(id(1,nsmp(2,l,1))),fric(46,l,1),us(1,l,1)
			   else
				consv(1,nsmp(1,l,1)) = fric(34,l,1)
				consv(2,nsmp(1,l,1)) =  fric(35,l,1)
				consv(3,nsmp(1,l,1)) =  fric(36,l,1)
				consv(1,nsmp(2,l,1)) =  fric(31,l,1)
				consv(2,nsmp(2,l,1)) =  fric(32,l,1)
				consv(3,nsmp(2,l,1)) =  fric(33,l,1)
			   endif 	
			enddo	
		endif

		itag=0 ! label the step needed written in files
	endif!END LOOP 2000 ![MUMPS]


	mumps_par%JOB = 4! 1 + 2
	CALL DMUMPS(mumps_par)
	
	stoptag = 0 ! set stoptag to FALSE.
	
   	do it = 1, nstep 
		if (stoptag == 1) exit ! exit EQquasi if stoptag is TRUE.

		if (me.eq.0) then
		
			resu_1 = resu
			
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
				write(*,'(X,A,40X,E15.7,4X,A)') '=',  maxsliprate, 'm/s'
			endif			
			
			call exit_criteria ! determine if to exit EQquasi.

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
			itag = 0
			
			!COMPUTE U*(t+1) FOR THE WHOLE VOLUME'
			call bound_load

			mumps_par%RHS = right
		endif
		call MPI_BCAST(stoptag, 1, MPI_INT, 0, MPI_COMM_WORLD, IERR)

		mumps_par%JOB = 3
		mumps_par%ICNTL(3) = -1 ! Goblal info, default 6 
		CALL DMUMPS( mumps_par )!MPI
		
		if (me.eq.0) then
			resu = mumps_par%RHS
		
			do i=1,numnp
				do j=1,ndof
					if (id(j,i)>0) then !Only nodes that have equation number have been updated here.
						constmp(j,i)=resu(id(j,i))					
					endif
				enddo 
			enddo 	
			!COMPUTE F*=KU*(t+1)=[KSTIFF].DOT.[CONSTRAINTMP]'
			!COMPUTE LUMPED MASS	FOR SPLIT NODES	

			call bound_ft_ku
			
			!GET FIRST PREDICTIONS OF V*(t+1).

			call faulting

			!V*(t+1) OBTAINED. 
			!UPDATE BOUNDARY U**(t+1),THEN DECLARE [CONSTRAIN] BY U**(t+1).	
			
			itag = 1
			
			call bound_load
				
			!COMPUTE U**(t+1) FOR THE WHOLE VOLUME. <-> REPEAT 3.5.2'
		
			mumps_par%RHS = right
		endif!MYID==0
		
		mumps_par%JOB = 3
		CALL DMUMPS( mumps_par )!MPI
		
		if (me.eq.0) then
			resu = mumps_par%RHS
		
			do i=1,numnp
				do j=1,ndof
					if (id(j,i)>0) then 
						cons(j,i) = resu(id(j,i))
					else
						cons(j,i) = constmp(j,i)
					endif
				enddo 
			enddo
			
			!COMPUTE F**=KU**(t+1)=[KSTIFF].DOT.[CONSTRAIN]'		
			call bound_ft_ku
			
			!GET SECOND PREDICTIONS OF V**(t+1), AND DECLARE V(t+1)=V**(t+1).'	
			call faulting

			status0=status1
			nev=int8(dtev/dt)

			dtev1=max(dt,nev*dt)	

			do i=1,numnp
				do j=1,ndof
					if (id(j,i)>0) then 
						consv(j,i) = (resu(id(j,i)) - resu_1(id(j,i)))/dtev1
					endif
				enddo 
			enddo						
			globaldat(1,it) = time
			globaldat(2,it) = maxsliprate
			globaldat(3,it) = totmomrate	
			globaldat(4,it) = tottaoruptarea
			globaldat(5,it) = totslipruptarea
			globaldat(6,it) = totruptarea
			
			! Write disp field. 
			if (mod(it,nt_output_stress) == 1) then 
				write(proc_str,'(I5.5)') it
				netcdf_outfile = 'disp.'//trim(proc_str)//'.nc'
				output_type = 'disp'
				call netcdf_write(netcdf_outfile, output_type)
				
				netcdf_outfile = 'fault.'//trim(proc_str)//'.nc'
				call netcdf_write_on_fault(netcdf_outfile)
			endif
		endif!MYID==0
	enddo!it

	!if (me.eq.0) then		
	!	DEALLOCATE( mumps_par%IRN )
	!	DEALLOCATE( mumps_par%JCN )
	!	DEALLOCATE( mumps_par%A )
	!	DEALLOCATE( mumps_par%RHS )
	!endif
	
	!Destroy the instance (deallocate internal data structures)
	!mumps_par%JOB = -2
	!CALL DMUMPS( mumps_par )	
	
	!CALL MPI_FINALIZE ( IERR )
	
end subroutine main

subroutine exit_criteria
! Subroutine to compute if to exit EQquasi. 
! If mode == 1 (quasi-dynamic/quasi-static), the exit happens AFTER 
!	a rupture when maxsliprate falls below slipr_thres.
! If mode == 2 (fully dynamic), the exit happens BEFORE a rupture when 
!	the maxsliprate rises above the slipr_thres.

! The system is determined by the combination of two status variables - status0/1 and t_end_status.
! status0/status1: code status in the last/current time step.
! 0: interseismic; 1: rupture. 
! 
	use globalvar
	implicit none
	
	if (eqquasi_mode == 1) then ! quasi-dynamic/quasi-static
		if (maxsliprate>slipr_thres)then
		! if maximum slip rate rises above the slipr_thres, entering/in the co-seismic phase.
			if (status0 == 0) then ! convert from inter-seismic to co-seismic.
				status1 = 1 ! current status is changed to co-seismic.
			elseif (status0 == 1) then ! convert from co-seismic to inter-seismic.
				if (t_end_status>0.0d0) then ! if t_end_status is TRUE, set it to FALSE and continue rupturing.
					t_end_status = -1.0d0 ! t_end_status is FALSE.
					tdynaend = -1000.0d0 ! don't need to record end time in the rupture phase.
				endif
			endif
		elseif (maxsliprate<=slipr_thres) then
		! if maximum slip rate falls below the slipr_thres, consider exiting the code. 
			if (status0 == 0) then 
				status1 = 0 ! if last step was in inter-seismic, no need to change.
			elseif (status0 == 1) then !if last step was in co-seismic, consider to exit.
				if (t_end_status<0.0d0) then 
					t_end_status = 1.0d0 ! change exit status to be TRUE.
					tdynaend = time ! record the end time of the rupture.
				else
					if ((time-tdynaend)>=100.0d0.and.tdynaend>0.0d0) then ! exit if the low slip rate status is kept for over 100 seconds.
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
		if (maxsliprate>slipr_thres)then
		! In this mode, if maximum slip rate rises above the slipr_thres, exit. 
			if (status0 == 0) then ! convert from inter- to co-seismic.
				status1 = 1
				stoptag = 1 ! stoptag is TRUE now. Direct the code to exit.
			endif
		endif
	endif 
end subroutine exit_criteria