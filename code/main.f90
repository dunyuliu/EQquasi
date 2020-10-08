subroutine main

	use globalvar
	implicit none
	include 'mpif.h'
	include 'dmumps_struc.h'
	
	TYPE (DMUMPS_STRUC) mumps_par
	integer(kind=8)::nev
	integer(kind=4)::i,j,inv,jnv,ntag,node_num,var,l,k
	real(kind=8)::tmparr(nftmx)
	
	![MUMPS]-2018/03
	!Define a communicator for the package
	mumps_par%COMM = MPI_COMM_WORLD
	!Initialize an instance of the package 
	!for LU factorization (sym=0, with working host)
	mumps_par%JOB = -1
	mumps_par%SYM = 1! 0: unsymmetric; 1: Positive definite symmetric ; 2: general symmetric
	mumps_par%PAR = 1
	CALL DMUMPS(mumps_par)

	if (icstart > 1) then 
		tmparr = 0.0d0
		do i = 1,nftnd(1)
		tmparr(i) = fric(47,i,1)
		enddo         
		maxsliprate = maxval(tmparr)
		!write(*,*) 'Maxsliprate from last cycle = ',maxsliprate
	elseif (icstart == 1) then 
		maxsliprate = 1.0d-2 
		!write(*,*) 'Maxsliprate from last cycle = ',maxsliprate
	endif 
	

	
	if (me.eq.0) then
		!LOOP 2000 ![MUMPS]
		write(*,*) 'LOOP 2000: CONSTRUCT KSTIFF in CRS; CONVERT CRS to MUMPS FORMAT'
		!Define problem on the host (processor 0)
		dtev = ksi*LL/maxsliprate
		dtev1 = ksi*LL/maxsliprate!Initial dtev1>dt to enter static state.
		cons = 0.0d0!Initialize cons, displacements.
		constmp = 0.0d0
		consvtmp = 0.0d0
		consv = 0.0d0
		consa = 0.0d0

		call crs

		!LOOP 2001 : Calculating element stiff, mass, damp and constructing the overall stiffness
		write(*,*) 'LOOP 2001'
		write(*,*) 'Calculating element stiff, mass, damp and constructing the overall stiffness.'

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

    do it=1,nstep 
		if (me.eq.0) then
		
			resu_1 = resu
			
			time=time+dtev1
			
			if (mod(it,nhplt) == 1) then
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
			
			if (maxsliprate>slipr_thres)then 
				if (status0 == 0) then        
					status1 = 1
				elseif (status0 == 1) then
					if (tcheckstatus>0.0d0) then
						tcheckstatus = -1.0d0
						tdynaend = -1000.0d0
					endif
				endif
			elseif (maxsliprate<=slipr_thres) then
				if (status0 == 0) then
					status1 = 0
				elseif (status0 == 1) then
					if (tcheckstatus<0.0d0) then
						tcheckstatus = 1.0d0
						tdynaend = time
					else
						if ((time-tdynaend)>=100.0d0.and.tdynaend>0.0d0) exit
						!if ((time-tdynaend)>=1.0d0*7*24*60*60.and.tdynaend>0.0d0) exit
					endif
				endif
			endif
			
			if ((status1-status0)==1) then
				tdynastart = time
			endif

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
			
		endif!MYID==0				
			
		if (me.eq.0) then

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
	
		endif!MYID==0
	enddo!it

	if (me.eq.0) then
		call output_onfault_st
		
		call output_offfault_st
		
		call output_onfault_transfer
		
		call output_timedy
		
		call output_globaldat
		
		DEALLOCATE( mumps_par%IRN )
		DEALLOCATE( mumps_par%JCN )
		DEALLOCATE( mumps_par%A )
		DEALLOCATE( mumps_par%RHS )
	endif
	
	!Destroy the instance (deallocate internal data structures)
	! mumps_par%JOB = -2
	! CALL DMUMPS( mumps_par )	
	
	! CALL MPI_FINALIZE ( IERR )
	
end subroutine main
