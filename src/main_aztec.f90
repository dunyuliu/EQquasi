subroutine main_aztec

	use globalvar
	implicit none
	include 'mpif.h'
	include 'az_aztecf.h'

	integer(kind=8)::nev
	integer(kind=4)::i,j,inv,jnv,ntag,node_num,var,l,k
	real(kind=8)::tmparr(nftmx)
	
	!Parameters for Aztec
	integer(kind=4)::proc_config(0:AZ_PROC_SIZE),options(0:AZ_OPTIONS_SIZE),ierror,dimtmp, AZ_find_index
	integer(kind=4)::data_org(0:dimestimate),external(0:dimestimate),update_index(0:dimestimate), external_index(0:dimestimate), btmp(0:dimestimate)
	real(kind=8)::b(0:dimestimate),azx(0:dimestimate),untransx(0:dimestimate),btmpx(0:dimestimate)
	real(kind=8)::div,params(0:AZ_PARAMS_SIZE),status(0:AZ_STATUS_SIZE),u0(0:neq-1)
	integer(kind=4)::N1(0:nprocs-1),N2(0:nprocs-1),N_update_count(0:nprocs-1),N_update_disp(0:nprocs-1),N_update4u0(0:neq-1)

	call AZ_set_proc_config(proc_config, MPI_COMM_WORLD)
	allocate(val(0:dimestimate*81),bindx(0:dimestimate*81), update(0:dimestimate))
	
	tmparr = 0.0d0
	do i = 1,nftnd(1)
	tmparr(i) = fric(47,i,1)
	enddo         
	maxSlipRate = maxval(tmparr)
	
	if (me.eq.0) then
		!LOOP 2000 ![MUMPS]
		write(*,*) '=                                                                   ='
		write(*,*) '=       Building KSTIFF in CRS format                               ='
		write(*,*) '=                                                                   ='
		write(*,*) '=       Converting CRS to MRS format                                ='
	endif
	

	!Define problem on the host (processor 0)
	dtev = ksi*minDc/maxSlipRate
	dtev1 = ksi*minDc/maxSlipRate!Initial dtev1>dt to enter static state.
	cons = 0.0d0!Initialize cons, displacements.
	constmp = 0.0d0
	consvtmp = 0.0d0
	consv = 0.0d0
	consa = 0.0d0

	call msr
	
	div=neq/nprocs
	dimtmp = floor(div)+nprocs		
	
	if (me == 0) then
		write(*,*) 'Arrays for AZTEC including update, external ... are allocated.'
	endif				
	b=0.0d0
	azx=0.0d0
	untransx = 0.0d0
	update = 0
	external = 0
	update_index = 0
	external_index =0
	data_org = 0
	
	call AZ_read_update(N_update,update,proc_config,neq,1,0)
	write(*,*) 'me,N_update',me,N_update
	N2(0)=0
	do i=0,nprocs-1
		N1(i)=1
		if (i>0)then
			N2(i)=N2(i-1)+N1(i-1)
		endif
	enddo
	call MPI_Allgatherv(N_update,1,MPI_INT,N_update_count, N1, &
					N2,MPI_INT,MPI_COMM_WORLD,ierror)
	N_update_disp(0)=0
	do i=0,nprocs-1
		if (i>0)then
			N_update_disp(i)=N_update_disp(i-1)+N_update_count(i-1)
		endif
	enddo					
	if (me==0) then
		do i = 0,nprocs-1
			write(*,*) 'N_update_count/disp,i',N_update_count(i),N_update_disp(i),i
		enddo
	endif		
	
	bindx = 0
	val   = 0.0d0
	mass  = 0.0d0
	f     = 0.0d0

	if (me == 0) then
		write(*,*) 'For AZTEC, bindx and val are allocated.'
	endif
	!Check update should be in ascending order
	do i = 0, N_update-2
		if (update(i+1)<=update(i)) then 
			write(*,*) 'me =',me
			stop 'update is in a reverse order in me='
		endif
	enddo
	minneq = update(0)
	maxneq = update(N_update-1)
	write(*,*) 'me,minneq,maxneq', me, minneq, maxneq
	!Create bindx from lie and num
	bindx(0) = N_update+1
	ntag = N_update
	do i = 0, N_update-1
		bindx(i+1) = bindx(i) + num(update(i)+1)-1 !To exclude the diagonal term, which is included in num and lie.
		do j = 1,num(update(i)+1)
			if (lie(update(i)+1,j)/=(update(i)+1)) then 
				ntag = ntag + 1
				bindx(ntag) = lie(update(i)+1,j)-1!Pay attention, in AZTEC, arrays start from 0 dimension.
			endif
		enddo 	
	enddo
	if (me == 0) then
		write(*,*) 'Bindx is successfully filled.'
	endif	
	if (me == 19) then
		write(*,*) 'bindx 19 65885',N_update,bindx(65884),bindx(65885),bindx(65886)
	endif
	!LOOP 2001 : Calculating element stiff, mass, damp and constructing the overall stiffness
	if (me == 0) then
		write(*,*) 'LOOP 2001'
		write(*,*) 'Calculating element stiff, mass, damp and constructing the distributed global stiffness in progress ...'
	endif		
	!LOOP 2001 : Calculating element stiff, mass, damp and constructing the overall stiffness
	write(*,*) 'LOOP 2001'
	write(*,*) 'Calculating element stiff, mass, damp and constructing the overall stiffness.'
	
	call elemcal_aztec
		
	if (me==0) then 
		write(*,*) 'Val is successfully filled'
		write(*,*) 'Before calling Az_transform.'
	endif
	call Az_transform(proc_config,external,bindx,val,update, &
						update_index,external_index,data_org,  &
						N_update,0,0,0,0,AZ_MSR_MATRIX)
	!write(*,*) 'me =',me, 'azx needs',data_org(AZ_N_internal)+data_org(AZ_N_external)+data_org(AZ_N_border)
	!write(*,*) 'me =',me, 'dimestimate = ',dimestimate
	!Initialize AZTEC options
	if (me==0) then 
		write(*,*) 'After calling Az_transform.'
		write(*,*) 'Before calling Az_defaults.'
	endif	
	call AZ_defaults(options, params)
	if (AZTEC_OPTIONS == 2) then
		options(AZ_max_iter) = azmaxiter
		options(AZ_output) = AZ_last
		params(AZ_tol) = aztol
		params(AZ_drop) = 0.0d0
		options(AZ_solver) = AZ_gmres
		options(AZ_precond) = AZ_dom_decomp
		options(AZ_subdomain_solve) = AZ_ilu
		params(AZ_ilut_fill) = 0
	elseif (AZTEC_OPTIONS == 3) then
		params(AZ_tol) = 1.0e-06
		options(AZ_solver) = AZ_bicgstab
		options(AZ_precond) = AZ_sym_GS
		options(AZ_poly_ord) =10
		options(AZ_scaling) = AZ_sym_diag
		options(AZ_subdomain_solve) = AZ_ilut
		options(AZ_max_iter) = 1500
		params(AZ_ilut_fill) = 3
		params(AZ_drop) = 0.0d0	
	elseif (AZTEC_OPTIONS == 4) then
		options(AZ_solver) = AZ_cg
		!options(AZ_scaling) = AZ_none
		!options(AZ_precond) = AZ_ls
		options(AZ_output) = AZ_last
		options(AZ_max_iter) = 2000
		!options(AZ_poly_ord) = 7
		params(AZ_tol) = 1.0d-7
		!params(AZ_drop) = 0.0d0
		!params(AZ_omega) = 1.0d0
	endif
	if (me==0) then 
		write(*,*) 'After calling Az_defaults.'
		write(*,*) 'AZTEC_OPTIONS == ', AZTEC_OPTIONS 		
	endif	

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


   	do it=1,nstep 
		if (stoptag == 1) exit

		
		resu_1 = resu
		
		time=time+dtev1
		
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
		
		if (maxSlipRate>slipr_thres)then 
			if (status0 == 0) then        
				status1 = 1
			elseif (status0 == 1) then
				if (t_end_status>0.0d0) then
					t_end_status = -1.0d0
					tdynaend = -1000.0d0
				endif
			endif
		elseif (maxSlipRate<=slipr_thres) then
			if (status0 == 0) then
				status1 = 0
			elseif (status0 == 1) then
				if (t_end_status<0.0d0) then
					t_end_status = 1.0d0
					tdynaend = time
				else
					if ((time-tdynaend)>=100.0d0.and.tdynaend>0.0d0) then !exit
						stoptag = 1
						!call MPI_Bcast(stoptag, 1, MPI_INT, 0, MPI_COMM_WORLD, IERR)
					endif
					!if ((time-tdynaend)>=1.0d0*7*24*60*60.and.tdynaend>0.0d0) exit
				endif
			endif
		endif
		
		if ((status1-status0)==1) then
			tdynastart = time
			t_start_status = 1
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
		do i =0, N_update-1
			b(i) = right(update(i)+1)
		enddo
		!azx = 0.0d0 
		call AZ_reorder_vec(b,data_org,update_index,0)
		call AZ_solve(azx,b,options,params,0,bindx,0,0,0,val,data_org,status,proc_config)
		call AZ_invorder_vec(azx,data_org,update_index,0,untransx)
		! do i =0, N_update-1
			! resu(update(i)+1) = untransx(i)
		! enddo
		! To get the solution for the whole model.
		call MPI_Allgatherv(untransx,N_update,MPI_DOUBLE_PRECISION, &
					u0,N_update_count,N_update_disp,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierror)
		call MPI_Allgatherv(update,N_update,MPI_INT, &
					N_update4u0,N_update_count,N_update_disp,MPI_INT,MPI_COMM_WORLD,ierror)	
		do i =0, neq-1
			azx(i) = u0(i) ! Use the solution as the initial guess for the next solve. 
			resu(N_update4u0(i)+1) = u0(i)
		enddo	
		
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
	
		do i =0, N_update-1
			b(i) = right(update(i)+1)
		enddo
		!azx = 0.0d0
		call AZ_reorder_vec(b,data_org,update_index,0)
		call AZ_solve(azx,b,options,params,0,bindx,0,0,0,val,data_org,status,proc_config)
		call AZ_invorder_vec(azx,data_org,update_index,0,untransx)

		! To get the solution for the whole model.
		call MPI_Allgatherv(untransx,N_update,MPI_DOUBLE_PRECISION, &
					u0,N_update_count,N_update_disp,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierror)
		call MPI_Allgatherv(update,N_update,MPI_INT, &
					N_update4u0,N_update_count,N_update_disp,MPI_INT,MPI_COMM_WORLD,ierror)		

		do i =0, neq-1
			azx(i) = u0(i) ! Use the solution as the initial guess for the next solve. 
			resu(N_update4u0(i)+1) = u0(i)
		enddo
				

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
		globaldat(2,it) = maxSlipRate
		globaldat(3,it) = totMomRate	
		globaldat(4,it) = totTaoRuptArea
		globaldat(5,it) = totSlipRuptArea
		globaldat(6,it) = totRuptArea
	
	enddo!it
	
end subroutine main_aztec
