subroutine CRSelem(nelement,nnode,neq,x,mat,ien,id,&
			nftnd,n4onf,xonfs,nonfs,&
			nftmx,nonmx,nsmp,un,us,ud,fric,arn,r4nuc,&
			anonfs,fnft,slp4fri,et,me)
use globalvar
implicit none

include 'mpif.h'
include 'dmumps_struc.h'
TYPE (DMUMPS_STRUC) mumps_par
integer(kind=4)::IERR,me

character(len=30)::sttmp,dptmp
integer(kind=4)::i,j,mm,m,nel,inv,jnv,&
	ntag,maxa,mloc,ki,kj,kloc,&
	temp,node_num,var,&
	it,itag,ift,l
integer(kind=4)::nelement,nnode,neq
integer(kind=4)::id(ndof,nnode),ien(nen,nelement),ia(neq+1),&
	num(neq),et(nelement)
integer(kind=4),allocatable::lie(:,:)
real(kind=8)::x(ndof,nnode),es(nee,nee),em(nee),eledet(nelement),&
	ef(nee),elemx(ndof,nen),elemu(nee),&
	mat(nelement,5),elemmat(5)!4:alpha,5:beta	
integer(kind=4)::elemvar(nee)
integer(kind=4),allocatable,dimension(:)::ja	
real(kind=8),allocatable,dimension(:)::kstiff
real(kind=8)::mass(neq),right(neq),right2(neq),&
	cuf(neq),cvf(neq),resu(neq),resv(neq),&
	u0(neq),u_1(neq),resa(neq)
real(kind=8),allocatable,dimension(:)::dump,left
real(kind=8)::f(neq),a0(neq),v0(neq),caf(neq),cRf(neq),cRfv(neq),cRfa(neq),cRfatrial(neq),a0trial(neq),Vol
real(kind=8),dimension(ndof,nnode)::constrain,constraintmp,constrainv,constrainvtmp,constraina,constrainf,constrainm
real(kind=8)::alfa,delt,aa0,aa1,aa2,aa3,aa4,aa5,aa6,aa7
!For faulting
integer(kind=4)::nftmx,nonmx,n4onf,nonfs
real(kind=8),dimension(2,nonmx,ntotft)::xonfs
integer(kind=4),dimension(ntotft)::nftnd
integer(kind=4),dimension(3,nonmx)::anonfs
integer(kind=4),dimension(2,nftmx,ntotft) :: nsmp
real(kind=8),dimension(nftmx,ntotft)::fnft,arn,r4nuc,slp4fri
real(kind=8),dimension(3,nftmx,ntotft)::un,us,ud,fltslp
real(kind=8),dimension(50,nftmx,ntotft)::fric
real(kind=8),dimension(10,nstep,nonmx)::fltsta
real(kind=8),dimension(10,nstep)::globaldat
real(kind=8)::dtev,dtev1,totmomrate,maxsliprate,totruptarea,tottaoruptarea, totslipruptarea,tdynastart,tdynaend,tcheck,tcheckstatus
integer(kind=8)::nev
!Parameters for pardiso
integer(kind=8)::pt(64)
integer(kind=4)::maxfct,mnum,mtype,phase,nrhs,error,msglvl,solver
integer(kind=4)::iparm(64),idum(1)
REAL*8  dparm(64)
real*8 ddum(1)
data nrhs/1/,maxfct/1/,mnum/1/
!--------------------------------------------------------------------c
integer(kind=4)::status0,status1
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

tmparr = 0.0d0
do i = 1,nftnd(1)
 tmparr(i) = fric(47,i,1)
enddo
	status0 = 0
	status1 = 0
        fnft = -1000.0d0
	fltsta = 0.0d0
        globaldat = 0.0d0
        totmomrate = 0.0d0
        maxsliprate = maxval(tmparr)
                write(*,*) 'Maxsliprate from last cycle = ',maxsliprate
        tdynastart = -1000.0d0
        tdynaend = -1000.0d0
        tcheck = 0.0d0
        tcheckstatus = -1.0d0
	if (me.eq.0) then
	!LOOP 2000 ![MUMPS]
	write(*,*) 'LOOP 2000: CONSTRUCT KSTIFF in CRS; CONVERT CRS to MUMPS FORMAT'
	!Define problem on the host (processor 0)
	dtev = ksi*LL/maxsliprate
	dtev1 = ksi*LL/maxsliprate!Initial dtev1>dt to enter static state.
	constrain = 0.0d0!Initialize constrain, displacements.
	constraintmp = 0.0d0
	constrainvtmp = 0.0d0
	constrainv = 0.0d0
	constraina = 0.0d0
	elemu=0.0d0

	allocate(lie(neq,100))
	do i=1,neq
		num(i)=0
		do j=1,100
			lie(i,j)=0
		enddo
	enddo

	!LOOP 1001: Get kloc(maxa),ja(maxa),ia(neq)
	write(*,*) 'Into subroutine CRSelem..'
	write(*,*) 'Starting LOOP 1001'

	do nel=1,nelement
		if (mod(nel,1000000)==0) then 
			write(*,*) 'Constructing overall stiffness'
			write(*,*) 'Progress..'
			write(*,*) int(nel/1000000),'M elems constructed'
		endif
	!--------------------------------------------------------------------c
		ntag=0
		do i=1,nen
			node_num=ien(i,nel)
			do j=1,ndof
				var=id(j,node_num)
				if (var.gt.0) then
					ntag=ntag+1
					elemvar(ntag)=var
				endif
			enddo
		enddo
	!--------------------------------------------------------------------c
		do i=1,ntag
			do j=i,ntag
				inv=elemvar(i)
				jnv=elemvar(j)
				if (jnv.lt.inv) then
					temp=jnv
					jnv=inv
					inv=temp
				endif

				if (num(inv).eq.0) then
					num(inv)=1
					lie(inv,1)=jnv
	!------------------------------------------------------------------
				elseif (num(inv).eq.1) then
					if (jnv.eq.lie(inv,1)) then
					elseif (jnv.lt.lie(inv,1)) then
						num(inv)=num(inv)+1
						lie(inv,2)=lie(inv,1)
						lie(inv,1)=jnv
					else
						num(inv)=num(inv)+1
						lie(inv,2)=jnv
					endif
	!------------------------------------------------------------------
				elseif (num(inv).gt.1) then
					if(jnv.lt.lie(inv,1)) then
						do mm=num(inv),1,-1
							lie(inv,mm+1)=lie(inv,mm)
						enddo
						lie(inv,1)=jnv
						num(inv)=num(inv)+1
					elseif (jnv.gt.lie(inv,num(inv))) then
						lie(inv,num(inv)+1)=jnv
						num(inv)=num(inv)+1
					else
					do m=1,num(inv)-1
						if ((jnv.eq.lie(inv,m)).or.  &
							(jnv.eq.lie(inv,m+1))) then
							elseif ((jnv.gt.lie(inv,m)).and.&
							(jnv.lt.lie(inv,m+1)))then
								do mm=num(inv),m+1,-1
									lie(inv,mm+1)=lie(inv,mm)
								enddo
								lie(inv,m+1)=jnv
								num(inv)=num(inv)+1
							endif
						enddo    
					endif        
				endif
	!-----------------------------------------------------------------
			enddo
		enddo
	enddo!nel loop 
	maxa=0
	do i=1,neq
		maxa=maxa+num(i)
	enddo
	write(*,*) 'maxa = ', maxa
	allocate(ja(maxa))
	mloc=0
	do i=1,neq
		do j=1,num(i)
			ja(mloc+j)=lie(i,j)
		enddo
		ia(i)=mloc+1
		mloc=mloc+num(i)
	enddo
	ia(neq+1)=mloc+1
	allocate(dump(maxa))
	allocate(kstiff(maxa))
	do i=1,maxa
		kstiff(i)=0.0d0
		dump(i)=0.0d0
	enddo
	do i=1,neq
		mass(i)=0.0d0
		f(i)=0.0d0
	enddo
	deallocate(lie)
	!LOOP 2001 : Calculating element stiff, mass, damp and constructing the overall stiffness
	write(*,*) 'LOOP 2001'
	write(*,*) 'Calculating element stiff, mass, damp and constructing the overall stiffness.'
	do  nel=1,nelement	
		if (mod(nel,1000000)==0) then 
			write(*,*) 'Constructing overall stiffness'
			write(*,*) 'Progress..'
			write(*,*) int(nel/1000000),'M elems constructed'
		endif
	!LOOP 130 ：Give the coordinate information of elements
		do j=1,nen
			node_num=ien(j,nel)
			do i=1,ndof
				elemx(i,j)=x(i,node_num)
			enddo
		enddo
		do j=1,5
			elemmat(j)=mat(nel,j)
		enddo
		! if (nel==5960775) then 
			! write(*,*) elemmat(1),elemmat(2),elemmat(3),elemmat(4),elemmat(5)
			! write(*,*) elemx(1,1),elemx(2,1),elemx(3,1)
			! write(*,*) elemx(1,2),elemx(2,2),elemx(3,2)
			! write(*,*) elemx(1,3),elemx(2,3),elemx(3,3)
			! write(*,*) elemx(1,4),elemx(2,4),elemx(3,4)
			! write(*,*) elemx(1,5),elemx(2,5),elemx(3,5)
			! write(*,*) elemx(1,6),elemx(2,6),elemx(3,6)
			! write(*,*) elemx(1,7),elemx(2,7),elemx(3,7)
			! write(*,*) elemx(1,8),elemx(2,8),elemx(3,8)
		! endif 
		
		call c8g2(elemx,elemmat,es,em,ef,Vol)
		eledet(nel)=Vol
		! if (nel==12) then 
			! write(*,*) es(1,1),es(1,2),es(1,3)
			! write(*,*) es(1,4),es(1,5),es(1,6)
			! write(*,*) es(1,7),es(1,8),es(1,9)
			! write(*,*) es(1,10),es(1,11),es(1,12)
			! write(*,*) es(1,13),es(1,14),es(1,15)
			! write(*,*) es(1,16),es(1,17),es(1,18)
			! write(*,*) es(1,19),es(1,20),es(1,21)
			! write(*,*) es(1,22),es(1,23),es(1,24)
			! write(*,*) '========em'
			! write(*,*) em(1),em(2),em(3)
			! write(*,*) em(4),em(5),em(6)
			! write(*,*) em(7),em(8),em(9)
			! write(*,*) em(10),em(11),em(12)
			! write(*,*) em(13),em(14),em(15)
			! write(*,*) em(16),em(17),em(18)
			! write(*,*) em(19),em(20),em(21)
			! write(*,*) em(22),em(23),em(24)
		! endif 	
		ntag=0
		do i=1,nen
			node_num=ien(i,nel)
			do j=1,ndof
			ntag=ntag+1
			elemvar(ntag)=id(j,node_num)
			!elemu(iii)=constrain(j,node_num)
			enddo
		enddo
	!   write(*,*)  (elemvar(i),i=1,iii)
	!------------------------------------------------------------------
		do i=1,nee
			inv=elemvar(i)
			if (inv.gt.0) then
				mass(inv)=mass(inv)+em(i)
				dump(ia(inv))=dump(ia(inv))+rdalfa*em(i)
				f(inv)=f(inv)+ef(i)
			endif
		enddo
		  
		do i=1,nee
			do j=i,nee
				inv=elemvar(i)          
				jnv=elemvar(j)
				if ((inv.gt.0).and.(jnv.gt.0)) then
					if (jnv.lt.inv) then
						temp=jnv
						jnv=inv
						inv=temp
					endif
	!...m为第inv行，对角线元向右非零元素在一维存储中的位置。
					do m=ia(inv),ia(inv)+num(inv)-1
						if (jnv.eq.ja(m)) then
							kstiff(m)=kstiff(m)+es(i,j)
							dump(m)=dump(m)+rdbeta*es(i,j)
						endif
					enddo
				endif
			enddo
		enddo

		do j=1,nee
			jnv=elemvar(j) 
			if (jnv.lt.0) then
				do i=1,nee                   
					inv=elemvar(i)
					if (inv.gt.0) then
						!f(inv)=f(inv)-es(i,j)*elemu(j)
					endif
				enddo
			endif
		enddo	
	enddo

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

	write(*,*) '--------------------------------------------------------'
	write(*,*) 'Begin the solver'
	write(*,*) 'Time integration scheme'
	if (IS==1) then 
		write(*,*) 'Central difference selected'
	elseif (IS==2) then 
		write(*,*) 'Newmark method selected & Solver Pardiso applied'
	endif
	if (nftnd(1) > 0) then !RSF
		do l=1,nftnd(1)
			constrainv(1,nsmp(1,l,1)) = fric(34,l,1)! -mass(id(1,nsmp(2,l,1)))/(mass(id(1,nsmp(1,l,1)))+mass(id(1,nsmp(2,l,1)))) * 1.0d-9 *us(1,l,1)!Slave
			constrainv(2,nsmp(1,l,1)) =  fric(35,l,1)!-mass(id(1,nsmp(2,l,1)))/(mass(id(1,nsmp(1,l,1)))+mass(id(1,nsmp(2,l,1)))) * 1.0d-9 *us(2,l,1)
			constrainv(3,nsmp(1,l,1)) =  fric(36,l,1)!-mass(id(1,nsmp(2,l,1)))/(mass(id(1,nsmp(1,l,1)))+mass(id(1,nsmp(2,l,1)))) * 1.0d-9 *us(3,l,1)
			constrainv(1,nsmp(2,l,1)) =  fric(31,l,1)!mass(id(1,nsmp(1,l,1)))/(mass(id(1,nsmp(1,l,1)))+mass(id(1,nsmp(2,l,1)))) * 1.0d-9 *us(1,l,1)!Master 
			constrainv(2,nsmp(2,l,1)) =  fric(32,l,1)!mass(id(1,nsmp(1,l,1)))/(mass(id(1,nsmp(1,l,1)))+mass(id(1,nsmp(2,l,1)))) * 1.0d-9 *us(2,l,1) 
			constrainv(3,nsmp(2,l,1)) =  fric(33,l,1)!mass(id(1,nsmp(1,l,1)))/(mass(id(1,nsmp(1,l,1)))+mass(id(1,nsmp(2,l,1)))) * 1.0d-9 *us(3,l,1) 	
		enddo
		
	endif

endif!END LOOP 2000 ![MUMPS]

if (me.eq.0) then

	allocate(left(maxa))
	left=0.0d0

	do i=1,neq
		u0(i)=0.0d0
		v0(i)=0.0d0
		a0(i)=0.0d0
		u_1(i)=0.0d0!u_1 is a dummy u0.
	enddo	  

!	NEWMARK Integration
!	MATRIX = [S]+[M]*a0+[C]*a1
!	FORC = [F]+[M*U1]*a0+[M*V1]*a2+[M]*[W1]*a3+[C*U1]*a1+[C]*[V1]*a4+[C]*[W1]*a5  
!	Form left matrix
!Analysis left matrix with pardiso   
    do i=1,maxa
        left(i)=kstiff(i)+aa1*dump(i)
    enddo
    do i=1,neq
        left(ia(i))=left(ia(i))+mass(i)*aa0
    enddo

    itag=0 ! label the step needed written in files


endif
!*****************************  III-v1  ********************************c
!---------------------     Parameters for MUMPS     --------------------c 
!	Call package for solution
mumps_par%JOB = 4! 1 + 2
CALL DMUMPS(mumps_par)
	! IF ( mumps par%INFOG ( 1 ) .LT. 0 ) THEN
	! WRITE( 6 , ’ (A, A, I6 , A, I 9 ) ’ ) ” ERROR RETURN: ” ,
	! & ” mumps par%INFOG ( 1 ) = ” , mumps par%INFOG ( 1 ) ,
	! & ” mumps par%INFOG ( 2 ) = ” , mumps par%INFOG ( 2 )
	! GOTO 500
	! END IF
	
    do it=1,nstep 
		if (me.eq.0) then
	!PART3===ENTER TIME LOOP===!
	!write(*,*) '!---3.1===CHECK CURRENT STATUS DYNAMIC/NON-DYNAMIC',me    
			!if (dtev1>dt) then!dtev1 was calculated in last timestep's faulting.
			!	status1=0
			!elseif (dtev1<=dt) then
		        !		dtev1 = dt
			!	status1=1
		!		exit
			!endif
	!write(*,*) '!---3.2===ADJUST LEFT_HAND MATRIX ACCORDING TO STATUS'   		
			!if ((status1-status0)==1) then
			!	write(*,*) 'CHANGE STATUS from STATIC to DYNAMIC'
			!elseif ((status1-status0)==-1) then
			!	write(*,*) 'CHANGE STATUS from DYNAMIC to STATIC'	
			!endif!
	!write(*,*) '!---3.3===UPDATE TIME & PRINT' 
			!time=(it)*dt
			time=time+dtev1
			write(*,*) it,'time step'
			write(*,*) time/60.0d0/60.0d0/24.0d0/365.0d0, 'year'
                        if (time/60/60/24/365>term) then
                                exit
                        endif
                        if (maxsliprate>1.0d-3)then 
                          if (status0 == 0) then        
                            status1 = 1
                          elseif (status0 == 1) then
                            if (tcheckstatus>0.0d0) then
                              tcheckstatus = -1.0d0
                              tdynaend = -1000.0d0
                            endif
                          endif
                        elseif (maxsliprate<=1.0d-3) then
                          if (status0 == 0) then
                            status1 = 0
                          elseif (status0 == 1) then
                            if (tcheckstatus<0.0d0) then
                                tcheckstatus = 1.0d0
                                tdynaend = time
                            else
                              if ((time-tdynaend)>=10.0d0.and.tdynaend>0.0d0) exit
                            endif
                          endif
                        endif
                        if ((status1-status0)==1) then
                                tdynastart = time
                        endif
	!write(*,*) '!---3.4===COSEISMIC DYNAMIC: STATUS == 1'
		endif!MYID==0	
		if (status1 == 1.or.status1 ==0) then		
		!elseif (status1==0) then
			
			if (me.eq.0) then
			
				right=0.0d0
				resu=0.0d0
				cRf=0.0d0
				constraintmp=0.0d0
	!write(*,*) '!---3.5.1===UPDATE BOUNDARY U*(t+1) INTO [CONSTRAINTMP]'	
	!--------===STORE EQUVILENT FORCES INTO [RIGHT]
				do nel=1,nelement
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
								if (elemvar(ntag)==-2) then
									constrainv(j,node_num)=-loadrate!-5.0d-10
									constrainvtmp(j,node_num) = constrainv(j,node_num)
									constraina(j,node_num)=0.0d0
								elseif (elemvar(ntag)==-3) then
									constrainv(j,node_num)=loadrate!5.0d-10
									constrainvtmp(j,node_num) = constrainv(j,node_num)
									constraina(j,node_num)=0.0d0					
								elseif (elemvar(ntag)==-1) then
									constrainv(j,node_num)=0.0d0
									constrainvtmp(j,node_num) = constrainv(j,node_num)
									constraina(j,node_num)=0.0d0
								endif 

								!First predictions of U*(t+1)[CONSTRAINTMP]
								constraintmp(j,node_num)=constrain(j,node_num)+constrainv(j,node_num)*dtev1
								elemu(ntag) = constraintmp(j,node_num)
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
					
				enddo!nel,3.5.1	
				!.. Computing Matrix(CRS)-vector multiplication      
				! do i=1,neq
					! do j=1,num(i)
						! ki=i
						! kj=ja(ia(i)+j-1)
						! kloc=ia(i)+j-1
						! if (ki.eq.kj) then
							! right(ki)=right(ki)-kstiff(kloc)*u0(ki)
						! else 
							! right(ki)=right(ki)-kstiff(kloc)*u0(kj)
							! right(kj)=right(kj)-kstiff(kloc)*u0(ki)							
						! endif
					! enddo
				! enddo 
	!write(*,*) '!---3.5.2===USE PARDISO TO COMPUTE U*(t+1) FOR THE WHOLE VOLUME'
	!--------===STORE ALL U*(t+1) INTO [CONSTRAINTMP]
				mumps_par%RHS = right
			endif!MYID==0
			mumps_par%JOB = 3
                        mumps_par%ICNTL(3) = -1 ! Goblal info, default 6 
			CALL DMUMPS( mumps_par )!MPI
			if (me.eq.0) then
				resu = mumps_par%RHS
				! call pardiso (pt,maxfct,mnum,mtype,phase,neq,kstiff,ia,ja,&
						! idum,nrhs,iparm,msglvl,right,resu,error,dparm) 
				do i=1,nnode
					do j=1,ndof
						if (id(j,i)>0) then !Only nodes that have equation number have been updated here.
							constraintmp(j,i)=resu(id(j,i))
						endif
					enddo 
				enddo 	
	!write(*,*) '!---3.5.3===COMPUTE F*=KU*(t+1)=[KSTIFF].DOT.[CONSTRAINTMP]'
	!--------===COMPUTE LUMPED MASS	FOR SPLIT NODES	
				constrainf=0.0d0
				constrainm=0.0d0!Updated only when itag == 0	
				do nel=1,nelement!Compute traction. 
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
								elemvar(ntag)=id(j,node_num)
								elemu(ntag)=constraintmp(j,node_num)
							enddo 
						enddo		
						do i=1,nen!-KU(t+1)
							node_num=ien(i,nel)
							do j=1,ndof 
								if (id(j,node_num)<0) then 
									do m=1,nee
										constrainf(j,node_num)=constrainf(j,node_num)-es((i-1)*3+j,m)*elemu(m)
									enddo
								constrainm(j,node_num)=constrainm(j,node_num)+em((i-1)*3+j)	
								endif
							enddo 
						enddo					
					endif
				enddo!nel
	!write(*,*) '!---3.5.4===USE FAULTING.F90 TO GET FIRST PREDICTIONS OF V*(t+1).'
	!--------===ITAG == 0, INDICATING THE FIRST TRY.
				!INPUT: - F* <-> [CONSTRAINF]
				!		- [CONSTRAINVTMP] TO STORE V*(t+1) AT THE END OF IT.
				itag = 0!
				!write(*,*) 'ITAG = 0'
				do ift=1,ntotft
					if (nftnd(ift)>0) then !only nonzero fault node, does faulting. B.D. 10/16/09
						call faulting(ift,nftnd(ift),nnode,neq,constrainm,right,constrain,constrainv,constraina,constrainf,id,&
							n4onf,fltsta,&
							nsmp(1,1,ift),fnft(1,ift),fltslp(1,1,ift),&
							un(1,1,ift),us(1,1,ift),ud(1,1,ift),fric(1,1,ift),arn(1,ift),r4nuc(1,ift),&
							slp4fri(1,ift),anonfs,nonmx,it,x,&
							kstiff,ia,maxa,right2,itag,dtev1,dtev,status1,constrainvtmp, &
                                                        maxsliprate,totmomrate,totruptarea,tottaoruptarea,totslipruptarea)
					endif
				enddo
	!write(*,*) '!---3.5.5===V*(t+1) OBTAINED. <-> REPEAT 3.5.1'
	!--------===UPDATE BOUNDARY U**(t+1),THEN DECLARE [CONSTRAIN] BY U**(t+1).	
	!--------===STORE EQUVILENT FORCES INTO [RIGHT]		
				right = 0.0d0
				constraintmp = 0.0d0
				do nel=1,nelement
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
							if (elemvar(ntag)<0)then
								!-2->-x; -3->+x; -5->fault
								if (elemvar(ntag)==-2) then
									constrainv(j,node_num)=-loadrate!-5.0d-10
									constraina(j,node_num)=0.0d0
								elseif (elemvar(ntag)==-3) then
									constrainv(j,node_num)=loadrate!5.0d-10
									constraina(j,node_num)=0.0d0					
								elseif (elemvar(ntag)==-1) then
									constrainv(j,node_num)=0.0d0
									constraina(j,node_num)=0.0d0
								endif 
	!write(*,*) '!---3.5.5.1-FINALL UPDATION OF [CONSTRAIN].'
	!-----------NO NEED TO WORRY NON-BOUNDARY NODES, SINCE THEIR CONSTRAINV&VTMP ARE ZERO.
	!-----------ERROR, CANNOT UPDATE CONTRAIN ITSELF. EACH D.O.F HAS BEEN VISITIED SEVERAL TIMES FOR IN THE LOOP.
								constraintmp(j,node_num) = constrain(j,node_num) + 0.5d0*(constrainv(j,node_num)+constrainvtmp(j,node_num))*dtev1
								elemu(ntag) = constraintmp(j,node_num)
							endif	
							enddo
						enddo			
						do j=1,nee
							jnv=elemvar(j) 
							if (jnv.lt.0) then
								do i=1,nee                   
									inv=elemvar(i)
									if (inv.gt.0) then
										right(inv) = right(inv) - es(i,j)*elemu(j)
									endif
								enddo
							endif
						enddo				
					endif
				enddo!nel
				! do i=1,neq
					! do j=1,num(i)
						! ki=i
						! kj=ja(ia(i)+j-1)
						! kloc=ia(i)+j-1
						! if (ki.eq.kj) then
							! right(ki)=right(ki)-kstiff(kloc)*resu(ki)
						! else 
							! right(ki)=right(ki)-kstiff(kloc)*resu(kj)
							! right(kj)=right(kj)-kstiff(kloc)*resu(ki)							
						! endif
					! enddo
				! enddo			
	!write(*,*) '!---3.5.6===USE PARDISO TO COMPUTE U**(t+1) FOR THE WHOLE VOLUME. <-> REPEAT 3.5.2'
	!--------===STORE ALL U**(t+1) INTO [CONSTRAIN]			
				mumps_par%RHS = right
			endif!MYID==0
			mumps_par%JOB = 3
			CALL DMUMPS( mumps_par )!MPI
			if (me.eq.0) then
				resu = mumps_par%RHS
				! call pardiso (pt,maxfct,mnum,mtype,phase,neq,kstiff,ia,ja,&
						! idum,nrhs,iparm,msglvl,right,resu,error,dparm) 
				open(1122,file='dispprofile.txt',position='append')
				open(1123,file='dispprofilecoor.txt')
				do i=1,nnode
					do j=1,ndof
						if (id(j,i)>0) then !Only nodes that have equation number have been updated here.
							constrain(j,i)=resu(id(j,i))
						else
							constrain(j,i)=constraintmp(j,i)
						endif
					enddo 
					! if (x(1,i)==0.0d0.and.x(2,i)>=0.0d0) then 
						! if (it==1) then 
							! write(1123,'(3e32.21e4)') x(1,i),x(2,i),x(3,i)
						! endif
						! write(1122,'(6e32.21e4)') constrain(1,i),constrain(2,i),constrain(3,i),right(id(1,i)),right(id(2,i)),right(id(3,i))
					! endif
				enddo
	!write(*,*) '!---3.5.7===COMPUTE F**=KU**(t+1)=[KSTIFF].DOT.[CONSTRAIN]'
	!--------===NO NEED TO COMPUTE LUMPED MASS	FOR SPLIT NODES			
				constrainf=0.0d0
				do nel=1,nelement!Compute traction. 
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
											ntag = 0
						call c8g2(elemx,elemmat,es,em,ef,Vol)		
						do i=1,nen!-KU**(t+1)
							node_num=ien(i,nel)
							do j=1,ndof
															ntag = ntag +1 
								elemu(ntag) = constrain(j,node_num)
							enddo 
						enddo	
						do i=1,nen!-KU(t+1)
							node_num=ien(i,nel)
							do j=1,ndof 
								if (id(j,node_num)<0) then 
									do m=1,nee
										constrainf(j,node_num)=constrainf(j,node_num)-es((i-1)*3+j,m)*elemu(m)
									enddo
								endif
							enddo 
						enddo					
					endif
				enddo!nel
	!write(*,*) '!---3.5.8===USE FAULTING.F90 TO GET SECOND PREDICTIONS OF V**(t+1), AND DECLARE V(t+1)=V**(t+1).'
	!--------===ITAG == 1, INDICATING THE SECOND AND FINAL TRY.
				!INPUT: - F* <-> [CONSTRAINF]
				!		- [CONSTRAINV] TO STORE V**(t+1) AT THE END OF IT.			
				itag = 1!
				!write(*,*) 'ITAG = 1'			
				do ift=1,ntotft
					if (nftnd(ift)>0) then !only nonzero fault node, does faulting. B.D. 10/16/09
						call faulting(ift,nftnd(ift),nnode,neq,constrainm,right,constrain,constrainv,constraina,constrainf,id,&
							n4onf,fltsta,&
							nsmp(1,1,ift),fnft(1,ift),fltslp(1,1,ift),&
							un(1,1,ift),us(1,1,ift),ud(1,1,ift),fric(1,1,ift),arn(1,ift),r4nuc(1,ift),&
							slp4fri(1,ift),anonfs,nonmx,it,x,&
							kstiff,ia,maxa,right2,itag,dtev1,dtev,status1,constrainvtmp,&
                                                        maxsliprate,totmomrate,totruptarea,tottaoruptarea,totslipruptarea)
					endif
				enddo
			endif!MYID==0				
		endif!status1 == 0/ ==1
		if (me.eq.0) then
			!PART4:FINISHING UP
			!------REFRESH STATUS,DTEV1		
			status0=status1
			nev=int8(dtev/dt)
			!if (status1==0) then
				dtev1=max(dt,nev*dt)	
			!else
			!	dtev1=dt
			!endif	
                        globaldat(1,it) = time
                        globaldat(2,it) = maxsliprate
                        globaldat(3,it) = totmomrate	
                        globaldat(4,it) = tottaoruptarea
						globaldat(5,it) = totslipruptarea
						globaldat(6,it) = totruptarea

			write(*,*) 'New time step dtev1 = ', dtev1, 'seconds'
                        write(*,*) 'Pma = ', pma
                        write(*,*) 'Sliprate max =', maxsliprate
			write(*,*) '******* Confirm Status 0',status1,'Confirmed*******'
		endif!MYID==0
		
	enddo!it
!endif!IS
!
if (me.eq.0) then
	write(*,*) 'Writing out results on fault stations'	
	!...on-fault stations. B.D. 8/11/10
	if(n4onf > 0) then
		do i=1,n4onf
			j=anonfs(3,i)
			if(j==1)  then  !main fault stations
				sttmp = '      '
				dptmp = '      '
				write(sttmp,'(i4.3)') int(xonfs(1,anonfs(2,i),j)/100.d0) 
				write(dptmp,'(i4.3)') int((-xonfs(2,anonfs(2,i),j)-60.0d3)/100.d0) 
				open(51,file='fltst_strk-'//trim(adjustl(sttmp))//'dp'//trim(adjustl(dptmp))//'.txt',status='unknown')

				! sttmp = '      '
				! dptmp = '      '
				! write(sttmp,'(f5.1)') xonfs(1,anonfs(2,i),j)/1000.d0 
				! write(dptmp,'(f5.1)') abs(xonfs(2,anonfs(2,i),j)/1000.d0) 
				! loca = '# location = on fault, '//trim(adjustl(sttmp))//' km along strike, '//trim(adjustl(dptmp))//' km down-dip'		
			endif
                        !write(51,*) '# This is the file header'
			!write(51,*) '# problem=SEAS Benchmark No.4'
			! write(51,*) '#Author=D.Liu'
		        !write(51,*) '# code = EQSimu'
			!write(51,*) '# version=1.0'
			!write(51,*) '# modeler = A.Modeler'
                        !write(51,*) '#element_size = 200 m x 200 m on fault'
			! write(51,*) '# Time series in 8 columns in format e15.7'
			! write(51,*) '# Column #1 = Time (s)'
			! write(51,*) '# Column #2 = horizontal slip (m)'
			! write(51,*) '# Column #3 = horizontal slip rate (m/s)'
			! write(51,*) '# Column #4 = horizontal shear stress (MPa)'
			! write(51,*) '# Column #3 = down-dip slip (m)'
			! write(51,*) '# Column #6 = down-dip slip rate (m/s)'
			! write(51,*) '# Column #7 = down-dip shear stress (MPa)'
			! !write(51,*) '# Column #8 = normal stress (MPa)'
			! write(51,*) '# The line below lists the names of the data fields:'
			! write(51,*) 't h-slip h-slip-rate h-shear-stress v-slip v-slip-rate v-shear-stress'
			! write(51,*) '#'		
			do j=1,it-1
			write(51,'( 8e32.21e4)') fltsta(1,j,i),fltsta(5,j,i),-fltsta(6,j,i), fltsta(2,j,i), &
                                fltsta(3,j,i),fltsta(8,j,i)/1.d6,-fltsta(9,j,i)/1.d6,dlog10(fltsta(4,j,i))
			enddo
			close(51)
		enddo
	endif
	write(*,*) 'Writing out cplot_EQquasi.txt'	
	if(nftnd(1) > 0) then
		open(unit=1111,file='cplot_EQquasi.txt',status='unknown')	!rupture time
		! write(1111,*) '#This is the file header:'
		! write(1111,*) '# problem=TPV8'
		! write(1111,*) '# author=D.Liu'
		! write(1111,*) '# date=2016/03/10'
		! write(1111,*) '# code=EQquasi3d'
		! write(1111,*) '# code_version=1.0'
		! write(1111,*) '# element_size=100 m'
		! write(1111,*) '# Column #1 = horizontal coordinate, distance along strike (m)'
		! write(1111,*) '# Column #2 = vertical coordinate, distance down-dip (m)'
		! write(1111,*) '# Column #3 = rupture time (s)'
		! write(1111,*) '# The line below lists the names of the data fields.'
		! write(1111,*) '# It indicates that the first column contains the horizontal'
		! write(1111,*) '# coordinate (j), the second column contains the vertical'
		! write(1111,*) '# coordinate (k), and the third column contains the time (t).'
		! write(1111,*) 'j k t'
		! write(1111,*) '# Here is the rupture history'	
		!write(1111,'(1x,7e32.21e4)') ((x(j,nsmp(1,i,1)),j=1,3),fnft(i,1),&
		!	(fltslp(j,i,1),j=1,3),i=1,nftnd(1))
		!!-FORMAT
		!!-1,    2,     3, fric(26,i),4,fric(20,i),5,fric(28,i),6,fric(29,i),7,fric(30,i)   
		!!-XCOOR,+ZCOOR, V_trial,     psi,        ,tstk,       ,tdip,       ,tnrm0
		write(1111,'(1x,16e32.21e4)') (x(1,nsmp(1,i,1)),-x(3,nsmp(1,i,1)),fric(26,i,1),fric(20,i,1),fric(28,i,1), &
				fric(29,i,1),fric(30,i,1),fric(31,i,1),fric(32,i,1),fric(33,i,1),fric(34,i,1),fric(35,i,1),fric(36,i,1), &
                                fric(44,i,1), fric(45,i,1), fnft(i,1), i=1,nftnd(1))	
		close(1111)
	endif
        
        write(*,*) 'Writing out tdynastart and tdynaend'
        open(1112,file='tdyna.txt',form = 'formatted', status = 'unknown')
                write(1112,'(2e32.21e4)') tdynastart, tdynaend
        close(1112)

        write(*,*) 'Writing out global.dat'
        open(1113,file='global.dat',form='formatted',status='unknown')
                write(1113,'(6e32.21e4)') (globaldat(1,i),globaldat(2,i),globaldat(3,i),globaldat(4,i),globaldat(5,i),globaldat(6,i),i=1,it-1)
        close(1113)
endif!MYID==0
! Deallocate user data
IF ( me.eq.0 )THEN
DEALLOCATE( mumps_par%IRN )
DEALLOCATE( mumps_par%JCN )
DEALLOCATE( mumps_par%A )
DEALLOCATE( mumps_par%RHS )
END IF
! Destroy the instance (deallocate internal data structures)
! mumps_par%JOB = -2
! CALL DMUMPS( mumps_par )	
	
end subroutine CRSelem
