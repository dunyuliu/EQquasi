subroutine meshgen

	use globalvar
	implicit none
	include 'mpif.h'	
		
	integer(kind=4)::nnode,nelement,neq0,ix,iy,iz,&
		edgex1,edgey1,edgez1,edgezn,&
		i,j,k,i1,j1,ift
	real(kind=dp)::dy,dz
	real(kind=dp)::tmp1,tmp2
	integer (kind=4)::nxuni,nyuni,nzuni	
	real(kind=dp)::tol,xcoor,ycoor,zcoor,xstep,ystep,zstep,&
			xmin1,xmax1,ymin1,ymax1,zmin1
	real (kind = dp), allocatable, dimension(:) :: xlinet,ylinet,zlinet,xline,yline,zline
	integer (kind = 4), allocatable, dimension(:,:) :: plane1,plane2
	integer (kind = 4), allocatable, dimension(:,:,:,:) :: fltrc	

	integer(kind=4)::n1,n2,n3,n4,m1,m2,m3,m4
	real(kind=dp)::a1,b1,c1,d1,p1,q1,area 
	integer (kind=4),dimension(ntotft) :: nftnd0,ixfi,izfi,ifs,ifd
	real (kind=dp) :: x1,x2,x3=0.d0,y1,y2,y3=0.d0,z1,z2,z3=0.d0 !3 points for fault plane (inc origin)

	real(kind=dp)::center(3),ycoort,strtmp,kinkx, pfx, pfz
	integer(kind=4)::kkk,nfx,nfz,ntmp 	
	real(kind=dp),allocatable::initial(:,:)
	
	allocate(n4yn(n4nds))
	n4yn = 0
	
	ntmp = (int((xmax-xmin)/dx)+1)*(int((zmax-zmin)/dx)+1)										
	allocate(initial(16,ntmp))
	initial = 0.0d0 
	
	if (icstart > 1) then 
		open(111,file = 'cplot_EQquasi.txt', form = 'formatted', status ='old')
		  do i = 1,ntmp 
			read(111,*) (initial(j,i),j=1,16)
		  enddo
		close(111)
	endif 
	dy=dx
	dz=dx
	tol=dx/100.d0
		
	nxuni=(fltxyz(2,1,1)-fltxyz(1,1,1)-2.0d0*dx)/dx+1
	xstep=dx
	xcoor=fltxyz(1,1,1)+dx
	do ix=1,np
		xstep=xstep
		xcoor=xcoor-xstep
		if(xcoor<=xmin) exit
	enddo
	edgex1=ix
	xstep=dx
	xcoor=fltxyz(2,1,1)-dx
	do ix=1,np
		xstep=xstep
		xcoor=xcoor+xstep
		if(xcoor>=xmax) exit
	enddo
	nxt=nxuni+edgex1+ix
	allocate(xlinet(nxt))
	!predetermine x-coor
	xlinet(edgex1+1)=fltxyz(1,1,1)+dx
	xstep=dx
	do ix=edgex1,1,-1
		xstep=xstep
		xlinet(ix)=xlinet(ix+1)-xstep
	enddo
	xmin1=xlinet(1)
	do ix=edgex1+2,edgex1+nxuni
		xlinet(ix)=xlinet(ix-1)+dx
	enddo
	xstep=dx
	do ix=edgex1+nxuni+1,nxt
		xstep=xstep
		xlinet(ix)=xlinet(ix-1)+xstep
	enddo
	xmax1=xlinet(nxt)
	!Y
	nyuni=dis4uniF+dis4uniB+1
	ystep=dy
	ycoor=-dy*(dis4uniF)
	do iy=1,np
		ystep=ystep*rat
		if (ystep>=dymax) ystep = dymax
		ycoor=ycoor-ystep
		if(ycoor<=ymin) exit
	enddo
	edgey1=iy
	ystep=dy
	ycoor=dy*(dis4uniB)
	do iy=1,np
		ystep=ystep*rat
		if (ystep>=dymax) ystep = dymax
		ycoor=ycoor+ystep
		if(ycoor>=ymax) exit
	enddo
	nyt=nyuni+edgey1+iy
	!...pre-determine y-coor
	allocate(ylinet(nyt))
	!...predetermine y-coor
	ylinet(edgey1+1)=-dy*(dis4uniF)
	ystep=dy
	do iy=edgey1,1,-1
		ystep=ystep*rat
		if (ystep>=dymax) ystep = dymax
		ylinet(iy)=ylinet(iy+1)-ystep
	enddo
	ymin1=ylinet(1)
	do iy=edgey1+2,edgey1+nyuni
		ylinet(iy)=ylinet(iy-1)+dy
	enddo
	ystep=dy
	do iy=edgey1+nyuni+1,nyt
		ystep=ystep*rat
		if (ystep>=dymax) ystep = dymax
		ylinet(iy)=ylinet(iy-1)+ystep
	enddo
	ymax1=ylinet(nyt)
	!Z
	zstep=dz
	zcoor=fltxyz(1,3,1)+dx
	do iz=1,np
		zstep=zstep
		zcoor=zcoor-zstep
		if(zcoor<=zmin) exit
	enddo
	edgezn=iz
	nzuni=(fltxyz(2,3,1)-fltxyz(1,3,1)-dx)/dx+1 
	nzt=edgezn+nzuni
	!...predetermine z-coor
	allocate(zlinet(nzt))
	zlinet(nzt)=zmax
	do iz=nzt-1,nzt-nzuni+1,-1
		zlinet(iz)=zlinet(iz+1)-dz
	enddo
	zstep=dz
	do iz=nzt-nzuni,1,-1
		zstep=zstep
		zlinet(iz)=zlinet(iz+1)-zstep
	enddo
	zmin1=zlinet(1)
	!...prepare for digitizing	
	allocate(plane1(nyt+ntotft,nzt),plane2(nyt+ntotft,nzt),fltrc(2,nxt,nzt,ntotft))
	nnode = 0
	nelement = 0
	neq0 = 0
	nftnd0 = 0 
	n4onf=0
	ixfi = 0
	izfi = 0
	!write(*,*) 'C1'
	!Strategy to assign boundary conditions. The fault is treated as boundaries.
	!1. Create a node. If the node is on the boundaries, do not allocate equation number. 
	!2. Loop over fault nodes. Split the node, b
	!write(*,*) xlinet(1),xlinet(2),xlinet(3)
	do ix = 1, nxt
		do iz = 1, nzt
			do iy = 1, nyt
				xcoor = xlinet(ix)
				ycoor = ylinet(iy)
				zcoor = zlinet(iz)
				nnode = nnode + 1
				x(1,nnode) = xcoor
				x(2,nnode) = ycoor
				x(3,nnode) = zcoor			
				strtmp=270.0d0/180.0d0*pi
				if (rough_fault == 1) then 
					call insert_rough_fault(xcoor, ycoor, zcoor, ycoort, pfx, pfz, ymax1, ymin1)
					x(2,nnode) = ycoort
				endif 
							!The normal vector at x = 0 belongs to the left segment.
				!if (xcoor>0.0d0) then 
	!
	!				tmp1 = (xcoor--0.0d3)*dtand(ttheta)
	!				tmp2=ttheta
	!				strtmp=(270.0d0-tmp2)/180.0d0*pi		
	!				if (ycoor>=0.0d3) then 	
	!					ycoort=(ylinet(nyt)-tmp1)*ycoor/ylinet(nyt) + tmp1
	!				elseif (ycoor<0.0d0) then 
	!					ycoort=ylinet(1)+(ycoor-ylinet(1))/(-ylinet(1))*(tmp1-ylinet(1))
	!				endif
	!				x(2,nnode)=ycoort							
	!			endif
				plane2(iy,iz) = nnode
				if(iy==1.or.iy==nyt.or.ycoor==0.0d0) then 
				!if(ycoor==0.0d0)then
					if (ycoor>0.0d0) then 
						id(1,nnode)=-3!+x movement
						id(2,nnode)=-1!fixed
						id(3,nnode)=-1!fixed
					elseif (ycoor<0.0d0) then 
						id(1,nnode)=-2!-x movement
						id(2,nnode)=-1!fixed
						id(3,nnode)=-1!fixed
					elseif (ycoor==0.0d0) then 
						id(1,nnode)=-5!fault boundaries
						id(2,nnode)=-5!
						id(3,nnode)=-5!
					endif 
				else 
					do i1=1,ndof
						neq0=neq0+1
						id(i1,nnode)=neq0!Only nodes inside two blocks have equation numbers.
					enddo
				endif
				if(ix>1.and.ix<nxt .and. iy>1.and.iy<nyt) then 
					do i=1,n4nds
						if(n4yn(i)==0) then
							if (abs(zcoor-x4nds(3,i))<tol) then
								if(abs(xcoor-x4nds(1,i))<tol .or.&
								(x4nds(1,i)>xlinet(ix-1).and.x4nds(1,i)<xcoor.and. &
								(xcoor-x4nds(1,i))<(x4nds(1,i)-xlinet(ix-1))) .or. &
								(x4nds(1,i)>xcoor.and.x4nds(1,i)<xlinet(ix+1).and. &
								(x4nds(1,i)-xcoor)<(xlinet(ix+1)-x4nds(1,i)))) then
									if(abs(ycoor-x4nds(2,i))<tol .or. &
									(x4nds(2,i)>ylinet(iy-1).and.x4nds(2,i)<ycoor.and. &
									(ycoor-x4nds(2,i))<(x4nds(2,i)-ylinet(iy-1))) .or. &
									(x4nds(2,i)>ycoor.and.x4nds(2,i)<ylinet(iy+1).and. &
									(x4nds(2,i)-ycoor)<(ylinet(iy+1)-x4nds(2,i)))) then
										n4yn(i) = 1
										n4out = n4out + 1
										an4nds(1,n4out) = i
										an4nds(2,n4out) = nnode
										exit 	!if node found, jump out the loop
									endif
								endif
							endif
						endif
					enddo
				endif				
				do ift=1,ntotft 
					if(xcoor>=(fltxyz(1,1,ift)-tol).and.xcoor<=(fltxyz(2,1,ift)+tol).and. &
					ycoor>=(fltxyz(1,2,ift)-tol).and.ycoor<=(fltxyz(2,2,ift)+tol).and. &
					zcoor>=(fltxyz(1,3,ift)-tol) .and. zcoor<=(fltxyz(2,3,ift)+tol)) then
						do i1=1,ndof
							id(i1,nnode)=-5
						enddo
						! if (xcoor==-18.0d3.and.zcoor==-21.0d3) then 
							! write(*,*) 'SLAVE:-18,-21,nnode',nnode,nftnd0(ift),ifs(ift),ifd(ift) 
							! write(*,*) '-18,-21,id',id(1,nnode),id(2,nnode),id(3,nnode)
						! endif					
						nftnd0(ift) = nftnd0(ift) + 1

						nsmp(1,nftnd0(ift),ift) = nnode !slave node						
						nnode = nnode + 1

						nsmp(2,nftnd0(ift),ift) = nnode !master node, on the +y side.
						plane2(nyt+ift,iz) = nnode
						x(1,nnode) = xcoor
						x(2,nnode) = ycoor
						x(3,nnode) = zcoor  
						if (rough_fault == 1) then 
							x(2,nnode) = ycoort
						endif
						!...identify output on-fault stations. B.D. 10/25/09
						!...fault-dependent for multiple faults. B.D. 1/8/
						do i1=1,nonfs(ift)
							if(abs(xcoor-xonfs(1,i1,ift))<tol .and. &
							abs(zcoor-xonfs(2,i1,ift))<tol) then
								n4onf = n4onf + 1
								anonfs(1,n4onf) = nftnd0(ift)
								anonfs(2,n4onf) = i1
								anonfs(3,n4onf) = ift
								exit         
							endif
						enddo  
						!...establish equation numbers for this master node
						do i1=1,ndof
							!neq0 = neq0 + 1
							id(i1,nnode)=-5!assume fault boundary nodes
						enddo
						!...assign unit vectors to the split node pair						
						un(1,nftnd0(ift),ift) = dcos(strtmp)*dsin(fltxyz(2,4,ift))
						un(2,nftnd0(ift),ift) = -dsin(strtmp)*dsin(fltxyz(2,4,ift))
						un(3,nftnd0(ift),ift) = dcos(fltxyz(2,4,ift))		
						us(1,nftnd0(ift),ift) = -dsin(strtmp)
						us(2,nftnd0(ift),ift) = -dcos(strtmp)
						us(3,nftnd0(ift),ift) = 0.0d0
						ud(1,nftnd0(ift),ift) = dcos(strtmp)*dcos(fltxyz(2,4,ift))
						ud(2,nftnd0(ift),ift) = dsin(strtmp)*dcos(fltxyz(2,4,ift))
						ud(3,nftnd0(ift),ift) = dsin(fltxyz(2,4,ift))
						if (rough_fault == 1) then
							un(1,nftnd0(ift),ift) = -pfx/(pfx**2 + 1.0d0 + pfz**2)**0.5
							un(2,nftnd0(ift),ift) = 1.0d0/(pfx**2 + 1.0d0 + pfz**2)**0.5
							un(3,nftnd0(ift),ift) = -pfz/(pfx**2 + 1.0d0 + pfz**2)**0.5
							us(1,nftnd0(ift),ift) = 1.0d0/(1.0d0 + pfx**2)**0.5
							us(2,nftnd0(ift),ift) = pfx/(1.0d0 + pfx**2)**0.5
							us(3,nftnd0(ift),ift) = 0.0d0
							ud(1,nftnd0(ift),ift) = us(2,nftnd0(ift),ift)*un(3,nftnd0(ift),ift) - us(3,nftnd0(ift),ift)*un(2,nftnd0(ift),ift)
							ud(2,nftnd0(ift),ift) = us(3,nftnd0(ift),ift)*un(1,nftnd0(ift),ift) - us(1,nftnd0(ift),ift)*un(3,nftnd0(ift),ift)
							ud(3,nftnd0(ift),ift) = us(1,nftnd0(ift),ift)*un(2,nftnd0(ift),ift) - us(2,nftnd0(ift),ift)*un(1,nftnd0(ift),ift)
						endif
						!...prepare for area calculation
						if(ixfi(ift)==0) ixfi(ift)=ix
						if(izfi(ift)==0) izfi(ift)=iz
						ifs(ift)=ix-ixfi(ift)+1
						ifd(ift)=iz-izfi(ift)+1
						fltrc(1,ifs(ift),ifd(ift),ift) = nnode	!master node
						fltrc(2,ifs(ift),ifd(ift),ift) = nftnd0(ift) !fault node num in sequence
						! if (xcoor==-18.0d3.and.zcoor==-21.0d3) then 
							! write(*,*) 'Master:-18,-21,nnode',nnode,nftnd0(ift),ifs(ift),ifd(ift) 
							! write(*,*) '-18,-21,id',id(1,nnode),id(2,nnode),id(3,nnode)
						! endif
						!...assign friction parameters for SCEC TPV19. B.D. 1/8/12
						!...now for Ma and Andrews (2010) model. B.D. 6/1/12
					!For TPV8	
						fric(1,nftnd0(ift),ift) = 0.18d0!mus
						fric(2,nftnd0(ift),ift) = 0.12d0!mud
						! if (xcoor<-20e3.and.xcoor>-40e3.and.zcoor<-3e3.and.zcoor>-8e3)then
							! fric(2,nftnd0(ift),ift) = 0.5   !mud
						! endif						
						fric(3,nftnd0(ift),ift) = 0.30d0    !D0
						fric(4,nftnd0(ift),ift) = 0.0d6  !cohesion
						fric(5,nftnd0(ift),ift) = 0.08d0  !Not used in TPV8 !T for forced rupture,i.e.,nucleation
						fric(6,nftnd0(ift),ift) = 0.0d0!rhow*grav*(-zcoor) !pore pressure
						! -------------------------------------------------------------------------
						! Case specific setup
						! if (friclaw == 3 .and. bp == 5) then
							! ! Setting a, b, Dc, v0, r0
							! if(abs(zcoor)>=18.0d3.or.abs(xcoor)>=32.0d3.or.abs(zcoor)<=2.0d3) then
								! fric(9,nftnd0(ift),ift) = fric_rsf_a + fric_rsf_deltaa0
							! elseif (abs(zcoor)<=16.0d3.and.abs(zcoor)>=4.0d3.and.abs(xcoor)<=30.0d3) then
								! fric(9,nftnd0(ift),ift) = fric_rsf_a
							! else
								! tmp1 = abs(abs(zcoor)-2.0d3-2.0d3-12.0d3/2.0d0) - 12.0d3/2.0d0
								! tmp1 = tmp1 / 2.0d3
								! tmp2 = (abs(xcoor)-30.0d3)/2.0d3		
								! fric(9,nftnd0(ift),ift) = fric_rsf_a + max(tmp1,tmp2)*fric_rsf_deltaa0					 
							! endif	
							! fric(10,nftnd0(ift),ift) = fric_rsf_b                                            
							! fric(11,nftnd0(ift),ift) = fric_rsf_Dc !RSF critical distance.
							! if (xcoor<-30.0d3+12.0d3+tol.and.xcoor>-30.0d3-tol.and.zcoor<-4.0d3+tol &
								 ! .and.zcoor>-16.0d3-tol) then
								! fric(11,nftnd0(ift),ift) = minDc
							! endif 
							! fric(12,nftnd0(ift),ift) = fric_rsf_v0 !RSF:V0
							! fric(13,nftnd0(ift),ift) = fric_rsf_r0 !RSF:miu0

							! fric(16,nftnd0(ift),ift) = 0.0d0 !RSF: initial normal slip rate 
							! fric(17,nftnd0(ift),ift) = fric_rsf_vinix !RSF:s 
							! fric(18,nftnd0(ift),ift) = fric_rsf_viniz !RSF:d
							! fric(19,nftnd0(ift),ift) = sqrt(fric_rsf_vinix**2+fric_rsf_viniz**2) !RSF:mag		
							! ! Setting theta
							! if (icstart == 1) then 
								! fric(7,nftnd0(ift),ift) = init_norm !-25.0d6 ! Initial normal stress 
								! fric(20,nftnd0(ift),ift )= fric(11,nftnd0(ift),ift)/load_slip_rate ! theta0 for steady state at v==load_slip_rate.
								! ! Setting initial shear stress. 
								! !fric(8,nftnd0(ift),ift) = - fric(7,nftnd0(ift),ift) * fric(9,nftnd0(ift),ift) *dasinh( &
								! !	1.0d-9/2.0d-6*dexp((0.6d0+0.03d0*dlog(1.0d-6/1.0d-9))/fric(9,nftnd0(ift),ift))) + 2670.0d0*3464.0d0/2.0d0*1.0d-9
								! fric(45, nftnd0(ift), ift) = 0.0d0
								! fric(46, nftnd0(ift), ift) = load_slip_rate ! initial slip rate.
								! call rsf_rd(fric(8,nftnd0(ift),ift), fric(7,nftnd0(ift),ift), fric(9,nftnd0(ift),ift), fric(10,nftnd0(ift),ift), fric(13,nftnd0(ift),ift), fric(12,nftnd0(ift),ift), mat0(1,2), mat0(1,3), fric(46, nftnd0(ift), ift))
								! fric(44, nftnd0(ift), ift) = fric(8, nftnd0(ift), ift) !tstk0
								! if (xcoor<-30.0d3+12.0d3+tol.and.xcoor>-30.0d3-tol.and.zcoor<-4.0d3+tol &
									 ! .and.zcoor>-16.0d3-tol) then
									! fric(46,nftnd0(ift),ift)=3.0d-2 ! In the special patch, increase initial slip rate. 
									! fric(8,nftnd0(ift),ift) = - fric(7,nftnd0(ift),ift) * fric(9,nftnd0(ift),ift) *dasinh( &
										! fric(46,nftnd0(ift),ift)/2.0d-6*dexp((0.6d0+0.03d0*dlog(1.0d-6/1.0d-9))/fric(9,nftnd0(ift),ift))) + 2670.0d0*3464.0d0/2.0d0*fric(46,nftnd0(ift),ift)
									! !call rsf_rd(fric(8,nftnd0(ift),ift), fric(7,nftnd0(ift),ift), fric(9,nftnd0(ift),ift), fric(10,nftnd0(ift),ift), fric(13,nftnd0(ift),ift), fric(12,nftnd0(ift),ift), mat0(1,2), mat0(1,3), fric(46, nftnd0(ift), ift))
									! fric(44, nftnd0(ift), ift) = fric(8, nftnd0(ift), ift)
									! fric(45, nftnd0(ift), ift) = 0.0d0
								! endif
								! fric(47,nftnd0(ift),ift) = fric(46,nftnd0(ift),ift)! Peak slip rate
							! elseif (icstart > 1) then
								! ! If icstart > 1, and use eqquasi to simulate 2nd and later cycles, adopts initials from previous models.
								! ! This is the elastic version that only transfers quantities on the fault. 
								! nfx = (xcoor - xmin)/dx+1
								! nfz = (zcoor - zmin)/dx+1
								! ntmp = (zmax-zmin)/dx + 1
								! fric(20,nftnd0(ift),ift) = initial(4,(nfx-1)*ntmp +nfz) !theta0
								! fric(7,nftnd0(ift),ift) = initial(7,(nfx-1)*ntmp +nfz) !tnrm0
								! fric(8,nftnd0(ift),ift) = initial(5,(nfx-1)*ntmp +nfz) !tstk
								! fric(49,nftnd0(ift),ift) = initial(6,(nfx-1)*ntmp +nfz) !tdip
								! fric(31,nftnd0(ift),ift) = initial(8,(nfx-1)*ntmp +nfz) ! vx0+
								! fric(32,nftnd0(ift),ift) = initial(9,(nfx-1)*ntmp +nfz) ! vy0+
								! fric(33,nftnd0(ift),ift) = initial(10,(nfx-1)*ntmp +nfz) ! vz0+
								! fric(34,nftnd0(ift),ift) = initial(11,(nfx-1)*ntmp +nfz) ! vx0-
								! fric(35,nftnd0(ift),ift) = initial(12,(nfx-1)*ntmp +nfz)
								! fric(36,nftnd0(ift),ift) = initial(13,(nfx-1)*ntmp +nfz)
								! fric(47,nftnd0(ift),ift) = initial(3,(nfx-1)*ntmp +nfz)! Peak slip rate
								! fric(44, nftnd0(ift), ift) = fric(8, nftnd0(ift), ift) !tstk0
								! fric(45, nftnd0(ift), ift) = 0.0d0
							 ! endif
						! endif !end case bp==5.
						if (friclaw == 3 .and. bp == 1001) then
							! Setting a, b, Dc, v0, r0
							if(abs(zcoor)>=18.0d3.or.abs(xcoor)>=20.0d3.or.abs(zcoor)<=2.0d3) then
								fric(9,nftnd0(ift),ift) = fric_rsf_a + fric_rsf_deltaa0
							elseif (abs(zcoor)<=16.0d3.and.abs(zcoor)>=4.0d3.and.abs(xcoor)<=18.0d3) then
								fric(9,nftnd0(ift),ift) = fric_rsf_a
							else
								tmp1 = abs(abs(zcoor)-2.0d3-2.0d3-12.0d3/2.0d0) - 12.0d3/2.0d0
								tmp1 = tmp1 / 2.0d3
								tmp2 = (abs(xcoor)-18.0d3)/2.0d3		
								fric(9,nftnd0(ift),ift) = fric_rsf_a + max(tmp1,tmp2)*fric_rsf_deltaa0					 
							endif	
							fric(10,nftnd0(ift),ift) = fric_rsf_b                                            
							fric(11,nftnd0(ift),ift) = fric_rsf_Dc !RSF critical distance.
							fric(12,nftnd0(ift),ift) = fric_rsf_v0 !RSF:V0
							fric(13,nftnd0(ift),ift) = fric_rsf_r0 !RSF:miu0
							fric(16,nftnd0(ift),ift) = 0.0d0 !RSF: initial normal slip rate 
							fric(17,nftnd0(ift),ift) = fric_rsf_vinix !RSF:s 
							fric(18,nftnd0(ift),ift) = fric_rsf_viniz !RSF:d
							fric(19,nftnd0(ift),ift) = sqrt(fric_rsf_vinix**2+fric_rsf_viniz**2) !RSF:mag		
							! Setting theta
							if (icstart == 1) then 
								fric(7,nftnd0(ift),ift) = init_norm !-25.0d6 ! Initial normal stress 
								fric(20,nftnd0(ift),ift )= fric(11,nftnd0(ift),ift)/load_slip_rate ! theta0 for steady state at v==load_slip_rate.
								! Setting initial shear stress. 
								!fric(8,nftnd0(ift),ift) = - fric(7,nftnd0(ift),ift) * fric(9,nftnd0(ift),ift) *dasinh( &
								!	1.0d-9/2.0d-6*dexp((0.6d0+0.03d0*dlog(1.0d-6/1.0d-9))/fric(9,nftnd0(ift),ift))) + 2670.0d0*3464.0d0/2.0d0*1.0d-9
								fric(45, nftnd0(ift), ift) = 0.0d0
								fric(46, nftnd0(ift), ift) = load_slip_rate ! initial slip rate.
								call rsf_rd(fric(8,nftnd0(ift),ift), fric(7,nftnd0(ift),ift), fric(9,nftnd0(ift),ift), fric(10,nftnd0(ift),ift), fric(13,nftnd0(ift),ift), fric(12,nftnd0(ift),ift), mat0(1,2), mat0(1,3), fric(46, nftnd0(ift), ift))
								fric(44, nftnd0(ift), ift) = fric(8, nftnd0(ift), ift) !tstk0
								
								! Specify initial high slip rate patch. 
								! if (xcoor<-8.0d3+tol.and.xcoor>-20.0d3-tol.and.zcoor<-4.0d3+tol &
									 ! .and.zcoor>-16.0d3-tol) then
									! fric(46,nftnd0(ift),ift)=3.0d-2 ! In the special patch, increase initial slip rate. 
									! fric(8,nftnd0(ift),ift) = - fric(7,nftnd0(ift),ift) * fric(9,nftnd0(ift),ift) *dasinh( &
										! fric(46,nftnd0(ift),ift)/2.0d-6*dexp((0.6d0+0.03d0*dlog(1.0d-6/1.0d-9))/fric(9,nftnd0(ift),ift))) + mat0(1,2)*mat0(1,3)/2.0d0*fric(46,nftnd0(ift),ift)
									! !call rsf_rd(fric(8,nftnd0(ift),ift), fric(7,nftnd0(ift),ift), fric(9,nftnd0(ift),ift), fric(10,nftnd0(ift),ift), fric(13,nftnd0(ift),ift), fric(12,nftnd0(ift),ift), mat0(1,2), mat0(1,3), fric(46, nftnd0(ift), ift))
									! fric(44, nftnd0(ift), ift) = fric(8, nftnd0(ift), ift)
									! fric(45, nftnd0(ift), ift) = 0.0d0
								! endif
								fric(47,nftnd0(ift),ift) = fric(46,nftnd0(ift),ift)! Peak slip rate
							elseif (icstart > 1) then
								! If icstart > 1, and use eqquasi to simulate 2nd and later cycles, adopts initials from previous models.
								! This is the elastic version that only transfers quantities on the fault. 
								nfx = (xcoor - xmin)/dx+1
								nfz = (zcoor - zmin)/dx+1
								ntmp = (zmax-zmin)/dx + 1
								fric(20,nftnd0(ift),ift) = initial(4,(nfx-1)*ntmp +nfz) !theta0
								fric(7,nftnd0(ift),ift) = initial(7,(nfx-1)*ntmp +nfz) !tnrm0
								fric(8,nftnd0(ift),ift) = initial(5,(nfx-1)*ntmp +nfz) !tstk
								fric(49,nftnd0(ift),ift) = initial(6,(nfx-1)*ntmp +nfz) !tdip
								fric(31,nftnd0(ift),ift) = initial(8,(nfx-1)*ntmp +nfz) ! vx0+
								fric(32,nftnd0(ift),ift) = initial(9,(nfx-1)*ntmp +nfz) ! vy0+
								fric(33,nftnd0(ift),ift) = initial(10,(nfx-1)*ntmp +nfz) ! vz0+
								fric(34,nftnd0(ift),ift) = initial(11,(nfx-1)*ntmp +nfz) ! vx0-
								fric(35,nftnd0(ift),ift) = initial(12,(nfx-1)*ntmp +nfz)
								fric(36,nftnd0(ift),ift) = initial(13,(nfx-1)*ntmp +nfz)
								fric(47,nftnd0(ift),ift) = initial(3,(nfx-1)*ntmp +nfz)! Peak slip rate
								fric(44, nftnd0(ift), ift) = fric(8, nftnd0(ift), ift) !tstk0
								fric(45, nftnd0(ift), ift) = 0.0d0
							 endif
						endif !end case bp==1001.
						exit !can only be on 1 fault, thus if ynft(ift), exit do loop       
					endif  !if flt range
				enddo  !do ift
				if(ix>=2 .and. iy>=2 .and. iz>=2) then
					! if (xcoor==0.and.ycoor==100.and.zcoor==-100)then 
						! write(*,*) 'nele',nelement
					! endif
					nelement=nelement+1
					if(nelement>numel) then
						write(*,*) 'more elements in meshgen than in mesh4num'
						write(*,*) 'x,y,z',xcoor,ycoor,zcoor
						pause
					endif	
					et(nelement)=1
					ien(1,nelement)=plane1(iy-1,iz-1)
					ien(2,nelement)=plane2(iy-1,iz-1)
					ien(3,nelement)=plane2(iy,iz-1)
					ien(4,nelement)=plane1(iy,iz-1)
					ien(5,nelement)=plane1(iy-1,iz)
					ien(6,nelement)=plane2(iy-1,iz)
					ien(7,nelement)=plane2(iy,iz)
					ien(8,nelement)=plane1(iy,iz)
					center=0.0d0
					do kkk=1,8
						center(1)=center(1)+x(1,ien(kkk,nelement))
						center(2)=center(2)+x(2,ien(kkk,nelement))
						center(3)=center(3)+x(3,ien(kkk,nelement))
					enddo
					center=center/8.0d0
					mat(nelement,1) = mat0(1,1) !Vp
					mat(nelement,2) = mat0(1,2) !Vs
					mat(nelement,3) = mat0(1,3) !Rou
					mat(nelement,4) = 0.0d0
					mat(nelement,5) = 0.0d-5				
					!special treatment for using master nodes above the main fault.
					! B.D. 1/7/12
					if (center(1)<0.0d0) then 
						kinkx=dy/2.0d0
					else
						kinkx=dy/2.0d0!center(1)/30.0d3*5289.8d0
					endif				
					
					!if((center(2)>0.0d0 + kinkx - dy/2.0d0).and.(center(2)<(dy/2.0d0 + tol +kinkx))) then
					! The fault resides on x-z plane and penetrates the whole model.
					! Use undistorted mesh ycoor to locate elements on the y+ side of the fault.
					if (ycoor>0 .and. abs(ycoor-dy)<tol) then
						do i=1,nftnd0(1)
							do k=1,8
								if(ien(k,nelement)==nsmp(1,i,1)) then
									ien(k,nelement) = nsmp(2,i,1)  !use master node for the node!
								endif
							enddo
						enddo
					endif
					!if((center(2)>(-2.0d0*dy-tol+kinkx)).and.(center(2)<(2.0d0*dy+tol+kinkx))) then
					! Use undistorted mesh ycoor to locate fault elements. 
					! Assign element et to 2.
					if (abs(ycoor)<2.0d0*dy) then  
						et(nelement)=2!Fault boundary, essential boundary treatment.
						! open(1002,file='faultelem.txt',form='formatted',status='unknown',position='append')
							! write(1002,'(1x,i10,3e18.7e4,32i10)') nelement,center(1),center(2),center(3),ien(1,nelement),id(1,ien(1,nelement)),id(2,ien(1,nelement)),id(3,ien(1,nelement)),&
							! ien(2,nelement),id(1,ien(2,nelement)),id(2,ien(2,nelement)),id(3,ien(2,nelement)),&
							! ien(3,nelement),id(1,ien(3,nelement)),id(2,ien(3,nelement)),id(3,ien(3,nelement)),&
							! ien(4,nelement),id(1,ien(4,nelement)),id(2,ien(4,nelement)),id(3,ien(4,nelement)),&	
							! ien(5,nelement),id(1,ien(5,nelement)),id(2,ien(5,nelement)),id(3,ien(5,nelement)),&
							! ien(6,nelement),id(1,ien(6,nelement)),id(2,ien(6,nelement)),id(3,ien(6,nelement)),&
							! ien(7,nelement),id(1,ien(7,nelement)),id(2,ien(7,nelement)),id(3,ien(7,nelement)),&
							! ien(8,nelement),id(1,ien(8,nelement)),id(2,ien(8,nelement)),id(3,ien(8,nelement))
					endif	
					!if((center(2)<(ylinet(1)+30.0d3)).or.(center(2)>(ylinet(nyt)-30.0d3))) then
					! Locate far field elements by ycoor. Assign et = 3. 
					if (ycoor < ymin1 + 2.0d0*dymax .or. ycoor > ymax1 - 2.0d0*dymax) then  
						et(nelement)=3!Regular boundaries, essential boundary treatment.
						! open(1003,file='boundelem.txt',form='formatted',status='unknown',position='append')
							! write(1003,'(1x,i10,3e18.7e4)') nelement,center(1),center(2),center(3)						
					endif
				endif !if element
			enddo!iy
		enddo!iz
		plane1 = plane2
	enddo!ix	

	!Checking phase
	!write(*,*) 'C2'
	if (nnode/=numnp.or.nelement/=numel.or.nftnd0(1)/=nftnd(1).or.neq0/=neq) then 
		write(*,*) 'Consistancy check'
		write(*,*) 'mesh4:meshgen(nnode)',numnp,nnode
		write(*,*) 'mesh4:meshgen(nelem)',numel,nelement
		write(*,*) 'mesh4:meshgen(nnode)',nftnd0(1),nftnd(1)
		write(*,*) 'mesh4:meshgen(nnode)',neq0,neq
		stop 2
	endif 
	!for multiple faults. B.D. 1/7/12
	do ift=1,ntotft
		if(nftnd0(ift)>0) then
		!...distance from the source point
		do i=1,ifd(ift)
			do j=1,ifs(ift)
				i1 = fltrc(1,j,i,ift)
				j1 = fltrc(2,j,i,ift)
				r4nuc(j1,ift) = sqrt((x(1,i1)-xsource)**2 + (x(2,i1)-ysource)**2 &
				+ (x(3,i1)-zsource)**2)
			enddo
		enddo
		!...element areas and distribute evenly to its four nodes
		do i=2,ifd(ift)
			do j=2,ifs(ift)
			!...4 nodes of quadrilateral
				n1 = fltrc(1,j,i,ift) !use nodal number in nnode
				n2 = fltrc(1,j-1,i,ift)
				n3 = fltrc(1,j-1,i-1,ift)
				n4 = fltrc(1,j,i-1,ift)
				m1 = fltrc(2,j,i,ift) !use nodal number in nftnd0
				m2 = fltrc(2,j-1,i,ift)
				m3 = fltrc(2,j-1,i-1,ift)
				m4 = fltrc(2,j,i-1,ift)
				!...calculate area of the quadrilateral
				!......if fault is not in coor axes plane
				a1=sqrt((x(1,n2)-x(1,n1))**2 + (x(2,n2)-x(2,n1))**2 &
				+ (x(3,n2)-x(3,n1))**2)
				b1=sqrt((x(1,n3)-x(1,n2))**2 + (x(2,n3)-x(2,n2))**2 &
				+ (x(3,n3)-x(3,n2))**2)
				c1=sqrt((x(1,n4)-x(1,n3))**2 + (x(2,n4)-x(2,n3))**2 &
				+ (x(3,n4)-x(3,n3))**2)
				d1=sqrt((x(1,n1)-x(1,n4))**2 + (x(2,n1)-x(2,n4))**2 &
				+ (x(3,n1)-x(3,n4))**2)
				p1=sqrt((x(1,n4)-x(1,n2))**2 + (x(2,n4)-x(2,n2))**2 &
				+ (x(3,n4)-x(3,n2))**2)
				q1=sqrt((x(1,n3)-x(1,n1))**2 + (x(2,n3)-x(2,n1))**2 &
				+ (x(3,n3)-x(3,n1))**2)
				area=0.25d0 * sqrt(4.0d0*p1*p1*q1*q1 - &
			   (b1*b1 + d1*d1 - a1*a1 -c1*c1)**2) 
				!...distribute above area to 4 nodes evenly
				area = 0.25d0 * area
				arn(m1,ift) = arn(m1,ift) + area
				arn(m2,ift) = arn(m2,ift) + area
				arn(m3,ift) = arn(m3,ift) + area
				arn(m4,ift) = arn(m4,ift) + area			
				enddo
			enddo
			! !...save area for moment (rate) calculation before MPI. B.D. 8/11/10
			! arn4m = arn
		endif!end if nftnd0=0/ not	
	enddo!ift
end subroutine meshgen
