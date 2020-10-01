subroutine meshgen(nnode0,nelement0,neq0,ien,mat,et,x,id,&
			nftnd1,n4onf,xonfs,nonfs,&
			nftmx,nonmx,nsmp,un,us,ud,fric,arn,r4nuc,&
			anonfs,me)
use globalvar
implicit none
!include 'mpif.h'	
	
integer(kind=4)::nnode0,nelement0,nnode,nelement,neq0,neq,nxt,nyt,nzt,ix,iy,iz,&
	edgex1,edgey1,edgez1,edgezn,&
	i,j,k,i1,j1,np=10000,ift,me	
real(kind=8)::dy,dz
real(kind=8)::tmp1,tmp2
integer (kind=4)::nxuni,nyuni,nzuni	
real(kind=8)::tol,xcoor,ycoor,zcoor,xstep,ystep,zstep,&
		xmin1,xmax1,ymin1,ymax1,zmin1
real(kind=8),allocatable,dimension(:)::xlinet,ylinet,zlinet,xline,yline,zline
integer(kind=4),allocatable,dimension(:,:) :: plane1,plane2
integer(kind=4),allocatable,dimension(:,:,:,:) :: fltrc	
integer(kind=4)::id(ndof,nnode0),ien(8,nelement0),et(nelement0)
real(kind=8)::mat(nelement0,5),x(ned,nnode0)
!...fault definition: 
integer(kind=4)::nftmx,nonmx,n1,n2,n3,n4,m1,m2,m3,m4
real(kind=8)::a1,b1,c1,d1,p1,q1,area 
integer (kind=4),dimension(ntotft) :: nftnd,nftnd1,ixfi,izfi,ifs,ifd
real (kind=8) :: x1,x2,x3=0.d0,y1,y2,y3=0.d0,z1,z2,z3=0.d0 !3 points for fault plane (inc origin)
integer (kind=4),dimension(3,nonmx) :: anonfs
integer (kind=4),dimension(2,nftmx,ntotft) :: nsmp
real (kind=8),dimension(nftmx,ntotft) :: fnft,arn,r4nuc,slp4fri
real (kind=8),dimension(3,nftmx,ntotft) :: un,us,ud
real (kind=8),dimension(50,nftmx,ntotft) :: fric

!...output on-fault stations. B.D. 10/25/09
integer (kind=4) :: n4onf
integer,dimension(ntotft)::nonfs
real (kind=8),dimension(2,9,ntotft) :: xonfs
!...output node coors. 6 Off-fault stations in this example. B.D. 10/25/11
integer (kind=4)::n4nds=6,n4out
integer (kind=4),dimension(6)::n4yn=0
integer (kind=4),dimension(2,6)::an4nds
!For TPV8
real (kind=8),dimension(2,6)::x4nds=reshape((/0.0,-3.0,&!8.0,6.0
												0.0,-2.0,&
												0.0,-1.0,&
												0.0,1.0,&
												12.0,-3.0,&
												12.0,3.0/),(/2,6/))
real(kind=8)::center(3),ycoort,strtmp,kinkx	
integer(kind=4)::kkk	

										
!on-fault station coordinates are given here 
!(along strike, z).
xonfs=reshape((/-18.0,-9.0,&
				-15.0,-9.0,&
				-12.0,-9.0,&
				-9.0,-9.0,&
				-6.0,-9.0,&
				6.0,-9.0,&
				12.0,-9.0,&
				18.0,-9.0,&
				0.0,-9.0/),(/2,9,1/))
xonfs = xonfs * 1000.0d0  !convert from km to m

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
	ycoor=ycoor-ystep
	if(ycoor<=ymin) exit
enddo
edgey1=iy
ystep=dy
ycoor=dy*(dis4uniB)
do iy=1,np
	ystep=ystep*rat
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
	ylinet(iy)=ylinet(iy+1)-ystep
enddo
ymin1=ylinet(1)
do iy=edgey1+2,edgey1+nyuni
	ylinet(iy)=ylinet(iy-1)+dy
enddo
ystep=dy
do iy=edgey1+nyuni+1,nyt
	ystep=ystep*rat
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
neq = 0
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
			ycoort=ycoor
			strtmp=270.0d0/180.0d0*pi
                        !The normal vector at x = 0 belongs to the left segment.
			if (xcoor>0.0d0) then 

				tmp1 = (xcoor--0.0d3)*dtand(ttheta)
				tmp2=ttheta
				strtmp=(270.0d0-tmp2)/180.0d0*pi		
				if (ycoor>=0.0d3) then 	
					ycoort=(ylinet(nyt)-tmp1)*ycoor/ylinet(nyt) + tmp1
				elseif (ycoor<0.0d0) then 
					ycoort=ylinet(1)+(ycoor-ylinet(1))/(-ylinet(1))*(tmp1-ylinet(1))
				endif
				x(2,nnode)=ycoort							
			endif
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
					neq=neq+1
					id(i1,nnode)=neq!Only nodes inside two blocks have equation numbers.
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
						! write(*,*) 'SLAVE:-18,-21,nnode',nnode,nftnd(ift),ifs(ift),ifd(ift) 
						! write(*,*) '-18,-21,id',id(1,nnode),id(2,nnode),id(3,nnode)
					! endif					
					nftnd(ift) = nftnd(ift) + 1
					if(nftnd1(ift)==0) then
					write(*,*) 'inconsistency between mesh4num and meshgen for nftnd'
					write(*,*) 'current node x,y,z:',xcoor,ycoor,zcoor
					stop 2001
					endif
					nsmp(1,nftnd(ift),ift) = nnode !slave node						
					nnode = nnode + 1

					nsmp(2,nftnd(ift),ift) = nnode !master node, on the +y side.
					plane2(nyt+ift,iz) = nnode
					x(1,nnode) = xcoor
					x(2,nnode) = ycoort
					x(3,nnode) = zcoor  
					!...identify output on-fault stations. B.D. 10/25/09
					!...fault-dependent for multiple faults. B.D. 1/8/
					do i1=1,nonfs(ift)
						if(abs(xcoor-xonfs(1,i1,ift))<tol .and. &
						abs(zcoor-xonfs(2,i1,ift))<tol) then
							n4onf = n4onf + 1
							anonfs(1,n4onf) = nftnd(ift)
							anonfs(2,n4onf) = i1
							anonfs(3,n4onf) = ift
							exit         
						endif
					enddo  
					!...establish equation numbers for this master node
					do i1=1,ndof
						!neq = neq + 1
						id(i1,nnode)=-5!assume fault boundary nodes
					enddo
					!...assign unit vectors to the split node pair						
					un(1,nftnd(ift),ift) = dcos(strtmp)*dsin(fltxyz(2,4,ift))
					un(2,nftnd(ift),ift) = -dsin(strtmp)*dsin(fltxyz(2,4,ift))
					un(3,nftnd(ift),ift) = dcos(fltxyz(2,4,ift))		
					us(1,nftnd(ift),ift) = -dsin(strtmp)
					us(2,nftnd(ift),ift) = -dcos(strtmp)
					us(3,nftnd(ift),ift) = 0.0d0
					ud(1,nftnd(ift),ift) = dcos(strtmp)*dcos(fltxyz(2,4,ift))
					ud(2,nftnd(ift),ift) = dsin(strtmp)*dcos(fltxyz(2,4,ift))
					ud(3,nftnd(ift),ift) = dsin(fltxyz(2,4,ift))
					!...prepare for area calculation
					if(ixfi(ift)==0) ixfi(ift)=ix
					if(izfi(ift)==0) izfi(ift)=iz
					ifs(ift)=ix-ixfi(ift)+1
					ifd(ift)=iz-izfi(ift)+1
					fltrc(1,ifs(ift),ifd(ift),ift) = nnode	!master node
					fltrc(2,ifs(ift),ifd(ift),ift) = nftnd(ift) !fault node num in sequence
					! if (xcoor==-18.0d3.and.zcoor==-21.0d3) then 
						! write(*,*) 'Master:-18,-21,nnode',nnode,nftnd(ift),ifs(ift),ifd(ift) 
						! write(*,*) '-18,-21,id',id(1,nnode),id(2,nnode),id(3,nnode)
					! endif
					!...assign friction parameters for SCEC TPV19. B.D. 1/8/12
					!...now for Ma and Andrews (2010) model. B.D. 6/1/12
				!For TPV8	
					fric(1,nftnd(ift),ift) = 0.18d0!mus
					fric(2,nftnd(ift),ift) = 0.12d0!mud
					! if (xcoor<-20e3.and.xcoor>-40e3.and.zcoor<-3e3.and.zcoor>-8e3)then
						! fric(2,nftnd(ift),ift) = 0.5   !mud
					! endif						
					fric(3,nftnd(ift),ift) = 0.30d0    !D0
					fric(4,nftnd(ift),ift) = 0.0d6  !cohesion
					fric(5,nftnd(ift),ift) = 0.08d0  !Not used in TPV8 !T for forced rupture,i.e.,nucleation
					fric(6,nftnd(ift),ift) = 0.0d0!rhow*grav*(-zcoor) !pore pressure						
					! fric(7,nftnd(ift),ift) = 7378.0*zcoor! Initial normal stress for ElasticVersion
					! if (abs(zcoor)<tol)then
						! fric(7,nftnd(ift),ift)=-7378.0*dx/4!
					! endif					
					! fric(8,nftnd(ift),ift)=0.55*abs(fric(7,nftnd(ift),ift))!Initial shear stress for ElasticVersion
					! if (abs(xcoor)<=1500.and.abs(zcoor--12e3)<=1500.)then
						! fric(8,nftnd(ift),ift) =1.0e6+1.005*0.760*abs(7378.0*zcoor)
					! endif		
					tmp2 = 270.0d0*pi/180.0d0 - strtmp
					
					fric(7,nftnd(ift),ift) = -50.0d6 - 30.0d6 * dcos(2.0d0*(45.0d0*pi/180.0d0 - tmp2))
					fric(8,nftnd(ift),ift) = 30.0d6 * dsin(2.0d0*(45.0d0*pi/180.0d0 - tmp2))
                                        fric(44, nftnd(ift), ift) = fric(8, nftnd(ift), ift) !tstk0
                                        fric(45, nftnd(ift), ift) = 0.0d0
					! if (xcoor<-10.0d3.and.xcoor>-15.0d3.and.zcoor<-2.67d3.and.zcoor>-14.34d3) then
						! fric(8,nftnd(ift),ift)=30.0d6*1.02d0
					! endif
!TPV104 2016.10.05 D.Liu
					if (friclaw==3.or.friclaw==4)then
						! if (abs(xcoor)<=12.d3) then 
							! tmp1=1.0d0
						! elseif ((abs(xcoor)<=30.d3).and.(abs(xcoor)>=12.d3)) then 
							! !tmp1=0.5d0*(1.d0+dtanh(5.d3/(abs(xcoor)-18.d3)+5.d3/(abs(xcoor)-12.d3)))
							! tmp1=(30.0d3-abs(xcoor))/18.0d3
						! else
							! tmp1=0.0d0
						! endif					
						! if (zcoor>=-12.0d3) then 
							 ! tmp2=1.0d0
						 ! elseif (zcoor<=-12.0d3.and.zcoor>=-30.0d3) then
							! tmp2=(30.0d3-abs(zcoor))/18.0d3
							! ! tmp2=0.5d0*(1.d0+dtanh(6.d3/(abs(zcoor--7.5d3)-13.5d3)+3.d3/(abs(zcoor--7.5d3)-7.5d3)))
						 ! else
							 ! tmp2=0.0d0
						 ! endif								
						! fric(9,nftnd(ift),ift)=0.008d0+0.024d0*(1.d0-tmp1*tmp2)!a
						! if (xcoor==17e3.and.zcoor==-17e3) then 
							! write(*,*) 'tmp1',fric(9,nftnd(ift),ift),tmp1,tmp2 
						! endif
						! fric(10,nftnd(ift),ift)=0.012d0!b 
						! fric(11,nftnd(ift),ift)=0.008d0!RSF critical distance.
						
						!a Lapusta&Liu(2009)
						if (abs(zcoor)<=4.0d3) then 
							fric(9,nftnd(ift),ift) = 0.007d0 + 0.004d0*(4.0d3-abs(zcoor))/4.0d3
						elseif (abs(zcoor)>4.0d3.and.abs(zcoor)<=17.5d3) then 
							fric(9,nftnd(ift),ift) = 0.007d0
						elseif (abs(zcoor)>17.5d3.and.abs(zcoor)<=24.0d3) then
							fric(9,nftnd(ift),ift) = 0.007d0 + 0.009d0*(abs(zcoor)-17.5d3)/6.5d3 
						elseif (abs(zcoor)>24.0d3) then
							fric(9,nftnd(ift),ift) = 0.016d0
						endif
						if (abs(zcoor)<=4.0d3) then 
							tmp1 = -0.004d0 + 0.008d0*(4.0d3-abs(zcoor))/4.0d3
							fric(10,nftnd(ift),ift) = fric(9,nftnd(ift),ift) - tmp1
						elseif (abs(zcoor)>4.0d3.and.abs(zcoor)<=13.5d3) then 
							tmp1 = -0.004d0 
							fric(10,nftnd(ift),ift) = fric(9,nftnd(ift),ift) - tmp1 
						elseif (abs(zcoor)>13.5d3.and.abs(zcoor)<=17.5d3) then
							tmp1 = -0.004d0 + 0.011d0*(abs(zcoor)-13.5d3)/4.0d3
							fric(10,nftnd(ift),ift) = fric(9,nftnd(ift),ift) - tmp1  
						elseif (abs(zcoor)>17.5d3.and.abs(zcoor)<=24.0d3) then
							tmp1 = 0.007d0 + 0.009d0*(abs(zcoor)-17.5d3)/6.5d3 
							fric(10,nftnd(ift),ift) = fric(9,nftnd(ift),ift) - tmp1
						elseif (abs(zcoor)>24.0d3) then
							tmp1 = 0.007d0 + 0.009d0
							fric(10,nftnd(ift),ift) = fric(9,nftnd(ift),ift) - tmp1
						endif
						if (abs(xcoor)>23.0d3.and.abs(xcoor)<28.0d3) then 
							!fric(9,nftnd(ift),ift) = fric(9,nftnd(ift),ift)+0.009d0*( abs(xcoor)-15.d3)/5.0d3
							fric(10,nftnd(ift),ift) = fric(10,nftnd(ift),ift)*(28.0d3 - abs(xcoor))/5.0d3
						elseif (abs(xcoor)>=28.0d3) then 
							!fric(9,nftnd(ift),ift) = fric(9,nftnd(ift),ift)+0.009d0*( abs(xcoor)-15.d3)/5.0d3
							fric(10,nftnd(ift),ift) = 0.000d0
						endif 
						!if (fric(10,nftnd(ift),ift)<0.001d0)then 
						!	fric(10,nftnd(ift),ift) = 0.001d0
						!endif 
						
						fric(11,nftnd(ift),ift)=0.011d0 + 0.016d0*max(0.0d0, (4.0d3-abs(zcoor))/4.0d3)!RSF critical distance.
						fric(12,nftnd(ift),ift)=1.0d-6!RSF:V0
						fric(13,nftnd(ift),ift)=0.6d0!RSF:miu0
						if(friclaw==4) then 
							fric(14,nftnd(ift),ift)  = 0.2d0 !RSF: fw for strong rate weakenging
							fric(15,nftnd(ift),ift)  = 0.1d0+0.9d0*(1.d0-tmp1*tmp2) !RSF: Vw for strong rate weakening
						endif	
						! fric(16,nftnd(ift),ift)=0.0d0 !RSF: initial normal slip rate 
						! fric(17,nftnd(ift),ift)=1.0d-9!RSF:s 
						! fric(18,nftnd(ift),ift)=0.0d0!RSF:d
						! fric(19,nftnd(ift),ift)=1.0d-9!RSF:mag	
						!if (zcoor>=-12.0d3.and.abs(xcoor)<=12.0d3) then
							fric(16,nftnd(ift),ift)=0.0d0 !RSF: initial normal slip rate 
							fric(17,nftnd(ift),ift)=1.0d-12!RSF:s 
							fric(18,nftnd(ift),ift)=0.0d0!RSF:d
							fric(19,nftnd(ift),ift)=1.0d-12!RSF:mag	
						!	if( fric(9,nftnd(ift),ift) - fric(10,nftnd(ift),ift) > 0.0001d0 ) then
						!		fric(16,nftnd(ift),ift)=0.0d0 !RSF: initial normal slip rate 
						!		fric(17,nftnd(ift),ift)=1.0d-9!RSF:s 
						!		fric(18,nftnd(ift),ift)=0.0d0!RSF:d
						!		fric(19,nftnd(ift),ift)=1.0d-9!RSF:mag	
						!	end if
						!endif		
						
					endif		
					!RSF: initalize state variable based on the inital friction and inital slip rate  B.L. 1/8/16
					if(friclaw == 3) then 
						!Theta=d0/v0*dexp(a*dlog(2*dsinh()))
						!fric-20th dimension stores the state variable. 2017.2.2. D.Liu
						 !fric(20,nftnd(ift),ift) = fric(11,nftnd(ift),ift)/fric(12,nftnd(ift),ift)*dexp((fric(9,nftnd(ift),ift)*dlog(2.0d0*dsinh &
						 !(sqrt(fric(8,nftnd(ift),ift)**2+0.0d0**2)/abs(fric(7,nftnd(ift),ift))/fric(9,nftnd(ift),ift)) &
						 !)-fric(13,nftnd(ift),ift)-fric(9,nftnd(ift),ift)*dlog(sqrt((fric(16,nftnd(ift),ift))**2+ &
						 !(fric(17,nftnd(ift),ift))**2+(fric(18,nftnd(ift),ift))**2)/fric(12,nftnd(ift),ift)))/fric(10,nftnd(ift),ift))
						 ! if (xcoor==-18.0d3.and.zcoor==-21.0d3) then 
							! write(*,*) 'state variable',fric(20,nftnd(ift),ift)
						 ! endif 
						fric(20,nftnd(ift),ift)=fric(11,nftnd(ift),ift)/1.0d-9!fric(19,nftnd(ift),ift)
						fric(8,nftnd(ift),ift) = - fric(7,nftnd(ift),ift) * (fric(13,nftnd(ift),ift) + (fric(9,nftnd(ift),ift) - fric(10,nftnd(ift),ift))*dlog(1.0d-9/fric(12,nftnd(ift),ift)))
						! if( fric(9,nftnd(ift),ift) - fric(10,nftnd(ift),ift) > 0.0001d0 ) then
							! fric(20,nftnd(ift),ift)=fric(11,nftnd(ift),ift)/1.0d-9!					
						! end if	
							
						! if (xcoor<-10.0d3.and.xcoor>-15.0d3.and.zcoor<-2.67d3.and.zcoor>-14.34d3) then
							! fric(20,nftnd(ift),ift)=fric(11,nftnd(ift),ift)/5.0d-8
						! endif						
					elseif(friclaw == 4) then
						!Theta0=a*log(2V0/Vini*sinh(TAOini/a/SigmaNini))
						! fric(20,nftnd(ift),ift)=fric(9,nftnd(ift),ift)*dlog(2.0d0*fric(12,nftnd(ift),ift)/sqrt((fric(16,nftnd(ift),ift))**2+(fric(17,nftnd(ift),ift))**2+(fric(18,nftnd(ift),ift))**2) &
						! *dsinh(sqrt(fric(8,nftnd(ift),ift)**2+0.0d0**2)/abs(fric(7,nftnd(ift),ift))/fric(9,nftnd(ift),ift)))
						fric(20,nftnd(ift),ift)=fric(11,nftnd(ift),ift)/fric(19,nftnd(ift),ift)
					endif		
						fric(21,nftnd(ift),ift) = 0.0d0 ! theta_dot = 0
						fric(22,nftnd(ift),ift) = 0.0d0 ! theta*(t+1) = 0
					!special values below.						
					if(abs(xcoor-fltxyz(1,1,ift))<tol .or. abs(xcoor-fltxyz(2,1,ift))<tol &
					.or. abs(zcoor-fltxyz(1,3,ift))<tol) then !-x,+x,-z for 1, +x,-z for 2
						fric(1,nftnd(ift),ift) = 10000.d0	!fault edge, pinned
					endif         
					exit !can only be on 1 fault, thus if ynft(ift), exit do loop       
				endif  !if flt range
			enddo  !do ift 			
			if(ix>=2 .and. iy>=2 .and. iz>=2) then
				! if (xcoor==0.and.ycoor==100.and.zcoor==-100)then 
					! write(*,*) 'nele',nelement
				! endif
				nelement=nelement+1
				if(nelement>nelement0) then
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
				mat(nelement,1)=6000.d0
				mat(nelement,2)=3464.d0
				mat(nelement,3)=2670.d0
				mat(nelement,4)=0.0d0
				mat(nelement,5)=0.0d-5				
				!special treatment for using master nodes above the main fault.
				! B.D. 1/7/12
				if (center(1)<0.0d0) then 
					kinkx=dy/2.0d0
				else
					kinkx=center(1)/30.0d3*5289.8d0
				endif				
				
				if((center(2)>0.0d0 + kinkx - dy/2.0d0).and.(center(2)<(dy/2.0d0 + tol +kinkx))) then
					! if (me==0) then
						! open(1000,file='Slave-Master-Elem-Replace.txt',form='formatted',status='unknown',position='append')
							! write(1000,'(1x,i10,3e18.7e4)') nelement,center(1),center(2),center(3)
					! endif
					do i=1,nftnd(1)
						do k=1,8
							if(ien(k,nelement)==nsmp(1,i,1)) then
								ien(k,nelement) = nsmp(2,i,1)  !use master node for the node!
							endif
						enddo
					enddo
				endif				
				if((center(2)>(-2.0d0*dy-tol+kinkx)).and.(center(2)<(2.0d0*dy+tol+kinkx))) then
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
				if((center(2)<(ylinet(1)+15.0d0*dx)).or.(center(2)>(ylinet(nyt)-15.0d0*dx))) then
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
if (nnode/=nnode0.or.nelement/=nelement0.or.nftnd(1)/=nftnd1(1).or.neq/=neq0) then 
	write(*,*) 'Consistancy check'
	write(*,*) 'mesh4:meshgen(nnode)',nnode0,nnode
	write(*,*) 'mesh4:meshgen(nelem)',nelement0,nelement
	write(*,*) 'mesh4:meshgen(nnode)',nftnd1(1),nftnd(1)
	write(*,*) 'mesh4:meshgen(nnode)',neq0,neq
	stop 2
endif 
!for multiple faults. B.D. 1/7/12
do ift=1,ntotft
	if(nftnd(ift)>0) then
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
			m1 = fltrc(2,j,i,ift) !use nodal number in nftnd
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
	endif!end if nftnd=0/ not	
enddo!ift	


end subroutine meshgen
