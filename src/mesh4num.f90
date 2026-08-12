subroutine mesh4num

	use globalvar
	implicit none
	include 'mpif.h'

	integer(kind=4)::nnode0,nelement0,neq0,ix,iy,iz,&
		edgex1,edgey1,edgez1,edgezn,&
		i,j,k,ift

	integer (kind=4),dimension(ntotft)::nftnd0
	real(kind=8)::dy,dz
	integer (kind=4)::nxuni,nyuni,nzuni

	real(kind=8)::tol,xcoor,ycoor,zcoor,xstep,ystep,zstep,&
			xmin1,xmax1,ymin1,ymax1,zmin1
	real(kind=8),allocatable,dimension(:)::xlinet,ylinet,zlinet,xline,yline,zline
	integer(kind=4)::i1
	real(kind=8)::fltx1,fltx2,fltz1,fltz2,yflt_lo,yflt_hi


	dy=dx
	dz=dx
	tol=dx/100.d0

	! Uniform-grid region must span the union of every fault's x/z extent, not
	! just fault 1's, so ntotft > 1 exercises the same statements ntotft == 1
	! does (fltx1/fltx2/fltz1/fltz2 collapse to fault 1's own box when
	! ntotft == 1).
	fltx1=minval(fltxyz(1,1,:))
	fltx2=maxval(fltxyz(2,1,:))
	fltz1=minval(fltxyz(1,3,:))
	fltz2=maxval(fltxyz(2,3,:))

	nxuni=(fltx2-fltx1-2.0d0*dx)/dx+1
	xstep=dx
	xcoor=fltx1+dx
	do ix=1,np
		xstep=xstep*1.0d0
		xcoor=xcoor-xstep
		if(xcoor<=xmin) exit
	enddo
	edgex1=ix
	xstep=dx
	xcoor=fltx2-dx
	do ix=1,np
		xstep=xstep*1.0d0
		xcoor=xcoor+xstep
		if(xcoor>=xmax) exit
	enddo
	nxt=nxuni+edgex1+ix
	allocate(xlinet(nxt))
	!predetermine x-coor
	xlinet(edgex1+1)=fltx1+dx
	xstep=dx
	do ix=edgex1,1,-1
		xstep=xstep*1.0d0
		xlinet(ix)=xlinet(ix+1)-xstep
	enddo
	xmin1=xlinet(1)
	do ix=edgex1+2,edgex1+nxuni
		xlinet(ix)=xlinet(ix-1)+dx
	enddo
	xstep=dx
	do ix=edgex1+nxuni+1,nxt
		xstep=xstep * 1.0d0
		xlinet(ix)=xlinet(ix-1)+xstep
	enddo
	xmax1=xlinet(nxt)
	!Y
	! Uniform-dy belt spans the union of every fault's y position, not just
	! fault 1's -- mirrors the x/z generalization above and must stay
	! consistent with meshgen.f90, which this mirrors. ntotft==1 collapses
	! yflt_lo==yflt_hi==0 and this reduces to the original single-fault belt.
	yflt_lo=minval(fltxyz(1,2,:))
	yflt_hi=maxval(fltxyz(2,2,:))
	nyuni=int((yflt_hi-yflt_lo)/dy+0.5d0)+dis4uniF+dis4uniB+1
	ystep=dy
	ycoor=yflt_lo-dy*(dis4uniF)
	do iy=1,np
		ystep=ystep*rat
		if (ystep>=dymax) ystep = dymax
		ycoor=ycoor-ystep
		if(ycoor<=ymin) exit
	enddo
	edgey1=iy
	ystep=dy
	ycoor=yflt_hi+dy*(dis4uniB)
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
	ylinet(edgey1+1)=yflt_lo-dy*(dis4uniF)
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
	zcoor=fltz1+dx
	do iz=1,np
		zstep=zstep
		zcoor=zcoor-zstep
		if(zcoor<=zmin) exit
	enddo
	edgezn=iz
	nzuni=(fltz2-fltz1-dx)/dx+1
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

	nnode0    = 0
	nelement0 = 0
	neq0      = 0
	nftnd0    = 0
	
	do ix=1,nxt
		do iz=1,nzt
			do iy=1,nyt
				xcoor=xlinet(ix)
				ycoor=ylinet(iy)
				zcoor=zlinet(iz)
				nnode0=nnode0 + 1
				! No equation number on a domain edge (iy==1/nyt) or on any
				! fault's own y-plane -- one check per fault, not just y==0.
				if(iy==1.or.iy==nyt.or.any(abs(ycoor-fltxyz(1,2,1:ntotft))<tol)) then

				else
					do i1=1,ndof
						neq0=neq0+1
					enddo
				endif
				do ift=1,ntotft 
					if(xcoor>=(fltxyz(1,1,ift)-tol).and.xcoor<=(fltxyz(2,1,ift)+tol).and. &
						ycoor>=(fltxyz(1,2,ift)-tol).and.ycoor<=(fltxyz(2,2,ift)+tol).and. &
						zcoor>=(fltxyz(1,3,ift)-tol) .and. zcoor<=(fltxyz(2,3,ift)+tol)) then			
						nftnd0(ift) = nftnd0(ift) + 1
						nnode0=nnode0 + 1
						!...establish equation numbers for this master node
						!do i1=1,ndof
							!neq0=neq0 + 1
						!enddo
						exit !can only be on 1 fault, thus if ynft(ift), exit do loop       
					endif
				enddo  !do ift						
				if(ix>=2.and.iy>=2.and.iz>=2) then
					nelement0=nelement0+1
				endif !if element
			enddo!iy
		enddo!iz
	enddo!ix

	numnp = nnode0
	numel = nelement0
	neq   = neq0 
	nftnd = nftnd0

end subroutine mesh4num
