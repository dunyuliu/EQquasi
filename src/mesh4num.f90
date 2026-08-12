subroutine mesh4num

	use globalvar
	implicit none
	include 'mpif.h'

	! Explicit interface: build_yline_belt's ylinet is intent(out),
	! allocatable. func_lib.f90 has no module wrapper (matches this file's
	! own external-subroutine style), so without this block the call below
	! has an IMPLICIT interface and gfortran does not pass the allocatable
	! descriptor correctly -- it silently corrupts memory instead of
	! erroring, and the crash surfaces as a segfault inside
	! build_yline_belt itself.
	interface
		subroutine build_yline_belt(ylinet, nyt)
			use globalvar, only: dp
			real (kind = dp), allocatable, intent(out) :: ylinet(:)
			integer (kind = 4), intent(out) :: nyt
		end subroutine build_yline_belt
	end interface

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
	! Uniform-dy belt spans [min fault y, max fault y], unchanged from
	! v1.5.0. build_yline_belt (func_lib.f90) additionally REFUSES at setup
	! if any declared fault's y-offset from the belt origin is not an
	! integer multiple of dy, rather than silently meshing that fault with
	! zero nodes. Shared with meshgen.f90 so the two stay consistent by
	! construction rather than by two hand-kept copies. ntotft==1 has one
	! fault, trivially an integer (zero) multiple of dy from itself, and
	! reduces to the original single-fault belt exactly.
	call build_yline_belt(ylinet, nyt)
	ymin1=ylinet(1)
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

	! Backstop: a declared fault with zero nodes is never valid. This is
	! belt-and-braces with build_yline_belt's up-front commensurability
	! refusal (func_lib.f90) -- that check should make nftnd(ift)==0
	! impossible by construction, but this one does not trust that and
	! catches ANY route to a zero-node fault, not only a non-commensurate
	! y-offset. It silently ran to completion once (v1.5.0's belt bug: fault
	! A landed off the mesh line, meshed with nftnd(1)==0, wrote an
	! all-zero slab to fault.*.nc, and ran 4240 steps to 11 m of slip on a
	! single fault with no partner -- the only hint was "Fault nodes = 0"
	! in the run summary, which read as cosmetic). Stop here, at the
	! earliest point nftnd is known, naming the fault, its y, and dy.
	do ift = 1, ntotft
		if (nftnd(ift) == 0) then
			if (me == 0) then
				write(*,*) '====================================================================='
				write(*,*) '=                          MESH ERROR                               ='
				write(*,*) '= Declared fault has zero fault nodes.                              ='
				write(*,*) '= A declared fault with no nodes is never valid -- refusing to run   ='
				write(*,*) '= rather than silently producing an unpartnered fault.               ='
				write(*,*) '=                                                                    ='
				write(*,'(X,A,I4,4X,A)')      '=   fault index = ', ift, '='
				write(*,'(X,A,E16.8,A,4X,A)') '=   fault y     = ', fltxyz(1,2,ift), ' m', '='
				write(*,'(X,A,E16.8,A,4X,A)') '=   dx = dy     = ', dx, ' m', '='
				write(*,*) '=                                                                    ='
				write(*,*) '= Check fltxyz (case_input compset faultgeom) against dx/dy and the  ='
				write(*,*) '= model domain (fymin/fymax).                                        ='
				write(*,*) '====================================================================='
			endif
			stop 6
		endif
	enddo

end subroutine mesh4num
