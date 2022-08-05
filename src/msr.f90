subroutine msr

	use globalvar
	implicit none
	include 'mpif.h'
	
	integer (kind = 4) :: i, j, nel, ntag, temp, var, node_num, inv, jnv, mm, m, mloc, elemvar(nee)
	
	
	allocate(lie(neq,100))
	do i=1,neq
		num(i)=0
		do j=1,100
			lie(i,j)=0
		enddo
	enddo

	do nel=1,numel
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
			! MSR storage required by AZTEC is for the whole submatrix. 
			do j=1,ntag !previous loop starts at i to ntag. Now, treate the whole matrix.
				inv=elemvar(i)
				jnv=elemvar(j)
				! if (jnv.lt.inv) then
					! temp=jnv
					! jnv=inv
					! inv=temp
				! endif

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
	
end subroutine msr