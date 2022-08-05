subroutine elemcal_aztec

	use globalvar
	implicit none
	include 'mpif.h'
	
	integer (kind = 4) :: tocalculate, nel, node_num, i, j, m, ntag, inv, jnv, temp, elemvar(nee),k
	real (kind = dp) :: elemmat(5), elemx(ndof,nen), es(nee,nee), em(nee), ef(nee), Vol
	
	
	do  nel=1,numel
		tocalculate = 0
		ntag=0
		do i=1,nen
			node_num=ien(i,nel)
			do j=1,ndof
			ntag=ntag+1
			elemvar(ntag)=id(j,node_num)
			!ERROR: I cannot simply loop over elements involving minneq-maxneq, because 
			! the columns inlove other elements.
			 if (((elemvar(ntag)-1)>=minneq-5 .and. (elemvar(ntag)-1)<=maxneq+5)) then
				 tocalculate = 1
			 endif
			!elemu(iii)=constrain(j,node_num)
			enddo
		enddo
		if (tocalculate == 1) then
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
			
			call c8g2(elemx,elemmat,es,em,ef,Vol)
			eledet(nel)=Vol

			do i=1,nee
				inv=elemvar(i)
				if (inv.gt.0) then
					mass(inv)=mass(inv)+em(i)
				endif
			enddo
			  
			do i=1,nee
				inv=elemvar(i)-1 
				do m = 0,N_update-1
					if (inv == update(m)) then
						val(m) = val(m) + es(i,i)
						do j=1,nee       
							jnv=elemvar(j)-1
							do k = bindx(m),bindx(m+1)-1
								if (jnv == bindx(k)) then
									val(k) = val(k) + es(i,j)
								endif
							enddo
						enddo
					endif
				enddo 
			enddo

		endif!tocalculate == 1
	enddo	
end subroutine elemcal_aztec