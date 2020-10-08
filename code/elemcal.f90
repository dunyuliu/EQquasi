subroutine elemcal

	use globalvar
	implicit none
	include 'mpif.h'
	
	integer (kind = 4) :: nel, node_num, i, j, m, ntag, inv, jnv, temp, elemvar(nee)
	real (kind = dp) :: elemmat(5), elemx(ndof,nen), es(nee,nee), em(nee), ef(nee), Vol
	
	do  nel=1,numel	
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
	
		call c8g2(elemx,elemmat,es,em,ef,Vol)
		
		eledet(nel)=Vol
	
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
	
end subroutine elemcal