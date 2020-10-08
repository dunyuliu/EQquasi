subroutine bound_load

	use globalvar
	implicit none
	
	integer (kind = 4) :: nel, node_num, i, j, ntag, inv, jnv, temp, elemvar(nee), Vol
	real (kind = dp) :: elemmat(5), elemx(ndof,nen), es(nee,nee), em(nee), ef(nee), elemu(nee)
	
	right=0.0d0
	constmp=0.0d0
!write(*,*) '!---3.5.1===UPDATE BOUNDARY U*(t+1) INTO [CONSTRAINTMP]'	
!--------===STORE EQUVILENT FORCES INTO [RIGHT]
	do nel=1,numel
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
					if (itag == 0) then 
						if (elemvar(ntag)==-2) then
							consv(j,node_num)=-loadrate!-5.0d-10
							consvtmp(j,node_num) = consv(j,node_num)
							consa(j,node_num)=0.0d0
						elseif (elemvar(ntag)==-3) then
							consv(j,node_num)=loadrate!5.0d-10
							consvtmp(j,node_num) = consv(j,node_num)
							consa(j,node_num)=0.0d0					
						elseif (elemvar(ntag)==-1) then
							consv(j,node_num)=0.0d0
							consvtmp(j,node_num) = consv(j,node_num)
							consa(j,node_num)=0.0d0
						endif 

						!First predictions of U*(t+1)[CONSTRAINTMP]
						constmp(j,node_num)=cons(j,node_num)+consv(j,node_num)*dtev1
					elseif (itag == 1) then 
						!-2->-x; -3->+x; -5->fault
						if (elemvar(ntag)==-2) then
							consv(j,node_num)=-loadrate!-5.0d-10
							consa(j,node_num)=0.0d0
						elseif (elemvar(ntag)==-3) then
							consv(j,node_num)=loadrate!5.0d-10
							consa(j,node_num)=0.0d0					
						elseif (elemvar(ntag)==-1) then
							consv(j,node_num)=0.0d0
							consa(j,node_num)=0.0d0
						endif 
!write(*,*) '!---3.5.5.1-FINALL UPDATION OF [CONSTRAIN].'
!-----------NO NEED TO WORRY NON-BOUNDARY NODES, SINCE THEIR CONSTRAINV&VTMP ARE ZERO.
!-----------ERROR, CANNOT UPDATE CONTRAIN ITSELF. EACH D.O.F HAS BEEN VISITIED SEVERAL TIMES FOR IN THE LOOP.
						constmp(j,node_num) = cons(j,node_num) + 0.5d0*(consv(j,node_num)+consvtmp(j,node_num))*dtev1					
					endif 
					elemu(ntag) = constmp(j,node_num)
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
		
	enddo

end subroutine bound_load


subroutine bound_ft_ku

	use globalvar
	implicit none
	
	integer (kind = 4) :: nel, node_num, i, j, m, ntag, inv, jnv, temp, elemvar(nee)
	real (kind = dp) :: elemmat(5), elemx(ndof,nen), es(nee,nee), em(nee), ef(nee), elemu(nee), Vol 
			
	consf=0.0d0
	if (itag == 0) then 
		consm=0.0d0!Updated only when itag == 0
	endif 
	
	do nel=1,numel!Compute traction. 
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
					if (itag == 0) then 
						elemu(ntag)=constmp(j,node_num)
					elseif (itag == 1) then 
						elemu(ntag)=cons(j,node_num)
					endif 
				enddo 
			enddo		
			do i=1,nen!-KU(t+1)
				node_num=ien(i,nel)
				do j=1,ndof 
					if (id(j,node_num)<0) then 
						do m=1,nee
							consf(j,node_num)=consf(j,node_num)-es((i-1)*3+j,m)*elemu(m)
						enddo
						if (itag == 0) then 
							consm(j,node_num)=consm(j,node_num)+em((i-1)*3+j)
						endif 
					endif
				enddo 
			enddo					
		endif
	enddo!nel
end subroutine bound_ft_ku