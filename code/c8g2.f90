subroutine c8g2(elemx,material,es,em,ef,Vol)

	use globalvar
	implicit none 
	
	integer(kind=4)::i,j,ii,jj,kk
	real(kind=8)::gl(3,3,8),elemx(3,8)
	real(kind=8)::shl(3,8,8),shg(3,8,8)!nrowsh*nen*nint
	real(kind=8)::shl0(8,8),w(8)
	real(kind=8)::coef(3,3,8),jacob(8),B(6,24),DD(24,6)
	real(kind=8)::es(24,24),em(24),ef(24)
	real(kind=8)::material(5),vp,vs,rou,D0(6,6),Vol

	call shlc8g2(shl,shl0,w)
	gl    = 0.0d0
	shg   = 0.0d0
	coef  = 0.0d0
	jacob = 0.0d0
	B     = 0.0d0
	Vol   = 0.0d0
	vp    = material(1)
	vs    = material(2)
	rou   = material(3)
	D0    = 0.0d0 
	D0(1,1)=vp**2*rou 
	D0(2,2)=D0(1,1)
	D0(3,3)=D0(1,1)
	D0(1,2)=vp**2*rou-2.0d0*vs**2*rou
	D0(1,3)=D0(1,2)
	D0(2,1)=D0(1,2)
	D0(3,1)=D0(1,2)
	D0(2,3)=D0(1,2)
	D0(3,2)=D0(1,2)
	D0(4,4)=vs**2*rou
	D0(5,5)=D0(4,4)
	D0(6,6)=D0(4,4)
	es   = 0.0d0
	em   = 0.0d0
	ef   = 0.0d0
	do i=1,nint
		do j=1,nen
			gl(1,1,i)=gl(1,1,i)+elemx(1,j)*shl(1,j,i)!x_kesai
			gl(1,2,i)=gl(1,2,i)+elemx(1,j)*shl(2,j,i)!x_yita
			gl(1,3,i)=gl(1,3,i)+elemx(1,j)*shl(3,j,i)!x_fai
			gl(2,1,i)=gl(2,1,i)+elemx(2,j)*shl(1,j,i)!y_
			gl(2,2,i)=gl(2,2,i)+elemx(2,j)*shl(2,j,i)
			gl(2,3,i)=gl(2,3,i)+elemx(2,j)*shl(3,j,i)
			gl(3,1,i)=gl(3,1,i)+elemx(3,j)*shl(1,j,i)
			gl(3,2,i)=gl(3,2,i)+elemx(3,j)*shl(2,j,i)
			gl(3,3,i)=gl(3,3,i)+elemx(3,j)*shl(3,j,i)		
		enddo
		coef(1,1,i)=gl(2,2,i)*gl(3,3,i)-gl(2,3,i)*gl(3,2,i)
		coef(1,2,i)=gl(2,3,i)*gl(3,1,i)-gl(2,1,i)*gl(3,3,i)
		coef(1,3,i)=gl(2,1,i)*gl(3,2,i)-gl(2,2,i)*gl(3,1,i)
		coef(2,1,i)=gl(3,2,i)*gl(1,3,i)-gl(3,3,i)*gl(1,2,i)
		coef(2,2,i)=gl(3,3,i)*gl(1,1,i)-gl(3,1,i)*gl(1,3,i)
		coef(2,3,i)=gl(3,1,i)*gl(1,2,i)-gl(3,2,i)*gl(1,1,i)
		coef(3,1,i)=gl(1,2,i)*gl(2,3,i)-gl(1,3,i)*gl(2,2,i)
		coef(3,2,i)=gl(1,3,i)*gl(2,1,i)-gl(1,1,i)*gl(2,3,i)
		coef(3,3,i)=gl(1,1,i)*gl(2,2,i)-gl(1,2,i)*gl(2,1,i)	
		jacob(i)=coef(1,1,i)*gl(1,1,i)+coef(1,2,i)*gl(1,2,i)+coef(1,3,i)*gl(1,3,i)
		do j=1,nen 
			shg(1,j,i)=(shl(1,j,i)*coef(1,1,i)+shl(2,j,i)*coef(1,2,i)+shl(3,j,i)*coef(1,3,i))/jacob(i)
			shg(2,j,i)=(shl(1,j,i)*coef(2,1,i)+shl(2,j,i)*coef(2,2,i)+shl(3,j,i)*coef(2,3,i))/jacob(i)
			shg(3,j,i)=(shl(1,j,i)*coef(3,1,i)+shl(2,j,i)*coef(3,2,i)+shl(3,j,i)*coef(3,3,i))/jacob(i)
			B(1,3*j-2)=shg(1,j,i)
			B(2,3*j-1)=shg(2,j,i)
			B(3,3*j)=shg(3,j,i)
			B(4,3*j-1)=shg(3,j,i)
			B(4,3*j)=shg(2,j,i)
			B(5,3*j-2)=shg(3,j,i)
			B(5,3*j)=shg(1,j,i)
			B(6,3*j-2)=shg(2,j,i)
			B(6,3*j-1)=shg(1,j,i)			
		enddo
		DD=0.0d0
		do ii=1,24
			do jj=1,6
				do kk=1,6
					DD(ii,jj)=DD(ii,jj)+B(kk,ii)*D0(kk,jj)
				enddo
			enddo
		enddo		
		do ii=1,24
			do jj=1,24
				do kk=1,6
					es(ii,jj)=es(ii,jj)+DD(ii,kk)*B(kk,jj)*w(i)*jacob(i)
				enddo 
			enddo 
		enddo 
		do ii=1,nen
			do jj=1,nen
				em(3*ii-2)=em(3*ii-2)+shl0(ii,i)*shl0(jj,i)*rou*w(i)*jacob(i)
				em(3*ii-1)=em(3*ii-1)+shl0(ii,i)*shl0(jj,i)*rou*w(i)*jacob(i)
				em(3*ii)=em(3*ii)+shl0(ii,i)*shl0(jj,i)*rou*w(i)*jacob(i)		
			enddo
		enddo 

		Vol=Vol+w(i)*jacob(i)

	enddo

end subroutine c8g2