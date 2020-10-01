SUBROUTINE hrglss(nelement,numnp,neq,ien,d,v,right,id,mat,eledet)
use globalvar
implicit none

!### program to calculate hourglass resistence and to add to
!	residual force vector for the 8-node hexahedral
!	element and assemble into the global right-hand-side 
!	vector. (one-gaussian quadrature).
!  Notice about an important assumption: infinitesimal 
!	deformation so that coordinates x() are not updated, 
!	so that lots of variables can be clculated once and 
!	save for time step loops to save time.
!  B.D. 8/21/05
integer (kind=4) :: nel,i,j,k,m,nelement,numnp,neq,non
integer (kind=4),dimension(nen,nelement)::ien
integer (kind=4),dimension(3,numnp)::id
real (kind=8),dimension(neq) :: right,d,v
real(kind=8)::f(24),Vol,coef,kapa_hg
real(kind=8),dimension(3,4)::q
real(kind=8),dimension(4,8)::fi
real(kind=8),dimension(nelement)::eledet
real(kind=8),dimension(nelement,5)::mat

kapa_hg=0.1

do nel=1,nelement
		!viscous hourglass control
		Vol=eledet(nel)
		coef=0.25*kapa_hg*mat(nel,3)*mat(nel,1)*Vol**(2./3.)
		fi(1,1)=1;fi(1,2)=1;fi(1,3)=-1;fi(1,4)=-1;fi(1,5)=-1;fi(1,6)=-1;fi(1,7)=1;fi(1,8)=1
		fi(2,1)=1;fi(2,2)=-1;fi(2,3)=-1;fi(2,4)=1;fi(2,5)=-1;fi(2,6)=1;fi(2,7)=1;fi(2,8)=-1
		fi(3,1)=1;fi(3,2)=-1;fi(3,3)=1;fi(3,4)=-1;fi(3,5)=1;fi(3,6)=-1;fi(3,7)=1;fi(3,8)=-1
		fi(4,1)=1;fi(4,2)=-1;fi(4,3)=1;fi(4,4)=-1;fi(4,5)=-1;fi(4,6)=1;fi(4,7)=-1;fi(4,8)=1
		!fi1=[1,1,-1,-1,-1,-1,1,1]
		!fi2=[1,-1,-1,1,-1,1,1,-1]
		!fi3=[1,-1,1,-1,1,-1,1,-1]
		!fi4=[1,-1,1,-1,-1,1,-1,1]
		q=0.0!q(3,4)
		f=0.0
		do i=1,3
			do j=1,4
				do k=1,8
					q(i,j)=q(i,j)+v(id(i,ien(k,nel)))*fi(j,k)
				enddo
			enddo
		enddo
		do k=1,8
			do i=1,3
				do j=1,4
					f((k-1)*3+i)=f((k-1)*3+i)-coef*q(i,j)*fi(j,k)
				enddo
			enddo
		enddo	
		do i=1,nen
			do j=1,ned
				non=ien(i,nel) 
				right(id(j,non))=right(id(j,non))+f((i-1)*3+j)
			enddo
		enddo			

enddo

end SUBROUTINE hrglss