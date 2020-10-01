subroutine shlc8g2(shl,shl0,w)
use globalvar
implicit none 
integer(kind=4)::i
real(kind=8)::shl(3,8,8)!nrowsh*nen*nint
real(kind=8)::shl0(8,8)!nrowsh*nen*nint
!1st-3rd:x,y,z derivatives of local shape function
!4th:shape function itself
!NO.1:-1,-1,-1; NO.5:-1,-1, 1;
!NO.2: 1,-1,-1; NO.6: 1,-1, 1;
!NO.3: 1, 1,-1; NO.7: 1, 1, 1; 
!NO.4:-1, 1,-1; NO.8:-1, 1, 1;
!-----------------------------
real(kind=8)::w(8),cst
real(kind=8)::gsloc(3,8)
!Gaussian points have the same order
cst=1.0d0/8.0d0!
gsloc(1,1)=-1.0d0/sqrt(3.0d0)!NO.1
gsloc(2,1)=-1.0d0/sqrt(3.0d0)
gsloc(3,1)=-1.0d0/sqrt(3.0d0)
w(1)=1.0d0
gsloc(1,2)=1.0d0/sqrt(3.0d0)!NO.2
gsloc(2,2)=-1.0d0/sqrt(3.0d0)
gsloc(3,2)=-1.0d0/sqrt(3.0d0)
w(2)=1.0d0
gsloc(1,3)=1.0d0/sqrt(3.0d0)!NO.3
gsloc(2,3)=1.0d0/sqrt(3.0d0)
gsloc(3,3)=-1.0d0/sqrt(3.0d0)
w(3)=1.0d0
gsloc(1,4)=-1.0d0/sqrt(3.0d0)!NO.4
gsloc(2,4)=1.0d0/sqrt(3.0d0)
gsloc(3,4)=-1.0d0/sqrt(3.0d0)
w(4)=1.0d0
gsloc(1,5)=-1.0d0/sqrt(3.0d0)!NO.5
gsloc(2,5)=-1.0d0/sqrt(3.0d0)
gsloc(3,5)=1.0d0/sqrt(3.0d0)
w(5)=1.0d0
gsloc(1,6)=1.0d0/sqrt(3.0d0)!NO.6
gsloc(2,6)=-1.0d0/sqrt(3.0d0)
gsloc(3,6)=1.0d0/sqrt(3.0d0)
w(6)=1.0d0
gsloc(1,7)=1.0d0/sqrt(3.0d0)!NO.7
gsloc(2,7)=1.0d0/sqrt(3.0d0)
gsloc(3,7)=1.0d0/sqrt(3.0d0)
w(7)=1.0d0
gsloc(1,8)=-1.0d0/sqrt(3.0d0)!NO.8
gsloc(2,8)=1.0d0/sqrt(3.0d0)
gsloc(3,8)=1.0d0/sqrt(3.0d0)
w(8)=1.0d0

shl=0.0d0
shl0=0.0d0
do i=1,nint
	shl(1,1,i)=cst*(-1.0d0)*(1.0d0-gsloc(2,i))*(1.0d0-gsloc(3,i))
	shl(2,1,i)=cst*(1.0d0-gsloc(1,i))*(-1.0d0)*(1.0d0-gsloc(3,i))
	shl(3,1,i)=cst*(1.0d0-gsloc(1,i))*(1.0d0-gsloc(2,i))*(-1.0d0)
	shl0(1,i)=cst*(1.0d0-gsloc(1,i))*(1.0d0-gsloc(2,i))*(1.0d0-gsloc(3,i))
!
	shl(1,2,i)=cst*(1.0d0)*(1.0d0-gsloc(2,i))*(1.0d0-gsloc(3,i))
	shl(2,2,i)=cst*(1.0d0+gsloc(1,i))*(-1.0d0)*(1.0d0-gsloc(3,i))
	shl(3,2,i)=cst*(1.0d0+gsloc(1,i))*(1.0d0-gsloc(2,i))*(-1.0d0)
	shl0(2,i)=cst*(1.0d0+gsloc(1,i))*(1.0d0-gsloc(2,i))*(1.0d0-gsloc(3,i))	
!
	shl(1,3,i)=cst*(1.0d0)*(1.0d0+gsloc(2,i))*(1.0d0-gsloc(3,i))
	shl(2,3,i)=cst*(1.0d0+gsloc(1,i))*(1.0d0)*(1.0d0-gsloc(3,i))
	shl(3,3,i)=cst*(1.0d0+gsloc(1,i))*(1.0d0+gsloc(2,i))*(-1.0d0)
	shl0(3,i)=cst*(1.0d0+gsloc(1,i))*(1.0d0+gsloc(2,i))*(1.0d0-gsloc(3,i))
!
	shl(1,4,i)=cst*(-1.0d0)*(1.0d0+gsloc(2,i))*(1.0d0-gsloc(3,i))
	shl(2,4,i)=cst*(1.0d0-gsloc(1,i))*(1.0d0)*(1.0d0-gsloc(3,i))
	shl(3,4,i)=cst*(1.0d0-gsloc(1,i))*(1.0d0+gsloc(2,i))*(-1.0d0)
	shl0(4,i)=cst*(1.0d0-gsloc(1,i))*(1.0d0+gsloc(2,i))*(1.0d0-gsloc(3,i))	
!
	shl(1,5,i)=cst*(-1.0d0)*(1.0d0-gsloc(2,i))*(1.0d0+gsloc(3,i))
	shl(2,5,i)=cst*(1.0d0-gsloc(1,i))*(-1.0d0)*(1.0d0+gsloc(3,i))
	shl(3,5,i)=cst*(1.0d0-gsloc(1,i))*(1.0d0-gsloc(2,i))*(1.0d0)
	shl0(5,i)=cst*(1.0d0-gsloc(1,i))*(1.0d0-gsloc(2,i))*(1.0d0+gsloc(3,i))	
!
	shl(1,6,i)=cst*(1.0d0)*(1.0d0-gsloc(2,i))*(1.0d0+gsloc(3,i))
	shl(2,6,i)=cst*(1.0d0+gsloc(1,i))*(-1.0d0)*(1.0d0+gsloc(3,i))
	shl(3,6,i)=cst*(1.0d0+gsloc(1,i))*(1.0d0-gsloc(2,i))*(1.0d0)
	shl0(6,i)=cst*(1.0d0+gsloc(1,i))*(1.0d0-gsloc(2,i))*(1.0d0+gsloc(3,i))
!
	shl(1,7,i)=cst*(1.0d0)*(1.0d0+gsloc(2,i))*(1.0d0+gsloc(3,i))
	shl(2,7,i)=cst*(1.0d0+gsloc(1,i))*(1.0d0)*(1.0d0+gsloc(3,i))
	shl(3,7,i)=cst*(1.0d0+gsloc(1,i))*(1.0d0+gsloc(2,i))*(1.0d0)
	shl0(7,i)=cst*(1.0d0+gsloc(1,i))*(1.0d0+gsloc(2,i))*(1.0d0+gsloc(3,i))	
!
	shl(1,8,i)=cst*(-1.0d0)*(1.0d0+gsloc(2,i))*(1.0d0+gsloc(3,i))
	shl(2,8,i)=cst*(1.0d0-gsloc(1,i))*(1.0d0)*(1.0d0+gsloc(3,i))
	shl(3,8,i)=cst*(1.0d0-gsloc(1,i))*(1.0d0+gsloc(2,i))*(1.0d0)
	shl0(8,i)=cst*(1.0d0-gsloc(1,i))*(1.0d0+gsloc(2,i))*(1.0d0+gsloc(3,i))
enddo
  
end subroutine shlc8g2
