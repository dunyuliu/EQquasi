MODULE globalvar

implicit none
save 
integer(kind=4),parameter::nen=8,ned=3,nee=24,ndof=3,nint=8
integer(kind=4),parameter::ntotft=1,IS=2
integer(kind=4)::ProType=1
integer(kind=4)::C_Nuclea=1
integer(kind=4)::dis4uniF=4,dis4uniB=4
real(kind=8)::dx=1000.d0,rat=1.3d0
real(kind=8)::srcrad0=3.0d3
real(kind=8)::xmin=-60.0d3,xmax=60.0d3,ymin=-100.0d3,ymax=100.0d3,zmin=-120.0d3,zmax=0.0d0
real(kind=8)::xsource=-9.0d3,ysource=0.0d0,zsource=-6.0d3
real(kind=8)::term=200.0d0,dt=0.1d0
integer(kind=4)::friclaw=3
real(kind=8)::fltxyz(2,4,1)
real(kind=8)::rdalfa=0.00d0,rdbeta=0.0d-3
integer(kind=4)::nplpts,nhplt=1,nstep=1000000
real(kind=8)::pi=3.14159265358979323846d0,time=0.0d0,critt0=3.0d3,k0=1.0d8
real(kind = 8) :: ksi=0.2d0, LL=0.04d0, ttheta=10.0d0, loadrate = 0.0d-10, pma, sliprmax
integer(kind=4)::OMPNUM=20,epsi = 1
end MODULE globalvar 
