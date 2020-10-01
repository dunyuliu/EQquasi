MODULE globalvar

implicit none
save 
integer(kind=4),parameter::nen=8,ned=3,nee=24,ndof=3,nint=8
integer(kind=4),parameter::ntotft=1,IS=2
integer(kind=4)::ProType=1
integer(kind=4)::C_Nuclea=1
integer(kind=4)::dis4uniF=4,dis4uniB=4
real(kind=8)::dx=300.d0,rat=1.3d0
real(kind=8)::srcrad0=3.0d3
real(kind=8)::xmin=-30.0d3,xmax=30.0d3,ymin=-30.0d3,ymax=30.0d3,zmin=-30.0d3,zmax=0.0d0
real(kind=8)::xsource=-9.0d3,ysource=0.0d0,zsource=-6.0d3
real(kind=8)::term=2.0d0,dt=0.025d0
integer(kind=4)::friclaw=3
real(kind=8)::fltxyz(2,4,1)
real(kind=8)::rdalfa=0.00d0,rdbeta=0.0d-3
integer(kind=4)::nplpts,nhplt=1,nstep=10000
real(kind=8)::pi=3.14159265358979323846d0,time=0.0d0,critt0=3.0d3,k0=1.0d8,ksi=0.2d0,LL=0.011d0,ttheta=10.0d0,normalcap=-10.0d6, pma, sliprmax
integer(kind=4)::OMPNUM=20,epsi = 9, sliprate_increase_count, sliprate_increase_count_VW
end MODULE globalvar 
