MODULE globalvar

	implicit none
	save 
	
	integer, parameter :: dp = selected_real_kind(15,307)
	character (len = 30) :: sttmp, dptmp
	!---------------------------------------------------------------!
	!------------------------- FE system ---------------------------!
	real (kind = dp)::timebegin, timeover
	integer(kind=4)::numnp, numel, neq, maxa, IERR, me, nprocs, np=1000000, status0, status1, it, itag, ivmax = 100, sol_op
	real (kind = dp), allocatable, dimension(:,:) :: x, mat, fnft,arn,r4nuc,arn4m,slp4fri,globaldat, &
										cons,constmp,consv,consvtmp,consa,consf,consm
	real (kind = dp), allocatable, dimension(:) :: brhs, d, v, mass, kstiff, dump, f, right, resu, resu_1, eledet
	integer(kind=4), allocatable, dimension(:,:) :: ien, id, anonfs,idhist	
	integer(kind=4), allocatable, dimension(:) :: et, nftnd, ia, num, ja
	integer(kind=4), allocatable,dimension(:,:,:) :: nsmp
	real (kind = dp), allocatable,dimension(:,:,:) :: fric,un,us,ud,fltslp, fltsta
	integer(kind=4),allocatable::lie(:,:), n4yn(:), an4nds(:,:), nonfs(:)
	integer(kind=4)::nftmx,nonmx,n4onf, n4nds, n4out=0, ndout
	real (kind = dp), allocatable :: x4nds(:,:), xonfs(:,:,:), dout(:,:)
	real(kind=8)::dtev,dtev1,totmomrate,maxsliprate,totruptarea,tottaoruptarea, totslipruptarea,tdynastart,tdynaend,tcheck,t_end_status, t_start_status	
	integer(kind=4),parameter::nen=8,ned=3,nee=24,ndof=3,nint=8, noid=2
	integer(kind=4),parameter::ntotft=1,IS=2
	integer(kind=4)::ProType=1
	integer(kind=4)::C_Nuclea=1
	integer(kind=4)::dis4uniF=5,dis4uniB=5
	real (kind = dp)::dx,rat=1.3d0
	real (kind = dp)::srcrad0=3.0d3
	real (kind = dp)::xmin=-60.0d3,xmax=60.0d3,ymin=-50.0d3,ymax=50.0d3,zmin=-60.0d3,zmax=0.0d3
	real (kind = dp)::xsource=-9.0d3,ysource=0.0d0,zsource=-6.0d3
	real (kind = dp)::term=200.0d0,dt=0.06d0
	integer(kind=4)::friclaw=3
	real (kind = dp)::fltxyz(2,4,1)
	real (kind = dp)::rdalfa=0.00d0,rdbeta=0.0d-3, slipr_thres = 1.0d-3, dymax = 3.0d3
	integer(kind=4)::nplpts, nhplt=20, nstep=1000000
	real (kind = dp)::pi=3.14159265358979323846d0,time=0.0d0,critt0=3.0d3,k0=1.0d8
	real(kind = 8) :: ksi=0.015d0, LL=0.13d0, ttheta=10.0d0, loadrate = 5.0d-10, pma, sliprmax
	integer(kind=4)::OMPNUM=20,epsi = 1, C_farfield = 0, icstart, stoptag
	
	!Parameters for pardiso
	integer(kind=8)::pt(64)
	integer(kind=4)::maxfct,mnum,mtype,phase,nrhs,error,msglvl,solver
	integer(kind=4)::iparm(64),idum(1)
	REAL*8  dparm(64)
	real*8 ddum(1)
	data nrhs/1/,maxfct/1/,mnum/1/	
	
	!Parameters for AZTEC
	
	integer(kind=4)::dimestimate, AZTEC_OPTIONS, minneq,maxneq, N_update
	real (kind = dp), allocatable :: val(:)
	integer (kind = 4), allocatable :: bindx(:), update(:)
	real (kind = dp) :: azmaxiter, aztol
	!AZTEC_OPTIONS == 1: default setup
	!AZTEC_OPTIONS == 2: options(AZ_max_iter) = 2000; options(AZ_output) = AZ_last; params(AZ_tol) = 1.0d-8
	!		params(AZ_drop) = 0.0d0; options(AZ_solver) = AZ_gmres; options(AZ_precond) = AZ_dom_decomp
	!		options(AZ_subdomain_solve) = AZ_ilu; params(AZ_ilut_fill) = 0;
	!AZTEC_OPTIONS == 3: 
end MODULE globalvar 
