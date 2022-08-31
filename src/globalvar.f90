MODULE globalvar
	implicit none
	!save 
	integer, parameter   :: dp = selected_real_kind(15,307)
	character (len = 30) :: sttmp, dptmp, proc_str
	real (kind = dp)::timebegin, timeover
	!---------------------------------------------------------------!
	!--------------------- system status variables -----------------!
	integer (kind = 4)   :: eqquasi_mode ! 1: quasi-dynamic/quasi-static; 2: fully-dynamic.
	integer (kind = 4)   :: status0, status1, stoptag, icstart
		! status0/status1: inter-seismic (0) or co-seismic (1) for the last/current time step of the code.
		! stoptag: exit (1) or continue running (0) the code.
		! icstart:
	!---------------------------------------------------------------!
	!------------------------- FE system ---------------------------!
	integer (kind = 4) :: numnp, numel, neq, maxa, IERR, me, nprocs, np=1000000, it, itag, ivmax = 100, sol_op, nxt,nyt,nzt
	real (kind = dp), allocatable, dimension(:,:) :: x, mat, fnft,arn,r4nuc,arn4m,slp4fri,globaldat, cons,constmp,consv,consvtmp,consa,consf,consm
	real (kind = dp), allocatable, dimension(:) :: brhs, d, v, mass, kstiff, dump, f, right, resu, resu_1, eledet
	integer (kind = 4), allocatable, dimension(:,:) :: ien, id, anonfs,idhist	
	integer (kind = 4), allocatable, dimension(:) :: et, nftnd, ia, num, ja
	integer (kind = 4), allocatable,dimension(:,:,:) :: nsmp
	real (kind = dp), allocatable,dimension(:,:,:) :: fric,un,us,ud,fltslp, fltsta
	!real (kind = dp), allocatable,dimension(:) :: xlinet, ylinet, zlinet
	integer(kind=4),allocatable::lie(:,:), n4yn(:), an4nds(:,:), nonfs(:)
	integer(kind=4)::nftmx,nonmx,n4onf, n4nds, n4out=0, ndout
	real (kind = dp), allocatable :: x4nds(:,:), xonfs(:,:,:), dout(:,:), rough_geo(:,:)
	real(kind=8)::dtev,dtev1,totmomrate,maxsliprate,totruptarea,tottaoruptarea, totslipruptarea,tdynastart,tdynaend,tcheck,t_end_status, t_start_status	
	integer (kind = 4) :: nen=8,ned=3,nee=24,ndof=3,nint=8, noid=2
	integer (kind = 4) :: C_elastic, C_Nuclea, friclaw, rough_fault, bp, ntotft, nnx, nnz, nres, nplpts, nhplt=20, nstep=1000000, dis4uniF, dis4uniB, C_farfield = 0
	integer (kind = 4) :: nt_output_stress
	real (kind = dp)::srcrad0=3.0d3, rdalfa=0.0d0, rdbeta=0.0d0, fltxyz(2,4,1), init_norm, critt0
	real (kind = dp)::xmin, xmax, ymin, ymax, zmin, zmax, dx, rat, dymax = 3.0d3, dxtmp
	real (kind = dp) :: xminc, xmaxc, zminc ! creeping zone boundaries. Outside of x>xmaxc, x<xminc, and z<zminc, the fault creeps at v = load_slip_rate assigned through input model.txt.
	real (kind = dp) :: rough_fx_min, rough_fx_max, rough_fz_min
	real (kind = dp)::xsource=-9.0d3,ysource=0.0d0,zsource=-6.0d3
	real (kind = dp)::term=200.0d0,dt=0.06d0
	real (kind = dp)::pi=3.14159265358979323846d0,time=0.0d0
	real (kind = dp) :: mat0(1,3), ksi, minDc, slipr_thres, ttheta=10.0d0, far_load_rate, load_slip_rate, pma, sliprmax
	real (kind = dp) :: fric_sw_fs, fric_sw_fd, fric_sw_D0, fric_rsf_a, fric_rsf_b, fric_rsf_Dc, fric_rsf_deltaa0, fric_rsf_vinix, fric_rsf_viniz, fric_rsf_r0, fric_rsf_v0, fric_rsf_fw, fric_rsf_vw, fric_rsf_deltavw0, fric_tp_a_th, fric_tp_rouc, fric_tp_lambda, fric_tp_h, fric_tp_a_hy, fric_tp_deltaa_hy0, fric_ww, fric_w, fric_ini_sliprate, fric_tp_Tini, fric_tp_pini
	! dymax: maximum element size along y/fault-normal direction 
	!Parameters for AZTEC
	integer (kind = 4) ::dimestimate, AZTEC_OPTIONS, minneq,maxneq, N_update
	real (kind = dp), allocatable :: val(:)
	integer (kind = 4), allocatable :: bindx(:), update(:)
	real (kind = dp) :: azmaxiter, aztol
	!AZTEC_OPTIONS == 1: default setup
	!AZTEC_OPTIONS == 2: options(AZ_max_iter) = 2000; options(AZ_output) = AZ_last; params(AZ_tol) = 1.0d-8
	!		params(AZ_drop) = 0.0d0; options(AZ_solver) = AZ_gmres; options(AZ_precond) = AZ_dom_decomp
	!		options(AZ_subdomain_solve) = AZ_ilu; params(AZ_ilut_fill) = 0;
	!AZTEC_OPTIONS == 3: 
end MODULE globalvar 
