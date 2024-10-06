! Copyright (C) 2020 Dunyu Liu <dliu@ig.utexas.edu>, Benchun Duan <bduan@tamu.edu>.
! MIT
MODULE globalvar
    implicit none

    integer, parameter :: dp = selected_real_kind(15,307) ! precision of double precision
    real (kind=dp) :: pi = 3.14159265358979323846d0
    character (len = 30) :: sttmp, dptmp, proc_str 
    
    ! System control parameters
    integer (kind=4) :: eqquasi_mode     ! 1: quasi-dynamic/quasi-static; 2: fully-dynamic.
    integer (kind=4) :: status0, status1 ! status0/status1: inter-seismic (0) or co-seismic (1) for the last/current time step of the code.
    integer (kind=4) :: stoptag          ! stoptag: exit (1) or continue running (0) the code.
    integer (kind=4) :: icstart          ! icstart:
    integer (kind=4) :: C_elastic        ! 1: elastic, default; 2: ?.
    integer (kind=4) :: C_Nuclea         ! 1: artificial nucleation; 0: none, default.   
    integer (kind=4) :: friclaw          ! 3: RSF with aging law, default; 4:?.
    integer (kind=4) :: rough_fault      ! 1: insert customized fault interface; 0: none, default.
    integer (kind=4) :: bp               ! SCEC benchmark problem ID or ID tag for customization. BP5/BP7.
    
    integer (kind=4) :: numnp, numel, neq, maxa, IERR, me, &
        nprocs, np=1000000, it, itag, ivmax = 100, sol_op, &
        nxt, nyt, nzt, nftmx, nonmx, n4onf, n4nds, n4out=0, &
        ndout, nen=8, ned=3, nee=24, ndof=3, nint=8, noid=2, &
        ntotft, nnx, nnz, nres, nplpts, nhplt=20, nstep, &
        dis4uniF, dis4uniB, C_farfield = 0, nt_output_stress
    real (kind=dp), allocatable :: x(:,:), mat(:,:), fnft(:,:), &
        arn(:,:), r4nuc(:,:), arn4m(:,:), slp4fri(:,:), &
        globaldat(:,:), cons(:,:), constmp(:,:), consv(:,:), &
        consvtmp(:,:), consa(:,:), consf(:,:), consm(:,:), &
        brhs(:), d(:), v(:), mass(:), kstiff(:), dump(:), &
        f(:), right(:), resu(:), resu_1(:), eledet(:), &
        fric(:,:,:), un(:,:,:), us(:,:,:), ud(:,:,:), &
        fltslp(:,:,:), fltsta(:,:,:), x4nds(:,:), xonfs(:,:,:), &
        dout(:,:), rough_geo(:,:)
    integer (kind=4), allocatable:: ien(:,:), id(:,:), anonfs(:,:), &
        idhist(:,:), et(:), nftnd(:), ia(:), num(:), ja(:), &
        nsmp(:,:,:), lie(:,:), n4yn(:), an4nds(:,:), nonfs(:)   
    real (kind=dp) :: dtev, dtev1 ! dtev/dtev1: time step defined by ksi*L/Max(slip-rate)
    ! dtev is calculated from the formula and dtev1 is round numbers of dt. 
    real (kind=dp) :: dt ! minimum time step.
    real (kind=dp) :: totMomRate, totMomRateVW, totRuptArea, totTaoRuptArea, &
        totSlipRuptArea, maxSlipRate, &
        tdynastart, tdynaend, tcheck, t_end_status, t_start_status, &
        srcrad0=3.0d3, rdalfa=0.0d0, rdbeta=0.0d0, fltxyz(2,4,1)
    real (kind=dp) :: xmin, xmax, ymin, ymax, zmin, zmax ! left/right/front/back/top/bot domain boundaries.
    real (kind=dp) :: dx ! grid cell sizes.
    real (kind=dp) :: rat ! geometrical enlarging ratio of cell size outside of uniform grid domain.
    ! Currently only along y direction outside of dis4uniF and dis4uniB.
    real (kind=dp) :: dymax, dxtmp
    real (kind=dp) :: xminc, xmaxc, zminc ! On-fault creeping zone left/right/top boundaries. 
    ! Outside of x>xmaxc, x<xminc, and z<zminc, the fault creeps at fixed slip rate of load_slip_rate.
    real (kind=dp) :: rough_fx_min, rough_fx_max, rough_fz_min ! left/right/top boundaries of the inserted fault interface.
    real (kind=dp) :: xsource=-9.0d3,ysource=0.0d0,zsource=-6.0d3 ! hypocenter location for artificial nucleation.
    real (kind=dp) :: term=200.0d0
    real (kind=dp) :: time=0.0d0
    real (kind=dp) :: mat0(1,3), ksi, minDc, slipr_thres, ttheta=10.0d0, far_load_rate, load_slip_rate, pma, sliprmax
    real (kind=dp) :: init_norm ! initial effective normal stress on the fault, MPa.
    real (kind=dp) :: critt0 ! critical time for the nucleation to occur.
    real (kind=dp) :: min_norm, max_norm ! minimum and maximum effective normal stress caps on the fault, MPa. 
    ! parameters for friction laws.
    real (kind=dp) :: fric_sw_fs,       fric_sw_fd,      fric_sw_D0, &
        fric_rsf_a,       fric_rsf_b,      fric_rsf_Dc, &
        fric_rsf_deltaa0, fric_rsf_vinix,  fric_rsf_viniz, &
        fric_rsf_r0,      fric_rsf_v0,     fric_rsf_fw, &
        fric_rsf_vw,      fric_rsf_deltavw0, fric_tp_a_th, &
        fric_tp_rouc,     fric_tp_lambda,  fric_tp_h, &
        fric_tp_a_hy,     fric_tp_deltaa_hy0, fric_ww, &
        fric_w,           fric_ini_sliprate, fric_tp_Tini, &
        fric_tp_pini
    
    ! parameters for nucleation in bp7. 
    real (kind=dp) :: nucx  = -50,      nucz = -50 ! hypocenter x and z coordinates in m. 
    real (kind=dp) :: nucdtao0 = 1.75d6             ! max amplitude of nucleation stress perturbation in Pa. 
    real (kind=dp) :: nucr  = 150.0d0           ! radius of the nucleation stress perturbation 
    real (kind=dp) :: nuct  = 1.0d0             ! duration of the nucleation stress perturbation.
    
    !Parameters for AZTEC
    integer (kind=4) ::dimestimate, AZTEC_OPTIONS, minneq,maxneq, N_update
    real (kind=dp), allocatable :: val(:)
    integer (kind=4), allocatable :: bindx(:), update(:)
    real (kind=dp) :: azmaxiter, aztol
    !AZTEC_OPTIONS == 1: default setup
    !AZTEC_OPTIONS == 2: options(AZ_max_iter) = 2000; options(AZ_output) = AZ_last; params(AZ_tol) = 1.0d-8
    !        params(AZ_drop) = 0.0d0; options(AZ_solver) = AZ_gmres; options(AZ_precond) = AZ_dom_decomp
    !        options(AZ_subdomain_solve) = AZ_ilu; params(AZ_ilut_fill) = 0;
    !AZTEC_OPTIONS == 3: 
end MODULE globalvar 
