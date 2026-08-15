#! /usr/bin/env python3

import numpy as np
from math import *

# this default parameters come from compset bp5.qdc.2000.

# ---------------------------------------------------------------------------
# fric()/on_fault_vars slot registry (rule 5). ONE authoritative list.
#
# The same numbers index three things that must agree: par.on_fault_vars'
# 4th axis here and in the compsets, case.setup's netcdf writer, and
# src/globalvar.f90's Fortran constants of the SAME NAMES. Slots 1-5 are
# friction-law-dependent legacy and stay numeric. Compsets may keep raw
# integers; this table is what they are checked against.
#
#  slot  name               meaning                                [unit] writer
#     6  FR_PORE_DP         pore pressure change                   [Pa]  porepressure.f90
#     7  FR_TNRM0           initial effective normal stress        [Pa]  netcdf_io (slot 8 of input)
#     8  FR_TSTK0           initial shear stress, strike           [Pa]  netcdf_io (slot 7)
#     9  FR_RSF_A           rate-and-state a                       [-]  netcdf_io (slot 1)
#    10  FR_RSF_B           rate-and-state b                       [-]  netcdf_io (slot 2)
#    11  FR_RSF_DC          critical slip distance Dc              [m]  netcdf_io (slot 3)
#    12  FR_RSF_V0          reference slip rate V0                 [m/s]  netcdf_io (slot 4)
#    13  FR_RSF_F0          reference friction f0/r0               [-]  netcdf_io (slot 5)
#    16  FR_VINI_N          initial slip rate, normal              [m/s]  meshgen
#    17  FR_VINI_S          initial slip rate, strike              [m/s]  meshgen
#    18  FR_VINI_D          initial slip rate, dip                 [m/s]  meshgen
#    19  FR_VINI_MAG        initial slip rate magnitude            [m/s]  meshgen
#    20  FR_STATE           state variable theta                   [s]  netcdf_io (slot 9) / faulting
#    22  FR_STATE_TRIAL     trial state update                     [s]  faulting
#    23  FR_THETA_PC        normal-stress state theta_pc           [Pa]  netcdf_io/faulting
#    24  FR_THETA_PC_TRIAL  trial theta_pc                         [Pa]  faulting
#    26  FR_V_TRIAL         trial slip rate magnitude              [m/s]  faulting
#    28  FR_TSTK            shear stress, strike                   [Pa]  faulting
#    29  FR_TDIP            shear stress, dip                      [Pa]  faulting
#    30  FR_TNRM            effective normal stress                [Pa]  faulting
#    31  FR_VXM             master node velocity x                 [m/s]  faulting
#    32  FR_VYM             master node velocity y                 [m/s]  faulting
#    33  FR_VZM             master node velocity z                 [m/s]  faulting
#    34  FR_VXS             slave node velocity x                  [m/s]  faulting
#    35  FR_VYS             slave node velocity y                  [m/s]  faulting
#    36  FR_VZS             slave node velocity z                  [m/s]  faulting
#    41  FR_ABS_KU          abs(KU) elastic force magnitude        [N]  faulting
#    42  FR_V_FINAL         final-step slip rate                   [m/s]  faulting
#    44  FR_SPARE44         spare (written zero)                   [-]  library_output
#    45  FR_SPARE45         spare (written zero)                   [-]  library_output
#    46  FR_VINIT           initial slip rate scalar               [m/s]  netcdf_io (slot 6)
#    47  FR_V_CURRENT       current slip rate (dt control)         [m/s]  faulting/solveTimeLoop
#    48  FR_STATE_INIT      initial state snapshot                 [s]  netcdf_io
#    49  FR_TDIP0           initial shear stress, dip              [Pa]  meshgen
#    51  FR_DARCY_S         Darcy velocity, strike                 [m/s]  porepressure
#    52  FR_DARCY_D         Darcy velocity, dip                    [m/s]  porepressure
#    71  FR_SLIP_S          accumulated slip, strike               [m]  faulting
#    72  FR_SLIP_D          accumulated slip, dip                  [m]  faulting
#    73  FR_SLIP_N          accumulated slip, normal               [m]  faulting
#    74  FR_VS_FINAL        final slip rate, strike                [m/s]  faulting
#    75  FR_VD_FINAL        final slip rate, dip                   [m/s]  faulting
#    76  FR_V_PEAK          running peak slip rate                 [m/s]  faulting
#    81  FR_RUPT_AREA       ruptured area                          [m^2]  faulting
#    82  FR_RUPT_SLIP       total slip in rupture                  [m]  faulting
#    83  FR_TRACT_START     traction at rupture start              [Pa]  faulting
#    84  FR_TRACT_END       traction at rupture end                [Pa]  faulting
# ---------------------------------------------------------------------------
FR_PORE_DP = 6
FR_TNRM0 = 7
FR_TSTK0 = 8
FR_RSF_A = 9
FR_RSF_B = 10
FR_RSF_DC = 11
FR_RSF_V0 = 12
FR_RSF_F0 = 13
FR_VINI_N = 16
FR_VINI_S = 17
FR_VINI_D = 18
FR_VINI_MAG = 19
FR_STATE = 20
FR_STATE_TRIAL = 22
FR_THETA_PC = 23
FR_THETA_PC_TRIAL = 24
FR_V_TRIAL = 26
FR_TSTK = 28
FR_TDIP = 29
FR_TNRM = 30
FR_VXM = 31
FR_VYM = 32
FR_VZM = 33
FR_VXS = 34
FR_VYS = 35
FR_VZS = 36
FR_ABS_KU = 41
FR_V_FINAL = 42
FR_SPARE44 = 44
FR_SPARE45 = 45
FR_VINIT = 46
FR_V_CURRENT = 47
FR_STATE_INIT = 48
FR_TDIP0 = 49
FR_DARCY_S = 51
FR_DARCY_D = 52
FR_SLIP_S = 71
FR_SLIP_D = 72
FR_SLIP_N = 73
FR_VS_FINAL = 74
FR_VD_FINAL = 75
FR_V_PEAK = 76
FR_RUPT_AREA = 81
FR_RUPT_SLIP = 82
FR_TRACT_START = 83
FR_TRACT_END = 84

class parameters:
    # cylce id. Simulate quasi-dynamic earthquake cycles from istart to iend.
    istart = 1
    iend = 1
    # mode of the code - quasi-dynamic (1) or fully-dynamic (2). 
    mode = 1

    dip = 90.

    # model_domain (in meters)
    fxmin, fxmax = -60.0e3, 60.0e3
    fymin, fymax = -50.0e3, 50.0e3
    fzmin, fzmax = -60.0e3, 0.0e3

    # creeping zone bounaries.
    # creeping zones are assinged on the lateral sides and bottom of 
    # the RSF controlled region and will slide at fixed loading slip rate.
    xminc, xmaxc, zminc = -50.0e3, 50.0e3, -40.0e3

    dx = 2000.0e0 # cell size, spatial resolution
    dy = dx
    dz = dx
    nuni_y_plus, nuni_y_minus = 5, 5 # along the fault-normal dimension, the number of cells share the dx cell size.
    enlarging_ratio = 1.3e0 # along the fault-normal dimension (y), cell size will be enlarged at this ratio compoundly.

    # Isotropic material propterty.
    # Vp, Vs, Rou
    vp, vs, rou = 6.0e3, 3.464e3, 2.67e3
    init_norm = -25.0e6 # initial normal stress in Pa. Negative compressive.

    # Controlling switches for EQquasi system
    insertFaultType = 0 # 0: no fault, 1: planar 2: rough.
    rough_fault = insertFaultType
    rheology    = 1 # elastic(1). 
    friclaw     = 3 # rsf_aging(3), rsf_slip(4).
    ntotft      = 1 # number of total faults.
    # Optional per-fault geometry override, one 5-tuple per fault:
    #   (xlo, xhi, ycoor, zlo, zhi) in meters.
    # None (default) means every fault uses the single domain-derived box
    # (ycoor = 0) that case.setup has always emitted; set this to a list of
    # len(ntotft) tuples to place faults at different y (e.g. a step-over).
    faultgeom   = None
    solver      = 1 # solver option. MUMPS(1, recommended). AZTEC(2).
    nstep       = 10000 # total num of time steps for exiting, if not exit via sliprate threshold
    nt_out      = 100 # Every nt_out time steps, disp of the whole model and on-fault variables will be written out in netCDF format.
    bp          = 5 
    # currently supported cases
    # 5 (SCEC-BP5)
    # 1001 (GM-cycle)

    # xi, minimum Dc
    # Lapusta et al. (2009) time-step factor, dtev = xi*Dc/Vmax. Measured for
    # BP8 at dx = 50 m: 0.05, 0.1 and 0.2 cost the same wall clock for a given
    # step count but reach 6.0, 17.0 and 22.3 days, agreeing to 0.003 log units
    # in peak slip rate. 0.2 is the default; tighten it per compset if a problem
    # needs it.
    xi = 0.2
    minDc = 0.13 # meters

    # loading 
    far_vel_load = 4e-10 # far field shear loading velocity (+x vel=far_vel_load on xz bound at y=ymax). A minus value is applied on the other side.
    far_norm_load_vel = 0.0 # far field fault-normal extensional loading velocity (+y vel=far_norm_load_vel on xz bound at y=ymax). 
    creep_slip_rate = 1.0e-9 # creeping slip rate outside of RSF controlled region.
    exit_slip_rate = 1.0e-3 # exiting slip rate for EQquasi [m/s].

    #################################
    ##### Frictional variables ######
    #################################
    # friclaw == 1, slip weakening
    fric_sw_fs = 0
    fric_sw_fd = 0
    fric_sw_D0 = 0
    # friclaw == 3, rate- and state- friction with aging law.
    fric_rsf_a, fric_rsf_b, fric_rsf_Dc = 0.004, 0.03, 0.14
    fric_rsf_deltaa = 0.036
    fric_rsf_r0 = 0.6
    fric_rsf_v0 = 1e-6
    # Creating the fault interface
    nfx = int((fxmax - fxmin)/dx + 1)
    nfz = int((fzmax - fzmin)/dz + 1)
    fx = np.linspace(fxmin,fxmax,nfx) # coordinates of fault grids along strike.
    fz = np.linspace(fzmin,fzmax,nfz) # coordinates of fault grids along dip.
    # Create on_fault_vars array for on_fault varialbes.
    # Always (ntotft, nfz, nfx, 100), with ntotft = 1 the degenerate case --
    # the same neutrality the Fortran requires (rule 13). A single-fault default
    # forced every multi-fault compset to allocate and fill its own 4-D array by
    # hand, which is two schemes for one thing (rule 19).
    #
    # For ntotft > 1 the strike/dip extents come from faultgeom, and the array
    # is sized by the LARGEST fault patch so every fault fits: the owner's
    # ruling, "keep the same nx nz ntotft structure, just arrays sized by the
    # largest fault patch".
    on_fault_vars = np.zeros((ntotft, nfz, nfx, 100))
    def shear_steady_state(a,b,v0,r0,load_rate,norm,slip_rate, rou, vs):
      # calculate shear stress at steady state
      res = -norm*a*asinh(slip_rate/2.0/v0*exp((r0+b*log(v0/load_rate))/a)) + rou*vs/2.0*slip_rate
      return res
      
    for ix, xcoor in enumerate(fx):
      for iz, zcoor in enumerate(fz):
      # assign a in RSF. a is a 2D distribution.
        if abs(zcoor)>=18e3 or abs(xcoor)>=32e3 or abs(zcoor)<=2e3: 
          on_fault_vars[0,iz,ix,FR_RSF_A] = fric_rsf_a + fric_rsf_deltaa
        elif abs(zcoor)<=16e3 and abs(zcoor)>=4e3 and abs(xcoor)<=30e3:
          on_fault_vars[0,iz,ix,FR_RSF_A] = fric_rsf_a
        else:
          tmp1 = (abs(abs(zcoor)-10e3) - 6e3)/2e3
          tmp2 = (abs(xcoor)-30e3)/2e3
          on_fault_vars[0,iz,ix,FR_RSF_A] = fric_rsf_a + max(tmp1,tmp2)*fric_rsf_deltaa
        on_fault_vars[0,iz,ix,FR_RSF_B] = fric_rsf_b # assign b in RSF 
        on_fault_vars[0,iz,ix,FR_RSF_DC] = fric_rsf_Dc # assign Dc in RSF.
        if (xcoor<=-18e3 and xcoor>=-30e3 and zcoor<=-4e3 and zcoor>=-16e3):
          on_fault_vars[0,iz,ix,FR_RSF_DC] = minDc # a special Dc zone.
        on_fault_vars[0,iz,ix,FR_RSF_V0] = fric_rsf_v0 # initial reference slip rate.
        on_fault_vars[0,iz,ix,FR_RSF_F0] = fric_rsf_r0 # initial reference friction.
        
        on_fault_vars[0,iz,ix,FR_VINIT] = creep_slip_rate # initial slip rates
        if (xcoor<=-18e3 and xcoor>=-30e3 and zcoor<=-4e3 and zcoor>=-16e3):
          on_fault_vars[0,iz,ix,FR_VINIT] = 0.03 # initial high slip rate patch.
        on_fault_vars[0,iz,ix,FR_STATE] = on_fault_vars[0,iz,ix,FR_RSF_DC]/creep_slip_rate # initial state var.
        on_fault_vars[0,iz,ix,FR_TNRM0] = init_norm # initial normal stress.
        on_fault_vars[0,iz,ix,FR_TSTK0] = shear_steady_state(on_fault_vars[0,iz,ix,FR_RSF_A], 
                                                    on_fault_vars[0,iz,ix,FR_RSF_B],
                                                    on_fault_vars[0,iz,ix,FR_RSF_V0],
                                                    on_fault_vars[0,iz,ix,FR_RSF_F0],
                                                    creep_slip_rate,
                                                    on_fault_vars[0,iz,ix,FR_TNRM0],
                                                    on_fault_vars[0,iz,ix,FR_VINIT],
                                                    rou,
                                                    vs)
    ###############################################
    ##### Domain boundaries for transferring ######
    ###############################################
    xmin_trans, xmax_trans = -25e3, 25e3
    zmin_trans = -25e3
    ymin_trans, ymax_trans = -5e3, 5e3
    dx_trans = 50 

    ###############################################################
    ##### Along-fault pore fluid diffusion, used when bp == 8 ######
    ###############################################################
    # Off by default, so every other compset is unaffected.
    fluid_src    = 0       # 0: off; 1: Gaussian source (GS); 2: Peaceman well (PW).
    fluid_q0     = 0.0     # total volume injection rate, m^3/s.
    fluid_toff   = 0.0     # injection turn-off time, s.
    fluid_tend   = 0.0     # final simulation time, s. 0 disables the time-based exit.
    fluid_Lgauss = 50.0    # characteristic size of the Gaussian source, m.
    fluid_Lfwid  = 1.0     # fault zone thickness, m.
    fluid_beta   = 1.0e-8  # pore and fluid compressibility, 1/Pa.
    fluid_phi    = 0.1     # porosity.
    fluid_perm   = 5.0e-14 # permeability, m^2.
    fluid_eta    = 1.0e-3  # fluid viscosity, Pa s.
    fluid_Swell  = 1.0e-7  # volumetric well storage, m^3/Pa.
    fluid_rwell  = 0.05    # true well radius, m. BP8 Table 1 (2026-08-12 revision).
    dtmax        = 0.0     # cap on the adaptive time step, s. 0 means no cap.
    # Prakash-Clifton relaxation distance for the normal-stress state variable, m.
    #   0.0 : default. No state; use the instantaneous effective normal stress.
    #  > 0  : relax the normal-stress state over this slip distance.
    fric_pc_L    = 0.0
    # Effective normal stress caps, Pa (negative = compressive). Applied only
    # when the fault is non-planar (rough_fault = 1) AND the material is
    # elastic (C_elastic = 1): on a bent or rough fault, normal stress can run
    # away at geometric asperities, and in the real Earth off-fault inelastic
    # deformation limits it. Previously these were the literals -10e6/-40e6 in
    # src/faulting.f90, invisible to every compset that they silently acted on.
    # Liu, Duan & Luo (2020) section 3.5 uses -10e6/-100e6.
    min_norm     = -10.0e6
    max_norm     = -40.0e6
    # ONE switch: min_norm/max_norm apply iff this is 1. Not coupled to
    # rough_fault. Kink/rough compsets set 1 (they historically had caps);
    # planar cases may set 1 to bound tip stress (BP1002 variants), the
    # regularisation Liu et al. 2020 section 3.5 applies at a bend.
    enforce_norm_caps = 0

    ####################################
    ##### HPC resource allocation ######
    ####################################
    casename = "bp5-qd-2000"
    HPC_nnode = 1 # Number of computing nodes. On LS6, one node has 128 CPUs.
    HPC_ncpu = 3 # Number of CPUs requested.
    HPC_queue = "normal" # q status. Depending on systems, job WALLTIME and Node requested.
    HPC_time = "00:30:00" # WALLTIME, in hh:mm:ss format.
    HPC_account = "EAR22013" # Project account to be charged SUs against.
    HPC_email = "dliu@ig.utexas.edu" # Email to receive job status.

    ##############################################
    ##### Single station time series output ######
    ##############################################

    # (x,z) coordinate pairs for on-fault stations (in km).
    st_coor_on_fault = [[-36.0, 0.0], [-16.0,0.0], [0.0,0.0], [16.0,0.0], \
       [36.0,0.0], [-24.0,0.0], [-16.0,0.0], [0.0,-10.0], [16.0,-10.0], [0.0,-22.0]]
       
    # (x,y,z) coordinates for off-fault stations (in km).
    st_coor_off_fault = [[0,8,0], [0,8,-10], [0,16,0], [0,16,-10], [0,32,0], \
       [0,32,-10], [0,48,0], [16,8,0], [-16,8,0]]
    n_on_fault = len(st_coor_on_fault)
    n_off_fault = len(st_coor_off_fault)

    # Additional solver options for AZTEC
    az_op = 2 # AZTEC options
    az_maxiter = 2000 # maximum iteration for AZTEC
    az_tol = 1.0e-7 # tolerance for solution in AZTEC.


