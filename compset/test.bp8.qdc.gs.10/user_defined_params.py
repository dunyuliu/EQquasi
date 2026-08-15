#! /usr/bin/env python3
# compset: test.bp8.qdc.gs.10
# gate   : e2e fast
# changed: 2026-08-12   (content; see git log for the full history)
# note   : Regression variant. Not for science.
# Status is mirrored in compset/README.md, which a contract test checks.
# ----------------------------------------------------------------------

# SCEC SEAS benchmark BP8-QD-GS (aging law), 50 m, +/-500 m cube.
#
# Regression smoke test: fast and deterministic, 8000 elements. Not a physics
# result -- 50 m cells put the Gaussian source L_inj = 50 m on a single cell,
# and the section 4.3 profiles carry 17 columns where the description asks for
# 81 at 10 m spacing. For BP8 physics use bp8.qdc.gs.10.
# Quasi-dynamic fluid injection in 3D: a planar fault in a homogeneous whole
# space, velocity-strengthening rate-and-state friction everywhere, driven
# purely by the pore pressure change from an injection at the fault centre.
# Parameter values follow Table 1 of the BP8 benchmark description (2026/07/31).

from defaultParameters import parameters
import numpy as np
from math import *

par = parameters()
# BP8 is aseismic by construction, so it is a single run rather than a cycle.
par.istart = 1
par.iend = 1
# mode of the code - quasi-dynamic (1) or fully-dynamic (2).
par.mode = 1

# model_domain (in meters)
par.fxmin, par.fxmax = -500.0, 500.0
par.fymin, par.fymax = -500.0, 500.0
par.fzmin, par.fzmax = -500.0, 500.0

# Frictional / fluid domain Omega_f. Half-length l_f = 400 m. Outside of it the
# fault is held at par.creep_slip_rate (the benchmark's V_zero, i.e. no slip)
# and the pore pressure solver imposes zero fluid flux.
par.xminc, par.xmaxc, par.zminc = -400.0, 400.0, -400.0

par.dx = 50.0e0 # cell size, spatial resolution. Coarse, for a fast smoke test.
par.dy = par.dx
par.dz = par.dx
par.nuni_y_plus, par.nuni_y_minus = 5, 5 # along the fault-normal dimension, the number of cells share the dx cell size.
par.enlarging_ratio = 1.0e0 # along the fault-normal dimension (y), cell size will be enlarged at this ratio compoundly.

# Isotropic material propterty.
# Vp, Vs, Rou. mu = rou*vs^2 = 32.04 GPa, matching Table 1.
par.vp, par.vs, par.rou = 6.0e3, 3.464e3, 2.67e3
par.init_norm = -25.0e6 # initial effective normal stress in Pa. Negative compressive.

# Controlling switches for EQquasi system
par.insertFaultType = 0
par.rough_fault = par.insertFaultType
par.rheology    = 1 # elastic(1).
par.friclaw     = 3 # rsf_aging(3) -> BP8-QD-A; rsf_slip(4) -> BP8-QD-S.
par.ntotft      = 1 # number of total faults.
par.solver      = 1 # solver option. MUMPS(1, recommended). AZTEC(2).
par.nstep       = 200 # short run for the regression test.
par.nt_out      = 100
par.bp          = 8
# currently supported cases
# 5 (SCEC-BP5)
# 7 (SCEC-BP7)
# 8 (SCEC-BP8)
# 1001 (GM-cycle)

# xi, minimum Dc
par.xi = 0.2 # Lapusta et al. (2009) time-step factor, dtev = xi*Dc/Vmax.
# Measured at dx = 50 m: xi = 0.05, 0.1 and 0.2 cost the same wall clock for
# 4000 steps but reach 6.0, 17.0 and 22.3 days of simulated time, and agree to
# 0.003 log units in peak slip rate, 0.3% in moment rate and 0.05% in slip.
par.minDc = 0.5e-3 # meters

# loading. BP8 has no tectonic loading; slip is driven only by fluid injection.
par.far_vel_load = 0.0
par.far_norm_load_vel = 0.0
par.creep_slip_rate = 1.0e-20 # V_zero outside Omega_f, i.e. no slip.
par.exit_slip_rate = 1.0e0 # bp8 exits on fluid_tend, not on slip rate.

#################################
##### Frictional variables ######
#################################
# friclaw == 1, slip weakening
par.fric_sw_fs = 0
par.fric_sw_fd = 0
par.fric_sw_D0 = 0
# friclaw == 3, rate- and state- friction with aging law.
# a > b everywhere: velocity strengthening, so the response is entirely aseismic.
par.fric_rsf_a, par.fric_rsf_b, par.fric_rsf_Dc = 0.016, 0.010, 0.5e-3
par.fric_rsf_r0 = 0.6
par.fric_rsf_v0 = 1e-6
par.init_slip_rate = 1e-12 # V_init
par.init_shear = 14.6e6 # Table 1 tau_init.

#####################################
##### Fluid injection (bp == 8) #####
#####################################
par.fluid_src    = 1        # 1: Gaussian source (GS); 2: Peaceman well (PW).
par.fluid_q0     = 0.003    # Q0, total volume injection rate, m^3/s.
par.fluid_toff   = 100*3600.0       # t_off, injection turn-off time, s.
par.fluid_tend   = 30*24*3600.0     # t_f, final simulation time, s.
par.fluid_Lgauss = 50.0     # L_inj, characteristic distance of the Gaussian source, m.
par.fluid_Lfwid  = 1.0      # L_fwid, fault thickness, m.
par.fluid_beta   = 1.0e-8   # beta, pore and fluid compressibility, 1/Pa.
par.fluid_phi    = 0.1      # phi, porosity.
par.fluid_perm   = 5.0e-14  # k, permeability, m^2.
par.fluid_eta    = 1.0e-3   # eta, fluid viscosity, Pa s.
par.fluid_Swell  = 1.0e-7   # S_well, volumetric well storage, m^3/Pa.
par.fluid_rwell  = 0.0508   # r_well, well radius, m. Not fixed by Table 1; 2 inches.
# Cap the adaptive time step. Without a cap, xi*Dc/V is ~1e7 s at V_init.
# 500 s also satisfies the explicit diffusion limit dx^2/(4*alpha) at dx = 10 m.
par.dtmax        = 500.0
# BP8 eq. (10) builds strength from the instantaneous effective normal stress;
# there is no normal-stress state variable. 0 selects that limit.
par.fric_pc_L    = 0.0

# Creating the fault interface
par.nfx = round((par.fxmax - par.fxmin)/par.dx + 1)
par.nfz = round((par.fzmax - par.fzmin)/par.dz + 1)
par.fx = np.linspace(par.fxmin,par.fxmax,par.nfx) # coordinates of fault grids along strike.
par.fz = np.linspace(par.fzmin,par.fzmax,par.nfz) # coordinates of fault grids along dip.
# Sized here, not in defaultParameters: the compset redefines nfx/nfz above,
# so the array must be allocated once the grid is known. Always 4-D --
# (ntotft, nfz, nfx, 100) -- with ntotft = 1 the degenerate case.
par.on_fault_vars = np.zeros((par.ntotft, par.nfz, par.nfx, 100))
# Create on_fault_vars array for on_fault varialbes.
def state_from_stress_and_rate(a, b, Dc, v0, r0, norm, tau, slip_rate):
  # BP8 reading C. Table 1 gives tau_init = 14.6 MPa and eq. (27) requires the
  # fault to start at a uniform V_init; eq. (12) then fixes theta_0, so eq. (29)
  # (theta_0 = D_RS/V_init) cannot also hold. The specification is
  # over-determined -- see README, "The initial condition is over-determined".
  #
  # This reading keeps Table 1 and eq. (27) and derives theta_0 = 4.0188e11 s.
  # It is the only one consistent with both things we can verify about
  # taehoKim_ref: its shear_stress_2 starts at 14.6 MPa, and it matches our
  # slip at both compared stations to within 2 % and 7 %.
  psi = a*log(2.0*v0/slip_rate*sinh(tau/(abs(norm)*a)))
  return (Dc/v0)*exp((psi - r0)/b)

for ix, xcoor in enumerate(par.fx):
  for iz, zcoor in enumerate(par.fz):
    par.on_fault_vars[0,iz,ix,9]  = par.fric_rsf_a  # a in RSF, uniform.
    par.on_fault_vars[0,iz,ix,10] = par.fric_rsf_b  # b in RSF, uniform.
    par.on_fault_vars[0,iz,ix,11] = par.fric_rsf_Dc # Dc in RSF.
    par.on_fault_vars[0,iz,ix,12] = par.fric_rsf_v0 # reference slip rate V*.
    par.on_fault_vars[0,iz,ix,13] = par.fric_rsf_r0 # reference friction f*.

    par.on_fault_vars[0,iz,ix,46] = par.init_slip_rate # initial slip rate V_init.
    # BP8 eq. (29): initial state at steady state with V_init. Prescribed by
    # the benchmark. See bp8.qdc.gs.10 for the tau_init note.
    par.on_fault_vars[0,iz,ix,20] = state_from_stress_and_rate(
        par.fric_rsf_a, par.fric_rsf_b, par.fric_rsf_Dc, par.fric_rsf_v0,
        par.fric_rsf_r0, par.init_norm, par.init_shear, par.init_slip_rate)
    par.on_fault_vars[0,iz,ix,7]  = par.init_norm  # initial effective normal stress.
    par.on_fault_vars[0,iz,ix,8]  = par.init_shear # tau^0, Table 1. # tau^0, derived. See above.

####################################
##### HPC resource allocation ######
####################################
par.casename = "test-bp8-qdc"
par.HPC_nnode = 1 # Number of computing nodes. On LS6, one node has 128 CPUs.
par.HPC_ncpu = 20 # Number of CPUs requested.
par.HPC_queue = "normal" # q status. Depending on systems, job WALLTIME and Node requested.
par.HPC_time = "24:00:00" # WALLTIME, in hh:mm:ss format.
par.HPC_account = "EAR22013" # Project account to be charged SUs against.
par.HPC_email = "dliu@ig.utexas.edu" # Email to receive job status.

##############################################
##### Single station time series output ######
##############################################
# BP8 asks for 9 on-fault stations at x2, x3 in {-200, 0, 200} m.
# Coordinates here are (x, z) in km; BP8's x3 is positive downward, so x3 = -z.
par.st_coor_on_fault = [[0,0], [-200,0], [0,-200], [200,0], [0,200], \
   [-200,200], [-200,-200], [200,200], [200,-200]]
par.st_coor_on_fault = (np.asarray(par.st_coor_on_fault)/1000.0).tolist()

# (x,y,z) coordinates for off-fault stations (in km).
par.st_coor_off_fault = [[0,200,0], [0,400,0]]
par.st_coor_off_fault = (np.asarray(par.st_coor_off_fault)/1000.0).tolist()
par.n_on_fault = len(par.st_coor_on_fault)
par.n_off_fault = len(par.st_coor_off_fault)

# Additional solver options for AZTEC
par.az_op = 2 # AZTEC options
par.az_maxiter = 2000 # maximum iteration for AZTEC
par.az_tol = 1.0e-7 # tolerance for solution in AZTEC.
