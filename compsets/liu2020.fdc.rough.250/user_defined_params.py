#! /usr/bin/env python3
# compset: liu2020.fdc.rough.250
# gate   : none -- ported to par schema 2026-08-15, UNVERIFIED (no reference)
# changed: 2026-08-12   (content; see git log for the full history)
# note   : Ported from the pre-par bare-variable format. No oracle exists;
#          runs but nothing checks the numbers.
# Status is mirrored in compsets/README.md, which a contract test checks.
# ----------------------------------------------------------------------

import numpy as np
from math import *
from defaultParameters import parameters
par = parameters()
from lib     import *

# cylce id. Simulate quasi-dynamic earthquake cycles from par.istart to par.iend.
par.istart     = 1  
par.iend       = 1
# par.mode of the code - quasi-dynamic (1) or fully-dynamic (2). 
par.mode       = 2

# model_domain (in meters)
par.fxmin, par.fxmax = -40e3, 40e3
par.fymin, par.fymax = -10e3, 10e3
par.fzmin, par.fzmax = -40e3, 0

# creeping zone bounaries.
# creeping zones are assinged on the lateral sides and bottom of 
# the RSF controlled region and will slide at fixed loading slip rate.
par.xminc, par.xmaxc, par.zminc = -35e3, 35e3, -35e3

par.dx = 250.0e0 # cell size, spatial resolution
par.dy = par.dx
par.dz = par.dx
par.nuni_y_plus, par.nuni_y_minus = 5, 5 # along the fault-normal dimension, the number of cells share the par.dx cell size.
par.enlarging_ratio = 1.3 # along the fault-normal dimension (y), cell size will be enlarged at this ratio compoundly.

# Isotropic material propterty.
# Vp, Vs, Rou
par.vp, par.vs, par.rou = 6.0e3, 3.464e3, 2.67e3
par.init_norm = -25.0e6 # initial normal stress in Pa. Negative compressive.

# Controlling switches for EQquasi system
par.rough_fault = 1 # include rough fault yes(1) or not(0).
par.rheology    = 1 # elastic(1). 
par.friclaw     = 3 # rsf_aging(3), rsf_slip(4).
par.ntotft      = 1 # number of total faults.
par.solver      = 1 # par.solver option. MUMPS(1, recommended). AZTEC(2).
par.nt_out      = 100 # Every par.nt_out time steps, disp of the whole model and on-fault variables will be written out in netCDF format.
par.bp          = 1000 # for description only
# currently supported cases
# 5 (SCEC-BP5)
# 7 (SCEC-BP7)
# 1000 (liu2020.fdc.planar, GM-cycle)

# par.xi, minimum Dc
par.xi              = 0.005 # 0.02 for planar; par.xi used to limit variable time step size. See Lapusta et al. (2009).
# par.xi=0.2/0.05 all show diverging results.
par.minDc           = 11e-3 # meters

# loading 
par.far_vel_load    = 5e-10 # far field loading velocity on xz planes. A minus value is applied on the other side.
par.creep_slip_rate = 1.0e-9 # creeping slip rate outside of RSF controlled region.
par.exit_slip_rate  = 1.5e-2 # exiting slip rate for EQquasi [m/s].
# setting par.exit_slip_rate to be 1.0e-2 m/s. Dynamic rupture slip rates may drop to the level 1e-2~1e-3 m/s after 100 seconds. It would be inefficient to run dynamic ruptures longer than a few hundreds seconds.

#################################
##### Frictional variables ######
#################################
# par.friclaw == 1, slip weakening
par.fric_sw_fs = 0
par.fric_sw_fd = 0
par.fric_sw_D0 = 0
# par.friclaw == 3, rate- and state- friction with aging law.
par.fric_rsf_a      = 0.007 
par.fric_rsf_b      = 0.011 
par.fric_rsf_Dc     = 11e-3
par.fric_rsf_deltaa = 0.01
par.fric_rsf_r0     = 0.6
par.fric_rsf_v0     = 1e-6
# Creating the fault interface
par.nfx = int((par.fxmax - par.fxmin)/par.dx + 1)
par.nfz = int((par.fzmax - par.fzmin)/par.dx + 1)
par.fx = np.linspace(par.fxmin,par.fxmax,par.nfx) # coordinates of fault grids along strike.
par.fz = np.linspace(par.fzmin,par.fzmax,par.nfz) # coordinates of fault grids along par.dip.
# Create par.on_fault_vars array for on_fault varialbes.
par.on_fault_vars = np.zeros((par.ntotft, par.nfz, par.nfx, 100))
def shear_steady_state(a,b,v0,r0,load_rate,norm,slip_rate):
  # calculate shear stress at steady state
  res = -norm*a*asinh(slip_rate/2.0/v0*exp((r0+b*log(v0/load_rate))/a)) + par.rou*par.vs/2.0*slip_rate
  return res

for ix, par.xcoor in enumerate(par.fx):
  for iz, par.zcoor in enumerate(par.fz):
  # assign a in RSF. a is a 2D distribution.
    tmp1  = linear1(par.xcoor, 20.e3, 3.e3)
    tmp2  = linear1(-par.zcoor-10.e3, 7.e3, 1.e3)
    par.on_fault_vars[0,iz,ix,9]  = par.fric_rsf_a + (1. - tmp1*tmp2)*par.fric_rsf_deltaa
    par.on_fault_vars[0,iz,ix,10] = par.fric_rsf_b # assign b in RSF 
    par.on_fault_vars[0,iz,ix,11] = par.fric_rsf_Dc # assign Dc in RSF.
    par.on_fault_vars[0,iz,ix,12] = par.fric_rsf_v0 # initial reference slip rate.
    par.on_fault_vars[0,iz,ix,13] = par.fric_rsf_r0 # initial reference friction.
    par.on_fault_vars[0,iz,ix,46] = par.creep_slip_rate # initial slip rates
    par.on_fault_vars[0,iz,ix,20] = par.on_fault_vars[0,iz,ix,11]/par.creep_slip_rate # initial state var.
    par.on_fault_vars[0,iz,ix,7] = par.init_norm # initial normal stress.
    par.on_fault_vars[0,iz,ix,8] = shear_steady_state(par.on_fault_vars[0,iz,ix,9], 
                                                par.on_fault_vars[0,iz,ix,10],
                                                par.on_fault_vars[0,iz,ix,12],
                                                par.on_fault_vars[0,iz,ix,13],
                                                par.creep_slip_rate,
                                                par.on_fault_vars[0,iz,ix,7],
                                                par.on_fault_vars[0,iz,ix,46])
###############################################
##### Domain boundaries for transferring ######
###############################################
#par.xmin_trans, par.xmax_trans = -25e3, 25e3
#par.zmin_trans = -25e3
#par.ymin_trans, par.ymax_trans = -5e3, 5e3
#par.dx_trans = 50 

####################################
##### HPC resource allocation ######
####################################
par.casename       = "liu2020.fdc.rough.250.quasi"
par.HPC_nnode      = 1 # Number of computing nodes. On LS6, one node has 128 CPUs.
par.HPC_ncpu       = 20 # Number of CPUs requested.
par.HPC_queue      = "normal" # q status. Depending on systems, job WALLTIME and Node requested.
par.HPC_time       = "48:00:00" # WALLTIME, in hh:mm:ss format.
par.HPC_account    = "EAR22012" # Project account to be charged SUs against.
par.HPC_email      = "dliu@ig.utexas.edu" # Email to receive job status.

##############################################
##### Single station time series output ######
##############################################

# (x,z) coordinate pairs for on-fault stations (in km).
par.st_coor_on_fault = [[0.0, -3.0], [0.0,-7.5], [0.0, -12.0], [9.0,-7.5], \
   [12.0, -3.0], [12.0,-12.0], [15.0, -7.5], [18.0,-7.5], [-9.0,-7.5], \
   [-12.0,-3.0], [-12.0,-12.0], [-15.0, -7.5], [-18.0, -7.5]]
   
# (x,y,z) coordinates for off-fault stations (in km).
par.st_coor_off_fault = [[0,9,0], [0,-9,0], [12,6,0], [12,-6,0], [-12,6,0], \
   [-12,-6,0]]

par.n_on_fault        = len(par.st_coor_on_fault)
par.n_off_fault       = len(par.st_coor_off_fault)

# Additional par.solver options for AZTEC
par.az_op = 2 # AZTEC options
par.az_maxiter = 2000 # maximum iteration for AZTEC
par.az_tol = 1.0e-7 # tolerance for solution in AZTEC.


