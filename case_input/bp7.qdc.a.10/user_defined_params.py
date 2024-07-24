#! /usr/bin/env python3

from defaultParameters import parameters
import numpy as np
from math import *

par = parameters()
# cylce id. Simulate quasi-dynamic earthquake cycles from istart to iend.

par.istart = 1
par.iend = 3
# mode of the code - quasi-dynamic (1) or fully-dynamic (2). 
par.mode = 1

# model_domain (in meters)
par.xmin, par.xmax = -500, 500
par.ymin, par.ymax = -500, 500
par.zmin, par.zmax = -500, 500

# creeping zone bounaries.
# creeping zones are assinged on the lateral sides and bottom of 
# the RSF controlled region and will slide at fixed loading slip rate.
par.xminc, par.xmaxc, par.zminc = -400, 400, -400

par.dx = 10.0e0 # cell size, spatial resolution
par.nuni_y_plus, par.nuni_y_minus = 5, 5 # along the fault-normal dimension, the number of cells share the dx cell size.
par.enlarging_ratio = 1.0e0 # along the fault-normal dimension (y), cell size will be enlarged at this ratio compoundly.

# Isotropic material propterty.
# Vp, Vs, Rou
par.vp, par.vs, par.rou = 6.0e3, 3.464e3, 2.67e3
par.init_norm = -25.0e6 # initial normal stress in Pa. Negative compressive.

# Controlling switches for EQquasi system
par.rough_fault = 0 # include rough fault yes(1) or not(0).
par.rheology    = 1 # elastic(1). 
par.friclaw     = 3 # rsf_aging(3), rsf_slip(4).
par.ntotft      = 1 # number of total faults.
par.solver      = 1 # solver option. MUMPS(1, recommended). AZTEC(2).
par.nstep       = 10000 # total num of time steps for exiting, if not exit via sliprate threshold
par.nt_out      = 100 # Every nt_out time steps, disp of the whole model and on-fault variables will be written out in netCDF format.
par.bp          = 7 
# currently supported cases
# 5 (SCEC-BP5)
# 1001 (GM-cycle)

# xi, minimum Dc
par.xi = 0.05 # xi used to limit variable time step size. See Lapusta et al. (2009).
par.minDc = 0.5e-3 # meters

# loading 
par.far_vel_load = 5e-10 # far field loading velocity on xz planes. A minus value is applied on the other side.
par.creep_slip_rate = 1.0e-9 # creeping slip rate outside of RSF controlled region.
par.exit_slip_rate = 1.0e-3 # exiting slip rate for EQquasi [m/s].

#################################
##### Frictional variables ######
#################################
# friclaw == 1, slip weakening
par.fric_sw_fs = 0
par.fric_sw_fd = 0
par.fric_sw_D0 = 0
# friclaw == 3, rate- and state- friction with aging law.
par.fric_rsf_a, par.fric_rsf_b, par.fric_rsf_Dc = 0.004, 0.01, 0.5e-3
par.fric_rsf_deltaa = 0.012
par.fric_rsf_r0 = 0.6
par.fric_rsf_v0 = 1e-6
# Creating the fault interface
par.nfx = round((par.xmax - par.xmin)/par.dx + 1)
par.nfz = round((par.zmax - par.zmin)/par.dx + 1)
par.fx = np.linspace(par.xmin,par.xmax,par.nfx) # coordinates of fault grids along strike.
par.fz = np.linspace(par.zmin,par.zmax,par.nfz) # coordinates of fault grids along dip.
# Create on_fault_vars array for on_fault varialbes.
par.on_fault_vars = np.zeros((par.nfz,par.nfx,100))
def shear_steady_state(a,b,v0,r0,load_rate,norm,slip_rate, rou, vs):
  # calculate shear stress at steady state
  res = -norm*a*asinh(slip_rate/2.0/v0*exp((r0+b*log(v0/load_rate))/a)) + rou*vs/2.0*slip_rate
  return res
  
rad = 200 # radius of the velocity weakening zone
for ix, xcoor in enumerate(par.fx):
  for iz, zcoor in enumerate(par.fz):
  # assign a in RSF. a is a 2D distribution.
    radii = sqrt(xcoor*xcoor + zcoor*zcoor)
    if radii <= rad:
      par.on_fault_vars[iz,ix,9] = par.fric_rsf_a
    else: 
      par.on_fault_vars[iz,ix,9] = par.fric_rsf_a + par.fric_rsf_deltaa
      
    par.on_fault_vars[iz,ix,10] = par.fric_rsf_b # assign b in RSF 
    par.on_fault_vars[iz,ix,11] = par.fric_rsf_Dc # assign Dc in RSF.
    par.on_fault_vars[iz,ix,12] = par.fric_rsf_v0 # initial reference slip rate.
    par.on_fault_vars[iz,ix,13] = par.fric_rsf_r0 # initial reference friction.
    
    par.on_fault_vars[iz,ix,46] = par.creep_slip_rate # initial slip rates
    par.on_fault_vars[iz,ix,20] = par.on_fault_vars[iz,ix,11]/par.creep_slip_rate # initial state var.
    par.on_fault_vars[iz,ix,7] = par.init_norm # initial normal stress.
    par.on_fault_vars[iz,ix,8] = shear_steady_state(par.on_fault_vars[iz,ix,9], 
                                                par.on_fault_vars[iz,ix,10],
                                                par.on_fault_vars[iz,ix,12],
                                                par.on_fault_vars[iz,ix,13],
                                                par.creep_slip_rate,
                                                par.on_fault_vars[iz,ix,7],
                                                par.on_fault_vars[iz,ix,46],
                                                par.rou,
                                                par.vs)

####################################
##### HPC resource allocation ######
####################################
par.casename  = "bp7.qdc.a.10"
par.HPC_nnode = 1 # Number of computing nodes. On LS6, one node has 128 CPUs.
par.HPC_ncpu  = 20 # Number of CPUs requested.
par.HPC_queue = "normal" # q status. Depending on systems, job WALLTIME and Node requested.
par.HPC_time  = "35:00:00" # WALLTIME, in hh:mm:ss format.
par.HPC_account = "EAR22012" # Project account to be charged SUs against.
par.HPC_email = "dliu@ig.utexas.edu" # Email to receive job status.

##############################################
##### Single station time series output ######
##############################################
# (x,z) coordinate pairs for on-fault stations (in km).
par.st_coor_on_fault  = [[0, 0], [-100,0], [0,100], [100,0], \
                    [0,-100], [-100,-100], [-100,100], [100,-100], [100,100], \
                    [-300,0], [0, 300], [300,0], [0,-300]]
par.st_coor_on_fault  = np.asarray(par.st_coor_on_fault)/1000.0
# (x,y,z) coordinates for off-fault stations (in km).
par.st_coor_off_fault = [[0,200,0], [0,400,0], [-300,400,0], [0,400,-300], [300,400,0], \
                     [0,400,300]]
par.st_coor_off_fault = np.asarray(par.st_coor_off_fault)/1000.0
par.n_on_fault = len(par.st_coor_on_fault)
par.n_off_fault = len(par.st_coor_off_fault)

# Additional solver options for AZTEC
par.az_op = 2 # AZTEC options
par.az_maxiter = 2000 # maximum iteration for AZTEC
par.az_tol = 1.0e-7 # tolerance for solution in AZTEC.


