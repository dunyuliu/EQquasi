#! /usr/bin/env python3

import numpy as np
from math import *

# cylce id. Simulate quasi-dynamic earthquake cycles from istart to iend.
istart = 1
iend = 1
# mode of the code - quasi-dynamic (1) or fully-dynamic (2). 
mode = 1

# model_domain (in meters)
xmin, xmax = -60.0e3, 60.0e3
ymin, ymax = -50.0e3, 50.0e3
zmin, zmax = -60.0e3, 0.0e3

# creeping zone bounaries.
# creeping zones are assinged on the lateral sides and bottom of 
# the RSF controlled region and will slide at fixed loading slip rate.
xminc, xmaxc, zminc = -50.0e3, 50.0e3, -40.0e3

dx = 2000.0e0 # cell size, spatial resolution
nuni_y_plus, nuni_y_minus = 5, 5 # along the fault-normal dimension, the number of cells share the dx cell size.
enlarging_ratio = 1.3e0 # along the fault-normal dimension (y), cell size will be enlarged at this ratio compoundly.

# Isotropic material propterty.
# Vp, Vs, Rou
vp, vs, rou = 6.0e3, 3.464e3, 2.67e3
init_norm = -25.0e6 # initial normal stress in Pa. Negative compressive.

# Controlling switches for EQquasi system
rough_fault = 0 # include rough fault yes(1) or not(0).
rheology    = 1 # elastic(1). 
friclaw     = 3 # rsf_aging(3), rsf_slip(4).
ntotft      = 1 # number of total faults.
solver      = 1 # solver option. MUMPS(1, recommended). AZTEC(2).
nt_out      = 2 # Every nt_out time steps, disp of the whole model and on-fault variables will be written out in netCDF format.
bp          = 5 
# currently supported cases
# 5 (SCEC-BP5)
# 1001 (GM-cycle)

# xi, minimum Dc
xi = 0.015 # xi used to limit variable time step size. See Lapusta et al. (2009).
minDc = 0.13 # meters

# loading 
far_vel_load = 4e-10 # far field loading velocity on xz planes. A minus value is applied on the other side.
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
nfx = int((xmax - xmin)/dx + 1)
nfz = int((zmax - zmin)/dx + 1)
fx = np.linspace(xmin,xmax,nfx) # coordinates of fault grids along strike.
fz = np.linspace(zmin,zmax,nfz) # coordinates of fault grids along dip.
# Create on_fault_vars array for on_fault varialbes.
on_fault_vars = np.zeros((nfz,nfx,100))
def shear_steady_state(a,b,v0,r0,load_rate,norm,slip_rate):
  # calculate shear stress at steady state
  res = -norm*a*asinh(slip_rate/2.0/v0*exp((r0+b*log(v0/load_rate))/a)) + rou*vs/2.0*slip_rate
  return res
  
for ix, xcoor in enumerate(fx):
  for iz, zcoor in enumerate(fz):
  # assign a in RSF. a is a 2D distribution.
    if abs(zcoor)>=18e3 or abs(xcoor)>=32e3 or abs(zcoor)<=2e3: 
      on_fault_vars[iz,ix,9] = fric_rsf_a + fric_rsf_deltaa
    elif abs(zcoor)<=16e3 and abs(zcoor)>=4e3 and abs(xcoor)<=30e3:
      on_fault_vars[iz,ix,9] = fric_rsf_a
    else:
      tmp1 = (abs(abs(zcoor)-10e3) - 6e3)/2e3
      tmp2 = (abs(xcoor)-30e3)/2e3
      on_fault_vars[iz,ix,9] = fric_rsf_a + max(tmp1,tmp2)*fric_rsf_deltaa
    on_fault_vars[iz,ix,10] = fric_rsf_b # assign b in RSF 
    on_fault_vars[iz,ix,11] = fric_rsf_Dc # assign Dc in RSF.
    if (xcoor<=-18e3 and xcoor>=-30e3 and zcoor<=-4e3 and zcoor>=-16e3):
      on_fault_vars[iz,ix,11] = minDc # a special Dc zone.
    on_fault_vars[iz,ix,12] = fric_rsf_v0 # initial reference slip rate.
    on_fault_vars[iz,ix,13] = fric_rsf_r0 # initial reference friction.
    
    on_fault_vars[iz,ix,46] = creep_slip_rate # initial slip rates
    if (xcoor<=-18e3 and xcoor>=-30e3 and zcoor<=-4e3 and zcoor>=-16e3):
      on_fault_vars[iz,ix,46] = 0.03 # initial high slip rate patch.
    on_fault_vars[iz,ix,20] = on_fault_vars[iz,ix,11]/creep_slip_rate # initial state var.
    on_fault_vars[iz,ix,7] = init_norm # initial normal stress.
    on_fault_vars[iz,ix,8] = shear_steady_state(on_fault_vars[iz,ix,9], 
                                                on_fault_vars[iz,ix,10],
                                                on_fault_vars[iz,ix,12],
                                                on_fault_vars[iz,ix,13],
                                                creep_slip_rate,
                                                on_fault_vars[iz,ix,7],
                                                on_fault_vars[iz,ix,46])
###############################################
##### Domain boundaries for transferring ######
###############################################
xmin_trans, xmax_trans = -25e3, 25e3
zmin_trans = -25e3
ymin_trans, ymax_trans = -5e3, 5e3
dx_trans = 50 

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


