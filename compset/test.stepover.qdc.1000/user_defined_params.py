#! /usr/bin/env python3
# compset: test.stepover.qdc.1000
# gate   : e2e fast
# changed: 2026-08-12   (content; see git log for the full history)
# note   : Fast-tier row; the only ntotft > 1 case in every-push CI.
# Status is mirrored in compset/README.md, which a contract test checks.
# ----------------------------------------------------------------------
"""Two-fault step-over, BP5 parameters, geometry changed only.

Segment A: x in [-42.5, -2.5] km, y = 0 km      (the "cuts through the model" side is x = -42.5 km)
Segment B: x in [-7.5, 32.5] km,  y = -5.0 km    (the "cuts through the model" side is x = +32.5 km)
Both z in [-20, 0] km, dx = 1000 m, right-lateral, releasing step-over.

Creeping (velocity-strengthening, a-b > 0) is assigned on:
  - the side of each segment that reaches the model's own x boundary (the
    fault "cuts through the model" there),
  - the bottom (z near -20 km) and a band at the top (z near 0) of both
    segments.
The interior tips that face each other across the step-over (segment A's
+x end at -2.5 km, segment B's -x end at -7.5 km) are left
velocity-weakening -- an idealization the project owner ruled acceptable
for this case; a real step-over would taper these too.

This is the deliverable for on-fault-input multi-fault support: proof that
netcdf_read_on_fault (src/netcdf_io.f90) and case.setup's
netcdf_write_on_fault_vars correctly route each fault's own (a, b, Dc, v0,
r0, init_slip_rate, init_shear_stress, init_normal_stress, init_state) onto
that fault's own nodes, not onto whichever fault happens to be first.
"""

from defaultParameters import parameters
import numpy as np
from math import *

par = parameters()

par.istart = 1
par.iend = 1
par.mode = 1  # quasi-dynamic

# Two-fault geometry: (xlo, xhi, ycoor, zlo, zhi) per fault, meters.
par.ntotft = 2
FAULT_A = (-42.5e3, -2.5e3, 0.0, -20.0e3, 0.0)
FAULT_B = (-7.5e3, 32.5e3, -5.0e3, -20.0e3, 0.0)
par.faultgeom = [FAULT_A, FAULT_B]

# Model domain: x, z bound the union of both faults; y spans both fault
# planes (0 and -5 km) with comparable far-field margin on each side.
par.fxmin, par.fxmax = -42.5e3, 32.5e3
par.fymin, par.fymax = -35.0e3, 30.0e3
par.fzmin, par.fzmax = -20.0e3, 0.0e3

# The hard RSF/non-RSF cutoff (faulting.f90) is a single global box, not
# per-fault, so it is set to the union of both faults' extents: every fault
# node is inside it, and the a-b sign (assigned below, per fault) alone
# controls which nodes creep.
par.xminc, par.xmaxc, par.zminc = par.fxmin, par.fxmax, par.fzmin

par.dx = 1000.0e0
par.dy = par.dx
par.dz = par.dx
par.nuni_y_plus, par.nuni_y_minus = 5, 5
par.enlarging_ratio = 1.3e0

par.vp, par.vs, par.rou = 6.0e3, 3.464e3, 2.67e3
par.init_norm = -25.0e6

par.insertFaultType = 0
par.rough_fault = par.insertFaultType
par.rheology    = 1
par.friclaw     = 3  # rsf_aging
par.solver      = 1
par.nstep       = 101
par.nt_out      = 100
par.bp          = 5  # BP5 parameters, step-over geometry

par.xi = 0.015
par.minDc = 0.13

par.far_vel_load = 4e-10  # right-lateral
par.creep_slip_rate = 1.0e-9
par.exit_slip_rate = 1.0e-3

par.fric_sw_fs = 0
par.fric_sw_fd = 0
par.fric_sw_D0 = 0
par.fric_rsf_a, par.fric_rsf_b, par.fric_rsf_Dc = 0.004, 0.03, 0.14
par.fric_rsf_deltaa = 0.036
par.fric_rsf_r0 = 0.6
par.fric_rsf_v0 = 1e-6


def shear_steady_state(a, b, v0, r0, load_rate, norm, slip_rate, rou, vs):
    res = -norm * a * asinh(slip_rate / 2.0 / v0 * exp((r0 + b * log(v0 / load_rate)) / a)) \
        + rou * vs / 2.0 * slip_rate
    return res


# Per-fault local (strike, dip) grids, padded to the largest fault's
# extent -- here both segments are 40 km x 20 km so nfxMax == nfx(ift) and
# nfzMax == nfz(ift) for both, but the code does not assume that.
nfx = [round((xhi - xlo) / par.dx + 1) for xlo, xhi, ycoor, zlo, zhi in par.faultgeom]
nfz = [round((zhi - zlo) / par.dz + 1) for xlo, xhi, ycoor, zlo, zhi in par.faultgeom]
nfxMax, nfzMax = max(nfx), max(nfz)

# (ntotft, nfzMax, nfxMax, 100): the multi-fault on_fault_vars layout that
# case.setup's netcdf_write_on_fault_vars detects from ndim == 4 and writes
# with an explicit leading fault dimension, matching what
# netcdf_read_on_fault (src/netcdf_io.f90) expects.
par.on_fault_vars = np.zeros((par.ntotft, nfzMax, nfxMax, 100))

# Column 99 (unused by the Fortran reader) carries fault*100 + local strike
# index, so on_fault_vars_input.nc can be inspected directly to confirm
# case.setup routed each fault's own slab to the right place.
#
# init_norm also differs by fault (a plausible per-segment choice, not a
# BP5 value) so the *Fortran-side* read can be checked independently: it
# lands in fric(7, n, ift) and is written back out as "effective_normal" in
# fault.00101.nc, so a swapped-fault or wrong-local-index read shows up as
# the wrong constant on the wrong fault's nodes rather than a value that
# happens to look plausible either way.
init_norm_by_fault = [par.init_norm, par.init_norm - 2.0e6]  # -25, -27 MPa... Pa

for ift, (xlo, xhi, ycoor, zlo, zhi) in enumerate(par.faultgeom):
    init_norm = init_norm_by_fault[ift]
    fx = np.linspace(xlo, xhi, nfx[ift])
    fz = np.linspace(zlo, zhi, nfz[ift])
    for ix, xcoor in enumerate(fx):
        # Distance from this segment's "cuts through the model" (creeping)
        # edge: fault A's is xlo, fault B's is xhi.
        d_creep_edge = (xcoor - xlo) if ift == 0 else (xhi - xcoor)
        for iz, zcoor in enumerate(fz):
            if abs(zcoor) >= 18e3 or abs(zcoor) <= 2e3 or d_creep_edge <= 2e3:
                a = par.fric_rsf_a + par.fric_rsf_deltaa  # creeping (VS)
            elif 4e3 <= abs(zcoor) <= 16e3 and d_creep_edge >= 4e3:
                a = par.fric_rsf_a  # locked (VW) -- includes the interior tip
            else:
                tmp1 = (abs(abs(zcoor) - 10e3) - 6e3) / 2e3
                tmp2 = (4e3 - d_creep_edge) / 2e3
                a = par.fric_rsf_a + max(tmp1, tmp2) * par.fric_rsf_deltaa

            par.on_fault_vars[ift, iz, ix, 9]  = a
            par.on_fault_vars[ift, iz, ix, 10] = par.fric_rsf_b
            par.on_fault_vars[ift, iz, ix, 11] = par.fric_rsf_Dc
            par.on_fault_vars[ift, iz, ix, 12] = par.fric_rsf_v0
            par.on_fault_vars[ift, iz, ix, 13] = par.fric_rsf_r0
            par.on_fault_vars[ift, iz, ix, 46] = par.creep_slip_rate
            par.on_fault_vars[ift, iz, ix, 20] = par.fric_rsf_Dc / par.creep_slip_rate
            par.on_fault_vars[ift, iz, ix, 7]  = init_norm
            # Steady state for *this node's own* initial slip rate, [..., 46],
            # not for creep_slip_rate. This case seeds no nucleation patch, so
            # the two are equal here and the value is unchanged; it is written
            # this way so that adding a patch later cannot silently pair a
            # 0.03 m/s node with creep-consistent shear. That inconsistency
            # collapses V to the creep rate on step 1, which sends the adaptive
            # step dtev = xi*Dc/V to ~2e6 s and flattens the model.
            par.on_fault_vars[ift, iz, ix, 8]  = shear_steady_state(
                a, par.fric_rsf_b, par.fric_rsf_v0, par.fric_rsf_r0,
                par.creep_slip_rate, init_norm,
                par.on_fault_vars[ift, iz, ix, 46],
                par.rou, par.vs)
            par.on_fault_vars[ift, iz, ix, 99] = ift * 100 + ix  # mapping marker, not physical

###############################################
##### Domain boundaries for transferring #######
###############################################
par.xmin_trans, par.xmax_trans = -25e3, 25e3
par.zmin_trans = -20e3
par.ymin_trans, par.ymax_trans = -10e3, 5e3
par.dx_trans = 50

####################################
##### HPC resource allocation ######
####################################
par.casename = "stepover-qd-1000"
par.HPC_nnode = 1
par.HPC_ncpu = 2
par.HPC_queue = "normal"
par.HPC_time = "00:30:00"
par.HPC_account = "EAR22013"
par.HPC_email = "dliu@ig.utexas.edu"

##############################################
##### Single station time series output ######
##############################################
# Same (x, z) list is written for every fault (case.setup does not yet
# support per-fault station lists), so pick points that fall inside BOTH
# segments' along-strike range: the intersection is x in [-7.5, -2.5] km.
# A station must land on a node measured from the fault CORNER: both corners
# sit at x = 500 mod 1000 m (-42.5 km, -7.5 km), so with dx = 1000 the x
# coordinates must be half-km values; integer-km stations match nothing
# (the solver now warns, eqquasi.f90).
# And a station must carry SIGNAL to be comparable across machines: a node
# in the locked VW core produces only MPI-reduction noise in this 101-step
# window, and its integrated slip columns accumulate that noise to ~1e-16,
# machine-dependent -- the first CI run of this row failed on exactly that.
# All three stations therefore sit in creeping/taper zones (|z| <= 3 km or
# |z| >= 18 km).
par.st_coor_on_fault = [[-6.5, -2.0], [-5.5, -3.0], [-3.5, -18.0]]
par.st_coor_off_fault = [[0, 5, 0], [0, 5, -10], [-10, 10, 0]]
par.n_on_fault = len(par.st_coor_on_fault)
par.n_off_fault = len(par.st_coor_off_fault)

par.az_op = 2
par.az_maxiter = 2000
par.az_tol = 1.0e-7
