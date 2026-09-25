#! /usr/bin/env python3
# compset: enechelon.qdc.gap4.1000
# gate   : none -- UNVERIFIED (no reference; collaborator's multi-cycle production run)
# changed: 2026-09-25   (content; see git log for the full history)
# note   : 4 km step-over, VW core, dx=1000. Parameters match TriBIE var-tf-500-v2.dat.
# Status is mirrored in compset/README.md, which a contract test checks.
# ----------------------------------------------------------------------
"""Two en-echelon faults 4 km apart, with a velocity-weakening core.

Derived from enechelon.qdc.2000.creep by two changes:

  1. Fault A moves from y = -12 km to y = -4 km, so the step-over GAP is
     4 km instead of 12 km. (The gap is the fault-NORMAL offset, the third
     entry of each faultgeom tuple; z is the depth range, unchanged at
     [-20, 0] km on both segments.)
  2. Friction is no longer velocity-strengthening everywhere. There is a
     velocity-weakening core at a - b = -0.026, tapering linearly over a
     5 km transition to a - b = +0.010 at the patch edge.
  3. dx is 1000 m, not 2000, and Dc is 0.07 m, not 0.14.

PARAMETER PARITY WITH TriBIE, AND WHERE IT ENDS. Columns 1-5 of
TriBIE_coding/twofaults/var-tf-500-v2.dat are seff, Dc, a, b, v_ini
(prepare_input.py's write order). This compset still matches TriBIE on
Dc = 0.07, b = 0.03, v_ini = 1e-9, f0 = 0.6, v0 = 1e-6 and seff = 25 MPa, and
on the endpoints of a (0.004 core, 0.04 edge).

`a` ITSELF NO LONGER MATCHES. TriBIE's bp5_taper is a 5 km ramp anchored at
the patch edge on both axes; the zoning below is BP5's -- an explicit
weakening box with 2 km transitions. The two give different a everywhere
except deep inside the box and hard against the rim. Restoring parity means
reverting the zoning to the 5 km edge-anchored ramp, not adjusting a constant.
That 5 km value was verified against TriBIE rather than assumed: rebuilding a
from twofaults.gts with bp5_taper at w = 5000 m reproduces var-tf-500-v2.dat
to 4e-13, and inverting w = e/(1 - frac) gives exactly 5000.0 m across all
12948 tapered elements (w = 2000 m is off by 2.2e-2). var-tf-500.dat is the
same.

This zoning supersedes the VW_WIDTH_Z = 10 km down-dip narrowing that stood
between 2026-09-10 and 2026-09-23.

It was briefly set to 2 km on 2026-09-08 and reverted the same day. Beware the
quantity being compared when eyeballing a plot: the VELOCITY-STRENGTHENING
RIM is 1.39 km wide, not 5 km, because a - b crosses zero at
e = 0.2778 * w. The taper width and the rim width are different numbers.

KNOWN LIMITATION at this width. In cycle 0 the rupture arrested at x = +11 km,
five kilometres short of the interior tip at +16 km, so nothing reached the
step-over: this configuration cannot be asked whether rupture crosses it.
With the 2 km BP5 zoning below, the flat a - b = -0.026 core is 56 x 12 km
(|x_local| <= 28 km, 4 <= depth <= 16 km; 57 x 13 nodes at dx = 1 km), i.e.
22.1 x 4.7 h* at h* = 2.53 km, and a - b stays negative out to 58.9 x 14.9 km
(the zero crossing sits 0.026/0.036 * 2 km = 1.44 km into the taper). The
"50 x 10 km, 19.7 x 3.9 h*" figure that stood here (corrected 2026-09-25) was
the 5 km edge-anchored taper's core, not this file's. The taper rule is the same function:
TriBIE combines the along-strike and down-dip ramps with a Chebyshev max,
which is algebraically the min-distance-to-nearest-edge used below (checked
numerically on this grid, max difference 7e-18).

Dc is 0.07 rather than BP5's 0.14 for the reason TriBIE's prepare_input.py
gives: effective normal stress here is 25 MPa, half BP5's 50 MPa, and
h* ~ mu* b Dc / ((b-a)^2 seff) is what has to stay put. At Dc = 0.14 and
25 MPa, h* = 5.5 km against a 12 km-tall VW core -- 2.2 h*, marginal for
nucleation. At 0.07, h* = 2.5 km and the core is 4.7 h* tall.

That is also why dx is 1000 m. At 2000 m, Dc = 0.07 puts h*/dx at 1.27 and
Lambda_0/dx at 1.32: the nucleation zone would span barely one cell. Note
that case.setup would NOT warn there -- its threshold is 1.0 -- and its
message cites BP5 as 1.32, which is wrong (the same code computes 2.64 for
BP5's own parameters). At dx = 1000 both ratios land near BP5's real values.

  Fault A (west): x in [-44, 16] km, y = -4 km, z in [-20, 0] km
  Fault B (east): x in [  0, 60] km, y =  0 km, z in [-20, 0] km

Length 60 km, width 20 km, gap 4 km, along-strike overlap 16 km, both
vertical, both breaching the free surface, all four along-strike tips buried.
Geometry otherwise follows twofaults/generate_gmsh.py; see that directory and
compset/README.md for why the overlap is 16 km rather than 15.

The model keeps the 180 deg rotational symmetry of the .creep pair: about the
vertical axis at (x, y) = (8, -2) km, x -> 16 - x and y -> -4 - y maps fault A
onto fault B, and the domain (-84 <-> 100 in x, -54 <-> 50 in y) and the
antisymmetric far-field load are both invariant under it. The .creep run
reproduced that symmetry bit-exactly in every output field, which is the
cheapest available check that each fault's on-fault input is routed to its own
nodes. Keep the domain symmetric and the check keeps working.

TAPER AT THE FREE SURFACE. BP5's depth zoning strengthens the top 2 km, so
rupture does not reach the surface. That is BP5's own choice and is kept.
"""

from defaultParameters import (parameters, FR_RSF_A, FR_RSF_B, FR_RSF_DC,
    FR_RSF_V0, FR_RSF_F0, FR_VINIT, FR_STATE, FR_TNRM0, FR_TSTK0)
import numpy as np
from math import *

par = parameters()

# One earthquake per cycle: in quasi-dynamic mode the solver runs THROUGH the
# rupture and exits once slip rate has stayed below exit_slip_rate for 1e5 s
# (exitCriteria, src/solveTimeLoopMUMPS.f90). There is no "stop at t = T"
# switch for bp /= 8, and the clock restarts at zero each cycle, so a target
# duration is reached by accumulating cycles: 3500 yr is roughly 40 of them at
# a BP5-like recurrence. iend is set generously; watch the accumulated time and
# stop the job once it passes 3500 yr -- each cycle is packaged into its own
# result/cycleN as it finishes, so nothing already done is lost.
par.istart = 1
par.iend = 100
par.mode = 1  # quasi-dynamic
par.bp = 5    # BP5 parameters; bp only gates the bp7/bp8 special paths

# Two-fault geometry: (xlo, xhi, ycoor, zlo, zhi) per fault, meters.
par.ntotft = 2
FAULT_A = (-44.0e3, 16.0e3, -4.0e3, -20.0e3, 0.0)
FAULT_B = (  0.0e3, 60.0e3,  0.0e3, -20.0e3, 0.0)
par.faultgeom = [FAULT_A, FAULT_B]

# Model domain. y is symmetric about the two fault planes (50 km beyond each),
# which both loads the segments equally and preserves the rotational symmetry
# described above.
par.fxmin, par.fxmax = -84.0e3, 100.0e3
par.fymin, par.fymax = -54.0e3,  50.0e3
par.fzmin, par.fzmax = -50.0e3,   0.0e3

# The RSF/fixed-creep cutoff (faulting.f90) is a single global box, not
# per-fault: set to the model domain so every fault node is rate-and-state
# governed and the a-b sign alone controls behaviour.
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
par.solver      = 1  # MUMPS

# nstep is a CEILING per cycle, not a target: the physics exit above is what
# ends a cycle, so this is set high enough that no cycle is truncated mid-event.
#
# It is ALSO a checkpoint interval, which matters if a cycle outlasts the queue.
# solveTimeLoopMUMPS.f90 writes disp.r.nc and fault.r.nc on `it == nstep` as
# well as on the physics exit, so a cycle that runs out of steps still ends
# cleanly, is packaged into result/cycleN, and the next cycle restarts from it.
# If a cycle cannot fit in one job, set nstep ~ (usable walltime) / (seconds per
# step, from the run summary) and the cycles become checkpoints rather than
# earthquakes. Two consequences if you do: `time` restarts at zero each cycle
# (it always does), and the interseismic/coseismic status bookkeeping in
# exitCriteria resets, so an event split across a boundary may be counted twice
# in tdyna.txt. The physics is continuous either way.
par.nstep  = 200000
par.nt_out = 2000

par.xi = 0.015
# minDc drives the adaptive step, dtev = ksi*minDc/maxSlipRate
# (solveTimeLoopMUMPS.f90). It is the SMALLEST Dc on the fault, not a label:
# left at BP5's 0.13 while Dc is 0.07 every step would be ~1.9x too long.
par.minDc = 0.07

# +4e-10 is RIGHT-LATERAL: solveTimeLoopMUMPS.f90 gives the ymax boundary an
# x-velocity of +far_load_rate and ymin -far_load_rate, so the y > 0 block moves
# +x relative to y < 0. Going +x along fault A, fault B sits at +y, i.e. to the
# LEFT: a left step. A left step under right-lateral shear is RESTRAINING
# (compressional across the 4 km step), not releasing. (Corrected 2026-09-25;
# the earlier comment here said "releasing".)
par.far_vel_load = 4e-10
par.creep_slip_rate = 1.0e-9
par.exit_slip_rate = 1.0e-3       # ends each cycle after its earthquake

par.fric_sw_fs = 0
par.fric_sw_fd = 0
par.fric_sw_D0 = 0
par.fric_rsf_b, par.fric_rsf_Dc = 0.03, 0.07
par.fric_rsf_r0 = 0.6
par.fric_rsf_v0 = 1e-6

# Frictional zoning: copied from compset/bp5.qdc.2000 and applied to BOTH
# faults, in each fault's own local coordinates. An explicit velocity-weakening
# box with 2 km transitions, not a ramp anchored at the patch edge.
#
# Down-dip, BP5's zoning fits this 20 km fault exactly and is used verbatim:
# strengthening for depth <= 2 km and >= 18 km, weakening for 4 to 16 km, 2 km
# transitions between.
#
# Along strike, BP5 wants weakening for |x| <= 30 km with a taper out to 32 km
# -- 64 km in all -- and these segments are 60 km. The taper is therefore moved
# INSIDE the patch (VW_HALF_X = 28 km, taper 28 -> 30 km), which keeps BP5's
# 2 km transition width and leaves the tips velocity-strengthening. The
# alternative, applying |x| <= 30 verbatim, would put the taper outside the
# fault and leave a - b = -0.026 all the way to both tips.
#
# The two axes combine with a max -- the more strengthening wins -- which is
# BP5's own rule, unchanged.
VW_HALF_X = 28.0e3       # half-length of the weakening box, local x
TAPER_M   = 2.0e3        # transition width, BP5's value, both axes
par.fric_rsf_a      = 0.004   # BP5 core value      -> a - b = -0.026
par.fric_rsf_deltaa = 0.036   # BP5 strengthening increment -> a - b = +0.010

# Normal-stress caps ON. The .creep pair leaves them off so the unclamping is
# visible, which is right for a 101-step diagnostic. This case is different:
# it runs tens of cycles with a velocity-weakening core and a step-over only
# 4 km wide, and the facing interior tips unclamp cycle after cycle.
# bp1002.qdc.2500 -- same mechanism, a 5 km gap -- reached zero effective
# normal stress at one of those tips and died on STOP 508 in cycle 3. Losing a
# multi-week run that way is not a risk worth taking to keep a default.
# Liu, Duan & Luo (2020) section 3.5 is the reference for the regularisation.
par.C_normal_stress_caps = 1
par.min_norm = -10.0e6
par.max_norm = -40.0e6


def shear_steady_state(a, b, v0, r0, load_rate, norm, slip_rate, rou, vs):
    res = -norm * a * asinh(slip_rate / 2.0 / v0 * exp((r0 + b * log(v0 / load_rate)) / a)) \
        + rou * vs / 2.0 * slip_rate
    return res


nfx = [round((xhi - xlo) / par.dx + 1) for xlo, xhi, ycoor, zlo, zhi in par.faultgeom]
nfz = [round((zhi - zlo) / par.dz + 1) for xlo, xhi, ycoor, zlo, zhi in par.faultgeom]
nfxMax, nfzMax = max(nfx), max(nfz)

# (ntotft, nfzMax, nfxMax, 100): the multi-fault layout case.setup detects from
# ndim == 4 and netcdf_read_on_fault (src/netcdf_io.f90) expects.
par.on_fault_vars = np.zeros((par.ntotft, nfzMax, nfxMax, 100))

for ift, (xlo, xhi, ycoor, zlo, zhi) in enumerate(par.faultgeom):
    fx = np.linspace(xlo, xhi, nfx[ift])
    fz = np.linspace(zlo, zhi, nfz[ift])
    for ix, xcoor in enumerate(fx):
        for iz, zcoor in enumerate(fz):
            # Distance to the nearest edge of this fault's rectangle. Both
            # along-strike tips are buried on both segments, so both count.
            # Local coordinates: xl from this segment's own centre, zd = depth.
            xl = xcoor - 0.5 * (xlo + xhi)
            zd = abs(zcoor)
            if zd >= 18e3 or zd <= 2e3 or abs(xl) >= VW_HALF_X + TAPER_M:
                a = par.fric_rsf_a + par.fric_rsf_deltaa
            elif 4e3 <= zd <= 16e3 and abs(xl) <= VW_HALF_X:
                a = par.fric_rsf_a
            else:
                tmp1 = (abs(zd - 10e3) - 6e3) / TAPER_M
                tmp2 = (abs(xl) - VW_HALF_X) / TAPER_M
                a = par.fric_rsf_a + max(tmp1, tmp2) * par.fric_rsf_deltaa

            par.on_fault_vars[ift, iz, ix, FR_RSF_A]  = a
            par.on_fault_vars[ift, iz, ix, FR_RSF_B]  = par.fric_rsf_b
            par.on_fault_vars[ift, iz, ix, FR_RSF_DC] = par.fric_rsf_Dc
            par.on_fault_vars[ift, iz, ix, FR_RSF_V0] = par.fric_rsf_v0
            par.on_fault_vars[ift, iz, ix, FR_RSF_F0] = par.fric_rsf_r0
            par.on_fault_vars[ift, iz, ix, FR_VINIT]  = par.creep_slip_rate
            par.on_fault_vars[ift, iz, ix, FR_STATE]  = par.fric_rsf_Dc / par.creep_slip_rate
            par.on_fault_vars[ift, iz, ix, FR_TNRM0]  = par.init_norm
            # Steady state for this node's own a and its own initial slip rate.
            # No nucleation patch is seeded: the first event grows out of the
            # relaxation from V_init = 1e-9 toward the plate rate 8e-10 and the
            # stress concentration at the tips and across the step. Cycle 0 is
            # therefore spin-up, not a representative earthquake.
            par.on_fault_vars[ift, iz, ix, FR_TSTK0]  = shear_steady_state(
                a, par.fric_rsf_b, par.fric_rsf_v0, par.fric_rsf_r0,
                par.creep_slip_rate, par.init_norm,
                par.on_fault_vars[ift, iz, ix, FR_VINIT],
                par.rou, par.vs)
            # Unread by the Fortran; a marker so on_fault_vars_input.nc can be
            # inspected to confirm each fault's slab went to the right place.
            par.on_fault_vars[ift, iz, ix, 99] = ift * 100 + ix

###############################################
##### Domain boundaries for transferring #######
###############################################
par.xmin_trans, par.xmax_trans = -50e3, 66e3
par.zmin_trans = -20e3
par.ymin_trans, par.ymax_trans = -10e3, 6e3
par.dx_trans = 50

####################################
##### HPC resource allocation ######
####################################
par.casename = "enechelon-qd-gap4-1000"
par.HPC_nnode = 1
par.HPC_ncpu = 20   # case.setup estimates 18 for this mesh; measure scaling
par.HPC_queue = "normal"
par.HPC_time = "48:00:00"
par.HPC_account = "EAR22013"
par.HPC_email = "dliu@ig.utexas.edu"

##############################################
##### Single station time series output ######
##############################################
# case.setup writes the same list to every fault, so every station must fall
# inside BOTH segments' along-strike range -- the intersection is x in [0, 16]
# km -- and must land on a node: x = -44 km + n*1 km and z = 0 - n*1 km at
# dx = 1000, so any integer kilometre. The two velocity-weakening cores overlap only in
# x in [5, 11] km at 5 <= |z| <= 15 km, so the core stations sit there; the
# other two watch the free surface at the step and the interior tip.
par.st_coor_on_fault = [[6.0, -10.0], [8.0, -6.0], [8.0, -10.0], [10.0, -14.0],
                        [0.0, 0.0], [16.0, -2.0]]
# Mid-step is y = -2 km, a node line (the belt spans [-4, 0] at dy = 1 km).
par.st_coor_off_fault = [[8, -2, 0], [8, -2, -10], [8, 6, 0], [-20, -4, 0], [40, 0, 0]]
par.n_on_fault = len(par.st_coor_on_fault)
par.n_off_fault = len(par.st_coor_off_fault)

par.az_op = 2
par.az_maxiter = 2000
par.az_tol = 1.0e-7
