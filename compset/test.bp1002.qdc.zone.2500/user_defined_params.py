#! /usr/bin/env python3
# compset: test.bp1002.qdc.zone.2500
# gate   : e2e full row registered, gate: none -- UNVERIFIED (skips until
#          reference/test.bp1002.qdc.zone.2500 exists; see compset/README.md)
# changed: 2026-09-25
# note   : 101-step twin of bp1002.qdc.zone.2500 (unequal segment lengths,
#          identical VW patches). Same mesh; nstep cut only.
# Status is mirrored in compset/README.md, which a contract test checks.
# ----------------------------------------------------------------------
"""BP1002 zoning experiment: unequal segment lengths, identical VW patches.

Segment A: x in [-55.0,   2.5] km, y =  0.0 km   (cuts the model at x = -55 km)
Segment B: x in [ -2.5, 100.0] km, y = -5.0 km   (cuts the model at x = +100 km)
Both z in [-20, 0] km, dx = 2500 m, right-lateral, releasing step-over.
A is 57.5 km long, B is 102.5 km; they overlap by 5 km along strike and are
offset by 5 km across strike, exactly as in bp1002.qdc.2500. The model spans
x in [-55, 100] km (155 km); each segment still runs from its own outer
boundary to 2.5 km past the step-over, so both keep cutting the mesh edge.

What is varied, and what is deliberately NOT
--------------------------------------------
Relative to bp1002.qdc.2500 the ONLY changes are FAULT_A's xlo (-60 -> -55 km),
FAULT_B's xhi (+60 -> +100 km), and fxmin/fxmax to match. The frictional
zoning below is written in the model's own |x| (38 / 40 / 50 km from the
step-over at x = 0) and par.xminc/xmaxc stay at -50 / +50 km, so every
segment carries the SAME absolute velocity-weakening patch (|x| <= 38 km,
4 <= |z| <= 16 km), the same 2 km transition, and the same 10 km
velocity-strengthening band (40 <= |x| <= 50 km), at the same distance from
the step-over. Only the outer non-rate-and-state fixed-creep margin
(|x| > 50 km, faulting.f90:150, driven at the load rate) changes: 5 km on A
(two nodes, -55 and -52.5 km) against 50 km on B, a factor of ten. This is
the owner's "hold VW AREA fixed" experiment (PATHWAY_FORWARD row 7,
2026-09-24), as opposed to holding the VW *fraction* of each segment fixed,
which would be a different experiment.

Why A is the short one: A carries the cycle-0 nucleation patch, so cycle 0
ruptures A by construction whatever its length. In bp1002.qdc.2500 the
rupturing segment then alternates in later cycles. If A keeps rupturing here
even though it is the shorter, less-loaded segment, the "same segment
regardless of length" reading survives; if B takes over earlier or A never
ruptures again, segment length (the creeping-margin loading) matters.

Why not shorter than 55 km: read_input.f90 refuses xminc < fxmin, and a
segment whose outer end is inside |x| <= 50 km would lose its fixed-creep
margin entirely (at 40 km it would lose the VS band too and the VW patch
would sit 2 km from the mesh edge) -- a different boundary condition, not a
length change. 55 km keeps two creep nodes on A, the minimum that still
represents "the fault continues beyond the mesh".

The experiment is run with par.iend = 3 overridden in the case (three
cycles); the compset itself keeps bp1002's iend = 1 so the e2e cycle0
convention holds.

dx = 2500 m, not BP5's production 2000 m: the fault-normal offset between the
two segments is 5000 m, and EQquasi's mesh (src/func_lib.f90,
build_yline_belt) requires that offset to be an integer multiple of dy = dx --
it refuses to mesh, rather than silently produce a fault with zero nodes,
otherwise. 5000 / 2000 = 2.5 (refused); 5000 / 2500 = 2 (meshes). The segment ends are
commensurate at 2500 m too: -60000, -2500, 2500 and 60000 are all exact
multiples of it.

Frictional zoning follows BP5's, in the model's own x:
  - |x| > 50 km   non-rate-and-state, fixed creep (par.xminc/xmaxc),
  - |x| >= 40 km  velocity-strengthening,
  - |x| <= 38 km  velocity-weakening where 4 <= |z| <= 16 km,
  - 2 km transitions, as BP5 uses.
The interior tips that face each other across the step-over (segment A's
+x end at +2.5 km, segment B's -x end at -2.5 km) are left
velocity-weakening -- an idealization the project owner ruled acceptable
for this case; a real step-over would taper these too. It is not free: in a
10-cycle run of the earlier geometry the effective normal stress reached zero
at one of these tips in cycle 3 and the run stopped there (STOP 508).

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
par.iend = 1   # one cycle, beginning to end
par.mode = 1  # quasi-dynamic

# Two-fault geometry: (xlo, xhi, ycoor, zlo, zhi) per fault, meters.
par.ntotft = 2
FAULT_A = (-55.0e3,   2.5e3,  0.0, -20.0e3, 0.0)   # 57.5 km, short
FAULT_B = ( -2.5e3, 100.0e3, -5.0e3, -20.0e3, 0.0)   # 102.5 km, long
par.faultgeom = [FAULT_A, FAULT_B]

# Model domain.
#
# x is not free: each segment cuts through the model at its outer end, by
# design, so the lateral bounds ARE the union of the two segments: here
# [-55, 100] km, 155 km along strike against bp1002's 120 km. The creeping
# band at each outer end is what represents the fault continuing beyond the
# mesh; its length is the one thing this compset varies between segments.
#
# y and z do follow BP5, and this is the part that was wrong before. The
# earlier version put fzmin at -20 km, the faults' own base, so each fault
# met the mesh boundary in two directions at once and its corner had no
# elastic volume to relax into. BP5 keeps 40 km of material below a fault
# that bottoms at 20 km; that margin is reproduced here.
par.fxmin, par.fxmax = -55.0e3, 100.0e3
par.fymin, par.fymax = -50.0e3, 50.0e3
par.fzmin, par.fzmax = -60.0e3, 0.0e3

# The hard RSF/non-RSF cutoff (faulting.f90) is a single global box, not
# per-fault, so it is set to the union of both faults' extents: every fault
# node is inside it, and the a-b sign (assigned below, per fault) alone
# controls which nodes creep.
# Tracks the faults, never the domain: tying it to fxmin/fzmin is what put
# the rate-and-state box on the mesh boundary before.
# BP5's non-rate-and-state margin, kept at BP5's |x| > 50 km and NOT moved
# with the segment ends: outside this box faulting.f90 drives the fault at a
# fixed creep rate rather than by rate and state, so the outer 5 km of A and
# the outer 50 km of B are kinematic boundaries. Moving xminc/xmaxc with the
# ends would change the VS band each segment sees, which is exactly what this
# compset holds fixed. zminc is the fault base here because these segments only
# reach -20 km, so there is nothing below them to exclude.
par.xminc, par.xmaxc, par.zminc = -50.0e3, 50.0e3, -20.0e3

par.dx = 2500.0e0
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
# extent -- here A has nfx = 24 and B nfx = 42, so A's slab is zero-padded
# beyond its own 24 strike nodes; case.setup writes nfx per fault.
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
init_norm_by_fault = [par.init_norm, par.init_norm]  # BP5: -25 MPa on both

for ift, (xlo, xhi, ycoor, zlo, zhi) in enumerate(par.faultgeom):
    init_norm = init_norm_by_fault[ift]
    fx = np.linspace(xlo, xhi, nfx[ift])
    fz = np.linspace(zlo, zhi, nfz[ift])
    for ix, xcoor in enumerate(fx):
        for iz, zcoor in enumerate(fz):
            # BP5's zoning, in the model's own x rather than per-segment.
            #
            # BP5 wraps an 868 km^2 velocity-weakening patch in 3416 km^2 of
            # velocity-strengthening and 3280 km^2 of non-rate-and-state
            # creep -- nearly 4:1 -- and that margin is what arrests its
            # rupture. The earlier bp1002 zoning measured everything from
            # each segment's own cutting edge and left ~51 % of every segment
            # velocity-weakening, so rupture ran end to end and slipped 9.5 m.
            # The same three zones are used here, at BP5's widths:
            #
            #   |x| >  50 km   non-RSF, fixed creep (set by xminc/xmaxc)
            #   |x| >= 40 km   velocity-strengthening
            #   38 < |x| < 40  transition, 2 km as in BP5
            #   |x| <= 38 km   velocity-weakening, if z also allows
            #
            # z is BP5's unchanged: VW for 4 <= |z| <= 16 km, VS beyond 18 km
            # and above 2 km, 2 km transitions between.
            #
            # The step-over sits at x = 0, so the interior tips facing each
            # other stay velocity-weakening -- that is the question the case
            # asks and it must not be tapered away.
            if abs(zcoor) >= 18e3 or abs(zcoor) <= 2e3 or abs(xcoor) >= 40e3:
                a = par.fric_rsf_a + par.fric_rsf_deltaa  # creeping (VS)
            elif 4e3 <= abs(zcoor) <= 16e3 and abs(xcoor) <= 38e3:
                a = par.fric_rsf_a  # locked (VW) -- includes the interior tips
            else:
                tmp1 = (abs(abs(zcoor) - 10e3) - 6e3) / 2e3
                tmp2 = (abs(xcoor) - 38e3) / 2e3
                a = par.fric_rsf_a + max(tmp1, tmp2) * par.fric_rsf_deltaa

            # Nucleation patch, BP5's 12 km, on fault A only and inside its
            # velocity-weakening region: low Dc plus a high initial slip rate.
            # Fault B has no patch -- whether it ruptures is the question this
            # case asks, so it must not be seeded.
            nucleate = (ift == 0 and -36.0e3 <= xcoor <= -24.0e3
                        and -16e3 <= zcoor <= -4e3)

            par.on_fault_vars[ift, iz, ix, 9]  = a
            par.on_fault_vars[ift, iz, ix, 10] = par.fric_rsf_b
            par.on_fault_vars[ift, iz, ix, 11] = par.minDc if nucleate else par.fric_rsf_Dc
            par.on_fault_vars[ift, iz, ix, 12] = par.fric_rsf_v0
            par.on_fault_vars[ift, iz, ix, 13] = par.fric_rsf_r0
            par.on_fault_vars[ift, iz, ix, 46] = 0.03 if nucleate else par.creep_slip_rate
            par.on_fault_vars[ift, iz, ix, 20] = (
                par.on_fault_vars[ift, iz, ix, 11] / par.creep_slip_rate)
            par.on_fault_vars[ift, iz, ix, 7]  = init_norm
            # The shear must be the steady state for *this node's own* initial
            # slip rate, not for the creep rate: inside the nucleation patch
            # that rate is 0.03 m/s, and pairing it with the creep-consistent
            # shear makes the initial condition self-inconsistent. The Newton
            # solve then resolves it in favour of the stress, dropping V to the
            # creep rate at step 1 -- which sends dtev = xi*Dc/V to ~2e6 s and
            # flattens the whole model. BP5 passes its own [..., 46] here.
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
par.casename = "bp1002-stepover-zone-qd-2500"
par.HPC_nnode = 1
par.HPC_ncpu = 3
par.HPC_queue = "normal"
par.HPC_time = "04:00:00"
par.HPC_account = "EAR22013"
par.HPC_email = "dliu@ig.utexas.edu"

##############################################
##### Single station time series output ######
##############################################
# Same (x, z) list is written for every fault (case.setup does not yet
# support per-fault station lists), so pick points that fall inside BOTH
# segments' along-strike range: the intersection is x in [-7.5, -2.5] km.
# The same (x, z) list is written for every fault, so the points must lie
# inside BOTH segments' along-strike range: the intersection is the 5 km
# overlap, x in [-2.5, 2.5] km.
#
# They must also land ON mesh nodes. dx = 2500 m and the segment ends are
# multiples of it, so the admissible x are -2.5, 0 and 2.5 km and the
# admissible z are multiples of 2.5 km. An earlier version asked for
# x = +-2.0 km and z = -2.0/-18.0 km, none of which is a node: those two
# stations were dropped with no message at all, leaving one station where
# three were requested.
#
# Note that only fault 1 gets station output -- src/library_output.f90
# writes them under `if (j==1)` and never reaches fault 2 -- so these are
# fault 1's stations regardless of the list being duplicated per fault.
par.st_coor_on_fault = [[-2.5, -2.5], [0.0, -10.0], [2.5, -17.5]]
par.st_coor_off_fault = [[0, 5, 0], [0, 5, -10], [-10, 10, 0]]
par.n_on_fault = len(par.st_coor_on_fault)
par.n_off_fault = len(par.st_coor_off_fault)

par.az_op = 2
par.az_maxiter = 2000
par.az_tol = 1.0e-7
