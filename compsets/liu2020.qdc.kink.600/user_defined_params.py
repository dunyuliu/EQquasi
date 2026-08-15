#! /usr/bin/env python3
# compset: liu2020.qdc.kink.600
# gate   : reference frozen (utilities-read); no e2e row -- a cycle is ~5 h
# changed: 2026-08-15
# note   : Paper params at dx=600 (Lambda0/dx=0.94, marginal). Coherent cycles; nstep=10000 truncates long ones.
# Status is mirrored in compsets/README.md, which a contract test checks.
# ----------------------------------------------------------------------
"""Liu, Duan & Luo (2020, GJI 220, 598-609): strike-slip fault with a kink.

Quasi-dynamic (EQquasi alone). The paper's model is fully dynamic, EQquasi
looping with EQdyna through EQsimu; this compset is the EQquasi-only half, so
it reproduces the geometry, the friction parameters and the interseismic and
nucleation phases, but not the coseismic dynamic rupture. Slip rate exits at
par.exit_slip_rate rather than being handed to EQdyna.

Parameters are the paper's table 1:

    fault length along strike   FS      60 km      -> x in [-30, 30] km
    fault depth                 FD      30 km      -> z in [-30, 0] km
    VW region length            Lseis   ~45 km
    VW region width             Wseis   11.4 km
    S-wave velocity             Vs      3464 m/s
    P-wave velocity             Vp      6000 m/s
    density                     rho     2670 kg/m3
    reference slip velocity     V0      1e-6 m/s
    reference friction          mu0     0.6
    RSF a (VW region)           a       0.007
    RSF b                       b       0.011
    initial normal stress       sigma_n -50 MPa
    initial shear stress        tau      30 MPa
    critical slip distance      L        11 mm
    element size                dx      300 m

The bend
--------
EQquasi cannot express a bent fault through par.faultgeom: that places each
fault on a constant-y plane, so a segment oblique to x has no representation.
insert_rough_fault (src/func_lib.f90) can -- it displaces each node's y by a
value read per (x, z) from bFault_Rough_Geometry.txt, which is how the rough
(fractal) faults work. A kink is just a different y(x): zero left of the bend,
(x - kinkX)*tan(kinkAngleDeg) right of it.

So the geometry comes from generateKinkGeometry.py in this directory, NOT
from scripts/generateFaultInterface, which only knows planar (1) and fractal
(2) surfaces. par.insertFaultType is therefore 0 -- case.setup must not call
generateFaultInterface and overwrite the kink -- while par.rough_fault is 1,
which is what the solver reads to decide whether to displace y at all.

    ./case.setup
    ./generateKinkGeometry.py .      # writes bFault_Rough_Geometry.txt

Set kinkAngleDeg to 0 to recover the paper's planar reference model without
changing anything else, which is the comparison figure 1(a) versus 1(b) makes.
"""

from defaultParameters import parameters
import numpy as np
from math import *

par = parameters()

par.dip = 90.          # vertical; the bend is in map view, not in dip

par.istart = 1
par.iend = 1
par.mode = 1           # quasi-dynamic. The paper's model is mode 2 via EQsimu.

# The kink, as module-level constants rather than attributes on par:
# par carries the solver's own parameter schema (scripts/defaultParameters.py)
# and nothing here is read by the solver. generateKinkGeometry.py imports this
# module and reads these directly, the way bp1002 exposes FAULT_A/FAULT_B.
KINK_ANGLE_DEG = 10.0   # paper's bend angle; 0 recovers the planar model
KINK_X         = 0.0    # bend at the fault centre, as in figure 1(b)

# Model domain. The fault is 60 km along strike and 30 km deep (table 1); the
# box extends well beyond it so the boundaries do not load the fault -- the
# paper notes EQquasi's model "is much smaller than that of EQdyna" but still
# applies tectonic loading on boundaries parallel to the fault.
par.fxmin, par.fxmax = -30.0e3, 30.0e3
par.fymin, par.fymax = -30.0e3, 30.0e3
par.fzmin, par.fzmax = -30.0e3, 0.0e3

# Creeping (non-rate-and-state) margin, inside the domain on both sides and
# below, as liu2020.fdc.planar had it.
par.xminc, par.xmaxc, par.zminc = -25.0e3, 25.0e3, -25.0e3

par.dx = 600.0e0       # 300 m (paper) needs 64-bit MUMPS; Lambda0/dx = 1.26 here vs the 2.3 the paper requires
par.dy = par.dx
par.dz = par.dx
par.nuni_y_plus, par.nuni_y_minus = 10, 10
par.enlarging_ratio = 1.3e0

par.vp, par.vs, par.rou = 6.0e3, 3.464e3, 2.67e3
par.init_norm = -50.0e6                       # table 1

# insertFaultType 0 keeps case.setup's hands off the geometry; rough_fault 1
# tells the solver to read bFault_Rough_Geometry.txt and displace y.
par.insertFaultType = 0
par.rough_fault = 1
par.rheology    = 1
par.friclaw     = 3    # rate and state, aging law -- the paper's eq. (7)
par.ntotft      = 1    # ONE fault: the kink is a shape, not a second fault
par.solver      = 1
par.nstep       = 10000
par.nt_out      = 1000
par.bp          = 1000  # liu2020, for description only

# The paper's epsilon, section 2.2 step (2): the time step is set from
# eps*L/max(sliprate), "we use eps = 0.2". The liu2020.fdc.* compsets in the
# repo carry 0.02, a factor of ten smaller; the paper's value is used here.
par.xi = 0.2
par.minDc = 11e-3

par.far_vel_load = 5e-10
par.creep_slip_rate = 1.0e-9
# 1.5e-2 m/s: in the coupled system this is where EQquasi hands over to
# EQdyna. Standing alone it is simply where the run stops, so a quasi-dynamic
# event here is truncated at nucleation rather than run to arrest.
par.exit_slip_rate = 1.5e-2

par.fric_sw_fs = 0
par.fric_sw_fd = 0
par.fric_sw_D0 = 0
par.fric_rsf_a, par.fric_rsf_b, par.fric_rsf_Dc = 0.007, 0.011, 11e-3
par.fric_rsf_deltaa = 0.01
par.fric_rsf_r0 = 0.6
par.fric_rsf_v0 = 1e-6

par.nfx = round((par.fxmax - par.fxmin)/par.dx + 1)
par.nfz = round((par.fzmax - par.fzmin)/par.dz + 1)
par.fx = np.linspace(par.fxmin, par.fxmax, par.nfx)
par.fz = np.linspace(par.fzmin, par.fzmax, par.nfz)
par.on_fault_vars = np.zeros((par.ntotft, par.nfz, par.nfx, 100))


def shear_steady_state(a, b, v0, r0, load_rate, norm, slip_rate, rou, vs):
    return -norm*a*asinh(slip_rate/2.0/v0*exp((r0 + b*log(v0/load_rate))/a)) \
        + rou*vs/2.0*slip_rate


# Velocity-weakening region: table 1 gives Lseis ~45 km along strike and
# Wseis 11.4 km down dip. Centred, that is |x| <= 22.5 km and z in
# [-16.7, -5.3] km, with a transition to velocity strengthening either side.
# Figure 1(d) and 1(e) show a and a-b tapering rather than stepping, which is
# what the transition widths below reproduce.
VW_HALF_LEN = 22.5e3
VW_TOP, VW_BOT = -5.3e3, -16.7e3
TRANS = 3.0e3          # taper width, along strike and down dip

for ix, xcoor in enumerate(par.fx):
    for iz, zcoor in enumerate(par.fz):
        dx_edge = VW_HALF_LEN - abs(xcoor)          # >0 inside, along strike
        dz_edge = min(zcoor - VW_BOT, VW_TOP - zcoor)  # >0 inside, down dip
        d = min(dx_edge, dz_edge)
        if d >= TRANS:
            a = par.fric_rsf_a                       # velocity weakening
        elif d <= 0.0:
            a = par.fric_rsf_a + par.fric_rsf_deltaa  # velocity strengthening
        else:
            a = par.fric_rsf_a + (1.0 - d/TRANS)*par.fric_rsf_deltaa

        par.on_fault_vars[0, iz, ix, 9]  = a
        par.on_fault_vars[0, iz, ix, 10] = par.fric_rsf_b
        par.on_fault_vars[0, iz, ix, 11] = par.fric_rsf_Dc
        par.on_fault_vars[0, iz, ix, 12] = par.fric_rsf_v0
        par.on_fault_vars[0, iz, ix, 13] = par.fric_rsf_r0
        par.on_fault_vars[0, iz, ix, 46] = par.creep_slip_rate
        par.on_fault_vars[0, iz, ix, 20] = \
            par.on_fault_vars[0, iz, ix, 11]/par.creep_slip_rate
        par.on_fault_vars[0, iz, ix, 7]  = par.init_norm
        # Steady state for THIS node's own initial slip rate. Passing the
        # creep rate for a node seeded faster leaves the patch short of the
        # stress its own rate needs, the Newton solve sides with the stress,
        # and the adaptive step blows up; that cost a day on the step-over.
        par.on_fault_vars[0, iz, ix, 8]  = shear_steady_state(
            a, par.fric_rsf_b, par.fric_rsf_v0, par.fric_rsf_r0,
            par.creep_slip_rate, par.init_norm,
            par.on_fault_vars[0, iz, ix, 46], par.rou, par.vs)

par.xmin_trans, par.xmax_trans = -25e3, 25e3
par.zmin_trans = -25e3
par.ymin_trans, par.ymax_trans = -5e3, 5e3
par.dx_trans = 50

par.casename = "liu2020-kink-qd-300"
par.HPC_nnode = 1
par.HPC_ncpu = 24
par.HPC_queue = "normal"
par.HPC_time = "04:00:00"
par.HPC_account = "EAR22013"
par.HPC_email = "dliu@ig.utexas.edu"

par.n_on_fault = 3
par.n_off_fault = 3
# Stations either side of the bend and on it, following the stars in
# figure 1(c).
par.st_coor_on_fault = [[-10.0, -9.0], [0.0, -9.0], [10.0, -9.0]]
par.st_coor_off_fault = [[0, 5, 0], [0, 5, -10], [-10, 10, 0]]

par.az_op = 2
par.az_maxiter = 2000
par.az_tol = 1.0e-7

# Liu, Duan & Luo (2020) section 3.5: normal stress caps of -100 to -10 MPa.
par.max_norm = -100.0e6

# Caps were historically implied by rough_fault; the switch is explicit now.
par.enforce_norm_caps = 1
