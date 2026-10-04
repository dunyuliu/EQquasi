# Parameters

A case is configured by one file, `user_defined_params.py`, in the case
root. It starts from the default `parameters` class and overrides what the
case needs:

```python
from defaultParameters import parameters
par = parameters()
par.dx = 1000.0e0
par.nstep = 5000
```

`case.setup` reads `par` and writes the solver's input files into `input/`.
The reference below is generated from the defaults class itself, so it
lists every parameter that exists and its current default.

## Things the reference cannot show

**Per-node fault properties.** Rate-and-state `a`, `b`, `Dc`, initial slip
rate, state, and initial normal and shear stress are set node by node in
the array `par.on_fault_vars`, shaped `(ntotft, nfz, nfx, 100)`. The
default class fills it with the BP5 distribution; a compset that needs
different properties fills its own. The fourth index selects the quantity;
the named slots (`FR_RSF_A`, `FR_RSF_B`, `FR_RSF_DC`, `FR_VINIT`,
`FR_STATE`, `FR_TNRM0`, `FR_TSTK0`, ...) are defined, with units, at the top
of `script/defaultParameters.py`. Import them from there rather than writing
raw numbers.

**Benchmark-specific behaviour.** `bp` selects behaviour built into the
solver. Values with special handling are `5` (SEAS BP5), `7` (SEAS BP7: the
prescribed nucleation perturbation), `8` (SEAS BP8: pore-fluid diffusion,
`.dat` output, profile files) and `1001`. Any other value runs the generic
path; some compsets use one as a label only.

**Multiple faults.** Set `ntotft` and give `faultgeom` one
`(xlo, xhi, ycoor, zlo, zhi)` tuple per fault, in metres. Each fault lies on
a constant-y plane. Every fault's y offset must be a whole number of `dy`
steps, and its x and z bounds whole numbers of `dx` steps from the mesh
origin. The solver refuses a geometry that breaks this (see
[Troubleshooting](troubleshooting.md)).

**Normal-stress caps.** `min_norm` and `max_norm` take effect only when
`C_normal_stress_caps = 1`. They are off by default. `case.setup` warns when
a non-planar or multi-fault case runs without them.

<!-- BEGIN PARAMETER REFERENCE (generated from script/defaultParameters.py by docs/user/gen_params.py; do not edit by hand) -->

Every entry below is an attribute of the `parameters` class in `script/defaultParameters.py`. A case overrides any of them in its own `user_defined_params.py` as `par.<name> = <value>`. The defaults are the class's own, taken from BP5 at 2000 m; each compset sets its own values on top.

* **`istart`** -- default `1`

  cylce id. Simulate quasi-dynamic earthquake cycles from istart to iend.

* **`iend`** -- default `1`

* **`mode`** -- default `1`

  mode of the code - quasi-dynamic (1) or fully-dynamic (2).

* **`dip`** -- default `90.0`

* **`fxmin`**, **`fxmax`** -- defaults `-60000.0`, `60000.0`

  model_domain (in meters)

* **`fymin`**, **`fymax`** -- defaults `-50000.0`, `50000.0`

* **`fzmin`**, **`fzmax`** -- defaults `-60000.0`, `0.0`

* **`xminc`**, **`xmaxc`**, **`zminc`** -- defaults `-50000.0`, `50000.0`, `-40000.0`

  creeping zone bounaries. creeping zones are assinged on the lateral sides and bottom of the RSF controlled region and will slide at fixed loading slip rate.

* **`dx`** -- default `2000.0`

  cell size, spatial resolution

* **`dy`** -- default `dx`

* **`dz`** -- default `dx`

* **`nuni_y_plus`**, **`nuni_y_minus`** -- defaults `5`, `5`

  along the fault-normal dimension, the number of cells share the dx cell size.

* **`enlarging_ratio`** -- default `1.3`

  along the fault-normal dimension (y), cell size will be enlarged at this ratio compoundly.

* **`enlarging_ratio_xz`** -- default `1.0`

  along strike (x) and dip (z), outside the fault box (par.faultgeom, else the domain): cell size grows at this ratio compoundly, capped at min(12*dx, 3 km). 1 = uniform (pre-1.21).

* **`vp`**, **`vs`**, **`rou`** -- defaults `6000.0`, `3464.0`, `2670.0`

  Isotropic material propterty. Vp, Vs, Rou

* **`init_norm`** -- default `-25000000.0`

  initial normal stress in Pa. Negative compressive.

* **`insertFaultType`** -- default `0`

  Controlling switches for EQquasi system 0: no fault, 1: planar 2: rough.

* **`rough_fault`** -- default `insertFaultType`

* **`rheology`** -- default `1`

  elastic(1).

* **`friclaw`** -- default `3`

  rsf_aging(3), rsf_slip(4).

* **`ntotft`** -- default `1`

  number of total faults.

* **`faultgeom`** -- default `None`

  Optional per-fault geometry override, one 5-tuple per fault: (xlo, xhi, ycoor, zlo, zhi) in meters. None (default) means every fault uses the single domain-derived box (ycoor = 0) that case.setup has always emitted; set this to a list of len(ntotft) tuples to place faults at different y (e.g. a step-over).

* **`solver`** -- default `1`

  solver option. MUMPS(1, recommended). PETSc(2).

* **`nstep`** -- default `10000`

  total num of time steps for exiting, if not exit via sliprate threshold

* **`nt_out`** -- default `100`

  Every nt_out time steps, disp of the whole model and on-fault variables will be written out in netCDF format.

* **`bp`** -- default `5`

  currently supported cases 5 (SCEC-BP5) 1001 (GM-cycle)

* **`xi`** -- default `0.2`

  xi, minimum Dc Lapusta et al. (2009) time-step factor, dtev = xi*Dc/Vmax. Measured for BP8 at dx = 50 m: 0.05, 0.1 and 0.2 cost the same wall clock for a given step count but reach 6.0, 17.0 and 22.3 days, agreeing to 0.003 log units in peak slip rate. 0.2 is the default; tighten it per compset if a problem needs it.

* **`minDc`** -- default `0.13`

  meters

* **`far_vel_load`** -- default `4e-10`

  loading far field shear loading velocity (+x vel=far_vel_load on xz bound at y=ymax). A minus value is applied on the other side.

* **`far_norm_load_vel`** -- default `0.0`

  far field fault-normal extensional loading velocity (+y vel=far_norm_load_vel on xz bound at y=ymax).

* **`creep_slip_rate`** -- default `1e-09`

  creeping slip rate outside of RSF controlled region.

* **`exit_slip_rate`** -- default `0.001`

  exiting slip rate for EQquasi [m/s].

## Frictional variables

* **`fric_sw_fs`** -- default `0`

  friclaw == 1, slip weakening

* **`fric_sw_fd`** -- default `0`

* **`fric_sw_D0`** -- default `0`

* **`fric_rsf_a`**, **`fric_rsf_b`**, **`fric_rsf_Dc`** -- defaults `0.004`, `0.03`, `0.14`

  friclaw == 3, rate- and state- friction with aging law.

* **`fric_rsf_deltaa`** -- default `0.036`

* **`fric_rsf_r0`** -- default `0.6`

* **`fric_rsf_v0`** -- default `1e-06`

## Domain boundaries for transferring

* **`xmin_trans`**, **`xmax_trans`** -- defaults `-25000.0`, `25000.0`

* **`zmin_trans`** -- default `-25000.0`

* **`ymin_trans`**, **`ymax_trans`** -- defaults `-5000.0`, `5000.0`

* **`dx_trans`** -- default `50`

## Along-fault pore fluid diffusion, used when bp == 8

* **`fluid_src`** -- default `0`

  Off by default, so every other compset is unaffected. 0: off; 1: Gaussian source (GS); 2: Peaceman well (PW).

* **`fluid_q0`** -- default `0.0`

  total volume injection rate, m^3/s.

* **`fluid_toff`** -- default `0.0`

  injection turn-off time, s.

* **`fluid_tend`** -- default `0.0`

  final simulation time, s. 0 disables the time-based exit.

* **`fluid_Lgauss`** -- default `50.0`

  characteristic size of the Gaussian source, m.

* **`fluid_Lfwid`** -- default `1.0`

  fault zone thickness, m.

* **`fluid_beta`** -- default `1e-08`

  pore and fluid compressibility, 1/Pa.

* **`fluid_phi`** -- default `0.1`

  porosity.

* **`fluid_perm`** -- default `5e-14`

  permeability, m^2.

* **`fluid_eta`** -- default `0.001`

  fluid viscosity, Pa s.

* **`fluid_Swell`** -- default `1e-07`

  volumetric well storage, m^3/Pa.

* **`fluid_rwell`** -- default `0.05`

  true well radius, m. BP8 Table 1 (2026-08-12 revision).

* **`dtmax`** -- default `0.0`

  cap on the adaptive time step, s. 0 means no cap.

* **`fric_pc_L`** -- default `0.0`

  Prakash-Clifton relaxation distance for the normal-stress state variable, m. 0.0 : default. No state; use the instantaneous effective normal stress. > 0  : relax the normal-stress state over this slip distance.

* **`min_norm`** -- default `-10000000.0`

* **`max_norm`** -- default `-40000000.0`

* **`C_normal_stress_caps`** -- default `0`

  ONE explicit switch, no gates: min_norm/max_norm apply iff this is 1. OFF by default -- most cases never need caps. A case that wants them sets this in its own user_defined_params.py; nothing else turns them on, and nothing else turns them off.

## HPC resource allocation

* **`casename`** -- default `'bp5-qd-2000'`

* **`HPC_nnode`** -- default `1`

  Number of computing nodes. On LS6, one node has 128 CPUs.

* **`HPC_ncpu`** -- default `3`

  Number of CPUs requested.

* **`HPC_queue`** -- default `'normal'`

  q status. Depending on systems, job WALLTIME and Node requested.

* **`HPC_time`** -- default `'00:30:00'`

  WALLTIME, in hh:mm:ss format.

* **`HPC_account`** -- default `'EAR22013'`

  Project account to be charged SUs against.

* **`HPC_email`** -- default `'dliu@ig.utexas.edu'`

  Email to receive job status.

## Single station time series output

* **`st_coor_on_fault`**

  Default:

  ```python
  [[-36.0, 0.0], [-16.0, 0.0], [0.0, 0.0], [16.0, 0.0], [36.0, 0.0], [-24.0, 0.0], [-16.0, 0.0], [0.0, -10.0], [16.0, -10.0], [0.0, -22.0]]
  ```

  (x,z) coordinate pairs for on-fault stations (in km).

* **`st_coor_off_fault`**

  Default:

  ```python
  [[0, 8, 0], [0, 8, -10], [0, 16, 0], [0, 16, -10], [0, 32, 0], [0, 32, -10], [0, 48, 0], [16, 8, 0], [-16, 8, 0]]
  ```

  (x,y,z) coordinates for off-fault stations (in km).

* **`az_op`** -- default `2`

  Additional solver options for AZTEC AZTEC options

* **`az_maxiter`** -- default `2000`

  maximum iteration for AZTEC

* **`az_tol`** -- default `1e-07`

  tolerance for solution in AZTEC.

<!-- END PARAMETER REFERENCE -->
