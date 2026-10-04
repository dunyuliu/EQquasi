# Output Files

Each completed cycle's output is in `result/cycle<N>/` of the case;
`log/cycle<N>.log` holds its screen output. Nothing is written to `input/`
after `case.setup`. This page describes each file family. For how the files
are compared against reference results, see [Benchmarks](benchmarks.md).

## Run record: `runInfo.json`

One JSON object describing the run: EQquasi `version`, `benchmark_id`
(`bp`), `run_timestamp`, `host`, `cpu_model`, `mpi_ranks`,
`omp_threads_per_rank`, mesh size (`num_nodes`, `num_elements`,
`num_fault_nodes`, `num_equations`, `element_size_m`), `steps_completed`,
`simulated_time_s`, timings (`time_loop_seconds`, `factorization_seconds`,
`seconds_per_step`), `max_slip_rate_final_m_s`, `solver`, `friclaw` and
`exit_reason`. Before comparing two runs, check that both were made by the
same `version`. `omp_threads_per_rank` is 0 when `OMP_NUM_THREADS` was not
set.

## Whole-fault time series: `global.dat`

One row per time step, 7 columns:

| column | quantity | unit |
|---|---|---|
| 1 | time | s |
| 2 | peak slip rate over all faults | m/s |
| 3 | moment rate, sum of mu V dx^2 | N m/s |
| 4 | shear traction times area, over nodes slipping faster than `exit_slip_rate` | N |
| 5 | slip times area, over the same nodes | m^3 |
| 6 | area of the same nodes | m^2 |
| 7 | moment rate within 200 m of the origin (BP7 only, otherwise 0) | N m/s |

BP8 writes the SEAS section 4.2 format instead: a commented header and
three columns, time, log10 of peak slip rate, and moment rate, starting
with the initial condition at t = 0.

`peak_sliprate_per_fault.dat` has the same time column followed by one
peak-slip-rate column per fault (1 + `ntotft` columns, m/s). The largest of
those columns equals `global.dat` column 2.

## On-fault stations: `fltst_strk<x>dp<z>.txt`

One file per on-fault station in `st_coor_on_fault`; `<x>` and `<z>` are
the station's along-strike position and depth in km. One row per time step,
9 columns:

| column | quantity | unit |
|---|---|---|
| 1 | time | s |
| 2 | slip along strike | m |
| 3 | slip along dip, positive downward | m |
| 4 | slip rate magnitude | m/s |
| 5 | slip rate along dip | m/s |
| 6 | shear stress along strike | MPa |
| 7 | shear stress along dip, sign reversed | MPa |
| 8 | effective normal stress, negative in compression | MPa |
| 9 | log10 of the state variable | log10 s |

For BP7 the station names are in metres rather than km. For BP8 the names
carry an explicit sign and use `.dat` (`fltst_strk+000dp-200.dat`), and
the file follows SEAS section 4.1: a commented header and 11 columns -- time,
slip along strike and dip, log10 slip rate along strike and dip, shear
stress along strike and dip (MPa), pore pressure (MPa), Darcy velocity
along strike and dip (m/s), and log10 state -- starting from t = 0.

A station is written only if its (x, z) lands on a fault node: both must be
multiples of `dx` measured from the fault's corner. The solver prints
`ON-FAULT STATIONS REQUESTED BUT NOT FOUND` when some do not.

## Off-fault stations: `srfst_strk<x>st<y>dp<z>.txt`

One file per station in `st_coor_off_fault`; `<x>`, `<y>` and `<z>` are
its along-strike, fault-normal and depth coordinates (km; metres for BP7 and
BP8). 7 columns: time (s), then displacement (m) and velocity (m/s) in turn
along strike, vertical, and fault-normal.

## Fault state at the end of the cycle: `cplot_EQquasi.txt`

One row per fault node, 16 columns: x (m), depth (m), slip rate (m/s),
state variable, shear stress along strike and dip (Pa), effective normal
stress (Pa), the x, y and z velocities of the node on the +y side and then
the -y side of the fault (m/s), two unused columns, and the node's rupture
time (s) -- the time its slip rate first reached 1 mm/s.
`plotRuptureTime.py` reads this file.

`cplot_ruptarea_trac_slip.txt` has one row per fault node: ruptured area
(m^2), slip during the rupture (m), and shear traction at the start and end
of the rupture (Pa). `plotMagnitudeTime.py` reads it to compute the moment
magnitude.

`tdyna.txt` holds two numbers: the times (s) at which the cycle entered and
left its coseismic phase.

## Snapshots: `fault.NNNNN.nc` and `disp.NNNNN.nc`

netCDF snapshots written every `nt_out` steps and at the last step, named
by step number. `fault.NNNNN.nc` holds the fault-plane fields on a
`(nid_fault, nid_dip, nid_strike)` grid: `shear_strike`, `shear_dip`,
`effective_normal` (Pa), `slip_rate` (m/s), `state_variable`,
`state_normal`, the nodal velocities `vxm vym vzm vxs vys vzs`, and
accumulated slip `slips`, `slipd`, `slipn` (m). `disp.NNNNN.nc` holds the
displacement of every node in the model. `plotOnFaultVars` maps every
`fault.NNNNN.nc`.

## Restart files: `fault.r.nc` and `disp.r.nc`

The fault and displacement state at the end of the cycle, in the same
formats. The next cycle starts from them; `plotAccumulated` reads slip from
`fault.r.nc`.

## Mesh: `eqquasi.mesh.coor.nc`, `eqquasi.mesh.ien.nc`, `eqquasi.mesh.nsmp.nc`

Node coordinates, element connectivity, and the fault's paired-node list,
written once per cycle.

## BP5 and BP7 profiles: `p1output.txt`, `p2output.txt`

For BP5 and BP7 only: time, slip along strike and dip, and shear stress
along strike and dip, at the nodes of one along-strike line (`p1output.txt`)
and one down-dip line (`p2output.txt`), appended every 20 steps between
earthquakes and every 10 steps during one.

## BP8 profiles: `<quantity>_strike.dat`, `<quantity>_depth.dat`

For BP8 only, the ten SEAS section 4.3 files, along the strike line and the
depth line through the injection point:

* `slip_2_strike.dat`, `slip_2_depth.dat` -- slip along strike (m)
* `slip_3_strike.dat`, `slip_3_depth.dat` -- slip along dip (m)
* `shear_stress_2_strike.dat`, `shear_stress_2_depth.dat` -- shear stress along strike (MPa)
* `shear_stress_3_strike.dat`, `shear_stress_3_depth.dat` -- shear stress along dip (MPa)
* `pore_pressure_strike.dat`, `pore_pressure_depth.dat` -- pore pressure (MPa)

 In each, row 1 holds the node coordinates; each later row holds the
time, log10 of the peak slip rate, and the value at every node.

## Cycle number: `currentcycle.txt`

The cycle number `run.sh` passed to the solver. Cycle 1 starts from the
initial conditions in `input/`; later cycles restart from the previous
cycle's `fault.r.nc` and `disp.r.nc`.

## Multiple faults

A case with `ntotft > 1` tags the per-fault files of fault 2 onward with
`ft<N>_`, and leaves fault 1's names as above:

* `fltst_ft2_strk<x>dp<z>.txt` -- on-fault stations, with coordinates local to that fault.
* `cplot_ft2_EQquasi.txt` and `cplot_ft2_ruptarea_trac_slip.txt`.

`global.dat` and the netCDF files are not split: the snapshots carry a
`nid_fault` dimension, and `peak_sliprate_per_fault.dat` has one column per
fault.
