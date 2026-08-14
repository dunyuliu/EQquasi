# Aztec solver path (dormant, not broken)

`main_aztec.f90`, `elemcal_aztec.f90` and `msr.f90` implement an iterative-solver
path through [AZTEC 2.1](https://trilinos.github.io/aztecoo.html), reached via
`sol_op == 2`.

**This worked.** It was a deliberate probe of an iterative solver, not an
abandoned or failed experiment. Production runs use MUMPS instead because MUMPS,
a direct sparse solver, is faster on problems of the size EQquasi actually runs --
iterative methods pay off at much larger scale than these meshes reach. So the
Aztec path was left in place rather than deleted, and gradually stopped being
built.

By v1.4.8 the state was:

  - `src/eqquasi.f90` printed "aztec is temporarily disabled" for `sol_op == 2`,
    with `call main_aztec` commented out;
  - the makefile rules for `main_aztec.o` and `elemcal_aztec.o` were commented
    out, so neither file had been compiled for some time;
  - `msr.f90` was still in `OBJ`, so it was compiled and linked into every
    binary even though its only caller, `main_aztec.f90:46`, was never built.

Moved here so `src/` contains only what ships. Nothing was reachable, so nothing
changed. If the iterative path is ever wanted again -- a much larger mesh, or a
machine where MUMPS memory becomes the binding constraint -- this is the starting
point, and it is known to have run rather than being a sketch.

`AZTEC_OPTIONS` deliberately stays in `src/globalvar.f90` and `src/read_input.f90`:
`model.txt` is a *positional* file written by `script/case.setup`, so removing
the field would shift every value after it and silently corrupt the parameters of
every case. It costs one integer.
