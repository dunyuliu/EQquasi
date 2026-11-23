"""MUMPS workspace ICNTL(14) = 100 is set in all three solver paths.

SOURCE-ORDER CHECK, NOT A RUN CHECK. It reads src/solveTimeLoopMUMPS.f90,
solveTimeLoopPETSc.f90 and solveTimeLoopPETScDD.f90 (comments stripped) and
asserts the setting exists, with value 100, at the place it can take effect:

* MUMPS path: `mumps_par%ICNTL(14) = 100` after the JOB=-1 `call DMUMPS`
  (init resets ICNTL to defaults) and before the JOB=4 factorization call.
* PETSc paths: `-mat_mumps_icntl_14` is set to '100' after
  PCFactorSetMatSolverType and before KSPSetFromOptions (which reads it),
  and KSPSetUp (the factorization) comes after.

It cannot prove the factorization used 100: a run check needs the solver to
print INFOG/ICNTL(14) (e.g. mumps_par%ICNTL(14) after JOB=4, or
-mat_mumps_icntl_14 echoed by -mat_view/-ksp_view) and a test to grep it.
That is a src/ change, so it is the owner's decision.
"""

import re

import pytest

from conftest import read, strip_fortran_comments


def _pos(code, pattern, after=0):
    m = re.compile(pattern, re.I).search(code, after)
    assert m, f"pattern not found after offset {after}: {pattern}"
    return m.start()


def test_mumps_path_sets_icntl14_after_init_and_before_factorization():
    code = strip_fortran_comments(read("src/solveTimeLoopMUMPS.f90"))
    init_job = _pos(code, r"mumps_par%JOB\s*=\s*-1\b")
    init_call = _pos(code, r"call\s+DMUMPS\s*\(\s*mumps_par\s*\)", init_job)
    setting = _pos(code, r"mumps_par%ICNTL\(14\)\s*=\s*100\b", init_call)
    fact_job = _pos(code, r"mumps_par%JOB\s*=\s*4\b", init_call)
    fact_call = _pos(code, r"call\s+DMUMPS\s*\(\s*mumps_par\s*\)", fact_job)
    assert init_call < setting < fact_job < fact_call, (
        f"order init {init_call}, ICNTL(14)=100 {setting}, JOB=4 {fact_job}, factorize {fact_call}")


@pytest.mark.parametrize("path,setup", [
    ("src/solveTimeLoopPETSc.f90", r"call\s+KSPSetUp\s*\(\s*ksp"),
    ("src/solveTimeLoopPETScDD.f90", r"call\s+KSPSetUp\s*\(\s*ksp"),
])
def test_petsc_path_sets_icntl14_100_between_solver_type_and_setup(path, setup):
    code = strip_fortran_comments(read(path))
    solver = _pos(code, r"call\s+PCFactorSetMatSolverType\s*\(\s*pc\s*,\s*MATSOLVERMUMPS")
    setting = _pos(code, r"call\s+PetscOptionsSetValue\s*\([^)]*'-mat_mumps_icntl_14'\s*,\s*'100'", solver)
    from_opts = _pos(code, r"call\s+KSPSetFromOptions\s*\(\s*ksp", setting)
    setup_at = _pos(code, setup, from_opts)
    assert solver < setting < from_opts < setup_at, (
        f"{path}: order solver type {solver}, icntl_14=100 {setting}, "
        f"KSPSetFromOptions {from_opts}, KSPSetUp {setup_at}")
