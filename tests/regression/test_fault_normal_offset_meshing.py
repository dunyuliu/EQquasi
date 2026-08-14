"""v1.5.0's uniform-dy belt landed only the MINIMUM fault y on a mesh line.

meshgen.f90 (and mesh4num.f90, which must count the same nodes) widened the
uniform-y belt in v1.5.0 to span every fault's y position, not just fault
1's -- but it still stepped a single dy from the minimum. A second fault
whose offset from that minimum was not an integer multiple of dy fell
between mesh lines and meshed with nftnd == 0: no fault nodes, an all-zero
slab in fault.*.nc, and the run completed anyway (4240 steps, 11 m of slip
on a single fault with no partner). The only symptom was "Fault nodes = 0"
in the run summary (fault 1's own count, not even the actual empty fault's
-- see the fix to that line too).

Owner's decision: do not reshape the mesh to fit an arbitrary offset.
build_yline_belt (src/func_lib.f90) refuses at setup, loudly, whenever a
fault's y-offset from the belt origin is not an integer multiple of dy --
before any expensive work. The user picks a commensurate dx/dy instead.

BP1002 (compsets/bp1002.qdc.2500) is exactly this: two faults offset by
5000 m. At dx = 2000 m (5000 / 2000 = 2.5) it must refuse; at dx = 2500 m
(5000 / 2500 = 2) it must mesh both faults.

This probes build_yline_belt directly -- the mesh-level unit, not a full
solver run -- so the guard is milliseconds and needs neither MPI nor MUMPS.
The compiled probe is driven twice (non-commensurate, then commensurate) so
both the refusal and the working case are guarded by the same test.
"""

import os
import shutil
import subprocess
import tempfile

import pytest

from conftest import ROOT

FIXTURES = ROOT / "tests" / "regression" / "fixtures"
FAULT_YS = [0.0, -5000.0]  # BP1002's two fault planes, meters.


def _mpif90():
    return shutil.which("mpif90") or shutil.which("mpif90.openmpi")


@pytest.fixture(scope="module")
def probe_exe():
    """Compile the probe once; tests below just re-run it with different stdin."""
    compiler = _mpif90()
    if not compiler:
        pytest.skip("no MPI Fortran wrapper on this host")

    build_dir = tempfile.mkdtemp(prefix="probe_yline_")

    def compile_one(src, obj):
        r = subprocess.run(
            [compiler, "-ffree-line-length-none", "-O0", "-g",
             "-c", str(src), "-o", str(obj), "-I", build_dir, "-J", build_dir],
            cwd=build_dir, capture_output=True, text=True)
        assert r.returncode == 0, f"compiling {src} failed:\n{r.stderr}"

    globalvar_o = os.path.join(build_dir, "globalvar.o")
    funclib_o = os.path.join(build_dir, "func_lib.o")
    driver_o = os.path.join(build_dir, "driver.o")
    exe = os.path.join(build_dir, "probe")

    compile_one(ROOT / "src" / "globalvar.f90", globalvar_o)
    compile_one(ROOT / "src" / "func_lib.f90", funclib_o)
    compile_one(FIXTURES / "probe_build_yline_belt.f90", driver_o)

    r = subprocess.run([compiler, globalvar_o, funclib_o, driver_o, "-o", exe],
                       cwd=build_dir, capture_output=True, text=True)
    assert r.returncode == 0, f"linking probe failed:\n{r.stderr}"

    yield exe
    shutil.rmtree(build_dir, ignore_errors=True)


def _run_probe(exe, dx, fault_ys):
    stdin = f"{len(fault_ys)} {dx}\n" + "\n".join(str(y) for y in fault_ys) + "\n"
    return subprocess.run([exe], input=stdin, capture_output=True, text=True)


def _parse_ylines(stdout):
    return [float(line.split()[-1]) for line in stdout.splitlines() if line.startswith("Y ")]


def test_noncommensurate_offset_refuses_loudly(probe_exe):
    """dx = 2000 m against a 5000 m offset (2.5x) must be REFUSED, not meshed."""
    r = _run_probe(probe_exe, 2000.0, FAULT_YS)
    assert r.returncode != 0, (
        f"dx = 2000 m against a 5000 m fault-normal offset (not an integer "
        f"multiple) meshed instead of refusing:\n{r.stdout}"
    )
    out = r.stdout + r.stderr
    assert "MESH ERROR" in out, f"refusal did not print a MESH ERROR box:\n{out}"
    # Named diagnostics: fault, its y, the belt origin, dy, and the ratio --
    # not just "something is wrong".
    assert "fault index" in out
    assert "fault y" in out
    assert "belt origin" in out
    assert "dy" in out
    assert "offset / dy" in out
    # No Y lines: build_yline_belt must stop before allocating/filling ylinet.
    assert not _parse_ylines(r.stdout), f"produced node lines despite refusing:\n{r.stdout}"


def test_commensurate_offset_meshes_both_faults(probe_exe):
    """dx = 2500 m against the same 5000 m offset (2x) must mesh cleanly."""
    r = _run_probe(probe_exe, 2500.0, FAULT_YS)
    assert r.returncode == 0, f"commensurate offset (5000 / 2500 = 2) was refused:\n{r.stdout}\n{r.stderr}"

    ys = _parse_ylines(r.stdout)
    assert ys, f"no node lines produced:\n{r.stdout}"

    tol = 2500.0 / 100.0
    for fy in FAULT_YS:
        hit = any(abs(y - fy) < tol for y in ys)
        assert hit, f"fault at y = {fy} m has no mesh line within {tol} m of it:\n{r.stdout}"

    # The belt is untouched (owner's design: refuse, don't reshape), so every
    # adjacent pair of lines inside the fault span is exactly dy apart --
    # not just close, exact, since nothing subdivides it.
    lo, hi = min(FAULT_YS), max(FAULT_YS)
    inside = sorted(y for y in ys if lo - tol <= y <= hi + tol)
    assert len(inside) >= 2
    steps = [round(b - a, 6) for a, b in zip(inside, inside[1:])]
    assert all(s == pytest.approx(2500.0, abs=1e-6) for s in steps), (
        f"belt is not uniform dy inside the fault span: {steps}"
    )
