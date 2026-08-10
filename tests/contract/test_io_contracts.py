"""Guards for the two undocumented contracts that couple Python to Fortran.

1. model.txt is positional and unlabelled: scripts/case.setup writes N lines
   and src/read_input.f90 reads them back in exactly that order. Inserting a
   line in one place and not the other silently shifts every later parameter.
   PROJECT_RULES.md rule 4.

2. The on-fault variables travel through on_fault_vars_input.nc by NAME on the
   netCDF side but land in fric(<magic index>, node, fault) on the Fortran
   side. The index mapping is duplicated with no single source of truth.
   PROJECT_RULES.md rule 5.
"""

import re

from conftest import read, strip_fortran_comments, strip_python_comments


def model_txt_write_count():
    """Number of lines scripts/case.setup writes into model.txt."""
    src = strip_python_comments(read("scripts/case.setup"))
    body = src.split("def create_model_input_file")[1].split("\ndef ")[0]
    return len(re.findall(r'f\.write\(', body))


def model_txt_read_count():
    """Number of records src/read_input.f90 reads from unit 1002."""
    src = strip_fortran_comments(read("src/read_input.f90"))
    body = src.split("file = 'model.txt'")[1].split("close(1002)")[0]
    return len(re.findall(r"read\(1002,\*", body))


def test_model_txt_write_and_read_counts_agree():
    w, r = model_txt_write_count(), model_txt_read_count()
    assert w == r, (
        f"scripts/case.setup writes {w} model.txt lines but "
        f"src/read_input.f90 reads {r}. model.txt is positional: a mismatch "
        "shifts every subsequent parameter silently."
    )


def test_netcdf_on_fault_variable_names_agree_between_writer_and_reader():
    """case.setup creates the variables; netcdf_io.f90 looks them up by name."""
    py = strip_python_comments(read("scripts/case.setup"))
    written = set(re.findall(r"createVariable\(\s*'([^']+)'", py))
    written -= {"dip", "strike"}  # coordinate variables, not on-fault fields

    f90 = strip_fortran_comments(read("src/netcdf_io.f90"))
    body = f90.split("subroutine netcdf_read_on_fault(")[1].split("end subroutine")[0]
    read_names = set(re.findall(r'nf90_inq_varid\(ncid,\s*"([^"]+)"', body))

    assert written == read_names, (
        "on_fault_vars_input.nc variable names disagree.\n"
        f"  written by case.setup but not read: {sorted(written - read_names)}\n"
        f"  read by netcdf_io.f90 but not written: {sorted(read_names - written)}"
    )


def test_fric_indices_used_by_the_reader_are_within_the_allocated_array():
    """fric is allocated as fric(100, nftmx, ntotft)."""
    alloc = read("src/eqquasi.f90")
    m = re.search(r"fric\((\d+),\s*nftmx", alloc)
    assert m, "could not find the fric allocation"
    nfric = int(m.group(1))

    used = set()
    for path in ("src/faulting.f90", "src/netcdf_io.f90", "src/library_output.f90",
                 "src/porepressure.f90", "src/meshgen.f90"):
        try:
            src = strip_fortran_comments(read(path))
        except FileNotFoundError:
            continue
        used.update(int(n) for n in re.findall(r"\bfric\((\d+)\s*,", src))

    over = sorted(i for i in used if i > nfric)
    assert not over, f"fric() indices {over} exceed the allocated first dimension {nfric}"


def test_fltsta_indices_are_within_the_allocated_array():
    """fltsta is allocated as fltsta(<n>, nstep, nonmx) and indexed by magic number."""
    alloc = read("src/eqquasi.f90")
    m = re.search(r"fltsta\((\d+),\s*nstep", alloc)
    assert m, "could not find the fltsta allocation"
    n = int(m.group(1))

    used = set()
    for path in ("src/faulting.f90", "src/library_output.f90"):
        src = strip_fortran_comments(read(path))
        used.update(int(x) for x in re.findall(r"\bfltsta\((\d+)\s*,", src))

    over = sorted(i for i in used if i > n)
    assert not over, (
        f"fltsta() indices {over} exceed the allocated first dimension {n}. "
        "Writing past it is out-of-bounds; reading past it returns garbage."
    )


def test_globaldat_indices_are_within_the_allocated_array():
    alloc = read("src/eqquasi.f90")
    m = re.search(r"globaldat\((\d+),\s*nstep", alloc)
    assert m, "could not find the globaldat allocation"
    n = int(m.group(1))

    used = set()
    for path in ("src/solveTimeLoopMUMPS.f90", "src/library_output.f90"):
        src = strip_fortran_comments(read(path))
        used.update(int(x) for x in re.findall(r"\bglobaldat\((\d+)\s*,", src))

    over = sorted(i for i in used if i > n)
    assert not over, \
        f"globaldat() indices {over} exceed the allocated first dimension {n}"


def test_all_faulting_accumulators_are_initialised():
    """Uninitialised stack memory leaked into global.dat and the pma printout.

    faulting declares seven automatic accumulator arrays but originally zeroed
    only three. momrate_arr, momRateVW and ma_bar_ku_arr are written only inside
    the frictional branch, so every fault node in the creeping / no-slip region
    contributed whatever was on the stack to sum(momrate_arr) and
    maxval(ma_bar_ku_arr). The symptom was a constant ~-3.1e304 in global.dat
    column 3, printed as '-0.3145234+305' because E15.7 cannot hold a three
    digit exponent.

    The regression oracle cannot catch this: it only compares fault.*.nc.
    """
    src = strip_fortran_comments(read("src/faulting.f90"))
    decl = re.search(r"real \(kind = dp\) :: (ma_bar_ku_arr.*?)\n\n", src, re.S)
    assert decl, "could not find the accumulator declaration block"
    declared = set(re.findall(r"(\w+)\(nftnd\(1\)\)", decl.group(1)))
    assert declared, "no accumulator arrays parsed"

    body = src.split("real (kind = dp) :: ma_bar_ku_arr")[1]
    initialised = set(re.findall(r"^\s*(\w+)\s*=\s*0\.0d0\s*$", body, re.M))

    missing = declared - initialised
    assert not missing, (
        f"these faulting accumulators are never zeroed: {sorted(missing)}. "
        "Any fault node that skips the frictional branch will contribute "
        "uninitialised memory to the totals."
    )


def test_benchmark_file_parsing_is_not_reimplemented_per_script():
    """One parser for benchmark output files, in scripts/seasio.py.

    Six scripts had each grown their own "strip, skip #, try float(tok[0])" loop.
    They agreed on the load-bearing subtlety -- a data row is one whose FIRST
    token parses as a number, because field names like slip_2 and x2 contain
    digits -- but six copies is six chances to get it wrong, and one of them
    (plotDomainSweep.py) had already drifted to globbing the wrong extension.
    """
    import glob
    import os
    import re
    from conftest import ROOT

    offenders = []
    for path in glob.glob(os.path.join(str(ROOT), "scripts", "*")):
        if os.path.isdir(path) or os.path.basename(path) == "seasio.py":
            continue
        try:
            src = open(path).read()
        except (UnicodeDecodeError, IsADirectoryError):
            continue
        if "float(tok[0])" in src or re.search(r"float\([a-z]+\.split\(\)\[0\]\)", src):
            offenders.append(os.path.basename(path))
    assert not offenders, (
        f"{offenders} reimplement benchmark-file parsing; import read_table, "
        f"read_rows, read_array or read_profile from seasio instead."
    )
