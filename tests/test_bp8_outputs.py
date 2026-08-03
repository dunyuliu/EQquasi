"""Guards for BP8's required output files and their format.

Every assertion corresponds to a way the BP8 submission was actually
non-conforming before it was checked against the benchmark description.
"""

import re

import pytest

from conftest import read, strip_fortran_comments

LIBOUT = "src/library_output.f90"

# Section 4.3 of the BP8 description: slip, shear stress (both components) and
# pore pressure, along each of the two cross-section lines.
PROFILE_FILES = [
    f"{q}_{line}.dat"
    for q in ("slip_2", "slip_3", "shear_stress_2", "shear_stress_3", "pore_pressure")
    for line in ("strike", "depth")
]

# Section 4.1: the 11 fields, in column order.
STATION_FIELDS = [
    "t", "slip_2", "slip_3", "slip_rate_2", "slip_rate_3", "shear_stress_2",
    "shear_stress_3", "pore_pressure", "darcy_vel_2", "darcy_vel_3", "state",
]


@pytest.mark.parametrize("fname", PROFILE_FILES)
def test_section_43_profile_file_is_written(fname):
    """All ten were missing entirely: BP8 would have been rejected at upload."""
    src = strip_fortran_comments(read(LIBOUT))
    assert f"'{fname}'" in src, \
        f"{fname} is required by section 4.3 of the BP8 description but is never written"


def test_station_file_lists_all_eleven_fields_in_order():
    src = read(LIBOUT)
    m = re.search(r"'(t slip_2[^']*)'\s*//\s*&\s*\n\s*'([^']*)'", src)
    assert m, "could not find the BP8 field-name line"
    fields = (m.group(1) + m.group(2)).split()
    assert fields == STATION_FIELDS, \
        f"BP8 field list wrong.\n expected {STATION_FIELDS}\n got      {fields}"


def test_offfault_station_names_do_not_collide_for_metre_scale_models():
    """BP8 stations at 200 m and 400 m both became srfst_...st000... .

    output_offfault_st divides coordinates by 1000 and truncates. For a
    metre-scale model that maps distinct stations onto one filename, and the
    second silently overwrites the first.
    """
    src = strip_fortran_comments(read(LIBOUT))
    body = src.split("subroutine output_offfault_st")[1].split("end subroutine")[0]
    branch = re.search(r"if \(([^)]*bp == 7[^)]*)\) then", body)
    assert branch, "could not find the metre-scale branch in output_offfault_st"
    assert "bp == 8" in branch.group(1), (
        "bp == 8 must take the metre-scale filename branch; dividing metre "
        "coordinates by 1000 collapses distinct stations onto one file"
    )


def test_bp8_headers_carry_a_date():
    """The description lists date as required, unlike version or step counts."""
    src = strip_fortran_comments(read(LIBOUT))
    for block, label in (
        (src.split("subroutine output_onfault_st")[1].split("end subroutine")[0],
         "station time series"),
        (src.split("subroutine output_globaldat")[1].split("end subroutine")[0],
         "global.dat"),
        (src.split("subroutine write_one_profile")[1].split("end subroutine")[0],
         "section 4.3 profiles"),
    ):
        assert "# date=" in block, f"{label} header is missing the required date field"


def test_profile_rows_use_the_recommended_fortran_formats():
    """E22.14 for time, E15.7 for data, per the description's Fortran note."""
    src = strip_fortran_comments(read(LIBOUT))
    body = src.split("subroutine write_one_profile")[1].split("end subroutine")[0]
    assert "E22.14" in body and "E15.7" in body, \
        "profile writer should use E22.14 for time and E15.7 for data fields"


def test_profile_first_row_is_two_zeros_then_coordinates():
    """The description specifies [0, 0, x2...] as the first row."""
    src = strip_fortran_comments(read(LIBOUT))
    body = src.split("subroutine write_one_profile")[1].split("end subroutine")[0]
    assert re.search(r"0\.0d0,\s*0\.0d0,\s*\(coord", body), \
        "first profile row must be two zeros followed by the node coordinates"


def test_profile_time_rows_carry_log10_max_slip_rate():
    src = strip_fortran_comments(read(LIBOUT))
    body = src.split("subroutine write_one_profile")[1].split("end subroutine")[0]
    assert "dlog10" in body and "profVmax" in body, \
        "each profile row must carry log10 of the max slip rate as column 2"


# --- run provenance -------------------------------------------------------

REQUIRED_RUNINFO_KEYS = [
    "code", "version", "benchmark_id", "run_timestamp", "host",
    "host_logical_cores", "cpu_model", "mpi_ranks", "omp_threads_per_rank",
    "num_nodes", "num_elements", "num_fault_nodes", "num_equations",
    "element_size_m", "steps_completed", "simulated_time_s",
    "time_loop_seconds", "factorization_seconds", "seconds_per_step",
]


@pytest.mark.parametrize("key", REQUIRED_RUNINFO_KEYS)
def test_runinfo_json_records_key(key):
    """A timing without its core count and machine is not reproducible."""
    src = strip_fortran_comments(read(LIBOUT))
    body = src.split("subroutine output_run_metadata")[1].split("end subroutine")[0]
    assert f'"{key}"' in body, f"runInfo.json must record {key}"


def test_runinfo_hostname_does_not_rely_on_the_environment_alone():
    """HOSTNAME is not exported to non-interactive shells; it read 'unknown'."""
    src = strip_fortran_comments(read(LIBOUT))
    body = src.split("subroutine output_run_metadata")[1].split("end subroutine")[0]
    assert "/proc/sys/kernel/hostname" in body, \
        "fall back to procfs when HOSTNAME is unset, or the host records as unknown"


def test_run_summary_is_printed_to_the_log():
    src = strip_fortran_comments(read("src/solveTimeLoopMUMPS.f90"))
    assert "RUN SUMMARY" in src, "the run should print a summary block to stdout"
    for token in ("MPI ranks", "Elements", "Seconds per step", "Time steps completed"):
        assert token in src, f"run summary should report {token}"
