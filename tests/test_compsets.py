"""Guards for compset parameter files.

The defect that motivated most of this file: case_input/bp7.qdc.a.10 sets
par.xmin/par.xmax/... but scripts/case.setup reads par.fxmin/par.fxmax/...
Assigning the wrong name silently leaves BP5's default 120x100x60 km domain in
place, so BP7 wrote a mesh many orders of magnitude too large and was
OOM-killed. Nothing in the repo noticed.
"""

import pytest

from conftest import compset_dirs, par_attributes_read_by, read, strip_python_comments

# Attributes case.setup reads that must resolve to something the compset
# actually meant to set, rather than an inherited default.
DOMAIN_ATTRS = ["fxmin", "fxmax", "fymin", "fymax", "fzmin", "fzmax"]

# Wrong-name variants that look right but are silently ignored.
SHADOW_ATTRS = {"xmin", "xmax", "ymin", "ymax", "zmin", "zmax"}


def compset_ids():
    return [d.name for d in compset_dirs()]


@pytest.mark.parametrize("name", compset_ids())
def test_compset_does_not_set_shadow_domain_attributes(name):
    """`par.xmin = ...` creates a new unused attribute; case.setup reads fxmin."""
    src = strip_python_comments(read(f"case_input/{name}/user_defined_params.py"))
    import re as _re
    assigned = set()
    for line in src.splitlines():
        if "=" not in line:
            continue
        lhs = line.split("=")[0]
        for attr in SHADOW_ATTRS:
            # Word boundary: par.xminc must not count as par.xmin.
            if _re.search(rf"\bpar\.{attr}\b", lhs):
                assigned.add(attr)
    assert not assigned, (
        f"{name} sets par.{sorted(assigned)}, which case.setup never reads. "
        "Use the f-prefixed names (fxmin, fxmax, fymin, fymax, fzmin, fzmax) "
        "or the compset silently inherits defaultParameters.py's domain."
    )


@pytest.mark.parametrize("name", compset_ids())
def test_compset_sets_dy_and_dz_when_it_overrides_dx(name):
    """defaultParameters evaluates dy = dx once, at class-definition time.

    A compset that overrides par.dx therefore leaves dy/dz at the default 2000,
    which corrupts case.setup's estimate_HPC_resource() (it divides by par.dz).
    """
    src = strip_python_comments(read(f"case_input/{name}/user_defined_params.py"))
    if "par.dx" not in src:
        return
    assert "par.dy" in src and "par.dz" in src, (
        f"{name} overrides par.dx but not par.dy/par.dz, which stay at the "
        "class-level default and desynchronise from dx"
    )


def test_case_setup_reads_only_attributes_that_exist_on_the_defaults():
    """Every par.<attr> case.setup touches must exist in defaultParameters."""
    import ast

    needed = par_attributes_read_by("scripts/case.setup")
    tree = ast.parse(read("scripts/defaultParameters.py"))
    cls = next(n for n in ast.walk(tree)
               if isinstance(n, ast.ClassDef) and n.name == "parameters")
    defined = set()
    for node in ast.walk(cls):
        if isinstance(node, ast.Assign):
            for t in node.targets:
                for sub in ast.walk(t):
                    if isinstance(sub, ast.Name):
                        defined.add(sub.id)

    missing = needed - defined
    assert not missing, (
        f"case.setup reads par.{sorted(missing)}, which defaultParameters.py "
        "does not define; any compset that forgets to set them crashes or, "
        "worse, silently uses a stale value"
    )


def test_compsets_txt_matches_case_input():
    """compsets.txt is documentation-only and can drift silently (rule 7)."""
    listed = {ln.strip() for ln in read("case_input/compsets.txt").splitlines()
              if ln.strip()}
    actual = {d.name for d in compset_dirs() if not d.name.startswith("test.")}
    assert listed == actual, (
        f"compsets.txt is out of date.\n"
        f"  listed but absent: {sorted(listed - actual)}\n"
        f"  present but unlisted: {sorted(actual - listed)}"
    )


def test_readme_compset_names_exist():
    """The README example used bp5.qd.2000; the real compset is bp5.qdc.2000."""
    readme = read("README.md")
    actual = {d.name for d in compset_dirs()}
    import re
    cited = set(re.findall(r"\b(bp\d[\w.]*|liu2020[\w.]*|test\.bp[\w.]*)\b", readme))
    # Only names shaped like compsets, not prose like "bp5" or "BP5".
    cited = {c for c in cited if "." in c}
    unknown = {c for c in cited if c not in actual}
    assert not unknown, f"README names compsets that do not exist: {sorted(unknown)}"


def test_readme_scripts_exist():
    """The README said create_newcase; the script is create.newcase."""
    readme = read("README.md")
    import re
    from conftest import ROOT
    for token in re.findall(r"\b(create[._]newcase|case\.setup|case\.submit)\b", readme):
        assert (ROOT / "scripts" / token).exists(), \
            f"README references scripts/{token}, which does not exist"


def test_readme_lists_netcdf4_dependency():
    """case.setup imports netCDF4 on line 6; the dependency list omitted it."""
    assert "netCDF4" in read("README.md"), \
        "scripts/case.setup imports netCDF4, so the README must list it"


def test_case_setup_checks_the_fault_interface_helper_exit_status():
    """os.system's return code was ignored, so a crashed helper looked like success.

    case_input/test.bp5.qdc.dip90 sets par.insertFaultType = 1, which makes
    case.setup shell out to generateFaultInterface. That script imports
    matplotlib. On CI the apt-installed matplotlib went to the system Python
    while actions/setup-python provided a different interpreter, so the import
    failed, the helper died, case.setup reported success anyway, and the run
    produced no output at all -- surfacing much later as a missing
    fault.00101.nc rather than as a setup error.
    """
    src = strip_python_comments(read("scripts/case.setup"))
    assert "generateFaultInterface" in src
    assert "waitstatus_to_exitcode" in src or "check_call" in src or "check_output" in src, (
        "case.setup must check generateFaultInterface's exit status; ignoring it "
        "turns a hard failure into a case that looks set up and then runs with "
        "no fault geometry"
    )


def test_ci_installs_matplotlib_into_the_pip_environment():
    """apt's python3-matplotlib is invisible to the setup-python interpreter."""
    wf = read(".github/workflows/test.yml")
    pip_lines = [l for l in wf.splitlines() if "pip install" in l]
    assert pip_lines, "no pip install line in the workflow"
    assert any("matplotlib" in l for l in pip_lines), (
        "generateFaultInterface imports matplotlib, so it must be pip-installed "
        "into the interpreter actions/setup-python provides, not only apt-installed"
    )
