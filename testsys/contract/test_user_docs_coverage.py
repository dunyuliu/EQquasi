"""The user docs site must keep up with what the code ships.

The style and command checks cannot tell that a page has fallen behind: a
new compset, output file, stop code, install target or launcher override can
land without the docs ever mentioning it. Each check here derives its list
from the code itself and requires every entry to be named on the page whose
job it is to cover it:

  - parameters.md's generated reference is current with
    script/defaultParameters.py (docs/user/gen_params.py --check);
  - every compset directory is named in benchmarks.md;
  - every `stop <code>` in src/*.f90 has a row in troubleshooting.md;
  - every output file the solver opens in its output directory is named in
    outputs.md;
  - every machine install.eqquasi.sh or the makefile builds for is named in
    getting-started.md;
  - every variable run.sh takes from the environment is named in
    running-a-case.md.
"""

import importlib.util
import re

from conftest import ROOT, compset_dirs, read, strip_fortran_comments

DOCS = ROOT / "docs" / "user"


def page(name):
    return (DOCS / name).read_text()


def missing(items, text):
    return sorted(i for i in items if i not in text)


def test_generated_parameter_page_is_current():
    spec = importlib.util.spec_from_file_location("gen_params", DOCS / "gen_params.py")
    gen = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(gen)
    problems = gen.check_or_update(update=False)
    assert not problems, "\n".join(problems)


def test_every_compset_is_in_benchmarks():
    names = [d.name for d in compset_dirs()]
    assert names
    gap = missing(names, page("benchmarks.md"))
    assert not gap, "compsets missing from docs/user/benchmarks.md: %s" % gap


def stop_codes():
    codes = set()
    for f in (ROOT / "src").glob("*.f90"):
        codes |= set(re.findall(r"\bstop\s+(\d+)\b",
                                strip_fortran_comments(f.read_text()), re.I))
    return codes


def test_every_stop_code_is_in_troubleshooting():
    codes = stop_codes()
    assert codes, "found no stop codes in src/ -- the scan is broken"
    text = page("troubleshooting.md")
    gap = sorted((c for c in codes if not re.search(r"^\|\s*%s\s*\|" % c, text, re.M)), key=int)
    assert not gap, "stop codes with no row in docs/user/troubleshooting.md: %s" % gap


def output_files():
    """Literal names (or name prefixes) the solver opens under its output
    directory, plus the BP8 profile file names."""
    names = set()
    for f in (ROOT / "src").glob("*.f90"):
        text = strip_fortran_comments(f.read_text())
        names |= set(re.findall(r"outDir\)\s*//\s*['\"]([^'\"]+)['\"]", text))
        names |= set(re.findall(r"call\s+write_one_profile\(\s*'([^']+)'", text))
    return names


def test_every_output_file_is_in_outputs():
    names = output_files()
    assert {"global.dat", "runInfo.json"} <= names, "output-file scan is broken: %s" % names
    gap = missing(names, page("outputs.md"))
    assert not gap, "solver output files missing from docs/user/outputs.md: %s" % gap


def build_machines():
    branches = set(re.findall(r'MACHINE"?\s*==\s*"([\w-]+)"', read("install.eqquasi.sh")))
    listed = re.search(r"Supported:\s*([\w ,-]+)", read("src/makefile"))
    return branches | set(listed.group(1).replace(",", " ").split() if listed else [])


def test_every_install_target_is_in_getting_started():
    machines = build_machines()
    assert {"ubuntu", "ls6"} <= machines, "machine scan is broken: %s" % machines
    text = page("getting-started.md")
    gap = sorted(m for m in machines if "-m %s" % m not in text)
    assert not gap, "install targets missing from docs/user/getting-started.md: %s" % gap


def test_every_run_sh_override_is_in_running_a_case():
    env = set(re.findall(r"\$\{(\w+):?-", read("script/case.setup")))
    assert "MPIRUN" in env, "run.sh override scan is broken: %s" % env
    gap = missing(env, page("running-a-case.md"))
    assert not gap, "run.sh overrides missing from docs/user/running-a-case.md: %s" % gap
