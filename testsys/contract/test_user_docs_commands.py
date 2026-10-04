"""Every command shown on the user docs site (docs/user/*.md) must resolve.

A page that tells a user to run a script that was renamed, a compset that was
removed, an option a tool no longer takes, or a launcher override run.sh no
longer reads, fails silently for the reader and for nobody else. This reads
every shell line in the pages' fenced code blocks and checks, against the
repository itself:

  - the program exists: a file in script/ or the repository root, a file
    create.newcase puts in a case, or a known system command;
  - the compset named to create.newcase is a directory under compset/;
  - the machine named to install.eqquasi.sh is one the installer handles;
  - every --option given to a script/ tool is one that tool defines;
  - every pytest -m marker is declared in pytest.ini;
  - every VAR=... override in front of `bash run.sh` is one the run.sh that
    case.setup writes actually reads;
  - every repository path quoted in inline code (script/..., compset/...,
    reference/..., testsys/..., src/...) exists.

Fences tagged with a language other than sh/bash/console (```python,
```text) are examples, not commands, and are skipped. Each check is also run
on synthetic pages that must fail, so none of them can pass vacuously.
"""

import re
import shlex

from conftest import ROOT, read

DOCS = ROOT / "docs" / "user"

# Commands a page may use that are not part of this repository.
SYSTEM = {"cd", "export", "git", "conda", "rm", "cp", "mkdir", "uptime",
          "source", "bash", "python3", "sudo", "pip"}

# Files create.newcase / case.setup put in a case root, run from there.
CASE_FILES = {"case.setup", "case.submit", "run.sh"}

ENV_ASSIGN = re.compile(r"^[A-Z_][A-Z0-9_]*=")
REPO_PATH = re.compile(r"`((?:script|compset|reference|testsys|src)/[^`\s<>*]+)`")


def doc_pages():
    return sorted(DOCS.glob("*.md"))


def shell_lines(text):
    """(lineno, line) for every command line in an untagged or shell fence."""
    out, fence, lang = [], False, ""
    for n, line in enumerate(text.splitlines(), 1):
        s = line.strip()
        if s.startswith("```"):
            if not fence:
                fence, lang = True, s[3:].strip().lower()
            else:
                fence = False
            continue
        if fence and lang in ("", "sh", "bash", "console") and s and not s.startswith("#"):
            out.append((n, s))
    return out


def commands(line):
    """Split one shell line into simple commands, comments dropped."""
    line = re.sub(r"\s+#.*$", "", line)
    return [c.strip() for c in re.split(r"&&|\|\||;|\|", line) if c.strip()]


def installer_machines():
    """Machines with a branch in the installer or the makefile."""
    branches = set(re.findall(r'MACHINE"?\s*==\s*"([\w-]+)"', read("install.eqquasi.sh")))
    listed = re.search(r"Supported:\s*([\w ,-]+)", read("src/makefile"))
    return branches | set(listed.group(1).replace(",", " ").split() if listed else [])


def pytest_markers():
    ini = read("pytest.ini")
    return set(re.findall(r"^\s{4}(\w+):", ini, re.M))


def run_sh_env_vars():
    """Variables the run.sh written by case.setup takes from the environment."""
    return set(re.findall(r"\$\{(\w+):?-", read("script/case.setup")))


def tool_defines(path, option):
    return re.search(r"""["']%s["']""" % re.escape(option), path.read_text()) is not None


def check_command(cmd):
    """Problems with one simple command, as strings."""
    try:
        toks = shlex.split(cmd)
    except ValueError as exc:
        return ["cannot parse %r (%s)" % (cmd, exc)]
    env = []
    while toks and ENV_ASSIGN.match(toks[0]):
        env.append(toks.pop(0).split("=", 1)[0])
    if not toks:
        return []
    if toks[0] == "sudo":
        toks = toks[1:]
    prog, args = toks[0], toks[1:]
    bad = []

    if prog == "python3" and args[:1] == ["-m"]:
        if args[1:2] != ["pytest"]:
            bad.append("python3 -m %s is not a module this project runs" % args[1:2])
        for i, a in enumerate(args):
            if a == "-m" and i + 1 < len(args) and i > 0:
                for word in re.findall(r"\w+", args[i + 1]):
                    if word not in ("not", "and", "or") and word not in pytest_markers():
                        bad.append("pytest marker %r is not declared in pytest.ini" % word)
        return bad
    if prog in ("bash", "source", "python3") and args:
        prog, args = args[0], args[1:]
    prog = prog[2:] if prog.startswith("./") else prog

    if prog in SYSTEM:
        return bad
    if env and prog != "run.sh":
        bad.append("environment override %s in front of %s, not run.sh" % (env, prog))
    for var in env:
        if var not in run_sh_env_vars():
            bad.append("run.sh does not read %s from the environment" % var)

    script = ROOT / "script" / prog
    if prog in CASE_FILES:
        if prog != "run.sh" and not (ROOT / "script" / prog).is_file():
            bad.append("%s is not in script/" % prog)
        if prog == "case.setup":
            for a in args:
                if a.startswith("--") and not tool_defines(script, a):
                    bad.append("case.setup does not take %s" % a)
        return bad
    if prog.startswith("input/"):
        name = prog.split("/", 1)[1]
        if not list((ROOT / "compset").glob("*/" + name)):
            bad.append("no compset ships %s" % name)
        return bad
    if script.is_file():
        for a in args:
            if a.startswith("--") and not tool_defines(script, a.split("=", 1)[0]):
                bad.append("%s does not define option %s" % (prog, a))
        if prog == "create.newcase":
            if len(args) != 2:
                bad.append("create.newcase takes <case dir> <compset>")
            elif not (ROOT / "compset" / args[1]).is_dir():
                bad.append("compset %s does not exist" % args[1])
        return bad
    if (ROOT / prog).is_file():
        if prog == "install.eqquasi.sh" and "-m" in args:
            m = args[args.index("-m") + 1]
            if m not in installer_machines():
                bad.append("install.eqquasi.sh does not support -m %s" % m)
        return bad
    return bad + ["%s is neither in script/, the repository root, a case, "
                  "nor a known system command" % prog]


def problems(name, text):
    bad = []
    for n, line in shell_lines(text):
        for cmd in commands(line):
            bad += ["%s:%d: %s" % (name, n, p) for p in check_command(cmd)]
    for n, line in enumerate(text.splitlines(), 1):
        for path in REPO_PATH.findall(line):
            if not (ROOT / path).exists():
                bad.append("%s:%d: %s does not exist" % (name, n, path))
    return bad


def test_checker_can_fail():
    clean = ("# T\n\n```\ncreate.newcase scratch/x test.bp5.qdc.2000\n"
             "./case.setup --force\nbash run.sh\nplotStations.py --time y\n"
             "python3 -m pytest testsys/ -m e2e_fast\n"
             "OMP_NUM_THREADS=1 bash run.sh\n```\n\n```python\nnot_a_command x\n```\n")
    assert problems("clean.md", clean) == []
    broken = {
        "unknown program": "```\nplotNothing.py\n```\n",
        "missing compset": "```\ncreate.newcase x no.such.compset\n```\n",
        "unknown option": "```\nplotStations.py --no-such-option\n```\n",
        "unknown marker": "```\npython3 -m pytest testsys/ -m e2e_nope\n```\n",
        "unknown machine": "```\nbash install.eqquasi.sh -m vax\n```\n",
        "unread override": "```\nNOT_READ=1 bash run.sh\n```\n",
        "missing path": "See `script/noSuchTool.py`.\n",
        "case.setup option": "```\n./case.setup --nope\n```\n",
    }
    for label, text in broken.items():
        assert problems("broken.md", "# T\n\n" + text), \
            "the %r check passed a page it should fail" % label


def test_pages_exist():
    assert doc_pages(), "no docs/user/*.md pages found"


def test_every_documented_command_resolves():
    bad = []
    for page in doc_pages():
        bad += problems("docs/user/" + page.name, page.read_text())
    assert not bad, "commands in the user docs that do not work:\n" + "\n".join(bad)
