"""Guard for PROJECT_RULES.md rule 22: the root is a whitelist, not a preference.

Two mechanical checks (fail the build):
  - the tracked top-level entries match rule 22's table exactly
  - no tracked file is over 5 MB

One report-only step (prints, never fails): `test_tidy_report` lists stale
worktrees under scratch/, merged local branches, and scratch/ run directories
PATHWAY_FORWARD.md does not mention -- tidying is judgment, not a gate.
"""

import subprocess

from conftest import ROOT

# Rule 22's table, root-relative. Directories are listed with a trailing "/"
# to distinguish them from files sharing a name (none here, but keeps intent
# explicit); git ls-tree gives bare top-level names for both, so we strip it
# before comparing.
ALLOWED_ROOT_ENTRIES = {
    "README.md",
    "PATHWAY_FORWARD.md",
    "PROJECT_RULES.md",
    "LICENSE",
    "environment.yml",
    "install.eqquasi.sh",
    "install.mumps.sh",
    "ubuntu.env.setup.sh",
    "make.scripts.executable.sh",
    "pytest.ini",
    "archive",
    "compset",
    "docs",
    "reference",
    "script",
    "src",
    "testsys",
    ".github",
    ".gitignore",
    # CLAUDE.md is in rule 22's table but not yet created (noted there as a
    # board item, PATHWAY_FORWARD.md row 18) -- NOT included here, so its
    # eventual addition does not need this test touched.
}

MAX_TRACKED_FILE_BYTES = 5 * 1024 * 1024  # 5 MB


def _git(args, report_only=False):
    r = subprocess.run(["git"] + args, cwd=str(ROOT), capture_output=True, text=True)
    if report_only:
        # The tidy report must never fail: CI's checkout has no local master
        # branch, and a missing ref is not a layout violation.
        return r.stdout if r.returncode == 0 else ""
    assert r.returncode == 0, f"git {args} failed: {r.stderr}"
    return r.stdout


def test_tracked_root_matches_rule_22():
    names = {
        line.split("\t", 1)[1] if "\t" in line else line
        for line in _git(["ls-tree", "--name-only", "HEAD"]).splitlines()
        if line.strip()
    }
    extra = names - ALLOWED_ROOT_ENTRIES
    missing = ALLOWED_ROOT_ENTRIES - names
    assert not extra, (
        f"tracked at repo root but not in rule 22's whitelist: {sorted(extra)} "
        "-- add a row to rule 22 in PROJECT_RULES.md, or fold the content "
        "into an existing entry"
    )
    assert not missing, (
        f"rule 22 whitelists these but they are no longer tracked at root: "
        f"{sorted(missing)} -- PROJECT_RULES.md rule 22 is now stale"
    )


def test_no_tracked_file_over_5mb():
    out = _git(["ls-tree", "-r", "-l", "HEAD"])
    offenders = []
    for line in out.splitlines():
        if not line.strip():
            continue
        meta, path = line.split("\t", 1)
        parts = meta.split()
        size = parts[3] if len(parts) > 3 else "-"
        if size == "-":
            continue
        if int(size) > MAX_TRACKED_FILE_BYTES:
            offenders.append(f"{path} ({int(size)} bytes)")
    assert not offenders, (
        f"tracked file(s) over {MAX_TRACKED_FILE_BYTES} bytes: {offenders} "
        "-- rule 22 caps committed file size"
    )


def test_tidy_report():
    """Report-only: never fails. Prints what a human tidy pass would remove."""
    lines = ["", "--- rule 22 tidy report (informational only) ---"]

    wt = _git(["worktree", "list", "--porcelain"], report_only=True)
    worktrees = [l.split(" ", 1)[1] for l in wt.splitlines() if l.startswith("worktree ")]
    stale = [w for w in worktrees if "/scratch/" in w and w.rstrip("/") != str(ROOT)]
    lines.append(f"worktrees under scratch/: {stale or 'none'}")

    branches = _git(["branch", "--merged", "master"], report_only=True).splitlines()
    merged = [b.strip().lstrip("* ").strip() for b in branches if b.strip() and "master" not in b]
    lines.append(f"local branches merged into master: {merged or 'none'}")

    scratch_dir = ROOT / "scratch"
    run_dirs = sorted(p.name for p in scratch_dir.iterdir()) if scratch_dir.is_dir() else []
    pathway = (ROOT / "PATHWAY_FORWARD.md").read_text() if (ROOT / "PATHWAY_FORWARD.md").is_file() else ""
    uncited = [d for d in run_dirs if d not in pathway]
    lines.append(f"scratch/ entries not mentioned in PATHWAY_FORWARD.md: {uncited or 'none'}")

    print("\n".join(lines))
