"""The user docs site (docs/user/*.md) is written for users.

These pages are read by people running EQquasi, not by people developing it.
Developer bookkeeping -- PR numbers, rule numbers, status-board references,
session logs, agent names, commit SHAs, source file:line citations -- means
nothing to that reader and goes stale. For each page:

  1. no internal markers in prose (outside fenced code, with inline code
     stripped, since commands and paths are user content); a source
     file:line reference or a backticked commit SHA is refused anywhere;
  2. exactly one H1, it is the first heading, and nothing deeper than H3;
  3. every bullet line is at most MAX_BULLET_CHARS long.

Each check first runs on synthetic pages, one per violation, that must be
flagged, and on clean and legitimate-content pages that must pass.
"""

import re

from conftest import ROOT

DOCS = ROOT / "docs" / "user"
MAX_BULLET_CHARS = 200
AGENTS = ("mira", "iris", "lars", "kai", "haruto", "nadia", "sophia", "zofia",
          "victor", "wei-lin", "wei lin", "dunyu-liu", "anya", "marta", "priya")

# Prose only, inline code stripped. Agent handles are matched in their
# lowercase form, so ordinary capitalised words (the IRIS data centre) pass.
PROSE_MARKERS = [
    ("PR number", re.compile(r"\bPR\s*#\d+|\(#\d+\)")),
    ("rule number", re.compile(r"\brules?\s+\d+[a-z]?\b", re.I)),
    ("board reference", re.compile(
        r"pathway_forward|project_rules|board row|status board|\bitem\s+\d+", re.I)),
    ("session log", re.compile(r"session[- ]log|owner decision|owner's ruling", re.I)),
    ("agent name", re.compile(r"\b(" + "|".join(re.escape(a) for a in AGENTS) + r")\b")),
    ("bare commit SHA", re.compile(
        r"(?<![\w/.-])(?=[0-9a-f]*[a-f])(?=[0-9a-f]*[0-9])(?:[0-9a-f]{7,12}|[0-9a-f]{40})(?![\w/.-])")),
]
# Whole line, inline code included.
ANYWHERE_MARKERS = [
    ("source line reference", re.compile(r"\.(py|f90|sh|md|yml|txt):\d+")),
]
SHA = re.compile(r"`[0-9a-f]{7,40}`")


def problems(name, text):
    bad = []
    m = SHA.search(text)
    if m:
        bad.append("%s: backticked commit SHA %s" % (name, m.group(0)))
    fence, headings = False, []
    for n, line in enumerate(text.split("\n"), 1):
        if line.lstrip().startswith("```"):
            fence = not fence
            continue
        if fence:
            continue
        prose = re.sub(r"`[^`]*`", "", line)
        for label, rx in PROSE_MARKERS:
            hit = rx.search(prose)
            if hit:
                bad.append("%s:%d: %s %r" % (name, n, label, hit.group(0)))
        for label, rx in ANYWHERE_MARKERS:
            hit = rx.search(line)
            if hit:
                bad.append("%s:%d: %s %r" % (name, n, label, hit.group(0)))
        h = re.match(r"(#{1,6})\s", line)
        if h:
            headings.append((n, len(h.group(1))))
        if re.match(r"\s*([*-]|\d+\.)\s", line) and len(line) > MAX_BULLET_CHARS:
            bad.append("%s:%d: bullet is %d chars (max %d)"
                       % (name, n, len(line), MAX_BULLET_CHARS))
    h1 = [n for n, lvl in headings if lvl == 1]
    if len(h1) != 1:
        bad.append("%s: %d H1 headings (want exactly 1)" % (name, len(h1)))
    elif headings[0][1] != 1:
        bad.append("%s:%d: first heading is not the H1" % (name, headings[0][0]))
    bad += ["%s:%d: heading deeper than H3" % (name, n) for n, lvl in headings if lvl > 3]
    return bad


def test_checker_can_fail():
    clean = ("# Title\n\n## Install\n\n```\npython3 -m pytest testsys/\n```\n\n"
             "* Short bullet with `script/case.setup`.\n")
    legit = clean + ("\nData can be fetched from IRIS; row 1 of the table holds "
                     "coordinates.\nRun `python3 -m pytest testsys/` from the root.\n")
    assert problems("clean.md", clean) == []
    assert problems("legit.md", legit) == [], problems("legit.md", legit)
    cases = {
        "PR number": "Fixed in PR #12.",
        "rule number": "See rule 14.",
        "board reference": "Tracked in PATHWAY_FORWARD.md row 15.",
        "rules file": "As PROJECT_RULES says.",
        "session log": "Recorded in the session log.",
        "agent name": "Ported by mira.",
        "bare sha": "Fixed in a9df42d last week.",
        "backticked sha": "Reference `a9df42d`.",
        "source line": "See `faulting.f90:623`.",
        "two H1": "# Second",
        "H4": "#### Deep",
        "long bullet": "* " + "x" * MAX_BULLET_CHARS,
    }
    for label, extra in cases.items():
        assert problems("bad.md", clean + "\n" + extra + "\n"), \
            "the %r check passed a page it should fail" % label


def test_pages_are_written_for_users():
    pages = sorted(DOCS.glob("*.md"))
    assert pages, "no docs/user/*.md pages found"
    bad = []
    for page in pages:
        bad += problems("docs/user/" + page.name, page.read_text())
    assert not bad, "developer-facing content in the user docs:\n" + "\n".join(bad)
