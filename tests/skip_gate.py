"""Fail a pytest log that holds a skip not named in tests/skip_allowlist.txt.

Rule 3a: in the gate a skip is a failure, except a skip whose reason names
absent reference data, which is reported and does not fail the run. CI used to
read only the passed/failed/error counts, so a skipped test could never turn it
red.

Usage:  python3 tests/skip_gate.py <pytest-log> [<pytest-log> ...]
The logs come from `pytest -rs`, whose lines read
`SKIPPED [N] path:line: reason`.  Allow-list lines are
`<path><TAB><reason substring>`; `#` starts a comment.  A skip is allowed only
when both its file and its reason match an entry.  Line numbers are not part of
the key, so editing a test does not silently disarm the gate.  Entries under
the `[data-absence]` header are reported (reason and count printed) and exit 0;
entries before it are tolerated silently-by-reason (a binary built without
PETSc).  Anything else fails.
"""
import re
import sys
from pathlib import Path

ALLOWLIST = Path(__file__).with_name("skip_allowlist.txt")
SKIP_LINE = re.compile(r"^SKIPPED \[(\d+)\] (\S+?):(\d+): (.*)$")
# A real pytest summary line: `=== 3 passed, 1 skipped in 1.00s ===` or the -q
# form `3 passed, 1 skipped in 1.00s (0:00:01)`; anchored, so a test name or a
# message that merely contains "1 skipped" is not one.
SUMMARY = re.compile(
    r"^(?:=+ )?(?:no tests ran|\d+ [a-z]+(?:, \d+ [a-z]+)*)"
    r" in \d+(?:\.\d+)?s(?: \(\d+:\d\d:\d\d\))?(?: =+)?$")
SECTION = "[data-absence]"


def load_allowlist(path=ALLOWLIST):
    """Return (tolerated, data_absence): two lists of (path, reason substring)."""
    tolerated, absence, cur = [], [], None
    cur = tolerated
    for n, raw in enumerate(Path(path).read_text().splitlines(), 1):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        if line == SECTION:
            cur = absence
            continue
        parts = line.split("\t")
        if len(parts) != 2 or not all(p.strip() for p in parts):
            raise ValueError(f"{path}:{n}: expected '<path>TAB<reason>', got {raw!r}")
        cur.append((parts[0].strip(), parts[1].strip()))
    return tolerated, absence


def classify_skips(log_text, allowlist):
    """Return (forbidden, reported): SKIPPED lines no entry covers, and those
    covered by a [data-absence] entry.

    Raises ValueError when the log has no pytest summary line (the run did not
    happen, or the log is truncated) or when the last summary line's skip count
    differs from the sum of the printed SKIPPED [n] reasons (the log was made
    without -rs, or lost lines).  Trailing output after the summary is fine.
    """
    tolerated, absence = allowlist
    lines = log_text.splitlines()
    summaries = [l.strip() for l in lines if SUMMARY.match(l.strip())]
    if not summaries:
        raise ValueError("no pytest summary line in the log; the run did not happen")
    bad, reported, seen = [], [], 0
    for line in lines:
        m = SKIP_LINE.match(line)
        if not m:
            continue
        seen += int(m.group(1))
        path, reason = m.group(2), m.group(4)
        if any(path == p and sub in reason for p, sub in tolerated):
            continue
        if any(path == p and sub in reason for p, sub in absence):
            reported.append(line)
        else:
            bad.append(line)
    m = re.search(r"(\d+) skipped", summaries[-1])
    n = int(m.group(1)) if m else 0
    if n != seen:
        raise ValueError(f"summary reports {n} skipped but {seen} reasons "
                         "were printed; run pytest with -rs")
    return bad, reported


def forbidden_skips(log_text, allowlist):
    return classify_skips(log_text, allowlist)[0]


def main(argv):
    if not argv:
        raise SystemExit(__doc__)
    allow = load_allowlist()
    failed = False
    for log in argv:
        bad, reported = classify_skips(Path(log).read_text(), allow)
        for line in reported:
            print(f"{log}: data-absence skip (reported, rule 3a): {line}")
        if reported:
            print(f"{log}: {sum(int(SKIP_LINE.match(l).group(1)) for l in reported)} "
                  "data-absence skip(s) reported")
        for line in bad:
            print(f"{log}: skip not on tests/skip_allowlist.txt: {line}")
        failed = failed or bool(bad)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
