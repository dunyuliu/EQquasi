"""Fail a pytest log that holds a skip not named in tests/skip_allowlist.txt.

Rule 3a: in the gate a skip is a failure. CI used to read only the
passed/failed/error counts, so a skipped test could never turn it red.

Usage:  python3 tests/skip_gate.py <pytest-log> [<pytest-log> ...]
The logs come from `pytest -rs`, whose lines read
`SKIPPED [N] path:line: reason`.  Allow-list lines are
`<path><TAB><reason substring>`; `#` starts a comment.  A skip is allowed only
when both its file and its reason match an entry.  Line numbers are not part of
the key, so editing a test does not silently disarm the gate.
"""
import re
import sys
from pathlib import Path

ALLOWLIST = Path(__file__).with_name("skip_allowlist.txt")
SKIP_LINE = re.compile(r"^SKIPPED \[(\d+)\] (\S+?):(\d+): (.*)$")
SUMMARY = re.compile(r"\d+ (passed|failed|error)")


def load_allowlist(path=ALLOWLIST):
    entries = []
    for n, raw in enumerate(Path(path).read_text().splitlines(), 1):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) != 2 or not all(p.strip() for p in parts):
            raise ValueError(f"{path}:{n}: expected '<path>TAB<reason>', got {raw!r}")
        entries.append((parts[0].strip(), parts[1].strip()))
    return entries


def forbidden_skips(log_text, allowlist):
    """Return the SKIPPED lines of log_text that no allow-list entry covers.

    Raises ValueError when the log has no pytest summary line (the run did not
    happen) or when the summary reports skips but no SKIPPED line was printed
    (the log was produced without -rs, so the reasons cannot be read).
    """
    if not SUMMARY.search(log_text):
        raise ValueError("no pytest summary line in the log; the run did not happen")
    bad, seen = [], 0
    for line in log_text.splitlines():
        m = SKIP_LINE.match(line)
        if not m:
            continue
        seen += int(m.group(1))
        path, reason = m.group(2), m.group(4)
        if not any(path == p and sub in reason for p, sub in allowlist):
            bad.append(line)
    m = re.search(r"(\d+) skipped", log_text.splitlines()[-1] if log_text.strip() else "")
    if m and int(m.group(1)) != seen:
        raise ValueError(f"summary reports {m.group(1)} skipped but {seen} reasons "
                         "were printed; run pytest with -rs")
    return bad


def main(argv):
    if not argv:
        raise SystemExit(__doc__)
    allow = load_allowlist()
    failed = False
    for log in argv:
        bad = forbidden_skips(Path(log).read_text(), allow)
        for line in bad:
            print(f"{log}: skip not on tests/skip_allowlist.txt: {line}")
        failed = failed or bool(bad)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
