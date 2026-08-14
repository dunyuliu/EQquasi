#! /usr/bin/env python3
"""Reading SEAS benchmark output files.

Every benchmark output file this project writes -- station time series, global.dat,
the section 4.3 profiles -- has the same three-part shape:

    # comment lines            the header
    t slip_2 slip_3 ...        one or more field-name lines
    0.0000E+00 1.234E-05 ...   numeric rows

and six scripts had each grown their own loop to take it apart. They agreed on
the important subtlety and it is worth stating once here rather than six times:

    A data row is one whose FIRST token parses as a number.

Testing for "contains a digit" misclassifies the field-name line, because names
like `slip_2` and `x2` contain digits. That distinction is the whole reason this
cannot be a two-line `numpy.loadtxt` call.

`read_table` is the single implementation. The two flags exist because the six
callers genuinely differed:

    strip     checkBP8Submission stores stripped header lines (it reports them);
              resampleBP8Profiles keeps them verbatim (it rewrites them to a new
              file, so leading whitespace must survive).
    Delimiter
              Whitespace or comma, detected per line. Gold references are
              comma-separated; run output is whitespace-separated.
    collect_bad
              checkBP8Submission must report malformed rows with line numbers
              rather than crash, since its job is to explain why a submission is
              unacceptable. Everyone else may fail loudly.
"""


def read_table(path, strip=True, collect_bad=False):
    """Split a benchmark output file into header, field-name lines and rows.

    Returns (header, names, rows), or (header, names, rows, bad) when
    collect_bad is set, where bad is a list of (line_number, truncated_text).
    With collect_bad False a malformed numeric row raises ValueError.
    """
    header, names, rows, bad = [], [], [], []
    for n, line in enumerate(open(path), 1):
        raw = line.rstrip("\n")
        text = raw.strip()
        if not text:
            continue
        if text.startswith("#"):
            header.append(text if strip else raw)
            continue
        # Run output is whitespace-separated; the frozen gold under reference/
        # is comma-separated with a one-line header. Detect per line rather than
        # requiring the caller to know which it has -- otherwise no
        # post-processing utility can be pointed at a gold directory, which is
        # how BP8's gold ended up unplottable.
        tok = text.split(",") if "," in text else text.split()
        try:
            float(tok[0])
        except ValueError:
            names.append(text if strip else raw)
            continue
        if collect_bad:
            try:
                rows.append([float(v) for v in tok])
            except ValueError:
                bad.append((n, raw[:70]))
        else:
            rows.append([float(v) for v in tok])
    if collect_bad:
        return header, names, rows, bad
    return header, names, rows


def read_rows(path):
    """Just the numeric rows, as a list of lists. The common case."""
    return read_table(path)[2]


def read_array(path):
    """Just the numeric rows, as a numpy array."""
    import numpy as np
    return np.array(read_table(path)[2])


def read_profile(path):
    """A section 4.3 profile as (x, t, field).

    Row 0 is `0 0 x2...` (or x3); every later row is
    `t max_slip_rate values...`.
    """
    d = read_array(path)
    return d[0, 2:], d[1:, 0], d[1:, 2:]
