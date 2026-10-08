#! /usr/bin/env python3
"""Generate docs/user/parameters.md's parameter reference from
scripts/defaultParameters.py's `parameters` class -- the class case.setup and
every compset's user_defined_params.py read their defaults from -- so the
page cannot list a knob that does not exist or show a stale default.

What is extracted:
  every assignment directly in the `parameters` class body
  (`name = value` or `a, b = 1, 2`), in source order, with its literal
  default and a note built from the comment on its own line plus the
  contiguous comment block directly above and below it. A
  separator / `##### Title #####` / separator banner becomes a heading.

What is left out:
  - derived working state, not something a user sets (SKIP_NAMES);
  - the per-node loop that fills on_fault_vars (it computes arrays from the
    parameters already listed);
  - any comment block containing a developer-only reference (a rule number,
    a source-file name, an agent name ...). The whole block is dropped, never
    one line of it, and each drop is printed so it is visible.

Usage:
  python3 docs/user/gen_params.py            # regenerate parameters.md
  python3 docs/user/gen_params.py --check    # exit 1 if it would change
"""
import ast
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
SOURCE = os.path.join(ROOT, 'scripts', 'defaultParameters.py')
TARGET = os.path.join(HERE, 'parameters.md')

BEGIN = ('<!-- BEGIN PARAMETER REFERENCE (generated from '
         'scripts/defaultParameters.py by docs/user/gen_params.py; do not '
         'edit by hand) -->')
END = '<!-- END PARAMETER REFERENCE -->'

# Working state computed from the parameters, not set by a user.
SKIP_NAMES = {'on_fault_vars', 'fx', 'fz', 'nfx', 'nfz',
              'n_on_fault', 'n_off_fault'}

# The same class of internal reference tests/contract/test_user_docs_style.py
# refuses on a user page.
INTERNAL_MARKERS = [
    re.compile(r'\bPR\s*#\d+|\(#\d+\)'),
    re.compile(r'\brules?\s+\d+[a-z]?\b', re.I),
    re.compile(r'pathway_forward|project_rules|board row|\bitem\s+\d+', re.I),
    re.compile(r'\b(mira|iris|lars|kai|haruto|nadia|sophia|zofia|victor|'
               r'wei-lin|wei lin|dunyu-liu|anya|marta|priya)\b'),
    re.compile(r'(?<![\w/.-])(?=[0-9a-f]*[a-f])(?=[0-9a-f]*[0-9])'
               r'(?:[0-9a-f]{7,12}|[0-9a-f]{40})(?![\w/.-])'),
    re.compile(r'\.(py|f90|sh|md|yml|txt):\d+'),
    # A Fortran source file is never something a user of this page edits.
    re.compile(r'\b[\w/]+\.f90\b'),
    re.compile(r"owner's ruling|owner decision", re.I),
]

SEP_RE = re.compile(r'^#{5,}\s*$')
BANNER_RE = re.compile(r'^#{2,}\s+(.+?)\s+#{2,}$')


def _is_internal(line):
    return any(rx.search(line) for rx in INTERNAL_MARKERS)


def _find_banners(lines):
    """({title line index: title}, set of every line a banner occupies)."""
    banners, consumed = {}, set()
    for i in range(1, len(lines) - 1):
        if SEP_RE.match(lines[i - 1].strip()) and SEP_RE.match(lines[i + 1].strip()):
            m = BANNER_RE.match(lines[i].strip())
            if m:
                banners[i] = m.group(1).strip()
                consumed |= {i - 1, i, i + 1}
    return banners, consumed


def _comment_block(lines, start, step, consumed, taken):
    """Contiguous comment-only lines from `start`, moving by `step`, stopping
    at code or at a line already claimed. Returned in reading order."""
    out, used, i = [], set(), start
    while 0 <= i < len(lines):
        raw = lines[i]
        if not raw.strip().startswith('#') or i in consumed or i in taken:
            break
        out.append(raw)
        used.add(i)
        i += step
    if step < 0:
        out.reverse()
    return out, used


def _clean_note(raw_lines, dropped, owner):
    """Join a comment block into one note, or drop it whole if any line is
    internal (splicing around one removed line leaves a broken sentence)."""
    texts = [raw.strip().lstrip('#').strip() for raw in raw_lines]
    texts = [t for t in texts if t]
    if not texts:
        return ''
    hit = next((t for t in texts if _is_internal(t)), None)
    if hit is not None:
        dropped.append('%s: dropped a %d-line comment block (contains %r)'
                       % (owner, len(texts), hit))
        return ''
    return ' '.join(texts)


def _assignments(cls):
    """[(node, names, values)] for every plain assignment in the class body."""
    out = []
    for node in cls.body:
        if not isinstance(node, ast.Assign):
            continue                      # loops, defs: not knobs
        target = node.targets[0]
        if isinstance(target, ast.Tuple):
            names = [t.id for t in target.elts if isinstance(t, ast.Name)]
        elif isinstance(target, ast.Name):
            names = [target.id]
        else:
            continue
        if isinstance(node.value, ast.Tuple) and len(node.value.elts) == len(names):
            values = node.value.elts
        else:
            values = [node.value] * len(names)
        out.append((node, names, values))
    return out


def extract(source_text=None):
    """Returns (rows, banners, dropped).
    rows    -- [(lineno, [names], [defaults], note)], one per assignment.
    banners -- [(lineno, title)].

    A comment block belongs to the assignment directly below it. A block
    directly below an assignment is attached to it only when no assignment
    claimed it from below (two passes, so the order is unambiguous).
    """
    text = source_text if source_text is not None else open(SOURCE).read()
    lines = text.splitlines()
    banners, consumed = _find_banners(lines)
    cls = next(n for n in ast.parse(text).body
               if isinstance(n, ast.ClassDef) and n.name == 'parameters')
    assigns = _assignments(cls)
    taken, lead, trail = set(), {}, {}
    for node, _, _ in assigns:
        lead[node], used = _comment_block(lines, node.lineno - 2, -1, consumed, taken)
        taken |= used
    for node, _, _ in assigns:
        trail[node], used = _comment_block(lines, node.end_lineno, +1, consumed, taken)
        taken |= used

    rows, dropped = [], []
    for node, names, values in assigns:
        keep = [(n, v) for n, v in zip(names, values) if n not in SKIP_NAMES]
        if not keep:
            continue
        defaults = []
        for _, value in keep:
            try:
                defaults.append(repr(ast.literal_eval(value)))
            except Exception:
                defaults.append(ast.get_source_segment(text, value))
        owner = '%s (scripts/defaultParameters.py line %d)' % (
            ', '.join(n for n, _ in keep), node.lineno)
        last = lines[node.end_lineno - 1]
        inline = last.split('#', 1)[1].strip() if '#' in last else ''
        parts = [_clean_note(lead[node], dropped, owner)]
        if inline:
            if _is_internal(inline):
                dropped.append('%s: dropped inline comment %r' % (owner, inline))
            else:
                parts.append(inline)
        parts.append(_clean_note(trail[node], dropped, owner))
        note = ' '.join(p for p in parts if p)
        if _is_internal(note):
            raise RuntimeError('gen_params: %s note is still internal: %r'
                               % (owner, note))
        rows.append((node.lineno, [n for n, _ in keep], defaults, note))
    return rows, sorted(banners.items()), dropped


def render(rows, banners):
    events = [(ln, 'banner', title) for ln, title in banners]
    events += [(ln, 'row', row) for ln, *row in rows]
    events.sort(key=lambda e: e[0])
    out = [BEGIN, '',
           'Every entry below is an attribute of the `parameters` class in '
           '`scripts/defaultParameters.py`. A case overrides any of them in its '
           'own `user_defined_params.py` as `par.<name> = <value>`. The '
           'defaults are the class\'s own, taken from BP5 at 2000 m; each '
           'compset sets its own values on top.', '']
    for _, kind, data in events:
        if kind == 'banner':
            out += ['## %s' % data, '']
            continue
        names, defaults, note = data
        label = ', '.join('**`%s`**' % n for n in names)
        if all('\n' not in d and len(d) <= 100 for d in defaults):
            word = 'default' if len(defaults) == 1 else 'defaults'
            out.append('* %s -- %s %s' % (label, word,
                                          ', '.join('`%s`' % d for d in defaults)))
        else:
            out += ['* %s' % label, '', '  Default:', '', '  ```python']
            for d in defaults:
                out += ['  ' + l for l in d.splitlines()]
            out.append('  ```')
        if note:
            out += ['', '  ' + note]
        out.append('')
    out.append(END)
    return '\n'.join(out)


def check_or_update(update):
    if not os.path.exists(TARGET):
        return ['docs/user/parameters.md does not exist']
    s = open(TARGET).read()
    if BEGIN not in s or END not in s:
        return ['docs/user/parameters.md has no generated parameter section']
    head, rest = s.split(BEGIN, 1)
    body, tail = rest.split(END, 1)
    rows, banners, dropped = extract()
    for d in dropped:
        print('  filtered: %s' % d)
    wanted = render(rows, banners)
    if BEGIN + body + END == wanted:
        return []
    if update:
        open(TARGET, 'w').write(head + wanted + tail)
        print('  docs/user/parameters.md regenerated (%d entries, %d sections)'
              % (len(rows), len(banners)))
        return []
    return ["docs/user/parameters.md no longer matches "
            "scripts/defaultParameters.py -- run 'python3 docs/user/gen_params.py'"]


def main():
    check = '--check' in sys.argv
    problems = check_or_update(update=not check)
    if problems:
        print('FAIL gen_params:')
        for p in problems:
            print(' -', p)
        return 1
    print('SUCCESS gen_params' + (' (check only, no write)' if check else ''))
    return 0


if __name__ == '__main__':
    sys.exit(main())
