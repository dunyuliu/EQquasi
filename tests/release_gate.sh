#!/bin/bash
# Full-tier release gate (rule 9a): a bare clone of <sha> under env -i with an
# empty HOME, install build, fast tier, then the e2e tier behind a shared flock.
# usage: MACHINE=<machine> CONDA_BASE=<conda root> GATE_LOCK=<lock file> \
#          tests/release_gate.sh <sha> <outdir>
# Writes <outdir>/status.txt with each stage's exit code; exits non-zero and
# names every failed stage if any stage failed. A gate run from any other
# script is not a gate.
set -u

die() { echo "release_gate: $*" >&2; exit 2; }

[ $# -eq 2 ] || die "usage: release_gate.sh <sha> <outdir> (got $# args)"
SHA=$1; OUT=$2
[ -n "$SHA" ] || die "<sha> is empty"
[ -n "$OUT" ] || die "<outdir> is empty"
[ -n "${MACHINE:-}" ] || die "MACHINE is unset; set it explicitly (e.g. MACHINE=conda-linux)"
[ -n "${CONDA_BASE:-}" ] || die "CONDA_BASE is unset; set it to the conda root holding eqquasi-petsc"
[ -x "$CONDA_BASE/bin/conda" ] || die "CONDA_BASE=$CONDA_BASE has no executable bin/conda"
[ -n "${GATE_LOCK:-}" ] || die "GATE_LOCK is unset; set it to the shared heavy-run lock file"
export MACHINE

HERE=$(cd "$(dirname "$0")" && pwd) || die "cannot resolve script directory"
URL=$(git -C "$HERE" remote get-url origin) || die "invoking checkout has no 'origin' remote"
[ -n "$URL" ] || die "origin URL is empty"

[ ! -e "$OUT/clone" ] || die "$OUT/clone already exists; use a fresh <outdir>"
mkdir -p "$OUT/home" || die "cannot create $OUT"
OUT=$(cd "$OUT" && pwd) || die "cannot resolve $OUT"
W=$OUT/clone
STATUS=$OUT/status.txt
# The e2e stage runs `run` inside `flock ... bash -c`: that child sees only exported variables.
export OUT W CONDA_BASE
echo "start $(date +%T) load $(cut -d' ' -f1-3 /proc/loadavg) machine=$MACHINE url=$URL" > "$STATUS"

FAILED=""
record() {  # record <stage> <rc>
  echo "$1 rc=$2 $(date +%T)" >> "$STATUS"
  if [ "$2" -ne 0 ]; then FAILED="$FAILED $1"; fi
}

git clone -q "$URL" "$W" 2> "$OUT/clone.err"; rc=$?
if [ $rc -eq 0 ]; then git -C "$W" checkout -q "$SHA" 2>> "$OUT/clone.err"; rc=$?; fi
record clone $rc
if [ $rc -ne 0 ]; then
  echo "done FAILED: clone" >> "$STATUS"
  echo "release_gate: FAILED stage(s): clone (see $OUT/clone.err)" >&2
  exit 1
fi
echo "clone head $(git -C "$W" rev-parse HEAD)" >> "$STATUS"

# run <command string>: env -i, empty HOME, conda env eqquasi-petsc, cwd = clone.
run() {
  env -i HOME="$OUT/home" PATH=/usr/bin:/bin OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
      EQQUASIROOT="$W" MACHINE="$MACHINE" CONDA_BASE="$CONDA_BASE" CLONE="$W" \
      bash -c 'eval "$($CONDA_BASE/bin/conda shell.bash hook)" && conda activate eqquasi-petsc && cd "$CLONE" && '"$1"
}

run 'bash install.eqquasi.sh -m "$MACHINE"' > "$OUT/build.log" 2>&1; record build $?
run 'python3 -m pytest -q -rs -m "not e2e" tests' > "$OUT/fast.log" 2>&1; record fast $?
flock "$GATE_LOCK" bash -c "$(declare -f run); run 'python3 -m pytest -q -rs -m e2e tests'" > "$OUT/e2e_full.log" 2>&1; record e2e_full $?

if [ -n "$FAILED" ]; then
  echo "done FAILED:$FAILED" >> "$STATUS"
  echo "release_gate: FAILED stage(s):$FAILED (see $STATUS)" >&2
  exit 1
fi
echo "done PASS" >> "$STATUS"
