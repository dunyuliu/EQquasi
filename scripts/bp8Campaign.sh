#! /bin/bash
# Autonomous BP8 study. Runs sequentially, appends one summary line per case to
# work/campaign.log, and never overwrites a directory another case is using.
#
# The problem is small, so everything runs on a single rank with a single
# OpenMP thread unless the element count makes that impractical.

set -u
R=/home/utig5/dliu/seas_bp10_eqquasi
EX=$R/bin/eqquasi
export EQQUASIROOT=$R
export PATH=$R/bin:$R/scripts:$PATH
export OMP_NUM_THREADS=1
LOG=$R/work/campaign.log

say() { echo "[$(date +%H:%M:%S)] $*" | tee -a "$LOG"; }

# ranks: 1 unless the mesh is large enough that a direct factorisation needs help
ranks_for() { [ "$1" -gt 150000 ] && echo 4 || echo 1; }

summarise() {  # $1 = case dir, $2 = label
  python3 - "$1" "$2" <<'PY' | tee -a "$LOG"
import sys, numpy as np, os, re
d, lab = sys.argv[1], sys.argv[2]
def rd(p):
    r=[]
    for l in open(p):
        l=l.strip()
        if not l or l.startswith('#'): continue
        t=l.split()
        try: float(t[0])
        except ValueError: continue
        r.append([float(v) for v in t])
    return np.array(r)
try:
    g = rd(os.path.join(d,'global.dat'))
    ne = ''
    log = os.path.join(d,'r.log')
    if os.path.exists(log):
        m = re.findall(r'Total elems =\s+(\d+)', open(log, errors='ignore').read())
        ne = m[-1] if m else ''
    s = rd(os.path.join(d,'fltst_strk+000dp+000.txt'))
    print(f"RESULT {lab:22s} elems={ne:>8s} t_end={g[-1,0]/86400:6.2f}d "
          f"peakV={g[:,1].max():7.3f}@{g[g[:,1].argmax(),0]/86400:5.2f}d "
          f"peakM={g[:,2].max():.3e} slip={s[-1,1]*1000:7.2f}mm")
except Exception as e:
    print(f"RESULT {lab:22s} FAILED: {e}")
PY
}

run_case() {   # $1 dir  $2 label  $3..= sed edits already applied by caller
  local D=$1 LAB=$2
  cd "$D" || return 1
  ./case.setup > setup.log 2>&1
  if [ ! -f model.txt ]; then say "SETUP FAILED $LAB"; return 1; fi
  echo 1 > currentcycle.txt
  local NE NP
  NE=$(grep -oiE "Estimated total cells in the model is [0-9]+" setup.log | grep -oE "[0-9]+" || echo 0)
  NP=$(ranks_for "${NE:-0}")
  say "start $LAB (est cells ${NE:-?}, np=$NP)"
  /usr/bin/time -v mpirun -np "$NP" "$EX" > r.log 2> t.log
  say "done  $LAB  wall=$(grep 'Elapsed (wall' t.log | awk '{print $NF}')"
  summarise "$D" "$LAB"
}

mkcase() {     # $1 dir  $2 compset
  rm -rf "$1"; cd $R && create.newcase "$1" "$2" > /dev/null
}

say "===== BP8 campaign start ====="

# ---------------------------------------------------------------- 1. scaling
say "--- phase 1: MPI/OpenMP scaling on the 8000-element case"
mkcase $R/work/camp.scal test.bp8.qdc
cd $R/work/camp.scal
sed -i "s/^par.nstep       = .*/par.nstep       = 200/" user_defined_params.py
./case.setup > setup.log 2>&1
for NP in 1 2 4 8; do
  echo 1 > currentcycle.txt
  /usr/bin/time -f "%e" mpirun -np $NP "$EX" > r.log 2> tt.log
  SPS=$(grep -oE "[0-9]+ +steps use +[0-9.]+" r.log | awk '{printf "%.4f", $3/$1}')
  say "SCALING np=$NP OMP=1  s/step=${SPS:-NA}  wall=$(tail -1 tt.log)"
done

# ------------------------------------------------------------ 2. domain sweep
say "--- phase 2: domain sweep, half-width in x, y and z"
for H in 500 1000 2000; do
  mkcase $R/work/camp.dom$H test.bp8.qdc
  cd $R/work/camp.dom$H
  sed -i "s/^par.fxmin, par.fxmax = .*/par.fxmin, par.fxmax = -$H.0, $H.0/" user_defined_params.py
  sed -i "s/^par.fymin, par.fymax = .*/par.fymin, par.fymax = -$H.0, $H.0/" user_defined_params.py
  sed -i "s/^par.fzmin, par.fzmax = .*/par.fzmin, par.fzmax = -$H.0, $H.0/" user_defined_params.py
  sed -i "s/^par.nstep       = .*/par.nstep       = 4000/" user_defined_params.py
  sed -i "s/^par.nt_out      = .*/par.nt_out      = 4000/" user_defined_params.py
  run_case $R/work/camp.dom$H "domain_${H}m"
done

# ------------------------------------------------------- 3. xi (time step) sweep
say "--- phase 3: xi sweep at the 1000 m box"
for XI in 0.05 0.1 0.2; do
  mkcase $R/work/camp.xi$XI test.bp8.qdc
  cd $R/work/camp.xi$XI
  sed -i "s/^par.fxmin, par.fxmax = .*/par.fxmin, par.fxmax = -1000.0, 1000.0/" user_defined_params.py
  sed -i "s/^par.fymin, par.fymax = .*/par.fymin, par.fymax = -1000.0, 1000.0/" user_defined_params.py
  sed -i "s/^par.fzmin, par.fzmax = .*/par.fzmin, par.fzmax = -1000.0, 1000.0/" user_defined_params.py
  sed -i "s/^par.xi = .*/par.xi = $XI/" user_defined_params.py
  sed -i "s/^par.nstep       = .*/par.nstep       = 4000/" user_defined_params.py
  sed -i "s/^par.nt_out      = .*/par.nt_out      = 4000/" user_defined_params.py
  run_case $R/work/camp.xi$XI "xi_${XI}"
done

# ------------------------------------------------------------- 4. comparison plot
say "--- phase 4: comparison plot"
cd $R && python3 scripts/plotDomainSweep.py >> "$LOG" 2>&1

say "===== BP8 campaign complete ====="
