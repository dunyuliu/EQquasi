#! /bin/bash
# xi sweep at the 500 m box, then a 2000 m domain confirmation.
#
# The xi sweep moved from the 1000 m box to the 500 m box because the domain
# sweep showed peak slip rate identical to three decimals and slip within 1 %
# between them: domain size does not affect the answer here, so there is no
# reason to pay 64k elements when 8k gives the same physics eight times faster.

set -u
R=/home/utig5/dliu/seas_bp10_eqquasi
EX=$R/bin/eqquasi
export EQQUASIROOT=$R
export PATH=$R/bin:$R/scripts:$PATH
export OMP_NUM_THREADS=1
LOG=$R/work/campaign.log

say() { echo "[$(date +%H:%M:%S)] $*" | tee -a "$LOG"; }

sums() {
  python3 - "$1" "$2" <<'PY' | tee -a "$LOG"
import sys, os, re, numpy as np
d, lab = sys.argv[1], sys.argv[2]
def rd(p):
    r = []
    for l in open(p):
        l = l.strip()
        if not l or l.startswith('#'): continue
        t = l.split()
        try: float(t[0])
        except ValueError: continue
        r.append([float(v) for v in t])
    return np.array(r)
try:
    g = rd(os.path.join(d, 'global.dat'))
    s = rd(os.path.join(d, 'fltst_strk+000dp+000.txt'))
    m = re.findall(r'Total elems =\s+(\d+)', open(os.path.join(d, 'r.log'), errors='ignore').read())
    print(f"RESULT {lab:16s} elems={(m[-1] if m else ''):>8s} steps={len(g)-1:6d} "
          f"t_end={g[-1,0]/86400:6.2f}d peakV={g[:,1].max():7.3f}@{g[g[:,1].argmax(),0]/86400:5.2f}d "
          f"peakM={g[:,2].max():.3e} slip={s[-1,1]*1000:7.2f}mm")
except Exception as e:
    print(f"RESULT {lab:16s} FAILED: {e}")
PY
}

runit() {   # dir label np
  cd "$1" || return 1
  ./case.setup > setup.log 2>&1
  if [ ! -f model.txt ]; then say "SETUP FAILED $2"; return 1; fi
  echo 1 > currentcycle.txt
  say "start $2 (np=$3)"
  /usr/bin/time -v mpirun -np "$3" "$EX" > r.log 2> t.log
  say "done  $2 wall=$(grep 'Elapsed (wall' t.log | awk '{print $NF}')"
  sums "$1" "$2"
}

say "--- phase 3b: xi sweep at the 500 m box"
for XI in 0.05 0.1 0.2; do
  D=$R/work/c3.xi$XI
  rm -rf "$D"; cd $R && create.newcase "$D" test.bp8.qdc.gs.10 > /dev/null
  cd "$D"
  sed -i "s/^par.xi = .*/par.xi = $XI # swept/" user_defined_params.py
  sed -i "s/^par.nstep       = .*/par.nstep       = 4000/" user_defined_params.py
  sed -i "s/^par.nt_out      = .*/par.nt_out      = 4000/" user_defined_params.py
  runit "$D" "xi_${XI}" 1
done

say "--- phase 2c: 2000 m box confirmation, nstep trimmed since the peak is at 0.84 d"
D=$R/work/c3.dom2000
rm -rf "$D"; cd $R && create.newcase "$D" test.bp8.qdc.gs.10 > /dev/null
cd "$D"
for v in fx fy fz; do
  sed -i "s/^par.${v}min, par.${v}max = .*/par.${v}min, par.${v}max = -2000.0, 2000.0/" user_defined_params.py
done
sed -i "s/^par.nstep       = .*/par.nstep       = 1500/" user_defined_params.py
sed -i "s/^par.nt_out      = .*/par.nt_out      = 1500/" user_defined_params.py
runit "$D" domain_2000m 8

cd $R && python3 script/plotDomainSweep.py >> "$LOG" 2>&1
say "===== phase 3b complete ====="
