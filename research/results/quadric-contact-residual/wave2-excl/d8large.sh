#!/bin/bash
# Job excl-d8large: the four remaining degree-8 pure G_aff sectors (26-29 unknowns, 52 quadrics): msolve -g 1, 16 threads, 3 h cap each.
cd /tmp/hodge/results/quadric-contact-residual/wave2-excl
export HODGE_W2LOG=/tmp/hodge/results/quadric-contact-residual/wave2-excl/excl-d8large.log.jsonl
mkdir -p local_d8
for tag in 1100 1101 0001 0000; do
  python3 - "$tag" <<'PY' 2>&1 | grep -v "engine v2" | tail -2
import sys, time, json
from w2lib import *
from family import *
from msverdict import verdict
tag = sys.argv[1]; c = tuple(int(x) for x in tag)
t0 = time.time()
info = sector_closure(8, c, f'local_d8/d8_chi{tag}', run_singular=False, verbose=False)
v = verdict(f'local_d8/d8_chi{tag}.ms', threads=16, timeout=10800, verbose=1, logfile=f'local_d8/d8_chi{tag}.msolve.log')
info['msolve'] = {k: val for k, val in v.items() if k not in ('log_tail', 'pure_powers')}
info['result'] = 'EXCLUDED (only c=0)' if v.get('verdict', '').startswith('ONLY-ORIGIN') else v.get('verdict')
info['seconds'] = round(time.time() - t0, 1)
log('closure-pure-sector', **info, note='remote run; gauge-fixed Hessian relaxation of the full E1-admissible family in this G_aff sector at d=8; msolve -g 1, 16 threads')
PY
done
echo D8LARGE DONE
