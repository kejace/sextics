import sys, time, json
from w2lib import *
from family import *
from msverdict import verdict
tags = sys.argv[1].split(',')
tmo = int(sys.argv[2]) if len(sys.argv) > 2 else 60
for tag in tags:
    c = tuple(int(x) for x in tag)
    t0 = time.time()
    info = sector_closure(8, c, f'local_d8/d8_chi{tag}', run_singular=False, verbose=False)
    if info.get('n_unknowns', 0) > 0 and info.get('n_quadrics', 0) > 0:
        v = verdict(f'local_d8/d8_chi{tag}.ms', threads=2, timeout=tmo, verbose=1)
        info['msolve'] = {k: val for k, val in v.items() if k not in ('log_tail', 'pure_powers')}
        info['result'] = 'EXCLUDED (only c=0)' if v.get('verdict', '').startswith('ONLY-ORIGIN') else v.get('verdict')
    info['seconds'] = round(time.time() - t0, 1)
    log('closure-pure-sector', **info, note='local run; gauge-fixed Hessian relaxation of the full E1-admissible family in this G_aff sector at d=8; msolve -g 1')
    print(tag, info.get('n_unknowns'), info.get('n_quadrics'), info.get('quadric_rank'), info.get('result'), info['seconds'], flush=True)
