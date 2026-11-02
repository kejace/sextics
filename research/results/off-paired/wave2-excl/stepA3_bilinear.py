"""A3: bilinear (D,N) sector systems, EXACT over Q, decided by Singular std over QQ (rigorous unit-ideal test)
and independently by msolve mod a new prime.  Usage: python3 stepA3_bilinear.py d [a b ...] [--mode qq|modp|both] [--threads k] [--timeout s]"""
import sys, os, json, time, subprocess
sys.path.insert(0, '/tmp/hodge/results/off-paired/wave2-excl')
from exlib import *
from fractions import Fraction

args = [a for a in sys.argv[1:] if not a.startswith('--')]
opts = {a.split('=')[0][2:]: a.split('=')[1] for a in sys.argv[1:] if a.startswith('--')}
mode = opts.get('mode', 'both'); threads = int(opts.get('threads', '1')); timeout = int(opts.get('timeout', '3600'))
PM = int(opts.get('prime', str(PRIMES[0])))
d = int(args[0])
chars = [(int(x[0]), int(x[1])) for x in args[1:]] if len(args) > 1 else [(a, b) for a in range(6) for b in range(6)]
F, A, q = load()
OUT = f'/tmp/hodge/results/off-paired/wave2-excl/bilinear_d{d}'
os.makedirs(OUT, exist_ok=True)
for c in chars:
    tag = f'd{d}_c{c[0]}{c[1]}'
    sysd, info = build_bilinear(F, A, q, d, c, p=None, verbose=True)
    if sysd is None:
        rec = {'step': 'A3', 'family': 'bilinear-sector', 'd': d, 'char': list(c), 'result': 'empty-colon-after-gauge', **info}
        print(json.dumps(rec)); log(rec); continue
    rec = {'step': 'A3', 'family': 'bilinear-sector', 'd': d, 'char': list(c), 'nt': info['nt'], 'ns': info['ns'], 'neq': info['neq'],
           'colon_dim': info['colon_dim'], 'gauge_dim': info['gauge_dim'], 'exact_build': True}
    json.dump({'nt': sysd['nt'], 'ns': sysd['ns'], 'eqs': [{f'{vname(a)}*{vname(b)}': str(v) for (a, b), v in eq.items()} for eq in sysd['eqs']],
               'D_basis': [core.to_str(D) for D in sysd['Ds']]}, open(f'{OUT}/{tag}_system_Q.json', 'w'))
    if mode in ('qq', 'both'):
        path = f'{OUT}/{tag}_QQ.sing'
        write_singular(sysd, path, p=0, charts=True)
        t0 = time.time()
        try:
            r = subprocess.run(['Singular', '-q', path], capture_output=True, text=True, timeout=timeout)
            res = r.stdout.strip().replace('\n', ' ')
        except subprocess.TimeoutExpired:
            res = 'TIMEOUT'
        rec['singular_QQ'] = res; rec['singular_QQ_seconds'] = round(time.time()-t0, 1)
        rec['rigorous_excluded'] = ('ALL-UNIT' in res)
    if mode in ('modp', 'both'):
        paths = write_msolve_charts(sysd, f'{OUT}/{tag}_p{PM}', PM)
        sts = []; tot = 0.0
        for path in paths:
            st, dt, res = run_msolve(path, path.replace('.ms', '.out'), threads=threads, timeout=timeout)
            sts.append(st); tot += dt
            if st not in ('no-solution', 'timeout'): rec.setdefault('msolve_output', []).append(res[:500])
        rec['msolve_charts'] = sts; rec['msolve_prime'] = PM; rec['msolve_seconds'] = round(tot, 1)
        rec['msolve_all_empty'] = all(s_ == 'no-solution' for s_ in sts)
    print(json.dumps(rec)); sys.stdout.flush(); log(rec)
