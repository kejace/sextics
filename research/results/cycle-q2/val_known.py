"""Validate the sector solver on accepted witnesses (must come out consistent) and on exact lift."""
import sys, time
sys.path.insert(0, '/tmp/hodge/results/cycle-q2')
import eq, modlu
from hodge import core
for iid in ['fermat-q2', 'cone-ladder-linear', 'polarized-contact-coupled-d2']:
    F, A, q, rec = core.load_instance(iid)
    w = core.load_witness(iid)
    D = w[1]
    t0 = time.time()
    S = eq.EqSystem(F, A, q, D, p=modlu.P22, verbose=True)
    r = eq.solve_fast(S, verbose=False)
    print(iid, 'q', q, 'D', core.to_str(D), '->', {k: v for k, v in r.items() if k in ('consistent', 'rank', 'nullity', 'rows', 'cols')}, f'{time.time()-t0:.1f}s', flush=True)
    # a wrong denominator for control
    Dw = core.add(D, core.parse('X0**' + str(core.degree(D))))
    S = eq.EqSystem(F, A, q, Dw, p=modlu.P22, verbose=False)
    r = eq.solve_fast(S, verbose=False)
    print('   control D + X0^d ->', {k: v for k, v in r.items() if k in ('consistent', 'rank', 'nullity')}, flush=True)
    # exact lift via big-prime path
    st, data = eq.solve_exact(F, A, q, D, max_primes=3, verbose=True)
    print('   exact lift:', st, {k: v for k, v in data.items() if k not in ('N', 'Q')}, flush=True)
