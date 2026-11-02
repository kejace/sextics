"""E1 filter for structured candidates with cached Jacobian reducers (one reducer per (degree, parity)).
usage: python3 e1filter.py cands.jsonl out.jsonl   (lines {"D": expr, "family": tag});  logs e1 pass/fail; passers also get
the full fixed-D solve (fixed_d.solve_fixed_D_fast) when d <= DMAX_FULL."""
import sys, json, time
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core, jacobian, fixed_d
LANE = '/tmp/hodge/results/diagonal-permutation'
DMAX_FULL = 10
F, A, q, rec = core.load_instance('diagonal-permutation')
p = core.PRIMES[0]
reducers = {}
def red(deg, par):
    if (deg, par) not in reducers:
        t0 = time.time()
        reducers[(deg, par)] = jacobian.jacobian_ring_basis(F, deg, p, par)[1]
        print(f'[reducer] deg {deg} par {par} built {time.time()-t0:.1f}s', file=sys.stderr, flush=True)
    return reducers[(deg, par)]
def log(e):
    e.setdefault('engine', 2); e.setdefault('ts', time.strftime('%Y-%m-%dT%H:%M:%S')); e['lens'] = 'search'
    with open(sys.argv[2], 'a') as f: f.write(json.dumps(e, default=str) + '\n')
    with open(f'{LANE}/wave2-search.log.jsonl', 'a') as f: f.write(json.dumps(e, default=str) + '\n')
done = set()
try:
    for l in open(sys.argv[2]): done.add(json.loads(l)['D'])
except FileNotFoundError: pass
for l in open(sys.argv[1]):
    c = json.loads(l)
    if c['D'] in done: continue
    t0 = time.time()
    try:
        D = core.parse(c['D'])
    except Exception as ex:
        log({'exp': 'e1filter', **c, 'error': repr(ex)}); continue
    d = core.degree(D); par = core.parity_of(D)
    e = {'exp': 'e1filter', 'D': c['D'], 'family': c.get('family'), 'd': d, 'parity': par, 'terms': len(D)}
    if par is None or not core.is_homogeneous(D) or not D:
        e['note'] = 'mixed parity / inhomogeneous / zero'; log(e); continue
    divF, _ = jacobian.divisible_by(F, D, p) if d >= 6 else (False, None)
    if divF:
        e['note'] = 'F | D'; log(e); continue
    AD = core.mul(core.reduce_mod(A, p), core.reduce_mod(D, p), p)
    nf = red(d + 6, (par + 1) % 2).normal_form(AD)
    e['e1'] = (len(nf) == 0)
    # also: is D itself in J (Jacobian type)?
    nfD = red(d, par).normal_form(core.reduce_mod(D, p))
    e['D_in_J'] = (len(nfD) == 0)
    if e['e1'] and d <= DMAX_FULL:
        out = fixed_d.solve_fixed_D_fast(F, A, q, D, verbose=False)
        e.update({'consistent': out['consistent'], 'rank': out['rank'], 'nullity': out['nullity'], 'rows': out['rows'], 'cols': out['cols']})
        if out['consistent'] and out.get('exact'):
            okv, vinfo = core.verify_witness(F, A, q, D, out['N'], out['Q'])
            e['verify'] = vinfo; e['note'] = 'SOLVED-EXACT' if okv else 'lift-bad'
            if okv:
                json.dump({'d': d, 'D': core.to_str(D), 'N': [core.to_str(n) for n in out['N']], 'Q': core.to_str(out['Q'])},
                          open(f'{LANE}/witness-candidate.json', 'w'), indent=1)
                print('WITNESS FOUND', file=sys.stderr, flush=True)
    e['seconds'] = round(time.time() - t0, 2)
    log(e)
print('E1FILTER DONE', flush=True)
