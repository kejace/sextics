"""Wave 2 job A (remote): (1) full fixed-D test of D8 (unique non-Jacobian E1-kernel element at d=8), compressed mod two
primes + exact lift attempt + served-target dimension; (2) E1 kernel (mult_kernel) at d=9, both Y-parities, bases saved.
Appends JSON lines to w2-d8d9.jsonl in the lane dir."""
import sys, json, time, os
sys.path.insert(0, '/home/kejace/devel/sextics-work')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hodge import core, fixed_d, jacobian
LANE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(LANE, 'w2-d8d9.jsonl')
def log(e):
    e.setdefault('engine', 2); e.setdefault('ts', time.strftime('%Y-%m-%dT%H:%M:%S'))
    with open(OUT, 'a') as f: f.write(json.dumps(e, default=str) + '\n')
    print(json.dumps(e, default=str)[:600], flush=True)
F, A, q, rec = core.load_instance('diagonal-permutation')
D8 = core.parse(json.load(open(os.path.join(LANE, 'D8_e1kernel.json')))['D'])
# (1) D8
for p in core.PRIMES[:2]:
    t0 = time.time()
    S = fixed_d.FixedDSystem(F, A, q, D8, p, verbose=True)
    r = fixed_d.solve_compressed(S, verbose=True)
    log({'exp': 'D8-fixedD-compressed', 'd': 8, 'p': p, 'consistent': r['consistent'], 'rank': r['rank'], 'nullity': r['nullity'],
         'rows': S.nrows, 'cols': S.ncols, 'seconds': round(time.time() - t0, 1)})
    if r['consistent']:
        t0 = time.time()
        sol = None
        for seed in range(3):
            sol = fixed_d.exact_solution_from_pivots(F, A, q, D8, r['pivots'], verbose=True, seed=seed + 7)
            if sol is not None: break
        if sol is not None:
            N, Q = sol
            ok, info = core.verify_witness(F, A, q, D8, N, Q)
            w = {'d': 8, 'D': core.to_str(D8), 'N': [core.to_str(n) for n in N], 'Q': core.to_str(Q)}
            json.dump(w, open(os.path.join(LANE, 'witness-candidate.json'), 'w'), indent=1)
            log({'exp': 'D8-exact', 'ok': ok, 'info': info, 'seconds': round(time.time() - t0, 1), 'note': 'WITNESS CANDIDATE WRITTEN'})
        else:
            log({'exp': 'D8-exact', 'ok': False, 'note': 'exact lift failed', 'seconds': round(time.time() - t0, 1)})
        break
# served-target dimension for D8 (compressed homogeneous system with target columns)
try:
    from hlib import target_image_fast
    t0 = time.time()
    dimA, basis, contains, info = target_image_fast(D8, verbose=True)
    log({'exp': 'D8-served-targets', 'served_dim': dimA, 'containsA': contains, 'info': info, 'seconds': round(time.time() - t0, 1)})
except Exception as ex:
    log({'exp': 'D8-served-targets', 'error': repr(ex)})
# (2) E1 kernel at d=9 both parities
for par in (0, 1):
    t0 = time.time()
    K = jacobian.mult_kernel(F, A, 9, parity=par, verbose=True)
    log({'exp': 'mult_kernel', 'd': 9, 'parity': par, 'dim_mod_J': K['dim'], 'kernel_dim': K['kernel_dim'], 'jacobian_dim': K['jacobian_dim'],
         'ring_dim': K['ring_dim'], 'rows': K['rows'], 'cols': K['cols'], 'seconds': round(time.time() - t0, 1)})
    if K['dim']:
        json.dump({'d': 9, 'parity': par, 'p': K['p'], 'basis_modp': [core.to_str({e: core.Fraction(v) for e, v in b.items()}) for b in K['basis']]},
                  open(os.path.join(LANE, f'w2-kernel-d9-par{par}.json'), 'w'))
print('JOB A DONE', flush=True)
