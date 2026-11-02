"""Full fixed-D solve for D = u * w (u = E_X(F)/6, w vertical).  Usage: python3 w4_fixedD.py <tag> <w-expression in L,M,eX,eY> [method]"""
import sys, json, time, os
sys.path.insert(0, '/tmp/hodge/results/short-power-presentation/wave4')
from w4_tlib import *
from w2c_lib import c, L, M, eX, eY, ppow
from hodge import fixed_d
tag, wexpr = sys.argv[1], sys.argv[2]
method = sys.argv[3] if len(sys.argv) > 3 else 'auto'
F, A, q, rec = core.load_instance(IID)
u = core.parse(json.load(open('/tmp/hodge/results/short-power-presentation/wave4/u_poly.json'))['u'])   # 7^7 u (integral)
env = {'L': L, 'M': M, 'eX': eX, 'eY': eY, 'mul': core.mul, 'add': core.add, 'pow': ppow, 'scale': core.scale}
w = eval(wexpr, env)
D = core.mul(u, w)
d = core.degree(D); par = core.parity_of(D)
print('tag', tag, 'deg D', d, 'parity', par, 'terms', len(D), 'engine', getattr(fixed_d, 'ENGINE_VERSION', None), flush=True)
p3 = core.PRIMES[2]
ok, info = fixed_d.e1_solvable(F, A, q, D, p=p3)
divF, _ = jacobian.divisible_by(reduce_mod(F, p3), reduce_mod(D, p3), p3)
print('E1(p3)', ok, 'F|D', bool(divF), flush=True)
t0 = time.time()
out = fixed_d.solve_fixed_D_fast(F, A, q, D, verbose=True, method=method)
rec_out = {'what': 'fixedD-u*w', 'tag': tag, 'w': wexpr, 'd': d, 'parD': par, 'prime': out.get('p'), 'engine': getattr(fixed_d, 'ENGINE_VERSION', None),
           'consistent': out.get('consistent'), 'rank': out.get('rank'), 'nullity': out.get('nullity'), 'rows': out.get('rows'), 'cols': out.get('cols'),
           'solver': out.get('solver'), 'exact': bool(out.get('exact')), 'seconds': round(time.time() - t0, 1), 'file': 'fixedD_%s.json' % tag,
           'result': 'solved' if out.get('exact') else ('inconsistent' if out.get('consistent') is False else ('lift-failed' if out.get('consistent') else 'unknown'))}
log(rec_out)
json.dump({k: (v if isinstance(v, (int, float, str, bool, type(None), list)) else str(v)) for k, v in out.items() if k not in ('N', 'Q', 'x')}, open('fixedD_%s.json' % tag, 'w'))
if out.get('exact'):
    okv = core.verify_witness(F, A, q, D, out['N'], out['Q'])
    print('VERIFY', okv, flush=True)
    json.dump({'D': core.to_str(D), 'N': [core.to_str(n) for n in out['N']], 'Q': core.to_str(out['Q'])},
              open('/tmp/hodge/results/short-power-presentation/witness-candidate-w4-%s.json' % tag, 'w'))
print('DONE', json.dumps(rec_out), flush=True)
