"""Fixed-D test for the discriminant components D_± = 27 R^2 ∓ 2 L^3 (degree 12, cones over surfaces in P^3_{X0X1Y0Y1}),
using the index-2 subgroup G' (D_± are G'-eigenvectors) and the compressed-row solver."""
import sys, time, json
sys.path.insert(0, '/tmp/hodge/results/off-paired')
import oplib
from hodge import core, fixed_d
from hodge.core import mul, add, sub, scale, power
F, A, q = oplib.load()
R = core.parse('X0**6 + X1**6 + Y0**6 + Y1**6 + X0**3*X1*Y0*Y1/2')
L = core.parse('X1**3*Y1/2 + X0**2*X1*Y0/3')
which = sys.argv[1] if len(sys.argv) > 1 else 'minus'
D = sub(scale(mul(R, R), 27), scale(power(L, 3), 2)) if which == 'minus' else add(scale(mul(R, R), 27), scale(power(L, 3), 2))
oplib.set_group('Gp')
print('D_%s: degree %d, terms %d, G\'-char %s' % (which, core.degree(D), len(D), oplib.char_of(D)), flush=True)
t0 = time.time()
ok, info = oplib.e1_reduced(F, A, q, D)
print('E1:', ok, info, f'{time.time()-t0:.1f}s', flush=True)
if ok:
    t0 = time.time()
    S = oplib.ReducedSystem(F, A, q, D, verbose=True)
    r = fixed_d.solve_compressed(S, k_extra=60, verbose=True)
    dt = time.time() - t0
    print('compressed solve:', {k: v for k, v in r.items() if k not in ('x', 'pivots')}, f'{dt:.1f}s', flush=True)
    oplib.log(f'27*R**2 {"-" if which=="minus" else "+"} 2*L**3 (R = X0^6+X1^6+Y0^6+Y1^6+X0^3X1Y0Y1/2, L = X1^3Y1/2+X0^2X1Y0/3)', 12, 'genus2-discriminant-component', True, r['consistent'], r['rank'], r['nullity'], dt, f"G'-reduced compressed solve rows={S.nrows} cols={S.ncols}")
    if r['consistent']:
        print('!!!!! CONSISTENT — attempting exact lift', flush=True)
        st, data = oplib.lift_solution(F, A, q, D, verbose=True)
        print(st, flush=True)
        if st == 'solved':
            json.dump({'D': core.to_str(D), 'N': [core.to_str(n) for n in data['N']], 'Q': core.to_str(data['Q'])}, open('witness-candidate.json', 'w'))
else:
    oplib.log(f'27*R**2 {"-" if which=="minus" else "+"} 2*L**3', 12, 'genus2-discriminant-component', False, None, info['rank'], 0, time.time()-t0, "E1 fails (G'-reduced)")
