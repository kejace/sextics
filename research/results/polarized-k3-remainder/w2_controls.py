"""Controls for the exact reduced solver at d=8 (D in J(F))."""
import sys, time, json
sys.path.insert(0, '/tmp/hodge/results/polarized-k3-remainder')
import w2k3 as W
from hodge import core
F, A, q = W.load()
fU = core.parse("(X0+Y0)**6 + 5*(X0+Y0)**5*(X1+Y1)/7 + 5*(X0+Y0)**2*(X1+Y1)**3*(X2+Y2)/7 + 15*(X0+Y0)*(X2+Y2)**5/7 + (X1+Y1)**6 + 10*(X1+Y1)**5*(X2+Y2)/7 + (X2+Y2)**6")
fV = core.parse("(X0-Y0)**6 + 5*(X0-Y0)**5*(X1-Y1)/7 + 5*(X0-Y0)**2*(X1-Y1)**3*(X2-Y2)/7 + 15*(X0-Y0)*(X2-Y2)**5/7 + (X1-Y1)**6 + 10*(X1-Y1)**5*(X2-Y2)/7 + (X2-Y2)**6")
chk = core.sub(core.scale(core.add(fU, fV), core.Fraction(1,2)), F)
print('F = (fU+fV)/2 check:', len(chk) == 0)
controls = [
  ('illegal F*X0^2', core.mul(F, core.parse('X0**2'))),
  ('fU*U0^2 + fV*V0^2', core.add(core.mul(fU, core.parse('(X0+Y0)**2')), core.mul(fV, core.parse('(X0-Y0)**2')))),
]
for name, D in controls:
    t0 = time.time()
    w = W.w_from_D_exact(F, 8, D)
    print(name, 'in J:', w is not None, f'({time.time()-t0:.1f}s)', flush=True)
    if w is None: continue
    st, data = W.exact_reduced_solve(F, A, 8, w, verbose=True, max_primes=4)
    print(name, '->', st, {k: v for k, v in data.items() if k not in ('D','N','Q')}, flush=True)
    rec = {'step': 'control', 'family': 'control', 'D': name, 'd': 8, 'status': st, 'rank': data.get('rank'), 'nullity': data.get('nullity'),
           'rows': data.get('rows'), 'cols': data.get('cols'), 'seconds': round(time.time()-t0,1), 'note': 'reduced Koszul-parametrized exact solver (w2k3.exact_reduced_solve)'}
    W.log_batch(rec)
