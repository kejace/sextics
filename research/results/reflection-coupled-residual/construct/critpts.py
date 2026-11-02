"""Exp F: critical points of the pencil map Y0/Y1 on X (singular fibres of the (Y0:Y1)-pencil):
solve {d0F = d1F = d2F = d5F = 0, F = 0} in P^5 over F_p with msolve (0-dim, Bezout 3750),
then classify solutions by Klein fixed-point type and by the value lambda = Y0/Y1.
Usage: python3 critpts.py <name> <prime>"""
import sys, json, time, subprocess, ast, re
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core
import flint
name, p = sys.argv[1], int(sys.argv[2])
t0 = time.time()
F6, A, q, rec = core.load_instance('reflection-coupled-residual')
Fm = core.reduce_mod(F6, p)
V = ['X0','X1','X2','Y0','Y1','Y2']
def ms(poly):
    terms = []
    for e, c in poly.items():
        t = [str(c)]
        for k, ex in enumerate(e):
            if ex: t.append(f'{V[k]}^{ex}' if ex > 1 else V[k])
        terms.append('*'.join(t))
    return '+'.join(terms) if terms else '0'
# generic affine chart: random linear form = 1 to avoid missing points
import random
random.seed(11)
lin = {tuple(1 if k == i else 0 for k in range(6)): random.randrange(1, p) for i in range(6)}
eqs = [core.pderiv(Fm, i, p) for i in (0, 1, 2, 5)] + [Fm]
eqs_str = [ms(e) for e in eqs] + [ms(lin) + '-1']
with open(f'{name}.ms', 'w') as f:
    f.write(','.join(V) + '\n' + str(p) + '\n' + ',\n'.join(eqs_str) + '\n')
r = subprocess.run(['msolve', '-v', '1', '-t', '6', '-P', '1', '-f', f'{name}.ms', '-o', f'{name}.param'], capture_output=True, text=True)
print(f'[{name}] msolve rc={r.returncode} {time.time()-t0:.1f}s', flush=True)
txt = open(f'{name}.param').read().replace('\n', ' ').rstrip(': ')
data = ast.literal_eval(txt.replace("'", '"'))
dim, body = data
print('dim', dim)
if dim != 0: sys.exit()
pp, nv, deg, vars_, linform, par = body
one, [[dw, wc], [dd, dc], params] = par
w = flint.nmod_poly(wc, p); dwp = w.derivative(); den = flint.nmod_poly(dc, p)
fac = w.factor()
degs = {}
for f, m in fac[1]: degs[f.degree()] = degs.get(f.degree(), 0) + m
print(f'[{name}] elim degree {w.degree()} (expected 3750), squarefree part {w.degree()-w.gcd(dwp).degree()}, factor degree multiset {sorted(degs.items())}', flush=True)
# F_p-rational points: decode and classify
out = {'name': name, 'p': p, 'degree': w.degree(), 'sqfree': w.degree()-w.gcd(dwp).degree(), 'factor_degrees': sorted(degs.items()), 'seconds': time.time()-t0}
rat = []
for f, m in fac[1]:
    if f.degree() != 1: continue
    th = -f[0] / f[1]
    wp = dwp(th); dn = den(th)
    vals = {}
    for i, pr in enumerate(params):
        Pi = flint.nmod_poly(pr[0][1], p)
        vals[vars_[i]] = Pi(th)
    vals[vars_[-1]] = th
    # fix convention by checking F=0 for candidates
    def evalF(pt):
        s = flint.nmod(0, p)
        for e, c in Fm.items():
            t = flint.nmod(c, p)
            for k, ex in enumerate(e):
                if ex: t = t * pt[k]**ex
            s += t
        return s
    cands = []
    for label, tr in [('raw', lambda x: x), ('-raw/wp', lambda x: -x / wp), ('raw/wp', lambda x: x / wp), ('-raw/dn', lambda x: -x / dn), ('raw/dn', lambda x: x / dn)]:
        pt = [tr(vals[v]) if v != vars_[-1] else vals[v] for v in V]
        if evalF(pt) == 0:
            cands.append((label, pt)); break
    if not cands:
        continue
    label, pt = cands[0]
    typ = []
    if pt[5] == 0: typ.append('Y2=0 (r-fixed)')
    if pt[3] == 0 and pt[4] == 0: typ.append('Y0=Y1=0 (rho-fixed)')
    if pt[0] == 0 and pt[1] == 0 and pt[2] == 0: typ.append('X=0 (tau-fixed)')
    lam = str(pt[3] / pt[4]) if pt[4] != 0 else 'inf'
    rat.append({'pt': [str(x) for x in pt], 'type': typ or ['free'], 'lambda': lam, 'conv': label, 'mult': m})
out['n_rational'] = len(rat)
from collections import Counter
out['rational_types'] = Counter(tuple(r['type']) for r in rat).most_common()
out['rational_points'] = rat[:400]
json.dump(out, open(f'{name}.json', 'w'), indent=1)
print(json.dumps({k: v for k, v in out.items() if k != 'rational_points'}), flush=True)
