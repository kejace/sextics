"""Exp B: find ALL reflection symmetries of a smooth hypersurface Phi (degree 6) over F_p via a
0-dimensional polynomial system solved by msolve.

A reflection s with (-1)-eigenvector v and mirror ell=0 preserves Phi iff D_v^k Phi vanishes on ker ell
for k=1,3,5, where D_v = sum v_i d_i.  The k=5 condition forces ell = ell_v := sum_i (d_iPhi)(v) z_i.
Unknowns: v=(1,v_1..v_{n-1}) (affine chart after a random linear coordinate change), lambda_i = d_iPhi(v).
Equations: lambda_i - d_iPhi(v) = 0; for charts j in {0,1,2}: coefficients in z of
  q_v(z_j -> -sum_{i!=j} lambda_i z_i, z_i -> lambda_j z_i)   for q_v = D_v Phi and D_v^3 Phi.
Usage: python3 reflsearch.py <name> <prime> <seed> <poly-expr> [subst 'Y0=l*Y1' ...]
Writes <name>.ms (msolve input), <name>.msout, <name>.json (summary).
"""
import sys, json, time, random, subprocess, re, itertools, shutil, os
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core
import flint

name, p, seed = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
expr = sys.argv[4]
substs = sys.argv[5:]
random.seed(seed)
t0 = time.time()
VARS6 = ['X0','X1','X2','Y0','Y1','Y2']
P6 = core.parse(expr)                         # dict poly over Q
# apply substitutions var = lam * var2 (pencil fibre), lam random in F_p
ctx6 = flint.nmod_mpoly_ctx.get(tuple(VARS6), ordering=flint.Ordering.lex, modulus=p)
def to_mpoly(poly, ctx, nvars):
    d = {}
    for m, cf in poly.items():
        cf = core.reduce_mod(cf, p) if not isinstance(cf, int) else cf % p
        d[tuple(m[:nvars])] = cf
    return ctx.from_dict(d)
# reduce coefficients mod p
Pm = {}
for m, cf in P6.items():
    num, den = cf.numerator, cf.denominator
    Pm[m] = (num * pow(den, -1, p)) % p
Phi = to_mpoly(Pm, ctx6, 6)
gens6 = ctx6.gens()
sub_info = []
for s in substs:
    a, b = s.split('=')
    lam = random.randrange(2, p - 1)
    i = VARS6.index(a.strip()); j = VARS6.index(b.strip().split('*')[-1].strip())
    Phi = Phi.compose(*[gens6[k] if k != i else lam * gens6[j] for k in range(6)])
    sub_info.append((a, b, lam))
# determine the active variables (those actually appearing)
active = [k for k in range(6) if Phi.degrees()[k] > 0]
n = len(active)
VARS = [VARS6[k] for k in active]
ctx = flint.nmod_mpoly_ctx.get(tuple(VARS), ordering=flint.Ordering.lex, modulus=p)
gens = ctx.gens()
Phi_n = ctx.from_dict({tuple(m[k] for k in active): cf for m, cf in Phi.to_dict().items()})
# random linear coordinate change z = G w (G random invertible mod p)
while True:
    G = [[random.randrange(p) for _ in range(n)] for _ in range(n)]
    if flint.nmod_mat(G, p).rank() == n:
        break
Phi_w = Phi_n.compose(*[sum((G[a][b] * gens[b] for b in range(n)), ctx.from_dict({})) for a in range(n)])
# unknowns: v_1..v_{n-1}, lambda_0..lambda_{n-1}  (v_0 = 1)
UN = [f'v{i}' for i in range(1, n)] + [f'l{i}' for i in range(n)]
ALL = VARS + UN
ctxA = flint.nmod_mpoly_ctx.get(tuple(ALL), ordering=flint.Ordering.lex, modulus=p)
gA = ctxA.gens()
z = gA[:n]; v = [ctxA.from_dict({(0,)*len(ALL): 1})] + list(gA[n:n + n - 1]); lam = list(gA[2 * n - 1:])
def lift(poly_ctx):  # Phi_w in ctx -> ctxA
    return ctxA.from_dict({m + (0,) * (len(ALL) - n): cf for m, cf in poly_ctx.to_dict().items()})
PhiA = lift(Phi_w)
dPhi = [PhiA.derivative(i) for i in range(n)]
# D_v Phi, D_v^3 Phi
def Dv(q):
    return sum((v[i] * q.derivative(i) for i in range(n)), ctxA.from_dict({}))
q1 = Dv(PhiA); q3 = Dv(Dv(q1))
eqs = []
# lambda_i = d_iPhi(v)
for i in range(n):
    e = dPhi[i].compose(*(v + [gA[k] for k in range(n, len(ALL))]))
    eqs.append(lam[i] - e)
zero = ctxA.from_dict({})
NCHARTS = int(os.environ.get('REFL_CHARTS', '3'))
K1_ONLY = os.environ.get('REFL_K1_ONLY', '0') == '1'
for j in range(NCHARTS):
    img = []
    for i in range(n):
        if i == j:
            img.append(-sum((lam[k] * z[k] for k in range(n) if k != j), zero))
        else:
            img.append(lam[j] * z[i])
    img += [gA[k] for k in range(n, len(ALL))]
    for q in ((q1,) if K1_ONLY else (q1, q3)):
        qs = q.compose(*img)
        # collect coefficients of z-monomials
        coeffs = {}
        for m, cf in qs.to_dict().items():
            zm = m[:n]; rest = m[n:]
            coeffs.setdefault(zm, {})[rest] = cf
        for zm, d in coeffs.items():
            eqs.append(ctxA.from_dict({(0,) * n + r: cf for r, cf in d.items()}))
# drop zero eqs, dedupe
eqs = [e for e in eqs if e != zero]
print(f'[{name}] n={n} unknowns={len(UN)} equations={len(eqs)} charts={NCHARTS} k1only={K1_ONLY} built in {time.time()-t0:.1f}s', flush=True)
# write msolve input (only unknown variables)
def ms_str(e):
    s = []
    for m, cf in e.to_dict().items():
        term = [str(cf)]
        for k, ex in enumerate(m):
            if ex:
                term.append(f'{ALL[k]}^{ex}' if ex > 1 else ALL[k])
        s.append('*'.join(term))
    return '+'.join(s) if s else '0'
with open(f'{name}.ms', 'w') as f:
    f.write(','.join(UN) + '\n' + str(p) + '\n')
    f.write(',\n'.join(ms_str(e) for e in eqs) + '\n')
t1 = time.time()
if shutil.which('msolve') is None:
    print(f'[{name}] msolve not found; wrote {name}.ms only'); sys.exit(0)
cmd = ['msolve', '-v', '1', '-t', os.environ.get('MSOLVE_THREADS', '8'), '-P', '1', '-f', f'{name}.ms', '-o', f'{name}.param']
r = subprocess.run(cmd, capture_output=True, text=True)
out = open(f'{name}.param').read() if r.returncode == 0 else ''
print(f'[{name}] msolve rc={r.returncode} in {time.time()-t1:.1f}s; stderr tail: {r.stderr[-400:]}', flush=True)
# parse: msolve prints [dim, [deg, nvars, ... ]] for 0-dim; dimension -1 means empty
summary = {'name': name, 'p': p, 'seed': seed, 'n': n, 'vars': VARS, 'subst': sub_info, 'G': G,
           'n_eqs': len(eqs), 'msolve_rc': r.returncode, 'seconds': time.time() - t0}
head = out[:200].replace('\n', ' ')
summary['msout_head'] = head
m = re.match(r'\s*\[\s*(-?\d+)\s*,', out)
if m:
    summary['dim'] = int(m.group(1))
    if summary['dim'] == 0:
        m2 = re.match(r'\s*\[\s*0\s*,\s*\[\s*(\d+)\s*,\s*(\d+)\s*,\s*(\d+)\s*,', out)
        if m2:
            summary['degree'] = int(m2.group(3))
summary['param_rc'] = r.returncode
print(json.dumps(summary), flush=True)
json.dump(summary, open(f'{name}.json', 'w'))
