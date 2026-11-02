"""M1 (wave 4): all reflections of a smooth degree-6 hypersurface Phi over F_p, lambda-eliminated formulation.
s_v = reflection with (-1)-eigenvector v and mirror ell_v = grad Phi(v) (forced by the degree-5 condition).
s_v preserves Phi  <=>  D_v Phi and D_v^3 Phi vanish on ker ell_v  (D_v^5 Phi = 5! ell_v automatically).
Unknowns: v = (1, v1..v_{n-1}) in an affine chart after a random linear change G.  No auxiliary lambdas:
lambda_i := d_i Phi(v) are substituted, so the equations are polynomials in v of degree <= 26 (k=1) and 18 (k=3),
in only n-1 unknowns (cheap Macaulay matrices).  Charts j in {0,1,2} (z_j -> -sum_{i!=j} lambda_i z_i, z_i -> lambda_j z_i).
Usage: python3 reflsearch2.py <name> <prime> <seed> <poly-expr> [subst 'Y0=l*Y1' ...]   env: REFL_CHARTS (3), REFL_K1_ONLY (0), MSOLVE_THREADS (4)
Writes <name>.ms, <name>.param, <name>.json (summary incl. elimination-polynomial factor degrees and verified reflections).
"""
import sys, json, time, random, subprocess, re, os, shutil, ast
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core
import flint

name, p, seed = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
expr = sys.argv[4]; substs = sys.argv[5:]
random.seed(seed); t0 = time.time()
VARS6 = ['X0', 'X1', 'X2', 'Y0', 'Y1', 'Y2']
P6 = core.parse(expr)
Pm = {m: (cf.numerator * pow(cf.denominator, -1, p)) % p for m, cf in P6.items()}
ctx6 = flint.nmod_mpoly_ctx.get(tuple(VARS6), ordering=flint.Ordering.lex, modulus=p)
Phi = ctx6.from_dict(Pm); gens6 = ctx6.gens()
sub_info = []
for s in substs:
    a, b = s.split('=')
    lam = random.randrange(2, p - 1)
    i = VARS6.index(a.strip()); j = VARS6.index(b.strip().split('*')[-1].strip())
    Phi = Phi.compose(*[gens6[k] if k != i else lam * gens6[j] for k in range(6)])
    sub_info.append((a, b, lam))
active = [k for k in range(6) if Phi.degrees()[k] > 0]
n = len(active); VARS = [VARS6[k] for k in active]
ctx = flint.nmod_mpoly_ctx.get(tuple(VARS), ordering=flint.Ordering.lex, modulus=p); gens = ctx.gens()
Phi_n = ctx.from_dict({tuple(m[k] for k in active): cf for m, cf in Phi.to_dict().items()})
while True:
    G = [[random.randrange(p) for _ in range(n)] for _ in range(n)]
    if flint.nmod_mat(G, p).rank() == n: break
Phi_w = Phi_n.compose(*[sum((G[a][b] * gens[b] for b in range(n)), ctx.from_dict({})) for a in range(n)])
UN = [f'v{i}' for i in range(1, n)]
ALL = VARS + UN
ctxA = flint.nmod_mpoly_ctx.get(tuple(ALL), ordering=flint.Ordering.lex, modulus=p); gA = ctxA.gens()
z = gA[:n]; v = [ctxA.from_dict({(0,) * len(ALL): 1})] + list(gA[n:])
PhiA = ctxA.from_dict({m + (0,) * (n - 1): cf for m, cf in Phi_w.to_dict().items()})
dPhi = [PhiA.derivative(i) for i in range(n)]
zero = ctxA.from_dict({})
# lambda_i(v) = d_i Phi(v): compose z -> v
lam = [dPhi[i].compose(*(v + list(gA[n:]))) for i in range(n)]
def Dv(q): return sum((v[i] * q.derivative(i) for i in range(n)), zero)
q1 = Dv(PhiA); q3 = Dv(Dv(q1))
NCHARTS = int(os.environ.get('REFL_CHARTS', '3')); K1_ONLY = os.environ.get('REFL_K1_ONLY', '0') == '1'
eqs = []
for j in range(NCHARTS):
    img = [(-sum((lam[k] * z[k] for k in range(n) if k != j), zero)) if i == j else lam[j] * z[i] for i in range(n)] + list(gA[n:])
    for q in ((q1,) if K1_ONLY else (q1, q3)):
        qs = q.compose(*img)
        coeffs = {}
        for m, cf in qs.to_dict().items():
            coeffs.setdefault(m[:n], {})[m[n:]] = cf
        for zm, d in coeffs.items():
            e = ctxA.from_dict({(0,) * n + r: cf for r, cf in d.items()})
            if e != zero: eqs.append(e)
# dedupe
seen = set(); eqs2 = []
for e in eqs:
    key = str(e)
    if key not in seen: seen.add(key); eqs2.append(e)
eqs = eqs2
maxdeg = int(max(e.total_degree() for e in eqs))
print(f'[{name}] n={n} unknowns={len(UN)} equations={len(eqs)} maxdeg={maxdeg} charts={NCHARTS} k1only={K1_ONLY} built in {time.time()-t0:.1f}s', flush=True)
def ms_str(e):
    s = []
    for m, cf in e.to_dict().items():
        term = [str(cf)]
        for k, ex in enumerate(m):
            if ex: term.append(f'{ALL[k]}^{ex}' if ex > 1 else ALL[k])
        s.append('*'.join(term))
    return '+'.join(s) if s else '0'
with open(f'{name}.ms', 'w') as f:
    f.write(','.join(UN) + '\n' + str(p) + '\n'); f.write(',\n'.join(ms_str(e) for e in eqs) + '\n')
summary = {'name': name, 'p': p, 'seed': seed, 'n': n, 'vars': VARS, 'subst': sub_info, 'G': G, 'n_eqs': len(eqs), 'maxdeg': maxdeg, 'charts': NCHARTS, 'k1only': K1_ONLY}
json.dump(summary, open(f'{name}.json', 'w'))
t1 = time.time()
cmd = ['msolve', '-v', '1', '-t', os.environ.get('MSOLVE_THREADS', '4'), '-P', '1', '-f', f'{name}.ms', '-o', f'{name}.param']
r = subprocess.run(cmd, capture_output=True, text=True)
summary['msolve_rc'] = r.returncode; summary['msolve_seconds'] = round(time.time() - t1, 1); summary['stderr_tail'] = r.stderr[-300:]
out = open(f'{name}.param').read() if os.path.exists(f'{name}.param') else ''
print(f'[{name}] msolve rc={r.returncode} in {time.time()-t1:.1f}s', flush=True)
m = re.match(r'\s*\[\s*(-?\d+)\s*,', out)
if m: summary['dim'] = int(m.group(1))
if summary.get('dim') == 0:
    data = ast.literal_eval(out.replace('\n', ' ').rstrip(': ').replace("'", '"'))
    pp, nv, deg, vars_, linform, par = data[1]
    summary['degree'] = int(deg)
    if deg > 0:
        one, [[dw, wc], [dd, dc], params] = par
        w = flint.nmod_poly(wc, p); dwp = w.derivative(); den = flint.nmod_poly(dc, p)
        fac = w.factor(); degs = sorted(int(f.degree()) for f, mm in fac[1]); mults = sorted(int(mm) for f, mm in fac[1])
        summary['elim_degree'] = int(w.degree()); summary['factor_degrees'] = degs; summary['factor_mults'] = mults
        # decode rational (F_p) solutions and verify reflections exactly over F_p
        Ginv = flint.nmod_mat(G, p).inv()
        def polyval(P, pt):
            val = 0
            for mono, cf in P.to_dict().items():
                t = int(cf)
                for a, e in enumerate(mono):
                    if e: t = t * pow(int(pt[a]), e, p) % p
                val = (val + t) % p
            return val
        genuine = []; spurious = 0; rational = 0
        for f, mm in fac[1]:
            if f.degree() != 1: continue
            theta = int(-f[0] / f[1]); rational += 1
            wp = int(dwp(theta)); dn = int(den(theta))
            raw = [int(sum(int(c) * pow(theta, kk, p) for kk, c in enumerate(pr[0][1])) % p) for pr in params]
            # msolve convention: x_i = -raw_i / w'(theta) (common); try the alternatives and keep the one satisfying all equations
            cand = []
            for lab, tr in [('-raw/wp', lambda x: (-x * pow(wp, -1, p)) % p if wp else None), ('raw/wp', lambda x: (x * pow(wp, -1, p)) % p if wp else None), ('raw', lambda x: x % p), ('-raw/dn', lambda x: (-x * pow(dn, -1, p)) % p if dn else None)]:
                vals = [tr(x) for x in raw]
                if any(x is None for x in vals): continue
                vals = vals + [theta]  # eliminating variable is the last one if linform = last var
                if list(linform) == [0] * (len(vars_) - 1) + [1]:
                    pt = vals
                else:
                    continue
                vv = [1] + pt[:n - 1]
                ok = all(polyval(e, [0] * n + vv[1:]) == 0 for e in eqs[:40])
                if ok: cand.append((lab, vv)); break
            if not cand: spurious += 1; continue
            lab, vv = cand[0]
            lamv = [polyval(dPhi[i], vv + [0] * (n - 1)) for i in range(n)]
            ellv = sum(lamv[i] * vv[i] for i in range(n)) % p
            if ellv == 0: spurious += 1; continue
            ok = True
            for trial in range(5):
                zz = [random.randrange(p) for _ in range(n)]
                ellz = sum(lamv[i] * zz[i] for i in range(n)) % p
                sz = [(zz[i] - 2 * ellz * pow(ellv, -1, p) * vv[i]) % p for i in range(n)]
                if polyval(PhiA, sz + [0] * (n - 1)) != polyval(PhiA, zz + [0] * (n - 1)): ok = False; break
            if ok:
                vorig = [int(sum(G[a][b] * vv[b] for b in range(n)) % p) for a in range(n)]
                genuine.append({'v_chart': vv, 'v_orig': vorig})
            else: spurious += 1
        summary['rational_solutions'] = rational; summary['genuine_Fp_reflections'] = len(genuine); summary['spurious_Fp'] = spurious
        summary['genuine'] = genuine
summary['seconds'] = round(time.time() - t0, 1)
print(json.dumps({k: vv for k, vv in summary.items() if k not in ('G', 'genuine')}), flush=True)
for g in summary.get('genuine', []): print('  reflection eigvec (orig coords):', g['v_orig'])
json.dump(summary, open(f'{name}.json', 'w'))
