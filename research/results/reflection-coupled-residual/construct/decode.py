"""Decode msolve -P 1 output of reflsearch.py, verify which solutions are genuine reflections (over F_p and F_{p^2}),
and map eigenvectors back to the original coordinates.  Usage: python3 decode.py <name> <poly-expr> [substs...]"""
import sys, json, ast, time, random
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core
import flint
name = sys.argv[1]; expr = sys.argv[2]; substs = sys.argv[3:]
S = json.load(open(f'{name}.json'))
p = S['p']; G = S['G']; n = S['n']; VARS = S['vars']
VARS6 = ['X0','X1','X2','Y0','Y1','Y2']
P6 = core.parse(expr)
Pm = {m: (cf.numerator * pow(cf.denominator, -1, p)) % p for m, cf in P6.items()}
ctx6 = flint.nmod_mpoly_ctx.get(tuple(VARS6), ordering=flint.Ordering.lex, modulus=p)
Phi = ctx6.from_dict(Pm); g6 = ctx6.gens()
for (a, b, lam) in S['subst']:
    i = VARS6.index(a.strip()); j = VARS6.index(b.strip().split('*')[-1].strip())
    Phi = Phi.compose(*[g6[k] if k != i else lam * g6[j] for k in range(6)])
active = [VARS6.index(v) for v in VARS]
ctx = flint.nmod_mpoly_ctx.get(tuple(VARS), ordering=flint.Ordering.lex, modulus=p)
gens = ctx.gens()
Phi_n = ctx.from_dict({tuple(m[k] for k in active): cf for m, cf in Phi.to_dict().items()})
Phi_w = Phi_n.compose(*[sum((G[a][b] * gens[b] for b in range(n)), ctx.from_dict({})) for a in range(n)])
dPhi_w = [Phi_w.derivative(i) for i in range(n)]
# parse param file
txt = open(f'{name}.param').read().replace('\n', ' ').rstrip(': ')
data = ast.literal_eval(txt.replace("'", '"'))
dim, body = data
if dim != 0:
    print('dim', dim); sys.exit()
pp, nv, deg, vars_, linform, par = body
one, [[dw, wc], [dd, dc], params] = par
w = flint.nmod_poly(wc, p); dwp = w.derivative(); den = flint.nmod_poly(dc, p)
fac = w.factor()
degs = sorted(f.degree() for f, m in fac[1])
print(f'[{name}] elim poly degree {w.degree()}, squarefree part {w.degree()-w.gcd(dwp).degree()}, factor degrees {degs}')
# work in F_{p^k} for k = each factor degree <= 2 (others reported only)
def eval_point(theta, Fctx):
    # returns dict var -> value in F_{p^k}; convention tested below
    vals = {}
    for i, pr in enumerate(params):
        Pi = pr[0][1]
        vals[vars_[i]] = sum(Fctx(c) * theta**k for k, c in enumerate(Pi))
    return vals
results = {'rational': 0, 'quadratic': 0, 'higher': 0, 'genuine': 0, 'spurious': 0, 'eigvecs': []}
Gm = flint.nmod_mat(G, p)
Ginv = Gm.inv()
for f, m in fac[1]:
    k = f.degree()
    if k > 2:
        results['higher'] += k; continue
    Fctx = flint.fq_default_ctx(p, k) if k > 1 else None
    if k == 1:
        thetas = [(-f[0] / f[1])]
        conv = lambda c: flint.nmod(int(c), p)
        to_f = lambda x: flint.nmod(int(x), p)
    else:
        Fctx = flint.fq_default_ctx(p, k)
        pctx = flint.fq_default_poly_ctx(Fctx)
        fx = pctx([Fctx(int(c)) for c in f.coeffs()])
        thetas = [r for r, mm in fx.roots()]
        conv = lambda c: Fctx(int(c))
    for theta in thetas:
        wp = sum(conv(int(c)) * theta**i for i, c in enumerate(dwp.coeffs()))
        dn = sum(conv(int(c)) * theta**i for i, c in enumerate(den.coeffs()))
        raw = {}
        for i, pr in enumerate(params):
            Pi = pr[0][1]
            raw[vars_[i]] = sum(conv(int(c)) * theta**kk for kk, c in enumerate(Pi))
        raw[vars_[-1]] = theta  # elimination variable = last var (linform picks it)
        ELIM = not any(vv.startswith('l') for vv in vars_)
        # try conventions: x = raw, x = -raw/wp, x = raw/wp, x = -raw/dn
        def dphi_at(v):
            lam = []
            for i in range(n):
                val = conv(0)
                for mono, cf in dPhi_w[i].to_dict().items():
                    t = conv(int(cf))
                    for a, e in enumerate(mono):
                        if e: t = t * v[a]**e
                    val = val + t
                lam.append(val)
            return lam
        def check(vals):
            v = [conv(1)] + [vals[f'v{i}'] for i in range(1, n)]
            if ELIM:
                return True
            lam = [vals[f'l{i}'] for i in range(n)]
            # lambda_i == dPhi_w(v)?
            for i in range(n):
                val = conv(0)
                for mono, cf in dPhi_w[i].to_dict().items():
                    t = conv(int(cf))
                    for a, e in enumerate(mono):
                        if e: t = t * v[a]**e
                    val = val + t
                if val != lam[i]: return False
            return True
        found = None
        for label, tr in [('raw', lambda x: x), ('-raw/wp', lambda x: -x / wp), ('raw/wp', lambda x: x / wp), ('-raw/dn', lambda x: -x / dn), ('raw/dn', lambda x: x / dn)]:
            vals = {kk: (tr(x) if kk != vars_[-1] else x) for kk, x in raw.items()}
            if check(vals):
                found = (label, vals); break
        if found is None:
            results['spurious'] += 1; results.setdefault('undecoded', 0); results['undecoded'] += 1
            continue
        label, vals = found
        v = [conv(1)] + [vals[f'v{i}'] for i in range(1, n)]
        lam = dphi_at(v) if ELIM else [vals[f'l{i}'] for i in range(n)]
        # genuine reflection test: Phi_w(s z) = Phi_w(z) with s z = z - 2 ell(z)/ell(v) v, ell = sum lam_i z_i
        ellv = sum((lam[i] * v[i] for i in range(n)), conv(0))
        genuine = ellv != 0
        if genuine:
            # test on random points z over the field
            for trial in range(6):
                z = [conv(random.randrange(p)) if k == 1 else Fctx.random_element() for _ in range(n)]
                ellz = sum((lam[i] * z[i] for i in range(n)), conv(0))
                sz = [z[i] - 2 * ellz / ellv * v[i] for i in range(n)]
                def ev(pt):
                    val = conv(0)
                    for mono, cf in Phi_w.to_dict().items():
                        t = conv(int(cf))
                        for a, e in enumerate(mono):
                            if e: t = t * pt[a]**e
                        val = val + t
                    return val
                if ev(sz) != ev(z):
                    genuine = False; break
        if genuine:
            results['genuine'] += 1
            # eigenvector in original coordinates: z = G w
            vz = [sum((flint.nmod(G[a][b], p) * (v[b] if k == 1 else conv(0)) for b in range(n)), flint.nmod(0, p)) if k == 1 else None for a in range(n)]
            results['eigvecs'].append({'field_deg': k, 'v_orig': [str(x) for x in vz] if k == 1 else 'F_p^2', 'conv': label})
        else:
            results['spurious'] += 1
        if k == 1: results['rational'] += 1
        else: results['quadratic'] += 1
print(json.dumps({kk: vv for kk, vv in results.items() if kk != 'eigvecs'}))
for e in results['eigvecs']:
    print('  eigvec', e)
S['decode'] = {kk: vv for kk, vv in results.items()}
json.dump(S, open(f'{name}.json', 'w'))
