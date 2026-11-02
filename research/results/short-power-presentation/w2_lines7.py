import sys, json
sys.path.insert(0, '/tmp/hodge/results/short-power-presentation')
from w2lib import *
from msparse import parse_ms, evaluate
from flint import nmod_poly
F, A, q, rec = load()
import sys as _s
base = _s.argv[1] if len(_s.argv) > 1 else 'relax_d7_0_w2'
vars_, p, eqs = parse_ms(base + '.ms'); meta = json.load(open(base + '_meta.json')); labels = meta['labels']
pos = {lab: i for i, lab in enumerate(labels)}
Fp = reduce_mod(F, p); Fi = [pderiv(Fp, i, p) for i in range(NV)]
E6 = json.load(open('E6.json'))
def field_for(target):
    """solve target = sum_j w_j F_j (deg w = 2) mod p -> dict (j, m) -> coeff"""
    part = core.parity_of(target)
    gens = []; meta_ = []
    for j in range(NV):
        pj = 0 if j < 3 else 1
        for m in monomials(2, (part - pj) % 2):
            gens.append(shift(Fi[j], m, 1, p)); meta_.append((j, m))
    rows = {}
    for g in gens:
        for e in g: rows.setdefault(e, len(rows))
    for e in target: rows.setdefault(e, len(rows))
    R, C = len(rows), len(gens); flat = [0] * (R * (C + 1))
    for jj, g in enumerate(gens):
        for e, v in g.items(): flat[rows[e] * (C + 1) + jj] = v
    for e, v in target.items(): flat[rows[e] * (C + 1) + C] = v
    M = nmod_mat(R, C + 1, flat, p); Mr, rank = M.rref()
    piv = []; j = 0
    for i in range(rank):
        while int(Mr[i, j]) == 0: j += 1
        piv.append(j); j += 1
    if piv and piv[-1] == C: return None
    w = {}
    for i, j in enumerate(piv):
        v = int(Mr[i, C])
        if v: w[meta_[j]] = v
    return w
def vec_of(w, s=0):
    x = [0] * len(labels)
    if s: x[pos['K0']] = s % p
    for (j, m), v in w.items():
        lab = f'J{j}:{"".join(map(str, m))}'
        if lab in pos: x[pos[lab]] = (x[pos[lab]] + v) % p
        # gauge-removed coordinate: the gauge direction was projected out; the vector w may have a component there:
        # subtract the gauge multiple h*E/6 to land in the chart (h from the (5, h*Y2) coefficient)
    return x
def fix_gauge(w):
    # remove the (5, h Y2) components by subtracting (6*coef) * (h/6) E  == coef * h * E
    w = dict(w)
    for h in monomials(1, int(base.split('_')[2])):
        key = (5, mon_add(h, (0, 0, 0, 0, 0, 1)))
        c = w.get(key, 0)
        if c:
            for j in range(NV):
                ej = tuple(1 if t == j else 0 for t in range(NV))
                k2 = (j, mon_add(h, ej)); w[k2] = (w.get(k2, 0) - c) % p
    return {k: v for k, v in w.items() if v}
Lp = reduce_mod(core.parse('X0+2*X1+3*X2'), p); Mp = reduce_mod(core.parse('Y0+3*Y1+5*Y2'), p)
cands = {}
E6L = reduce_mod(core.parse(E6['L']), p); E6M = reduce_mod(core.parse(E6['M']), p)
for name, tgt in [('L*E6', mul(Lp, E6L, p)), ('M*E6M', mul(Mp, E6M, p)), ('L*E6M', mul(Lp, E6M, p)), ('M*E6', mul(Mp, E6L, p))]:
    if core.parity_of(tgt) != meta.get('parity', int(base.split('_')[2])):
        continue
    w = field_for(tgt)
    cands['w(F)=' + name] = fix_gauge(w) if w is not None else None
for name, w in cands.items():
    if w is None: print(name, ': target not in J_7 (no field)', flush=True); continue
    x = vec_of(w)
    nz = sum(1 for eq in eqs if evaluate(eq, x, p) != 0)
    print(f'{name}: nonzero equations at (s=0, w) = {nz} / {len(eqs)}', flush=True)
    # line (s, t): x(s,t) = s*e_K0 + t*w : each equation is a binary quadratic form in (s,t); projective common roots
    xs = vec_of({}, 1); xt = vec_of(w)
    def q_coeffs(eq):
        # evaluate at (1,0),(0,1),(1,1) to get a s^2 + b t^2 + c s t
        a_ = evaluate(eq, xs, p); b_ = evaluate(eq, xt, p); c_ = (evaluate(eq, [(u + v) % p for u, v in zip(xs, xt)], p) - a_ - b_) % p
        return a_, b_, c_
    g = None
    for eq in eqs:
        a_, b_, c_ = q_coeffs(eq)
        # dehomogenize s=1: b t^2 + c t + a
        poly = nmod_poly([a_, c_, b_], p)
        g = poly if g is None else g.gcd(poly)
        if g.degree() == 0: break
    print(f'   line span(D7, {name}): gcd degree over t (s=1) = {g.degree()}; roots mod p: {g.roots() if g.degree() > 0 else []}', flush=True)
