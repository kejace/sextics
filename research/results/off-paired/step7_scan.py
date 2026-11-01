"""Unified scan of structured eigenvector denominators (degrees 7..10): E1 test, fixed-D consistency for A, and
served-target dimension (nontrivial G-invariant (3,1) classes carried by V(D))."""
import oplib, time, sys, itertools, random, json
from hodge import core
from hodge.core import mul, pderiv, add, sub, monomials, shift
from flint import nmod_mat
F, A, q = oplib.load()
p = oplib.P0
Ap = core.reduce_mod(A, p)
Fi = [pderiv(F, i) for i in range(6)]
Fij = [[pderiv(Fi[i], j) for j in range(6)] for i in range(6)]
D0 = {}
for e, c in Fi[0].items():
    e2 = list(e); e2[0] -= 1; D0[tuple(e2)] = c / 6
V = core.VARS
def mstr(e): return '*'.join(f'{V[i]}**{k}' if k > 1 else V[i] for i, k in enumerate(e) if k) or '1'
rng = random.Random(11)
cands = []
# --- degree 8: Hessian 2x2 minors (all), D0^2
for (i, j) in itertools.combinations(range(6), 2):
    for (k, l) in itertools.combinations_with_replacement(range(6), 2):
        Dm = sub(mul(Fij[i][k], Fij[j][l]), mul(Fij[i][l], Fij[j][k]))
        if Dm: cands.append((f'Hess[{i}{j}|{k}{l}]', 8, 'hessian-2x2-minor', Dm))
cands.append(('D0**2', 8, 'quartic-colon-square', mul(D0, D0)))
# --- degree 8: sampled partial*cubic, D0*quartic, E_I(F)*quadratic
ms3 = monomials(3); rng.shuffle(ms3)
for m in ms3[:12]:
    for i in range(6): cands.append((f'F_{V[i]}*{mstr(m)}', 8, 'partial*monomial(sample)', shift(Fi[i], m)))
ms4 = monomials(4); rng.shuffle(ms4)
for m in ms4[:30]: cands.append((f'D0*{mstr(m)}', 8, 'quartic-colon*monomial(sample)', shift(D0, m)))
subsets = [I for r in (1, 2, 3) for I in itertools.combinations(range(6), r)]
ms2 = monomials(2)
for I in subsets:
    EI = {}
    for i in I: EI = add(EI, shift(Fi[i], tuple(1 if j == i else 0 for j in range(6))))
    for m in rng.sample(ms2, 2):
        cands.append((f'E_{"".join(V[i] for i in I)}(F)*{mstr(m)}', 8, 'subset-euler*monomial(sample)', shift(EI, m)))
# --- degree 9
for i in range(6): cands.append((f'D0*F_{V[i]}', 9, 'quartic-colon*partial', mul(D0, Fi[i])))
ms4b = monomials(4); rng.shuffle(ms4b)
for m in ms4b[:10]:
    for i in range(6): cands.append((f'F_{V[i]}*{mstr(m)}', 9, 'partial*monomial(sample)', shift(Fi[i], m)))
# --- degree 10
for (i, j) in itertools.combinations_with_replacement(range(6), 2):
    cands.append((f'F_{V[i]}*F_{V[j]}', 10, 'partial*partial', mul(Fi[i], Fi[j])))
for m in monomials(2): cands.append((f'D0**2*{mstr(m)}', 10, 'quartic-colon-square*monomial', shift(mul(D0, D0), m)))
ms5 = monomials(5); rng.shuffle(ms5)
for m in ms5[:6]:
    for i in range(6): cands.append((f'F_{V[i]}*{mstr(m)}', 10, 'partial*monomial(sample)', shift(Fi[i], m)))
print('candidates', len(cands), file=sys.stderr)
seen = set(); out = open('served_scan.jsonl', 'a')
Amons_idx = None
for name, d, fam, D in cands:
    Dp = core.reduce_mod(D, p)
    key = tuple(sorted(Dp.items()))
    if key in seen or not Dp: continue
    seen.add(key)
    cD = oplib.char_of(Dp)
    if cD is None:
        print('skip non-eigen', name, file=sys.stderr); continue
    t0 = time.time()
    ok, info = oplib.e1_reduced(F, A, q, D, p)
    r = oplib.served_targets(F, q, D, p)
    # is A served?  A in span(served) mod p  (J^odd_6 is always served, so this is exactly fixed-D consistency)
    Amons = oplib.monos_char(6, (3, 0)); idx = {m: k for k, m in enumerate(Amons)}
    def vec(P):
        v = [0]*len(Amons)
        for e, val in P.items(): v[idx[e]] = val % p
        return v
    rows = [vec(b) for b in r['basis']]
    rk = nmod_mat(len(rows), len(Amons), [x for rr in rows for x in rr], p).rank() if rows else 0
    rows2 = rows + [vec(Ap)]
    rk2 = nmod_mat(len(rows2), len(Amons), [x for rr in rows2 for x in rr], p).rank()
    consistent = (rk2 == rk)
    dt = time.time() - t0
    rec = {'name': name, 'd': d, 'family': fam, 'char': cD, 'e1': ok, 'consistent': consistent, 'dim_served': r['dim_served'],
           'dim_nontrivial': r['dim_nontrivial'], 'rows': r['rows'], 'cols': r['cols'], 'seconds': round(dt, 1)}
    out.write(json.dumps(rec) + '\n'); out.flush()
    oplib.log(core.to_str(D) if len(D) < 30 else name, d, fam, ok, consistent, r['cols'] - r['dim_served'], r['dim_served'], dt,
              f"{name}: char={cD}; served (3,0)-targets dim {r['dim_served']} (nontrivial {r['dim_nontrivial']}); rows={r['rows']} cols={r['cols']}")
    flag = ' !!! NONTRIVIAL' if r['dim_nontrivial'] else ''
    if consistent:
        flag += ' !!!!! A SERVED -> CERTIFICATE'
        st, data = oplib.lift_solution(F, A, q, D)
        if st == 'solved':
            json.dump({'D': core.to_str(D), 'N': [core.to_str(n) for n in data['N']], 'Q': core.to_str(data['Q'])}, open('witness-candidate.json', 'w'))
    print(f"{name:34s} d={d} char={cD} e1={ok} consistent={consistent} served={r['dim_served']} nontriv={r['dim_nontrivial']} {dt:.1f}s{flag}", file=sys.stderr)
    if r['dim_nontrivial']:
        for b in r['basis']:
            lb = oplib.lift_poly_small(b, p)
            print('    served:', core.to_str(lb) if lb else '(no small lift)', file=sys.stderr)
print('done', file=sys.stderr)
