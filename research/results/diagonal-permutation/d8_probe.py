"""Probe the coset D8 + J_8^+ (mod p): symmetry, sparse representatives, membership of natural monomials, vanishing loci."""
import sys, json, time, itertools; sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core
from hodge.core import monomials, shift, pderiv, parse, to_str, mul, add, sub, scale
from flint import nmod_mat
F, A, q, rec = core.load_instance('diagonal-permutation')
p = core.PRIMES[0]
Fp = core.reduce_mod(F, p)
Fi = [pderiv(Fp, i, p) for i in range(6)]
d = 8
Dmons = monomials(d, 0); idx = {m: j for j, m in enumerate(Dmons)}; CA = len(Dmons)
def vec(P):
    Pp = core.reduce_mod(P, p); v = [0] * CA
    for e, c in Pp.items(): v[idx[e]] = c
    return v
jd = []
for i in range(6):
    par_needed = (0 - (0 if i < 3 else 1)) % 2
    for m in monomials(d - 5, par_needed): jd.append(shift(Fi[i], m, 1, p))
Jrows = [vec(P) for P in jd]
J = nmod_mat(len(Jrows), CA, sum(Jrows, []), p); rkJ = J.rank()
D8 = parse(json.load(open('D8_e1kernel.json'))['D']); v8 = vec(D8)
def rank_with(rows):
    return nmod_mat(len(Jrows) + len(rows), CA, sum(Jrows, []) + sum(rows, []), p).rank()
assert rank_with([v8]) == rkJ + 1
print('rkJ', rkJ, 'D8 not in J: ok')
# symmetry
def perm(P, sig):
    out = {}
    for e, v in P.items():
        ne = [0] * 6
        for i in range(3): ne[sig[i]] = e[i]; ne[3 + sig[i]] = e[3 + i]
        out[tuple(ne)] = (out.get(tuple(ne), 0) + v)
    return {k: v for k, v in out.items() if v}
D01 = perm(D8, (1, 0, 2)); Dc = perm(D8, (1, 2, 0))
for nm, P in [('(01)D8', D01), ('cycD8', Dc)]:
    r = rank_with([v8, vec(P)])
    print(nm, 'in span(J, D8)?', r == rkJ + 1)
    if r == rkJ + 1:
        for s in (1, -1):
            print('   ', nm, '=', s, '* D8 mod J ?', rank_with([vec(sub(P, scale(D8, s)))]) == rkJ)
# natural monomials/forms in the coset?
cands = {}
X = [parse(v) for v in ['X0','X1','X2','Y0','Y1','Y2']]
T = [mul(core.power(X[i],3), core.power(X[3+i],3)) for i in range(3)]
cands['X0^2X1^2Y0^2Y1^2'] = parse('X0**2*X1**2*Y0**2*Y1**2')
cands['X0^2Y0^2 X1^2Y1^2 + X0^2Y0^2 X2^2Y2^2 + X1^2Y1^2 X2^2Y2^2'] = parse('X0**2*Y0**2*X1**2*Y1**2 + X0**2*Y0**2*X2**2*Y2**2 + X1**2*Y1**2*X2**2*Y2**2')
cands['(X0Y0)^2 (X1Y1)^2 - ... anti'] = parse('X0**2*Y0**2*X2**2*Y2**2 - X1**2*Y1**2*X2**2*Y2**2')
cands['T0*X0^2 - T1*X1^2'] = parse('X0**5*Y0**3 - X1**5*Y1**3')
cands['X0^4Y0^4 + X1^4Y1^4'] = parse('X0**4*Y0**4 + X1**4*Y1**4'); cands['X0^4Y0^4 - X1^4Y1^4'] = parse('X0**4*Y0**4 - X1**4*Y1**4')
cands['X2^4 Y2^4'] = parse('X2**4*Y2**4'); cands['X2^8'] = parse('X2**8'); cands['Y2^8'] = parse('Y2**8')
cands['A^2 / (X0^3Y0^3 - X1^3Y1^3)... A*(T0+T1)?'] = mul(A, add(T[0], T[1]))
cands['A*T2'] = mul(A, T[2]); cands['A*(X0^2-X1^2)?'] = None
cands['(T0-T1)(T0+T1-2T2)'] = mul(A, parse('X0**3*Y0**3 + X1**3*Y1**3 - 2*X2**3*Y2**3'))
cands['T0^2 - T1^2 + ... (X0^3Y0^3)^2 - (X1^3Y1^3)^2'] = parse('X0**6*Y0**6 - X1**6*Y1**6')
cands['F_2-part squared? (sumX^2)(sumY^2)... '] = None
for nm, P in cands.items():
    if P is None or core.degree(P) != 8 or core.parity_of(P) != 0: continue
    inJ = rank_with([vec(P)]) == rkJ
    incoset = rank_with([v8, vec(P)]) == rkJ + 1
    print(f'{nm:70s} inJ={inJ} in span(J,D8)={incoset}')
# vanishing loci of the coset: exists c != 0, j in J with c*D8 + j in I^k  <=>  D8 in I^k + J
def ideal_power_rows(lin_forms, k):
    """basis (mod p) of degree-8 Y-even part of I^k, I = (lin_forms): all products of k forms times monomials of degree 8-k"""
    rows = []
    for combo in itertools.combinations_with_replacement(range(len(lin_forms)), k):
        P = {(0,)*6: 1}
        for c in combo: P = mul(P, lin_forms[c], p)
        for m in monomials(8 - k, None):
            Q = shift(P, m, 1, p)
            if core.parity_of(Q) == 0: rows.append(vec(Q))
    return rows
H = core.Fraction(1, 2) if hasattr(core, 'Fraction') else None
from fractions import Fraction as Fr
lin = {'X0': X[0], 'X1': X[1], 'X2': X[2], 'Y0': X[3], 'Y1': X[4], 'Y2': X[5],
       'p': add(X[0], X[1]), 'r': sub(X[0], X[1]), 'q': add(X[3], X[4]), 's': sub(X[3], X[4]),
       'r1': sub(X[0], X[1]), 'r2': sub(X[1], X[2]), 's1': sub(X[3], X[4]), 's2': sub(X[4], X[5])}
tests = [(['X0','X1','X2'], 3), (['X0','X1','X2'], 4), (['Y0','Y1','Y2'], 3), (['Y0','Y1','Y2'], 4),
         (['X0','Y0'], 2), (['X0','Y0'], 3), (['X2','Y2'], 2), (['X2','Y2'], 3), (['X0','X1'], 2), (['Y0','Y1'], 2),
         (['r','s'], 2), (['r','s'], 3), (['p','q'], 2), (['p','q'], 3),
         (['X0','X1','X2','Y0'], 4), (['X0','X1','Y0','Y1'], 4), (['X0','Y0','X1'], 3), (['X0','Y0','Y1'], 3),
         (['r','s','X2'], 3), (['r','s','Y2'], 3), (['p','q','X2','Y2'], 4), (['r','s','X2','Y2'], 4),
         (['r1','r2','s1','s2'], 4), (['r1','r2','s1','s2'], 3), (['r1','r2'], 2), (['s1','s2'], 2), (['X0','X1','X2','Y0','Y1','Y2'], 5)]
for names, k in tests:
    t0 = time.time()
    rows = ideal_power_rows([lin[n] for n in names], k)
    base = nmod_mat(len(Jrows) + len(rows), CA, sum(Jrows, []) + sum(rows, []), p).rank()
    withD = nmod_mat(len(Jrows) + len(rows) + 1, CA, sum(Jrows, []) + sum(rows, []) + v8, p).rank()
    print(f'I=({",".join(names)})^{k}: dim(I^k+J)_8={base} (J={rkJ}, I^k rows={len(rows)}), D8 in I^k+J: {withD == base}  [{time.time()-t0:.1f}s]', flush=True)
