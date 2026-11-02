import sys, json, time; sys.path.insert(0,'/tmp/hodge/results/eight-plane-configuration/verify')
from v_lib import *
import sympy
from sympy import Matrix, I, Rational, factorint
from flint import fmpq_mat, fmpq, fmpz_mat
F, A, q, fields = load_instance()
res = {}
# 1. support
allY_even = set(monomials(6, 0))
res['support_is_all_236_Yeven'] = (set(F) == allY_even) and all(c != 0 for c in F.values())
# 2. denominators
dens = sorted({c.denominator for c in F.values()})
primes = set()
for dn in dens: primes |= set(factorint(dn).keys())
res['denominator_primes'] = sorted(primes)
res['all_5_integral'] = all(c.denominator % 5 != 0 for c in F.values())
Fmod5 = reduce_mod(F, 5)
fermat5 = {tuple(6 if j==i else 0 for j in range(6)): 1 for i in range(6)}
res['F_mod5_is_Fermat'] = (Fmod5 == fermat5)
# G = (F - Fermat)/5 is 5-integral?
G = scale(add(F, scale(fermat5, -1)), Fraction(1,5))
res['G_5_integral'] = all(c.denominator % 5 != 0 for c in G.values())
# 3. planes and F in (Q0,Q1,Q2)
X = sympy.symbols('X0 X1 X2 Y0 Y1 Y2')
Fs = sympy.sympify(fields['F'])
planes_ok = True
for e0 in (1,-1):
  for e1 in (1,-1):
    for e2 in (1,-1):
      sub = {X[3]: e0*I*X[0], X[4]: e1*I*X[1], X[5]: e2*I*X[2]}
      v = sympy.expand(Fs.subs(sub, simultaneous=True))
      if v != 0: planes_ok = False
res['eight_planes_in_X'] = planes_ok
K = json.load(open('/tmp/hodge/results/eight-plane-configuration/K_decomp.json'))
Qs = [parse_poly('X0**2+Y0**2'), parse_poly('X1**2+Y1**2'), parse_poly('X2**2+Y2**2')]
Ks = [parse_poly(s) for s in K['K']]
S = {}
for Qi, Ki in zip(Qs, Ks): S = add(S, mul(Qi, Ki))
res['F_eq_sum_QiKi_exact'] = (S == F)
Hs = [parse_poly(s) for s in K['H']]
ok = True
for i, (Ki, Hi, pi) in enumerate(zip(Ks, Hs, (7,11,13))):
    diag = parse_poly(f'X{i}**4-X{i}**2*Y{i}**2+Y{i}**4')
    ok &= (Ki == add(diag, scale(Hi, Fraction(5, pi))))
res['K_i_eq_diag_plus_5Hi_over_pi'] = ok
res['H_sizes'] = [len(h) for h in Hs]
res['H_all_Yeven'] = [parity_of(h) for h in Hs]
# 4. SNF of exponent-difference lattice
exps = sorted(F)
e0 = exps[0]
M = Matrix([[a-b for a, b in zip(e, e0)] for e in exps[1:]])
from sympy.matrices.normalforms import smith_normal_form
snf = smith_normal_form(M, domain=sympy.ZZ)
inv = [snf[i,i] for i in range(min(snf.shape)) if snf[i,i] != 0]
res['snf_invariants'] = [int(x) for x in inv]; res['lattice_rank'] = len(inv)
# diagonal symmetries: characters t with t^v=1 for all v in lattice; torsion part Z/2 -> candidates: sign vectors
signs = []
import itertools
for s in itertools.product((1,-1), repeat=6):
    if s[0] == 1 and all((__import__('math').prod(si**ei for si, ei in zip(s, e)) == 1) for e in exps): signs.append(s)
res['sign_symmetries_mod_scalar'] = signs
# 5. A not in J_6 exactly over Q
Fi = [pderiv(F, i) for i in range(6)]
gens = []
for i in range(6):
    for j in range(6):
        m = tuple(1 if k==j else 0 for k in range(6))
        g = shift(Fi[i], m)
        if parity_of(g) == 1: gens.append(g)
rows = sorted({e for g in gens for e in g} | set(A))
def fq(P, cols):
    M = fmpq_mat(len(rows), len(cols))
    for j, g in enumerate(cols):
        for e, c in g.items():
            M[rows.index(e), j] = fmpq(c.numerator, c.denominator)
    return M
MJ = fq(None, gens); MJA = fq(None, gens + [A])
res['J6odd_rank_Q'] = MJ.rank(); res['J6odd_plus_A_rank_Q'] = MJA.rank()
res['A_notin_J_exact_Q'] = MJA.rank() == MJ.rank() + 1
# 6. Hessian symmetrizer dimension exactly over Q: M (6x6) with M*Hess symmetric
Hess = [[pderiv(Fi[i], j) for j in range(6)] for i in range(6)]
# unknowns M_ab (36). (M H)_{ij} = sum_a M_ia H_aj ; condition (MH)_{ij} - (MH)_{ji} = 0 for i<j
rowsH = sorted(set(monomials(4)))
eqs = []
for i in range(6):
    for j in range(i+1, 6):
        # coefficient of unknown M_ia is H_aj ; of M_ja is -H_ai
        for e in rowsH:
            row = [Fraction(0)]*36
            for a in range(6):
                row[i*6+a] += Hess[a][j].get(e, 0)
                row[j*6+a] -= Hess[a][i].get(e, 0)
            if any(row): eqs.append(row)
ME = fmpq_mat(len(eqs), 36)
for r, row in enumerate(eqs):
    for c, v in enumerate(row):
        if v: ME[r, c] = fmpq(v.numerator, v.denominator)
res['hessian_symmetrizer_dim_Q'] = 36 - ME.rank()
# 7. independence of quartics F_ij
def rankQ(polys):
    rws = sorted({e for P in polys for e in P})
    M = fmpq_mat(len(rws), len(polys))
    for j, P in enumerate(polys):
        for e, c in P.items(): M[rws.index(e), j] = fmpq(c.numerator, c.denominator)
    return M.rank()
XX = [Hess[i][j] for i in range(3) for j in range(i, 3)]
YY = [Hess[i][j] for i in range(3, 6) for j in range(i, 6)]
XY = [Hess[i][j] for i in range(3) for j in range(3, 6)]
res['rank_FXX_6'] = rankQ(XX); res['rank_FYY_6'] = rankQ(YY); res['rank_even12'] = rankQ(XX+YY); res['rank_FXY_9'] = rankQ(XY); res['rank_all21'] = rankQ(XX+YY+XY)
# 8. plane rigidity over Q(i): map (L_0,L_1,L_2) in (S_1(X))^3 -> sum_i F_{Y_i}|_Pi * L_i, rank must be 9
rig = {}
Xs = X[:3]
for eps in itertools.product((1,-1), repeat=3):
    sub = {X[3]: eps[0]*I*X[0], X[4]: eps[1]*I*X[1], X[5]: eps[2]*I*X[2]}
    FY = [sympy.expand(sympy.diff(Fs, X[3+i]).subs(sub, simultaneous=True)) for i in range(3)]
    cols = []
    for i in range(3):
        for xj in Xs:
            cols.append(sympy.Poly(sympy.expand(FY[i]*xj), *Xs))
    mons = sorted({m for c in cols for m in c.monoms()})
    Mr = Matrix([[c.coeff_monomial(sympy.Mul(*[x**k for x, k in zip(Xs, m)])) for c in cols] for m in mons])
    rig[str(eps)] = Mr.rank(simplify=True)
res['plane_rigidity_rank9_over_Qi'] = rig
print(json.dumps(res, indent=1, default=str))
for k, v in res.items():
    log({'check': 'structural:'+k, 'value': v if not isinstance(v, (list, dict)) or len(str(v)) < 300 else str(v)[:300], 'method': 'exact over Q / Q(i), own code (v_lib.py + sympy + flint fmpq)'})
