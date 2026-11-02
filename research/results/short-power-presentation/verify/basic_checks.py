"""Cheap exact checks of structural claims (independent of hodge.*)."""
import sys, json, itertools, subprocess
sys.path.insert(0, '/tmp/hodge/results/short-power-presentation/verify')
from vlib import *
from flint import nmod_mat, fmpq_mat, fmpq
X0, X1, X2, Y0, Y1, Y2 = VARS
Fs, As, q = sympy_F()
L = X0 + 2*X1 + 3*X2; M = Y0 + 3*Y1 + 5*Y2
c = sp.Rational(5, 7**7)
res = {}
# C1: presentation
pres = sum(v**6 for v in VARS) + c*((L+M)**6 + (L-M)**6)
res['C1_presentation'] = sp.expand(Fs - pres) == 0
print('C1 F == eX+eY+c((L+M)^6+(L-M)^6):', res['C1_presentation'])
print('   q =', q, ' A =', As, ' #terms F =', len(sp.Poly(Fs, *VARS).terms()))
# C: coefficient sanity
Pf = sp.Poly(Fs, *VARS)
print('   X0^6 coeff:', Pf.coeff_monomial(X0**6), '=', 1 + 2*c, ';  Y2^6:', Pf.coeff_monomial(Y2**6), '=', 1 + 2*c*5**6)
# quartic power rank: X_i^4, Y_i^4, (L+M)^4, (L-M)^4 independent over Q?
forms = [v for v in VARS] + [L + M, L - M]
quartics = [sp.Poly(f**4, *VARS) for f in forms]
mons4 = monomials(4)
Mq = sp.Matrix([[Q.coeff_monomial(sp.prod([v**e for v, e in zip(VARS, m)])) for m in mons4] for Q in quartics])
res['quartic_power_rank'] = Mq.rank()
print('quartic power rank (exact over Q):', res['quartic_power_rank'])
sextics = [sp.Poly(f**6, *VARS) for f in forms]
mons6 = monomials(6)
Ms = sp.Matrix([[Q.coeff_monomial(sp.prod([v**e for v, e in zip(VARS, m)])) for m in mons6] for Q in sextics])
res['sextic_power_rank'] = Ms.rank()
print('sextic power rank (l_k^6 independent?):', res['sextic_power_rank'])
# Hessian symmetrizer: {C in M_6(Q): Hess(F) C symmetric} -- exact over Q via sympy linear algebra on coefficients
Hs = sp.hessian(Fs, VARS)
Cs = sp.Matrix(6, 6, sp.symbols('c0:36'))
Sym = Hs * Cs - (Hs * Cs).T
eqs = []
for i in range(6):
    for j in range(i + 1, 6):
        P = sp.Poly(sp.expand(Sym[i, j]), *VARS)
        eqs.extend(P.coeffs())
Aeq, _ = sp.linear_eq_to_matrix(eqs, list(Cs))
res['hessian_symmetrizer_dim_Q'] = 36 - Aeq.rank()
print('Hessian symmetrizer dim (exact over Q):', res['hessian_symmetrizer_dim_Q'])
# diagonal group: exponent-difference lattice of supp F
exps = [tuple(m) for m, _ in Pf.terms()]
e0 = exps[0]
diffs = sp.Matrix([[a - b for a, b in zip(e, e0)] for e in exps[1:]])
# the group is dual to Z^6/L; projectively: L contains (1,...,1)? no: differences have sum 0.  Compute Smith normal form of diffs
from sympy.matrices.normalforms import smith_normal_form
snf = smith_normal_form(diffs, domain=sp.ZZ)
diag = [snf[i, i] for i in range(min(snf.shape))]
print('SNF diag of exponent-difference lattice:', diag, ' (rank', diffs.rank(), ')')
# lattice index in {sum=0}: product of nonzero invariant factors
res['diag_group_order_projective'] = int(sp.prod([d for d in diag if d != 0]))
print('=> diagonal projective group order (index of L in sum-zero lattice):', res['diag_group_order_projective'])
# E_I^2(F) = 6 E_I(F)?  for proper subsets I
bad = []
for r in range(1, 6):
    for I in itertools.combinations(range(6), r):
        EI = lambda G: sum(VARS[i]*sp.diff(G, VARS[i]) for i in I)
        if sp.expand(EI(EI(Fs)) - 6*EI(Fs)) == 0: bad.append(I)
res['E_I_squared_subsets_passing'] = bad
print('proper subsets I with E_I^2 F = 6 E_I F:', bad)
# apolar quadrics (q(d) F = 0) and their base locus, exact via Singular over Q
D = [sp.diff]  
quad_mons = monomials(2)
rows = []
for m in quad_mons:
    G = Fs
    for v, e in zip(VARS, m):
        G = sp.diff(G, v, e) if e else G
    rows.append(sp.Poly(sp.expand(G), *VARS))
mons4 = monomials(4)
Mcat = sp.Matrix([[R.coeff_monomial(sp.prod([v**e for v, e in zip(VARS, m)])) for m in mons4] for R in rows])
ns = Mcat.T.nullspace()  # vectors a with sum a_m d^m F = 0
res['apolar_quadrics'] = len(ns)
print('number of apolar quadrics (exact):', len(ns), ' catalecticant Cat(2,4) rank:', Mcat.rank())
# middle catalecticant rank
mons3 = monomials(3)
rows3 = []
for m in mons3:
    G = Fs
    for v, e in zip(VARS, m):
        G = sp.diff(G, v, e) if e else G
    rows3.append(sp.Poly(sp.expand(G), *VARS))
Mcat3 = sp.Matrix([[R.coeff_monomial(sp.prod([v**e for v, e in zip(VARS, m)])) for m in mons3] for R in rows3])
res['middle_catalecticant_rank_Q'] = Mcat3.rank()
print('middle catalecticant rank (exact):', res['middle_catalecticant_rank_Q'])
# base locus of apolar quadrics: Singular over Q
u = sp.symbols('u0:6')
quads = []
for vec in ns:
    Qd = sum(vec[i] * sp.prod([u[k]**e for k, e in zip(range(6), m)]) for i, m in enumerate(quad_mons))
    quads.append(sp.expand(Qd))
sing = 'ring R = 0, (u0,u1,u2,u3,u4,u5), dp;\nideal I = ' + ',\n'.join(str(qd).replace('**', '^') for qd in quads) + ';\n'
sing += 'ideal G = std(I); print("dim = " + string(dim(G))); print("degree = " + string(mult(G)));\n'
sing += 'LIB "primdec.lib"; list Lp = minAssGTZ(G); print("ncomp = " + string(size(Lp))); Lp;\n'
# check the 8 expected points lie on all quadrics
sing += 'quit;\n'
open('apolar_baselocus.sing', 'w').write(sing)
out = subprocess.run(['Singular', '-q', 'apolar_baselocus.sing'], capture_output=True, text=True, timeout=600).stdout
print('--- Singular apolar base locus ---'); print(out[:3000])
res['apolar_baselocus_singular'] = out[:3000]
pts = [[1 if k == i else 0 for k in range(6)] for i in range(6)] + [[1, 2, 3, 1, 3, 5], [1, 2, 3, -1, -3, -5]]
onall = all(all(qd.subs(dict(zip(u, pt))) == 0 for qd in quads) for pt in pts)
res['8_points_on_all_apolar_quadrics'] = onall
print('8 power points lie on all apolar quadrics:', onall)
json.dump(res, open('basic_checks.json', 'w'), indent=1, default=str)
for k, v in res.items():
    if k != 'apolar_baselocus_singular': logjson(claim='basic:' + k, verdict='DATA', value=v if not isinstance(v, list) else str(v))
