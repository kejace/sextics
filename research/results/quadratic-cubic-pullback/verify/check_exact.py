"""Cheap exact checks (own code, sympy exact / mod p=1000000007). Logs each check."""
import sys, json, itertools, time
sys.path.insert(0, '/tmp/hodge/results/quadratic-cubic-pullback/verify')
import vpoly as V
from vlog import log
from fractions import Fraction
import sympy
p = 1000000007
F, A, q = V.load_instance_F_A()
chk = json.load(open('/home/kejace/devel/sextics/data/quadratic-cubic-pullback/checks.json'))
Q = [V.parse(s) for s in chk['quadrics']]

# C1: F = sum Q_i^3 exactly; parity; term count; pure sixth powers
S = {}
for Qi in Q: S = V.add(S, V.power(Qi, 3))
diff = V.add(S, V.scale(F, -1))
pure = all(F.get(tuple(6 if j == i else 0 for j in range(6))) == 1 for i in range(6))
log('C1 F = sum_i Q_i^3 exactly; F Y-even, 50 terms, six pure sixth powers with coefficient 1',
    'CONFIRMED' if (not diff and V.parity(F) == 0 and len(F) == 50 and pure) else 'REFUTED',
    f'F - sum Q^3 has {len(diff)} terms; parity {V.parity(F)}; {len(F)} terms; pure sixth powers coeff 1: {pure}')

# C3: mod-5 reduction: denominators coprime to 5, F mod 5 = Fermat, partials mod 5 = Z_i^5
dens_ok = all(c.denominator % 5 != 0 for c in F.values())
F5 = V.modp(F, 5)
ferm5 = {tuple(6 if j == i else 0 for j in range(6)): 1 for i in range(6)}
log('C3 good reduction mod 5: F == Fermat mod 5 (hence X smooth over Q by openness of smoothness)',
    'CONFIRMED' if dens_ok and F5 == ferm5 else 'REFUTED', f'denominators coprime to 5: {dens_ok}; F mod 5 == Fermat: {F5 == ferm5}')

# C4: projective diagonal symmetry group = Z/2 via SNF of exponent differences inside the sum-zero lattice
sup = sorted(F)
e0 = sup[0]
diffs = [tuple(e[i]-e0[i] for i in range(6)) for e in sup[1:]]
M = sympy.Matrix(diffs)
from sympy.matrices.normalforms import smith_normal_form
Sm = smith_normal_form(M, domain=sympy.ZZ)
dg = [abs(Sm[i, i]) for i in range(min(Sm.shape))]
# the quotient Z^6/L(F) = Z (scalars) + torsion; projective group = torsion part = product of nontrivial invariant factors
tors = [x for x in dg if x not in (0, 1)]
log('C4 projective diagonal symmetry group of F is exactly Z/2 = {1,tau} (SNF of exponent-difference lattice)',
    'CONFIRMED' if tors == [2] and dg.count(0) == 1 else 'REFUTED', f'SNF diag = {dg}; torsion = {tors}')

# C5: monomial symmetries (permutation * sign): only identity and tau
syms = []
for perm in itertools.permutations(range(6)):
    for signs in itertools.product((1, -1), repeat=6):
        G = {}
        ok = True
        for e, c in F.items():
            e2 = [0]*6
            s = 1
            for i in range(6):
                e2[perm[i]] = e[i]
                if signs[i] == -1 and e[i] % 2: s = -s
            e2 = tuple(e2)
            if F.get(e2) != c*s: ok = False; break
        if ok: syms.append((perm, signs))
ident = tuple(range(6))
expected = {(ident, (1,)*6), (ident, (1, 1, 1, -1, -1, -1)), (ident, (-1, -1, -1, 1, 1, 1)), (ident, (-1,)*6)}
got = set(syms)
log('C5 monomial symmetries (perm x sign, 46080 maps) of F: only 1 and tau (projectively)',
    'CONFIRMED' if got == expected else 'REFUTED', f'{len(syms)} affine maps found = {sorted(got)}; projectively {{1,tau}}: {got == expected}')

# C6: linear symmetrizer: dim {(M, lam): sum_ij M_ij Z_j d_i F = lam F} == 1 (Euler only), mod p
Fp = V.modp(F, p)
Fi = [V.pderiv(Fp, i, p) for i in range(6)]
cols = []
for i in range(6):
    for j in range(6):
        m = [0]*6; m[j] = 1
        cols.append(V.shift(Fi[i], tuple(m), 1, p))
cols.append(V.scale(Fp, p-1, p))
mons6 = V.monomials(6)
r = V.rank_mod_p(cols, mons6, p)
log('C6 Lie algebra of linear symmetries of F (up to scalar) is 1-dim (Euler only): hessian_symmetrizer_dimension = 1',
    'CONFIRMED' if 37 - r == 1 else 'REFUTED', f'37 unknowns, rank {r}, nullity {37-r} (mod {p})')

# C7: Jacobian relations in degree 6: 36 products Z_j F_i independent (exact), parity split 18/18; R_6 = 426 = 218+208
Fi_ex = [V.pderiv(F, i) for i in range(6)]
prods = []
par = []
for i in range(6):
    for j in range(6):
        m = [0]*6; m[j] = 1
        P = V.shift(Fi_ex[i], tuple(m))
        prods.append(P); par.append(V.parity(P))
odd = [P for P, s in zip(prods, par) if s == 1]
even = [P for P, s in zip(prods, par) if s == 0]
mons6o = V.monomials(6, 1); mons6e = V.monomials(6, 0)
ro = V.rank_exact(odd, mons6o); re_ = V.rank_exact(even, mons6e)
log('C7 h^{3,1}(X)^+ = 226 - 18 = 208, h^{3,1}(X)^- = 236 - 18 = 218, h^{3,1} = 426 (18+18 independent degree-6 Jacobian relations, exact over Q)',
    'CONFIRMED' if (len(odd), ro, len(even), re_, len(mons6o), len(mons6e)) == (18, 18, 18, 18, 226, 236) else 'REFUTED',
    f'odd relations {len(odd)} rank {ro}; even {len(even)} rank {re_}; S_6 odd/even {len(mons6o)}/{len(mons6e)}')

# C8: delta = det(dQ_i/dZ_j) exact; 64 terms, Y-odd; delta mod 5 = -X0X1X2Y0Y1Y2; A not in J_6; A, delta independent mod J_6 (exact)
X = sympy.symbols('X0 X1 X2 Y0 Y1 Y2')
def topoly(P): return sum(sympy.Rational(c.numerator, c.denominator)*sympy.prod([X[i]**e[i] for i in range(6)]) for e, c in P.items())
Qs = [topoly(Qi) for Qi in Q]
Jm = sympy.Matrix(6, 6, lambda i, j: sympy.diff(Qs[i], X[j]))
delta_s = sympy.expand(Jm.det())
delta = V.parse(str(delta_s))
d5 = V.modp(delta, 5)
rA = V.rank_exact(odd + [A], mons6o)
rAd = V.rank_exact(odd + [A, delta], mons6o)
rd = V.rank_exact(odd + [delta], mons6o)
log('C8 delta = det(dQ/dZ): 64 terms, Y-odd, == -X0X1X2Y0Y1Y2 mod 5; A not in J_6; delta not in J_6; A, delta independent mod J_6 (exact over Q)',
    'CONFIRMED' if (len(delta) == 64 and V.parity(delta) == 1 and d5 == {(1,1,1,1,1,1): 4} and rA == 19 and rd == 19 and rAd == 20) else 'REFUTED',
    f'terms {len(delta)}, parity {V.parity(delta)}, delta mod 5 = {d5}; ranks J6odd+A={rA}, +delta={rd}, +A+delta={rAd}')
open('/tmp/hodge/results/quadratic-cubic-pullback/verify/delta_exact.txt', 'w').write(str(delta_s))

# C9: Fermat control witness: strict identity E = 0 exactly
F0 = V.parse('X0**6+X1**6+X2**6+Y0**6+Y1**6+Y2**6')
D0 = V.parse('X0**6+X1**6+X2**6')
N = [V.parse('X0**4*Y0**3/6'), V.parse('X0**3*X1*Y0**3/6'), V.parse('X0**3*X2*Y0**3/6'), {}, {}, {}]
E = V.mul(A, V.mul(D0, D0))
for i in range(6):
    E = V.add(E, V.scale(V.mul(D0, V.mul(N[i], V.pderiv(F0, i))), -1))
    E = V.add(E, V.mul(F0, V.add(V.mul(D0, V.pderiv(N[i], i)), V.scale(V.mul(N[i], V.pderiv(D0, i)), -1))))
log('C9 Fermat-limit strict witness v = A E_X/(6(X0^6+X1^6+X2^6)) for (Fermat, 1, X0^3Y0^3): cleared identity E = 0 exactly',
    'CONFIRMED' if not E else 'REFUTED', f'E has {len(E)} terms')

# C10: F is not a join: rank of 21 second partials = 21 exactly
second = [V.pderiv(Fi_ex[i], j) for i in range(6) for j in range(i, 6)]
r2 = V.rank_exact(second, V.monomials(4))
log('C10 F is not a join G(l)+H(l\') in any coordinates (any split k+(6-k)): the 21 second partials are linearly independent (exact)',
    'CONFIRMED' if r2 == 21 else 'REFUTED', f'rank of 21 second partials over Q = {r2}; a k+(6-k) join has rank <= 21 - k(6-k)')

# C11: Q_0..Q_5 regular sequence: Hilbert function of S/(Q) = (1+t)^6 mod p up to degree 7, parity split 1 | 6+9 | 9+6 | 0+1
Qp = [V.modp(Qi, p) for Qi in Q]
hf = []
for k in range(0, 8):
    for par_ in (0, 1):
        mons = V.monomials(k, par_)
        gens = []
        for Qi in Qp:
            for m in V.monomials(k-2, par_):
                gens.append(V.shift(Qi, m, 1, p))
        r = V.rank_mod_p(gens, mons, p) if k >= 2 else 0
        hf.append((k, par_, len(mons) - r))
hfd = {(k, s): v for k, s, v in hf}
exp = {(0,0):1,(0,1):0,(1,0):3,(1,1):3,(2,0):6,(2,1):9,(3,0):10,(3,1):10,(4,0):9,(4,1):6,(5,0):3,(5,1):3,(6,0):0,(6,1):1,(7,0):0,(7,1):0}
log('C11 Q_0..Q_5 regular sequence, B = S/(Q) Hilbert series (1+t)^6, tau-trace (1-t^2)^3: (phi_*O)^+ = O + 6O(-1) + 9O(-2)',
    'CONFIRMED' if hfd == exp else 'REFUTED', f'bigraded Hilbert function of S/(Q) mod {p}: {hfd}')

# C12: cubics in Q (56-dim) meet J_6^even exactly in C*F
cubQ = []
for c in itertools.combinations_with_replacement(range(6), 3):
    P = {(0,)*6: 1}
    for i in c: P = V.mul(P, Qp[i], p)
    cubQ.append(P)
evenp = [V.modp(P, p) for P in even]
rc = V.rank_mod_p(cubQ, mons6e, p); rj = V.rank_mod_p(evenp, mons6e, p); ru = V.rank_mod_p(cubQ + evenp, mons6e, p)
log('C12 the 56-dim space of cubics in Q_0..Q_5 meets J_6^even exactly in the line C*F',
    'CONFIRMED' if (rc, rj, rc + rj - ru) == (56, 18, 1) else 'REFUTED', f'rank cubics {rc}, rank J6even {rj}, intersection dim {rc + rj - ru} (mod {p})')
print('done')
