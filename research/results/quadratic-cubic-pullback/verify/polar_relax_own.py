"""Own-code Hessian relaxation for polar denominators D = w.grad F (deg w = d-5).
Necessary condition for ANY certificate with denominator D (derived independently):
  N(F) = A D - cF  =>  N = A w - (c/6) Z + K.grad F (Koszul)  =>  E2: cD + D div N - N(D) = F Q  =>  A * phi(w) in J_{2d},
  phi(w) := w^T Hess(F) w  (all other terms of N(D) lie in J).  With K_{2d-6} = J_{2d-6} (true for 2d-6 <= 9): phi(w) in J_{2d-6}.
d = 5: phi(w) in J_4 = 0 => phi(w) == 0; hand proof via the six pure-power coefficients. d = 6: injectivity test of the pair map."""
import sys, json, itertools
sys.path.insert(0, '/tmp/hodge/results/quadratic-cubic-pullback/verify')
import vpoly as V
from vlog import log
from fractions import Fraction
import sympy
F, A, q = V.load_instance_F_A()
Fi = [V.pderiv(F, i) for i in range(6)]
H = [[V.pderiv(Fi[i], j) for j in range(6)] for i in range(6)]
w = sympy.symbols('w0:6')
# --- d = 5: phi(w) = sum_ij w_i w_j H_ij, coefficients of the pure fourth powers
def coeff_pure(k):
    e = tuple(4 if t == k else 0 for t in range(6))
    return sympy.expand(sum(sympy.Rational(H[i][j].get(e, 0).numerator, H[i][j].get(e, 0).denominator)*w[i]*w[j] for i in range(6) for j in range(6)))
pure = [coeff_pure(k) for k in range(6)]
print('pure-power coefficients of phi(w):'); [print(' ', k, pure[k]) for k in range(6)]
# Y-even w = (w0,w1,w2,0,0,0): first three equations give w_i^2 = -(1/p_i) * (product of the other two) => (w0 w1 w2)^2 = -(w0w1w2)^2/(7*11*13) => w = 0
even_sys = [pure[k].subs({w[3]: 0, w[4]: 0, w[5]: 0}) for k in range(3)]
odd_sys = [pure[k].subs({w[0]: 0, w[1]: 0, w[2]: 0}) for k in range(3, 6)]
G_even = sympy.groebner(even_sys, w[0], w[1], w[2], order='grevlex')
G_odd = sympy.groebner(odd_sys, w[3], w[4], w[5], order='grevlex')
def only_origin(G, vars_):
    # homogeneous ideal; zero-dimensional affine <=> projectively empty <=> only w = 0
    return G.is_zero_dimensional
ok_e, ok_o = only_origin(G_even, w[:3]), only_origin(G_odd, w[3:])
log('C27 no Y-even degree-5 denominator exists (own derivation of the Hessian relaxation + hand proof from the X_i^4 coefficients: 30 w_i^2 + (30/p_i) w_j w_k = 0 => w=0; uses K_4=0)',
    'CONFIRMED' if ok_e else 'REFUTED', f'even eqs: {even_sys}; zero-dim: {ok_e}')
log('C27b no Y-odd degree-5 denominator exists (same argument with the Y_i^4 coefficients, p_i = 17,19,23)',
    'CONFIRMED' if ok_o else 'REFUTED', f'odd eqs: {odd_sys}; zero-dim: {ok_o}')
# mixed parity d=5: all coefficients of phi(w) (full system) -- Groebner over Q
mons4 = V.monomials(4)
eqs = []
for e in mons4:
    expr = sympy.expand(sum(sympy.Rational(H[i][j].get(e, 0).numerator, H[i][j].get(e, 0).denominator)*w[i]*w[j] for i in range(6) for j in range(6)))
    if expr != 0: eqs.append(expr)
G_mixed = sympy.groebner(eqs, *w, order='grevlex')
log('C27c no degree-5 denominator of mixed parity either (full phi(w)=0 system in 6 unknowns, Groebner over Q is zero-dimensional)',
    'CONFIRMED' if G_mixed.is_zero_dimensional else 'REFUTED', f'{len(eqs)} quadratic equations; zero-dim: {G_mixed.is_zero_dimensional}; GB size {len(G_mixed.exprs)}')
# --- d = 6: w_i linear. pure parity modes. pair map Sym^2(W_gauge-fixed) -> R_6^even must be injective.
p = 1000000007
Hp = [[V.modp(H[i][j], p) for j in range(6)] for i in range(6)]
mons6e = V.monomials(6, 0)
J6e = []
for i in range(6):
    for m in V.monomials(1, (0 - (0 if i < 3 else 1)) % 2):
        J6e.append(V.shift(V.modp(Fi[i], p), m, 1, p))
assert V.rank_mod_p(J6e, mons6e, p) == 18
for mode, pD in (('even', 0), ('odd', 1)):
    basis = []
    for i in range(6):
        par_i = (pD - (0 if i < 3 else 1)) % 2
        for m in V.monomials(1, par_i): basis.append((i, m))
    # gauge: w = l*Z with l constant of parity pD-? : l*X_i must have parity par_0 = pD -> l constant exists only if pD = 0
    if pD == 0:
        # drop the basis element (0, X0) and replace: gauge vector has all (i, Z_i) entries 1; fix gauge by removing (0, X0)
        basis_fixed = [b for b in basis if b != (0, (1, 0, 0, 0, 0, 0))]
    else:
        basis_fixed = basis
    pairs = []
    for a in range(len(basis_fixed)):
        for b in range(a, len(basis_fixed)):
            (i, mi), (j, mj) = basis_fixed[a], basis_fixed[b]
            P_ = V.shift(Hp[i][j], tuple(x + y for x, y in zip(mi, mj)), 1, p)
            pairs.append(P_)
    r_all = V.rank_mod_p(pairs + J6e, mons6e, p)
    inj = (r_all - 18 == len(pairs))
    log(f'C28 no {mode} degree-6 polar denominator passes E2: pair map Sym^2(W) -> R_6^even injective (own code, p={p}); w-space dim {len(basis)}, gauge-fixed {len(basis_fixed)}, pairs {len(pairs)}',
        'CONFIRMED' if inj else 'UNVERIFIED', f'rank(pairs+J6even)={r_all}, minus 18 = {r_all-18} vs pairs {len(pairs)}; injective: {inj}')
