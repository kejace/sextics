"""Formal check of the q-raising identity (dossier §5 P2 / A.7) as an identity of rational functions, with sympy, in a reduced
setting (3 variables suffice for the formal algebra; the identity only uses the Leibniz/Euler-free manipulations)."""
import sympy as sp
x, y, z, lam = sp.symbols('x y z lam')
vars_ = [x, y, z]
# generic polynomials with symbolic coefficients would be heavy; use random explicit polynomials in Q[x,y,z,lam]
import random
random.seed(3)
def rpoly(deg, with_lam=True):
    terms = []
    for i in range(deg + 1):
        for j in range(deg + 1 - i):
            k = deg - i - j
            c = random.randint(-3, 3) + (random.randint(-2, 2) * lam if with_lam else 0)
            terms.append(c * x**i * y**j * z**k)
    return sp.expand(sum(terms))
M1 = rpoly(2, False); B = rpoly(3, False); F0 = rpoly(2, False)
F = sp.expand(F0 + (lam - sp.Rational(1, 2)) * M1)        # dF/dlam = M1
A = sp.expand(M1 * B)
D = rpoly(2); N = [rpoly(3) for _ in vars_]
v = [n / D for n in N]
div = lambda w: sum(sp.diff(w[i], vars_[i]) for i in range(3))
app = lambda w, P: sum(w[i] * sp.diff(P, vars_[i]) for i in range(3))
K1 = B - app(v, F) + F * div(v)                 # q=1 cleared source; define Q := K1 / F^2 (rational)
Q = K1 / F**2
vdot = [sp.diff(vi, lam) for vi in v]
w = [sp.Rational(1, 2) * (M1 * v[i] - F * vdot[i]) for i in range(3)]
K2 = A - 2 * app(w, F) + F * div(w)
lhs = sp.together(K2 + sp.Rational(1, 2) * F**3 * sp.diff(Q, lam))
num = sp.expand(sp.numer(lhs))
print('K2(w) + 1/2 F^3 dQ/dlam == 0 identically:', num == 0)
import sys; sys.path.insert(0, '/tmp/hodge/results/cycle-q2/verify'); import vcheck as V
V.log(claim='q-raising identity: if v is a q=1 certificate for (F_lam, B) with F_lam = F + (lam-1/2) M1, A = M1 B, then w = 1/2 (M1 v - F v_dot) satisfies A - 2 w(F) + F div w = -1/2 F^3 Q_dot (q=2 certificate with denominator D^2)',
      source='dossier §5 P2 / A.7', verdict='CONFIRMED' if num == 0 else 'REFUTED', evidence={'method': 'sympy rational-function identity with random polynomials in 3 variables; Q defined as K1/F^2', 'zero': bool(num == 0)},
      notes='also checked by hand: K2(w) = A + c(2A + F^3 Q_dot) for w = c(F v_dot - M1 v); degrees deg N = d+1 at q=1 give deg w-numerators = 2d+7 = deg(D^2)+7 as required; needs D_{1/2} != 0 and gcd(D_{1/2},F)=1')
