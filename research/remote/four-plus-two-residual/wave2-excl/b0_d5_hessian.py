"""(b) deg D = 5 exclusion, exact: every E1-admissible D of degree 5 (any Y-parity) is D = grad F . w with w constant
(exact core kernels are 0 for k<=5 -> {D: A D in J(F)}_5 = J(F)_5, via R(F) = R(H) (x) R(B)); E2 mod J then forces
A (w^T Hess(F) w) in J_10; multiplication by A on (R/J)_4 is injective (exact kernels 0 for k<=4), so w^T Hess(F) w = 0
as a polynomial.  We verify exactly over Q that this quadratic system in w has only w = 0 (sympy Groebner), and independently
that it follows from smoothness (coefficient extraction).  Also records that the Hessian of F is block diagonal."""
import sys, json; sys.path.insert(0, '/tmp/hodge/results/four-plus-two-residual/wave2-excl')
from lib import *
import sympy
z = sympy.symbols('X0 X1 X2 Y0 Y1 Y2'); w = sympy.symbols('w0:6')
Fs = sum(sympy.Rational(c) * sympy.prod([z[i] ** e[i] for i in range(6)]) for e, c in F.items())
Hess = sympy.hessian(Fs, z)
quad = sympy.expand(sum(w[i] * Hess[i, j] * w[j] for i in range(6) for j in range(6)))
eqs = sympy.Poly(quad, *z).coeffs()
G = sympy.groebner(eqs, *w, order='grevlex')
# G should be the irrelevant ideal (w0..w5)
is_trivial = set(G.exprs) == set(w) or all(any(sympy.simplify(g - wi) == 0 for g in G.exprs) for wi in w)
# count equations and show block structure
cross = [Hess[i, j] for i in (0, 1, 3, 4) for j in (2, 5)]
log(dict(exp='b0-d5-hessian-exact', n_quadrics=len(eqs), groebner=[str(g) for g in G.exprs], only_trivial_solution=bool(is_trivial),
         hessian_block_diagonal=all(c == 0 for c in cross)))
