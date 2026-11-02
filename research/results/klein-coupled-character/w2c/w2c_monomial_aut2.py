"""Monomial automorphisms (permutation * diagonal, projective) of a polynomial with rational coefficients.
For each support-preserving permutation pi, we need d in (C*)^n with  c_e prod_i d_i^{e_i} = lam c_{pi(e)} for all e.
Taking ratios between monomials e and a base e0 kills lam:  prod d^{e-e0} = c_{pi e} c_{e0} / (c_e c_{pi e0}).
Let L = lattice spanned by {e-e0}; the equation is a homomorphism L -> C*, extended to Z^n; solutions exist iff
the assignment is consistent on relations of L (checked via Smith normal form), and then they form a torsor
under Hom(Z^n/L, C*) = diagonal symmetry group (finite modulo the Euler scalars if F is 'nondegenerate').
We compute: order of the projective diagonal group, and for each pi whether a compatible diagonal exists.
Consistency check: for integer relation vectors r (kernel of the exponent-difference matrix), the product of the
RHS ratios^r must equal 1 (exactly, in Q).  Output counts.
"""
import itertools, sys
from fractions import Fraction as Fr
import sympy as sp
def smith_order(rows, n):
    M = sp.Matrix(rows)  # rows = exponent differences, plus the all-ones row (projective)
    M = M.col_join(sp.Matrix([[1]*n]))
    from sympy.matrices.normalforms import smith_normal_form
    D = smith_normal_form(M, domain=sp.ZZ)
    diag = [abs(D[i,i]) for i in range(min(D.shape)) if D[i,i] != 0]
    order = 1
    for d in diag: order *= d
    return order, diag
def monomial_auts(terms, n):
    F = {tuple(e): Fr(c) for e,c in terms}
    supp = list(F.keys()); e0 = supp[0]
    rows = [[e[i]-e0[i] for i in range(n)] for e in supp[1:]]
    order, diag = smith_order(rows, n)
    out = {'diag_order': order, 'snf': diag, 'perms': []}
    A = sp.Matrix(rows + [[1]*n])          # relations among the difference vectors: left kernel of A
    K = A.T.nullspace()                   # integer relations r with r^T A = 0  (rational basis)
    for pi in itertools.permutations(range(n)):
        def img(e):
            e2=[0]*n
            for i in range(n): e2[pi[i]] = e[i]
            return tuple(e2)
        if set(img(e) for e in supp) != set(supp): continue
        # ratio for monomial e (vs e0): rho_e = c_{pi e} c_{e0} / (c_e c_{pi e0})
        rho = [F[img(e)]*F[e0]/(F[e]*F[img(e0)]) for e in supp[1:]] + [Fr(1)]   # last row = Euler, ratio 1 (lam absorbed)
        ok = True
        for r in K:
            r = r * sp.lcm([sp.fraction(x)[1] for x in r])
            val = Fr(1)
            for rk, rh in zip(r, rho):
                val *= rh ** int(rk)
            if val != 1: ok = False; break
        out['perms'].append((pi, ok))
    return out
S = [((5,1,0,0),1),((0,5,1,0),1),((1,0,5,0),1),((1,3,0,2),Fr(10,7)),((3,0,1,2),Fr(60,7)),((0,1,3,2),Fr(45,7)),((0,0,0,6),794)]
V = [((5,1,0,0,0),1),((0,5,1,0,0),1),((1,0,5,0,0),1),((1,3,0,1,1),Fr(5,7)),((3,0,1,0,2),Fr(20,7)),((0,1,3,1,1),Fr(30,7)),((0,0,0,6,0),1),((0,0,0,0,6),65)]
X = [((5,1,0,0,0,0),1),((0,5,1,0,0,0),1),((1,0,5,0,0,0),1),((1,3,0,1,1,0),Fr(5,7)),((3,0,1,0,1,1),Fr(10,7)),((0,1,3,1,0,1),Fr(15,7)),((0,0,0,6,0,0),1),((0,0,0,0,6,0),1),((0,0,0,0,0,6),1)]
for name, terms, n in (('S_(1,2,3)',S,4),('V_(1,2)',V,5),('X',X,6)):
    r = monomial_auts(terms, n)
    print(name, 'projective diagonal group order', r['diag_order'], 'SNF', r['snf'])
    print('   support-preserving permutations (pi, diagonal-compatible):', r['perms'])
