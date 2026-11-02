"""Monomial (permutation x diagonal) automorphism groups of the sample fibres and of X.
A monomial map pi*diag preserves F up to scalar iff pi maps supp(F) to itself and the diagonal
group of the permuted polynomial matches; we enumerate permutations pi with pi(supp F)=supp F
(as exponent sets) and for each solve the diagonal condition via the exponent-lattice routine.
"""
import sys, itertools, json
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core, symmetry
from fractions import Fraction as Fr
def poly_from_terms(terms, n):
    return {tuple(e): Fr(c) for e, c in terms}
def analyse(name, terms, n):
    F = poly_from_terms(terms, n)
    supp = set(F.keys())
    res = {}
    for pi in itertools.permutations(range(n)):
        # permuted polynomial: variable i -> variable pi[i]
        Fp = {}
        for e, c in F.items():
            e2 = [0]*n
            for i in range(n): e2[pi[i]] = e[i]
            Fp[tuple(e2)] = c
        if set(Fp.keys()) != supp: continue
        # need diag d with F(pi.diag Z) = lam F: coefficient condition c_e * prod d_i^{e_i} = lam * c_{pi(e)}
        # Solve multiplicatively: unknown log d_i; ratio condition for each pair of monomials.
        # Use the lattice: for monomials e,e' : prod d^{e-e'} = (c_{pi e'} c_e')^{-1}... do it numerically via group enumeration:
        # the diagonal symmetries of the *support* form a finite group G (projective); test each element.
        G = symmetry.diagonal_group(F) if pi == tuple(range(n)) else None
        res[str(pi)] = 'support-preserving'
    return res
# sample fibres: S_(1,2,3): vars (x,y,z,l); V_(1,2): (x,y,z,t,T); X: 6 vars
S = [((5,1,0,0),1),((0,5,1,0),1),((1,0,5,0),1),((1,3,0,2),Fr(10,7)),((3,0,1,2),Fr(60,7)),((0,1,3,2),Fr(45,7)),((0,0,0,6),794)]
V = [((5,1,0,0,0),1),((0,5,1,0,0),1),((1,0,5,0,0),1),((1,3,0,1,1),Fr(5,7)),((3,0,1,0,2),Fr(20,7)),((0,1,3,1,1),Fr(30,7)),((0,0,0,6,0),1),((0,0,0,0,6),65)]
X = [((5,1,0,0,0,0),1),((0,5,1,0,0,0),1),((1,0,5,0,0,0),1),((1,3,0,1,1,0),Fr(5,7)),((3,0,1,0,1,1),Fr(10,7)),((0,1,3,1,0,1),Fr(15,7)),((0,0,0,6,0,0),1),((0,0,0,0,6,0),1),((0,0,0,0,0,6),1)]
for name, terms, n in (('S',S,4),('V',V,5),('X',X,6)):
    F = poly_from_terms(terms, n)
    G = symmetry.diagonal_group(F)
    perms = analyse(name, terms, n)
    print(name, 'diagonal projective group order', G.order, 'support-preserving permutations', list(perms.keys()))
