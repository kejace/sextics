"""Shared helpers: identity size statistics and Lean emission in the accepted template."""
import sys, json
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core, leangen
from hodge.core import NV

def sumprod(P, N):
    s = {}
    for i in range(NV):
        s = core.add(s, core.mul(P[i], N[i]))
    return s

def stats(F, A, q, D, N, Q):
    dF = [core.pderiv(F, i) for i in range(NV)]
    dD = [core.pderiv(D, i) for i in range(NV)]
    e1_pairs = sum(len(dF[i]) * len(N[i]) for i in range(NV)) + len(A) * len(D)
    lhs1 = core.scale(sumprod(dF, N), q)
    rhs1 = core.mul(A, D)
    divN = {}
    for i in range(NV):
        divN = core.add(divN, core.pderiv(N[i], i))
    Fq = core.power(F, q)
    e2_pairs = len(D) * len(divN) + sum(len(dD[i]) * len(N[i]) for i in range(NV)) + len(Fq) * len(Q)
    lhs2 = core.sub(core.mul(D, divN), sumprod(dD, N))
    rhs2 = core.mul(Fq, Q)
    # raw (pre-cancellation) product monomials: count of all pairwise product terms = pairs; distinct monomials:
    return {'F': len(F), 'A': len(A), 'q': q, 'D': len(D), 'N': [len(n) for n in N], 'Q': len(Q), 'd': core.degree(D),
            'e1_pairs': e1_pairs, 'e1_monos': len(lhs1), 'e2_pairs': e2_pairs, 'e2_monos': len(lhs2),
            'e1_ok': core.sub(lhs1, rhs1) == {}, 'e2_ok': core.sub(lhs2, rhs2) == {},
            'pairs_total': e1_pairs + e2_pairs, 'monos_total': len(lhs1) + len(lhs2)}

HEADER = '''import HodgeComputer.Verify

set_option Elab.async false
set_option maxHeartbeats {hb}
set_option maxRecDepth {rd}
set_option linter.unusedTactic false
set_option linter.unusedSimpArgs false
set_option linter.unnecessarySeqFocus false
set_option linter.unreachableTactic false

noncomputable section

namespace {ns}

open MvPolynomial HodgeTestbed
'''

def nest(terms):
    return terms[0] if len(terms) == 1 else f'{terms[0]} + ({nest(terms[1:])})'

def e1_theorem(q):
    e1_sum = nest([f'pderiv ({i} : Fin 6) F * N{i}' for i in range(NV)])
    return (f'theorem e1 : ({q} : Poly) * directional F N = A*D := by\n'
            f'  simp only [directional, Fin.sum_univ_succ, Fin.sum_univ_zero, add_zero]\n'
            f'  change ({q}:Poly) * ({e1_sum}) = A*D\n'
            f'  simp only [F,A,D,N0,N1,N2,N3,N4,N5,map_zero,map_add,pderiv_C,pderiv_mul,pderiv_pow,{leangen.DERIV_NAMES}]\n'
            f'  apply MvPolynomial.funext; intro x\n'
            f'  norm_num <;> ring\n')

def e2_theorem(q):
    div_sum = nest([f'pderiv ({i} : Fin 6) N{i}' for i in range(NV)])
    dir_sum = nest([f'pderiv ({i} : Fin 6) D * N{i}' for i in range(NV)])
    return (f'theorem e2 : divergenceNumerator D N = F^{q}*Q := by\n'
            f'  simp only [divergenceNumerator, divergence, directional, Fin.sum_univ_succ, Fin.sum_univ_zero, add_zero]\n'
            f'  change D * ({div_sum}) - ({dir_sum}) = F^{q}*Q\n'
            f'  simp only [F,D,Q,N0,N1,N2,N3,N4,N5,map_zero,map_add,pderiv_C,pderiv_mul,pderiv_pow,{leangen.DERIV_NAMES}]\n'
            f'  apply MvPolynomial.funext; intro x\n'
            f'  norm_num <;> ring\n')

def defs_block(F, A, D, Q, N):
    L = [f'def F : Poly := {leangen._poly_lean(F)}\n', f'def A : Poly := {leangen._poly_lean(A)}\n',
         f'def D : Poly := {leangen._poly_lean(D)}\n', f'def Q : Poly := {leangen._poly_lean(Q)}\n']
    for i in range(NV):
        L.append(f'def N{i} : Poly := {leangen._poly_lean(N[i])}\n')
    L.append('def N : Fin 6 → Poly := fun i => ' + ''.join(f'if i.val = {i} then N{i} else ' for i in range(5)) + 'N5\n')
    return '\n'.join(L)

def module_text(ns, F, A, q, D, N, Q, which=('e1', 'e2'), hb=2000000, rd=4000):
    L = [HEADER.format(hb=hb, rd=rd, ns=ns), defs_block(F, A, D, Q, N), leangen._deriv_lemmas()]
    if 'e1' in which: L.append(e1_theorem(q))
    if 'e2' in which: L.append(e2_theorem(q))
    L.append(f'end {ns}\n')
    return '\n'.join(L)
