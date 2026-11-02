"""Emit a reflection-checked module Cap/R_<label>.lean for the same synthetic data as synth.py.

The identities e1/e2 are stated exactly as in the accepted template (so the rest of the template can be reused),
but proved via `HodgeComputer.Reflect`: F := eval rF etc., and the core equalities by `decide +kernel`.

Usage: python3 reflect_gen.py --F m3-coupled --d 4 --nD 40 [--q 1] [--seed 0] [--which e1,e2]
"""
import sys, json, argparse
sys.path.insert(0, '/tmp/hodge/lean/capacity/gen'); sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core
from hodge.core import NV
from fractions import Fraction
import ident, synth
CAP = '/tmp/hodge/lean/capacity'

import math
def denom(P):
    k = 1
    for c in P.values():
        k = math.lcm(k, Fraction(c).denominator)
    return k

def spoly_lean(P, k=1):
    """k*P as an integer list in ascending lexicographic order on exponent tuples (matches Mono.lt)."""
    items = []
    for e, c in sorted(P.items()):
        c = Fraction(c) * k
        assert c.denominator == 1, 'scaling by the denominator LCM must give integers'
        if c == 0:
            continue
        items.append(f'(⟨{e[0]}, {e[1]}, {e[2]}, {e[3]}, {e[4]}, {e[5]}⟩, {c.numerator})')
    if not items:
        return '[]'
    return '[' + ',\n  '.join(items) + ']'

def nest(terms):
    return terms[0] if len(terms) == 1 else f'{terms[0]} + ({nest(terms[1:])})'

def nest_add(terms):
    return terms[0] if len(terms) == 1 else f'(add {terms[0]} {nest_add(terms[1:])})'

def module_text(ns, F, A, q, D, N, Q, which=('e1', 'e2'), hb=2000000, rd=8000):
    L = []
    L.append('import Cap.Reflect\n')
    L.append(f'set_option Elab.async false\nset_option maxHeartbeats {hb}\nset_option maxRecDepth {rd}\n'
             'set_option linter.unusedTactic false\nset_option linter.unusedSimpArgs false\n'
             'set_option linter.unnecessarySeqFocus false\nset_option linter.unreachableTactic false\n')
    L.append(f'noncomputable section\n\nnamespace {ns}\n\nopen MvPolynomial HodgeTestbed HodgeComputer.Reflect\n')
    kF, kA, kD, kQ = denom(F), denom(A), denom(D), denom(Q)
    kN = 1
    for i in range(NV):
        kN = math.lcm(kN, denom(N[i]))
    L.append(f'def rF : SPoly := {spoly_lean(F, kF)}\n')
    L.append(f'def rA : SPoly := {spoly_lean(A, kA)}\n')
    L.append(f'def rD : SPoly := {spoly_lean(D, kD)}\n')
    L.append(f'def rQ : SPoly := {spoly_lean(Q, kQ)}\n')
    for i in range(NV):
        L.append(f'def rN{i} : SPoly := {spoly_lean(N[i], kN)}\n')
    L.append(f'def F : Poly := evalQ rF {kF}\ndef A : Poly := evalQ rA {kA}\ndef D : Poly := evalQ rD {kD}\ndef Q : Poly := evalQ rQ {kQ}\n')
    for i in range(NV):
        L.append(f'def N{i} : Poly := evalQ rN{i} {kN}\n')
    L.append('def N : Fin 6 → Poly := fun i => ' + ''.join(f'if i.val = {i} then N{i} else ' for i in range(5)) + 'N5\n')
    qcast = '' if q == 1 else f'  rw [show ({q} : Poly) = (({q} : ℕ) : Poly) by norm_num]\n'
    if 'e1' in which:
        e1_sum = nest([f'pderiv ({i} : Fin 6) F * N{i}' for i in range(NV)])
        L.append(f'theorem e1 : ({q} : Poly) * directional F N = A*D := by\n'
                 f'  simp only [directional, Fin.sum_univ_succ, Fin.sum_univ_zero, add_zero]\n'
                 f'  change ({q}:Poly) * ({e1_sum}) = A*D\n' + qcast +
                 f'  simp only [F, A, D, N0, N1, N2, N3, N4, N5, evalQ_pderiv, evalQ_mul, evalQ_add, evalQ_nsmul, one_mul]\n'
                 f'  exact evalQ_eq_of_eq (by decide) (by decide) (by decide +kernel)\n')
    if 'e2' in which:
        div_sum = nest([f'pderiv ({i} : Fin 6) N{i}' for i in range(NV)])
        dir_sum = nest([f'pderiv ({i} : Fin 6) D * N{i}' for i in range(NV)])
        L.append(f'theorem e2 : divergenceNumerator D N = F^{q}*Q := by\n'
                 f'  simp only [divergenceNumerator, divergence, directional, Fin.sum_univ_succ, Fin.sum_univ_zero, add_zero]\n'
                 f'  change D * ({div_sum}) - ({dir_sum}) = F^{q}*Q\n'
                 f'  simp only [F, D, Q, N0, N1, N2, N3, N4, N5, evalQ_pderiv, evalQ_mul, evalQ_add, evalQ_sub, evalQ_pow]\n'
                 f'  exact evalQ_eq_of_eq (by decide) (by decide) (by decide +kernel)\n')
    L.append(f'end {ns}\n')
    return '\n'.join(L)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--F', default='m3-coupled'); ap.add_argument('--d', type=int, default=4)
    ap.add_argument('--nD', type=int, default=40); ap.add_argument('--q', type=int, default=1)
    ap.add_argument('--seed', type=int, default=0); ap.add_argument('--which', default='e1,e2')
    ap.add_argument('--yeven', action='store_true'); ap.add_argument('--label', default=None)
    a = ap.parse_args()
    F, A, q, D, N, Q = synth.build(a.F, a.d, a.nD, a.q, a.seed, a.yeven)
    which = tuple(a.which.split(','))
    label = a.label or f'{synth.FSHORT.get(a.F, a.F)}_d{a.d}_n{len(D)}' + (f'_q{q}' if q != 1 else '') + ('' if which == ('e1', 'e2') else '_' + ''.join(which))
    ns = f'HodgeComputer.Instances.R_{label}'
    st = ident.stats(F, A, q, D, N, Q)
    st.update({'label': 'R_' + label, 'Fid': a.F, 'which': list(which), 'seed': a.seed, 'route': 'reflect'})
    txt = module_text(ns, F, A, q, D, N, Q, which)
    open(f'{CAP}/src/Cap/R_{label}.lean', 'w').write(txt)
    st['bytes'] = len(txt.encode())
    json.dump(st, open(f'{CAP}/src/Cap/R_{label}.json', 'w'), indent=1)
    print('R_' + label, json.dumps(st))

if __name__ == '__main__':
    main()
