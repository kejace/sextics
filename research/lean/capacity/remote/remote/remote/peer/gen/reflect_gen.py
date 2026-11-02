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

def spoly_lean(P):
    """Ascending lexicographic order on exponent tuples (matches Mono.lt), integer coefficients only."""
    items = []
    for e, c in sorted(P.items()):
        c = Fraction(c)
        assert c.denominator == 1, 'integral coefficients required in this prototype'
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
    L.append(f'def rF : SPoly := {spoly_lean(F)}\n')
    L.append(f'def rA : SPoly := {spoly_lean(A)}\n')
    L.append(f'def rD : SPoly := {spoly_lean(D)}\n')
    L.append(f'def rQ : SPoly := {spoly_lean(Q)}\n')
    for i in range(NV):
        L.append(f'def rN{i} : SPoly := {spoly_lean(N[i])}\n')
    L.append('def F : Poly := eval rF\ndef A : Poly := eval rA\ndef D : Poly := eval rD\ndef Q : Poly := eval rQ\n')
    for i in range(NV):
        L.append(f'def N{i} : Poly := eval rN{i}\n')
    L.append('def N : Fin 6 → Poly := fun i => ' + ''.join(f'if i.val = {i} then N{i} else ' for i in range(5)) + 'N5\n')
    if 'e1' in which:
        lhs_r = f'(smul {q} {nest_add([f"(mul (pderivP {i} rF) rN{i})" for i in range(NV)])})'
        rhs_r = '(mul rA rD)'
        L.append(f'theorem e1_core : eval {lhs_r} = eval {rhs_r} := eval_eq_of_eq (by decide +kernel)\n')
        e1_sum = nest([f'pderiv ({i} : Fin 6) F * N{i}' for i in range(NV)])
        L.append(f'theorem e1 : ({q} : Poly) * directional F N = A*D := by\n'
                 f'  simp only [directional, Fin.sum_univ_succ, Fin.sum_univ_zero, add_zero]\n'
                 f'  change ({q}:Poly) * ({e1_sum}) = A*D\n'
                 f'  have h := e1_core\n'
                 f'  simp only [eval_smul, eval_add, eval_mul, eval_pderivP, Int.cast_one, Int.cast_ofNat, map_one, map_ofNat] at h\n'
                 f'  exact h\n')
    if 'e2' in which:
        div_r = nest_add([f'(pderivP {i} rN{i})' for i in range(NV)])
        dir_r = nest_add([f'(mul (pderivP {i} rD) rN{i})' for i in range(NV)])
        lhs_r = f'(sub (mul rD {div_r}) {dir_r})'
        rhs_r = f'(mul (pow rF {q}) rQ)'
        L.append(f'theorem e2_core : eval {lhs_r} = eval {rhs_r} := eval_eq_of_eq (by decide +kernel)\n')
        div_sum = nest([f'pderiv ({i} : Fin 6) N{i}' for i in range(NV)])
        dir_sum = nest([f'pderiv ({i} : Fin 6) D * N{i}' for i in range(NV)])
        L.append(f'theorem e2 : divergenceNumerator D N = F^{q}*Q := by\n'
                 f'  simp only [divergenceNumerator, divergence, directional, Fin.sum_univ_succ, Fin.sum_univ_zero, add_zero]\n'
                 f'  change D * ({div_sum}) - ({dir_sum}) = F^{q}*Q\n'
                 f'  have h := e2_core\n'
                 f'  simp only [eval_sub, eval_mul, eval_add, eval_pderivP, eval_pow] at h\n'
                 f'  exact h\n')
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
