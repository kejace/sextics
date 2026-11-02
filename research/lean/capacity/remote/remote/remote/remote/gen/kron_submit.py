#!/usr/bin/env python3
"""Produce a complete solve proposal (Kron.lean, Witness.lean, Solution.lean, proof.json, witness.json) for an
instance from an explicit witness, using the Kronecker checker for the identities and the W-AST binding for the
protected witness terms.  The output folder can be checked with `python3 -m hodge.leancheck --instance <id> --dir <out>`
(local replica of the verifier) or `python3 -m hodge.sandboxcheck <id> <out>` (official Docker sandbox).

Usage: python3 kron_submit.py --instance m3-coupled --out /tmp/hodge/lean/capacity/submit/m3-coupled
       (witness taken from data/<instance>/witness.json unless --witness <json> is given)
"""
from __future__ import annotations
import argparse, ast, json, math, os, sys
from fractions import Fraction
sys.path.insert(0, '/tmp/hodge/lean/capacity/gen'); sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core, leangen
from hodge.core import NV, VARS
import kron_gen as K
import bind_test as B

KRON_SRC = '/tmp/hodge/lean/capacity/src/Reflect/Kron.lean'


def lcm_den(*polys):
    s = 1
    for P in polys:
        for c in P.values():
            den = Fraction(c).denominator
            s = s * den // math.gcd(s, den)
    return s


def nest(terms):
    return terms[0] if len(terms) == 1 else f'{terms[0]} + ({nest(terms[1:])})'


def nestE(op, terms):
    return terms[0] if len(terms) == 1 else f'(.{op} {terms[0]} {nestE(op, terms[1:])})'


def generate(instance_id, witness, out, name='KronWitness', heartbeats=4000000, recdepth=8000):
    F, A, q, rec = core.load_instance(instance_id)
    D = core.parse(witness['D']); N = [core.parse(s) for s in witness['N']]; Q = core.parse(witness['Q'])
    d = int(witness['d'])
    ok, info = core.verify_witness(F, A, q, D, N, Q)
    assert ok, info
    # integral witness: scale D, N by s and Q by s^2
    s = lcm_den(D, Q, *N)
    D = core.scale(D, s); N = [core.scale(n, s) for n in N]; Q = core.scale(Q, s * s)
    ok, info = core.verify_witness(F, A, q, D, N, Q)
    assert ok, info
    Di = {e: int(c) for e, c in D.items()}; Ni = [{e: int(c) for e, c in n.items()} for n in N]
    Qi = {e: int(c) for e, c in Q.items()}
    # F, A scaling
    sFA = lcm_den(F, A)
    Fi = {e: int(Fraction(c) * sFA) for e, c in F.items()}; Ai = {e: int(Fraction(c) * sFA) for e, c in A.items()}
    degN = d + 6 * q - 5
    deg1 = d + 6 * q; deg2 = 2 * d + 6 * q - 6
    m = max(deg1, deg2) + 1
    l1 = K.l1
    dF = [core.pderiv(Fi, i) for i in range(NV)]; dD = [core.pderiv(Di, i) for i in range(NV)]
    dN = [core.pderiv(Ni[i], i) for i in range(NV)]
    b1 = q * sum(l1(dF[i]) * l1(Ni[i]) for i in range(NV)) + l1(Ai) * l1(Di)
    b2 = (sFA ** q) * (l1(Di) * sum(l1(dN[i]) for i in range(NV)) + sum(l1(dD[i]) * l1(Ni[i]) for i in range(NV))) + l1(Fi) ** q * l1(Qi)
    W = max(max(b1, b2).bit_length() + 2, 8)
    ns = f'HodgeComputer.Instances.{name}'
    os.makedirs(out, exist_ok=True)
    # ---------------- Kron.lean (library, flat helper module)
    kron = open(KRON_SRC).read()
    open(os.path.join(out, 'Kron.lean'), 'w').write(kron)
    # ---------------- witness.json (balanced strings, integral coefficients)
    def wstr(P, deg):
        return B.witness_string(P, m) if P else '0'
    wjson = {'d': d, 'D': wstr(Di, d), 'N': [wstr(Ni[i], degN) for i in range(NV)], 'Q': wstr(Qi, 2 * d - 6)}
    json.dump(wjson, open(os.path.join(out, 'witness.json'), 'w'), indent=1)
    # ---------------- Witness.lean
    L = ['import Proposal.Kron\n', f'set_option maxHeartbeats {heartbeats}\nset_option maxRecDepth {recdepth}\n',
         'set_option linter.all false\n', f'namespace {ns}\n', 'open MvPolynomial HodgeTestbed HodgeComputer.Kron\n']
    L.append(f'noncomputable def tF : PT := {K.pt_lean(Fi, m, 6)}\n')
    L.append(f'noncomputable def tA : PT := {K.pt_lean(Ai, m, 6 * q)}\n')
    L.append(f'noncomputable def tD : PT := {K.pt_lean(Di, m, d)}\n')
    L.append(f'noncomputable def tQ : PT := {K.pt_lean(Qi, m, 2 * d - 6)}\n')
    for i in range(NV):
        L.append(f'noncomputable def tN{i} : PT := {K.pt_lean(Ni[i], m, degN)}\n')
    L.append('noncomputable section\n')
    if sFA == 1:
        L.append('def F : Poly := toC tF\ndef A : Poly := toC tA\n')
    else:
        L.append(f'def F : Poly := C (1 / {sFA} : ℂ) * toC tF\ndef A : Poly := C (1 / {sFA} : ℂ) * toC tA\n')
    L.append('def D : Poly := toC tD\ndef Q : Poly := toC tQ\n')
    for i in range(NV):
        L.append(f'def N{i} : Poly := toC tN{i}\n')
    L.append('def N : Fin 6 → Poly := fun i => ' + ''.join(f'if i.val = {i} then N{i} else ' for i in range(5)) + 'N5\n')
    side = '(by decide +kernel) (by decide +kernel) (by decide +kernel) (by decide +kernel) (by decide +kernel) (by decide +kernel)'
    # e1
    lhs = f'(.smulN {q} {nestE("add", [f"(.mul (.pd {i} tF) (.poly tN{i}))" for i in range(NV)])})'
    rhs = '(.mul (.poly tA) (.poly tD))'
    L.append(f'noncomputable def L1 : E := {lhs}\nnoncomputable def R1 : E := {rhs}\n')
    L.append(f'theorem e1_core : denC L1 = denC R1 := soundC (m := {m}) (W := {W}) L1 R1 {side}\n')
    e1_sum = nest([f'pderiv ({i} : Fin 6) F * N{i}' for i in range(NV)])
    if sFA == 1:
        L.append(f'theorem e1 : ({q} : Poly) * directional F N = A*D := by\n'
                 f'  simp only [directional, Fin.sum_univ_succ, Fin.sum_univ_zero, add_zero]\n'
                 f'  change ({q}:Poly) * ({e1_sum}) = A*D\n'
                 f'  have h := e1_core\n'
                 f'  simp only [L1, R1, denC_smulN, denC_add, denC_mul, denC_poly, denC_pd, Nat.cast_one, Nat.cast_ofNat] at h\n'
                 f'  simpa only [F, A, D, N0, N1, N2, N3, N4, N5] using h\n')
    else:
        L.append(f'theorem e1 : ({q} : Poly) * directional F N = A*D := by\n'
                 f'  simp only [directional, Fin.sum_univ_succ, Fin.sum_univ_zero, add_zero]\n'
                 f'  change ({q}:Poly) * ({e1_sum}) = A*D\n'
                 f'  have h := congrArg (fun p => C (1 / {sFA} : ℂ) * p) e1_core\n'
                 f'  simp only [L1, R1, denC_smulN, denC_add, denC_mul, denC_poly, denC_pd, Nat.cast_one, Nat.cast_ofNat] at h\n'
                 f'  simp only [F, A, D, N0, N1, N2, N3, N4, N5, pderiv_C_mul]\n'
                 f'  linear_combination h\n')
    # e2
    div_r = nestE('add', [f'(.pd {i} tN{i})' for i in range(NV)])
    dir_r = nestE('add', [f'(.mul (.pd {i} tD) (.poly tN{i}))' for i in range(NV)])
    lhs = f'(.sub (.mul (.poly tD) {div_r}) {dir_r})'
    if sFA != 1:
        lhs = f'(.smulN {sFA ** q} {lhs})'
    rhs = f'(.mul (.pow (.poly tF) {q}) (.poly tQ))'
    L.append(f'noncomputable def L2 : E := {lhs}\nnoncomputable def R2 : E := {rhs}\n')
    L.append(f'theorem e2_core : denC L2 = denC R2 := soundC (m := {m}) (W := {W}) L2 R2 {side}\n')
    div_sum = nest([f'pderiv ({i} : Fin 6) N{i}' for i in range(NV)])
    dir_sum = nest([f'pderiv ({i} : Fin 6) D * N{i}' for i in range(NV)])
    if sFA == 1:
        L.append(f'theorem e2 : divergenceNumerator D N = F^{q}*Q := by\n'
                 f'  simp only [divergenceNumerator, divergence, directional, Fin.sum_univ_succ, Fin.sum_univ_zero, add_zero]\n'
                 f'  change D * ({div_sum}) - ({dir_sum}) = F^{q}*Q\n'
                 f'  have h := e2_core\n'
                 f'  simp only [L2, R2, denC_sub, denC_add, denC_mul, denC_pow, denC_poly, denC_pd] at h\n'
                 f'  simpa only [F, D, Q, N0, N1, N2, N3, N4, N5] using h\n')
    else:
        sq = sFA ** q
        div_u = nest([f'pderiv ({i} : Fin 6) (toC tN{i})' for i in range(NV)])
        dir_u = nest([f'pderiv ({i} : Fin 6) (toC tD) * toC tN{i}' for i in range(NV)])
        L.append(f'theorem hF{q} : F ^ {q} = C (1 / {sq} : ℂ) * toC tF ^ {q} := by\n  rw [F, mul_pow, ← C_pow]; norm_num\n')
        L.append(f'theorem hs : (C (1 / {sq} : ℂ) : Poly) * ({sq} : Poly) = 1 := by\n  rw [← map_ofNat C {sq}, ← C_mul]; norm_num\n')
        L.append(f'theorem e2 : divergenceNumerator D N = F^{q}*Q := by\n'
                 f'  simp only [divergenceNumerator, divergence, directional, Fin.sum_univ_succ, Fin.sum_univ_zero, add_zero]\n'
                 f'  change D * ({div_sum}) - ({dir_sum}) = F^{q}*Q\n'
                 f'  have h := congrArg (fun p => C (1 / {sq} : ℂ) * p) e2_core\n'
                 f'  simp only [L2, R2, denC_smulN, denC_sub, denC_add, denC_mul, denC_pow, denC_poly, denC_pd, Nat.cast_ofNat] at h\n'
                 f'  rw [hF{q}]\n'
                 f'  simp only [D, Q, N0, N1, N2, N3, N4, N5]\n'
                 f'  linear_combination h - (toC tD * ({div_u}) - ({dir_u})) * hs\n')
    # homogeneity, nonvanishing
    L.append(f'theorem denominator_homogeneous : D.IsHomogeneous {d} := toC_isHomogeneous {d} tD (by decide +kernel)\n')
    for i in range(NV):
        L.append(f'theorem N{i}_homogeneous : N{i}.IsHomogeneous {degN} := toC_isHomogeneous {degN} tN{i} (by decide +kernel)\n')
    L.append(f'theorem numerator_homogeneous : ∀ i, (N i).IsHomogeneous {degN} := by\n  intro i; fin_cases i\n' +
             ''.join(f'  · exact N{i}_homogeneous\n' for i in range(NV)))
    nzpt = leangen.choose_nonzero_point(D)
    vec = '![' + ','.join(str(v) for v in nzpt) + ']'
    L.append(f'theorem denominator_ne_zero : D ≠ 0 := toC_ne_zero_of_evalPT {vec} tD (by decide +kernel)\n')
    # legality: template machinery on the explicit (unfolded) polynomials
    jj, pt, (bu, bv), fF, fD = leangen.choose_specialization(F, D)
    others = [i for i in range(NV) if i != jj]
    lc = fF[6]
    to_branches = []
    for i in range(NV):
        to_branches.append('(Polynomial.X)' if i == jj else f'(Polynomial.C (MvPolynomial.X {others.index(i)}))')
    to_fun = ''.join(f'if i.val = {i} then {to_branches[i]} else ' for i in range(5)) + to_branches[5]
    from_branches = [f'(MvPolynomial.X {i})' for i in others]
    from_fun = ''.join(f'if i.val = {t} then {from_branches[t]} else ' for t in range(4)) + from_branches[4]
    L.append(f'def toUni : Poly →+* Polynomial Coeff := MvPolynomial.eval₂Hom (Polynomial.C.comp MvPolynomial.C) (fun i => {to_fun})\n')
    L.append(f'def fromUni : Polynomial Coeff →+* Poly := Polynomial.eval₂RingHom (MvPolynomial.eval₂Hom MvPolynomial.C (fun i => {from_fun})) (MvPolynomial.X {jj})\n')
    for i in range(NV):
        L.append(f'@[simp] theorem toUni_X{i} : toUni (MvPolynomial.X {i}) = {to_branches[i]} := by simp only [toUni,MvPolynomial.eval₂Hom_X\']; rfl\n')
    L.append('@[simp] theorem toUni_C (c : ℂ) : toUni (MvPolynomial.C c) = Polynomial.C (MvPolynomial.C c) := by simp [toUni]\n')
    for t, i in enumerate(others):
        L.append(f'@[simp] theorem fromUni_w{t} : fromUni (Polynomial.C (MvPolynomial.X {t})) = (MvPolynomial.X {i}) := by simp only [fromUni,Polynomial.coe_eval₂RingHom,Polynomial.eval₂_C,MvPolynomial.eval₂Hom_X\']; rfl\n')
    L.append(f'@[simp] theorem fromUni_t : fromUni Polynomial.X = MvPolynomial.X {jj} := by simp [fromUni]\n')
    simp_from = ','.join(['fromUni_t'] + [f'fromUni_w{t}' for t in range(5)])
    L.append('theorem retract : fromUni.comp toUni = RingHom.id Poly := by\n'
             '  apply MvPolynomial.ringHom_ext\n'
             '  · intro c; simp [toUni,fromUni]\n'
             '  · intro i; fin_cases i <;> simp only [RingHom.comp_apply,RingHom.id_apply]\n' +
             ''.join(f'    · change fromUni (toUni (X {i})) = X {i}; simp only [toUni_X{i},map_add,{simp_from}] <;> ring\n' for i in range(NV)))
    by_deg = {}
    for e, c in F.items():
        k = e[jj]; e2 = list(e); e2[jj] = 0
        by_deg.setdefault(k, {})[tuple(e2)] = c
    funi_terms = [f'Polynomial.C ({leangen._coeff_poly_lean(by_deg[k], jj)}) * Polynomial.X^{k}' for k in sorted(by_deg, reverse=True)]
    L.append('def fUni : Polynomial Coeff := ' + ' +\n    '.join(funi_terms) + '\n')
    tos = ','.join(f'toUni_X{i}' for i in range(NV))
    # explicit forms of F and D for the univariate specialization (unfold the trees)
    L.append(f'theorem toUni_F : toUni F = fUni := by\n  simp only [F, tF, toC, fUni, map_add, map_mul, map_pow, toUni_C, {tos}]; polynomial\n')
    L.append('theorem fUni_degree : fUni.natDegree = 6 := by\n  unfold fUni; compute_degree!\n')
    L.append(f'theorem fUni_lc : fUni.leadingCoeff = MvPolynomial.C ({leangen._coef_str(lc)} : ℂ) := by\n'
             '  rw [Polynomial.leadingCoeff, fUni_degree]\n'
             '  simp only [fUni, Polynomial.coeff_add, Polynomial.coeff_C_mul_X_pow]\n'
             '  norm_num\n')
    L.append('def point : Fin 5 → ℂ := fun i => ' + ''.join(f'if i.val = {t} then ({pt[t]} : ℂ) else ' for t in range(4)) + f'({pt[4]} : ℂ)\n')
    L.append('def specialize : Coeff →+* ℂ := MvPolynomial.eval point\n')
    L.append(f'def fFiber : Polynomial ℂ := {leangen._uni_lean(fF)}\n')
    L.append(f'def dFiber : Polynomial ℂ := {leangen._uni_lean(fD)}\n')
    L.append(f'def bezoutF : Polynomial ℂ := {leangen._sympy_uni_lean(bu)}\n')
    L.append(f'def bezoutD : Polynomial ℂ := {leangen._sympy_uni_lean(bv)}\n')
    L.append('theorem fiber_F : (toUni F).map specialize = fFiber := by\n  rw [toUni_F]\n  norm_num [fUni,specialize,point,fFiber] <;> polynomial\n')
    L.append(f'theorem fiber_D : (toUni D).map specialize = dFiber := by\n  simp only [D,tD,toC,map_add,map_mul,map_pow,toUni_C,{tos}]\n  norm_num [specialize,point,dFiber] <;> polynomial\n')
    L.append('theorem fiber_coprime : IsCoprime fFiber dFiber := by\n  refine ⟨bezoutF, bezoutD, ?_⟩\n  unfold bezoutF bezoutD fFiber dFiber\n  polynomial\n')
    L.append('theorem denominator_legal : IsRelPrime F D := by\n'
             '  apply relPrime_of_retract toUni fromUni retract\n'
             '  apply relPrime_of_specialization specialize\n'
             '  · rw [toUni_F, fUni_lc]\n'
             f'    exact (isUnit_iff_ne_zero.mpr (by norm_num : ({leangen._coef_str(lc)} : ℂ) ≠ 0)).map MvPolynomial.C\n'
             '  · rw [fiber_F, fiber_D]\n'
             '    exact fiber_coprime.isRelPrime\n')
    L.append(f'theorem certificate : Verify F A {q} {d} D N Q where\n'
             '  denominator_ne_zero := denominator_ne_zero\n'
             '  denominator_homogeneous := denominator_homogeneous\n'
             '  numerator_homogeneous := numerator_homogeneous\n'
             '  denominator_legal := denominator_legal\n'
             '  identity := normalized_identity (by simpa using e1) e2\n')
    L.append(f'theorem solved : HasCertificate F A {q} := ⟨{d}, D, N, Q, certificate⟩\n')
    L.append(f'end\nend {ns}\n')
    open(os.path.join(out, 'Witness.lean'), 'w').write('\n'.join(L))
    # ---------------- Solution.lean: bind protected terms
    S = ['import HodgeWitness\nimport Proposal.Kron\nimport Proposal.Witness\n',
         f'set_option maxHeartbeats {heartbeats}\nset_option maxRecDepth {recdepth}\nset_option linter.all false\n',
         'noncomputable section\nnamespace Submission\nopen HodgeComputer.Kron\n']
    for nm, src in (('F', 'HodgeSubmission.F'), ('A', 'HodgeSubmission.A')):
        S.append(f'theorem h{nm} : {src} = {ns}.{nm} := by\n  apply MvPolynomial.funext; intro x\n'
                 f'  norm_num [{src}, {ns}.{nm}, {ns}.t{nm}, toC] <;> ring\n')
    def bind(nm, P, hwname, tname):
        if not P:
            S.append(f'theorem h{nm} : {hwname} = {ns}.{nm} := by\n  simp [{hwname}, {ns}.{nm}, {ns}.{tname}, toC]\n')
        else:
            S.append(f'noncomputable def w{nm} : W := {B.w_ast(wjson_str[nm])}\n')
            S.append(f'theorem h{nm}_rfl : {hwname} = denW w{nm} := rfl\n')
            S.append(f'theorem h{nm} : {hwname} = {ns}.{nm} := h{nm}_rfl.trans (matchW_sound w{nm} {ns}.{tname} (by decide +kernel))\n')
    wjson_str = {'D': wjson['D'], 'Q': wjson['Q']}
    for i in range(NV):
        wjson_str[f'N{i}'] = wjson['N'][i]
    bind('D', Di, 'HodgeWitness.D', 'tD')
    bind('Q', Qi, 'HodgeWitness.Q', 'tQ')
    for i in range(NV):
        bind(f'N{i}', Ni[i], f'HodgeWitness.N{i}', f'tN{i}')
    S.append(f'theorem hN : HodgeWitness.N = {ns}.N := by\n  funext i; fin_cases i\n' +
             ''.join(f'  · exact hN{i}\n' for i in range(NV)))
    S.append(f'theorem witness : HodgeComputer.Verify HodgeSubmission.F HodgeSubmission.A HodgeSubmission.q HodgeWitness.d HodgeWitness.D HodgeWitness.N HodgeWitness.Q := by\n'
             f'  simpa only [hF,hA,hD,hQ,hN,HodgeSubmission.q,HodgeWitness.d] using {ns}.certificate\n')
    S.append('theorem certificate : HodgeComputer.HasCertificate HodgeSubmission.F HodgeSubmission.A HodgeSubmission.q := ⟨HodgeWitness.d,HodgeWitness.D,HodgeWitness.N,HodgeWitness.Q,witness⟩\n')
    S.append('end Submission\nend\n')
    open(os.path.join(out, 'Solution.lean'), 'w').write('\n'.join(S))
    json.dump({'schemaVersion': 2, 'instanceId': instance_id, 'theorem': 'Submission.certificate',
               'witnessTheorem': 'Submission.witness'}, open(os.path.join(out, 'proof.json'), 'w'), indent=1)
    open(os.path.join(out, 'solution.md'), 'w').write(f'# Kronecker-checked witness for {instance_id} (capacity prototype)\n')
    meta = {'instance': instance_id, 'd': d, 'q': q, 'scale_witness': s, 'scale_FA': sFA, 'm': m, 'W': W,
            'terms': {'F': len(Fi), 'A': len(Ai), 'D': len(Di), 'N': [len(n) for n in Ni], 'Q': len(Qi)},
            'bytes': {f: os.path.getsize(os.path.join(out, f)) for f in os.listdir(out)}}
    json.dump(meta, open(os.path.join(out, 'meta.json'), 'w'), indent=1)
    return meta


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--instance', required=True)
    ap.add_argument('--witness')
    ap.add_argument('--out', required=True)
    ap.add_argument('--name', default='KronWitness')
    a = ap.parse_args()
    if a.witness:
        witness = json.load(open(a.witness))
    else:
        witness = json.load(open(f'/home/kejace/devel/sextics/data/{a.instance}/witness.json'))
    print(json.dumps(generate(a.instance, witness, a.out, a.name)))


if __name__ == '__main__':
    main()
