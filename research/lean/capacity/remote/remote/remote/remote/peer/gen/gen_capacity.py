#!/usr/bin/env python3
"""Synthetic, exactly-true certificate identities of adjustable size, emitted in the accepted Lean template.

Construction (q = pole order, d = deg D):
    M_i := sum_j d_j H_ij        with H antisymmetric of degree 6q-4   => div M = 0
    U_i := sum_j d_jF T_ij       with T antisymmetric of degree d-10   => sum_i d_iF U_i = 0 (Koszul)
    N_i := D*M_i + F^q*U_i,  A := q*sum_i d_iF*M_i,  Q := D*div U - sum_i d_iD*U_i
Then  q*sum_i d_iF N_i = A*D  and  D*div N - sum_i d_iD N_i = F^q*Q  hold exactly (T = 0 when d < 10, then Q = 0).
The identities have genuine certificate shape (LHS products of pderivs with N, RHS A*D and F^q*Q) and sizes are
controlled by |D| (k_D random monomials or dense), |H| (k_H), |T| (k_T) and the choice of F (sparse accepted lane
or the dense-rational sextic).

Usage:
  python3 gen_capacity.py --F m3 --q 1 --d 6 --kD 20 --kH 2 --kT 0 --parts e1 --name S1 [--seed 0] [--stats-only]
Writes src/Cap/<name>.lean and src/Cap/<name>.json (size statistics) under /tmp/hodge/lean/capacity.
"""
from __future__ import annotations
import argparse, json, random, sys, os
from fractions import Fraction
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core, leangen
from hodge.core import NV

CAP = '/tmp/hodge/lean/capacity'
FSRC = {'m3': 'm3-coupled', 'dense': 'dense-rational', 'hesse': 'hesse-triangle', 'contact': 'contact46006',
        'cycle': 'cycle-q2', 'klein': 'klein-rational'}


def rand_poly(rng, deg, k, lo=-9, hi=9, par=None):
    mons = core.monomials(deg, par)
    if k is None or k >= len(mons):
        chosen = mons
    else:
        chosen = rng.sample(mons, k)
    P = {}
    for e in chosen:
        c = 0
        while c == 0:
            c = rng.randint(lo, hi)
        P[e] = Fraction(c)
    return P


def antisym_family(rng, deg, k):
    """H_ij antisymmetric: dict (i,j)->poly for i<j (H_ji = -H_ij)."""
    H = {}
    if deg < 0 or k == 0:
        return H
    for i in range(NV):
        for j in range(i + 1, NV):
            P = rand_poly(rng, deg, k)
            if P:
                H[(i, j)] = P
    return H


def contract(H, G):
    """V_i = sum_j G_j * H_ij with H antisymmetric (G_j given polys)."""
    V = [{} for _ in range(NV)]
    for (i, j), P in H.items():
        V[i] = core.add(V[i], core.mul(G[j], P))
        V[j] = core.sub(V[j], core.mul(G[i], P))
    return V


def build(Fkey, q, d, kD, kH, kT, seed=0, dense_D=False):
    rng = random.Random(seed)
    F, _A, _q, _rec = core.load_instance(FSRC[Fkey])
    D = rand_poly(rng, d, None if dense_D else kD)
    # M_i = sum_j d_j H_ij, deg H = 6q-4
    H = antisym_family(rng, 6 * q - 4, kH)
    # d_j H_ij: contract with "derivative operator" -- do it directly
    M = [{} for _ in range(NV)]
    for (i, j), P in H.items():
        M[i] = core.add(M[i], core.pderiv(P, j))
        M[j] = core.sub(M[j], core.pderiv(P, i))
    dF = [core.pderiv(F, i) for i in range(NV)]
    T = antisym_family(rng, d - 10, kT) if d >= 10 else {}
    U = contract(T, dF)
    Fq = core.power(F, q)
    N = [core.add(core.mul(D, M[i]), core.mul(Fq, U[i])) for i in range(NV)]
    A = {}
    for i in range(NV):
        A = core.add(A, core.mul(dF[i], M[i]))
    A = core.scale(A, Fraction(q))
    divU = {}
    dirDU = {}
    dD = [core.pderiv(D, i) for i in range(NV)]
    for i in range(NV):
        divU = core.add(divU, core.pderiv(U[i], i))
        dirDU = core.add(dirDU, core.mul(dD[i], U[i]))
    Q = core.sub(core.mul(D, divU), dirDU)
    ok, info = core.verify_witness(F, A, q, D, N, Q)
    assert ok, info
    return F, A, D, N, Q, dF, dD


def stats(F, A, D, N, Q, q, dF, dD):
    Fq = core.power(F, q)
    e1_lhs = sum(len(dF[i]) * len(N[i]) for i in range(NV))
    e1_rhs = len(A) * len(D)
    e2_lhs = len(D) * sum(len(core.pderiv(N[i], i)) for i in range(NV)) + sum(len(dD[i]) * len(N[i]) for i in range(NV))
    e2_rhs = len(Fq) * len(Q)
    e1_col = len(core.mul(A, D))
    e2_col = len(core.mul(Fq, Q))
    return {'F': len(F), 'A': len(A), 'D': len(D), 'N': [len(n) for n in N], 'Q': len(Q), 'q': q,
            'e1_raw': e1_lhs + e1_rhs, 'e1_lhs_raw': e1_lhs, 'e1_rhs_raw': e1_rhs, 'e1_collected': e1_col,
            'e2_raw': e2_lhs + e2_rhs, 'e2_lhs_raw': e2_lhs, 'e2_rhs_raw': e2_rhs, 'e2_collected': e2_col,
            'defs_terms': len(F) + len(A) + len(D) + sum(len(n) for n in N) + len(Q)}


def nest(terms):
    return terms[0] if len(terms) == 1 else f'{terms[0]} + ({nest(terms[1:])})'


def module(name, F, A, D, N, Q, q, d, parts, heartbeats=2000000, recdepth=4000, async_elab=False):
    L = leangen
    ns = f'HodgeComputer.Cap.{name}'
    degN = d + 6 * q - 5
    out = ['import HodgeComputer.Verify\n']
    if not async_elab:
        out.append('set_option Elab.async false\n')
    out += [f'set_option maxHeartbeats {heartbeats}\n', f'set_option maxRecDepth {recdepth}\n',
            'set_option linter.unusedTactic false\nset_option linter.unusedSimpArgs false\n'
            'set_option linter.unnecessarySeqFocus false\nset_option linter.unreachableTactic false\n',
            'noncomputable section\n', f'namespace {ns}\n', 'open MvPolynomial HodgeTestbed\n']
    out.append(f'def F : Poly := {L._poly_lean(F)}\n')
    out.append(f'def A : Poly := {L._poly_lean(A)}\n')
    out.append(f'def D : Poly := {L._poly_lean(D)}\n')
    out.append(f'def Q : Poly := {L._poly_lean(Q)}\n')
    for i in range(NV):
        out.append(f'def N{i} : Poly := {L._poly_lean(N[i])}\n')
    out.append('def N : Fin 6 → Poly := fun i => ' + ''.join(f'if i.val = {i} then N{i} else ' for i in range(5)) + 'N5\n')
    out.append(L._deriv_lemmas())
    if 'e1' in parts:
        e1_sum = nest([f'pderiv ({i} : Fin 6) F * N{i}' for i in range(NV)])
        out.append(f'theorem e1 : ({q} : Poly) * directional F N = A*D := by\n'
                   f'  simp only [directional, Fin.sum_univ_succ, Fin.sum_univ_zero, add_zero]\n'
                   f'  change ({q}:Poly) * ({e1_sum}) = A*D\n'
                   f'  simp only [F,A,D,N0,N1,N2,N3,N4,N5,map_zero,map_add,pderiv_C,pderiv_mul,pderiv_pow,{L.DERIV_NAMES}]\n'
                   f'  apply MvPolynomial.funext; intro x\n'
                   f'  norm_num <;> ring\n')
    if 'e2' in parts:
        div_sum = nest([f'pderiv ({i} : Fin 6) N{i}' for i in range(NV)])
        dir_sum = nest([f'pderiv ({i} : Fin 6) D * N{i}' for i in range(NV)])
        out.append(f'theorem e2 : divergenceNumerator D N = F^{q}*Q := by\n'
                   f'  simp only [divergenceNumerator, divergence, directional, Fin.sum_univ_succ, Fin.sum_univ_zero, add_zero]\n'
                   f'  change D * ({div_sum}) - ({dir_sum}) = F^{q}*Q\n'
                   f'  simp only [F,D,Q,N0,N1,N2,N3,N4,N5,map_zero,map_add,pderiv_C,pderiv_mul,pderiv_pow,{L.DERIV_NAMES}]\n'
                   f'  apply MvPolynomial.funext; intro x\n'
                   f'  norm_num <;> ring\n')
    if 'homog' in parts:
        out.append(L._homog_proof(D, 'D', d).replace('theorem D_homogeneous', 'theorem denominator_homogeneous'))
        for i in range(NV):
            out.append(L._homog_proof(N[i], f'N{i}', degN))
        out.append(f'theorem numerator_homogeneous : ∀ i, (N i).IsHomogeneous {degN} := by\n  intro i; fin_cases i\n' +
                   ''.join(f'  · exact N{i}_homogeneous\n' for i in range(NV)))
    if 'legal' in parts:
        jj, pt, (bu, bv), fF, fD = L.choose_specialization(F, D)
        others = [i for i in range(NV) if i != jj]
        nzpt = L.choose_nonzero_point(D)
        lc = fF[6]
        to_branches = []
        for i in range(NV):
            if i == jj:
                to_branches.append('(Polynomial.X)')
            else:
                to_branches.append(f'(Polynomial.C (MvPolynomial.X {others.index(i)}))')
        to_fun = ''.join(f'if i.val = {i} then {to_branches[i]} else ' for i in range(5)) + to_branches[5]
        from_branches = [f'(MvPolynomial.X {i})' for i in others]
        from_fun = ''.join(f'if i.val = {t} then {from_branches[t]} else ' for t in range(4)) + from_branches[4]
        out.append(f'def toUni : Poly →+* Polynomial Coeff := MvPolynomial.eval₂Hom (Polynomial.C.comp MvPolynomial.C) (fun i => {to_fun})\n')
        out.append(f'def fromUni : Polynomial Coeff →+* Poly := Polynomial.eval₂RingHom (MvPolynomial.eval₂Hom MvPolynomial.C (fun i => {from_fun})) (MvPolynomial.X {jj})\n')
        for i in range(NV):
            out.append(f'@[simp] theorem toUni_X{i} : toUni (MvPolynomial.X {i}) = {to_branches[i]} := by simp only [toUni,MvPolynomial.eval₂Hom_X\']; rfl\n')
        out.append('@[simp] theorem toUni_C (c : ℂ) : toUni (MvPolynomial.C c) = Polynomial.C (MvPolynomial.C c) := by simp [toUni]\n')
        for t, i in enumerate(others):
            out.append(f'@[simp] theorem fromUni_w{t} : fromUni (Polynomial.C (MvPolynomial.X {t})) = (MvPolynomial.X {i}) := by simp only [fromUni,Polynomial.coe_eval₂RingHom,Polynomial.eval₂_C,MvPolynomial.eval₂Hom_X\']; rfl\n')
        out.append(f'@[simp] theorem fromUni_t : fromUni Polynomial.X = MvPolynomial.X {jj} := by simp [fromUni]\n')
        simp_from = ','.join(['fromUni_t'] + [f'fromUni_w{t}' for t in range(5)])
        out.append('theorem retract : fromUni.comp toUni = RingHom.id Poly := by\n'
                   '  apply MvPolynomial.ringHom_ext\n'
                   '  · intro c; simp [toUni,fromUni]\n'
                   '  · intro i; fin_cases i <;> simp only [RingHom.comp_apply,RingHom.id_apply]\n' +
                   ''.join(f'    · change fromUni (toUni (X {i})) = X {i}; simp only [toUni_X{i},map_add,{simp_from}] <;> ring\n' for i in range(NV)))
        by_deg = {}
        for e, c in F.items():
            k = e[jj]
            e2 = list(e); e2[jj] = 0
            by_deg.setdefault(k, {})[tuple(e2)] = c
        funi_terms = []
        for k in sorted(by_deg, reverse=True):
            funi_terms.append(f'Polynomial.C ({L._coeff_poly_lean(by_deg[k], jj)}) * Polynomial.X^{k}')
        out.append('def fUni : Polynomial Coeff := ' + ' +\n    '.join(funi_terms) + '\n')
        tos = ','.join(f'toUni_X{i}' for i in range(NV))
        out.append(f'theorem toUni_F : toUni F = fUni := by\n  simp only [F, fUni, map_add, map_mul, map_pow, toUni_C, {tos}]; polynomial\n')
        out.append('theorem fUni_degree : fUni.natDegree = 6 := by\n  unfold fUni; compute_degree!\n')
        out.append(f'theorem fUni_lc : fUni.leadingCoeff = MvPolynomial.C ({L._coef_str(lc)} : ℂ) := by\n'
                   '  rw [Polynomial.leadingCoeff, fUni_degree]\n'
                   '  simp only [fUni, Polynomial.coeff_add, Polynomial.coeff_C_mul_X_pow]\n'
                   '  norm_num\n')
        out.append('def point : Fin 5 → ℂ := fun i => ' + ''.join(f'if i.val = {t} then ({pt[t]} : ℂ) else ' for t in range(4)) + f'({pt[4]} : ℂ)\n')
        out.append('def specialize : Coeff →+* ℂ := MvPolynomial.eval point\n')
        out.append(f'def fFiber : Polynomial ℂ := {L._uni_lean(fF)}\n')
        out.append(f'def dFiber : Polynomial ℂ := {L._uni_lean(fD)}\n')
        out.append(f'def bezoutF : Polynomial ℂ := {L._sympy_uni_lean(bu)}\n')
        out.append(f'def bezoutD : Polynomial ℂ := {L._sympy_uni_lean(bv)}\n')
        out.append('theorem fiber_F : (toUni F).map specialize = fFiber := by\n  rw [toUni_F]\n  norm_num [fUni,specialize,point,fFiber] <;> polynomial\n')
        out.append(f'theorem fiber_D : (toUni D).map specialize = dFiber := by\n  simp only [D,map_add,map_mul,map_pow,toUni_C,{tos}]\n  norm_num [specialize,point,dFiber] <;> polynomial\n')
        out.append('theorem fiber_coprime : IsCoprime fFiber dFiber := by\n  refine ⟨bezoutF, bezoutD, ?_⟩\n  unfold bezoutF bezoutD fFiber dFiber\n  polynomial\n')
        out.append('theorem denominator_legal : IsRelPrime F D := by\n'
                   '  apply relPrime_of_retract toUni fromUni retract\n'
                   '  apply relPrime_of_specialization specialize\n'
                   '  · rw [toUni_F, fUni_lc]\n'
                   f'    exact (isUnit_iff_ne_zero.mpr (by norm_num : ({L._coef_str(lc)} : ℂ) ≠ 0)).map MvPolynomial.C\n'
                   '  · rw [fiber_F, fiber_D]\n'
                   '    exact fiber_coprime.isRelPrime\n')
        vec = '![' + ','.join(str(v) for v in nzpt) + ']'
        shows = ', '.join(f'show ({vec} : Fin 6 → ℂ) {i} = {nzpt[i]} from rfl' for i in range(NV))
        out.append('theorem denominator_ne_zero : D ≠ 0 := by\n'
                   '  intro h\n'
                   f'  have he := congrArg (MvPolynomial.eval ({vec} : Fin 6 → ℂ)) h\n'
                   f'  norm_num [D, {shows}] at he\n')
    if all(p in parts for p in ('e1', 'e2', 'homog', 'legal')):
        out.append(f'theorem certificate : Verify F A {q} {d} D N Q where\n'
                   '  denominator_ne_zero := denominator_ne_zero\n'
                   '  denominator_homogeneous := denominator_homogeneous\n'
                   '  numerator_homogeneous := numerator_homogeneous\n'
                   '  denominator_legal := denominator_legal\n'
                   '  identity := normalized_identity (by simpa using e1) e2\n')
        out.append(f'theorem solved : HasCertificate F A {q} := ⟨{d}, D, N, Q, certificate⟩\n')
    out.append(f'end {ns}\n')
    return '\n'.join(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--F', default='m3', choices=sorted(FSRC))
    ap.add_argument('--q', type=int, default=1)
    ap.add_argument('--d', type=int, default=6)
    ap.add_argument('--kD', type=int, default=20)
    ap.add_argument('--denseD', action='store_true')
    ap.add_argument('--kH', type=int, default=2)
    ap.add_argument('--kT', type=int, default=0)
    ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--parts', default='e1', help='comma list of e1,e2,homog,legal or full')
    ap.add_argument('--name', required=False)
    ap.add_argument('--heartbeats', type=int, default=2000000)
    ap.add_argument('--recdepth', type=int, default=4000)
    ap.add_argument('--async', dest='async_elab', action='store_true')
    ap.add_argument('--stats-only', action='store_true')
    ap.add_argument('--dump', help='also dump the polynomials as JSON (witness syntax) to this path')
    a = ap.parse_args()
    parts = ['e1', 'e2', 'homog', 'legal'] if a.parts == 'full' else a.parts.split(',')
    F, A, D, N, Q, dF, dD = build(a.F, a.q, a.d, a.kD, a.kH, a.kT, a.seed, a.denseD)
    st = stats(F, A, D, N, Q, a.q, dF, dD)
    st.update({'Fsrc': a.F, 'd': a.d, 'kD': a.kD, 'denseD': a.denseD, 'kH': a.kH, 'kT': a.kT, 'seed': a.seed, 'parts': parts})
    print(json.dumps(st))
    if a.stats_only:
        return
    name = a.name or f'{a.F}_q{a.q}_d{a.d}_D{"dense" if a.denseD else a.kD}_H{a.kH}_T{a.kT}_s{a.seed}_{"".join(parts)}'
    src = module(name, F, A, D, N, Q, a.q, a.d, parts, a.heartbeats, a.recdepth, a.async_elab)
    path = os.path.join(CAP, 'src', 'Cap', name + '.lean')
    with open(path, 'w') as f:
        f.write(src)
    st['bytes'] = len(src.encode())
    st['name'] = name
    with open(os.path.join(CAP, 'src', 'Cap', name + '.json'), 'w') as f:
        json.dump(st, f)
    if a.dump:
        with open(a.dump, 'w') as f:
            json.dump({'F': core.to_str(F), 'A': core.to_str(A), 'q': a.q, 'd': a.d, 'D': core.to_str(D),
                       'N': [core.to_str(n) for n in N], 'Q': core.to_str(Q)}, f)
    print(path, st['bytes'], 'bytes')


if __name__ == '__main__':
    main()
