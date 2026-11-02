#!/usr/bin/env python3
"""Kronecker-checked modules with slicing by the degree in X4 (`Kron.sound4`): every slice identity is an
identity between X4-free polynomials and is checked with a Kronecker point on four variables (m^4 positions), so
kernel memory per declaration shrinks by a factor ~m.  Same synthetic data as gen_capacity.build.

Usage: python3 kron_slice_gen.py --F m3 --q 1 --d 10 --kD 300 --kH 2 --kT 2 --parts e1,e2,homog --name KS_m3k300
"""
from __future__ import annotations
import argparse, json, math, os, sys
from fractions import Fraction
sys.path.insert(0, '/tmp/hodge/lean/capacity/gen'); sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core
from hodge.core import NV
import gen_capacity as G
import kron_gen as K

CAP = '/tmp/hodge/lean/capacity'

SIDE8 = ' '.join(['(by decide +kernel)'] * 8)
PD4_SIMP = ('pderiv_X4_pow_zero, pderiv_X4_pow_one, pderiv_X4_pow_two, pderiv_X4_pow_three, pderiv_X4_pow_four, '
            'pderiv_X4_pow_five, pderiv_four_toC_sliceK')


def nest(terms):
    return terms[0] if len(terms) == 1 else f'{terms[0]} + ({nest(terms[1:])})'


def nestE(op, terms):
    return terms[0] if len(terms) == 1 else f'(.{op} {terms[0]} {nestE(op, terms[1:])})'


def maxe4(P):
    return max((e[4] for e in P), default=0)


def slice_nonempty(P, k):
    return any(e[4] == k for e in P)


def l1(P):
    return sum(abs(int(c)) for c in P.values())


def integralize(D, N, Q):
    """Scale the witness to integer coefficients: D, N by s, Q by s^2 (the identities are homogeneous of degree 2)."""
    s = 1
    for P in (D, Q, *N):
        for c in P.values():
            den = Fraction(c).denominator
            s = s * den // math.gcd(s, den)
    return (core.scale(D, Fraction(s)), [core.scale(n, Fraction(s)) for n in N], core.scale(Q, Fraction(s * s)), s)


def module(name, F, A, D, N, Q, q, d, parts, W=None, heartbeats=2000000, recdepth=4000, ns=None, import_line='import Reflect.Kron', extra=''):
    Fi0, sF = K.clear_denominators(F); Ai0, sA = K.clear_denominators(A)
    s_ = sF * sA // math.gcd(sF, sA)
    Fi = {e: int(Fraction(c) * s_) for e, c in F.items()}; Ai = {e: int(Fraction(c) * s_) for e, c in A.items()}
    Di = {e: int(c) for e, c in D.items()}; Qi = {e: int(c) for e, c in Q.items()}
    Ni = [{e: int(c) for e, c in n.items()} for n in N]
    degN = d + 6 * q - 5
    deg1 = d + 6 * q; deg2 = 2 * d + 6 * q - 6
    m = max(deg1, deg2) + 1            # position base for the trees (sorted order) and for sound4
    polys = {'F': Fi, 'A': Ai, 'D': Di, 'Q': Qi}
    degs = {'F': 6, 'A': 6 * q, 'D': d, 'Q': 2 * d - 6}
    for i in range(NV):
        polys[f'N{i}'] = Ni[i]; degs[f'N{i}'] = degN
    Kmax = {k: maxe4(P) for k, P in polys.items()}
    # ℓ1 bound (same as kron_gen; slices only lower it)
    dF = [core.pderiv(Fi, i) for i in range(NV)]; dD = [core.pderiv(Di, i) for i in range(NV)]
    dN = [core.pderiv(Ni[i], i) for i in range(NV)]
    b1 = q * sum(l1(dF[i]) * l1(Ni[i]) for i in range(NV)) + l1(Ai) * l1(Di)
    b2 = (s_ ** q) * (l1(Di) * sum(l1(dN[i]) for i in range(NV)) + sum(l1(dD[i]) * l1(Ni[i]) for i in range(NV))) + l1(Fi) ** q * l1(Qi)
    bound = max(b1, b2)
    if W is None:
        W = max(bound.bit_length() + 2, 8)
    ns = ns or f'HodgeComputer.Cap.{name}'
    L = [import_line + '\n', f'set_option maxHeartbeats {heartbeats}\nset_option maxRecDepth {recdepth}\n',
         'set_option linter.all false\n', f'namespace {ns}\n', 'open MvPolynomial HodgeTestbed HodgeComputer.Kron\n']
    for k, P in polys.items():
        L.append(f'noncomputable def t{k} : PT := {K.pt_lean(P, m, degs[k])}\n')
    L.append('noncomputable section\n')
    if s_ == 1:
        L.append('def F : Poly := toC tF\ndef A : Poly := toC tA\n')
    else:
        L.append(f'def F : Poly := C (1 / {s_} : ℂ) * toC tF\ndef A : Poly := C (1 / {s_} : ℂ) * toC tA\n')
    L.append('def D : Poly := toC tD\ndef Q : Poly := toC tQ\n')
    for i in range(NV):
        L.append(f'def N{i} : Poly := toC tN{i}\n')
    L.append('def N : Fin 6 → Poly := fun i => ' + ''.join(f'if i.val = {i} then N{i} else ' for i in range(5)) + 'N5\n')
    # slice bridges
    for k in polys:
        L.append(f'theorem s{k} : toC t{k} = ∑ j ∈ Finset.range ({Kmax[k]} + 1), X 4 ^ j * toC (sliceK j t{k}) := '
                 f'toC_eq_sum_sliceK t{k} {Kmax[k]} (by decide +kernel)\n')
    bridges = ', '.join(f's{k}' for k in polys)
    simp_expand = (f'  simp only [Finset.sum_range_succ, Finset.sum_range_zero, zero_add, map_add, map_sum, pderiv_mul, '
                   f'{PD4_SIMP}, mul_zero, zero_mul, add_zero, Nat.cast_ofNat, Nat.cast_one, Nat.cast_zero]\n')

    def sl(k, a):   # E for slice a of polynomial k (always the slice: atoms must match the expanded goal)
        return f'(.poly (sliceK {a} t{k}))'

    def pd(i, k, a):  # E for pderiv i of slice a of polynomial k (i ≠ 4)
        return f'(.pd {i} (sliceK {a} t{k}))'

    def zero_of_deg(dg):   # syntactically-zero expression (0 • leaf) so that `ring` sees 0
        return f'(.smulN 0 (.poly (.leaf ⟨{max(dg, 0)}, 0, 0, 0, 0, 0⟩ 0)))'

    def sumE(terms, dg):
        return nestE('add', terms) if terms else zero_of_deg(dg)

    slice_thms = []
    if 'e1' in parts:
        M1 = Kmax['F'] + max(Kmax[f'N{i}'] for i in range(NV))
        M1 = max(M1, Kmax['A'] + Kmax['D'])
        for mm in range(M1 + 1):
            lhs_terms = []
            for i in range(NV):
                for a in range(0, Kmax['F'] + 1):
                    b = mm - a
                    if b < 0 or b > Kmax[f'N{i}']:
                        continue
                    if i == 4:
                        if a + 1 <= Kmax['F']:
                            lhs_terms.append(f'(.smulN {a + 1} (.mul {sl("F", a + 1)} {sl(f"N{i}", b)}))')
                    else:
                        lhs_terms.append(f'(.mul {pd(i, "F", a)} {sl(f"N{i}", b)})')
            rhs_terms = []
            for a in range(0, Kmax['A'] + 1):
                b = mm - a
                if 0 <= b <= Kmax['D']:
                    rhs_terms.append(f'(.mul {sl("A", a)} {sl("D", b)})')
            dg = deg1 - mm
            L.append(f'noncomputable def L1_{mm} : E := (.smulN {q} {sumE(lhs_terms, dg)})\nnoncomputable def R1_{mm} : E := {sumE(rhs_terms, dg)}\n')
            L.append(f'theorem e1s_{mm} : denC L1_{mm} = denC R1_{mm} := sound4 (m := {m}) (W := {W}) L1_{mm} R1_{mm} {SIDE8}\n')
            slice_thms.append(f'X 4 ^ {mm} * e1s_{mm}')
        e1_int_lhs = f'({q} : Poly) * ({nest([f"pderiv ({i} : Fin 6) (toC tF) * toC tN{i}" for i in range(NV)])})'
        e1_sum = nest([f'pderiv ({i} : Fin 6) F * N{i}' for i in range(NV)])
        slice_simp = ', '.join([f'L1_{mm}, R1_{mm}' for mm in range(M1 + 1)])
        haves = ''.join(f'  have h{mm} := e1s_{mm}\n' for mm in range(M1 + 1))
        hnames = ' '.join(f'h{mm}' for mm in range(M1 + 1))
        L.append(f'theorem e1_int : {e1_int_lhs} = toC tA * toC tD := by\n'
                 f'  simp only [{bridges}]\n' + simp_expand + haves +
                 f'  simp only [{slice_simp}, denC_smulN, denC_add, denC_mul, denC_poly, denC_pd, Nat.cast_ofNat, Nat.cast_one, Nat.cast_zero, zero_mul] at {hnames}\n'
                 f'  linear_combination {" + ".join(f"X 4 ^ {mm} * h{mm}" for mm in range(M1 + 1))}\n')
        if s_ == 1:
            L.append(f'theorem e1 : ({q} : Poly) * directional F N = A*D := by\n'
                     f'  simp only [directional, Fin.sum_univ_succ, Fin.sum_univ_zero, add_zero]\n'
                     f'  change ({q}:Poly) * ({e1_sum}) = A*D\n'
                     f'  simp only [F, A, D, N0, N1, N2, N3, N4, N5]\n'
                     f'  exact e1_int\n')
        else:
            L.append(f'theorem e1 : ({q} : Poly) * directional F N = A*D := by\n'
                     f'  simp only [directional, Fin.sum_univ_succ, Fin.sum_univ_zero, add_zero]\n'
                     f'  change ({q}:Poly) * ({e1_sum}) = A*D\n'
                     f'  simp only [F, A, D, N0, N1, N2, N3, N4, N5, pderiv_C_mul]\n'
                     f'  linear_combination (C (1 / {s_} : ℂ)) * e1_int\n')
    slice_thms2 = []
    if 'e2' in parts:
        KN = max(Kmax[f'N{i}'] for i in range(NV))
        M2 = max(Kmax['D'] + KN, q * Kmax['F'] + Kmax['Q'])
        for mm in range(M2 + 1):
            pos_terms, neg_terms = [], []
            # D * div N : i≠4: D_a * pd_i N_i_b ; i=4: (b+1) D_a * N4_{b+1}
            for a in range(0, Kmax['D'] + 1):
                b = mm - a
                if b < 0:
                    continue
                for i in range(NV):
                    if i == 4:
                        if b + 1 <= Kmax['N4']:
                            pos_terms.append(f'(.smulN {b + 1} (.mul {sl("D", a)} {sl("N4", b + 1)}))')
                    else:
                        if b <= Kmax[f'N{i}']:
                            pos_terms.append(f'(.mul {sl("D", a)} {pd(i, f"N{i}", b)})')
            # Σ pderiv i D * N_i : i≠4: pd_i D_a * N_i_b ; i=4: (a+1) D_{a+1} * N4_b
            for a in range(0, Kmax['D'] + 1):
                b = mm - a
                if b < 0:
                    continue
                for i in range(NV):
                    if i == 4:
                        if a + 1 <= Kmax['D'] and b <= Kmax['N4']:
                            neg_terms.append(f'(.smulN {a + 1} (.mul {sl("D", a + 1)} {sl("N4", b)}))')
                    else:
                        if b <= Kmax[f'N{i}']:
                            neg_terms.append(f'(.mul {pd(i, "D", a)} {sl(f"N{i}", b)})')
            rhs_terms = []
            if q == 1:
                for a in range(0, Kmax['F'] + 1):
                    c = mm - a
                    if 0 <= c <= Kmax['Q']:
                        rhs_terms.append(f'(.mul {sl("F", a)} {sl("Q", c)})')
            else:
                for a1 in range(0, Kmax['F'] + 1):
                    for a2 in range(0, Kmax['F'] + 1):
                        c = mm - a1 - a2
                        if 0 <= c <= Kmax['Q']:
                            rhs_terms.append(f'(.mul (.mul {sl("F", a1)} {sl("F", a2)}) {sl("Q", c)})')
            dg = deg2 - mm
            lhs = f'(.sub {sumE(pos_terms, dg)} {sumE(neg_terms, dg)})'
            if s_ != 1:
                lhs = f'(.smulN {s_ ** q} {lhs})'
            L.append(f'noncomputable def L2_{mm} : E := {lhs}\nnoncomputable def R2_{mm} : E := {sumE(rhs_terms, dg)}\n')
            L.append(f'theorem e2s_{mm} : denC L2_{mm} = denC R2_{mm} := sound4 (m := {m}) (W := {W}) L2_{mm} R2_{mm} {SIDE8}\n')
            slice_thms2.append(f'X 4 ^ {mm} * e2s_{mm}')
        div_u = nest([f'pderiv ({i} : Fin 6) (toC tN{i})' for i in range(NV)])
        dir_u = nest([f'pderiv ({i} : Fin 6) (toC tD) * toC tN{i}' for i in range(NV)])
        sq = s_ ** q
        scal = f'({sq} : Poly) * ' if s_ != 1 else ''
        slice_simp2 = ', '.join([f'L2_{mm}, R2_{mm}' for mm in range(M2 + 1)])
        haves2 = ''.join(f'  have g{mm} := e2s_{mm}\n' for mm in range(M2 + 1))
        gnames = ' '.join(f'g{mm}' for mm in range(M2 + 1))
        L.append(f'theorem e2_int : {scal}(toC tD * ({div_u}) - ({dir_u})) = toC tF ^ {q} * toC tQ := by\n'
                 f'  simp only [{bridges}]\n' + simp_expand + haves2 +
                 f'  simp only [{slice_simp2}, denC_smulN, denC_sub, denC_add, denC_mul, denC_pow, denC_poly, denC_pd, Nat.cast_ofNat, Nat.cast_one, Nat.cast_zero, zero_mul] at {gnames}\n'
                 f'  linear_combination {" + ".join(f"X 4 ^ {mm} * g{mm}" for mm in range(M2 + 1))}\n')
        div_sum = nest([f'pderiv ({i} : Fin 6) N{i}' for i in range(NV)])
        dir_sum = nest([f'pderiv ({i} : Fin 6) D * N{i}' for i in range(NV)])
        if s_ == 1:
            L.append(f'theorem e2 : divergenceNumerator D N = F^{q}*Q := by\n'
                     f'  simp only [divergenceNumerator, divergence, directional, Fin.sum_univ_succ, Fin.sum_univ_zero, add_zero]\n'
                     f'  change D * ({div_sum}) - ({dir_sum}) = F^{q}*Q\n'
                     f'  simp only [F, D, Q, N0, N1, N2, N3, N4, N5]\n'
                     f'  exact e2_int\n')
        else:
            L.append(f'theorem hF{q} : F ^ {q} = C (1 / {sq} : ℂ) * toC tF ^ {q} := by\n  rw [F, mul_pow, ← C_pow]; norm_num\n')
            L.append(f'theorem hs : (C (1 / {sq} : ℂ) : Poly) * ({sq} : Poly) = 1 := by\n  rw [← map_ofNat C {sq}, ← C_mul]; norm_num\n')
            L.append(f'theorem e2 : divergenceNumerator D N = F^{q}*Q := by\n'
                     f'  simp only [divergenceNumerator, divergence, directional, Fin.sum_univ_succ, Fin.sum_univ_zero, add_zero]\n'
                     f'  change D * ({div_sum}) - ({dir_sum}) = F^{q}*Q\n'
                     f'  simp only [D, Q, N0, N1, N2, N3, N4, N5]\n'
                     f'  rw [hF{q}]\n'
                     f'  have h := congrArg (fun p => C (1 / {sq} : ℂ) * p) e2_int\n'
                     f'  simp only [← mul_assoc, hs, one_mul] at h\n'
                     f'  rw [h]; ring\n')
    if 'homog' in parts:
        L.append(f'theorem denominator_homogeneous : D.IsHomogeneous {d} := toC_isHomogeneous {d} tD (by decide +kernel)\n')
        for i in range(NV):
            L.append(f'theorem N{i}_homogeneous : N{i}.IsHomogeneous {degN} := toC_isHomogeneous {degN} tN{i} (by decide +kernel)\n')
        L.append(f'theorem numerator_homogeneous : ∀ i, (N i).IsHomogeneous {degN} := by\n  intro i; fin_cases i\n' +
                 ''.join(f'  · exact N{i}_homogeneous\n' for i in range(NV)))
        pt = G.leangen.choose_nonzero_point(D)
        vec = '![' + ','.join(str(v) for v in pt) + ']'
        L.append(f'theorem denominator_ne_zero : D ≠ 0 := toC_ne_zero_of_evalPT {vec} tD (by decide +kernel)\n')
    if extra:
        L.append(extra)
    L.append(f'end\nend {ns}\n')
    meta = {'m': m, 'W': W, 'bound_bits': bound.bit_length(), 'deg1': deg1, 'deg2': deg2, 'scale': s_,
            'slices_e1': (M1 + 1) if 'e1' in parts else 0, 'slices_e2': (M2 + 1) if 'e2' in parts else 0, 'Kmax': Kmax}
    return '\n'.join(L), meta


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--F', default='m3', choices=sorted(G.FSRC))
    ap.add_argument('--q', type=int, default=1)
    ap.add_argument('--d', type=int, default=6)
    ap.add_argument('--kD', type=int, default=20)
    ap.add_argument('--denseD', action='store_true')
    ap.add_argument('--kH', type=int, default=1)
    ap.add_argument('--kT', type=int, default=0)
    ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--parts', default='e1,e2')
    ap.add_argument('--name', required=True)
    ap.add_argument('--W', type=int)
    ap.add_argument('--heartbeats', type=int, default=2000000)
    ap.add_argument('--recdepth', type=int, default=4000)
    a = ap.parse_args()
    parts = ['e1', 'e2', 'homog'] if a.parts == 'full' else a.parts.split(',')
    F, A, D, N, Q, dF, dD = G.build(a.F, a.q, a.d, a.kD, a.kH, a.kT, a.seed, a.denseD)
    D, N, Q, sW = integralize(D, N, Q)
    st = G.stats(F, A, D, N, Q, a.q, dF, dD)
    st['scale_witness'] = sW
    st.update({'Fsrc': a.F, 'd': a.d, 'kD': a.kD, 'denseD': a.denseD, 'kH': a.kH, 'kT': a.kT, 'seed': a.seed,
               'parts': parts, 'route': 'kron-slice'})
    src, meta = module(a.name, F, A, D, N, Q, a.q, a.d, parts, a.W, a.heartbeats, a.recdepth)
    st.update(meta)
    path = os.path.join(CAP, 'src', 'Cap', a.name + '.lean')
    open(path, 'w').write(src)
    st['bytes'] = len(src.encode()); st['name'] = a.name
    json.dump(st, open(os.path.join(CAP, 'src', 'Cap', a.name + '.json'), 'w'))
    print(json.dumps(st))


if __name__ == '__main__':
    main()
