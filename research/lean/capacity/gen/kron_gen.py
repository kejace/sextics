#!/usr/bin/env python3
"""Emit a Kronecker-checked module Cap/K_<name>.lean for synthetic data from gen_capacity.build.

Polynomials are emitted as balanced `PT` trees (leaves sorted by Kronecker position); the identities e1/e2 are
stated as in the accepted template and proved from `HodgeComputer.Kron.soundC` with `decide +kernel` side goals.

Usage: python3 kron_gen.py --F m3 --q 1 --d 6 --kD 40 --kH 1 --kT 0 --parts e1,e2 --name K1 [--seed 0] [--W 64]
"""
from __future__ import annotations
import argparse, json, math, os, sys
from fractions import Fraction
sys.path.insert(0, '/tmp/hodge/lean/capacity/gen'); sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core
from hodge.core import NV
import gen_capacity as G

CAP = '/tmp/hodge/lean/capacity'


def clear_denominators(P):
    """Return (integer poly, scale s) with s*P integral."""
    s = 1
    for c in P.values():
        s = s * Fraction(c).denominator // math.gcd(s, Fraction(c).denominator)
    return {e: int(Fraction(c) * s) for e, c in P.items()}, s


def pos(e, m):
    return e[0] + m * (e[1] + m * (e[2] + m * (e[3] + m * e[4])))


def pt_lean(P, m, deg=0):
    """Balanced tree over leaves sorted by Kronecker position (zero polynomial: one zero leaf of degree `deg`)."""
    items = sorted(((pos(e, m), e, c) for e, c in P.items() if c), key=lambda t: t[0])
    if not items:
        return f'.leaf ⟨{deg}, 0, 0, 0, 0, 0⟩ 0'
    def build(lo, hi):
        if hi - lo == 1:
            _, e, c = items[lo]
            return f'.leaf ⟨{e[0]}, {e[1]}, {e[2]}, {e[3]}, {e[4]}, {e[5]}⟩ ({c})'
        mid = (lo + hi) // 2
        return f'.node ({build(lo, mid)}) ({build(mid, hi)})'
    return build(0, len(items))


def l1(P):
    return sum(abs(int(c)) for c in P.values())


def nest(terms):
    return terms[0] if len(terms) == 1 else f'{terms[0]} + ({nest(terms[1:])})'


def nestE(op, terms):
    return terms[0] if len(terms) == 1 else f'(.{op} {terms[0]} {nestE(op, terms[1:])})'


def integralize(D, N, Q):
    """Scale the witness to integer coefficients: D, N by s, Q by s^2 (the identities are homogeneous of degree 2)."""
    s = 1
    for P in (D, Q, *N):
        for c in P.values():
            den = Fraction(c).denominator
            s = s * den // math.gcd(s, den)
    return (core.scale(D, Fraction(s)), [core.scale(n, Fraction(s)) for n in N], core.scale(Q, Fraction(s * s)), s)


def module(name, F, A, D, N, Q, q, d, parts, W=None, heartbeats=2000000, recdepth=4000):
    Fi, sF = clear_denominators(F)
    Ai, sA = clear_denominators(A)
    # common scale s for F and A: F = (1/s) F', A = (1/s) A' with integral F', A'
    s_ = sF * sA // math.gcd(sF, sA)
    Fi = {e: int(Fraction(c) * s_) for e, c in F.items()}
    Ai = {e: int(Fraction(c) * s_) for e, c in A.items()}
    for P in (D, Q, *N):
        for c in P.values():
            assert Fraction(c).denominator == 1
    Di = {e: int(c) for e, c in D.items()}; Qi = {e: int(c) for e, c in Q.items()}
    Ni = [{e: int(c) for e, c in n.items()} for n in N]
    degN = d + 6 * q - 5
    deg1 = d + 6 * q            # degree of e1 identity
    deg2 = 2 * d + 6 * q - 6    # degree of e2 identity
    m = max(deg1, deg2) + 1
    # ℓ1 bounds (match Kron.l1E)
    l1F = l1(Fi); l1A = l1(Ai); l1D = l1(Di); l1Q = l1(Qi); l1N = [l1(n) for n in Ni]
    dF = [core.pderiv(Fi, i) for i in range(NV)]; dD = [core.pderiv(Di, i) for i in range(NV)]
    dN = [core.pderiv(Ni[i], i) for i in range(NV)]
    b1 = q * sum(l1(dF[i]) * l1N[i] for i in range(NV)) + l1A * l1D
    b2 = (s_ ** q) * (l1D * sum(l1(dN[i]) for i in range(NV)) + sum(l1(dD[i]) * l1N[i] for i in range(NV))) + l1F ** q * l1Q
    bound = max(b1, b2)
    if W is None:
        W = max(bound.bit_length() + 2, 8)
    ns = f'HodgeComputer.Cap.{name}'
    L = ['import Reflect.Kron\n',
         f'set_option maxHeartbeats {heartbeats}\nset_option maxRecDepth {recdepth}\n'
         'set_option linter.all false\n',
         f'namespace {ns}\n', 'open MvPolynomial HodgeTestbed HodgeComputer.Kron\n']
    L.append(f'noncomputable def tF : PT := {pt_lean(Fi, m, 6)}\n')
    L.append(f'noncomputable def tA : PT := {pt_lean(Ai, m, 6 * q)}\n')
    L.append(f'noncomputable def tD : PT := {pt_lean(Di, m, d)}\n')
    L.append(f'noncomputable def tQ : PT := {pt_lean(Qi, m, 2 * d - 6)}\n')
    for i in range(NV):
        L.append(f'noncomputable def tN{i} : PT := {pt_lean(Ni[i], m, degN)}\n')
    L.append('noncomputable section\n')
    if s_ == 1:
        L.append('def F : Poly := toC tF\ndef A : Poly := toC tA\n')
    else:
        L.append(f'def F : Poly := C (1 / {s_} : ℂ) * toC tF\ndef A : Poly := C (1 / {s_} : ℂ) * toC tA\n')
    L.append('def D : Poly := toC tD\ndef Q : Poly := toC tQ\n')
    for i in range(NV):
        L.append(f'def N{i} : Poly := toC tN{i}\n')
    L.append('def N : Fin 6 → Poly := fun i => ' + ''.join(f'if i.val = {i} then N{i} else ' for i in range(5)) + 'N5\n')
    side = f'(by decide +kernel) (by decide +kernel) (by decide +kernel) (by decide +kernel) (by decide +kernel) (by decide +kernel)'
    if 'e1' in parts:
        lhs = f'(.smulN {q} {nestE("add", [f"(.mul (.pd {i} tF) (.poly tN{i}))" for i in range(NV)])})'
        rhs = '(.mul (.poly tA) (.poly tD))'
        L.append(f'noncomputable def L1 : E := {lhs}\nnoncomputable def R1 : E := {rhs}\n')
        L.append(f'theorem e1_core : denC L1 = denC R1 := soundN (m := {m}) (W := {W}) L1 R1 {side}\n')
        e1_sum = nest([f'pderiv ({i} : Fin 6) F * N{i}' for i in range(NV)])
        if s_ == 1:
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
                     f'  have h := congrArg (fun p => C (1 / {s_} : ℂ) * p) e1_core\n'
                     f'  simp only [L1, R1, denC_smulN, denC_add, denC_mul, denC_poly, denC_pd, Nat.cast_one, Nat.cast_ofNat] at h\n'
                     f'  simp only [F, A, D, N0, N1, N2, N3, N4, N5, pderiv_C_mul]\n'
                     f'  linear_combination h\n')
    if 'e2' in parts:
        div_r = nestE('add', [f'(.pd {i} tN{i})' for i in range(NV)])
        dir_r = nestE('add', [f'(.mul (.pd {i} tD) (.poly tN{i}))' for i in range(NV)])
        lhs = f'(.sub (.mul (.poly tD) {div_r}) {dir_r})'
        if s_ != 1:
            lhs = f'(.smulN {s_ ** q} {lhs})'
        rhs = f'(.mul (.pow (.poly tF) {q}) (.poly tQ))'
        L.append(f'noncomputable def L2 : E := {lhs}\nnoncomputable def R2 : E := {rhs}\n')
        L.append(f'theorem e2_core : denC L2 = denC R2 := soundN (m := {m}) (W := {W}) L2 R2 {side}\n')
        div_sum = nest([f'pderiv ({i} : Fin 6) N{i}' for i in range(NV)])
        dir_sum = nest([f'pderiv ({i} : Fin 6) D * N{i}' for i in range(NV)])
        if s_ == 1:
            L.append(f'theorem e2 : divergenceNumerator D N = F^{q}*Q := by\n'
                     f'  simp only [divergenceNumerator, divergence, directional, Fin.sum_univ_succ, Fin.sum_univ_zero, add_zero]\n'
                     f'  change D * ({div_sum}) - ({dir_sum}) = F^{q}*Q\n'
                     f'  have h := e2_core\n'
                     f'  simp only [L2, R2, denC_sub, denC_add, denC_mul, denC_pow, denC_poly, denC_pd] at h\n'
                     f'  simpa only [F, D, Q, N0, N1, N2, N3, N4, N5] using h\n')
        else:
            sq = s_ ** q
            div_u = nest([f'pderiv ({i} : Fin 6) (toC tN{i})' for i in range(NV)])
            dir_u = nest([f'pderiv ({i} : Fin 6) (toC tD) * toC tN{i}' for i in range(NV)])
            L.append(f'theorem hF{q} : F ^ {q} = C (1 / {sq} : ℂ) * toC tF ^ {q} := by\n'
                     f'  rw [F, mul_pow, ← C_pow]; norm_num\n')
            L.append(f'theorem hs : (C (1 / {sq} : ℂ) : Poly) * ({sq} : Poly) = 1 := by\n'
                     f'  rw [← map_ofNat C {sq}, ← C_mul]; norm_num\n')
            L.append(f'theorem e2 : divergenceNumerator D N = F^{q}*Q := by\n'
                     f'  simp only [divergenceNumerator, divergence, directional, Fin.sum_univ_succ, Fin.sum_univ_zero, add_zero]\n'
                     f'  change D * ({div_sum}) - ({dir_sum}) = F^{q}*Q\n'
                     f'  have h := congrArg (fun p => C (1 / {sq} : ℂ) * p) e2_core\n'
                     f'  simp only [L2, R2, denC_smulN, denC_sub, denC_add, denC_mul, denC_pow, denC_poly, denC_pd, Nat.cast_ofNat] at h\n'
                     f'  rw [hF{q}]\n'
                     f'  simp only [D, Q, N0, N1, N2, N3, N4, N5]\n'
                     f'  linear_combination h - (toC tD * ({div_u}) - ({dir_u})) * hs\n')
    if 'homog' in parts:
        L.append(f'theorem denominator_homogeneous : D.IsHomogeneous {d} := toC_isHomogeneous {d} tD (by decide +kernel)\n')
        for i in range(NV):
            L.append(f'theorem N{i}_homogeneous : N{i}.IsHomogeneous {degN} := toC_isHomogeneous {degN} tN{i} (by decide +kernel)\n')
        L.append(f'theorem numerator_homogeneous : ∀ i, (N i).IsHomogeneous {degN} := by\n  intro i; fin_cases i\n' +
                 ''.join(f'  · exact N{i}_homogeneous\n' for i in range(NV)))
        pt = G.leangen.choose_nonzero_point(D)
        vec = '![' + ','.join(str(v) for v in pt) + ']'
        L.append(f'theorem denominator_ne_zero : D ≠ 0 := toC_ne_zero_of_evalPT {vec} tD (by decide +kernel)\n')
    L.append(f'end\nend {ns}\n')
    meta = {'m': m, 'W': W, 'bound_bits': bound.bit_length(), 'deg1': deg1, 'deg2': deg2, 'scale': s_}
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
               'parts': parts, 'route': 'kron'})
    src, meta = module(a.name, F, A, D, N, Q, a.q, a.d, parts, a.W, a.heartbeats, a.recdepth)
    st.update(meta)
    path = os.path.join(CAP, 'src', 'Cap', a.name + '.lean')
    open(path, 'w').write(src)
    st['bytes'] = len(src.encode()); st['name'] = a.name
    json.dump(st, open(os.path.join(CAP, 'src', 'Cap', a.name + '.json'), 'w'))
    print(json.dumps(st))


if __name__ == '__main__':
    main()
