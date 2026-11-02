#!/usr/bin/env python3
"""Degree-sliced variant of the accepted template (q = 1): every explicit polynomial P is written as
P = Σ_j X0^j * P_j with X0-free slices P_j; the identities e1/e2 are proved slice by slice with the
template tactic (`MvPolynomial.funext; norm_num <;> ring`) and assembled with `linear_combination`
at the level of opaque slice atoms.  Same synthetic data as gen_capacity.build.

Usage: python3 slice_gen.py --F m3 --d 6 --kD 40 --kH 1 --parts e1,e2 --name SL_m3k40
"""
from __future__ import annotations
import argparse, json, os, sys
from fractions import Fraction
sys.path.insert(0, '/tmp/hodge/lean/capacity/gen'); sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core, leangen
from hodge.core import NV
import gen_capacity as G

CAP = '/tmp/hodge/lean/capacity'


def slices(P, deg):
    """list of X0-free slice polys indexed by X0-degree 0..deg"""
    out = [dict() for _ in range(deg + 1)]
    for e, c in P.items():
        e2 = (0,) + tuple(e[1:])
        out[e[0]][e2] = out[e[0]].get(e2, 0) + c
    return out


def nest(terms):
    return terms[0] if len(terms) == 1 else f'{terms[0]} + ({nest(terms[1:])})'


def module(name, F, A, D, N, Q, d, parts, heartbeats=2000000, recdepth=4000):
    q = 1
    degN = d + 1
    degs = {'F': 6, 'A': 6, 'D': d, 'Q': 2 * d - 6}
    polys = {'F': F, 'A': A, 'D': D, 'Q': Q}
    for i in range(NV):
        polys[f'N{i}'] = N[i]; degs[f'N{i}'] = degN
    sl = {k: slices(P, degs[k]) for k, P in polys.items()}
    ns = f'HodgeComputer.Cap.{name}'
    L = ['import HodgeComputer.Verify\n', 'set_option Elab.async false\n',
         f'set_option maxHeartbeats {heartbeats}\nset_option maxRecDepth {recdepth}\nset_option linter.all false\n',
         'noncomputable section\n', f'namespace {ns}\n', 'open MvPolynomial HodgeTestbed\n']
    # slice definitions and assembled polynomials
    for k, P in polys.items():
        for j, Pj in enumerate(sl[k]):
            L.append(f'def {k}_{j} : Poly := {leangen._poly_lean(Pj)}\n')
        L.append(f'def {k} : Poly := ' + ' + '.join(f'X 0 ^ {j} * {k}_{j}' for j in range(degs[k] + 1)) + '\n')
    L.append('def N : Fin 6 → Poly := fun i => ' + ''.join(f'if i.val = {i} then N{i} else ' for i in range(5)) + 'N5\n')
    L.append(leangen._deriv_lemmas())
    # X0-freeness of slices
    for k in polys:
        for j in range(degs[k] + 1):
            L.append(f'theorem h0_{k}_{j} : pderiv (0 : Fin 6) {k}_{j} = 0 := by\n'
                     f'  simp only [{k}_{j}, map_zero, map_add, pderiv_C, pderiv_mul, pderiv_pow, {leangen.DERIV_NAMES}, '
                     f'mul_zero, zero_mul, add_zero, zero_add]\n')
    tactic = (f'  simp only [{{defs}}, map_zero, map_add, pderiv_C, pderiv_mul, pderiv_pow, {leangen.DERIV_NAMES}]\n'
              f'  apply MvPolynomial.funext; intro x\n  norm_num <;> ring\n')
    def pd(i, k, j):  # symbolic derivative of a slice (as it appears after simp of pderiv on the assembled poly)
        return f'pderiv ({i} : Fin 6) {k}_{j}'
    if 'e1' in parts:
        # slice m of Σ_i pderiv i F * N_i : i=0: Σ_{j+k=m} (j+1) F_{j+1} N0_k ; i≥1: Σ_{j+k=m} pderiv i F_j * Ni_k
        deg1 = 6 + degN  # symbolic X0-degree bound (top slices are identically zero but must be stated)
        slice_names = []
        for m in range(deg1 + 1):
            lhs_terms, defs = [], set()
            for j in range(0, 7):
                k = m - j
                if 0 <= k <= degN:
                    if j + 1 <= 6:
                        lhs_terms.append(f'({j + 1} : Poly) * F_{j + 1} * N0_{k}'); defs |= {f'F_{j + 1}', f'N0_{k}'}
                    for i in range(1, NV):
                        lhs_terms.append(f'{pd(i, "F", j)} * N{i}_{k}'); defs |= {f'F_{j}', f'N{i}_{k}'}
            rhs_terms = []
            for j in range(0, 7):
                k = m - j
                if 0 <= k <= d:
                    rhs_terms.append(f'A_{j} * D_{k}'); defs |= {f'A_{j}', f'D_{k}'}
            lhs = ' + '.join(lhs_terms) if lhs_terms else '0'
            rhs = ' + '.join(rhs_terms) if rhs_terms else '0'
            L.append(f'theorem e1_slice_{m} : {lhs} = {rhs} := by\n' + tactic.replace('{defs}', ','.join(sorted(defs))))
            slice_names.append(m)
        e1_sum = nest([f'pderiv ({i} : Fin 6) F * N{i}' for i in range(NV)])
        h0s = ','.join(f'h0_F_{j}' for j in range(7))
        lc = ' + '.join(f'X 0 ^ {m} * e1_slice_{m}' for m in slice_names)
        L.append(f'theorem e1 : (1 : Poly) * directional F N = A*D := by\n'
                 f'  simp only [directional, Fin.sum_univ_succ, Fin.sum_univ_zero, add_zero]\n'
                 f'  change (1:Poly) * ({e1_sum}) = A*D\n'
                 f'  simp only [F, A, D, N0, N1, N2, N3, N4, N5, map_add, pderiv_mul, pderiv_pow, {h0s}, '
                 f'derivative_0_0, derivative_1_0, derivative_2_0, derivative_3_0, derivative_4_0, derivative_5_0]\n'
                 f'  linear_combination {lc}\n')
    if 'e2' in parts:
        deg2 = d + degN  # symbolic X0-degree bound
        slice_names = []
        for m in range(deg2 + 1):
            lhs_terms, defs = [], set()
            # D * div N:  i=0: Σ_{a+k=m} (k+1) D_a N0_{k+1} ; i≥1: Σ_{a+k=m} D_a * pderiv i Ni_k
            for a in range(0, d + 1):
                k = m - a
                if 0 <= k <= degN:
                    if k + 1 <= degN:
                        lhs_terms.append(f'({k + 1} : Poly) * D_{a} * N0_{k + 1}'); defs |= {f'D_{a}', f'N0_{k + 1}'}
                    for i in range(1, NV):
                        lhs_terms.append(f'D_{a} * {pd(i, f"N{i}", k)}'); defs |= {f'D_{a}', f'N{i}_{k}'}
            # - Σ_i pderiv i D * N_i : i=0: Σ_{a+k=m} (a+1) D_{a+1} N0_k ; i≥1: pderiv i D_a * Ni_k
            neg_terms = []
            for a in range(0, d + 1):
                k = m - a
                if 0 <= k <= degN:
                    if a + 1 <= d:
                        neg_terms.append(f'({a + 1} : Poly) * D_{a + 1} * N0_{k}'); defs |= {f'D_{a + 1}', f'N0_{k}'}
                    for i in range(1, NV):
                        neg_terms.append(f'{pd(i, "D", a)} * N{i}_{k}'); defs |= {f'D_{a}', f'N{i}_{k}'}
            rhs_terms = []
            for b in range(0, 7):
                c = m - b
                if 0 <= c <= 2 * d - 6:
                    rhs_terms.append(f'F_{b} * Q_{c}'); defs |= {f'F_{b}', f'Q_{c}'}
            lhs = (' + '.join(lhs_terms) if lhs_terms else '0') + ' - (' + (' + '.join(neg_terms) if neg_terms else '0') + ')'
            rhs = ' + '.join(rhs_terms) if rhs_terms else '0'
            L.append(f'theorem e2_slice_{m} : {lhs} = {rhs} := by\n' + tactic.replace('{defs}', ','.join(sorted(defs))))
            slice_names.append(m)
        div_sum = nest([f'pderiv ({i} : Fin 6) N{i}' for i in range(NV)])
        dir_sum = nest([f'pderiv ({i} : Fin 6) D * N{i}' for i in range(NV)])
        h0s = ','.join([f'h0_D_{a}' for a in range(d + 1)] + [f'h0_N{i}_{k}' for i in range(NV) for k in range(degN + 1)])
        lc = ' + '.join(f'X 0 ^ {m} * e2_slice_{m}' for m in slice_names)
        L.append(f'theorem e2 : divergenceNumerator D N = F^1*Q := by\n'
                 f'  simp only [divergenceNumerator, divergence, directional, Fin.sum_univ_succ, Fin.sum_univ_zero, add_zero]\n'
                 f'  change D * ({div_sum}) - ({dir_sum}) = F^1*Q\n'
                 f'  simp only [F, D, Q, N0, N1, N2, N3, N4, N5, map_add, pderiv_mul, pderiv_pow, {h0s}, '
                 f'derivative_0_0, derivative_1_0, derivative_2_0, derivative_3_0, derivative_4_0, derivative_5_0]\n'
                 f'  linear_combination {lc}\n')
    L.append(f'end {ns}\n')
    return '\n'.join(L)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--F', default='m3', choices=sorted(G.FSRC))
    ap.add_argument('--d', type=int, default=6)
    ap.add_argument('--kD', type=int, default=20)
    ap.add_argument('--denseD', action='store_true')
    ap.add_argument('--kH', type=int, default=1)
    ap.add_argument('--kT', type=int, default=0)
    ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--parts', default='e1,e2')
    ap.add_argument('--name', required=True)
    ap.add_argument('--heartbeats', type=int, default=2000000)
    ap.add_argument('--recdepth', type=int, default=4000)
    a = ap.parse_args()
    parts = a.parts.split(',')
    F, A, D, N, Q, dF, dD = G.build(a.F, 1, a.d, a.kD, a.kH, a.kT, a.seed, a.denseD)
    st = G.stats(F, A, D, N, Q, 1, dF, dD)
    st.update({'Fsrc': a.F, 'd': a.d, 'kD': a.kD, 'denseD': a.denseD, 'kH': a.kH, 'kT': a.kT, 'seed': a.seed,
               'parts': parts, 'route': 'slice'})
    src = module(a.name, F, A, D, N, Q, a.d, parts, a.heartbeats, a.recdepth)
    path = os.path.join(CAP, 'src', 'Cap', a.name + '.lean')
    open(path, 'w').write(src)
    st['bytes'] = len(src.encode()); st['name'] = a.name
    json.dump(st, open(os.path.join(CAP, 'src', 'Cap', a.name + '.json'), 'w'))
    print(json.dumps(st))


if __name__ == '__main__':
    main()
