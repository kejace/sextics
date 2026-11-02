"""Exact (rational) fixed-denominator test in the Koszul normalization, (U,V) frame, for D in J(G).

Normalized equations (thread REPORT.md / DEGREE_SEVEN.md):  M(G) = T D,  B_D(M) := D div M - M(D) in (G).
With D = sum_i P_i a_i (a_i given, degree d-5) every solution of the first equation is
   M = T a + S_q,   S_q = sum_{j<k} q_jk (P_k e_j - P_j e_k),  deg q_jk = d-4   (Koszul exactness; no parity ansatz).
So a certificate exists iff   B_D(T a) + B_D(S_q) = G * Q   is solvable in (q, Q), deg Q = 2d-6 -- a LINEAR system.
If D and T are homogeneous for the mu6 grading (V-degree mod 6) the system splits into character blocks and only the
blocks hit by B_D(Ta) matter.  We compute rank and augmented rank per block mod p AND over Q (exact fmpz/fmpq).
An exact rational rank gap (rank_Q < aug_rank_Q) is a rigorous exclusion of D for the given target (over C).
usage: fixedD_koszul.py <label> <d> <a-spec> [--target=full|k3] [--prime=P] [--modular-only] [--out=FILE]
  a-spec: python expression evaluating to a list of 6 polys in the variables U0,U1,U2,V0,V1,V2 (as dict polys via
          helper functions provided: z(i), mon(m), f_U, f_V, P[i], U[i], V[i], add, mul, scale, power), e.g.
          "[scale(U[0],Fr(1,42)), scale(U[1],Fr(1,42)), scale(U[2],Fr(1,42)), {}, {}, {}]"  gives D = f(U) (7f/7).
"""
import sys, json, time
sys.path.insert(0, '/home/kejace/devel/sextics-work'); sys.path.insert(0, '/tmp/hodge/results/polarized-k3-remainder/wave2-excl')
from fractions import Fraction as Fr
from flint import nmod_mat, fmpz_mat, fmpq_mat
from hodge import core
from hodge.core import monomials, mul, add, sub, scale, pderiv, shift, power
from relax import G, P, target, vdeg, z, mon, poly7f
NV = 6
U = [z(i) for i in range(3)]; V = [z(3+i) for i in range(3)]
f_U = scale(poly7f(0), Fr(1, 7)); f_V = scale(poly7f(3), Fr(1, 7))

def run(label, d, a, kind='full', p=65521, modular_only=False, verbose=True):
    t0 = time.time()
    T = target(kind)
    D = {}
    for i in range(NV):
        if a[i]: D = add(D, mul(P[i], a[i]))
    assert core.degree(D) == d, (core.degree(D), d)
    dD = [pderiv(D, i) for i in range(NV)]
    Ta = [mul(T, a[i]) if a[i] else {} for i in range(NV)]
    # rhs = -B_D(Ta) = -(D div(Ta) - sum Ta_i d_i D)
    divTa = {}
    for i in range(NV):
        if Ta[i]: divTa = add(divTa, pderiv(Ta[i], i))
    B0 = mul(D, divTa)
    for i in range(NV):
        if Ta[i]: B0 = sub(B0, mul(Ta[i], dD[i]))
    rhs = scale(B0, -1)
    deg2 = 2 * d
    # columns: Koszul coefficients q_jk z^m (deg m = d-4): B_D(S) for S = z^m (P_k e_j - P_j e_k)
    cols = []   # (kind, info, poly)
    for j in range(NV):
        for k in range(j + 1, NV):
            for m in monomials(d - 4):
                mz = mon(m)
                Sj = mul(mz, P[k]); Sk = scale(mul(mz, P[j]), -1)
                # div S = d_j(Sj) + d_k(Sk);  S(D) = Sj dD_j + Sk dD_k
                divS = add(pderiv(Sj, j), pderiv(Sk, k))
                col = sub(mul(D, divS), add(mul(Sj, dD[j]), mul(Sk, dD[k])))
                if col: cols.append(('q', (j, k, m), col))
    for m in monomials(deg2 - 6):
        cols.append(('Q', m, scale(shift(G, m), -1)))
    # character blocks
    def chars(Pp): return {vdeg(e) % 6 for e in Pp}
    rhs_chars = chars(rhs)
    # every column is mu6-homogeneous when D is (G, P_i, R are graded); T's inhomogeneity only spreads the rhs over
    # several characters, and the system decomposes block-wise regardless.
    homogeneous = all(len(chars(c)) <= 1 for _, _, c in cols) and len(chars(D)) == 1
    res = {'label': label, 'd': d, 'target': kind, 'prime': p, 'D_terms': len(D), 'homogeneous_mu6': homogeneous,
           'rhs_chars': sorted(rhs_chars), 'n_q_cols': sum(1 for c in cols if c[0] == 'q'), 'n_Q_cols': sum(1 for c in cols if c[0] == 'Q'),
           'blocks': [], 'D_str': core.to_str(D) if len(D) <= 40 else f'{len(D)} terms'}
    blocks = sorted(rhs_chars) if homogeneous else [None]
    consistent_p = True; consistent_Q = True
    for r in blocks:
        tb = time.time()
        if r is None:
            mons = monomials(deg2); cs = cols; rb = rhs
        else:
            mons = [m for m in monomials(deg2) if vdeg(m) % 6 == r]
            cs = [c for c in cols if chars(c[2]) == {r}]
            rb = {e: v for e, v in rhs.items() if vdeg(e) % 6 == r}
        idx = {m: i for i, m in enumerate(mons)}
        R, C = len(mons), len(cs)
        # integer matrix: clear denominators per column
        import math
        def intcol(Pp):
            den = 1
            for v in Pp.values():
                fv = Fr(v); den = den * fv.denominator // math.gcd(den, fv.denominator)
            return {e: int(Fr(v) * den) for e, v in Pp.items()}
        flatA = [0] * (R * (C + 1))
        for j, (_, _, c) in enumerate(cs):
            ic = intcol(c)
            for e, v in ic.items(): flatA[idx[e] * (C + 1) + j] = v
        irb = intcol(rb)
        for e, v in irb.items(): flatA[idx[e] * (C + 1) + C] = v
        Mp = nmod_mat(R, C + 1, [v % p for v in flatA], p)
        rank_aug_p = Mp.rank()
        M0p = nmod_mat(R, C, [flatA[i * (C + 1) + j] % p for i in range(R) for j in range(C)], p)
        rank_p = M0p.rank()
        blk = {'char': r, 'rows': R, 'cols': C, 'rank_p': rank_p, 'aug_rank_p': rank_aug_p, 'consistent_p': rank_p == rank_aug_p}
        consistent_p &= blk['consistent_p']
        if not modular_only:
            Mz = fmpz_mat(R, C + 1, flatA)
            rank_aug_Q = Mz.rank()
            M0z = fmpz_mat(R, C, [flatA[i * (C + 1) + j] for i in range(R) for j in range(C)])
            rank_Q = M0z.rank()
            blk.update({'rank_Q': rank_Q, 'aug_rank_Q': rank_aug_Q, 'consistent_Q': rank_Q == rank_aug_Q})
            consistent_Q &= blk['consistent_Q']
        blk['seconds'] = round(time.time() - tb, 1)
        res['blocks'].append(blk)
        if verbose: print('[block]', json.dumps(blk), flush=True)
    res['consistent_p'] = consistent_p
    if not modular_only: res['consistent_Q'] = consistent_Q
    res['seconds'] = round(time.time() - t0, 1)
    return res

if __name__ == '__main__':
    label = sys.argv[1]; d = int(sys.argv[2]); aspec = sys.argv[3]
    kind = 'full'; p = 65521; mo = '--modular-only' in sys.argv; out = None
    for s in sys.argv[4:]:
        if s.startswith('--target='): kind = s.split('=')[1]
        if s.startswith('--prime='): p = int(s.split('=')[1])
        if s.startswith('--out='): out = s.split('=')[1]
    a = eval(aspec)
    res = run(label, d, a, kind, p, mo)
    print('[result]', json.dumps({k: v for k, v in res.items() if k != 'blocks'}), flush=True)
    if out:
        with open(out, 'a') as fh: fh.write(json.dumps(res) + '\n')
