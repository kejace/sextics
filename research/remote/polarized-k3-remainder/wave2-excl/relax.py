"""Independent re-implementation of the L-relaxation ("Hessian relaxation") for polarized-k3-remainder.

(U,V) frame, slots 0..2 = U, 3..5 = V.  G = 7 f(U) + 7 f(V) (integer coefficients), P_i = d_i G, H_ij = d_i P_j.
Target T = prod_i (U_i^2 - V_i^2)  (= 64 * A-hat, full instance target)   or
       T = U0^2 U1^2 U2^2 - V0^2 V1^2 V2^2  (= 64 * K3-sector projection of the target, --target k3).
In degrees d where [T] is injective on the Jacobian ring, E1 forces D = sum_i P_i a_i (deg a_i = d-5);
the Euler gauge D ~ D + G h lets us drop the coefficients of a_5 divisible by z_5 (gauge-first).
Necessary condition (thread, DEGREE_SIX.md):  B_D(T a) := D div(T a) - (T a)(D)  in  L_{2d},
   L = (G) + J(G)^2 + (R_ijk),  R_ijk = P_k H_ij - P_j H_ik.
B_D(Ta) is a quadratic form in the parameters c; L is graded by V-degree mod 6 ("character" r).
We compute, per character block r: rank of the L-generator matrix mod p (and over Q with --exact-ranks),
and the constraint rows = quotient coordinates (S_{2d}/L_{2d})_r of every quadratic column.
Outputs JSON receipts + .npz of the constraint matrix.  Pure necessary condition; no HasCertificate claim.
usage: relax.py d p [--target full|k3] [--exact-ranks] [--exact-constraints] [--out PREFIX] [--no-gauge]
"""
import sys, json, time, itertools, os, math
sys.path.insert(0, '/home/kejace/devel/sextics-work')
import numpy as np
from flint import nmod_mat, fmpq_mat, fmpz_mat
from hodge import core
from hodge.core import monomials, mul, add, sub, scale, pderiv, shift

NV = 6
def z(i):
    e = [0]*6; e[i] = 1; return {tuple(e): 1}
def mon(m): return {tuple(m): 1}
def poly7f(off):
    # 7 f(z) = 7(z0^6+z1^6+z2^6) + 5(z0^5 z1 + 2 z1^5 z2 + 3 z2^5 z0 + z0^2 z1^3 z2), on slots off..off+2
    def e(*ab):
        v = [0]*6
        for i, k in enumerate(ab): v[off+i] = k
        return tuple(v)
    return {e(6,0,0): 7, e(0,6,0): 7, e(0,0,6): 7, e(5,1,0): 5, e(0,5,1): 10, e(1,0,5): 15, e(2,3,1): 5}
G = add(poly7f(0), poly7f(3))
P = [pderiv(G, i) for i in range(NV)]
H = [[pderiv(P[j], i) for j in range(NV)] for i in range(NV)]
def target(kind):
    if kind == 'full':
        T = {(0,)*6: 1}
        for i in range(3):
            T = mul(T, sub(mul(z(i), z(i)), mul(z(i+3), z(i+3))))
        return T
    if kind == 'k3':
        return {(2,2,2,0,0,0): 1, (0,0,0,2,2,2): -1}
    raise ValueError(kind)
def vdeg(e): return e[3]+e[4]+e[5]
def char_of(Pp):
    cs = {vdeg(e) % 6 for e in Pp}
    assert len(cs) == 1, cs
    return cs.pop()

def L_generators(deg2):
    """Generators of L_{deg2} as (char, poly)."""
    gens = []
    for m in monomials(deg2 - 6):
        gens.append(shift(G, m))
    if deg2 >= 10:
        PP = {}
        for i in range(NV):
            for j in range(i, NV):
                q = mul(P[i], P[j])
                if q: PP[(i, j)] = q
        for q in PP.values():
            for m in monomials(deg2 - 10):
                gens.append(shift(q, m))
    if deg2 >= 9:
        Rs, seen = [], set()
        for i in range(NV):
            for j in range(NV):
                for k in range(j+1, NV):
                    r = sub(mul(P[k], H[i][j]), mul(P[j], H[i][k]))
                    if not r: continue
                    # primitive, sign-normalized key
                    g = 0
                    for v in r.values(): g = math.gcd(g, abs(int(v)))
                    r = {e: int(v)//g for e, v in r.items()}
                    key = tuple(sorted(r.items()))
                    keyneg = tuple(sorted({e: -v for e, v in r.items()}.items()))
                    if key in seen or keyneg in seen: continue
                    seen.add(key); Rs.append(r)
        for r in Rs:
            for m in monomials(deg2 - 9):
                gens.append(shift(r, m))
    out = [(char_of(g), g) for g in gens if g]
    return out

def parameters(d, gauge=True, char=None):
    """Parameters (i, m): D = sum_i P_i a_i with a_i = sum_m c_{i,m} z^m.  char=b restricts to mu6-homogeneous D of
    character b (V-degree of P_i z^m == b mod 6; char(P_i) = 0 for i<3, 5 for i>=3); the gauge-first rule then removes
    exactly the quadrics h of character b (the directions a -> a + h z/6)."""
    pars = []
    for i in range(NV):
        for m in monomials(d - 5):
            if gauge and i == 5 and m[5] > 0: continue
            if char is not None and ((0 if i < 3 else 5) + vdeg(m)) % 6 != char % 6: continue
            pars.append((i, m))
    return pars

def quadratic_columns(d, T, pars, verbose=True):
    """Return list of (alpha, beta, poly) for alpha<=beta: coefficient polynomial of c_a c_b in B_D(Ta)."""
    t0 = time.time()
    Dp = [shift(P[i], m) for (i, m) in pars]                  # D_alpha = P_i z^m
    dT = [pderiv(T, i) for i in range(NV)]
    Vp = [add(shift(dT[i], m), mul(T, pderiv(mon(m), i))) for (i, m) in pars]   # d_i (T z^m)
    dD = [[pderiv(Dp[b], j) for j in range(NV)] for b in range(len(pars))]      # d_j D_beta
    Tz = [mul(T, mon(m)) for (i, m) in pars]                 # T z^m
    cols = []
    n = len(pars)
    for a in range(n):
        ia = pars[a][0]
        for b in range(a, n):
            ib = pars[b][0]
            if a == b:
                q = sub(mul(Dp[a], Vp[a]), mul(Tz[a], dD[a][ia]))
            else:
                q = add(mul(Dp[a], Vp[b]), mul(Dp[b], Vp[a]))
                q = sub(q, add(mul(Tz[a], dD[b][ia]), mul(Tz[b], dD[a][ib])))
            cols.append((a, b, q))
    if verbose:
        print(f'[quad] d={d} params={n} pairs={len(cols)} in {time.time()-t0:.1f}s', flush=True)
    return cols

def build(d, p, kind='full', gauge=True, exact_ranks=False, exact_constraints=False, out=None, verbose=True, char=None):
    t0 = time.time()
    T = target(kind)
    deg2 = 2 * d
    pars = parameters(d, gauge, char)
    gens = L_generators(deg2)
    cols = quadratic_columns(d, T, pars, verbose)
    npairs = len(cols)
    result = {'d': d, 'p': p, 'target': kind, 'gauge': gauge, 'char': char, 'params': len(pars), 'pairs': npairs,
              'L_generators': len(gens), 'blocks': []}
    allmons = monomials(deg2)
    constraints = []     # list of (block r, nmod rows as python lists)  -> assembled below
    C_blocks = []
    exact_rows = []
    for r in range(6):
        tb = time.time()
        mons = [m for m in allmons if vdeg(m) % 6 == r]
        idx = {m: i for i, m in enumerate(mons)}
        gr = [g for (c, g) in gens if c == r]
        mr, nr = len(mons), len(gr)
        flat = [0] * (mr * nr)
        for j, g in enumerate(gr):
            for e, v in g.items():
                flat[idx[e] * nr + j] = v % p
        Lam = nmod_mat(mr, nr, flat, p)
        rank_p = Lam.rank()
        blk = {'r': r, 'rows': mr, 'L_cols': nr, 'L_rank_p': rank_p}
        if exact_ranks or exact_constraints:
            Lq = fmpz_mat(mr, nr, [0]*(mr*nr))
            for j, g in enumerate(gr):
                for e, v in g.items():
                    Lq[idx[e], j] = int(v)
            rank_q = Lq.rank()
            blk['L_rank_Q'] = rank_q
        # quotient map mod p: Y = nullspace(Lam^T)  (mr x (mr-rank))
        Yn, nul = Lam.transpose().nullspace()
        assert nul == mr - rank_p
        Y = np.array([[int(Yn[i, k]) for k in range(nul)] for i in range(mr)], dtype=np.int64)  # mr x nul
        # quadratic columns restricted to block r, as sparse
        ri, ci, vals = [], [], []
        for j, (a, b, q) in enumerate(cols):
            for e, v in q.items():
                if vdeg(e) % 6 == r:
                    ri.append(idx[e]); ci.append(j); vals.append(int(v) % p)
        import scipy.sparse as sp
        Qm = sp.csr_matrix((np.array(vals, dtype=np.int64), (np.array(ri), np.array(ci))), shape=(mr, npairs))
        # C = Y^T Q  (nul x npairs) mod p with 16-bit splitting of Y
        QT = Qm.T.tocsr()
        Ylo = Y & 0xFFFF; Yhi = Y >> 16
        lo = (QT @ Ylo) % p; hi = (QT @ Yhi) % p
        C = ((lo + (hi * (1 << 16)) % p) % p).T   # nul x npairs
        C = np.ascontiguousarray(C)
        Cn = nmod_mat(nul, npairs, C.ravel().tolist(), p)
        crank = Cn.rank()
        blk.update({'quotient_dim_p': nul, 'constraint_rank_p': crank, 'seconds': round(time.time()-tb, 1)})
        C_blocks.append(C)
        if exact_constraints:
            # exact rational augmented rref  [Lambda | Q]  -> rows beyond L-rank give exact constraints
            aug = fmpq_mat(mr, nr + npairs)
            for j, g in enumerate(gr):
                for e, v in g.items():
                    aug[idx[e], j] = int(v)
            for j, (a, b, q) in enumerate(cols):
                for e, v in q.items():
                    if vdeg(e) % 6 == r:
                        aug[idx[e], nr + j] = int(v)
            rr, rk = aug.rref()
            assert rk >= blk['L_rank_Q']
            for i in range(blk['L_rank_Q'], rk):
                assert all(rr[i, j] == 0 for j in range(nr))
                exact_rows.append([rr[i, nr + j] for j in range(npairs)])
            blk['aug_rank_Q'] = rk
        result['blocks'].append(blk)
        if verbose:
            print('[block]', json.dumps(blk), flush=True)
    Call = np.vstack(C_blocks) if C_blocks else np.zeros((0, npairs), dtype=np.int64)
    Cn = nmod_mat(Call.shape[0], npairs, Call.ravel().tolist(), p)
    result['constraint_rows'] = int(Call.shape[0])
    result['constraint_rank_p'] = Cn.rank()
    result['seconds'] = round(time.time() - t0, 1)
    result['pairs_index'] = [(a, b) for (a, b, q) in cols]
    result['parameter_basis'] = [(i, list(m)) for (i, m) in pars]
    if out:
        np.savez_compressed(out + '.npz', C=Call, pairs=np.array([(a, b) for (a, b, q) in cols], dtype=np.int32))
        json.dump(result, open(out + '.json', 'w'))
        if exact_constraints:
            json.dump([[str(x) for x in row] for row in exact_rows], open(out + '-exactQ.json', 'w'))
    if verbose:
        print('[result]', json.dumps({k: v for k, v in result.items() if k not in ('pairs_index', 'parameter_basis')}), flush=True)
    return result, Call, cols, pars, exact_rows

def square_elimination(C, pairs, nparams, p=None, exact_rows=None, verbose=True):
    """Iteratively: if c_a^2 is in the row span of the constraints (restricted to surviving params), force c_a = 0.
    Modular (p given, C numpy) or exact (exact_rows: list of lists of fmpq). Returns (killed set, steps)."""
    killed = set(); steps = []
    while True:
        keep = [j for j, (a, b) in enumerate(pairs) if a not in killed and b not in killed]
        if not keep:
            steps.append({'remaining': 0, 'note': 'all parameters forced to zero'}); break
        if exact_rows is not None:
            M = fmpq_mat([[row[j] for j in keep] for row in exact_rows])
            rr, rank = M.rref()
            new = set()
            for i in range(rank):
                nz = [j for j in range(len(keep)) if rr[i, j] != 0]
                if len(nz) == 1 and pairs[keep[nz[0]]][0] == pairs[keep[nz[0]]][1]:
                    new.add(pairs[keep[nz[0]]][0])
        else:
            sub_ = np.ascontiguousarray(C[:, keep])
            M = nmod_mat(sub_.shape[0], sub_.shape[1], sub_.ravel().tolist(), p)
            rr, rank = M.rref()
            new = set()
            for i in range(rank):
                nz = [j for j in range(len(keep)) if int(rr[i, j]) != 0]
                if len(nz) == 1 and pairs[keep[nz[0]]][0] == pairs[keep[nz[0]]][1]:
                    new.add(pairs[keep[nz[0]]][0])
        steps.append({'remaining': nparams - len(killed), 'rank': int(rank), 'forced': sorted(new)})
        if verbose: print('[squares]', json.dumps(steps[-1]), flush=True)
        if not new: break
        killed |= new
    return killed, steps

if __name__ == '__main__':
    d = int(sys.argv[1]); p = int(sys.argv[2])
    kind = 'full'; exact_ranks = '--exact-ranks' in sys.argv; exact_constraints = '--exact-constraints' in sys.argv
    gauge = '--no-gauge' not in sys.argv
    for a in sys.argv:
        if a.startswith('--target='): kind = a.split('=', 1)[1]
    out = None; char = None
    for a in sys.argv:
        if a.startswith('--out='): out = a.split('=', 1)[1]
        if a.startswith('--char='): char = int(a.split('=', 1)[1])
    res, C, cols, pars, exact_rows = build(d, p, kind, gauge, exact_ranks, exact_constraints, out, char=char)
    pairs = [(a, b) for (a, b, q) in cols]
    if '--squares' in sys.argv:
        killed, steps = square_elimination(C, pairs, len(pars), p=p)
        res['square_elimination_mod_p'] = {'killed': sorted(killed), 'steps': steps}
        if exact_constraints:
            killedq, stepsq = square_elimination(None, pairs, len(pars), exact_rows=exact_rows)
            res['square_elimination_Q'] = {'killed': sorted(killedq), 'steps': stepsq}
        if out: json.dump(res, open(out + '.json', 'w'))
        print('[squares-summary]', json.dumps({'killed_mod_p': len(killed), 'params': len(pars),
              'killed_Q': (len(killedq) if exact_constraints else None)}), flush=True)
