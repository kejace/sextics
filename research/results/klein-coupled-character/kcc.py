"""klein-coupled-character attack helpers: sigma-weight (Z/7) decomposed fixed-D solver.

sigma acts with weights w = (1,2,4,0,0,0) mod 7; F has weight 0, A = X0^3 Y0^3 has weight 3.
If D is a sigma-eigenvector of weight k, every term of the cleared identity lives in weight 3+2k,
and the numerator N_i may be restricted to weight 3+k+w_i, Q to weight 3+2k (averaging argument).
"""
import sys, time, json, os
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core, fixed_d
from hodge.core import NV, monomials, shift, mul, add, sub, scale, pderiv, power, yparity
import flint
from flint import nmod_mat
import numpy as np

W = (1, 2, 4, 0, 0, 0)
P0 = core.PRIMES[0]
LOG = '/tmp/hodge/results/klein-coupled-character/log.jsonl'

def wt(e):
    return (e[0] + 2*e[1] + 4*e[2]) % 7

def sigma_weight(P):
    ws = {wt(e) for e in P}
    return ws.pop() if len(ws) == 1 else None

def wmonomials(deg, par, k):
    """monomials of given degree, Y-parity (None=any), sigma-weight k (None=any)."""
    return [e for e in monomials(deg, par) if k is None or wt(e) == k]

def log(entry):
    entry.setdefault('time', time.strftime('%H:%M:%S'))
    with open(LOG, 'a') as f:
        f.write(json.dumps(entry, default=str) + '\n')

class WSystem:
    """Fixed-D linear system restricted to the sigma-weight block (when D is an eigenvector)."""
    def __init__(self, F, A, q, D, p=P0, target_weight=None, with_target_space=False, verbose=False, weights=True):
        t0 = time.time()
        self.p, self.q = p, q
        Fp, Ap, Dp = core.reduce_mod(F, p), core.reduce_mod(A, p), core.reduce_mod(D, p)
        self.F, self.A, self.D = Fp, Ap, Dp
        d = core.degree(D); self.d = d
        pars = fixed_d.numerator_parities(D, A)
        k = sigma_weight(D); self.k = k
        a = sigma_weight(A) if target_weight is None else target_weight
        self.a = a
        restrict = (weights and k is not None and a is not None)
        self.restrict = restrict
        degN, degQ = d + 6*q - 5, 2*d - 6
        Fi = [pderiv(Fp, i, p) for i in range(NV)]
        Di = [pderiv(Dp, i, p) for i in range(NV)]
        FD = mul(Fp, Dp, p)
        G = [scale(add(scale(mul(Dp, Fi[i], p), q, p), mul(Fp, Di[i], p), p), p-1, p) for i in range(NV)]  # -(q D F_i + F D_i)  [engine v2 sign]
        Fq1 = scale(power(Fp, q+1, p), p-1, p)
        cols, colpolys = [], []
        for i in range(NV):
            kN = (a + k + W[i]) % 7 if restrict else None
            for m in wmonomials(degN, pars[i], kN):
                col = shift(G[i], m, 1, p)
                if m[i]:
                    m2 = list(m); m2[i] -= 1
                    col = add(col, shift(FD, tuple(m2), m[i], p), p)
                cols.append(('N', i, m)); colpolys.append(col)
        self.nN = len(cols)
        if degQ >= 0:
            kQ = (a + 2*k) % 7 if restrict else None
            for m in wmonomials(degQ, core.parity_of(A), kQ):
                cols.append(('Q', -1, m)); colpolys.append(shift(Fq1, m, 1, p))
        self.nQ = len(cols) - self.nN
        D2 = mul(Dp, Dp, p)
        if with_target_space:
            kA = a if restrict else None
            for m in wmonomials(6*q, 1, kA):
                cols.append(('A', -1, m)); colpolys.append(shift(D2, m, 1, p))
            rhs = {}
        else:
            rhs = scale(mul(Ap, D2, p), p-1, p)
        rows = {}
        for col in colpolys:
            for e in col:
                if e not in rows: rows[e] = len(rows)
        for e in rhs:
            if e not in rows: rows[e] = len(rows)
        self.cols, self.colpolys, self.rows, self.rhs = cols, colpolys, rows, rhs
        self.nrows, self.ncols = len(rows), len(cols)
        self.nnz = sum(len(c) for c in colpolys)
        self.build_time = time.time() - t0
        if verbose:
            print(f'[wsys] d={d} k={k} a={a} restrict={restrict} rows={self.nrows} cols={self.ncols} '
                  f'(N {self.nN}, Q {self.nQ}) nnz={self.nnz} built {self.build_time:.1f}s', file=sys.stderr)

    def dense(self, with_rhs=True):
        R, C, p = self.nrows, self.ncols, self.p
        Wd = C + (1 if with_rhs else 0)
        flat = [0] * (R * Wd)
        rows = self.rows
        for j, col in enumerate(self.colpolys):
            for e, v in col.items():
                flat[rows[e] * Wd + j] = v
        if with_rhs:
            for e, v in self.rhs.items():
                flat[rows[e] * Wd + C] = v
        return nmod_mat(R, Wd, flat, p)

    def solve(self, verbose=False, compress=None):
        """consistency / rank / nullity (mod p). compress: None=auto."""
        t0 = time.time()
        R, C = self.nrows, self.ncols
        if compress is None:
            compress = R > C + 200 and R * (C+1) > 4e8
        if compress:
            r = fixed_d.solve_compressed(self, k_extra=60, verbose=verbose)
            r['compressed'] = True
            r['seconds'] = round(time.time() - t0, 2)
            return r
        M = self.dense()
        Mr, rank = M.rref()
        pivots = []
        for i in range(rank):
            j = 0
            while j <= C and int(Mr[i, j]) == 0:
                j += 1
            pivots.append(j)
        consistent = not (pivots and pivots[-1] == C)
        x = None
        if consistent:
            x = [0] * C
            for i, j in enumerate(pivots):
                x[j] = int(Mr[i, C])
        res = {'consistent': consistent, 'rank': rank, 'nullity': C - rank, 'x': x, 'pivots': pivots,
               'rows': R, 'cols': C, 'compressed': False, 'seconds': round(time.time() - t0, 2)}
        if verbose:
            print(f'[wsolve p={self.p}] consistent={consistent} rank={rank} nullity={C-rank} {res["seconds"]}s', file=sys.stderr)
        return res

def e1_block(F, A, q, D, p=P0, weights=True):
    """E1 (q N(F) = A D) restricted to the weight block. Returns (ok, info)."""
    Fp, Ap, Dp = core.reduce_mod(F, p), core.reduce_mod(A, p), core.reduce_mod(D, p)
    d = core.degree(D)
    pars = fixed_d.numerator_parities(D, A)
    k = sigma_weight(D); a = sigma_weight(A)
    restrict = weights and k is not None and a is not None
    degN = d + 6*q - 5
    Fi = [pderiv(Fp, i, p) for i in range(NV)]
    cols = []
    for i in range(NV):
        kN = (a + k + W[i]) % 7 if restrict else None
        for m in wmonomials(degN, pars[i], kN):
            cols.append(shift(Fi[i], m, q, p))
    rhs = mul(Ap, Dp, p)
    rows = {}
    for col in cols:
        for e in col: rows.setdefault(e, len(rows))
    for e in rhs: rows.setdefault(e, len(rows))
    R, C = len(rows), len(cols)
    flat = [0] * (R * (C+1))
    for j, col in enumerate(cols):
        for e, v in col.items(): flat[rows[e]*(C+1)+j] = v
    for e, v in rhs.items(): flat[rows[e]*(C+1)+C] = v
    M = nmod_mat(R, C+1, flat, p)
    ra = M.rank()
    M0 = nmod_mat(R, C, [flat[r*(C+1)+j] for r in range(R) for j in range(C)], p)
    r0 = M0.rank()
    return ra == r0, {'rank': r0, 'rank_aug': ra, 'rows': R, 'cols': C, 'k': k}

def test_D(F, A, q, D, name, family, p=P0, verbose=False, full=True, note='', compress=None, weights=True):
    """E1 block test, then (if E1 ok and full) the full fixed-D system. Logs and returns a dict."""
    t0 = time.time()
    d = core.degree(D)
    k = sigma_weight(D)
    if not D or not core.is_homogeneous(D):
        raise ValueError('D must be nonzero homogeneous')
    if core.parity_of(D) is None:
        raise ValueError('D must have pure Y-parity')
    ok, info = e1_block(F, A, q, D, p, weights=weights)
    entry = {'D': name if len(name) < 200 else name[:200]+'...', 'd': d, 'family': family, 'k': k, 'e1': ok, 'engine': 2,
             'e1_rank': info['rank'], 'e1_rows': info['rows'], 'e1_cols': info['cols'],
             'consistent': None, 'rank': None, 'nullity': None, 'note': note}
    if ok and full:
        S = WSystem(F, A, q, D, p, verbose=verbose, weights=weights)
        r = S.solve(verbose=verbose, compress=compress)
        entry.update({'consistent': r['consistent'], 'rank': r['rank'], 'nullity': r['nullity'],
                      'rows': S.nrows, 'cols': S.ncols, 'compressed': r['compressed']})
        entry['_system'] = S; entry['_result'] = r
    entry['seconds'] = round(time.time() - t0, 2)
    log({kk: v for kk, v in entry.items() if not kk.startswith('_')})
    return entry

def served_targets(F, q, D, p=P0, target_weight=None, verbose=False):
    """dimension and basis of weight-restricted Y-odd degree-6q targets served by D (mod p)."""
    A0 = {monomials(6*q, 1)[0]: 1}
    S = WSystem(F, A0, q, D, p, target_weight=target_weight, with_target_space=True, verbose=verbose)
    M = S.dense(with_rhs=False)
    X, nullity = M.nullspace()
    a_idx = [j for j, c in enumerate(S.cols) if c[0] == 'A']
    if nullity == 0 or not a_idx:
        return 0, [], S
    proj = nmod_mat(len(a_idx), nullity, [int(X[j, kk]) for j in a_idx for kk in range(nullity)], p)
    T = proj.transpose()
    Tr, rk = T.rref()
    basis = []
    for i in range(rk):
        poly = {}
        for jj, j in enumerate(a_idx):
            v = int(Tr[i, jj])
            if v: poly[S.cols[j][2]] = v
        basis.append(poly)
    return rk, basis, S

def lift_solution(F, A, q, D, S, r, verbose=False, weights=True):
    """Try to lift a consistent mod-p solution (free vars 0) to Q using several primes + CRT."""
    primes = core.PRIMES
    acc = None; pivots0 = r['pivots']
    xs_list = []
    for kk, p in enumerate(primes[:8]):
        Sp = WSystem(F, A, q, D, p, weights=weights)
        rp = Sp.solve(compress=False)
        if not rp['consistent']:
            return None
        if rp['pivots'] != pivots0:
            if verbose: print('[lift] pivot mismatch', file=sys.stderr)
            acc = None; pivots0 = rp['pivots']
        x = rp['x']
        if acc is None: acc = (x, p)
        else:
            xs, m = acc
            xs = [core.crt_pair(aa, m, b, p)[0] for aa, b in zip(xs, x)]
            acc = (xs, m*p)
        xs, m = acc
        rat = [core.ratrecon(aa, m) for aa in xs]
        if all(v is not None for v in rat):
            N = [{} for _ in range(NV)]; Q = {}
            for (kind, i, mon), v in zip(Sp.cols, rat):
                if v == 0: continue
                if kind == 'N': N[i][mon] = v
                else: Q[mon] = v
            ok, info = core.verify_witness(F, A, q, D, N, Q)
            if ok:
                return {'N': N, 'Q': Q, 'primes': kk+1}
    return None

def load():
    F, A, q, rec = core.load_instance('klein-coupled-character')
    return F, A, q

if __name__ == '__main__':
    F, A, q = load()
    print('F weight', sigma_weight(F), 'A weight', sigma_weight(A))
