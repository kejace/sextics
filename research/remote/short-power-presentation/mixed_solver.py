"""Parity-free fixed-D system (for denominators of mixed Y-parity).  Consistency mod p via compressed rows,
optionally dense rref.  Mirrors fixed_d.FixedDSystem with no parity restriction on N and Q."""
import sys, time
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core, fixed_d
from hodge.core import NV, monomials, shift, mul, add, scale, pderiv, power

class MixedSystem:
    def __init__(self, F, A, q, D, p, verbose=False):
        t0 = time.time()
        self.p = p; self.q = q
        Fp = core.reduce_mod(F, p); Ap = core.reduce_mod(A, p); Dp = core.reduce_mod(D, p)
        d = core.degree(D); self.d = d
        degN = d + 6 * q - 5; degQ = 2 * d - 6
        Fi = [pderiv(Fp, i, p) for i in range(NV)]
        Di = [pderiv(Dp, i, p) for i in range(NV)]
        FD = mul(Fp, Dp, p)
        G = [scale(add(scale(mul(Dp, Fi[i], p), q, p), mul(Fp, Di[i], p), p), p - 1, p) for i in range(NV)]
        Fq1 = scale(power(Fp, q + 1, p), p - 1, p)
        cols, colpolys = [], []
        for i in range(NV):
            for m in monomials(degN, None):
                col = shift(G[i], m, 1, p)
                if m[i]:
                    m2 = list(m); m2[i] -= 1
                    col = add(col, shift(FD, tuple(m2), m[i], p), p)
                cols.append(('N', i, m)); colpolys.append(col)
        for m in monomials(degQ, None):
            cols.append(('Q', -1, m)); colpolys.append(shift(Fq1, m, 1, p))
        rhs = scale(mul(Ap, mul(Dp, Dp, p), p), p - 1, p)
        rows = {}
        for col in colpolys:
            for e in col: rows.setdefault(e, len(rows))
        for e in rhs: rows.setdefault(e, len(rows))
        self.cols, self.colpolys, self.rhs, self.rows = cols, colpolys, rhs, rows
        self.nrows, self.ncols = len(rows), len(cols)
        if verbose:
            print(f'[mixed system] d={d} rows={self.nrows} cols={self.ncols} built {time.time()-t0:.1f}s', file=sys.stderr, flush=True)

def mixed_consistent(F, A, q, D, p=None, verbose=True, k_extra=60):
    p = p or core.PRIMES[0]
    S = MixedSystem(F, A, q, D, p, verbose=verbose)
    r = fixed_d.solve_compressed(S, k_extra=k_extra, verbose=verbose)
    return {'consistent': r['consistent'], 'rank': r['rank'], 'nullity': r['nullity'], 'rows': S.nrows, 'cols': S.ncols, 'd': S.d}
