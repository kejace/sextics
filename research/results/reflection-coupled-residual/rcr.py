"""Infrastructure for the reflection-coupled-residual attack.
Jacobian-ring normal forms mod p, r2-sector split, served-space accumulation, logging."""
import sys, time, json, os, itertools
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core, fixed_d
from hodge.core import NV, monomials, mul, add, sub, scale, pderiv, shift
from flint import nmod_mat
from fractions import Fraction

IID = 'reflection-coupled-residual'
WORK = f'/tmp/hodge/results/{IID}'
LOG = f'{WORK}/log.jsonl'
P = core.PRIMES[0]
F, A, q, rec = core.load_instance(IID)
Fp = core.reduce_mod(F, P); Ap = core.reduce_mod(A, P)
Fi = [pderiv(Fp, i, P) for i in range(NV)]

def log(D_expr, d, family, e1, consistent, rank, nullity, seconds, note='', **extra):
    ent = {"D": D_expr, "d": d, "family": family, "e1": e1, "consistent": consistent,
           "rank": rank, "nullity": nullity, "seconds": round(seconds, 2), "note": note, "engine": 2}
    ent.update(extra)
    with open(LOG, 'a') as f:
        f.write(json.dumps(ent) + '\n')

def y2par(e):
    return e[5] & 1

class JacobianDegree:
    """Normal forms modulo the Jacobian ideal in a fixed degree and Y-parity (mod P)."""
    def __init__(self, deg, par):
        self.deg, self.par = deg, par
        self.mons = monomials(deg, par)
        self.idx = {m: i for i, m in enumerate(self.mons)}
        rows = []
        for i in range(NV):
            # multiplier monomials of degree deg-5 with parity so that m*F_i has parity par
            pi = 0 if i < 3 else 1
            for m in monomials(deg - 5, (par - pi) % 2):
                rows.append(shift(Fi[i], m, 1, P))
        self.nJ = len(rows)
        R, C = len(rows), len(self.mons)
        if R == 0:
            self.rank = 0; self.pivots = []; self.Mr = None
        else:
            flat = [0] * (R * C)
            for r, poly in enumerate(rows):
                for e, v in poly.items():
                    flat[r * C + self.idx[e]] = v
            M = nmod_mat(R, C, flat, P)
            Mr, rank = M.rref()
            self.rank = rank
            piv = []
            for i in range(rank):
                j = 0
                while int(Mr[i, j]) == 0: j += 1
                piv.append(j)
            self.pivots = piv
            self.Mr = Mr
        self.nonpiv = [j for j in range(C) if j not in set(self.pivots)]
        self.quot_dim = C - self.rank
        self.npidx = {j: k for k, j in enumerate(self.nonpiv)}

    def normal_form(self, poly):
        """Return coordinates (list length quot_dim, mod P) of poly's class in S/J."""
        v = [0] * len(self.mons)
        for e, c in poly.items():
            v[self.idx[e]] = (v[self.idx[e]] + int(c)) % P
        for i, j in enumerate(self.pivots):
            if v[j]:
                c = v[j]
                for k in range(len(self.mons)):
                    a = int(self.Mr[i, k])
                    if a:
                        v[k] = (v[k] - c * a) % P
        return [v[j] for j in self.nonpiv]

    def sector_mask(self, y2parity):
        return [k for k, j in enumerate(self.nonpiv) if y2par(self.mons[j]) == y2parity]

class Span:
    """Incremental span mod P with membership tests (row-reduced basis)."""
    def __init__(self, n):
        self.n = n; self.rows = []  # list of (pivot, vector) with vector[pivot]=1, reduced
    def reduce(self, v):
        v = list(v)
        for piv, r in self.rows:
            if v[piv]:
                c = v[piv]
                for k in range(self.n):
                    if r[k]:
                        v[k] = (v[k] - c * r[k]) % P
        return v
    def add(self, v):
        v = self.reduce(v)
        piv = next((k for k in range(self.n) if v[k]), None)
        if piv is None:
            return False
        inv = pow(v[piv], -1, P)
        v = [(x * inv) % P for x in v]
        # keep other rows reduced w.r.t. new pivot
        for i, (pv, r) in enumerate(self.rows):
            if r[piv]:
                c = r[piv]
                self.rows[i] = (pv, [(r[k] - c * v[k]) % P for k in range(self.n)])
        self.rows.append((piv, v))
        return True
    def contains(self, v):
        return not any(self.reduce(v))
    def dim(self):
        return len(self.rows)
    def save(self, path):
        json.dump({'n': self.n, 'rows': self.rows}, open(path, 'w'))
    @classmethod
    def load(cls, path):
        d = json.load(open(path)); s = cls(d['n']); s.rows = [(p, v) for p, v in d['rows']]; return s

J6 = None
def get_J6():
    global J6
    if J6 is None:
        J6 = JacobianDegree(6 * q, 1)
    return J6

def served(D, verbose=False):
    """Served space of D: returns dict with dim, basis (in R_6^- coords), sector dims, timing."""
    t0 = time.time()
    dim, basis, S = fixed_d.target_image(F, q, D, p=P, verbose=verbose)
    t1 = time.time()
    J = get_J6()
    vecs = [J.normal_form(b) for b in basis]
    # rank of reduced vectors
    sp = Span(J.quot_dim)
    for v in vecs:
        sp.add(v)
    even_mask = J.sector_mask(0)  # targets even in Y2 -> r2-anti-invariant residues (our sector)
    odd_mask = J.sector_mask(1)
    # project (valid if D is r2-invariant)
    sp_even = Span(J.quot_dim); sp_odd = Span(J.quot_dim)
    for v in vecs:
        ve = [v[k] if k in set(even_mask) else 0 for k in range(J.quot_dim)]
        vo = [v[k] if k in set(odd_mask) else 0 for k in range(J.quot_dim)]
        sp_even.add(ve); sp_odd.add(vo)
    return {'dim_raw': dim, 'dim_modJ': sp.dim(), 'dim_even': sp_even.dim(), 'dim_odd': sp_odd.dim(),
            'span': sp, 'span_even': sp_even, 'span_odd': sp_odd, 'rows': S.nrows, 'cols': S.ncols,
            't_image': round(t1 - t0, 1), 't_total': round(time.time() - t0, 1)}

def A_vec():
    return get_J6().normal_form(Ap)

def Dstr(D):
    return core.to_str(D)
