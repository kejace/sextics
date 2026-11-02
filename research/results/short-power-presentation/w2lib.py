"""Wave-2 helpers for short-power-presentation: logging, kernel bases, reducers (mod p)."""
import sys, os, json, time
sys.path.insert(0, '/home/kejace/devel/sextics-work')
sys.path.insert(0, '/tmp/hodge/results/short-power-presentation')
from hodge import core, fixed_d, jacobian
from fractions import Fraction
from flint import nmod_mat
import numpy as np

IID = 'short-power-presentation'
RES = '/tmp/hodge/results/short-power-presentation'
LOG = os.path.join(RES, 'wave2-search.log.jsonl')
PRIMES = list(core.PRIMES)

def wlog(**rec):
    rec.setdefault('t', time.strftime('%Y-%m-%dT%H:%M:%S'))
    rec.setdefault('engine', 2)
    with open(LOG, 'a') as f:
        f.write(json.dumps(rec, default=str) + '\n')
    return rec

def load():
    return core.load_instance(IID)

def jac_cols(Fp, deg, par, p):
    """Generators m*F_i of J(F)_deg with Y-parity par (mod p)."""
    Fi = [core.pderiv(Fp, i, p) for i in range(6)]
    cols, meta = [], []
    for i in range(6):
        pi = 0 if i < 3 else 1
        for m in core.monomials(deg - 5, (par - pi) % 2):
            cols.append(core.shift(Fi[i], m, 1, p)); meta.append((i, m))
    return cols, meta

def mat_from_polys(polys, p, mons):
    idx = {m: i for i, m in enumerate(mons)}
    R = len(polys); C = len(mons)
    flat = [0] * (R * C)
    for r, P in enumerate(polys):
        for e, v in P.items(): flat[r * C + idx[e]] = v
    return nmod_mat(R, C, flat, p)

class Reducer:
    """Normal form modulo J(F)_deg^par over GF(p): rref of the generator matrix (rows = generators)."""
    def __init__(self, Fp, deg, par, p, verbose=False):
        t0 = time.time()
        self.p, self.deg, self.par = p, deg, par
        self.mons = core.monomials(deg, par)
        self.idx = {m: i for i, m in enumerate(self.mons)}
        gens, _ = jac_cols(Fp, deg, par, p)
        M = mat_from_polys(gens, p, self.mons)
        Rr, r = M.rref()
        n = len(self.mons)
        self.rank = r
        self.R = nmod_mat(r, n, [int(Rr[i, j]) for i in range(r) for j in range(n)], p) if r else None
        piv = []
        for i in range(r):
            j = 0
            while int(Rr[i, j]) == 0: j += 1
            piv.append(j)
        self.pivots = piv
        pset = set(piv)
        self.standard = [self.mons[j] for j in range(n) if j not in pset]
        if verbose:
            print(f'[reducer] deg {deg} par {par}: dim S={n} rank J={r} dim R={n-r} ({time.time()-t0:.0f}s)', flush=True)

    def nf_rows(self, polys):
        """Normal forms of a list of polys (mod p) as an nmod_mat (rows)."""
        U = mat_from_polys(polys, self.p, self.mons)
        if self.R is None: return U
        k = len(polys)
        Up = nmod_mat(k, self.rank, [int(U[i, pj]) for i in range(k) for pj in self.pivots], self.p)
        return U - Up * self.R

    def nf(self, P):
        N = self.nf_rows([P])
        return {self.mons[j]: int(N[0, j]) for j in range(len(self.mons)) if int(N[0, j])}

def poly_from_row(M, i, mons):
    return {mons[j]: int(M[i, j]) for j in range(len(mons)) if int(M[i, j])}

def kernel_reps(Fp, Ap, d, parD, p, red_d=None, verbose=False):
    """ker([A]: R_d^parD -> R_{d+6}) reps reduced mod J_d (rref-normal forms).  Returns (reps, info, reducer_d)."""
    t0 = time.time()
    monsD = core.monomials(d, parD); nD = len(monsD)
    Jcols, _ = jac_cols(Fp, d + 6, (parD + core.parity_of(Ap)) % 2, p)
    AD = [core.shift(Ap, m, 1, p) for m in monsD]
    rows = {}
    for P in Jcols + AD:
        for e in P: rows.setdefault(e, len(rows))
    R = len(rows); C = len(Jcols) + nD
    flat = [0] * (R * C)
    for j, P in enumerate(Jcols + AD):
        for e, v in P.items(): flat[rows[e] * C + j] = v
    M = nmod_mat(R, C, flat, p)
    del flat
    X, nul = M.nullspace()
    nJ = len(Jcols)
    Dvecs = [{monsD[j]: int(X[nJ + j, k]) for j in range(nD) if int(X[nJ + j, k])} for k in range(nul)]
    del X, M
    if red_d is None:
        red_d = Reducer(Fp, d, parD, p)
    NF = red_d.nf_rows(Dvecs) if Dvecs else None
    reps = []
    if NF is not None:
        NFr, rk = NF.rref()
        for i in range(rk):
            reps.append(poly_from_row(NFr, i, monsD))
    info = {'d': d, 'parD': parD, 'p': p, 'nullity_full': nul, 'dim J_d': red_d.rank, 'dim S_d': nD,
            'dim R_d': nD - red_d.rank, 'kernel_mod_J': len(reps), 'rows': R, 'cols': C, 'seconds': round(time.time() - t0, 1)}
    if verbose: print(info, flush=True)
    return reps, info, red_d

def crt_reconstruct(reps_by_prime, primes):
    """CRT + rational reconstruction across primes for lists of normal-form reps (same order). None if fails."""
    out = []
    for k in range(len(reps_by_prime[0])):
        keys = set()
        for reps in reps_by_prime: keys |= set(reps[k])
        D = {}
        for m in keys:
            x, mod = reps_by_prime[0][k].get(m, 0), primes[0]
            for reps, pr in zip(reps_by_prime[1:], primes[1:]):
                x, mod = core.crt_pair(x, mod, reps[k].get(m, 0), pr)
            v = core.ratrecon(x, mod)
            if v is None: return None, k
            if v: D[m] = v
        out.append(D)
    return out, None

def save_poly_json(path, items):
    with open(path, 'w') as f: json.dump(items, f, indent=1)
