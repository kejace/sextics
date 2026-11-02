"""Common helpers for the wave-2 exclusion lens (eight-plane-configuration)."""
import sys, os, json, time
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core, fixed_d, jacobian
from hodge.core import monomials, mul, add, sub, scale, shift, pderiv, reduce_mod, parity_of, degree
from flint import nmod_mat
from math import comb

IID = 'eight-plane-configuration'
WD = f'/tmp/hodge/results/{IID}'
PK = f'{WD}/wave2-excl'
F, A, q, rec = core.load_instance(IID)
NV = 6
P_NEW = 2147483587   # core.PRIMES[2]; wave 1 used PRIMES[0]=2147483647 (local) and PRIMES[1]=2147483629 (remote)

def wlog(entry, fn=None):
    entry.setdefault('ts', time.strftime('%Y-%m-%dT%H:%M:%S'))
    fn = fn or f'{WD}/wave2-excl.log.jsonl'
    with open(fn, 'a') as f:
        f.write(json.dumps(entry) + '\n')

def S_dim(k, par=None):
    tot = 0
    for j in range(k + 1):
        if par is not None and j % 2 != par: continue
        tot += comb(j + 2, 2) * comb(k - j + 2, 2)
    return tot

def R_dim(k, par):
    v = 0
    for a in range(4):
        for b in range(4):
            kk = k - 5 * (a + b)
            if kk < 0: continue
            v += (-1) ** (a + b) * comb(3, a) * comb(3, b) * S_dim(kk, (par - b) % 2)
    return v

def J_dim(k, par):
    return S_dim(k, par) - R_dim(k, par)

def partials(p):
    Fp = reduce_mod(F, p)
    return Fp, [pderiv(Fp, i, p) for i in range(NV)]

def J_generators(Fi, deg, par, p):
    """list of polys m*F_i spanning J_deg^par (mod p)."""
    out = []
    for i in range(NV):
        parM = (par - (0 if i < 3 else 1)) % 2
        for m in monomials(deg - 5, parM):
            out.append(shift(Fi[i], m, 1, p))
    return out

def vectorize(polys, mons=None):
    if mons is None:
        mons = {}
        for P_ in polys:
            for e in P_: mons.setdefault(e, len(mons))
    return mons

class Reducer:
    """Normal forms modulo J_deg^par (mod p) via rref of the generator matrix (rows = generators)."""
    def __init__(self, Fi, deg, par, p):
        self.p, self.deg, self.par = p, deg, par
        self.mons = monomials(deg, par)
        self.idx = {m: i for i, m in enumerate(self.mons)}
        gens = J_generators(Fi, deg, par, p)
        R, C = len(gens), len(self.mons)
        flat = [0] * (R * C)
        for r, g in enumerate(gens):
            for e, v in g.items(): flat[r * C + self.idx[e]] = v
        M = nmod_mat(R, C, flat, p)
        Mr, rank = M.rref()
        self.rank = rank
        piv = []
        for i in range(rank):
            j = 0
            while j < C and int(Mr[i, j]) == 0: j += 1
            piv.append(j)
        self.pivots = piv
        self.pivset = set(piv)
        self.standard = [j for j in range(C) if j not in self.pivset]
        self.std_idx = {j: k for k, j in enumerate(self.standard)}
        # store rref rows (only the rank rows) as dense python lists restricted to standard columns
        self.rows_std = [[int(Mr[i, j]) for j in self.standard] for i in range(rank)]
    def nf_vec(self, poly):
        """normal form as a list over standard monomials (mod p)."""
        p = self.p
        v = [0] * len(self.standard)
        pivcoef = {}
        for e, c in poly.items():
            j = self.idx[e]
            if j in self.pivset:
                pivcoef[j] = (pivcoef.get(j, 0) + c) % p
            else:
                k = self.std_idx[j]; v[k] = (v[k] + c) % p
        for i, j in enumerate(self.pivots):
            c = pivcoef.get(j, 0)
            if c:
                row = self.rows_std[i]
                for k in range(len(v)):
                    if row[k]: v[k] = (v[k] - c * row[k]) % p
        return v
