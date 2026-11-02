"""Helpers for the wave-2 exclusion/rigor lens on symmetric-sextic-residual.

Independent (engine-free where it matters) builders for Jacobian matrices mod p, complete-intersection
dimension targets, logging.  Polys are hodge.core dicts {exp tuple: coeff}.
"""
import sys, os, json, time, itertools
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from fractions import Fraction
import numpy as np
from flint import nmod_mat, fmpz_mat, fmpq_mat
from hodge import core

IID = 'symmetric-sextic-residual'
RES = f'/tmp/hodge/results/{IID}'
PK = f'{RES}/wave2-excl'
LOG = f'{RES}/wave2-excl.log.jsonl'
P1 = 2147483647          # wave-1 prime
P2 = 2147483629          # second prime (engine PRIMES[1])
P3 = 1000000007          # 30-bit prime, different size class

def load():
    F, A, q, rec = core.load_instance(IID)
    return F, A, q

def log(**rec):
    rec.setdefault('lens', 'excl')
    rec.setdefault('engine', 2)
    rec['ts'] = time.strftime('%Y-%m-%dT%H:%M:%S')
    with open(LOG, 'a') as f:
        f.write(json.dumps(rec, default=str) + '\n')
    return rec

def receipt(name, obj):
    with open(f'{PK}/{name}.json', 'w') as f:
        json.dump(obj, f, indent=1, default=str)

# ---------- symmetric functions / group ----------
def var(i):
    e = [0] * 6; e[i] = 1
    return {tuple(e): Fraction(1)}
X = [var(0), var(1), var(2)]
Y = [var(3), var(4), var(5)]
def e_k(vs, k):
    out = {}
    for c in itertools.combinations(vs, k):
        t = {(0,) * 6: Fraction(1)}
        for v in c:
            t = core.mul(t, v)
        out = core.add(out, t)
    return out
def p_k(vs, k):
    out = {}
    for v in vs:
        out = core.add(out, core.power(v, k))
    return out
def permute(Pol, perm):
    """variable i -> variable perm[i]"""
    out = {}
    for e, c in Pol.items():
        e2 = [0] * 6
        for i in range(6):
            e2[perm[i]] += e[i]
        out[tuple(e2)] = out.get(tuple(e2), 0) + c
    return {e: c for e, c in out.items() if c}
def signs(Pol, sg):
    out = {}
    for e, c in Pol.items():
        s = 1
        for i in range(6):
            if sg[i] < 0 and e[i] % 2: s = -s
        out[e] = c * s
    return out
def perm_sign(p):
    s = 1
    p = list(p)
    for i in range(len(p)):
        for j in range(i + 1, len(p)):
            if p[i] > p[j]: s = -s
    return s
def G_elements():
    """S3(X) x S3(Y) as 6-permutations, with (sgnX, sgnY)."""
    out = []
    for px in itertools.permutations(range(3)):
        for py in itertools.permutations(range(3)):
            perm = tuple(list(px) + [3 + j for j in py])
            out.append((perm, perm_sign(px), perm_sign(py)))
    return out

# ---------- complete-intersection targets ----------
def ci_ring_dim(k, par=None):
    """dim R_k (Y-parity par) for a smooth sextic: monomials with all exponents <= 4."""
    cnt = 0
    for e in itertools.product(range(5), repeat=6):
        if sum(e) != k: continue
        if par is not None and (e[3] + e[4] + e[5]) % 2 != par: continue
        cnt += 1
    return cnt
def n_mons(k, par=None):
    return len(core.monomials(k, par))
def ci_jac_rank(k, par=None):
    return n_mons(k, par) - ci_ring_dim(k, par)

# ---------- vectorised Jacobian matrix mod p ----------
def mono_codes(mons):
    arr = np.array(mons, dtype=np.int64)
    return (arr * (32 ** np.arange(6, dtype=np.int64))).sum(axis=1)

def poly_to_arrays(P, p):
    E = np.array(list(P.keys()), dtype=np.int64).reshape(-1, 6)
    C = np.array([int(Fraction(c).numerator * pow(Fraction(c).denominator, -1, p)) % p if not isinstance(c, int) else c % p for c in P.values()], dtype=np.int64)
    return E, C

def jac_matrix(F, k, par, p, extra_cols=(), rows=None):
    """nmod_mat with rows = monomials of degree k, parity par; columns = m*F_i (deg m = k-5) then extra_cols.
    Returns (M, mons, ncolsJ)."""
    Fp = core.reduce_mod(F, p)
    Fi = [core.pderiv(Fp, i, p) for i in range(6)]
    mons = core.monomials(k, par) if rows is None else rows
    codes = mono_codes(mons)
    order = np.argsort(codes); scodes = codes[order]
    cols = []   # list of (row_idx array, val array)
    if k >= 5:
        for i in range(6):
            E, C = poly_to_arrays(Fi[i], p)
            if len(E) == 0: continue
            par_m = par if i < 3 else (par + 1) % 2
            for m in core.monomials(k - 5, par_m):
                ec = mono_codes([m])[0]
                cc = (E * (32 ** np.arange(6, dtype=np.int64))).sum(axis=1) + ec
                pos = np.searchsorted(scodes, cc)
                assert np.all(scodes[np.minimum(pos, len(scodes) - 1)] == cc), 'monomial not in row set'
                cols.append((order[pos], C))
    nJ = len(cols)
    for Pc in extra_cols:
        Pp = core.reduce_mod(Pc, p)
        E, C = poly_to_arrays(Pp, p)
        cc = (E * (32 ** np.arange(6, dtype=np.int64))).sum(axis=1)
        pos = np.searchsorted(scodes, cc)
        assert np.all(scodes[np.minimum(pos, len(scodes) - 1)] == cc), 'extra col monomial not in row set'
        cols.append((order[pos], C))
    R, Cn = len(mons), len(cols)
    M = nmod_mat(R, Cn, p)
    for j, (ri, vals) in enumerate(cols):
        for r, v in zip(ri.tolist(), vals.tolist()):
            M[int(r), j] = int(v)
    return M, mons, nJ

def rank_mod(M):
    return M.rank()

def sym_int(v, p):
    v = int(v) % p
    return v if v <= p // 2 else v - p
