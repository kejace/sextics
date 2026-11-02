"""Wave-2 exclusion-lens helpers for one-reflection-skew (independent of wave-1 scripts).
Logging: one JSON line per experiment to wave2-excl.log.jsonl."""
import sys, time, json, os, itertools
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from fractions import Fraction
from hodge import core
from hodge.core import mul, add, sub, scale, pderiv, shift, monomials
from flint import nmod_mat, fmpq_mat, fmpq, fmpz_mat
INST = 'one-reflection-skew'
WD = f'/tmp/hodge/results/{INST}'
LOG = f'{WD}/wave2-excl.log.jsonl'
PK = f'{WD}/wave2-excl'
F, A, q, rec = core.load_instance(INST)
P1, P2, P3 = core.PRIMES[0], core.PRIMES[1], core.PRIMES[2]
P_SMALL = 32003
X = [core.parse(v) for v in core.VARS]
Fi = [pderiv(F, i) for i in range(6)]
f5 = {e: c for e, c in F.items() if e[5] == 0}     # f = F - Y2^6
fi5 = [pderiv(f5, i) for i in range(5)]

def log(name, **kw):
    rec_ = {'name': name, 'time': time.strftime('%Y-%m-%dT%H:%M:%S'), 'lens': 'excl'}
    rec_.update(kw)
    with open(LOG, 'a') as fh:
        fh.write(json.dumps(rec_, default=str) + '\n')
    return rec_

def receipt(fname, obj):
    with open(os.path.join(PK, fname), 'w') as fh:
        json.dump(obj, fh, indent=1, default=str)

def hilbert_ci(nvars, gens_deg, ngens, maxdeg):
    """Hilbert function of a complete intersection of ngens forms of degree gens_deg in nvars variables."""
    # series ((1 - t^g)^n) / (1-t)^nvars
    import numpy as np
    ser = [0] * (maxdeg + 1); ser[0] = 1
    for _ in range(ngens):
        new = ser[:]
        for k in range(gens_deg, maxdeg + 1):
            new[k] -= ser[k - gens_deg]
        ser = new
    for _ in range(nvars):
        for k in range(1, maxdeg + 1):
            ser[k] += ser[k - 1]
    return ser

def vec_index(polys):
    rows = {}
    for P in polys:
        for e in P:
            rows.setdefault(e, len(rows))
    return rows

def nmod_from_cols(cols, p, rows=None):
    rows = rows if rows is not None else vec_index(cols)
    R, C = len(rows), len(cols)
    flat = [0] * (R * C)
    for j, col in enumerate(cols):
        for e, v in col.items():
            flat[rows[e] * C + j] = int(v) % p
    return nmod_mat(R, C, flat, p), rows

def modp(P, p):
    return core.reduce_mod(P, p)

def polys_from_rref_rows(Mr, rank, mons):
    out = []
    for i in range(rank):
        poly = {}
        for j, m in enumerate(mons):
            v = int(Mr[i, j])
            if v: poly[m] = v
        out.append(poly)
    return out
