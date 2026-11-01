"""E1-kernel sweep for both parities of D: kernel of [A]: R_d^{par} -> R_{d+6}^{par+1}, mod p."""
import sys, time, json
sys.path.insert(0, '/tmp/hodge/results/short-power-presentation')
from hlib import *
from flint import nmod_mat
p = core.PRIMES[0]

def jac_cols(F, deg, par, p):
    Fi = [core.pderiv(F, i, p) for i in range(6)]
    cols = []
    for i in range(6):
        pi = 0 if i < 3 else 1
        for m in core.monomials(deg - 5, (par - pi) % 2):
            cols.append(core.shift(Fi[i], m, 1, p))
    return cols

def mat(polys, p, rows):
    R = len(rows); C = len(polys)
    flat = [0] * (R * C)
    for j, P in enumerate(polys):
        for e, v in P.items(): flat[rows[e] * C + j] = v
    return nmod_mat(R, C, flat, p)

def e1_kernel(F, A, d, parD):
    t0 = time.time()
    Fp = core.reduce_mod(F, p); Ap = core.reduce_mod(A, p)
    monsD = core.monomials(d, parD)
    Jcols = jac_cols(Fp, d + 6, (parD + 1) % 2, p)
    rows = {}
    for P in Jcols:
        for e in P: rows.setdefault(e, len(rows))
    AD = [core.shift(Ap, m, 1, p) for m in monsD]
    for P in AD:
        for e in P: rows.setdefault(e, len(rows))
    rJ = mat(Jcols, p, rows).rank()
    rall = mat(Jcols + AD, p, rows).rank()
    kdim = len(monsD) - (rall - rJ)
    JD = jac_cols(Fp, d, parD, p)
    rowsD = {}
    for P in JD:
        for e in P: rowsD.setdefault(e, len(rowsD))
    rJD = mat(JD, p, rowsD).rank() if JD else 0
    out = {'d': d, 'parD': parD, 'dim S_d': len(monsD), 'dim J_d': rJD, 'dim R_d': len(monsD) - rJD,
           'rank J_{d+6}': rJ, 'kernel dim (incl J_d)': kdim, 'kernel mod J': kdim - rJD,
           'rows': len(rows), 'cols': len(Jcols) + len(monsD), 'p': p, 'engine': 2,
           'seconds': round(time.time() - t0, 1)}
    print(json.dumps(out), flush=True)
    with open('/tmp/hodge/results/short-power-presentation/e1kernel.jsonl', 'a') as f:
        f.write(json.dumps(out) + '\n')
    return out

if __name__ == '__main__':
    F, A, q, rec = load()
    jobs = sys.argv[1:]  # like 7e 7o 8e
    for j in jobs:
        d = int(j[:-1]); par = 0 if j[-1] == 'e' else 1
        e1_kernel(F, A, d, par)
