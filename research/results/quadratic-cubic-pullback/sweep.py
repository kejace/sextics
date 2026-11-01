"""E1 kernel sweep: K_d = {D in S_d^even : A*D in J(F)}, compared with J_d^even.

usage: python3 sweep.py dmin dmax [which]   (which = ours|fermat)
Writes kernel_sweep.jsonl and kernel_basis_d<d>.txt (basis of K_d mod J_d, mod p).
"""
import sys, time, json
sys.path.insert(0, '/tmp/hodge/results/quadratic-cubic-pullback')
import qcp
from qcp import core
from flint import nmod_mat

p = 2147483647
dmin, dmax = int(sys.argv[1]), int(sys.argv[2])
which = sys.argv[3] if len(sys.argv) > 3 else 'ours'
F_ = qcp.FERMAT if which == 'fermat' else qcp.F
A_ = qcp.A
Fp = core.reduce_mod(F_, p); Ap = core.reduce_mod(A_, p)
Fi = [core.pderiv(Fp, i, p) for i in range(6)]
WD = qcp.WD


def sweep(d):
    t0 = time.time()
    cols_list = core.monomials(d + 6, 1)           # Y-odd monomials of degree d+6
    col = {e: j for j, e in enumerate(cols_list)}
    C = len(cols_list)
    # J rows
    Jrows = []
    for i in range(6):
        par = 1 if i < 3 else 0                     # m*F_i must be Y-odd
        for m in core.monomials(d + 1, par):
            Jrows.append(core.shift(Fi[i], m, 1, p))
    Amons = core.monomials(d, 0)
    Arows = [core.shift(Ap, m, 1, p) for m in Amons]
    RJ, RA = len(Jrows), len(Arows)
    flat = [0] * ((RJ + RA) * C)
    for r, poly in enumerate(Jrows + Arows):
        base = r * C
        for e, v in poly.items():
            flat[base + col[e]] = v
    M = nmod_mat(RJ + RA, C, flat, p)
    MJ = nmod_mat(RJ, C, flat[:RJ * C], p)
    t1 = time.time()
    rJ = MJ.rank()
    rM = M.rank()
    t2 = time.time()
    dimK = RA - (rM - rJ)
    dimJd = 3 * len(core.monomials(d - 5)) if d >= 5 else 0
    out = {'d': d, 'which': which, 'cols': C, 'rowsJ': RJ, 'rowsA': RA, 'rankJ': rJ, 'rank_stack': rM,
           'dimK': dimK, 'dimJd_even': dimJd, 'dimK_mod_J': dimK - dimJd,
           'dimR_odd': C - rJ, 'seconds': round(t2 - t0, 1)}
    print(json.dumps(out), flush=True)
    with open(f'{WD}/kernel_sweep.jsonl', 'a') as f:
        f.write(json.dumps(out) + '\n')
    # basis of K_d if nontrivial beyond J_d: left nullspace of stacked matrix restricted to A-part
    if dimK > 0 and dimK - dimJd > 0 and d <= 9:
        # transpose nullspace: vectors (u, c) with u MJ + c MA = 0
        Mt = M.transpose()
        X, nul = Mt.nullspace()
        # extract c-parts (rows RJ..RJ+RA of X), reduce to a basis
        cvecs = []
        for k in range(nul):
            v = [int(X[RJ + r, k]) for r in range(RA)]
            if any(v):
                cvecs.append(v)
        if cvecs:
            Cm = nmod_mat(len(cvecs), RA, [x for v in cvecs for x in v], p)
            Cr, rk = Cm.rref()
            with open(f'{WD}/kernel_basis_d{d}_{which}.txt', 'w') as f:
                for i in range(rk):
                    poly = {}
                    for r in range(RA):
                        v = int(Cr[i, r])
                        if v:
                            poly[Amons[r]] = v if v < p // 2 else v - p
                    f.write(core.to_str(poly) + '\n')
            print(f'  wrote {rk} kernel basis elements (mod p, symmetric lift) to kernel_basis_d{d}_{which}.txt', flush=True)
    return out


for d in range(dmin, dmax + 1):
    sweep(d)
