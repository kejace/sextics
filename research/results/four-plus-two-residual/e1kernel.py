"""E1 kernel sweep: for each degree d compute dim{D Y-even of degree d : A*D in J(F)} and dim J(F)_d^{even}.
Also extract a basis of the quotient (kernel mod J) when small, saved to e1kernel-basis-d{d}.json."""
import sys, json, time
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core
from hodge.core import monomials, shift, mul, pderiv, reduce_mod
from flint import nmod_mat
F, A, q, rec = core.load_instance('four-plus-two-residual')
p = core.PRIMES[0]
Fp = reduce_mod(F, p); Ap = reduce_mod(A, p)
Fi = [pderiv(Fp, i, p) for i in range(6)]
RES = '/tmp/hodge/results/four-plus-two-residual'
dmax = int(sys.argv[1]) if len(sys.argv) > 1 else 8
dmin = int(sys.argv[2]) if len(sys.argv) > 2 else 1
parD = int(sys.argv[3]) if len(sys.argv) > 3 else 0   # parity of D (0 even, 1 odd)

def build(cols_polys, rowidx=None):
    rows = {} if rowidx is None else rowidx
    for col in cols_polys:
        for e in col:
            rows.setdefault(e, len(rows))
    R, C = len(rows), len(cols_polys)
    flat = [0] * (R * C)
    for j, col in enumerate(cols_polys):
        for e, v in col.items():
            flat[rows[e] * C + j] = v
    return nmod_mat(R, C, flat, p), rows

for d in range(dmin, dmax + 1):
    t0 = time.time()
    monsD = monomials(d, parD)
    # A*D has parity (1 + parD) % 2 ; J columns F_i * m' must have that parity: parity(F_i)= 0 (i<3), 1 (i>=3)
    tpar = (1 + parD) % 2
    Jcols = []
    Jinfo = []
    for i in range(6):
        for m2 in monomials(d + 1, (tpar - (0 if i < 3 else 1)) % 2):
            Jcols.append(shift(Fi[i], m2, 1, p)); Jinfo.append((i, m2))
    Acols = [shift(Ap, m, 1, p) for m in monsD]
    MJ, rows = build(Jcols)
    rJ = MJ.rank()
    MAJ, rows2 = build(Acols + Jcols)
    rAJ = MAJ.rank()
    dimKer = len(monsD) - (rAJ - rJ)
    # dim J_d^{parD}: columns F_i * m'' with deg m'' = d-5, parity so that product has parity parD
    if d >= 5:
        Jd = []
        for i in range(6):
            for m3 in monomials(d - 5, (parD - (0 if i < 3 else 1)) % 2):
                Jd.append(shift(Fi[i], m3, 1, p))
        MJd, _ = build(Jd)
        dimJd = MJd.rank()
    else:
        dimJd = 0
    out = {'d': d, 'parD': parD, 'dim_monsD': len(monsD), 'dim_E1kernel': dimKer, 'dim_Jd': dimJd,
           'dim_kernel_mod_J': dimKer - dimJd, 'rank_J_d+6': rJ, 'rank_AJ': rAJ, 'rows': MAJ.nrows(), 'cols': MAJ.ncols(),
           'seconds': round(time.time() - t0, 1)}
    print(json.dumps(out), flush=True)
    with open(f'{RES}/e1kernel.jsonl', 'a') as f:
        f.write(json.dumps(out) + '\n')
    if dimKer - dimJd > 0 and dimKer - dimJd <= 400:
        # extract kernel basis: nullspace of [A-cols | J-cols], project to A-coordinates, reduce mod J_d
        X, nullity = MAJ.nullspace()
        nA = len(monsD)
        vecs = []
        for k in range(nullity):
            v = [int(X[j, k]) for j in range(nA)]
            if any(v):
                vecs.append(v)
        # row-reduce the projected vectors
        if vecs:
            Mv = nmod_mat(len(vecs), nA, [x for v in vecs for x in v], p)
            Mr, rk = Mv.rref()
            basis = []
            for i in range(rk):
                poly = {monsD[j]: int(Mr[i, j]) for j in range(nA) if int(Mr[i, j])}
                basis.append(poly)
            # now reduce modulo J_d: compute which basis elements are independent mod J_d
            if d >= 5:
                MJd2, rowsJ = build(Jd)
                rJd = MJd2.rank()
                indep = []
                cur = list(Jd)
                currank = rJd
                for b in basis:
                    Mtest, _ = build(cur + [b])
                    r = Mtest.rank()
                    if r > currank:
                        indep.append(b); cur.append(b); currank = r
                basis = indep
            with open(f'{RES}/e1kernel-basis-d{d}-par{parD}.json', 'w') as f:
                json.dump([{str(k): v for k, v in b.items()} for b in basis], f)
            print(f'  saved {len(basis)} kernel-mod-J basis elements for d={d}', flush=True)
