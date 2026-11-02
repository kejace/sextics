"""Threefold reduction: Y_H : F' = H(X0,X1,Y0,Y1) + U^6 (U = X2 slot, Y2 unused), target A' = a*U, a = X0^3 Y0^3,
core-only denominator D(X0,X1,Y0,Y1). Fixed-D linear system mod p in 5 variables, graded by U-degree mod 6.
Also: 4-variable E1 kernel (a*D in J(H)) basis at degree d, with multi-prime rational lift."""
import sys, json, time, random
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core
from hodge.core import monomials, shift, pderiv, reduce_mod, mul, add, scale, ratrecon, crt_pair
from flint import nmod_mat
from fractions import Fraction
F, A, q, rec = core.load_instance('four-plus-two-residual')
H = {e: c for e, c in F.items() if e[2] == 0 and e[5] == 0}
U6 = {(0, 0, 6, 0, 0, 0): Fraction(1)}
Fp_ = core.add(H, U6)                     # F' = H + U^6
Ap_ = {(3, 0, 1, 3, 0, 0): Fraction(1)}   # a*U
CORE = [0, 1, 3, 4]
def mons5(deg, umod=None):
    out = []
    for e in monomials(deg):
        if e[5]: continue
        if umod is not None and e[2] % 6 != umod: continue
        out.append(e)
    return out
def cmon(deg, par=None):
    return [e for e in monomials(deg, par) if e[2] == 0 and e[5] == 0]
def build(cols, p, rows=None):
    rows = {} if rows is None else rows
    for col in cols:
        for e in col: rows.setdefault(e, len(rows))
    R, C = len(rows), len(cols)
    flat = [0] * (R * C)
    for j, col in enumerate(cols):
        for e, v in col.items(): flat[rows[e] * C + j] = v % p
    return nmod_mat(R, C, flat, p), rows
def reduced_test(Hcore, a, core_idx, D, p, uidx=None):
    """Reduced model: F' = Hcore + U^6 (U at index uidx, not in core), A' = a U^(c-1), c = -(deg a + n_core) mod 6.
    Fixed-D linear system mod p restricted to the U-degree sector. Returns dict."""
    t0 = time.time()
    n_core = len(core_idx)
    if uidx is None:
        uidx = [i for i in range(6) if i not in core_idx][0]
    c = (-(core.degree(a) + n_core)) % 6
    assert 1 <= c <= 5
    allowed = set(core_idx) | {uidx}
    def mons_sec(deg, umod):
        return [e for e in monomials(deg) if all(e[j] == 0 for j in range(6) if j not in allowed) and e[uidx] % 6 == umod]
    U6 = {tuple(6 if j == uidx else 0 for j in range(6)): 1}
    Uc1 = {tuple(c - 1 if j == uidx else 0 for j in range(6)): 1}
    Fp = add(reduce_mod(Hcore, p), U6, p); Ap = mul(reduce_mod(a, p), Uc1, p); Dp = reduce_mod(D, p)
    d = core.degree(D); degA = core.degree(Ap)
    degN = degA + d - 5; degQ = degA + 2 * d - 12
    Fi = [pderiv(Fp, i, p) for i in range(6)]; Di = [pderiv(Dp, i, p) for i in range(6)]
    FD = mul(Fp, Dp, p); FF = mul(Fp, Fp, p)
    cols = []
    for i in list(core_idx) + [uidx]:
        umod = (c - 1) % 6 if i != uidx else c % 6
        G = add(mul(Dp, Fi[i], p), mul(Fp, Di[i], p), p)
        for m in mons_sec(degN, umod):
            dm = pderiv({m: 1}, i, p)
            cols.append(add(scale(shift(G, m, 1, p), p - 1, p), mul(FD, dm, p), p))
    if degQ >= 0:
        for m in mons_sec(degQ, (c - 1) % 6):
            cols.append(scale(shift(FF, m, 1, p), p - 1, p))
    rhs = scale(mul(Ap, mul(Dp, Dp, p), p), p - 1, p)
    M, rows = build(cols, p); r1 = M.rank()
    M2, _ = build(cols + [rhs], p, dict(rows)); r2 = M2.rank()
    return dict(consistent=(r1 == r2), rank=r1, nullity=len(cols) - r1, rows=M2.nrows(), cols=len(cols), seconds=round(time.time() - t0, 2), d=d, c=c, degN=degN, degQ=degQ)
def threefold_test(D, p):
    return reduced_test(H, A, CORE, D, p, uidx=2)
def e1_kernel_core(k, par, p):
    """basis (rref, mod p) of {D core deg k parity par : a D in J(H)} modulo J(H)_k."""
    Hp = reduce_mod(H, p); ap = reduce_mod(A, p)
    Hi = [pderiv(Hp, i, p) for i in CORE]
    mons = cmon(k, par)
    tpar = (1 + par) % 2
    Jcols = []
    for idx, i in enumerate(CORE):
        for m2 in cmon(k + 1, (tpar - (0 if i < 3 else 1)) % 2):
            Jcols.append(shift(Hi[idx], m2, 1, p))
    Acols = [shift(ap, m, 1, p) for m in mons]
    MAJ, _ = build(Acols + Jcols, p)
    X, nullity = MAJ.nullspace()
    nA = len(mons)
    vecs = [[int(X[j, kk]) for j in range(nA)] for kk in range(nullity)]
    vecs = [v for v in vecs if any(v)]
    if not vecs: return [], mons
    Mv = nmod_mat(len(vecs), nA, [x for v in vecs for x in v], p)
    Mr, rk = Mv.rref()
    basis = [{mons[j]: int(Mr[i, j]) for j in range(nA) if int(Mr[i, j])} for i in range(rk)]
    Jk = []
    if k >= 5:
        for idx, i in enumerate(CORE):
            for m3 in cmon(k - 5, (par - (0 if i < 3 else 1)) % 2):
                Jk.append(shift(Hi[idx], m3, 1, p))
    rJ = build(Jk, p)[0].rank() if Jk else 0
    indep = []; cur = list(Jk); currank = rJ
    for b in basis:
        Mt, _ = build(cur + [b], p)
        if Mt.rank() > currank:
            indep.append(b); cur.append(b); currank += 1
    return indep, mons
def lift_basis(k, par, nprimes=6):
    """rational reconstruction of the rref kernel basis using CRT over several primes."""
    primes = core.PRIMES[:nprimes]
    bases = [e1_kernel_core(k, par, p)[0] for p in primes]
    n = len(bases[0])
    if any(len(b) != n for b in bases): return None
    out = []
    for i in range(n):
        keys = set(bases[0][i])
        if any(set(b[i]) != keys for b in bases): return None
        poly = {}
        for e in keys:
            a, m = bases[0][i][e] % primes[0], primes[0]
            for b, p in zip(bases[1:], primes[1:]):
                a, m = crt_pair(a, m, b[i][e] % p, p)
            r = ratrecon(a, m)
            if r is None: return None
            poly[e] = r
        out.append(poly)
    return out
LOG = '/tmp/hodge/results/four-plus-two-residual/log.jsonl'
def log(ent):
    with open(LOG, 'a') as f: f.write(json.dumps(ent) + '\n')
def run3(D, note, family, p=None, extra=''):
    p = p or core.PRIMES[0]
    r = threefold_test(D, p)
    ent = {'D': note, 'd': r['d'], 'family': family, 'e1': True, 'consistent': r['consistent'],
           'rank': r['rank'], 'nullity': r['nullity'], 'rows': r['rows'], 'cols': r['cols'], 'seconds': r['seconds'],
           'note': 'threefold model Y_H: U^6+H, target aU, mod p only' + extra, 'engine': 'threefold-v1', 'p': p}
    log(ent)
    print({k: ent[k] for k in ('D', 'consistent', 'rank', 'nullity', 'rows', 'cols', 'seconds')}, flush=True)
    if r['consistent']:
        print('!!!!! CONSISTENT in threefold model:', note, flush=True)
        with open('/tmp/hodge/results/four-plus-two-residual/threefold-hits.jsonl', 'a') as f:
            f.write(json.dumps({'note': note, 'D': {str(k): str(v) for k, v in D.items()}}) + '\n')
    return ent
if __name__ == '__main__':
    p = core.PRIMES[0]
    random.seed(7)
    dmin = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    dmax = int(sys.argv[2]) if len(sys.argv) > 2 else 8
    nbasis = int(sys.argv[3]) if len(sys.argv) > 3 else 1000
    for d in range(dmin, dmax + 1):
        for par in (0, 1):
            t0 = time.time()
            basis, mons = e1_kernel_core(d, par, p)
            print(f'd={d} par={par}: core E1-kernel mod J(H) dim {len(basis)} ({round(time.time()-t0,1)}s)', flush=True)
            if not basis: continue
            tests = []
            for r in range(3):
                D = {}
                for b in basis: D = add(D, scale(b, random.randrange(1, p), p), p)
                tests.append((D, f'random combo #{r} of {len(basis)} core E1-kernel elements, d={d} par={par}'))
            for i, b in enumerate(basis[:nbasis]):
                tests.append((b, f'core E1-kernel rref basis #{i}, d={d} par={par}, {len(b)} terms'))
            for D, note in tests:
                run3(D, note, 'threefold-reduction-core-kernel', p)
