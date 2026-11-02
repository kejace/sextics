"""wave-4 T-reduction E1 library (short-power-presentation).
J(F) = I' + (F0, F3), I' = (X1^5-2X0^5, X2^5-3X0^5, Y1^5-3Y0^5, Y2^5-5Y0^5)  (F1-2F0 = 6(X1^5-2X0^5), ...).
I' is a Groebner basis (pure-power leading terms), T = S/I' is free over C[X0,Y0] with basis X1^a X2^b Y1^c Y2^e, a..e<=4.
A*D in J  <=>  reduce_T(A*D) in span{F0*m, F3*m : m in T_{n-5}} inside T_n (n = d+6), per Y-parity sector."""
import sys, json, time, itertools
sys.path.insert(0, '/home/kejace/devel/sextics-work')
sys.path.insert(0, '/tmp/hodge/results/short-power-presentation')
from hodge import core, jacobian
from hodge.core import monomials, reduce_mod, pderiv, shift
from flint import nmod_mat
from fractions import Fraction
from math import comb

IID = 'short-power-presentation'
WAVE4 = '/tmp/hodge/results/short-power-presentation/wave4'
RED = {1: (0, 2), 2: (0, 3), 4: (3, 3), 5: (3, 5)}   # var -> (target var, coefficient): X1^5 -> 2 X0^5, ...

def log(rec):
    rec = dict(rec); rec.setdefault('time', time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()))
    with open(WAVE4 + '/log.jsonl', 'a') as f: f.write(json.dumps(rec) + '\n')

def reduce_T(P, p):
    out = {}
    for e, v in P.items():
        e = list(e); coef = v % p
        for i, (t, cf) in RED.items():
            k, r = divmod(e[i], 5)
            if k:
                coef = coef * pow(cf, k, p) % p
                e[t] += 5 * k; e[i] = r
        e = tuple(e)
        out[e] = (out.get(e, 0) + coef) % p
    return {e: v for e, v in out.items() if v}

def tbasis(n, par):
    out = []
    if n < 0: return out
    for e1, e2, e4, e5 in itertools.product(range(5), repeat=4):
        rem = n - e1 - e2 - e4 - e5
        if rem < 0: continue
        for e3 in range(rem + 1):
            if ((e3 + e4 + e5) & 1) != par: continue
            out.append((rem - e3, e1, e2, e3, e4, e5))
    return out

def dimS(k, par):
    if k < 0: return 0
    return sum(comb(j + 2, 2) * comb(k - j + 2, 2) for j in range(k + 1) if (j & 1) == par)

def dimR(n, par):
    """complete-intersection prediction for dim (S/J)_n^par (F0,F1,F2 Y-even, F3,F4,F5 Y-odd generators of degree 5)."""
    return sum((-1) ** (i + j) * comb(3, i) * comb(3, j) * dimS(n - 5 * (i + j), (par - j) % 2)
               for i in range(4) for j in range(4))

def dimT(n, par):
    return sum((-1) ** (i + j) * comb(2, i) * comb(2, j) * dimS(n - 5 * (i + j), (par - j) % 2)
               for i in range(3) for j in range(3))

class TBlock:
    """rref of [ (F0,F3)-block of T_n^par | family columns ] mod p; E1 kernel = family vectors landing in the block span."""
    def __init__(self, F, n, par, fam_polys, p, verbose=False):
        t0 = time.time()
        Fp = reduce_mod(F, p)
        F0 = reduce_T(pderiv(Fp, 0, p), p); F3 = reduce_T(pderiv(Fp, 3, p), p)
        rows = tbasis(n, par); idx = {e: i for i, e in enumerate(rows)}
        cols = [reduce_T(shift(F0, m, 1, p), p) for m in tbasis(n - 5, par)]
        cols += [reduce_T(shift(F3, m, 1, p), p) for m in tbasis(n - 5, 1 - par)]
        CJ = len(cols)
        for P in fam_polys:
            Q = reduce_T(reduce_mod(P, p), p)
            for e in Q:
                if e not in idx: raise ValueError('family polynomial outside T_%d^%d: %s' % (n, par, e))
            cols.append(Q)
        R, C = len(rows), len(cols)
        flat = [0] * (R * C)
        for j, col in enumerate(cols):
            for e, v in col.items(): flat[idx[e] * C + j] = v
        t1 = time.time()
        Mr, rank = nmod_mat(R, C, flat, p).rref()
        t2 = time.time()
        pivots = []; j = 0
        for i in range(rank):
            while j < C and int(Mr[i, j]) == 0: j += 1
            pivots.append(j); j += 1
        rank_J = sum(1 for j in pivots if j < CJ)
        fam_piv = [j for j in pivots if j >= CJ]; fps = set(fam_piv)
        pivrow = {j: i for i, j in enumerate(pivots)}
        vecs = []
        for j in range(CJ, C):
            if j in fps: continue
            v = {j - CJ: 1}
            for jp in fam_piv:
                x = int(Mr[pivrow[jp], j])
                if x: v[jp - CJ] = (-x) % p
            vecs.append(v)
        self.vecs = vecs; self.p = p
        K = len(fam_polys)
        self.info = dict(n=n, par=par, p=p, rows=R, cols=C, CJ=CJ, K=K, rank_total=rank, rank_J=rank_J,
                         dimT=R, dimT_pred=dimT(n, par), ring_dim=R - rank_J, ring_dim_pred=dimR(n, par),
                         hilbert_ok=bool(R - rank_J == dimR(n, par) and R == dimT(n, par)),
                         kernel_dim=len(vecs), t_build=round(t1 - t0, 1), t_rref=round(t2 - t1, 1))
        if verbose: print(json.dumps(self.info), flush=True)

def e1_kernel(F, A, fam, d, parD, p, verbose=False):
    """E1 kernel of {D in span(fam)} (fam: list of (name, poly) of degree d, Y-parity parD) mod p.
    Returns (vectors as {name: coeff}, info)."""
    parA = core.parity_of(A); par_t = (parA + parD) % 2
    polys = [core.mul(reduce_mod(A, p), reduce_mod(P, p), p) for _, P in fam]
    B = TBlock(F, d + 6, par_t, polys, p, verbose)
    names = [nm for nm, _ in fam]
    vecs = [{names[k]: c for k, c in v.items()} for v in B.vecs]
    return vecs, B.info

def in_J_family(F, fam, d, parD, p):
    """{D in span(fam) : D in J_d} mod p (same machinery at degree d)."""
    polys = [reduce_mod(P, p) for _, P in fam]
    B = TBlock(F, d, parD, polys, p)
    names = [nm for nm, _ in fam]
    return [{names[k]: c for k, c in v.items()} for v in B.vecs], B.info

def vec_to_poly(vec, fam, p=None):
    name2P = {n: P for n, P in fam}
    out = {}
    for n, cf in vec.items():
        out = core.add(out, core.scale(name2P[n], cf, p), p)
    return out

def sym(v, p):
    return {k: (x if x <= p // 2 else x - p) for k, x in v.items()}
