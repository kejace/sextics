"""Frobenius cycle types of the finite cover f: X -> C (degree 32, projective) and of
X/tau -> C (degree 16) over random F_p-points u of the Fermat cubic C: sum u_i^3 = 0.

Fibre of phi over u: affine scheme Q_i(z) = u_i (64 points), homogenized as
g_i = Q_i(z) - u_i w^2 in F_p[z0..z5, w]; (g_i) is a CI of 6 quadrics in 7 variables with
Hilbert function 1,7,22,42,57,63,64,64,...  Multiplication by a form g of degree 2 on the
64-dim space S_6/I_6 -> S_8/I_8, composed with (mult by w^2)^{-1}, has eigenvalues g(z)/1 at
the 64 affine points.  For g tau-even and (-1)-even (Y-even quadric) the squarefree part of the
characteristic polynomial has degree 16 and its irreducible factor degrees over F_p are the
Frobenius orbit sizes on the 16 points of the fibre of X/tau -> C.  For g Y-odd (tau-odd,
(-1)-even) the squarefree part has degree 32: orbits on the 32 projective points.
Writes one JSON line per sample to w2c_monodromy.jsonl.
"""
import sys, json, time, random, itertools
from flint import nmod_mat, nmod_poly, nmod
P = 1000000007          # p = 2 mod 3: cubing is a bijection, cube roots unique
NV = 7                  # z0..z5, w
PRIMES = [7, 11, 13, 17, 19, 23]
CROSS = [((1,2),(3,4)), ((0,2),(4,5)), ((0,1),(3,5)), ((0,1),(4,5)), ((1,2),(3,5)), ((0,2),(3,4))]

def inv(a): return pow(a, P-2, P)
def monos(deg, nv):
    out = []
    def rec(i, rem, cur):
        if i == nv-1:
            out.append(tuple(cur+[rem])); return
        for e in range(rem, -1, -1):
            rec(i+1, rem-e, cur+[e])
    rec(0, deg, [])
    return out
M6 = monos(6, NV); M8 = monos(8, NV); M4 = monos(4, NV)
IDX8 = {m:i for i,m in enumerate(M8)}; IDX6 = {m:i for i,m in enumerate(M6)}
def quadrics():
    Q = []
    for i in range(6):
        c = (5 * inv(PRIMES[i])) % P
        q = {}
        e = [0]*7; e[i] = 2; q[tuple(e)] = 1
        for (a,b) in CROSS[i]:
            e = [0]*7; e[a] += 1; e[b] += 1; q[tuple(e)] = (q.get(tuple(e),0) + c) % P
        Q.append(q)
    return Q
QS = quadrics()
def polymul_mono(poly, m):
    return {tuple(a+b for a,b in zip(k, m)): v for k,v in poly.items()}
def macaulay(gens, mlist_mult, mlist_cols, idx):
    rows = []
    for g in gens:
        for m in mlist_mult:
            r = [0]*len(mlist_cols)
            for k,v in polymul_mono(g, m).items():
                r[idx[k]] = v % P
            rows.append(r)
    return nmod_mat(rows, P)
def rref_info(M):
    R, rank = M.rref()
    ncols = M.ncols()
    pivots = []; c = 0
    for r in range(rank):
        while c < ncols and int(R[r, c]) == 0: c += 1
        pivots.append(c); c += 1
    piv_rows = nmod_mat([[int(R[r, j]) for j in range(ncols)] for r in range(rank)], P)
    return piv_rows, pivots, rank
def reduce_vectors(V, piv_rows, pivots, ncols):
    # V: nmod_mat (k x ncols); subtract V[:,pivots] * piv_rows  -> supported on non-pivot cols
    k = V.nrows()
    Vp = nmod_mat([[int(V[i, c]) for c in pivots] for i in range(k)], P)
    W = V - Vp * piv_rows
    nonpiv = [c for c in range(ncols) if c not in set(pivots)]
    return nmod_mat([[int(W[i, c]) for c in nonpiv] for i in range(k)], P), nonpiv
def cube_root(w):
    return pow(w, (2*P-1)//3, P)
def sample_point_on_C(rng):
    u = [rng.randrange(P) for _ in range(5)]
    s = (-sum(pow(x,3,P) for x in u)) % P
    u.append(cube_root(s))
    assert sum(pow(x,3,P) for x in u) % P == 0
    return u
def random_quadric(rng, parity):
    # parity 0: Y-even monomials (X_iX_j, Y_iY_j); parity 1: X_i Y_j
    q = {}
    for i in range(6):
        for j in range(i, 6):
            yi = (i>=3)+(j>=3)
            if yi % 2 == parity:
                e = [0]*7; e[i] += 1; e[j] += 1
                q[tuple(e)] = rng.randrange(1, P)
    return q
def cycle_type(M):
    cp = M.charpoly()
    # squarefree part
    d = cp.derivative()
    g = cp.gcd(d)
    sf = cp // g
    fac = sf.factor()
    degs = sorted([f.degree() for f, e in fac[1]], reverse=True)
    mult = sorted(set(e for f, e in fac[1]))
    return degs, sf.degree(), mult
def main(nsamp, seed, outpath, on_C=True):
    rng = random.Random(seed)
    for s in range(nsamp):
        t0 = time.time()
        u = sample_point_on_C(rng) if on_C else [rng.randrange(P) for _ in range(6)]
        gens = []
        for i in range(6):
            g = dict(QS[i]); e = [0]*7; e[6] = 2
            g[tuple(e)] = (g.get(tuple(e), 0) - u[i]) % P
            gens.append(g)
        A8 = macaulay(gens, M6, M8, IDX8); piv8, pivots8, rank8 = rref_info(A8)
        A6 = macaulay(gens, M4, M6, IDX6); piv6, pivots6, rank6 = rref_info(A6)
        nonpiv6 = [c for c in range(len(M6)) if c not in set(pivots6)]
        if rank8 != len(M8) - 64 or rank6 != len(M6) - 64:
            rec = {"sample": s, "u": u, "error": "rank", "rank8": rank8, "rank6": rank6}
            open(outpath, "a").write(json.dumps(rec)+"\n"); continue
        basis6 = [M6[c] for c in nonpiv6]           # 64 standard monomials of degree 6
        w2 = {(0,0,0,0,0,0,2): 1}
        def mult_matrix(gq):
            rows = []
            for m in basis6:
                r = [0]*len(M8)
                for k,v in polymul_mono(gq, m).items(): r[IDX8[k]] = (r[IDX8[k]] + v) % P
                rows.append(r)
            Red, _ = reduce_vectors(nmod_mat(rows, P), piv8, pivots8, len(M8))
            return Red
        Mw = mult_matrix(w2)
        Mwinv = Mw.inv()
        out = {"sample": s, "u": u, "on_C": on_C}
        for parity, name in ((0, "ct16"), (1, "ct32")):
            gq = random_quadric(rng, parity)
            Mg = mult_matrix(gq)
            M = Mwinv * Mg
            degs, sfdeg, mult = cycle_type(M)
            out[name] = degs; out[name+"_sfdeg"] = sfdeg; out[name+"_mult"] = mult
        out["seconds"] = round(time.time()-t0, 2)
        open(outpath, "a").write(json.dumps(out)+"\n")
        print(out, flush=True)
if __name__ == "__main__":
    nsamp = int(sys.argv[1]); seed = int(sys.argv[2]); outpath = sys.argv[3]
    on_C = (sys.argv[4] != "P5") if len(sys.argv) > 4 else True
    main(nsamp, seed, outpath, on_C)
