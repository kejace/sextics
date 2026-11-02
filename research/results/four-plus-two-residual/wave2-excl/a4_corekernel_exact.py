"""(a/b) EXACT (characteristic 0) core E1 kernels: ker(a : R(H)_k -> R(H)_{k+6}), both Y-parities, k = 0..kmax.
Integer matrix [A-cols | J(H)-cols] (H scaled by 7), exact nullspace with flint fmpz_mat; D-parts reduced modulo J(H)_k exactly.
Writes corekernel-exact.json: for each (k,par): dim kernel, dim J_k, dim kernel/J, and an exact rational basis of a complement
of J(H)_k inside the kernel (polys with Fraction coefficients, as strings)."""
import sys, json; sys.path.insert(0, '/tmp/hodge/results/four-plus-two-residual/wave2-excl')
from lib import *
from flint import fmpz_mat, fmpq_mat, fmpq
kmax = int(sys.argv[1]) if len(sys.argv) > 1 else 12
H7 = {e: int(c * 7) for e, c in H.items()}
assert all(Fraction(c) * 7 == int(c * 7) for c in H.values())
Hi = [pderiv(H7, i) for i in CORE]
a = {e: int(c) for e, c in A.items()}
def ibuild(cols):
    rows = {}
    for col in cols:
        for e in col: rows.setdefault(e, len(rows))
    R, C = len(rows), len(cols)
    M = [[0] * C for _ in range(R)]
    for j, col in enumerate(cols):
        for e, v in col.items(): M[rows[e]][j] = int(v)
    return fmpz_mat(M) if R and C else None
out = {}
for k in range(0, kmax + 1):
    for par in (0, 1):
        t0 = time.time()
        mons = cmon(k, par); tpar = (1 + par) % 2
        Jcols = []
        for idx, i in enumerate(CORE):
            for m2 in cmon(k + 1, (tpar - (0 if i < 3 else 1)) % 2):
                Jcols.append(shift(Hi[idx], m2, 1))
        Acols = [shift(a, m, 1) for m in mons]
        M = ibuild(Acols + Jcols)
        rAJ = M.rank()
        MJ = ibuild(Jcols); rJ = MJ.rank() if MJ is not None else 0
        X, nullity = M.nullspace()
        nA = len(mons)
        vecs = []
        for kk in range(nullity):
            v = [int(X[j, kk]) for j in range(nA)]
            if any(v): vecs.append(v)
        kernel_dim = len(mons) - (rAJ - rJ)
        # kernel as a Q-space: row space of vecs (projected nullspace); rank should equal kernel_dim
        Kmat = fmpz_mat(vecs) if vecs else None
        rk_proj = Kmat.rank() if Kmat is not None else 0
        assert rk_proj == kernel_dim, (k, par, rk_proj, kernel_dim)
        # J(H)_k basis
        Jk = []
        if k >= 5:
            for idx, i in enumerate(CORE):
                for m3 in cmon(k - 5, (par - (0 if i < 3 else 1)) % 2):
                    Jk.append(shift(Hi[idx], m3, 1))
        idx = {m: i for i, m in enumerate(mons)}
        Jrows = [[0] * nA for _ in Jk]
        for r, P in zip(Jrows, Jk):
            for e, v in P.items(): r[idx[e]] = int(v)
        rJk = fmpz_mat(Jrows).rank() if Jrows else 0
        # complement of J_k inside kernel: greedy over rref rows of kernel
        indep = []
        cur = list(Jrows); currank = rJk
        if vecs:
            Kq = fmpq_mat(Kmat); Kr, rk = Kq.rref()
            for i in range(rk):
                row = [Kr[i, j] for j in range(nA)]
                # scale to integers
                from math import lcm
                den = 1
                for x in row: den = lcm(den, int(x.q))
                irow = [int(x.p) * (den // int(x.q)) for x in row]
                test = fmpz_mat(cur + [irow]).rank()
                if test > currank:
                    indep.append(irow); cur.append(irow); currank += 1
        rec_ = dict(k=k, par=par, n_mons=nA, kernel_dim=kernel_dim, dim_J_k=rJk, kernel_mod_J=kernel_dim - rJk,
                    rank_J_k6=rJ, dim_J_k6_hilbert=dim_Jcore(k + 6, tpar), dim_J_k_hilbert=dim_Jcore(k, par),
                    ring_dim=dim_RH(k, par), seconds=round(time.time() - t0, 2), exact=True)
        assert rJ == rec_['dim_J_k6_hilbert'] and rJk == rec_['dim_J_k_hilbert']
        assert len(indep) == kernel_dim - rJk
        basis = []
        for irow in indep:
            P = {mons[j]: irow[j] for j in range(nA) if irow[j]}
            basis.append(core.to_str(P))
        rec_['basis'] = basis
        rec_['max_abs_coeff'] = max([abs(x) for irow in indep for x in irow] + [0])
        out[f'{k},{par}'] = rec_
        log(dict(exp='a4-corekernel-exact', **{kk: vv for kk, vv in rec_.items() if kk != 'basis'}))
json.dump(out, open('corekernel-exact.json', 'w'), indent=1)
