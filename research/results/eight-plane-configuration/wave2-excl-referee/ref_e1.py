"""E1 kernel {D in S_d^par : A*D in J(F)_{d+6}} at a NEW prime, own code (ranks).
dim_p V = nD - (rank_p[A m | J] - rank_p J).  Rigor over Q: dim_Q V <= dim_p V when rank_p J_{d+6} = dim_Q J_{d+6} (exact Koszul value);
V contains J_d (dim 3|S_{d-5}| for d <= 9), so kermodJ_p = 0 proves V = J_d over Q (hence over C)."""
import sys, time, os
sys.path.insert(0, '/tmp/hodge/results/eight-plane-configuration/wave2-excl-referee')
import ref_lib as L
from flint import nmod_mat

p = int(sys.argv[1]) if len(sys.argv) > 1 else 1000000009
ds = [int(x) for x in sys.argv[2:]] or [4, 5, 6, 7, 8, 9]
F, A, q = L.load()
Fp = L.red(F, p); Ap = L.red(A, p)
Fi = [L.pderiv(Fp, i, p) for i in range(L.NV)]

def avail_gb():
    for line in open('/proc/meminfo'):
        if line.startswith('MemAvailable'): return int(line.split()[1]) / 1048576

for d in ds:
    for par in (0, 1):
        parp = (par + 1) % 2
        mons = L.monomials(d + 6, parp); idx = {m: i for i, m in enumerate(mons)}
        gens = L.J_generators(Fi, d + 6, parp, p)
        Dm = L.monomials(d, par)
        need_gb = (len(gens) + len(Dm)) * len(mons) * 8 / 1e9 * 2.2
        if avail_gb() < need_gb + 1.5:
            L.log({'check': f'e1kernel d={d} par={par}', 'p': p, 'skipped': True, 'reason': f'memory: need ~{need_gb:.1f} GB, available {avail_gb():.1f} GB'}); continue
        t0 = time.time()
        M = nmod_mat(len(gens) + len(Dm), len(mons), p)
        for r, P in enumerate(gens):
            for e, v in P.items(): M[r, idx[e]] = v
        for r, m in enumerate(Dm):
            for e, v in L.shift(Ap, m, 1, p).items(): M[len(gens) + r, idx[e]] = v
        t1 = time.time()
        rank_full = M.rank()
        for r, m in enumerate(Dm):
            for e, v in L.shift(Ap, m, 1, p).items(): M[len(gens) + r, idx[e]] = 0
        rank_J = M.rank()
        del M
        nD = len(Dm)
        dimV = nD - (rank_full - rank_J)
        exactJ = L.J_dim(d + 6, parp); dimJd = L.J_dim(d, par)
        forced = max(0, L.R_dim(d, par) - L.R_dim(d + 6, parp))
        L.log({'check': f'e1kernel d={d} par={par}', 'p': p, 'rows': len(gens) + nD, 'cols': len(mons), 'nD': nD, 'rank_full': rank_full,
               'rank_J_p': rank_J, 'exact_J_dim': exactJ, 'nodrop': rank_J == exactJ, 'dimV_p': dimV, 'dimJ_d_exact': dimJd,
               'kermodJ_p': dimV - dimJd, 'forced_lower_bound': forced, 't_build': round(t1 - t0, 1), 't_ranks': round(time.time() - t1, 1),
               'verdict': ('V = J_d over Q (rigorous)' if (dimV == dimJd and rank_J == exactJ) else ('kernel mod J = %d (forced >= %d)%s' % (dimV - dimJd, forced, '' if rank_J == exactJ else ' NO-DROP FAILED')))})
