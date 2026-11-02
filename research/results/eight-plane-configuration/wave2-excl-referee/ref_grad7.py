"""d = 7 gradient subfamily: w = grad G, G cubic of parity par (28 monomials). Hess_F(grad m_a, grad m_b) reduced mod J_8^+ (R_8^+ = 483).
If the 406 symmetric pairs are independent in R_8^+ (mod p, with rank_p J_8^+ = 168 exact), then Hess(w,w) in J_8 forces g_a g_b = 0 => G = 0 over Q."""
import sys, time
sys.path.insert(0, '/tmp/hodge/results/eight-plane-configuration/wave2-excl-referee')
import ref_lib as L
from flint import nmod_mat
p = int(sys.argv[1]) if len(sys.argv) > 1 else 1000000009
F, A, q = L.load(); Fp = L.red(F, p)
Fi = [L.pderiv(Fp, i, p) for i in range(L.NV)]
Fij = [[L.pderiv(Fi[i], j, p) for j in range(L.NV)] for i in range(L.NV)]
mons8 = L.monomials(8, 0); idx = {m: i for i, m in enumerate(mons8)}
gens = L.J_generators(Fi, 8, 0, p)
for par in (0, 1):
    t0 = time.time()
    cub = L.monomials(3, par)
    grads = [[L.pderiv({m: 1}, i, p) for i in range(L.NV)] for m in cub]
    pairs = [(a, b) for a in range(len(cub)) for b in range(a, len(cub))]
    M = nmod_mat(len(gens) + len(pairs), len(mons8), p)
    for r, P in enumerate(gens):
        for e, v in P.items(): M[r, idx[e]] = v
    for r, (a, b) in enumerate(pairs):
        h = {}
        for i in range(L.NV):
            for j in range(L.NV):
                if grads[a][i] and grads[b][j]:
                    h = L.add(h, L.mul(L.mul(grads[a][i], grads[b][j], p), Fij[i][j], p), p)
        for e, v in h.items(): M[len(gens) + r, idx[e]] = v
    rfull = M.rank()
    for r in range(len(gens), len(gens) + len(pairs)):
        for c in range(len(mons8)): pass
    # rank of J part: recompute on a J-only matrix (small)
    MJ = nmod_mat(len(gens), len(mons8), p)
    for r, P in enumerate(gens):
        for e, v in P.items(): MJ[r, idx[e]] = v
    rJ = MJ.rank()
    L.log({'check': f'd=7 gradient subfamily, G parity {par}', 'p': p, 'n_cubics': len(cub), 'dimSym2': len(pairs), 'rank_J8_p': rJ, 'exact_J8': L.J_dim(8, 0),
           'rank_pairs_mod_J8': rfull - rJ, 'full_rank': rfull - rJ == len(pairs), 'seconds': round(time.time() - t0, 1),
           'verdict': 'only G = 0: no degree-7 denominator grad(F).grad(G) (rigorous over C)' if (rfull - rJ == len(pairs) and rJ == L.J_dim(8, 0)) else 'FAIL'})
