"""Fast builders for large bilinear systems: C-level multiplication (flint nmod_mpoly) and matrix-based normal forms."""
import sys, time
sys.path.insert(0, '/tmp/hodge/results/off-paired/wave2-excl')
import numpy as np
from flint import nmod_mat, nmod_mpoly_ctx, Ordering
from hodge.core import NV, VARS

def ctx_for(p):
    return nmod_mpoly_ctx.get(VARS, ordering=Ordering.lex, modulus=p)

def to_mp(ctx, P):
    return ctx.from_dict({tuple(int(x) for x in e): int(v) for e, v in P.items()}) if P else ctx.from_dict({})

class MatReducer:
    """Normal forms modulo the row space of `polys` within the monomial list `mons`, applied to many vectors at once:
    coords(V) = V[:, std] - V[:, piv] @ R[:, std]  where R is the rref of the generator matrix."""
    def __init__(self, polys, mons, p, verbose=False):
        t0 = time.time()
        self.p = p; self.mons = mons; self.idx = {m: i for i, m in enumerate(mons)}
        C = len(mons); R = len(polys)
        if R == 0:
            self.pivots = []; self.rank = 0; self.std = list(range(C)); self.Rstd = None
        else:
            flat = [0] * (R * C)
            for r, P in enumerate(polys):
                for e, v in P.items(): flat[r * C + self.idx[e]] = v
            M = nmod_mat(R, C, flat, p)
            Mr, rank = M.rref()
            self.rank = rank
            piv = []
            for i in range(rank):
                j = 0
                while j < C and int(Mr[i, j]) == 0: j += 1
                piv.append(j)
            self.pivots = piv
            ps = set(piv)
            self.std = [j for j in range(C) if j not in ps]
            # R[:, std] as nmod_mat (rank x s)
            self.Rstd = nmod_mat(rank, len(self.std), [int(Mr[i, j]) for i in range(rank) for j in self.std], p) if self.std else None
        self.nstd = len(self.std)
        if verbose: print(f'[MatReducer] gens {R} mons {C} rank {self.rank} std {self.nstd} ({time.time()-t0:.1f}s)', file=sys.stderr, flush=True)

    def coords_block(self, vecs):
        """vecs: list of dict poly (mod p).  Returns list of dict {std_index: value} (normal-form coordinates)."""
        n = len(vecs); p = self.p
        if n == 0: return []
        if self.rank == 0:
            out = []
            sidx = {m: k for k, m in enumerate(self.std)}
            for P in vecs:
                out.append({sidx[self.idx[e]]: v for e, v in P.items() if v})
            return out
        r = self.rank; s = self.nstd
        pivpos = {j: k for k, j in enumerate(self.pivots)}
        stdpos = {j: k for k, j in enumerate(self.std)}
        Vp = [0] * (n * r); Vs = {}
        for i, P in enumerate(vecs):
            for e, v in P.items():
                j = self.idx[e]
                if j in pivpos: Vp[i * r + pivpos[j]] = v
                else: Vs[(i, stdpos[j])] = v
        Mp = nmod_mat(n, r, Vp, p)
        Prod = Mp * self.Rstd   # n x s
        out = []
        for i in range(n):
            row = {}
            for k in range(s):
                v = (Vs.get((i, k), 0) - int(Prod[i, k])) % p
                if v: row[k] = v
            out.append(row)
        return out

def B_products(ctx, Ds, Ns, syz, p, verbose=False):
    """All B_{D_a}(N_b) (a,b) and B_{D_a}(S_k) as mod-p dict polys, via nmod_mpoly.  Returns list of (key, poly)."""
    t0 = time.time()
    Dm = [to_mp(ctx, D) for D in Ds]
    dD = [[D.derivative(i) for i in range(NV)] for D in Dm]
    def prep(N):
        Nm = [to_mp(ctx, N[i]) for i in range(NV)]
        div = sum((Nm[i].derivative(i) for i in range(NV)), ctx.from_dict({}))
        return Nm, div
    NN = [prep(N) for N in Ns]; SS = [prep(S) for S in syz]
    out = []
    nt = len(Ds)
    for a in range(nt):
        for b in range(nt):
            Nm, div = NN[b]
            B = Dm[a] * div
            for i in range(NV):
                if Nm[i]: B = B - Nm[i] * dD[a][i]
            out.append((('t', min(a, b), 't', max(a, b)), B.to_dict()))
        for k in range(len(syz)):
            Nm, div = SS[k]
            B = Dm[a] * div
            for i in range(NV):
                if Nm[i]: B = B - Nm[i] * dD[a][i]
            out.append((('t', a, 's', k), B.to_dict()))
    if verbose: print(f'[B_products] {len(out)} products in {time.time()-t0:.1f}s', file=sys.stderr, flush=True)
    return out
