"""Shared helpers for the wave-2 exclusion lens on reflection-coupled-residual."""
import sys, json, time, os, itertools
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core, jacobian
from hodge.core import NV, monomials, shift, pderiv, mul, add, sub, scale, reduce_mod
from flint import nmod_mat
LANE = 'reflection-coupled-residual'
RES = f'/tmp/hodge/results/{LANE}'
PK = f'{RES}/wave2-excl'
LOG = f'{RES}/wave2-excl.log.jsonl'
NEWPRIMES = [1000000007, 998244353, 1000000009, 1073741827]

def log(rec):
    rec = dict(rec); rec['ts'] = time.strftime('%Y-%m-%dT%H:%M:%S'); rec['lens'] = 'excl'; rec['engine'] = 3
    with open(LOG, 'a') as f: f.write(json.dumps(rec) + '\n')
    print(json.dumps(rec)[:400], flush=True)

def load():
    return core.load_instance(LANE)

def expected_ring_dim(k, parity=None):
    """dim (S/J)_k for a smooth sextic (Hilbert series (1+t+..+t^4)^6), optionally by Y-parity."""
    n = 0
    for e in itertools.product(range(5), repeat=6):
        if sum(e) == k and (parity is None or (e[3] + e[4] + e[5]) % 2 == parity):
            n += 1
    return n

def expected_ring_dim_klein(k):
    """dict {(y2par, y01par): dim (S/J)_k in that monomial character} (Fermat basis; valid for smooth F with Klein symmetry)."""
    out = {}
    for e in itertools.product(range(5), repeat=6):
        if sum(e) == k:
            key = (e[5] % 2, (e[3] + e[4]) % 2)
            out[key] = out.get(key, 0) + 1
    return out

def klein(e):
    return (e[5] % 2, (e[3] + e[4]) % 2)

def poly_to_vec(P, idx, C, p):
    v = [0] * C
    for e, c in P.items():
        v[idx[e]] = c % p
    return v

class Reducer:
    """rref of J(F)_deg (one Y-parity, or all if parity None) mod p; vectorized normal forms."""
    def __init__(self, F, deg, p, parity=None, verbose=False):
        t0 = time.time()
        self.p, self.deg, self.parity = p, deg, parity
        self.mons = monomials(deg, parity)
        self.idx = {e: i for i, e in enumerate(self.mons)}
        C = len(self.mons)
        Fp = reduce_mod(F, p)
        Fi = [pderiv(Fp, i, p) for i in range(NV)]
        rows = []
        if deg >= 5:
            for i in range(NV):
                pi = 0 if i < 3 else 1
                for m in monomials(deg - 5, None if parity is None else (parity - pi) % 2):
                    rows.append(shift(Fi[i], m, 1, p))
        R = len(rows)
        if R == 0:
            self.rank = 0; self.pivots = []; self.Jr = None
        else:
            flat = [0] * (R * C)
            for r, poly in enumerate(rows):
                base = r * C
                for e, v in poly.items(): flat[base + self.idx[e]] = v
            M = nmod_mat(R, C, flat, p)
            Jr, rank = M.rref()
            self.rank = rank
            piv = []
            for i in range(rank):
                j = 0
                while int(Jr[i, j]) == 0: j += 1
                piv.append(j)
            self.pivots = piv
            # keep only the rank rows
            self.Jr = nmod_mat(rank, C, [int(Jr[i, j]) for i in range(rank) for j in range(C)], p) if rank else None
        self.pivset = set(self.pivots)
        self.standard = [j for j in range(C) if j not in self.pivset]
        self.ring_dim = C - self.rank
        self.t = time.time() - t0
        if verbose:
            print(f'[Reducer] deg={deg} par={parity} mons={C} Jrows={R} rank={self.rank} ring={self.ring_dim} ({self.t:.1f}s)', flush=True)

    def nf_rows(self, polys):
        """Normal forms of a list of polys (same degree/parity) -> nmod_mat (len(polys) x C)."""
        p, C = self.p, len(self.mons)
        n = len(polys)
        flat = [0] * (n * C)
        for r, P in enumerate(polys):
            base = r * C
            for e, v in P.items(): flat[base + self.idx[e]] = v % p
        V = nmod_mat(n, C, flat, p)
        if self.rank == 0:
            return V
        k = self.rank
        flatp = [0] * (n * k)
        pivmons = [self.mons[j] for j in self.pivots]
        for r, P in enumerate(polys):
            base = r * k
            for t, m in enumerate(pivmons):
                c = P.get(m)
                if c: flatp[base + t] = c % p
        Vp = nmod_mat(n, k, flatp, p)
        return V - Vp * self.Jr

    def standard_coeffs(self, NF):
        """Rows of NF restricted to standard columns -> list of lists of ints."""
        n = NF.nrows()
        out = []
        for r in range(n):
            out.append([int(NF[r, j]) for j in self.standard])
        return out

def e1_solve(F, A, K, p, verbose=False):
    """Particular solution N (mod p) of sum N_i F_i = A*K with deg N_i = deg K + 1 and pure parities; None if inconsistent."""
    Fp = reduce_mod(F, p); Ap = reduce_mod(A, p); Kp = reduce_mod(K, p)
    d = core.degree(K)
    parA = core.parity_of(Ap); parK = core.parity_of(Kp)
    Fi = [pderiv(Fp, i, p) for i in range(NV)]
    cols = []; polys = []
    for i in range(NV):
        pi = 0 if i < 3 else 1
        par = (parA + parK - pi) % 2
        for m in monomials(d + 1, par):
            cols.append((i, m)); polys.append(shift(Fi[i], m, 1, p))
    rhs = mul(Ap, Kp, p)
    rows = {}
    for P in polys:
        for e in P: rows.setdefault(e, len(rows))
    for e in rhs: rows.setdefault(e, len(rows))
    R, C = len(rows), len(cols)
    flat = [0] * (R * (C + 1))
    for j, P in enumerate(polys):
        for e, v in P.items(): flat[rows[e] * (C + 1) + j] = v
    for e, v in rhs.items(): flat[rows[e] * (C + 1) + C] = v
    M = nmod_mat(R, C + 1, flat, p)
    Mr, rank = M.rref()
    piv = []
    for i in range(rank):
        j = 0
        while j <= C and int(Mr[i, j]) == 0: j += 1
        piv.append(j)
    if piv and piv[-1] == C:
        return None, {'rank': rank, 'rows': R, 'cols': C, 'consistent': False}
    N = [{} for _ in range(NV)]
    for i, j in enumerate(piv):
        v = int(Mr[i, C])
        if v:
            ii, m = cols[j]; N[ii][m] = v
    return N, {'rank': rank, 'rows': R, 'cols': C, 'consistent': True}

def B(D, N, p):
    """B_D(N) = D div N - N(D)  (mod p)."""
    divN = {}
    for i in range(NV):
        divN = add(divN, pderiv(N[i], i, p), p)
    out = mul(D, divN, p)
    for i in range(NV):
        out = sub(out, mul(N[i], pderiv(D, i, p), p), p)
    return out
