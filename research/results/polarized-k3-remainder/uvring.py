"""Jacobian ring of the polarized split F = (f(U)+f(V))/2 as R_f (x) R_f, mod p.
Ternary polys: dict {(a,b,c): int mod p}.  Six-variable polys use core's 6-tuples with
slots 0..2 = U, 3..5 = V (when in the (U,V) frame)."""
import sys, itertools, json, time
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core
from flint import nmod_mat
P = 2147483647

def tern_monomials(k):
    return [(a, b, k - a - b) for a in range(k + 1) for b in range(k + 1 - a)]

def tmul(Pp, Qq, p=P):
    out = {}
    for e1, c1 in Pp.items():
        for e2, c2 in Qq.items():
            e = (e1[0]+e2[0], e1[1]+e2[1], e1[2]+e2[2])
            out[e] = (out.get(e, 0) + c1 * c2) % p
    return {e: v for e, v in out.items() if v}

def tderiv(Pp, i, p=P):
    out = {}
    for e, c in Pp.items():
        if e[i]:
            e2 = list(e); e2[i] -= 1
            out[tuple(e2)] = (c * e[i]) % p
    return {e: v for e, v in out.items() if v}

class TernaryJacobianRing:
    def __init__(self, f, p=P, maxdeg=18):
        self.p = p
        self.f = f
        self.partials = [tderiv(f, i, p) for i in range(3)]
        self.data = {}
        for k in range(maxdeg + 1):
            self._build(k)

    def _build(self, k):
        p = self.p
        mons = tern_monomials(k)
        idx = {m: i for i, m in enumerate(mons)}
        gens = []
        if k >= 5:
            for m in tern_monomials(k - 5):
                for Pi in self.partials:
                    g = tmul({m: 1}, Pi, p)
                    gens.append(g)
        if gens:
            flat = [0] * (len(gens) * len(mons))
            for r, g in enumerate(gens):
                for e, v in g.items():
                    flat[r * len(mons) + idx[e]] = v
            M = nmod_mat(len(gens), len(mons), flat, p)
            Mr, rank = M.rref()
            pivots = []
            rows = []
            for i in range(rank):
                j = 0
                while int(Mr[i, j]) == 0:
                    j += 1
                pivots.append(j)
                rows.append([int(Mr[i, jj]) for jj in range(len(mons))])
        else:
            rank, pivots, rows = 0, [], []
        std = [m for j, m in enumerate(mons) if j not in set(pivots)]
        self.data[k] = dict(mons=mons, idx=idx, pivots=pivots, rows=rows, std=std,
                            stdidx={m: i for i, m in enumerate(std)}, dim=len(std))

    def dim(self, k):
        return self.data[k]['dim'] if k in self.data else 0

    def nf(self, Pp, k):
        """Normal form of a homogeneous degree-k ternary poly, as a vector over std monomials."""
        d = self.data[k]
        p = self.p
        vec = [0] * len(d['mons'])
        for e, v in Pp.items():
            vec[d['idx'][e]] = (vec[d['idx'][e]] + v) % p
        for row, pc in zip(d['rows'], d['pivots']):
            c = vec[pc]
            if c:
                for j, rv in enumerate(row):
                    if rv:
                        vec[j] = (vec[j] - c * rv) % p
        return [vec[d['idx'][m]] for m in d['std']]

def load_f(p=P):
    uv = json.load(open('/tmp/hodge/results/polarized-k3-remainder/uv_forms.json'))
    f6 = core.parse(uv['f'], p)
    return {(e[0], e[1], e[2]): v for e, v in f6.items()}

class TensorRing:
    """R_F = R_f (x) R_f in the (U,V) frame."""
    def __init__(self, R, maxdeg=24):
        self.R = R
        self.p = R.p
        self.bases = {}
        for d in range(maxdeg + 1):
            basis = []
            for a in range(d + 1):
                b = d - a
                if R.dim(a) == 0 or R.dim(b) == 0:
                    continue
                for i in range(R.dim(a)):
                    for j in range(R.dim(b)):
                        basis.append((a, b, i, j))
            self.bases[d] = basis
            # index map
        self.index = {d: {t: i for i, t in enumerate(self.bases[d])} for d in self.bases}

    def dim(self, d):
        return len(self.bases.get(d, []))

    def vector_of(self, Q6, d):
        """Class of a 6-variable poly (U slots 0..2, V slots 3..5) homogeneous of degree d."""
        R, p = self.R, self.p
        vec = [0] * self.dim(d)
        # group monomials by (U-degree)
        by = {}
        for e, v in Q6.items():
            a = e[0] + e[1] + e[2]
            by.setdefault(a, {})[e] = v
        for a, part in by.items():
            b = d - a
            if R.dim(a) == 0 or R.dim(b) == 0:
                continue
            # expand: sum over monomials m_U * m_V: nf(m_U) (x) nf(m_V)
            for e, v in part.items():
                nu = R.nf({(e[0], e[1], e[2]): 1}, a)
                nv = R.nf({(e[3], e[4], e[5]): 1}, b)
                for i, x in enumerate(nu):
                    if not x: continue
                    for j, y in enumerate(nv):
                        if not y: continue
                        k = self.index[d][(a, b, i, j)]
                        vec[k] = (vec[k] + v * x * y) % p
        return vec

    def mult_matrix(self, Q6, dq, d):
        """Matrix (dim R_{d+dq} x dim R_d) of multiplication by Q6 (degree dq)."""
        R, p = self.R, self.p
        rows_n, cols_n = self.dim(d + dq), self.dim(d)
        flat = [0] * (rows_n * cols_n)
        # precompute per term of Q6: U-monomial and V-monomial
        terms = []
        for e, v in Q6.items():
            terms.append((v, (e[0], e[1], e[2]), e[0]+e[1]+e[2], (e[3], e[4], e[5]), e[3]+e[4]+e[5]))
        cache = {}
        def nfm(mon, k):
            key = (mon, k)
            if key not in cache:
                cache[key] = R.nf({mon: 1}, k)
            return cache[key]
        for col, (a, b, i, j) in enumerate(self.bases[d]):
            mu = R.data[a]['std'][i]
            mv = R.data[b]['std'][j]
            for v, eu, du, ev, dv in terms:
                a2, b2 = a + du, b + dv
                if R.dim(a2) == 0 or R.dim(b2) == 0:
                    continue
                pu = (mu[0]+eu[0], mu[1]+eu[1], mu[2]+eu[2])
                pv = (mv[0]+ev[0], mv[1]+ev[1], mv[2]+ev[2])
                nu = nfm(pu, a2); nv = nfm(pv, b2)
                for ii, x in enumerate(nu):
                    if not x: continue
                    for jj, y in enumerate(nv):
                        if not y: continue
                        r = self.index[d + dq][(a2, b2, ii, jj)]
                        flat[r * cols_n + col] = (flat[r * cols_n + col] + v * x * y) % p
        return nmod_mat(rows_n, cols_n, flat, p)

def swap_matrix(T, d):
    """Matrix of the block swap U<->V on R_d in the tensor basis (it permutes basis elements)."""
    n = T.dim(d)
    flat = [0] * (n * n)
    for col, (a, b, i, j) in enumerate(T.bases[d]):
        r = T.index[d][(b, a, j, i)]
        flat[r * n + col] = 1
    return nmod_mat(n, n, flat, T.p)
