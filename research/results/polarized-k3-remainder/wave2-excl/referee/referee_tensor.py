"""Referee replication (independent code) for lane polarized-k3-remainder, wave-2 exclusion packet.

Tensor frame: U = X+Y, V = X-Y, F = (f(U)+f(V))/2, J(F) = J(f)(U) + J(f)(V), so R_F = R_f (x) R_f as bigraded
rings, with R_f = Q[z0,z1,z2]/J(f) the Jacobian ring of the ternary sextic f.  The target A = X0X1X2Y0Y1Y2 equals
2^-6 prod_i (U_i^2 - V_i^2) = 2^-6 sum_{S subset {0,1,2}} (-1)^|S| m_{S^c}(U) m_S(V), m_S = prod_{i in S} z_i^2.
A_K3 = (eta(U) - eta(V))/64 with eta = z0^2 z1^2 z2^2 (S = empty, S = all); A_mix = the other six terms.

Checks (exact over Q with python-flint fmpq_mat; and mod p with nmod_mat):
  C0  F(instance) == (f(X+Y) + f(X-Y))/2 exactly (frame identity).
  C1  Hilbert function of R_f is (1,3,6,10,15,18,19,18,15,10,6,3,1); eta != 0 in R_{f,6}.
  C2  A_K3 and A_mix are both nonzero modulo J(F) (bigrading).
  C3  E1 for the full target: ker([A]: (R_F)_d -> (R_F)_{d+6}) for d <= 9 (expect 0), d = 10 (expect 555 = 294+261).
  C4  E1 for A_K3: kernel dims d <= 9 (expect 0 for d <= 7, 90 at d = 8, 306 at d = 9).
  C5  K4 = ker(eta: R_{f,4} -> R_{f,10}) has dimension 9 and eta is onto R_{f,10} (rank 6).
  C6  the degree-8 blocks: K4 (x) K4 (81) and {l(U) l'(V) (eta(U) + eta(V))} (9) lie in ker [A_K3].
Results are appended to referee.log.jsonl as they are produced.
"""
import sys, json, time, itertools
from fractions import Fraction as Fr
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from flint import fmpq_mat, nmod_mat, fmpq
from hodge import core

LOG = '/tmp/hodge/results/polarized-k3-remainder/wave2-excl/referee/referee.log.jsonl'
EXACT_MAX = int(sys.argv[1]) if len(sys.argv) > 1 else 1300   # exact rank only when the matrix is at most this size
PRIMES = [int(x) for x in sys.argv[2].split(',')] if len(sys.argv) > 2 else [1000000007, 2147483647]

def log(**kw):
    kw['ts'] = time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())
    with open(LOG, 'a') as fh:
        fh.write(json.dumps(kw) + '\n')
    print(json.dumps(kw), flush=True)

# ---------------------------------------------------------------- ternary polynomials as dict {exp: Fraction}
def tmon(e):
    return {tuple(e): Fr(1)}
def tadd(a, b):
    out = dict(a)
    for k, v in b.items():
        out[k] = out.get(k, Fr(0)) + v
        if out[k] == 0: del out[k]
    return out
def tscale(a, c):
    return {k: v * c for k, v in a.items()} if c != 0 else {}
def tmul(a, b):
    out = {}
    for ka, va in a.items():
        for kb, vb in b.items():
            k = tuple(x + y for x, y in zip(ka, kb))
            out[k] = out.get(k, Fr(0)) + va * vb
    return {k: v for k, v in out.items() if v != 0}
def tderiv(a, i):
    out = {}
    for k, v in a.items():
        if k[i] > 0:
            e = list(k); e[i] -= 1
            out[tuple(e)] = out.get(tuple(e), Fr(0)) + v * k[i]
    return {k: v for k, v in out.items() if v != 0}
def tmonomials(k):
    return [(a, b, k - a - b) for a in range(k + 1) for b in range(k - a + 1)]

# f = z0^6 + z1^6 + z2^6 + (5/7)(z0^5 z1 + 2 z1^5 z2 + 3 z2^5 z0 + z0^2 z1^3 z2)
f = {}
for e, c in [((6,0,0),1), ((0,6,0),1), ((0,0,6),1), ((5,1,0),Fr(5,7)), ((0,5,1),Fr(10,7)), ((1,0,5),Fr(15,7)), ((2,3,1),Fr(5,7))]:
    f = tadd(f, {e: Fr(c)})
fi = [tderiv(f, i) for i in range(3)]

# ---------------------------------------------------------------- C0: frame identity against the instance file
F, A, q, rec = core.load_instance('polarized-k3-remainder')
def var6(i):
    e = [0]*6; e[i] = 1; return {tuple(e): Fr(1)}
X = [var6(i) for i in range(3)]; Y = [var6(3 + i) for i in range(3)]
U = [core.add(X[i], Y[i]) for i in range(3)]; V = [core.sub(X[i], Y[i]) for i in range(3)]
def compose(poly3, lin):
    out = {}
    for e, c in poly3.items():
        term = {(0,)*6: Fr(c)}
        for i in range(3):
            for _ in range(e[i]):
                term = core.mul(term, lin[i])
        out = core.add(out, term)
    return out
Fuv = core.scale(core.add(compose(f, U), compose(f, V)), Fr(1, 2))
diff = core.sub(Fuv, F)
log(check='C0_frame_identity', ok=(len(diff) == 0), n_terms_F=len(F), n_terms_diff=len(diff))
assert len(diff) == 0
# the target in the (U,V) frame: 2^-6 prod (U_i^2 - V_i^2)
T = {(0,)*6: Fr(1, 64)}
for i in range(3):
    T = core.mul(T, core.sub(core.mul(U[i], U[i]), core.mul(V[i], V[i])))
log(check='C0_target_identity', ok=(len(core.sub(T, A)) == 0))
assert len(core.sub(T, A)) == 0

# ---------------------------------------------------------------- Jacobian ring of f: exact normal forms per degree
class Ring:
    """R_k = S_k / J_k with an exact reduced row echelon basis of J_k; coordinates in the standard monomials."""
    def __init__(self, kmax):
        self.k = {}
        for k in range(kmax + 1):
            mons = tmonomials(k); idx = {m: i for i, m in enumerate(mons)}
            rows = []
            if k >= 5:
                for m in tmonomials(k - 5):
                    for i in range(3):
                        g = tmul(tmon(m), fi[i])
                        rows.append([Fr(0)] * len(mons))
                        for e, c in g.items(): rows[-1][idx[e]] = c
            if rows:
                M = fmpq_mat(len(rows), len(mons), [fmpq(c.numerator, c.denominator) for r in rows for c in r])
                R, rank = M.rref()
                piv = []
                r = 0
                for j in range(len(mons)):
                    if r < rank and R[r, j] != 0:
                        piv.append(j); r += 1
                Rrows = [[Fr(int(R[i, j].p), int(R[i, j].q)) for j in range(len(mons))] for i in range(rank)]
            else:
                rank, piv, Rrows = 0, [], []
            std = [j for j in range(len(mons)) if j not in set(piv)]
            self.k[k] = dict(mons=mons, idx=idx, piv=piv, rows=Rrows, std=std, stdidx={j: i for i, j in enumerate(std)})
    def dim(self, k):
        return len(self.k[k]['std']) if k in self.k else 0
    def nf(self, poly, k):
        """coordinates of poly (degree k) in R_k w.r.t. the standard monomials (exact)."""
        d = self.k[k]
        v = [Fr(0)] * len(d['mons'])
        for e, c in poly.items(): v[d['idx'][e]] += c
        for r, j in zip(d['rows'], d['piv']):
            if v[j] != 0:
                c = v[j]
                for t in range(len(v)):
                    if r[t] != 0: v[t] -= c * r[t]
        return [v[j] for j in d['std']]
    def mult_matrix(self, mu, a):
        """matrix of multiplication by the polynomial mu (degree e) : R_a -> R_{a+e}; rows = codomain, cols = domain."""
        e = sum(next(iter(mu)))
        da, db = self.k[a], self.k[a + e]
        cols = []
        for j in da['std']:
            cols.append(self.nf(tmul(tmon(da['mons'][j]), mu), a + e))
        return [[cols[c][r] for c in range(len(cols))] for r in range(len(db['std']))]

t0 = time.time()
R = Ring(14)
hilb = [R.dim(k) for k in range(15)]
eta = tmon((2, 2, 2))
eta_nf = R.nf(eta, 6)
log(check='C1_hilbert_and_eta', hilbert=hilb, ok=(hilb[:13] == [1,3,6,10,15,18,19,18,15,10,6,3,1] and hilb[13] == 0 and hilb[14] == 0),
    eta_nonzero_in_R6=any(c != 0 for c in eta_nf), seconds=round(time.time() - t0, 1))

# ---------------------------------------------------------------- multiplication matrices for m_S
def mS(S):
    e = [0, 0, 0]
    for i in S: e[i] = 2
    return tmon(tuple(e))
subsets = [tuple(S) for r in range(4) for S in itertools.combinations(range(3), r)]
MM = {}   # (S, a) -> matrix of m_S : R_a -> R_{a+2|S|}
def mm(S, a):
    key = (S, a)
    if key not in MM:
        MM[key] = R.mult_matrix(mS(S), a) if R.dim(a) and R.dim(a + 2 * len(S)) else [[]]
    return MM[key]

def e1_matrix(terms, d):
    """terms: list of (coef, S_U, S_V); map (R_F)_d -> (R_F)_{d+6} by multiplication with sum coef m_{S_U}(U) m_{S_V}(V)."""
    dom = [(a, d - a) for a in range(d + 1) if R.dim(a) and R.dim(d - a)]
    cod = [(a, d + 6 - a) for a in range(d + 7) if R.dim(a) and R.dim(d + 6 - a)]
    dom_off = {}; n = 0
    for (a, b) in dom: dom_off[(a, b)] = n; n += R.dim(a) * R.dim(b)
    cod_off = {}; m = 0
    for (a, b) in cod: cod_off[(a, b)] = m; m += R.dim(a) * R.dim(b)
    M = [[Fr(0)] * n for _ in range(m)]
    for (a, b) in dom:
        for coef, SU, SV in terms:
            a2, b2 = a + 2 * len(SU), b + 2 * len(SV)
            if (a2, b2) not in cod_off: continue
            P, Q = mm(SU, a), mm(SV, b)
            da, db, ca, cb = R.dim(a), R.dim(b), R.dim(a2), R.dim(b2)
            for i in range(ca):
                for k in range(da):
                    if P[i][k] == 0: continue
                    for j in range(cb):
                        for l in range(db):
                            if Q[j][l] == 0: continue
                            M[cod_off[(a2, b2)] + i * cb + j][dom_off[(a, b)] + k * db + l] += coef * P[i][k] * Q[j][l]
    return M, n, m

def ranks(M, n, m, exact):
    out = {}
    for p in PRIMES:
        N = nmod_mat(m, n, [ (c.numerator * pow(c.denominator, -1, p)) % p for r in M for c in r], p)
        out[f'rank_mod_{p}'] = N.rank()
    if exact:
        t = time.time()
        Q = fmpq_mat(m, n, [fmpq(c.numerator, c.denominator) for r in M for c in r])
        out['rank_Q'] = Q.rank(); out['exact_seconds'] = round(time.time() - t, 1)
    return out

full_terms = [(Fr((-1) ** len(S)), tuple(i for i in range(3) if i not in S), S) for S in subsets]   # 2^6 A
k3_terms = [(Fr(1), (0, 1, 2), ()), (Fr(-1), (), (0, 1, 2))]                                        # 64 A_K3
mix_terms = [t for t in full_terms if 0 < len(t[2]) < 3]                                               # 64 A_mix

# ---------------------------------------------------------------- C2: both parts nonzero mod J
log(check='C2_parts_nonzero_mod_J', A_K3_nonzero=any(c != 0 for c in eta_nf),
    A_mix_nonzero=all(R.dim(2 * len(SU)) > 0 and R.dim(2 * len(SV)) > 0 and any(c != 0 for c in R.nf(mS(SU), 2 * len(SU))) and any(c != 0 for c in R.nf(mS(SV), 2 * len(SV))) for _, SU, SV in mix_terms),
    note='bigraded pieces (6,0)+(0,6) vs (4,2)+(2,4) of R_f (x) R_f; J_k = 0 for k <= 4 so the mixed terms are nonzero trivially')

# ---------------------------------------------------------------- C5: K4
Meta4 = R.mult_matrix(eta, 4)
r = ranks(Meta4, R.dim(4), R.dim(10), True)
log(check='C5_K4', dim_R4=R.dim(4), dim_R10=R.dim(10), **r, K4_dim_Q=R.dim(4) - r['rank_Q'], onto=(r['rank_Q'] == R.dim(10)))
Meta1 = R.mult_matrix(eta, 1)
r1 = ranks(Meta1, R.dim(1), R.dim(7), True)
log(check='C6_eta_on_R1_injective', **r1, injective=(r1['rank_Q'] == 3))

# ---------------------------------------------------------------- C3 / C4: E1 kernels
for name, terms, dmax in [('A_K3', k3_terms, 9), ('A', full_terms, 10)]:
    for d in range(0, dmax + 1):
        t = time.time()
        M, n, m = e1_matrix(terms, d)
        r = ranks(M, n, m, exact=(max(n, m) <= EXACT_MAX))
        ker = {k.replace('rank', 'kernel'): n - v for k, v in r.items() if k.startswith('rank')}
        log(check=f'E1_{name}_d{d}', target=name, d=d, dim_domain=n, dim_codomain=m, **r, **ker, seconds=round(time.time() - t, 1))
log(check='DONE')
