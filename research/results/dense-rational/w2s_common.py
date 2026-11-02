"""w2s_common.py -- numerical toolkit for the wave-2 dense-rational joint search (lane-local; the engine is not modified).

Polynomials of a fixed degree k are dense float64 vectors over mons(k) = core.monomials(k) (all Y-parities, fixed order).
Products / adjoints use precomputed index tables.  Everything is numerical (float64); exactness is restored only by the
rational reconstruction + exact re-solve step (hodge.fixed_d / core.verify_witness)."""
import sys, os, time, json, math
sys.path.insert(0, '/home/kejace/devel/sextics-work')
import numpy as np
import scipy.linalg as sla
import scipy.sparse as sp
from fractions import Fraction
from hodge import core, fixed_d, jacobian
from hodge.core import NV, monomials, pderiv, mul, add, sub, scale, shift, power

IID = 'dense-rational'
LANE = '/tmp/hodge/results/dense-rational'
LOG = os.path.join(LANE, 'wave2-search.log.jsonl')
CACHE = os.path.join(LANE, 'w2s_cache')
os.makedirs(CACHE, exist_ok=True)
PVEC = [0, 0, 0, 1, 1, 1]
BASE = 32


def log_event(**rec):
    rec.setdefault('ts', time.strftime('%Y-%m-%dT%H:%M:%S'))
    rec.setdefault('engine', 2)
    rec.setdefault('host', os.uname().nodename)
    rec.setdefault('lane', IID)
    with open(LOG, 'a') as f:
        f.write(json.dumps(rec, default=float) + '\n')
    return rec


def load():
    F, A, q, rec = core.load_instance(IID)
    return F, A, q


class Mons:
    _cache = {}

    @classmethod
    def get(cls, k):
        if k not in cls._cache:
            ms = monomials(k) if k >= 0 else []
            E = np.array(ms, dtype=np.int64).reshape(len(ms), NV)
            keys = (E * (BASE ** np.arange(NV))).sum(1)
            order = np.argsort(keys)
            par = (E[:, 3:].sum(1) % 2).astype(np.int8)
            cls._cache[k] = {'mons': ms, 'E': E, 'keys': keys, 'order': order, 'par': par,
                             'idx': {m: i for i, m in enumerate(ms)}, 'n': len(ms)}
        return cls._cache[k]


def mons(k): return Mons.get(k)['mons']
def nmons(k): return Mons.get(k)['n']
def par_of(k): return Mons.get(k)['par']
def E_of(k): return Mons.get(k)['E']
def sel(k, par): return np.nonzero(Mons.get(k)['par'] == par)[0]


def lookup(k, keys):
    M = Mons.get(k)
    pos = np.searchsorted(M['keys'], keys, sorter=M['order'])
    pos = np.minimum(pos, M['n'] - 1)
    out = M['order'][pos]
    if not np.all(M['keys'][out] == keys):
        raise KeyError('monomial not found in degree %d' % k)
    return out


_T = {}


def mult_table(a, b):
    """int32 (n_a x n_b): index in mons(a+b) of e_a + e_b."""
    if (a, b) not in _T:
        Ea, Eb = E_of(a), E_of(b)
        ka = (Ea * (BASE ** np.arange(NV))).sum(1)
        kb = (Eb * (BASE ** np.arange(NV))).sum(1)
        _T[(a, b)] = lookup(a + b, (ka[:, None] + kb[None, :]).ravel()).reshape(len(ka), len(kb)).astype(np.int32)
    return _T[(a, b)]


def vec(P, k):
    v = np.zeros(nmons(k))
    idx = Mons.get(k)['idx']
    for e, c in P.items():
        v[idx[e]] += float(c)
    return v


def unvec(v, k, tol=0.0):
    ms = mons(k)
    return {ms[i]: float(v[i]) for i in np.nonzero(np.abs(v) > tol)[0]}


def polymul(u, a, v, b):
    iu = np.nonzero(u)[0]; iv = np.nonzero(v)[0]
    n = nmons(a + b)
    if len(iu) == 0 or len(iv) == 0:
        return np.zeros(n)
    T = mult_table(a, b)[np.ix_(iu, iv)]
    return np.bincount(T.ravel(), weights=np.outer(u[iu], v[iv]).ravel(), minlength=n)


def mult_adjoint(u, a, r, c):
    """y in degree c-a with y[e] = sum_{e'} u[e'] r[e'+e]  (adjoint of multiplication by u)."""
    b = c - a
    iu = np.nonzero(u)[0]
    if len(iu) == 0:
        return np.zeros(nmons(b))
    T = mult_table(a, b)[iu, :]
    return (u[iu, None] * r[T]).sum(0)


_D = {}


def deriv_matrix(k, i):
    """sparse (n_{k-1} x n_k): d/dz_i."""
    if (k, i) not in _D:
        E = E_of(k)
        nz = np.nonzero(E[:, i])[0]
        E2 = E[nz].copy(); E2[:, i] -= 1
        rows = lookup(k - 1, (E2 * (BASE ** np.arange(NV))).sum(1))
        _D[(k, i)] = sp.csr_matrix((E[nz, i].astype(float), (rows, nz)), shape=(nmons(k - 1), nmons(k)))
    return _D[(k, i)]


def gen_matrix(Fi, k, par):
    """Dense generator matrix of J(F)_k^par: rows = mons(k) restricted to parity par, columns = m*F_i.
    Returns (G, cols) with cols = list of (i, monomial index in mons(k-5))."""
    rows = sel(k, par)
    rowmap = -np.ones(nmons(k), dtype=np.int64); rowmap[rows] = np.arange(len(rows))
    cols = []
    for i in range(NV):
        for c in sel(k - 5, (par - PVEC[i]) % 2):
            cols.append((i, c))
    G = np.zeros((len(rows), len(cols)))
    T = mult_table(5, k - 5)
    for j, (i, c) in enumerate(cols):
        f = Fi[i]; nz = np.nonzero(f)[0]
        r = rowmap[T[nz, c]]
        G[r, j] = f[nz]
    return G, cols, rows


def column_normalize(G):
    nrm = np.linalg.norm(G, axis=0); nrm[nrm == 0] = 1.0
    return G / nrm, nrm


def mult_matrix(u, a, b, rows_sel=None, cols_sel=None):
    """Dense matrix of multiplication by u (degree a) from degree b to degree a+b, restricted to row/col index lists."""
    if cols_sel is None:
        cols_sel = np.arange(nmons(b))
    T = mult_table(a, b)
    nz = np.nonzero(u)[0]
    if rows_sel is None:
        rows_sel = np.arange(nmons(a + b))
    rowmap = -np.ones(nmons(a + b), dtype=np.int64); rowmap[rows_sel] = np.arange(len(rows_sel))
    M = np.zeros((len(rows_sel), len(cols_sel)))
    R = rowmap[T[np.ix_(nz, cols_sel)]]
    C = np.broadcast_to(np.arange(len(cols_sel))[None, :], R.shape)
    V = np.broadcast_to(u[nz][:, None], R.shape)
    ok = R >= 0
    M[R[ok], C[ok]] = V[ok]
    return M


def svd_complement(G, tol=1e-9, verbose=True, label=''):
    """Orthonormal basis of the orthogonal complement of col(G) (left null space), numerical rank, singular values."""
    t0 = time.time()
    Gn, _ = column_normalize(G)
    U, s, Vt = sla.svd(Gn, full_matrices=True, lapack_driver='gesdd', check_finite=False, overwrite_a=True)
    r = int(np.sum(s > tol * s[0]))
    gap = (s[r - 1] / s[0], (s[r] / s[0] if r < len(s) else 0.0))
    if verbose:
        print(f'[svd {label}] shape {G.shape} rank {r} (sv[r-1]/s0={gap[0]:.2e}, sv[r]/s0={gap[1]:.2e}) {time.time()-t0:.1f}s',
              file=sys.stderr, flush=True)
    return U[:, r:], r, s, gap


def nullspace(Mx, tol=1e-9):
    """Orthonormal basis of the right null space of a dense matrix (columns), with rank and gap info."""
    U, s, Vt = sla.svd(Mx, full_matrices=True, check_finite=False)
    r = int(np.sum(s > tol * s[0])) if len(s) else 0
    gap = (s[r - 1] / s[0] if r else 0.0, (s[r] / s[0] if r < len(s) else 0.0))
    return Vt[r:].T, r, gap


# --------------------------------------------------------------------------------------------------
# cached heavy objects (J complements, K_full bases, lifts)
# --------------------------------------------------------------------------------------------------

def cache_path(name):
    return os.path.join(CACHE, name + '.npz')


def J_complement(Fi, k, par, verbose=True):
    """U_perp (n_rows x r): orthonormal basis of the complement of J_k^par inside S_k^par; cached."""
    p = cache_path(f'Jcomp_k{k}_par{par}')
    if os.path.exists(p):
        z = np.load(p)
        return z['U'], int(z['rank']), z['rows']
    G, cols, rows = gen_matrix(Fi, k, par)
    U, r, s, gap = svd_complement(G, label=f'J_{k}^{par}', verbose=verbose)
    np.savez(p, U=U, rank=r, rows=rows, s_head=s[:5], s_tail=s[max(0, r - 3):r + 3], gap=np.array(gap))
    log_event(event='J-complement', k=k, par=par, rows=int(G.shape[0]), gens=int(G.shape[1]), rank=r,
              complement=int(U.shape[1]), gap_last_nonzero=float(gap[0]), gap_first_zero=float(gap[1]))
    return U, r, rows


def K_full(F, A, Fi, d, pD, verbose=True):
    """Orthonormal basis (n_d x dimK) of K_full(d)^pD = {D in S_d^pD : A D in J_{d+6}}, embedded in mons(d) (zeros off-parity).
    Also returns the numerical rank of [A]: S_d^pD -> S_{d+6}/J (= dim of the image in the Jacobian ring)."""
    p = cache_path(f'Kfull_d{d}_par{pD}')
    if os.path.exists(p):
        z = np.load(p)
        return z['K'], int(z['rankA']), z['cols']
    parT = (pD + 1) % 2   # A is Y-odd
    U, rJ, rows = J_complement(Fi, d + 6, parT, verbose=verbose)
    Av = vec(A, 6)
    cols = sel(d, pD)
    MA = mult_matrix(Av, 6, d, rows_sel=rows, cols_sel=cols)        # (n_rows x n_cols)
    Mred = U.T @ MA                                                  # (r x n_cols)
    Nb, rankA, gap = nullspace(Mred)
    K = np.zeros((nmons(d), Nb.shape[1]))
    K[cols, :] = Nb
    np.savez(p, K=K, rankA=rankA, cols=cols, gap=np.array(gap))
    log_event(event='K_full', d=d, pD=pD, dimS=int(len(cols)), rankA=int(rankA), dimKfull=int(Nb.shape[1]),
              dimR_target=int(U.shape[1]), gap_last_nonzero=float(gap[0]), gap_first_zero=float(gap[1]))
    if verbose:
        print(f'[K_full d={d} par={pD}] dim S={len(cols)} rank[A]={rankA} (target R dim {U.shape[1]}) dim K_full={Nb.shape[1]} gap {gap}',
              file=sys.stderr, flush=True)
    return K, rankA, cols


def J_basis(Fi, d, pD):
    """Orthonormal basis of J_d^pD inside mons(d) (zeros off-parity) and its dimension."""
    G, cols, rows = gen_matrix(Fi, d, pD)
    Gn, _ = column_normalize(G)
    U, s, Vt = sla.svd(Gn, full_matrices=False, check_finite=False)
    r = int(np.sum(s > 1e-9 * s[0]))
    B = np.zeros((nmons(d), r)); B[rows, :] = U[:, :r]
    return B, r


def e1_lifts(F, A, Fi, d, pD, Dmat, verbose=True):
    """Minimum-norm numerical solutions N0 of N0.grad F = A D for each row D of Dmat (K x n_d).
    Returns N0 (K x 6 x n_{d+1}) and the relative E1 residuals."""
    parT = (pD + 1) % 2
    G, cols, rows = gen_matrix(Fi, d + 6, parT)
    Gn, nrm = column_normalize(G)
    Av = vec(A, 6)
    rhs = np.stack([polymul(Av, 6, Dmat[k], d)[rows] for k in range(Dmat.shape[0])], axis=1)   # n_rows x K
    X, res, rk, s = sla.lstsq(Gn, rhs, cond=1e-10, lapack_driver='gelsd', check_finite=False)
    X = X / nrm[:, None]
    resid = np.linalg.norm(G @ X - rhs, axis=0) / np.maximum(np.linalg.norm(rhs, axis=0), 1e-300)
    N0 = np.zeros((Dmat.shape[0], NV, nmons(d + 1)))
    for j, (i, c) in enumerate(cols):
        N0[:, i, c] = X[j, :]
    if verbose:
        print(f'[e1_lifts d={d} par={pD}] K={Dmat.shape[0]} max rel E1 residual {resid.max():.2e}', file=sys.stderr, flush=True)
    return N0, resid


def jac_family(F, A, Fi, d, pD):
    """Jacobian-type family at degree d: w-basis (i, m) with m in S_{d-5}^{(pD+p_i)%2}; D_k = m F_i, N0_k = A m e_i (exact).
    Returns Dmat (K x n_d), N0mat (K x 6 x n_{d+1}), labels, gauge vectors (rows, in w-coordinates)."""
    Av = vec(A, 6)
    labels = []
    for i in range(NV):
        for c in sel(d - 5, (pD + PVEC[i]) % 2):
            labels.append((i, c))
    K = len(labels)
    Dmat = np.zeros((K, nmons(d)))
    N0 = np.zeros((K, NV, nmons(d + 1)))
    T5 = mult_table(5, d - 5); T6 = mult_table(6, d - 5)
    for k, (i, c) in enumerate(labels):
        f = Fi[i]; nz = np.nonzero(f)[0]
        Dmat[k, T5[nz, c]] = f[nz]
        nzA = np.nonzero(Av)[0]
        N0[k, i, T6[nzA, c]] = Av[nzA]
    # gauge: w = l * Euler, l in S_{d-6}^{pD}: coefficient 1 on (i, m = l z_i)
    gauge = []
    idx = {lab: k for k, lab in enumerate(labels)}
    Ed6 = E_of(d - 6) if d >= 6 else None
    if d >= 6:
        for c in sel(d - 6, pD):
            v = np.zeros(K)
            for i in range(NV):
                e = Ed6[c].copy(); e[i] += 1
                cc = lookup(d - 5, np.array([(e * (BASE ** np.arange(NV))).sum()]))[0]
                v[idx[(i, cc)]] = 1.0
            gauge.append(v)
    return Dmat, N0, labels, np.array(gauge)


def complement_projector(K, vectors):
    """Orthonormal basis (K x dimP) of the orthogonal complement of the span of the given vectors (rows)."""
    if vectors is None or len(vectors) == 0:
        return np.eye(K)
    Q, R = np.linalg.qr(np.asarray(vectors, float).T)   # K x g
    r = int(np.sum(np.abs(np.diag(R)) > 1e-10))
    U, s, Vt = np.linalg.svd(np.asarray(vectors, float).T, full_matrices=True)
    return U[:, r:]


# --------------------------------------------------------------------------------------------------
# exact side: rational reconstruction + exact re-solve
# --------------------------------------------------------------------------------------------------

def rational_candidates(Dv, d, dens=(1, 2, 6, 12, 60, 120, 1000, 10000, 100000, 1000000), tol_rel=1e-7, max_terms=None):
    """Rational approximations of the numerical D (normalized by its largest coefficient)."""
    Dv = np.asarray(Dv, float)
    j = int(np.argmax(np.abs(Dv)))
    Dv = Dv / Dv[j]
    out = []
    nz = np.nonzero(np.abs(Dv) > 1e-12)[0]
    seen = set()
    for den in dens:
        fr = [Fraction(float(Dv[i])).limit_denominator(den) for i in nz]
        err = max(abs(float(f) - Dv[i]) for f, i in zip(fr, nz))
        key = tuple(fr)
        if key in seen:
            continue
        seen.add(key)
        out.append({'den': den, 'err': err, 'poly': {mons(d)[i]: f for f, i in zip(fr, nz) if f != 0}})
    return out


def exact_confirm(F, A, q, Dpoly, verbose=True):
    """E1 (mod p) then exact fixed-D solve; returns (status, data).  status in {'e1-fails','inconsistent','solved','lift-failed'}."""
    if not Dpoly:
        return 'inconsistent', {'note': 'D = 0'}
    ok, info = fixed_d.e1_solvable(F, A, q, Dpoly)
    if not ok:
        return 'e1-fails', {'e1': info}
    try:
        out = fixed_d.solve_fixed_D_fast(F, A, q, Dpoly, verbose=verbose)
        if not out['consistent']:
            return 'inconsistent', {'rank': out['rank'], 'nullity': out['nullity']}
        if out.get('exact'):
            return 'solved', {'N': out['N'], 'Q': out['Q'], 'rank': out['rank'], 'nullity': out['nullity']}
    except Exception as ex:
        if verbose:
            print(f'[exact_confirm] fast solver failed ({ex}); dense fallback', file=sys.stderr)
    status, data = fixed_d.solve_fixed_D(F, A, q, Dpoly, verbose=verbose)
    return status, data


def save_witness_candidate(F, A, q, D, N, Q, path=os.path.join(LANE, 'witness-candidate.json')):
    ok, info = core.verify_witness(F, A, q, D, N, Q)
    divF, _ = jacobian.divisible_by(F, D)
    rec = {'instance': IID, 'd': core.degree(D), 'D': core.to_str(D), 'N': [core.to_str(n) for n in N],
           'Q': core.to_str(Q), 'verified': bool(ok), 'F_divides_D': bool(divF), 'engine': 2}
    with open(path, 'w') as f:
        json.dump(rec, f, indent=1)
    return rec
