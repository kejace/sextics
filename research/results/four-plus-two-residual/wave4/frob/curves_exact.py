"""Exact L-polynomials over K = Q(zeta3) of the chi-sectors of the curves C_B (w^3=B), E' (sigma-quotient), P (Prym), C_f (w^3=f), C_g (w^3=g)
from twisted cubic-character sums over F_{p^k}, k <= kmax.  Convention as in sums_k1.py (chi = cubic residue symbol mod pi, pi primary, zeta3 -> omega_p).
Over F_{p^k} the character is chi(N_{F_{p^k}/F_p}(.)).  F_{p^k} = F_p[x]/(x^k - m), m a non-k-th power (k=2: non-square; k=3: non-cube; k=4: x^4-m irreducible <=> m non-square, checked).
Power sums s_j = -(sum over P^1(F_{p^j}) of chi_j(curve poly)) = Tr(Frob^j | sector); sigma-twist separates E' and P.
Output: jsonl records with L-polynomial coefficients (as [u,v] = u+v*zeta3) and checks.
"""
import numpy as np, sys, json, time, math, itertools
sys.path.insert(0, '/tmp/hodge/results/four-plus-two-residual/wave4/frob')
from sums_k1 import chi_table, nonresidue, primary

def ext_elems(p, k, m):
    """all elements of F_{p^k} as (p^k, k) int64 array of coefficients (basis 1, x, .., x^{k-1}), x^k = m."""
    grids = np.meshgrid(*[np.arange(p, dtype=np.int64)]*k, indexing='ij')
    return np.stack([g.ravel() for g in grids], axis=1)

def mul(a, b, p, k, m):
    """multiply arrays of F_{p^k} elements (coefficient arrays shape (N,k)), x^k = m."""
    N = a.shape[0]; res = np.zeros((N, 2*k-1), dtype=np.int64)
    for i in range(k):
        for j in range(k):
            res[:, i+j] = (res[:, i+j] + a[:, i]*b[:, j]) % p
    out = res[:, :k].copy()
    for d in range(k, 2*k-1):
        out[:, d-k] = (out[:, d-k] + m*res[:, d]) % p
    return out

def norm(a, p, k, m):
    """N_{F_{p^k}/F_p}(a) = a * Frob(a) * ... ; Frob(x) = x^p = m^((p-1)/k) x  (since x^k=m, x^p = x * (x^k)^((p-1)/k))."""
    c = pow(m, (p-1)//k, p)   # x^p = c x
    cur = a.copy(); prod = a.copy()
    for j in range(1, k):
        # apply Frobenius to cur: coefficient of x^i gets multiplied by c^i
        cur = cur * np.array([pow(c, i, p) for i in range(k)], dtype=np.int64) % p
        prod = mul(prod, cur, p, k, m)
    assert (prod[:, 1:] == 0).all()
    return prod[:, 0]

def poly_eval(coeffs_fp, t, p, k, m):
    """evaluate sum c_i t^i (c_i in F_p) at t in F_{p^k} (array (N,k))."""
    N = t.shape[0]; val = np.zeros((N, k), dtype=np.int64); pw = np.zeros((N, k), dtype=np.int64); pw[:, 0] = 1
    for c in coeffs_fp:
        val = (val + c*pw) % p; pw = mul(pw, t, p, k, m)
    return val

def find_m(p, k):
    for m in range(2, p):
        if k == 2 and pow(m, (p-1)//2, p) == p-1: return m
        if k == 3 and pow(m, (p-1)//3, p) != 1: return m
        if k == 4 and pow(m, (p-1)//2, p) == p-1 and pow(m, (p-1)//4, p) not in (1, p-1): return m
    raise ValueError

def twisted_sum(p, k, tab, coeffs, m, n_k=None):
    """sum over P^1(F_{p^k}) of chi(N(poly(t))) with poly = sum coeffs[i] t^i, coeffs in F_p, plus the point at infinity (leading coeff)."""
    t = ext_elems(p, k, m)
    if n_k is None:
        val = poly_eval(coeffs, t, p, k, m)
    else:
        # sigma-twisted: evaluate B(1, c t) with c^2 = n_k (n_k non-square in F_{p^k}): coefficient i gets n_k^(i/2) for even i
        val = np.zeros((t.shape[0], k), dtype=np.int64); pw = np.zeros((t.shape[0], k), dtype=np.int64); pw[:, 0] = 1
        nk = np.zeros((1, k), dtype=np.int64); nk[0] = n_k; nkpow = np.zeros((1, k), dtype=np.int64); nkpow[0, 0] = 1
        for i, c in enumerate(coeffs):
            if i % 2 == 0 and i > 0: nkpow = mul(nkpow, nk, p, k, m)
            if c: val = (val + c*mul(pw, np.repeat(nkpow, t.shape[0], axis=0), p, k, m)) % p
            pw = mul(pw, t, p, k, m)
    nv = norm(val, p, k, m)
    cnt = np.bincount(tab[nv], minlength=4)
    inf_val = coeffs[-1] % p if n_k is None else None
    if n_k is None:
        cnt = cnt + np.bincount(tab[np.array([inf_val])], minlength=4)
    else:
        # point at infinity (0:1) -> (0 : c) ~ (0:1): B(0,1)*c^6 = n_k^3 ; chi(N(n_k^3))
        nk3 = mul(mul(nk, nk, p, k, m), nk, p, k, m)
        cnt = cnt + np.bincount(tab[norm(nk3, p, k, m)], minlength=4)
    return [int(cnt[0]-cnt[2]), int(cnt[1]-cnt[2])]

def nonsquare_ext(p, k, m):
    # find element of F_{p^k} with y^((p^k-1)/2) = -1: use x itself if k odd? x^((p^k-1)/2): just search small elements
    q = p**k
    for coeffs in itertools.product(range(p), repeat=k):
        y = np.array([coeffs], dtype=np.int64)
        if not y.any(): continue
        # compute y^((q-1)/2) by square-and-multiply
        r = np.zeros((1, k), dtype=np.int64); r[0, 0] = 1; b = y.copy(); e = (q-1)//2
        while e:
            if e & 1: r = mul(r, b, p, k, m)
            b = mul(b, b, p, k, m); e >>= 1
        if r[0, 0] == p-1 and not r[0, 1:].any(): return list(coeffs)
    raise ValueError

class KInt:
    """elements u + v zeta3 of Z[zeta3]"""
    def __init__(s, u, v=0): s.u, s.v = int(u), int(v)
    def __add__(s, o): return KInt(s.u+o.u, s.v+o.v)
    def __sub__(s, o): return KInt(s.u-o.u, s.v-o.v)
    def __mul__(s, o):
        if isinstance(o, int): return KInt(s.u*o, s.v*o)
        return KInt(s.u*o.u - s.v*o.v, s.u*o.v + s.v*o.u - s.v*o.v)   # zeta^2 = -1 - zeta
    def conj(s): return KInt(s.u - s.v, -s.v)
    def norm(s): return s.u*s.u - s.u*s.v + s.v*s.v
    def div_int(s, d):
        assert s.u % d == 0 and s.v % d == 0, (s.u, s.v, d); return KInt(s.u//d, s.v//d)
    def lst(s): return [s.u, s.v]
    def __repr__(s): return f'({s.u}+{s.v}z)'

def newton(power_sums, r):
    """elementary symmetric e_1..e_r from power sums s_1..s_r (KInt), L(T) = 1 - e1 T + e2 T^2 - ... """
    e = [KInt(1)]
    for k in range(1, r+1):
        acc = KInt(0)
        for i in range(1, k+1):
            term = e[k-i] * power_sums[i-1]
            acc = acc + (term if i % 2 == 1 else KInt(-term.u, -term.v))
        e.append(acc.div_int(k))
    return e

def one_prime(p, kmax):
    t0 = time.time(); tab, om, a, b = chi_table(p); n = nonresidue(p); i7 = pow(7, -1, p)
    polys = {'f': [1, 5*i7 % p, 0, 0, 0, 10*i7 % p, 1], 'g': [1, 5*i7 % p, 0, 0, 0, 15*i7 % p, 1], 'B': [1, 0, 15*i7 % p, 0, 10*i7 % p, 0, 1]}
    S = {name: [] for name in polys}; SBn = []
    for k in range(1, kmax+1):
        m = find_m(p, k) if k > 1 else 0
        if k == 1:
            from sums_k1 import curve_sum
            for name, co in polys.items(): S[name].append(curve_sum(p, tab, co))
            SBn.append(curve_sum(p, tab, [1, 0, 15*i7*n, 0, 10*i7*n*n, 0, pow(n, 3, p)]))
            S = {name: [[int(c[0]-c[2]), int(c[1]-c[2])] for c in v] for name, v in S.items()}; SBn = [[int(c[0]-c[2]), int(c[1]-c[2])] for c in SBn]
        else:
            for name, co in polys.items(): S[name].append(twisted_sum(p, k, tab, co, m))
            SBn.append(twisted_sum(p, k, tab, polys['B'], m, n_k=nonsquare_ext(p, k, m)))
    # power sums of Frobenius on sectors: s_j = -(sum)
    def neg(x): return KInt(-x[0], -x[1])
    sE = [(neg(S['B'][j]) + neg(SBn[j])).div_int(2) for j in range(kmax)]
    sP = [(neg(S['B'][j]) - neg(SBn[j])).div_int(2) for j in range(kmax)]
    sf = [neg(x) for x in S['f']]; sg = [neg(x) for x in S['g']]
    out = {'p': p, 'kmax': kmax, 'omega': int(om), 'pi': [a, b], 'secs': round(time.time()-t0, 1)}
    out['E_trace'] = sE[0].lst(); out['E_norm_check'] = sE[0].norm() == p
    out['E_power_sums'] = [x.lst() for x in sE]
    # E' rank 1: s_j must equal s_1^j
    pw = KInt(1); ok = True
    for j in range(kmax):
        pw = pw * sE[0]; ok = ok and (pw.lst() == sE[j].lst())
    out['E_rank1_check'] = ok
    if kmax >= 3:
        eP = newton(sP, 3); out['LP_chi'] = [x.lst() for x in eP]   # L = 1 - e1 T + e2 T^2 - e3 T^3
        out['LP_det_norm_check'] = eP[3].norm() == p**3
    out['P_power_sums'] = [x.lst() for x in sP]
    out['f_power_sums'] = [x.lst() for x in sf]; out['g_power_sums'] = [x.lst() for x in sg]
    if kmax >= 4:
        ef = newton(sf, 4); eg = newton(sg, 4)
        out['Lf_chi'] = [x.lst() for x in ef]; out['Lg_chi'] = [x.lst() for x in eg]
        out['f_det_norm_check'] = ef[4].norm() == p**4; out['g_det_norm_check'] = eg[4].norm() == p**4
    return out

if __name__ == '__main__':
    lo, hi, kmax, outf = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
    from sympy import primerange
    with open(outf, 'a') as fh:
        for p in primerange(lo, hi+1):
            if p % 3 != 1 or p == 7: continue
            rec = one_prime(p, kmax); fh.write(json.dumps(rec) + '\n'); fh.flush(); print(p, rec['secs'], rec.get('E_norm_check'), rec.get('E_rank1_check'), rec.get('LP_det_norm_check'), rec.get('f_det_norm_check'), flush=True)
