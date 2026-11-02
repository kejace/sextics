"""Frobenius traces (k=1) on the chi-sectors of T_H: v^3 = H (P(1,1,1,1,2)) and of the curves C_B, E', P, C_f, C_g,
from cubic-character sums over F_p, p = 1 mod 3.  H = f(X0,X1) + g(Y0,Y1) + M,  M = (5/7)X0^3X1Y0Y1 + (10/7)X0X1Y0^3Y1.
S0 = sum_{P^3(F_p)} chi(H), S1 = sum chi(H_n) (s-twist, s=(X0,X1)->(-X0,-X1), n a non-residue):
  Tr(Frob | H^3_c(P^3-S_H, L_chi)^{s=+1}) = -(S0+S1)/2,   Tr(... ^{s=-1}) = -(S0-S1)/2.
Convention: p = pi*pibar, pi = a + b*zeta3 primary (a=2 mod 3, b=0 mod 3, b>0), omega_p = -a/b mod p; chi(x) = omega^k with x^((p-1)/3) = omega^k.
Sums are recorded as [u, v] meaning u + v*zeta3.
"""
import numpy as np, sys, json, time, math
from multiprocessing import Pool

def primary(p):
    for b in range(3, int(2*math.sqrt(p))+2, 3):
        for a in range(-int(2*math.sqrt(p))-1, int(2*math.sqrt(p))+2):
            if a*a - a*b + b*b == p and a % 3 == 2:
                return a, b
    raise ValueError(p)

def chi_table(p):
    a, b = primary(p)
    om = (-a * pow(b, -1, p)) % p
    assert (om*om + om + 1) % p == 0 and om != 1
    x = np.arange(p, dtype=np.int64)
    e = np.ones(p, dtype=np.int64)
    base = x.copy(); n = (p-1)//3
    while n:
        if n & 1: e = (e*base) % p
        base = (base*base) % p; n >>= 1
    tab = np.full(p, -1, dtype=np.int8)
    tab[e == 1] = 0; tab[e == om] = 1; tab[e == (om*om) % p] = 2; tab[0] = 3
    assert (tab >= 0).all()
    return tab, om, a, b

def nonresidue(p):
    for n in range(2, p):
        if pow(n, (p-1)//2, p) == p-1: return n
    raise ValueError

def acc(tab, vals, cnt):
    c = np.bincount(tab[vals % tab.shape[0]], minlength=4)
    cnt += c[:4]

def th_sums(p, tab, nf, nM1, nM2):
    """sum over P^3(F_p) of chi(nf*f + nM1*c1*x0^3x1y0y1 + nM2*c2*x0x1y0^3y1 + g) ; returns counts [n_1, n_om, n_om2, n_zero]."""
    c1 = 5*pow(7,-1,p) % p; c2 = 10*pow(7,-1,p) % p; c3 = 15*pow(7,-1,p) % p
    y0 = np.repeat(np.arange(p, dtype=np.int64), p); y1 = np.tile(np.arange(p, dtype=np.int64), p)
    y02 = y0*y0 % p; y03 = y02*y0 % p; y05 = y03*y02 % p; y12 = y1*y1 % p; y15 = y12*y12 % p * y1 % p
    g = (y05*y0 + c1*y05 % p*y1 + c3*y0 % p*y15 + y15*y1) % p
    m1 = (nM1*c1 % p) * (y0*y1 % p) % p
    m2 = (nM2*c2 % p) * (y03*y1 % p) % p
    cnt = np.zeros(4, dtype=np.int64)
    # chart x0 = 1
    for x1 in range(p):
        fval = (1 + c1*x1 + c2*pow(x1,5,p) + pow(x1,6,p)) % p
        vals = (nf*fval % p + x1*m1 + x1*m2 + g) % p
        acc(tab, vals, cnt)
    # chart x0 = 0, x1 = 1: f(0,1) = 1, M = 0
    acc(tab, (nf % p + g) % p, cnt)
    # x0 = x1 = 0: P^1 in (y0:y1): (1:y1) and (0:1)
    t = np.arange(p, dtype=np.int64); t2 = t*t % p; t5 = t2*t2 % p*t % p
    acc(tab, (1 + c1*t + c3*t5 + t5*t) % p, cnt)
    acc(tab, np.array([1], dtype=np.int64), cnt)
    return cnt

def curve_sum(p, tab, coeffs):
    """sum over P^1(F_p) of chi(sum coeffs[i] t^i) with leading coeff for (0:1)."""
    t = np.arange(p, dtype=np.int64); vals = np.zeros(p, dtype=np.int64); pw = np.ones(p, dtype=np.int64)
    for c in coeffs:
        vals = (vals + c*pw) % p; pw = pw*t % p
    cnt = np.zeros(4, dtype=np.int64); acc(tab, vals, cnt)
    acc(tab, np.array([coeffs[-1] % p], dtype=np.int64), cnt)
    return cnt

def uv(cnt):
    return [int(cnt[0]-cnt[2]), int(cnt[1]-cnt[2])]

def one_prime(p):
    t0 = time.time()
    tab, om, a, b = chi_table(p); n = nonresidue(p)
    i7 = pow(7,-1,p)
    S0 = th_sums(p, tab, 1, 1, 1)
    S1 = th_sums(p, tab, pow(n,3,p), n*n % p, n)
    # curves: f(1,t) = 1 + 5/7 t + 10/7 t^5 + t^6 ; g(1,t) = 1 + 5/7 t + 15/7 t^5 + t^6 ; B(1,t) = 1 + 15/7 t^2 + 10/7 t^4 + t^6
    Sf = curve_sum(p, tab, [1, 5*i7, 0, 0, 0, 10*i7, 1])
    Sg = curve_sum(p, tab, [1, 5*i7, 0, 0, 0, 15*i7, 1])
    SB = curve_sum(p, tab, [1, 0, 15*i7, 0, 10*i7, 0, 1])
    SBn = curve_sum(p, tab, [1, 0, 15*i7*n, 0, 10*i7*n*n, 0, pow(n,3,p)])
    rec = {'p': p, 'pi': [a, b], 'omega': int(om), 'n': n,
           'S0': uv(S0), 'S1': uv(S1), 'NSH': int(S0[3]), 'NSH_tw': int(S1[3]),
           'Sf': uv(Sf), 'Sg': uv(Sg), 'SB': uv(SB), 'SBn': uv(SBn), 'zf': int(Sf[3]), 'zg': int(Sg[3]), 'zB': int(SB[3]),
           'secs': round(time.time()-t0, 2)}
    return rec

if __name__ == '__main__':
    lo, hi, nproc, out = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
    from sympy import primerange
    ps = [p for p in primerange(lo, hi+1) if p % 3 == 1 and p != 7]
    with Pool(nproc) as pool, open(out, 'a') as fh:
        for rec in pool.imap_unordered(one_prime, ps):
            fh.write(json.dumps(rec) + '\n'); fh.flush()
            print(rec['p'], rec['secs'], flush=True)
