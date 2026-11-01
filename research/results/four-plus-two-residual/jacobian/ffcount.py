"""Twisted cubic-character sums for T_B: y^3 = f(x), f(x)=B(1,x)=1+(15/7)x^2+(10/7)x^4+x^6 = gt(x^2),
and E': y^3 = gt(V), gt(V)=1+(15/7)V+(10/7)V^2+V^3, over F_q, q=p^k, via numpy vectorised F_q arithmetic.
Elements of Z[zeta] are pairs (a,b) = a + b*zeta, zeta^2 = -1 - zeta."""
import numpy as np, sys, itertools
from flint import nmod_poly

def irreducible_monic(p, k):
    if k == 1: return [0, 1]
    for coeffs in itertools.product(range(p), repeat=k):
        m = nmod_poly(list(coeffs) + [1], p)
        lead, facs = m.factor()
        if len(facs) == 1 and facs[0][1] == 1 and facs[0][0].degree() == k:
            return list(coeffs) + [1]
    raise RuntimeError

class Fq:
    def __init__(self, p, k):
        self.p, self.k = p, k
        self.m = irreducible_monic(p, k)   # monic, coefficients low->high, length k+1
        self.q = p**k
    def all_elements(self):
        p, k = self.p, self.k
        idx = np.arange(self.q, dtype=np.int64)
        E = np.zeros((self.q, k), dtype=np.int64)
        for j in range(k):
            E[:, j] = (idx // p**j) % p
        return E
    def mul(self, A, Bm):
        p, k = self.p, self.k
        n = A.shape[0]
        C = np.zeros((n, 2*k-1), dtype=np.int64)
        for i in range(k):
            for j in range(k):
                C[:, i+j] += A[:, i] * Bm[:, j]
                if (i*j) % 3 == 2: C %= p   # occasional reduction to avoid overflow (p<=400 so products < 1.6e5; fine)
        C %= p
        m = self.m
        for d in range(2*k-2, k-1, -1):
            c = C[:, d].copy()
            for j in range(k):
                C[:, d-k+j] = (C[:, d-k+j] - c * m[j]) % p
            C[:, d] = 0
        return C[:, :k] % p
    def const(self, c, n):
        A = np.zeros((n, self.k), dtype=np.int64); A[:, 0] = c % self.p; return A
    def power(self, A, e):
        n = A.shape[0]
        R = self.const(1, n); base = A.copy()
        while e:
            if e & 1: R = self.mul(R, base)
            e >>= 1
            if e: base = self.mul(base, base)
        return R
    def eq(self, A, Bm): return np.all(A == Bm, axis=1)
    def is_zero(self, A): return np.all(A == 0, axis=1)

def inv_mod(a, p): return pow(a, -1, p)

def cube_root_of_unity(F, E):
    """Return omega in F_q (as 1xk array) with omega^2+omega+1=0; for p=1 mod 3 take the one in F_p (smallest residue)."""
    p, q = F.p, F.q
    assert q % 3 == 1
    if p % 3 == 1:
        for w in range(2, p):
            if (w*w + w + 1) % p == 0:
                return F.const(w, 1)
    # p = 2 mod 3, k even: find a non-cube z and omega = z^((q-1)/3)
    for z in E[1:]:
        om = F.power(z[None, :], (q-1)//3)
        if not F.eq(om, F.const(1,1))[0]:
            return om
    raise RuntimeError

def twisted_sums(p, k, poly_coeffs_rat):
    """poly in x with rational coefficients (list low->high, as (num,den)), leading coeff 1.
    Returns dict with counts n0,n1,nw,nw2 over affine x in F_q, and S = sum_{x in P^1} iota(chi(f(x))) as (a,b)."""
    F = Fq(p, k); E = F.all_elements(); q = F.q; n = q
    coeffs = [(num * inv_mod(den, p)) % p for (num, den) in poly_coeffs_rat]
    # Horner: f = sum c_i x^i
    val = F.const(coeffs[-1], n)
    for c in coeffs[-2::-1]:
        val = F.mul(val, E)
        val[:, 0] = (val[:, 0] + c) % p
    om = cube_root_of_unity(F, E)
    om2 = F.mul(om, om)
    chi = F.power(val, (q-1)//3)
    z = F.is_zero(val)
    one = F.eq(chi, F.const(1, n)) & ~z
    w1 = F.eq(chi, np.repeat(om, n, axis=0))
    w2 = F.eq(chi, np.repeat(om2, n, axis=0))
    n0, n1, nw, nw2 = int(z.sum()), int(one.sum()), int(w1.sum()), int(w2.sum())
    assert n0 + n1 + nw + nw2 == q, (n0, n1, nw, nw2, q)
    # infinity: leading coeff 1 -> chi=1, contributes +1
    a = n1 + 1 - nw2; b = nw - nw2
    npts = n0 + 3*n1 + 3
    return dict(p=p, k=k, q=q, n0=n0, n1=n1, nw=nw, nw2=nw2, S=(a, b), npts=npts, omega=int(om[0,0]) if k==1 or p%3==1 else om.tolist())

f_B = [(1,1),(0,1),(15,7),(0,1),(10,7),(0,1),(1,1)]      # f(x) = 1 + 15/7 x^2 + 10/7 x^4 + x^6
g_E = [(1,1),(15,7),(10,7),(1,1)]                           # gt(V) = 1 + 15/7 V + 10/7 V^2 + V^3

if __name__ == '__main__':
    p = int(sys.argv[1]); k = int(sys.argv[2])
    print(twisted_sums(p, k, f_B)); print(twisted_sums(p, k, g_E))
