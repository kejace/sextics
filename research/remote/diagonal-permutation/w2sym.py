"""(01)-transposition sector machinery for w2pt (lane diagonal-permutation).

sigma = (01): swaps X0<->X1 and Y0<->Y1 (coordinates 0<->1, 3<->4).  F is sigma-invariant, A is sigma-anti-invariant.
If D is sigma-invariant and (N,Q) certifies A then averaging gives a certificate with N sigma-ANTI-equivariant
(sigma_* N = -N), so for invariant D one may restrict the Koszul unknowns to the anti-equivariant sector and the
E1 lifts to anti-equivariant ones.  (For target B = T0+T1-2T2 the roles flip: B invariant, N equivariant.)

Actions.  On monomials: sigma(e) permutes the exponents.  On fields: (sigma_* w)_i(z) = w_{sigma(i)}(sigma z).
On Koszul coefficients c_ij (K_c = sum_{i<j} c_ij (F_j d_i - F_i d_j)): (sigma.c)_{ij}(z) = c_{sigma(i) sigma(j)}(sigma z)
with c_{ji} = -c_{ij}.
"""
import numpy as np
from fractions import Fraction
import sys
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core
from hodge.core import NV, monomials, pderiv, add, scale

SIG = [1, 0, 2, 4, 3, 5]


def sig_exp(e):
    return (e[1], e[0], e[2], e[4], e[3], e[5])


def sig_poly(P):
    out = {}
    for e, c in P.items():
        e2 = sig_exp(e); out[e2] = out.get(e2, 0) + c
    return {e: v for e, v in out.items() if v}


def sig_field(N):
    """sigma_* N : (sigma_* N)_i (z) = N_{sigma(i)} (sigma z)."""
    return [sig_poly(N[SIG[i]]) for i in range(NV)]


def w_sector_matrix(wbasis, sign):
    """Columns spanning the sign-eigenspace (sign=+1 equivariant, -1 anti) of the sigma action on the w-basis
    (list of (i, m)).  (sigma_* w)_i(z) = w_{sigma i}(sigma z): basis element (i, m) -> (sigma i, sigma m)."""
    idx = {wb: k for k, wb in enumerate(wbasis)}
    seen = set(); cols = []
    for k, (i, m) in enumerate(wbasis):
        if k in seen: continue
        k2 = idx[(SIG[i], sig_exp(m))]
        seen.add(k); seen.add(k2)
        v = np.zeros(len(wbasis))
        if k2 == k:
            if sign == 1: v[k] = 1.0
            else: continue
        else:
            v[k] = 1.0; v[k2] = sign * 1.0
        cols.append(v)
    return np.array(cols).T if cols else np.zeros((len(wbasis), 0))


def c_sector_matrix(cols, sign):
    """cols: list of (i, j, m) Koszul monomial columns (i<j).  Returns the (nK x nK') matrix of the sign-eigenspace.
    (sigma.c)_{ij} = c_{sigma i, sigma j}(sigma z); if sigma i > sigma j the entry picks a minus sign."""
    idx = {c: k for k, c in enumerate(cols)}
    seen = set(); out = []
    for k, (i, j, m) in enumerate(cols):
        if k in seen: continue
        a, b = SIG[i], SIG[j]
        s = 1.0
        if a > b: a, b = b, a; s = -1.0
        k2 = idx[(a, b, sig_exp(m))]
        seen.add(k); seen.add(k2)
        v = np.zeros(len(cols))
        if k2 == k:
            # self-paired: sigma.c = s * c ; belongs to sector sign iff s == sign
            if s == sign: v[k] = 1.0
            else: continue
        else:
            v[k] = 1.0; v[k2] = sign * s
        out.append(v)
    return np.array(out).T if out else np.zeros((len(cols), 0))


def symmetrize_pair(D, N, sign_D=+1):
    """Given an E1 pair (D, N) with N(F) = A D (A anti-invariant), return (D', N') with D' = D + sign_D*sigma(D)
    and N' = N - sign_D * sigma_* N, so that N'(F) = A D'  (uses sigma(A) = -A)."""
    D2 = add(D, scale(sig_poly(D), sign_D))
    sN = sig_field(N)
    N2 = [add(N[i], scale(sN[i], -sign_D)) for i in range(NV)]
    return D2, N2
