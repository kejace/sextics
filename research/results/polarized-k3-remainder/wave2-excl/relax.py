"""Exact Koszul-parametrized relaxation builder (UV frame, integer-scaled G = 14 Fhat, T = 64 Ahat).
D = sum_i a_i P_i (P_i = d_i G, deg a_i = d-5), normalized numerator N' = T a (+ Koszul K),
E2(D,N') = D div N' - N'(D)  (must lie in (G) + B_D(Koszul) subset L = (G) + J(G)^2 + (R_ijk)).
Parameter-linear polys: {zmon: {param: int}}; parameter-quadratic: {zmon: {(p1<=p2): int}}."""
import sys, time, itertools
sys.path.insert(0, '.')
from k3lib import *

def integer_forms():
    Fh, Ah, q = load_uv()
    G = {e: int(c * 14) for e, c in Fh.items()}
    Tt = {e: int(c * 64) for e, c in Ah.items()}
    assert all(Fr(c * 14).denominator == 1 for c in Fh.values()) and all(Fr(c * 64).denominator == 1 for c in Ah.values())
    return G, Tt

def ideriv(P, i):
    out = {}
    for e, c in P.items():
        if e[i]:
            e2 = list(e); e2[i] -= 1; out[tuple(e2)] = c * e[i]
    return out

def imul(P, Q):
    out = {}
    for e1, c1 in P.items():
        for e2, c2 in Q.items():
            e = tuple(x + y for x, y in zip(e1, e2)); out[e] = out.get(e, 0) + c1 * c2
    return {e: v for e, v in out.items() if v}

def lin_scale_mon(poly, m, pidx, coef=1):
    out = {}
    for e, c in poly.items():
        e2 = tuple(e[t] + m[t] for t in range(NV)); out.setdefault(e2, {})[pidx] = out.get(e2, {}).get(pidx, 0) + c * coef
    return out

def lin_add(A_, B_):
    out = {e: dict(v) for e, v in A_.items()}
    for e, v in B_.items():
        dst = out.setdefault(e, {})
        for pi, c in v.items():
            dst[pi] = dst.get(pi, 0) + c
            if dst[pi] == 0: del dst[pi]
        if not dst: del out[e]
    return out

def lin_deriv(A_, i):
    out = {}
    for e, v in A_.items():
        if e[i]:
            e2 = list(e); e2[i] -= 1; out[tuple(e2)] = {pi: c * e[i] for pi, c in v.items()}
    return out

def quad_mul(A_, B_):
    out = {}
    for e1, v1 in A_.items():
        for e2, v2 in B_.items():
            e = tuple(x + y for x, y in zip(e1, e2)); dst = out.setdefault(e, {})
            for p1, c1 in v1.items():
                for p2, c2 in v2.items():
                    key = (p1, p2) if p1 <= p2 else (p2, p1); dst[key] = dst.get(key, 0) + c1 * c2
    return out

def quad_add(A_, B_, sign=1):
    out = {e: dict(v) for e, v in A_.items()}
    for e, v in B_.items():
        dst = out.setdefault(e, {})
        for key, c in v.items():
            dst[key] = dst.get(key, 0) + sign * c
            if dst[key] == 0: del dst[key]
        if not dst: del out[e]
    return out

class Relaxation:
    def __init__(self, d, with_koszul=False, verbose=True):
        self.d = d
        G, Tt = integer_forms(); self.G, self.T = G, Tt
        self.P = [ideriv(G, i) for i in range(NV)]
        self.H = [[ideriv(self.P[i], j) for j in range(NV)] for i in range(NV)]
        mons_a = monomials(d - 5) if d >= 5 else []
        self.params = []; self.aidx = {}; self.kidx = {}
        for i in range(NV):
            for m in mons_a:
                self.aidx[(i, m)] = len(self.params); self.params.append(('a', i, m))
        self.na = len(self.params)
        if with_koszul:
            mons_k = monomials(d - 4)
            for i in range(NV):
                for j in range(i + 1, NV):
                    for m in mons_k:
                        self.kidx[(i, j, m)] = len(self.params); self.params.append(('k', i, j, m))
        self.npar = len(self.params)
        t0 = time.time()
        Dl = {}
        for (i, m), pid in self.aidx.items():
            Dl = lin_add(Dl, lin_scale_mon(self.P[i], m, pid))
        Nl = [{} for _ in range(NV)]
        for (i, m), pid in self.aidx.items():
            Nl[i] = lin_add(Nl[i], lin_scale_mon(Tt, m, pid))
        for (i, j, m), pid in self.kidx.items():
            Nl[i] = lin_add(Nl[i], lin_scale_mon(self.P[j], m, pid))
            Nl[j] = lin_add(Nl[j], lin_scale_mon(self.P[i], m, pid, -1))
        divN = {}
        for i in range(NV):
            divN = lin_add(divN, lin_deriv(Nl[i], i))
        E2 = quad_mul(Dl, divN)
        for i in range(NV):
            E2 = quad_add(E2, quad_mul(Nl[i], lin_deriv(Dl, i)), -1)
        self.Dl, self.Nl, self.E2 = Dl, Nl, E2
        self.deg_out = 2 * d
        if verbose:
            print(f'[relax d={d}] params={self.npar} (a {self.na}) E2 monomials={len(E2)} built {time.time()-t0:.1f}s', flush=True)

    def quad_keys(self):
        keys = set()
        for v in self.E2.values():
            keys.update(v.keys())
        return sorted(keys)

    def L_generators(self, deg, with_R=True, with_J2=True):
        """Integer polys spanning L_deg = (G)_deg + (J(G)^2)_deg + (R_ijk)_deg."""
        gens = []
        if deg >= 6:
            for m in monomials(deg - 6):
                gens.append(('G', m, imul(self.G, {m: 1})))
        if with_J2 and deg >= 10:
            PP = {}
            for i in range(NV):
                for j in range(i, NV):
                    PP[(i, j)] = imul(self.P[i], self.P[j])
            for (i, j), pp in PP.items():
                for m in monomials(deg - 10):
                    gens.append(('PP', (i, j, m), imul(pp, {m: 1})))
        if with_R and deg >= 9:
            Rs = {}
            for i in range(NV):
                for j in range(NV):
                    for k in range(j + 1, NV):
                        r = {}
                        for e, c in imul(self.P[k], self.H[i][j]).items(): r[e] = r.get(e, 0) + c
                        for e, c in imul(self.P[j], self.H[i][k]).items(): r[e] = r.get(e, 0) - c
                        r = {e: c for e, c in r.items() if c}
                        if r: Rs[(i, j, k)] = r
            for key, r in Rs.items():
                for m in monomials(deg - 9):
                    gens.append(('R', key + (m,), imul(r, {m: 1})))
        return gens
