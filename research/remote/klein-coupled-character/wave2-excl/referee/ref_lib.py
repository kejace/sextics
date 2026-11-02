"""Referee replay library for the klein-coupled-character wave-2 exclusion packet (independent of w2lib.py).
Polynomials: dict {exponent 6-tuple: coeff}; coeff = Fraction (exact) or int mod p.
Variables X0,X1,X2,Y0,Y1,Y2; sigma weights W=(1,2,4,0,0,0); Y-parity = (e3+e4+e5) mod 2.
"""
import itertools, sys, json, time
from fractions import Fraction
from math import comb
from flint import nmod_mat, fmpq_mat, fmpq, fmpz_mat

NV = 6
W = (1, 2, 4, 0, 0, 0)
PI = (0, 0, 0, 1, 1, 1)
P1 = 1000000009   # referee prime (not used by the excl / verify lenses)

# 7F with integer coefficients (instance.rec): F = X0^5 X1 + 10/7 X0^3 X2 Y1 Y2 + 5/7 X0 X1^3 Y0 Y1 + X0 X2^5 + X1^5 X2
#   + 15/7 X1 X2^3 Y0 Y2 + Y0^6 + Y1^6 + Y2^6
F7 = {(5,1,0,0,0,0): 7, (3,0,1,0,1,1): 10, (1,3,0,1,1,0): 5, (1,0,5,0,0,0): 7, (0,5,1,0,0,0): 7,
      (0,1,3,1,0,1): 15, (0,0,0,6,0,0): 7, (0,0,0,0,6,0): 7, (0,0,0,0,0,6): 7}
A = {(3,0,0,3,0,0): 1}

def wt(e): return (e[0] + 2*e[1] + 4*e[2]) % 7
def par(e): return (e[3] + e[4] + e[5]) % 2

def monomials(deg, p=None, k=None):
    """exponent tuples of total degree deg (optionally Y-parity p, sigma-weight k), fixed canonical order"""
    if deg < 0: return []
    out = []
    for c in itertools.combinations(range(deg + NV - 1), NV - 1):
        e, prev = [], -1
        for x in c:
            e.append(x - prev - 1); prev = x
        e.append(deg + NV - 2 - prev); e = tuple(e)
        if (p is None or par(e) == p) and (k is None or wt(e) == k): out.append(e)
    return out

def red(P, p):
    if p is None: return {e: Fraction(c) for e, c in P.items() if c}
    out = {}
    for e, c in P.items():
        c = Fraction(c); v = (c.numerator * pow(c.denominator, -1, p)) % p
        if v: out[e] = v
    return out

def add(P, Q, p=None):
    R = dict(P)
    for e, c in Q.items():
        v = R.get(e, 0) + c
        if p: v %= p
        if v: R[e] = v
        else: R.pop(e, None)
    return R

def scale(P, c, p=None):
    if p: c %= p
    R = {}
    for e, v in P.items():
        w = (v*c) % p if p else v*c
        if w: R[e] = w
    return R

def mul(P, Q, p=None):
    R = {}
    for e1, c1 in P.items():
        for e2, c2 in Q.items():
            e = tuple(a+b for a, b in zip(e1, e2)); v = R.get(e, 0) + c1*c2
            if p: v %= p
            R[e] = v
    return {e: v for e, v in R.items() if v}

def shift(P, m, c=1, p=None):
    R = {}
    for e, v in P.items():
        w = (v*c) % p if p else v*c
        if w: R[tuple(a+b for a, b in zip(e, m))] = w
    return R

def pderiv(P, i, p=None):
    R = {}
    for e, c in P.items():
        if e[i]:
            e2 = list(e); e2[i] -= 1; v = c*e[i]
            if p: v %= p
            if v: R[tuple(e2)] = v
    return R

def degree(P): return max(sum(e) for e in P) if P else -1

def Fpolys(p=None):
    F = red(F7, p)                      # 7F (scalar multiple: same ideal, same certificates up to scaling)
    Fi = [pderiv(F, i, p) for i in range(NV)]
    H = [[pderiv(Fi[i], j, p) for j in range(NV)] for i in range(NV)]
    return F, Fi, H

def mat(rows, cols, entries, p):
    """entries: dict {(r,c): v}"""
    if p:
        flat = [0]*(rows*cols)
        for (r, c), v in entries.items(): flat[r*cols + c] = v % p
        return nmod_mat(rows, cols, flat, p)
    M = fmpq_mat(rows, cols)
    for (r, c), v in entries.items(): M[r, c] = fmpq(v.numerator, v.denominator) if isinstance(v, Fraction) else fmpq(v)
    return M

def rank(M): return M.rank()

def J_generators(Fi, d, tp, tk, p=None):
    """generators F_i * m of J_d in the (parity tp, weight tk) sector"""
    gens = []
    for i in range(NV):
        for m in monomials(d-5):
            if (par(m) + PI[i]) % 2 == tp and (wt(m) - W[i]) % 7 == tk:
                gens.append(((i, m), shift(Fi[i], m, 1, p)))
    return gens

def sector_dims(d, tp, tk, Fi, p):
    mons = monomials(d, tp, tk); idx = {e: r for r, e in enumerate(mons)}
    gens = J_generators(Fi, d, tp, tk, p)
    ent = {}
    for c, (_, g) in enumerate(gens):
        for e, v in g.items(): ent[(idx[e], c)] = v
    M = mat(len(mons), len(gens), ent, p) if gens else None
    rJ = rank(M) if gens else 0
    return len(mons), rJ, mons, idx, gens, M

def colon_dim(d, tp, tk, mult, Fi, p):
    """dim (J : mult)_d in sector (tp,tk) and dim J_d there; mult homogeneous of parity pm, weight km"""
    pm = (par(next(iter(mult))) ) ; km = wt(next(iter(mult)))
    dm = degree(mult)
    nS, rJ, mons, idx, gens, _ = sector_dims(d, tp, tk, Fi, p)
    # target sector of mult * S_d
    tp2, tk2 = (tp + pm) % 2, (tk + km) % 7
    mons2 = monomials(d + dm, tp2, tk2); idx2 = {e: r for r, e in enumerate(mons2)}
    gens2 = J_generators(Fi, d + dm, tp2, tk2, p)
    ent = {}
    for c, (_, g) in enumerate(gens2):
        for e, v in g.items(): ent[(idx2[e], c)] = v
    nJ = len(gens2)
    for c, m in enumerate(mons):
        g = shift(mult, m, 1, p)
        for e, v in g.items(): ent[(idx2[e], nJ + c)] = v
    M = mat(len(mons2), nJ + len(mons), ent, p)
    rAll = rank(M)
    rJ2 = rank(mat(len(mons2), nJ, {k: v for k, v in ent.items() if k[1] < nJ}, p)) if nJ else 0
    colon = len(mons) - (rAll - rJ2)
    return {'d': d, 'par': tp, 'k': tk, 'dimS': len(mons), 'dimJ': rJ, 'dim_colon': colon, 'kernel_mod_J': colon - rJ}

class Reducer:
    """normal form modulo the subspace spanned by given polynomials (rows = monomials), exact or mod p"""
    def __init__(self, mons, gens, p):
        self.mons, self.p = mons, p; self.idx = {e: r for r, e in enumerate(mons)}; nm = len(mons)
        ent = {}
        for r, g in enumerate(gens):
            for e, v in g.items(): ent[(r, self.idx[e])] = v
        G = mat(len(gens), nm, ent, p) if gens else None
        if gens:
            R, rk = G.rref() if p else (G.rref()[0], G.rank())
            if p: R, rk = R, rk
            self.R = R; self.rk = rk
            piv = []
            for i in range(rk):
                for j in range(nm):
                    if R[i, j] != 0: piv.append(j); break
            self.piv = piv
        else:
            self.R, self.rk, self.piv = None, 0, []
        self.nonpiv = [j for j in range(nm) if j not in set(self.piv)]
    def reduce_rows(self, V):
        """V: matrix (rows x nm); returns V reduced (pivot coords zero), as matrix"""
        if self.rk == 0: return V
        nm = len(self.mons)
        # V - V[:,piv] * R_rows(normalized so pivot entries = 1 -- rref already)
        if self.p:
            Vp = nmod_mat(V.nrows(), self.rk, [int(V[r, j]) for r in range(V.nrows()) for j in self.piv], self.p)
            Rr = nmod_mat(self.rk, nm, [int(self.R[i, j]) for i in range(self.rk) for j in range(nm)], self.p)
            return V - Vp * Rr
        Vp = fmpq_mat(V.nrows(), self.rk)
        for r in range(V.nrows()):
            for c, j in enumerate(self.piv): Vp[r, c] = V[r, j]
        Rr = fmpq_mat(self.rk, nm)
        for i in range(self.rk):
            for j in range(nm): Rr[i, j] = self.R[i, j]
        return V - Vp * Rr

def hess_relaxation(d, pars, chars, Fi, H, p, gauge=True):
    """Gauge-fixed Hessian relaxation for D = sum w_i F_i (deg w_i = d-5), D-parity in pars, D-weight in chars.
    Equations: coordinates of w^T Hess(7F) w modulo J_{2d-6} (valid for 2d-6 <= 8 where (J:A) = J).
    Returns vars (list of (i,m)), eqs (list of dict {(a,b): coeff}, a<=b)."""
    assert 2*d - 6 <= 8
    pars, chars = set(pars), set(chars)
    unknowns = [(i, m) for i in range(NV) for m in monomials(d-5) if (par(m)+PI[i]) % 2 in pars and (wt(m)-W[i]) % 7 in chars]
    removed = set()
    if gauge and d >= 6:
        for mp in monomials(d-6):
            if par(mp) in pars and wt(mp) in chars:
                m0 = list(mp); m0[0] += 1; removed.add((0, tuple(m0)))
    vars_ = [u for u in unknowns if u not in removed]
    assert len(removed) == len(unknowns) - len(vars_)
    n = len(vars_); e2 = 2*d - 6
    def dsec(u): i, m = u; return ((par(m)+PI[i]) % 2, (wt(m)-W[i]) % 7)
    groups = {}
    for a in range(n):
        for b in range(a, n):
            sa, sb = dsec(vars_[a]), dsec(vars_[b])
            groups.setdefault(((sa[0]+sb[0]) % 2, (sa[1]+sb[1]) % 7), []).append((a, b))
    eqs = []
    for (tp, tk), pairs in groups.items():
        mons = monomials(e2, tp, tk); idx = {e: r for r, e in enumerate(mons)}
        gens = [g for _, g in J_generators(Fi, e2, tp, tk, p)]
        redu = Reducer(mons, gens, p)
        ent = {}
        for r, (a, b) in enumerate(pairs):
            ia, ma = vars_[a]; ib, mb = vars_[b]
            Pab = H[ia][ib]
            if not Pab: continue
            mm = tuple(x+y for x, y in zip(ma, mb)); c = 1 if a == b else 2
            for e, v in Pab.items():
                ee = tuple(x+y for x, y in zip(e, mm)); key = (r, idx[ee])
                ent[key] = ent.get(key, 0) + c*v
        V = mat(len(pairs), len(mons), ent, p)
        NF = redu.reduce_rows(V)
        for j in redu.nonpiv:
            eq = {}
            for r in range(len(pairs)):
                v = NF[r, j]
                if v != 0: eq[pairs[r]] = (int(v) if p else Fraction(int(v.p), int(v.q)))
            if eq: eqs.append(eq)
    return {'d': d, 'pars': sorted(pars), 'chars': sorted(chars), 'vars': vars_, 'n': n, 'eqs': eqs, 'n_eqs': len(eqs), 'n_gauge': len(removed)}

def poly_str(eq, names, p):
    terms = []
    for (a, b), c in sorted(eq.items()):
        if p: cs = str(c % p)
        else:
            cs = str(c.numerator) if c.denominator == 1 else f'({c.numerator}/{c.denominator})'
        terms.append(f'{cs}*{names[a]}^2' if a == b else f'{cs}*{names[a]}*{names[b]}')
    return '+'.join(terms).replace('+-', '-')

def write_singular(sysd, path, p, extra=''):
    names = [f'u{i}' for i in range(sysd['n'])]
    body = f"ring R = {p if p else 0}, ({','.join(names)}), dp;\nideal I = " + ',\n'.join(poly_str(e, names, p) for e in sysd['eqs']) + ';\n' + extra
    open(path, 'w').write(body)

def write_msolve(sysd, path, p):
    names = [f'u{i}' for i in range(sysd['n'])]
    # msolve: integer coefficients; clear denominators for char 0
    lines = []
    for e in sysd['eqs']:
        if p: lines.append(poly_str(e, names, p))
        else:
            den = 1
            for c in e.values():
                from math import gcd; den = den * c.denominator // gcd(den, c.denominator)
            ee = {k: Fraction(c*den) for k, c in e.items()}
            lines.append(poly_str(ee, names, None))
    open(path, 'w').write(','.join(names) + f'\n{p if p else 0}\n' + ',\n'.join(lines) + '\n')

def ci_hilbert(d):
    """dim (S/J)_d for a regular sequence of six quintics"""
    tot = 0
    for j in range(0, 7):
        if d - 5*j < 0: break
        tot += (-1)**j * comb(6, j) * comb(d - 5*j + 5, 5)
    return tot
