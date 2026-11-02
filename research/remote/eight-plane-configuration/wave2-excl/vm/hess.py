"""Hessian relaxation for Jacobian denominators D = w(F), w in a linear space W of vector fields of degree d-5.

Necessary condition (exact, see README): A*Hess_F(w,w) in J_{2d}; for d <= 7 (Ann_R(A)_{2d-6} = 0) this is
Hess_F(w,w) in J_{2d-6}.  Modulo T = {h*Euler} the condition depends only on the class of w.
This module builds, for a user-given list of basis fields (exact integer coefficients), the quadratic forms
  q_l(x) = l( NF_J( Hess_F(sum x_k w_k, sum x_k w_k) ) ),  l running over a basis of (S_{2d-6}/J_{2d-6})^vee  (mod p),
writes them as an msolve input, and (optionally) computes Macaulay-matrix ranks with flint as an independent check.
'A-twisted' variant (d >= 8): coordinates of A*Hess(w,w) in S_{2d}/J_{2d}.
"""
import sys, time, json, itertools, subprocess, os
sys.path.insert(0, '/tmp/hodge/results/eight-plane-configuration/wave2-excl/vm')
from w2lib import *
from kred import field_apply, unit_field
from flint import nmod_mat

Fij_int = [[pderiv(Fi_int[i], j) for j in range(6)] for i in range(6)]

def hess_bilinear(w1, w2):
    """Hess_F(w1,w2) = sum_{i,j} w1_i w2_j F_ij (exact)."""
    out = {}
    for i in range(6):
        if not w1[i]: continue
        for j in range(6):
            if not w2[j]: continue
            out = add(out, mul(mul(w1[i], w2[j]), Fij_int[i][j]))
    return out

def euler_multiples(deg):
    """basis of T = {h * Euler : h in S_deg} as fields."""
    out = []
    for m in monomials(deg):
        w = [{} for _ in range(6)]
        for i in range(6):
            e = list(m); e[i] += 1
            w[i] = {tuple(e): 1}
        out.append(w)
    return out

def field_vec(w, degw):
    """flatten a field of coefficient degree degw to a vector over the basis (i, monomial)."""
    mons = monomials(degw)
    idx = {m: k for k, m in enumerate(mons)}
    v = [0] * (6 * len(mons))
    for i in range(6):
        for m, c in w[i].items():
            v[i * len(mons) + idx[m]] = c
    return v

class HessianSystem:
    def __init__(self, basis, d, p, twisted=False, verbose=False):
        """basis: list of exact fields (degree d-5). Quadrics in len(basis) variables."""
        t0 = time.time()
        self.p, self.d, self.basis = p, d, basis
        n = len(basis)
        self.n = n
        deg = 2 * d - 6 + (6 if twisted else 0)
        self.deg = deg
        red = jacobian._Reducer(F, deg, p, None)       # all parities
        self.red = red
        self.rank_J_p = red.jacobian_dim
        self.J_exact = len(monomials(deg)) - hilbert_CI(deg)
        std = red.standard
        sidx = {m: k for k, m in enumerate(std)}
        # bilinear forms H_kl = NF(Hess(w_k, w_l)) (symmetric)
        self.pairs = [(k, l) for k in range(n) for l in range(k, n)]
        H = {}
        for (k, l) in self.pairs:
            h = hess_bilinear(basis[k], basis[l])
            if twisted: h = mul(Aint, h)
            nf = red.normal_form(reduce_mod(h, p))
            H[(k, l)] = nf
        self.H = H
        # quadric coefficient matrix: rows = standard monomials (coordinates), cols = pairs; entry = coeff of x_k x_l
        # (for k<l the monomial x_k x_l appears from H_kl + H_lk = 2 H_kl; we store 2*H_kl for k<l, H_kk for k=l)
        R, C = len(std), len(self.pairs)
        flat = [0] * (R * C)
        for j, (k, l) in enumerate(self.pairs):
            mult = 1 if k == l else 2
            for m, c in H[(k, l)].items():
                flat[sidx[m] * C + j] = (c * mult) % p
        self.Qmat = nmod_mat(R, C, flat, p)      # each row = one quadric (coefficients in the x_k x_l basis)
        self.rank_lin = self.Qmat.rank()
        self.dimSym2 = C
        self.t_build = round(time.time() - t0, 1)
        if verbose:
            print(f'[hess] d={d} n={n} deg={deg} R_deg={len(std)} rank_p J={self.rank_J_p} exact={self.J_exact} '
                  f'quadrics(indep)={self.rank_lin} dimSym2={C} ({self.t_build}s)', file=sys.stderr, flush=True)

    def independent_quadrics(self):
        """rref rows of Qmat -> list of dict {(k,l): coeff}."""
        Qr, rk = self.Qmat.rref()
        out = []
        for i in range(rk):
            q = {}
            for j, (k, l) in enumerate(self.pairs):
                v = int(Qr[i, j])
                if v: q[(k, l)] = v
            out.append(q)
        return out

    def write_msolve(self, path):
        quads = self.independent_quadrics()
        names = [f'x{k}' for k in range(self.n)]
        with open(path, 'w') as f:
            f.write(','.join(names) + '\n')
            f.write(f'{self.p}\n')
            lines = []
            for q in quads:
                terms = []
                for (k, l), v in q.items():
                    mon = f'{names[k]}^2' if k == l else f'{names[k]}*{names[l]}'
                    terms.append(f'{v}*{mon}')
                lines.append('+'.join(terms))
            f.write(',\n'.join(lines) + '\n')
        return len(quads)

    def macaulay_rank(self, k):
        """rank mod p of the degree-k Macaulay matrix of the independent quadrics (flint). Returns (rank, dim S_k)."""
        quads = self.independent_quadrics()
        n = self.n
        monk = list(itertools.combinations_with_replacement(range(n), k))
        monk2 = list(itertools.combinations_with_replacement(range(n), k - 2))
        idx = {m: i for i, m in enumerate(monk)}
        R, C = len(quads) * len(monk2), len(monk)
        flat = [0] * (R * C)
        r = 0
        for q in quads:
            for m2 in monk2:
                base = r * C
                for (a, b), v in q.items():
                    mon = tuple(sorted(m2 + (a, b)))
                    flat[base + idx[mon]] = (flat[base + idx[mon]] + v) % self.p
                r += 1
        M = nmod_mat(R, C, flat, self.p)
        return M.rank(), C, R

def run_msolve(path, out, threads=1, timeout=None):
    cmd = ['msolve', '-g', '1', '-t', str(threads), '-f', path, '-o', out]
    t0 = time.time()
    subprocess.run(cmd, check=False, timeout=timeout, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return round(time.time() - t0, 1)

def hilbert_from_leading(out_path, n, maxdeg=12):
    """Hilbert function of S/in(I) from msolve -g 1 output (leading monomials)."""
    txt = open(out_path).read()
    # msolve -g 1 prints e.g. "#Leading ideal ... \n[x0^2, x0*x1, ...]:" ; parse monomials
    import re
    body = txt[txt.find('['):txt.rfind(']') + 1] if '[' in txt else ''
    mons = []
    for term in body.strip('[]').split(','):
        term = term.strip()
        if not term: continue
        e = [0] * n
        if term != '1':
            for fac in term.split('*'):
                fac = fac.strip()
                if '^' in fac:
                    v, k = fac.split('^'); e[int(v[1:])] += int(k)
                else:
                    e[int(fac[1:])] += 1
        else:
            return [0] * (maxdeg + 1), mons
        mons.append(tuple(e))
    # count standard monomials of each degree
    hf = []
    for k in range(maxdeg + 1):
        cnt = 0
        for m in itertools.combinations_with_replacement(range(n), k):
            e = [0] * n
            for v in m: e[v] += 1
            if not any(all(e[i] >= lm[i] for i in range(n)) for lm in mons):
                cnt += 1
        hf.append(cnt)
        if cnt == 0: break
    return hf, mons
