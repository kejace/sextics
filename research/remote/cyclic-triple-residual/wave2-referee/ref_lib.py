"""Referee replay library for lane cyclic-triple-residual (q = 1, A = X0^3 Y0^3).

Written from scratch for the wave-2 referee pass: no imports from the hodge engine, from w2lib.py or from lib.py.
Polynomials are dicts {exponent 6-tuple: int}.  Variables X0,X1,X2,Y0,Y1,Y2 = indices 0..5.
Fz = 7*F (integer form); J(Fz) = J(F) as ideals over Q and over F_p for p != 7.
"""
import itertools, json, time, os, re, subprocess
from fractions import Fraction
from flint import fmpz_mat, nmod_mat, fmpq_mat, fmpq, fmpz

# --- the exact target (data/cyclic-triple-residual/instance.rec), cleared by 7 ---------------------------------
FZ = {(6,0,0,0,0,0): 7, (0,6,0,0,0,0): 7, (0,0,6,0,0,0): 7, (0,0,0,6,0,0): 7, (0,0,0,0,6,0): 7, (0,0,0,0,0,6): 7,
      (3,1,0,1,1,0): 5, (3,0,3,0,0,0): 5, (2,0,0,1,1,2): 5, (1,1,0,0,3,1): 5, (1,0,3,0,1,1): 5, (0,1,3,1,1,0): 5}
AZ = {(3,0,0,3,0,0): 1}
Q = 1
# diagonal symmetry group G (order 18): weights mod 6 of the generators mu3 and s
MU3 = (0, 0, 2, 0, 0, 0)
SW = (1, 5, 1, 2, 2, 0)
W = [(MU3[i], SW[i]) for i in range(6)]          # character of the coordinate Z_i

def chi(m):
    return (sum(a * b for a, b in zip(m, MU3)) % 6, sum(a * b for a, b in zip(m, SW)) % 6)

def cadd(a, b): return ((a[0] + b[0]) % 6, (a[1] + b[1]) % 6)
def csub(a, b): return ((a[0] - b[0]) % 6, (a[1] - b[1]) % 6)
def cneg(a): return ((-a[0]) % 6, (-a[1]) % 6)

CHI_A = chi(next(iter(AZ)))                       # (0, 3)
SECTORS = sorted({(a, b) for a in (0, 2, 4) for b in range(6)})

def check_setup():
    assert all(chi(m) == (0, 0) for m in FZ), 'F is not G-invariant'
    assert CHI_A == (0, 3)
    assert len(SECTORS) == 18

# --- polynomial utilities -------------------------------------------------------------------------------------
def mons(n):
    """All exponent 6-tuples of total degree n (n < 0 -> empty)."""
    if n < 0: return []
    out = []
    for c in itertools.combinations_with_replacement(range(6), n):
        e = [0] * 6
        for i in c: e[i] += 1
        out.append(tuple(e))
    return out

_MONS = {}
def mons_c(n):
    if n not in _MONS: _MONS[n] = mons(n)
    return _MONS[n]

def mons_sector(n, c):
    return [m for m in mons_c(n) if chi(m) == c]

def padd(P, Qp, s=1):
    R = dict(P)
    for e, v in Qp.items():
        R[e] = R.get(e, 0) + s * v
        if R[e] == 0: del R[e]
    return R

def pmul(P, Qp):
    R = {}
    for e1, v1 in P.items():
        for e2, v2 in Qp.items():
            e = tuple(a + b for a, b in zip(e1, e2))
            R[e] = R.get(e, 0) + v1 * v2
    return {e: v for e, v in R.items() if v}

def pscale(P, s): return {e: v * s for e, v in P.items()} if s else {}

def pderiv(P, i):
    R = {}
    for e, v in P.items():
        if e[i]:
            e2 = list(e); e2[i] -= 1
            R[tuple(e2)] = v * e[i]
    return R

def mono(m, c=1): return {tuple(m): c}

def pdeg(P): return sum(next(iter(P))) if P else -1

DF = [pderiv(FZ, i) for i in range(6)]
DDF = [[pderiv(DF[i], k) for k in range(6)] for i in range(6)]

def ci_hilbert(n):
    """dim R_n for R = S/J(F), F smooth sextic in 6 variables: coefficient of t^n in (1+t+..+t^4)^6."""
    c = [1]
    for _ in range(6):
        new = [0] * (len(c) + 4)
        for i, v in enumerate(c):
            for j in range(5): new[i + j] += v
        c = new
    return c[n] if 0 <= n < len(c) else 0

def binom(n, k):
    from math import comb
    return comb(n, k)

def dimS(n): return binom(n + 5, 5)

# --- J(F) in one sector -----------------------------------------------------------------------------------------
def jac_gens(n, c):
    """Integer polynomials m * dF_i spanning J_n in sector c (chi(m) = c + w_i)."""
    out = []
    if n < 5: return out
    for i in range(6):
        for m in mons_sector(n - 5, cadd(c, W[i])):
            out.append(pmul(mono(m), DF[i]))
    return out

def to_rows(polys, idx, p=None):
    C = len(idx)
    flat = [0] * (len(polys) * C)
    for r, P in enumerate(polys):
        for e, v in P.items():
            flat[r * C + idx[e]] = v if p is None else v % p
    return flat

def rank_rows(polys, idx, p=None):
    if not polys: return 0
    C = len(idx)
    flat = to_rows(polys, idx, p)
    if p is None:
        return fmpz_mat(len(polys), C, flat).rank()
    return nmod_mat(len(polys), C, flat, p).rank()

class Reducer:
    """Row-reduced J_n in sector c (mod p, or exact over Q with p=None); normal forms on standard monomials."""
    def __init__(self, n, c, p):
        self.n, self.c, self.p = n, c, p
        self.mons = mons_sector(n, c)
        self.idx = {m: i for i, m in enumerate(self.mons)}
        gens = jac_gens(n, c)
        C = len(self.mons)
        self.rank = 0; self.piv = []; self.Rr = None
        if gens:
            flat = to_rows(gens, self.idx, p)
            if p is None:
                M = fmpq_mat(len(gens), C, flat)
                Mr, rk = M.rref()
            else:
                M = nmod_mat(len(gens), C, flat, p)
                Mr, rk = M.rref()
            self.rank = rk
            piv = []
            for i in range(rk):
                j = 0
                while j < C and Mr[i, j] == 0: j += 1
                piv.append(j)
            self.piv = piv
            if rk:
                if p is None:
                    self.Rr = fmpq_mat(rk, C, [Mr[i, j] for i in range(rk) for j in range(C)])
                else:
                    self.Rr = nmod_mat(rk, C, [int(Mr[i, j]) for i in range(rk) for j in range(C)], p)
        ps = set(self.piv)
        self.standard = [j for j in range(C) if j not in ps]
        self.dimR = len(self.standard)
    def nf(self, polys):
        """Normal forms of polys (all in this sector): list of dicts {standard-monomial position: value}."""
        if not polys: return []
        k, C, p = len(polys), len(self.mons), self.p
        flat = to_rows(polys, self.idx, p)
        if p is None:
            V = fmpq_mat(k, C, flat)
            if self.rank:
                Vp = fmpq_mat(k, self.rank, [V[i, j] for i in range(k) for j in self.piv])
                V = V - Vp * self.Rr
            out = []
            for i in range(k):
                out.append({s: V[i, j] for s, j in enumerate(self.standard) if V[i, j] != 0})
            return out
        V = nmod_mat(k, C, flat, p)
        if self.rank:
            Vp = nmod_mat(k, self.rank, [int(V[i, j]) for i in range(k) for j in self.piv], p)
            V = V - Vp * self.Rr
        out = []
        for i in range(k):
            out.append({s: int(V[i, j]) for s, j in enumerate(self.standard) if int(V[i, j]) != 0})
        return out

# --- E1 kernel: Ann([A])_d mod J_d, by sector --------------------------------------------------------------------
def e1_kernel(d, p=None, verbose=False):
    """Returns dict with per-sector kernel dims, totals, and the J-rank checks against the c.i. values."""
    tot_ker = 0; tot_ann = 0; rankJd = 0; rankJd6 = 0; per = {}
    for c in SECTORS:
        Sd = mons_sector(d, c)
        c2 = cadd(c, CHI_A)
        tgt = mons_sector(d + 6, c2)
        idx = {m: i for i, m in enumerate(tgt)}
        gens = jac_gens(d + 6, c2)
        rJ = rank_rows(gens, idx, p)
        cols = [pmul(AZ, mono(m)) for m in Sd]
        rAug = rank_rows(gens + cols, idx, p)
        ker = len(Sd) - (rAug - rJ)
        # dim J_d in sector c
        idxd = {m: i for i, m in enumerate(Sd)}
        rJd = rank_rows(jac_gens(d, c), idxd, p)
        ann = ker - rJd
        per[c] = {'dimS': len(Sd), 'kerA': ker, 'dimJ': rJd, 'ann': ann}
        tot_ker += ker; tot_ann += ann; rankJd += rJd; rankJd6 += rJ
        if verbose: print(f'  d={d} sector {c}: dimS={len(Sd)} kerA={ker} dimJ={rJd} ann={ann}', flush=True)
    ciJd = dimS(d) - ci_hilbert(d); ciJd6 = dimS(d + 6) - ci_hilbert(d + 6)
    return {'d': d, 'mode': 'exact' if p is None else f'mod {p}', 'kerA_total': tot_ker, 'ann_total': tot_ann,
            'dimJ_d': rankJd, 'ci_dimJ_d': ciJd, 'dimJ_d6': rankJd6, 'ci_dimJ_d6': ciJd6,
            'ci_ok': (rankJd == ciJd and rankJd6 == ciJd6), 'sectors': {str(k): v for k, v in per.items()}}

# --- membership tests -----------------------------------------------------------------------------------------
def in_J_exact(P):
    """Is the (sector-pure) integer polynomial P in J_n over Q?  Returns (in_J, rankJ, rankAug)."""
    n = pdeg(P); c = chi(next(iter(P)))
    assert all(chi(e) == c for e in P), 'P must be sector-pure'
    tgt = mons_sector(n, c); idx = {m: i for i, m in enumerate(tgt)}
    gens = jac_gens(n, c)
    r1 = rank_rows(gens, idx); r2 = rank_rows(gens + [P], idx)
    return (r1 == r2, r1, r2)

# --- fixed-D single cleared identity (formal contract, q = 1), G-sector reduced --------------------------------
def fixedD_columns(Dz, sector_only=True):
    """Columns of  E(N,Q) = A D^2 - D sum N_i F_i + F (D div N - N(D)) - F^2 Q = 0  written as  sum col*unknown = -A D^2.
    N_i in S_{d+1}, Q in S_{2d-6}.  If sector_only, keep the unknowns whose column lies in the character of A D^2."""
    d = pdeg(Dz); cD = chi(next(iter(Dz)))
    assert all(chi(e) == cD for e in Dz), 'D must be G-semi-invariant for the sector reduction'
    cR = cadd(CHI_A, cadd(cD, cD))
    dD = [pderiv(Dz, i) for i in range(6)]
    F2 = pmul(FZ, FZ)
    cols, meta = [], []
    for i in range(6):
        for m in mons_c(d + 1):
            if sector_only and cadd(chi(m), cneg(W[i])) != cadd(CHI_A, cD): continue
            col = pscale(pmul(Dz, pmul(mono(m), DF[i])), -1)
            t2 = padd(pmul(Dz, pderiv(mono(m), i)), pmul(mono(m), dD[i]), -1)
            col = padd(col, pmul(FZ, t2))
            cols.append(col); meta.append(('N', i, m))
    for m in mons_c(2 * d - 6):
        if sector_only and chi(m) != cR: continue
        cols.append(pscale(pmul(mono(m), F2), -1)); meta.append(('Q', -1, m))
    rhs = pscale(pmul(AZ, pmul(Dz, Dz)), -1)
    return cols, meta, rhs

def exact_ranks(cols, rhs):
    rows = {}
    for col in cols:
        for e in col: rows.setdefault(e, len(rows))
    for e in rhs: rows.setdefault(e, len(rows))
    R, C = len(rows), len(cols)
    flat = [0] * (R * (C + 1))
    for j, col in enumerate(cols):
        for e, v in col.items(): flat[rows[e] * (C + 1) + j] = v
    for e, v in rhs.items(): flat[rows[e] * (C + 1) + C] = v
    M = fmpz_mat(R, C + 1, flat)
    r_aug = M.rank()
    M0 = fmpz_mat(R, C, [flat[r * (C + 1) + j] for r in range(R) for j in range(C)])
    r0 = M0.rank()
    return {'rows': R, 'cols': C, 'rank': r0, 'rank_aug': r_aug, 'consistent': r0 == r_aug}

# --- Hessian relaxation quadrics ---------------------------------------------------------------------------------
def ypar(m): return (m[3] + m[4] + m[5]) % 2

def theta_vars(d, sector=None, parity=None):
    """Coordinates (i, m) of G_i; sector: G-character of D (= chi(m F_i)); parity: Y-parity of D (0 even, 1 odd)."""
    th = []
    for i in range(6):
        for m in mons_c(d - 5):
            if sector is not None and csub(chi(m), W[i]) != sector: continue
            if parity is not None and (ypar(m) + (1 if i >= 3 else 0)) % 2 != parity: continue
            th.append((i, m))
    return th

def relaxation_quadrics(d, p, sector=None, verbose=True, parity=None):
    """Quadrics q_s(theta) = coordinate s of NF( A * sum_{a,b} theta_a theta_b m_a m_b F_{i_a i_b} ) in R_{2d}, mod p
    (p=None: exact over Q).  Returns (theta list, quadrics as dicts {(a,b): coeff} with a<=b, info)."""
    th = theta_vars(d, sector, parity)
    n = len(th)
    pairs = {}
    for a in range(n):
        ia, ma = th[a]
        for b in range(a, n):
            ib, mb = th[b]
            P = pmul(AZ, pmul(mono(ma), pmul(mono(mb), DDF[ia][ib])))
            if not P: continue
            c = chi(next(iter(P)))
            pairs.setdefault(c, []).append(((a, b), P))
    quads = {}
    info = {}; rankJ = 0
    for c, lst in pairs.items():
        red = Reducer(2 * d, c, p)
        rankJ += red.rank
        nfs = red.nf([P for _, P in lst])
        for ((a, b), _), nf in zip(lst, nfs):
            f = 1 if a == b else 2
            for s, v in nf.items():
                key = (c, s)
                quads.setdefault(key, {})[(a, b)] = (v * f) % p if p is not None else v * f
        info[str(c)] = {'pairs': len(lst), 'dimR': red.dimR, 'rankJ': red.rank}
    # sectors of S_{2d} not hit by any pair still count toward the c.i. check
    allJ = 0
    for c in SECTORS:
        if str(c) in info: allJ += info[str(c)]['rankJ']
        else:
            allJ += rank_rows(jac_gens(2 * d, c), {m: i for i, m in enumerate(mons_sector(2 * d, c))}, p)
    qraw = [q for q in quads.values() if q]
    qlist = dedup_quadrics(qraw, n, p)
    ci = dimS(2 * d) - ci_hilbert(2 * d)
    if verbose: print(f'  d={d} sector={sector}: {n} unknowns, {len(qraw)} coordinate quadrics spanning {len(qlist)}; rank J_{2*d} = {allJ} (c.i. {ci})', flush=True)
    return th, qlist, {'n_theta': n, 'n_quadrics_raw': len(qraw), 'n_quadrics': len(qlist), 'rankJ_2d': allJ, 'ci_dimJ_2d': ci, 'sectors': info}

def dedup_quadrics(qraw, n, p):
    """Basis of the span of the quadrics (same ideal, fewer generators)."""
    if not qraw: return []
    pairs = [(a, b) for a in range(n) for b in range(a, n)]
    pidx = {pr: k for k, pr in enumerate(pairs)}
    C = len(pairs)
    flat = [0] * (len(qraw) * C)
    for r, q in enumerate(qraw):
        for pr, v in q.items(): flat[r * C + pidx[pr]] = (v % p) if p is not None else v
    if p is None:
        M, rk = fmpq_mat(len(qraw), C, flat).rref()
        rows = [{pairs[j]: M[i, j] for j in range(C) if M[i, j] != 0} for i in range(rk)]
    else:
        M, rk = nmod_mat(len(qraw), C, flat, p).rref()
        rows = [{pairs[j]: int(M[i, j]) for j in range(C) if int(M[i, j]) != 0} for i in range(rk)]
    return rows

def gauge_slice(th, which='X0'):
    """Indices to remove: component j, monomials divisible by Z_j (complement of the Euler line {h E})."""
    j = {'X0': 0, 'X1': 1, 'X2': 2, 'Y0': 3, 'Y1': 4, 'Y2': 5}[which]
    return {k for k, (i, m) in enumerate(th) if i == j and m[j] >= 1}

def restrict(quads, n, removed):
    keep = [k for k in range(n) if k not in removed]
    kidx = {k: i for i, k in enumerate(keep)}
    out = []
    for q in quads:
        q2 = {(kidx[a], kidx[b]): v for (a, b), v in q.items() if a not in removed and b not in removed}
        if q2: out.append(q2)
    return keep, out

def quad_str(q, names, p=None):
    terms = []
    for (a, b), v in q.items():
        if p is not None: v = v % p
        if v == 0: continue
        mon = f'{names[a]}^2' if a == b else f'{names[a]}*{names[b]}'
        terms.append(f'{v}*{mon}' if p is not None else f'({v})*{mon}')
    return '+'.join(terms) if terms else '0'

def write_msolve(path, names, p, polys):
    with open(path, 'w') as f:
        f.write(','.join(names) + '\n' + str(p) + '\n' + ',\n'.join(polys) + '\n')

def msolve_lead(path_in, path_out, threads=2, timeout=None):
    t0 = time.time()
    r = subprocess.run(['msolve', '-g', '1', '-t', str(threads), '-f', path_in, '-o', path_out],
                       capture_output=True, text=True, timeout=timeout)
    txt = open(path_out).read() if os.path.exists(path_out) else ''
    body = txt[txt.find('['):] if '[' in txt else ''
    mons_ = re.findall(r'u\d+(?:\^\d+)?(?:\*u\d+(?:\^\d+)?)*', body)
    return r.returncode, time.time() - t0, mons_

def pure_powers(lead_monomials, names):
    pure = set()
    for mon in lead_monomials:
        vs = set(re.findall(r'[A-Za-z_]\w*', mon))
        if len(vs) == 1: pure |= vs
    return pure, all(nm in pure for nm in names)

def write_singular(path, names, p, polys, extra='print(dim(I)); print(vdim(I));'):
    with open(path, 'w') as f:
        f.write(f'ring r = {p},({",".join(names)}),dp;\n')
        f.write('ideal I = ' + ',\n'.join(polys) + ';\n')
        f.write('option(redSB); I = std(I);\n' + extra + '\nquit;\n')

def run_singular(path, timeout=None):
    t0 = time.time()
    r = subprocess.run(['Singular', '-q', path], capture_output=True, text=True, timeout=timeout)
    return r.returncode, time.time() - t0, r.stdout.strip().splitlines(), r.stderr[-300:]
