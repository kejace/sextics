"""Strict five-variable (Y2-free) reduction, exact.

Lemma (proved in the wave-2 report): for a Y2-free homogeneous D != 0 of degree d, a certificate (D, N, Q) for
(F = f + Y2^6, A, q = 1) exists  <=>  there is N' in Q[X0,X1,X2,Y0,Y1]^5 of degree d+1 with
        E1':  N'(f) = A D          E2':  D div N' - N'(D) = 0        (no Q, no F; everything Y2-free).
Here F7 = 7F = f7 + 7 Y2^6, f7 integer; scaling f by a constant is harmless.
This module runs the same relaxation pipeline as w2x in that setting: J_f-duals (mod J_f Hessian relaxation) and the
strict relaxed duals (left kernel of {B_{D_a}(K)} over 5-variable Koszul corrections K, no F columns)."""
import sys, time, math, itertools
sys.path.insert(0, '/tmp/hodge/results/one-reflection-skew/wave2-excl')
import w2x
from w2x import ctx, Z, NV, chr_of, chr_add, chr_neg, CHI_A, W, monom, vec_of, mat_from_cols, poly_deg, poly_chr, exact_solve, BD
from flint import fmpz_mat, nmod_mat, fmpq

f7 = w2x.F7 - 7 * Z[5]**6
assert all(e[5] == 0 for e in f7.to_dict())
fi = [f7.derivative(i) for i in range(5)]
CHI_f = poly_chr(f7)
assert CHI_f == w2x.CHI_F

_m5 = {}
def mons5(n):
    if n not in _m5:
        _m5[n] = [e for e in w2x.mons(n) if e[5] == 0]
    return _m5[n]
_s5 = {}
def secmons5(n, chi):
    key = (n, chi)
    if key not in _s5:
        _s5[key] = [e for e in mons5(n) if chr_of(e) == chi]
    return _s5[key]
_si5 = {}
def secidx5(n, chi):
    key = (n, chi)
    if key not in _si5:
        _si5[key] = {e: k for k, e in enumerate(secmons5(n, chi))}
    return _si5[key]
def vec5(P, n, chi):
    idx = secidx5(n, chi)
    return {idx[e]: int(c) for e, c in P.to_dict().items()}

def jgens5(n, chi):
    out = []
    if n < 5: return out
    for i in range(5):
        ci = poly_chr(fi[i]); need = chr_add(chi, chr_neg(ci))
        for m in secmons5(n - 5, need):
            out.append((i, m, monom(m) * fi[i]))
    return out

_jd5 = {}
def jac_dual5(n, chi):
    key = (n, chi)
    if key in _jd5: return _jd5[key]
    rows = secmons5(n, chi)
    gens = jgens5(n, chi)
    M = mat_from_cols([vec5(g[2], n, chi) for g in gens], len(rows))
    if M.ncols() == 0:
        Y = fmpz_mat(len(rows), len(rows))
        for i in range(len(rows)): Y[i, i] = 1
        _jd5[key] = (Y, 0); return _jd5[key]
    Yfull, k = M.transpose().nullspace()
    Y = fmpz_mat(len(rows), k)
    for i in range(len(rows)):
        for j in range(k): Y[i, j] = Yfull[i, j]
    _jd5[key] = (Y, len(rows) - k)
    return _jd5[key]

def e1_kernel5(d, chi):
    Dm = secmons5(d, chi)
    tchi = chr_add(chi, CHI_A)
    Y, rankJ = jac_dual5(d + 6, tchi)
    cols = [vec5(w2x.A * monom(m), d + 6, tchi) for m in Dm]
    Mc = mat_from_cols(cols, len(secmons5(d + 6, tchi)))
    C = Y.transpose() * Mc
    if C.nrows() == 0:
        basis = [monom(m) for m in Dm]
    else:
        K, kd = C.nullspace()
        basis = []
        for j in range(kd):
            dct = {Dm[i]: int(K[i, j]) for i in range(len(Dm)) if int(K[i, j])}
            g = 0
            for v in dct.values(): g = math.gcd(g, abs(v))
            basis.append(ctx.from_dict({e: v // g for e, v in dct.items()}))
    rJ = jac_dual5(d, chi)[1] if d >= 5 else 0
    return {'d': d, 'chi': chi, 'nD': len(Dm), 'dimK': len(basis), 'dimJ': rJ, 'extra': len(basis) - rJ, 'basis': basis}

def e1_lift5(D):
    d = poly_deg(D); chi = poly_chr(D); tchi = chr_add(chi, CHI_A); n = d + 6
    gens = jgens5(n, tchi); rows = secmons5(n, tchi)
    M = mat_from_cols([vec5(g[2], n, tchi) for g in gens], len(rows))
    b = fmpz_mat(len(rows), 1)
    for r, v in vec5(w2x.A * D, n, tchi).items(): b[r, 0] = v
    x = exact_solve(M, b)
    if x is None: return None
    den = 1
    for v in x:
        dd = int(fmpq(v).q); den = den * dd // math.gcd(den, dd)
    N = [ctx.from_dict({}) for _ in range(NV)]
    for (i, m, _), v in zip(gens, x):
        if v != 0: N[i] += monom(m, int(v * den))
    assert sum((N[i] * fi[i] for i in range(5)), ctx.from_dict({})) == w2x.A * D * den
    return N, den

def in_jf_exact(P):
    n = poly_deg(P); chi = poly_chr(P)
    if n < 0: return True
    Y, _ = jac_dual5(n, chi); v = vec5(P, n, chi)
    for j in range(Y.ncols()):
        if sum(int(Y[r, j]) * c for r, c in v.items()): return False
    return True

class Family5(w2x.Family):
    @staticmethod
    def from_e1(d, sectors=None, verbose=False):
        items = []
        for chi in (sectors or w2x.ALLCHR):
            r = e1_kernel5(d, chi)
            for k, D in enumerate(r['basis']):
                N, den = e1_lift5(D)
                items.append((f'5v:d{d}:{chi}:{k}', D, N, den, chi, 'J' if in_jf_exact(D) else 'extra'))
            if verbose and r['dimK']: print('5v sector', chi, 'dimK', r['dimK'], 'dimJ', r['dimJ'], 'extra', r['extra'], file=sys.stderr, flush=True)
        return Family5(d, items)

def koszul5(d, chiD):
    out = []
    if d - 4 < 0: return out
    for i in range(5):
        for j in range(i + 1, 5):
            need = chr_add(CHI_A, chiD, W[i], W[j])
            for m in secmons5(d - 4, need):
                out.append((i, j, m))
    return out

def BD_koszul5(D, i, j, m):
    c = monom(m)
    return fi[j] * (D * c.derivative(i) - c * D.derivative(i)) - fi[i] * (D * c.derivative(j) - c * D.derivative(j))

def jac_duals_for5(fam):
    return {psi: jac_dual5(2 * fam.d, psi)[0] for psi in w2x.ALLCHR}

def quadrics_from_duals5(fam, duals, which_pairs=None):
    """Same as w2x.quadrics_from_duals but rows are the Y2-free monomials of S_{2d}[psi]."""
    n2 = 2 * fam.d
    bysec = {}
    for a, b in (which_pairs or fam.pairs()):
        psi = chr_add(fam.chi[a], fam.chi[b], CHI_A); bysec.setdefault(psi, []).append((a, b))
    quads = []
    for psi, prs in bysec.items():
        if psi not in duals: continue
        Y = duals[psi]; k = Y.ncols()
        if k == 0: continue
        nrows = len(secmons5(n2, psi))
        cols = []; scales = []
        for a, b in prs:
            P, sc = w2x.pair_B(fam, a, b)
            cols.append(vec5(P, n2, psi) if not P.is_zero() else {}); scales.append(sc)
        Qm = Y.transpose() * mat_from_cols(cols, nrows)
        L = 1
        for sc in scales: L = L * sc // math.gcd(L, sc)
        for i in range(k):
            qd = {}
            for j, (a, b) in enumerate(prs):
                v = int(Qm[i, j])
                if v: qd[(a, b)] = v * (L // scales[j])
            if qd:
                g = 0
                for v in qd.values(): g = math.gcd(g, abs(v))
                quads.append({'psi': psi, 'coef': {kk: v // g for kk, v in qd.items()}})
    return quads

def relaxed_duals5(fam, subset=None, bank=None, verbose=False):
    """Exact left kernel of the strict relaxed column space {B_{D_a}(K)} (5-var Koszul K of all characters), per row sector."""
    d = fam.d; n2 = 2 * d
    subset = list(range(fam.n)) if subset is None else sorted(subset)
    need = {chr_add(fam.chi[a], fam.chi[b], CHI_A) for a in subset for b in subset if a <= b}
    cols_by = {psi: [] for psi in need}
    for a in subset:
        _, Da, _, _, chia, _ = fam.items[a]
        for psi in need:
            if bank is not None:
                cols_by[psi].extend(bank[psi].get(a, []))
                continue
            chiK = chr_add(psi, chr_neg(chia), chr_neg(CHI_A))
            for (i, j, m) in koszul5(d, chiK):
                P = BD_koszul5(Da, i, j, m)
                if P.is_zero(): continue
                cols_by[psi].append(vec5(P, n2, psi))
    out = {}
    for psi in need:
        rows = secmons5(n2, psi)
        M = mat_from_cols(cols_by[psi], len(rows))
        if M.ncols() == 0:
            Y = fmpz_mat(len(rows), len(rows))
            for i in range(len(rows)): Y[i, i] = 1
            out[psi] = Y; continue
        Yfull, k = M.transpose().nullspace()
        Y = fmpz_mat(len(rows), k)
        for i in range(len(rows)):
            for j in range(k): Y[i, j] = Yfull[i, j]
        out[psi] = Y
        if verbose: print('5v relaxed sector', psi, 'rows', len(rows), 'cols', M.ncols(), 'duals', k, file=sys.stderr, flush=True)
    return out

def build_bank5(fam):
    """{psi: {a: [sparse cols]}} of all strict correction columns."""
    d = fam.d; n2 = 2 * d
    bank = {psi: {} for psi in w2x.ALLCHR}
    for a in range(fam.n):
        _, Da, _, _, chia, _ = fam.items[a]
        for psi in w2x.ALLCHR:
            chiK = chr_add(psi, chr_neg(chia), chr_neg(CHI_A))
            cols = []
            for (i, j, m) in koszul5(d, chiK):
                P = BD_koszul5(Da, i, j, m)
                if not P.is_zero(): cols.append(vec5(P, n2, psi))
            if cols: bank[psi][a] = cols
    return bank

def fixed_D_system5(D):
    """Strict 5-var system for a Y2-free eigenvector D: unknowns N_i (i<5), E1': N(f)=AD, E2': B_D(N)=0."""
    d = poly_deg(D); chi = poly_chr(D)
    Ncols = []
    for i in range(5):
        need = chr_add(CHI_A, chi, chr_neg(CHI_f), W[i])
        for m in secmons5(d + 1, need): Ncols.append((i, m))
    e1chi = chr_add(CHI_A, chi); e2chi = chr_add(CHI_A, chi, chi)
    r1 = secmons5(d + 6, e1chi); r2 = secmons5(2 * d, e2chi); n1 = len(r1)
    Di = [D.derivative(i) for i in range(5)]
    cols = []
    for (i, m) in Ncols:
        mm = monom(m); col = {}
        for r, v in vec5(mm * fi[i], d + 6, e1chi).items(): col[r] = v
        for r, v in vec5(D * mm.derivative(i) - mm * Di[i], 2 * d, e2chi).items(): col[n1 + r] = v
        cols.append(col)
    M = mat_from_cols(cols, n1 + len(r2))
    b = fmpz_mat(n1 + len(r2), 1)
    for r, v in vec5(w2x.A * D, d + 6, e1chi).items(): b[r, 0] = v
    return M, b, {'d': d, 'chi': chi, 'nN': len(Ncols), 'rows': n1 + len(r2), 'rowsE1': n1, 'rowsE2': len(r2)}, Ncols

def test_fixed_D5(D, name, family, p=2147483647, want_exact_dual=True, verbose=False):
    t0 = time.time()
    M, b, meta, Ncols = fixed_D_system5(D)
    r, c = M.nrows(), M.ncols()
    Mp = nmod_mat(r, c + 1, p); Mp0 = nmod_mat(r, c, p)
    for i in range(r):
        for j in range(c):
            v = int(M[i, j])
            if v: Mp[i, j] = v % p; Mp0[i, j] = v % p
        v = int(b[i, 0])
        if v: Mp[i, c] = v % p
    rank = int(Mp0.rank()); rank_aug = int(Mp.rank())
    out = dict(meta); out.update({'name': name, 'family': family, 'p': p, 'rank': rank, 'rank_aug': rank_aug, 'nullity': c - rank})
    if rank_aug == rank:
        x = exact_solve(M, b, p)
        if x is not None:
            den = 1
            for v in x:
                dd = int(fmpq(v).q); den = den * dd // math.gcd(den, dd)
            N = [ctx.from_dict({}) for _ in range(NV)]
            for (i, m), v in zip(Ncols, x):
                if v != 0: N[i] += monom(m, int(v * den))
            ok = (sum((N[i] * fi[i] for i in range(5)), ctx.from_dict({})) == w2x.A * D * den) and BD(D, N).is_zero()
            # 6-variable certificate: (D, N/den, Q=0) with N_5 = 0 for F7 = f7 + 7 Y2^6
            lhs = w2x.A * D * D * den - D * sum((N[i] * w2x.Fi[i] for i in range(NV)), ctx.from_dict({})) + w2x.F7 * BD(D, N)
            ok6 = lhs.is_zero()
            out['status'] = 'WITNESS' if (ok and ok6) else 'consistent-modp-lift-failed'
            out['witness'] = {'den': den, 'D': str(D), 'N': [str(n) for n in N], 'Q': '0', 'verified5': ok, 'verified6': ok6}
        else:
            out['status'] = 'consistent-modp-lift-failed'
    else:
        out['status'] = 'inconsistent-modp'
        if want_exact_dual:
            cert = w2x.exact_dual_certificate(M, b, p)
            if cert is not None: out['status'] = 'EXCLUDED-exact'
    out['seconds'] = round(time.time() - t0, 2)
    w2x.logexp(exp='fixed-D-5var', D=name, family=family, d=meta['d'], chi=list(meta['chi']), prime=p, result=out['status'], rank=rank, rank_aug=rank_aug,
               cols=c, rows=r, seconds=out['seconds'], note=f"strict Y2-free reduction E1'+E2' (no Q, no F); nN={meta['nN']}")
    if verbose: print(name, {k: v for k, v in out.items() if k != 'witness'}, file=sys.stderr, flush=True)
    return out

# ---------------------------------------------------------------- second-level (lifted) relaxation, strict 5-var
def koszul5_all(d):
    """All 5-var Koszul basis elements for deg N = d+1: list of (i, j, mon, charN)."""
    out = []
    for i in range(5):
        for j in range(i + 1, 5):
            for mon in mons5(d - 4):
                charN = chr_add(chr_of(mon), CHI_f, chr_neg(W[i]), chr_neg(W[j]))
                out.append((i, j, mon, charN))
    return out

def build_bank5k(fam):
    """bank[a] = {kidx: (psi, sparse col)} of B_{D_a}(K_kidx) (nonzero ones), psi = chi_a + charN(K)."""
    d = fam.d; n2 = 2 * d
    Ks = koszul5_all(d)
    bank = {}
    for a in range(fam.n):
        Da = fam.items[a][1]; chia = fam.chi[a]
        cols = {}
        for kidx, (i, j, mon, charN) in enumerate(Ks):
            P = BD_koszul5(Da, i, j, mon)
            if P.is_zero(): continue
            psi = chr_add(chia, charN)
            cols[kidx] = (psi, vec5(P, n2, psi))
        bank[a] = cols
    return bank, Ks

def lift2_cubics5(fam, alive, bank, Ks, verbose=False, max_cols=60000):
    """Lifted relaxation: multiply E2' by lambda_m, relax Z_{ {a,m}, K } = lambda_a lambda_m c_K (symmetric).
    Exact left kernels per total character Psi; returns list of cubics {sorted (m,a,b): int} (den-cleared)."""
    from fractions import Fraction
    d = fam.d; n2 = 2 * d
    alive = sorted(alive)
    Psis = {chr_add(fam.chi[m], fam.chi[a], Ks[k][3]) for m in alive for a in alive for k in bank[a] if a <= m or True}
    cubics = []
    for Psi in sorted(Psis):
        # row blocks
        off = {}; nrows = 0
        for m in alive:
            psi_m = chr_add(Psi, chr_neg(fam.chi[m]))
            off[m] = (nrows, psi_m); nrows += len(secmons5(n2, psi_m))
        cols = []
        for ia, a in enumerate(alive):
            for m in alive[ia:]:
                for kidx, (psi_a_k, col_a) in bank[a].items():
                    # need chi_a + chi_m + charN(K) = Psi  <=> psi_a_k (= chi_a + charN) = Psi - chi_m
                    if psi_a_k != off[m][1]: continue
                    v = {off[m][0] + r: c for r, c in col_a.items()}
                    if a != m:
                        cm = bank[m].get(kidx)
                        if cm is not None:
                            psi_m_k, col_m = cm
                            assert psi_m_k == off[a][1]
                            for r, c in col_m.items():
                                v[off[a][0] + r] = v.get(off[a][0] + r, 0) + c
                    cols.append(v)
        if nrows == 0: continue
        if len(cols) > max_cols:
            if verbose: print('lift2: skipping Psi', Psi, 'cols', len(cols), file=sys.stderr, flush=True)
            continue
        M = mat_from_cols(cols, nrows)
        if M.ncols() == 0:
            Yk = None
        else:
            Yfull, k = M.transpose().nullspace()
        if verbose: print('lift2 Psi', Psi, 'rows', nrows, 'cols', len(cols), 'duals', (nrows if M.ncols() == 0 else k), file=sys.stderr, flush=True)
        # cubic coefficients: for each dual y and each block m: y^{(m)} . b_ab over pairs in sector psi_m
        pairs_by_psi = {}
        for ia, a in enumerate(alive):
            for b in alive[ia:]:
                pairs_by_psi.setdefault(chr_add(fam.chi[a], fam.chi[b], CHI_A), []).append((a, b))
        # precompute pair vectors per sector
        pv = {}
        for psi, prs in pairs_by_psi.items():
            for (a, b) in prs:
                P, sc = w2x.pair_B(fam, a, b)
                pv[(a, b)] = (sc, vec5(P, n2, psi) if not P.is_zero() else {})
        nd = nrows if M.ncols() == 0 else k
        for jdx in range(nd):
            cub = {}
            for m in alive:
                o, psi_m = off[m]
                for (a, b) in pairs_by_psi.get(psi_m, []):
                    sc, vec = pv[(a, b)]
                    s = 0
                    for r, c in vec.items():
                        yv = 1 if (M.ncols() == 0 and o + r == jdx) else (int(Yfull[o + r, jdx]) if M.ncols() else 0)
                        if yv: s += yv * c
                    if s:
                        key = tuple(sorted((m, a, b)))
                        cub[key] = cub.get(key, Fraction(0)) + Fraction(s, sc)
            cub = {kk: v for kk, v in cub.items() if v}
            if cub:
                L = 1
                for v in cub.values(): L = L * v.denominator // math.gcd(L, v.denominator)
                ic = {kk: int(v * L) for kk, v in cub.items()}
                g = 0
                for v in ic.values(): g = math.gcd(g, abs(v))
                cubics.append({'Psi': Psi, 'coef': {kk: v // g for kk, v in ic.items()}})
    return cubics

def polys_to_msolve(quads, cubics, alive, path, p=1073741827):
    alive = sorted(alive); names = [f'l{a}' for a in alive]
    lines = []
    for qd in w2x.restrict_quads(quads, alive):
        terms = [f'{v % p}*l{a}*l{b}' for (a, b), v in qd['coef'].items() if v % p]
        if terms: lines.append('+'.join(terms))
    aset = set(alive)
    for cb in cubics:
        terms = []
        for key, v in cb['coef'].items():
            if all(x in aset for x in key) and v % p:
                terms.append(f'{v % p}*' + '*'.join(f'l{x}' for x in key))
        if terms: lines.append('+'.join(terms))
    with open(path, 'w') as fh:
        fh.write(','.join(names) + '\n' + str(p) + '\n' + ',\n'.join(lines) + '\n')
    return len(lines), len(names)
