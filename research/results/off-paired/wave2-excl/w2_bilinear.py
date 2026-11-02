"""Bilinear E1+E2 systems for denominators free in a character sector (or unrestricted), chart decomposition, msolve.

D = sum_a t_a D_a  (D_a: basis of the colon sector (J:A)_{d,c} modulo the Euler gauge F*S_{d-6}, or of the
non-Jacobian part only: --tsub nonjac), N = sum_a t_a N_a + sum_k s_k S_k  (N_a particular lifts, S_k Koszul kernel).
Normalized certificate <=> B_D(N) = D div N - N(D) == 0 (mod F)   [Euler normalization N(F) = A D is lossless].
Modes:
  sector  : full bilinear system in (t, s): equations = coefficients of B_D(N) modulo F*S_{2d-6} (lossless).
  relaxJ  : necessary pure-t quadrics: B_D(S_k) lies in J(F), so B_{D(t)}(N_{D(t)}) == 0 modulo J(F)_{2d}.
Charts: 'earliest nonzero coordinate' t_b = 0 (b < a), t_a = 1 -- a complete cover of {t != 0} (system is homogeneous
in (t, s) jointly).  msolve [-1] on every chart = no solution with D != 0 modulo p (SCREEN; rigor needs certificates).
usage: python3 w2_bilinear.py --d 8 --chars 00,20,40 [--mode sector|relaxJ] [--tsub all|nonjac] [--prime P]
       [--charts all|random|a,b,c] [--threads T] [--timeout S] --job NAME [--no-msolve] [--addJ]
"""
import sys, os, json, time, argparse, random, subprocess
sys.path.insert(0, '/tmp/hodge/results/off-paired/wave2-excl')
import w2lib
from w2lib import char, char_of, cadd, cneg, VARCHAR, monos_char, N6
from hodge import core, jacobian
from hodge.core import NV, monomials, mul, add, sub, scale, pderiv, shift, reduce_mod
from flint import nmod_mat

PK = w2lib.PK
PM = 2147483647   # msolve prime (< 2^31); wave 1 used 1073741827

def vec_of(P, idx):
    v = [0] * len(idx)
    for e, val in P.items(): v[idx[e]] = val
    return v

def nullspace_rows(M, p):
    X, n = M.nullspace()
    if n == 0: return []
    C = M.ncols()
    T = nmod_mat(n, C, [int(X[j, k]) for k in range(n) for j in range(C)], p)
    Tr, rk = T.rref()
    return [[int(Tr[i, j]) for j in range(C)] for i in range(rk)]

class Reducer:
    """Normal forms modulo the row space of a list of polynomials (rref'd) within a fixed monomial list."""
    def __init__(self, polys, mons, p):
        self.p = p; self.mons = mons; self.idx = {m: i for i, m in enumerate(mons)}
        R, C = len(polys), len(mons)
        if R == 0:
            self.pivots = []; self.rows = []; self.rank = 0
        else:
            flat = [0] * (R * C)
            for r, P in enumerate(polys):
                for e, v in P.items(): flat[r * C + self.idx[e]] = v
            M = nmod_mat(R, C, flat, p)
            Mr, rank = M.rref()
            self.rank = rank
            piv = []; rows = []
            for i in range(rank):
                j = 0
                while j < C and int(Mr[i, j]) == 0: j += 1
                piv.append(j)
                rows.append({mons[k]: int(Mr[i, k]) for k in range(j, C) if int(Mr[i, k])})
            self.pivots = piv; self.rows = rows
        ps = set(self.pivots)
        self.standard = [j for j in range(len(mons)) if j not in ps]
        self.std_index = {j: i for i, j in enumerate(self.standard)}
    def normal_form(self, P):
        P = dict(P); p = self.p
        for piv, row in zip(self.pivots, self.rows):
            m = self.mons[piv]
            c = P.get(m, 0)
            if c:
                for e, v in row.items():
                    nv = (P.get(e, 0) - c * v) % p
                    if nv: P[e] = nv
                    else: P.pop(e, None)
        return P
    def coords(self, P):
        """Coordinates of the normal form on the standard monomials (dict std-index -> value)."""
        nf = self.normal_form(P)
        return {self.std_index[self.idx[e]]: v for e, v in nf.items()}

def colon_and_lifts(F, A, q, d, c, p, tsub='all'):
    """Basis of the colon sector modulo gauge (or the non-Jacobian part) with lifts, and the Koszul kernel.
    c = None: unrestricted (all monomials).  Returns (Ds, Ns, syz, info)."""
    Fp = reduce_mod(F, p); Ap = reduce_mod(A, p)
    cA = char_of(Ap)
    degN = d + 6 * q - 5
    Fi = [pderiv(Fp, i, p) for i in range(NV)]
    Dmons = monos_char(d, c) if c is not None else monomials(d)
    kinds = []; cols = []
    for m in Dmons:
        cols.append(scale(shift(Ap, m, 1, p), p - 1, p)); kinds.append(('D', -1, m))
    for i in range(NV):
        ms = monos_char(degN, cadd(cadd(c, VARCHAR[i]), cA)) if c is not None else monomials(degN)
        for m in ms:
            cols.append(shift(Fi[i], m, q, p)); kinds.append(('N', i, m))
    rows = {}
    for col in cols:
        for e in col: rows.setdefault(e, len(rows))
    R, C = len(rows), len(cols)
    flat = [0] * (R * C)
    for j, col in enumerate(cols):
        for e, v in col.items(): flat[rows[e] * C + j] = v
    M = nmod_mat(R, C, flat, p)
    vecs = nullspace_rows(M, p)
    nD = len(Dmons)
    part = []; syz = []
    for v in vecs:
        D = {}; N = [{} for _ in range(NV)]
        for j, val in enumerate(v):
            if not val: continue
            kind, comp, m = kinds[j]
            if kind == 'D': D[m] = val
            else: N[comp][m] = val
        (part if D else syz).append((D, N))
    syz = [N for D, N in syz]
    info = {'colon_dim': len(part), 'syz_dim': len(syz), 'nD': nD, 'e1_rows': R, 'e1_cols': C}
    idx = {m: i for i, m in enumerate(Dmons)}
    if tsub == 'all':
        gauge = [mul(Fp, {m: 1}, p) for m in (monos_char(d - 6, c) if c is not None else monomials(d - 6))] if d >= 6 else []
        base = [vec_of(g, idx) for g in gauge]
        cur = nmod_mat(len(base), nD, [x for r in base for x in r], p).rank() if base else 0
        chosen = []
        for (D, N) in part:
            rows_ = base + [vec_of(D_, idx) for (D_, N_) in chosen] + [vec_of(D, idx)]
            rk = nmod_mat(len(rows_), nD, [x for r in rows_ for x in r], p).rank()
            if rk > cur + len(chosen):
                chosen.append((D, N))
        info['gauge_dim'] = cur
    else:
        # non-Jacobian part: colon elements independent modulo J(F)_d (sector); representatives = normal forms mod J
        Jgens = []
        for i in range(NV):
            for m in (monos_char(d - 5, cadd(c, VARCHAR[i])) if c is not None else monomials(d - 5)):
                Jgens.append(shift(Fi[i], m, 1, p))
        red = Reducer(Jgens, Dmons, p)
        info['dimJ_sector'] = red.rank
        chosen = []
        nfs = []
        for (D, N) in part:
            nf = red.normal_form(D)
            if not nf: continue
            rows_ = nfs + [vec_of(nf, idx)]
            rk = nmod_mat(len(rows_), nD, [x for r in rows_ for x in r], p).rank()
            if rk > len(nfs):
                nfs.append(vec_of(nf, idx))
                # lift of nf: N - lift of (D - nf) where D - nf = sum m_i F_i in J  -> N_J = A * (m_i) e_i; we recompute the
                # J-part by solving: easier to keep the original lift N and use D (same class); but we want D_a = nf itself.
                chosen.append((nf, None, D, N))
        # lifts for nf: solve N(F) = A nf via the same linear system restricted: reuse nullspace trick with rhs
        out = []
        for (nf, _, D, N) in chosen:
            diff = sub(D, nf, p)           # in J_d
            # express diff = sum_i g_i F_i (g_i of degree d-5) by linear algebra, then lift N_J = A g e_i
            gcols = []; gk = []
            for i in range(NV):
                for m in (monos_char(d - 5, cadd(c, VARCHAR[i])) if c is not None else monomials(d - 5)):
                    gcols.append(shift(Fi[i], m, 1, p)); gk.append((i, m))
            rws = {}
            for col in gcols:
                for e in col: rws.setdefault(e, len(rws))
            for e in diff: rws.setdefault(e, len(rws))
            Rg, Cg = len(rws), len(gcols)
            fl = [0] * (Rg * (Cg + 1))
            for j, col in enumerate(gcols):
                for e, v in col.items(): fl[rws[e] * (Cg + 1) + j] = v
            for e, v in diff.items(): fl[rws[e] * (Cg + 1) + Cg] = v
            Mg = nmod_mat(Rg, Cg + 1, fl, p)
            Mr, rk = Mg.rref()
            piv = []
            for i in range(rk):
                j = 0
                while j <= Cg and int(Mr[i, j]) == 0: j += 1
                piv.append(j)
            assert piv[-1] != Cg, 'diff not in J?'
            g = [{} for _ in range(NV)]
            for i, j in enumerate(piv):
                v = int(Mr[i, Cg])
                if v: g[gk[j][0]][gk[j][1]] = v
            NJ = [mul(Ap, g[i], p) for i in range(NV)]
            Nnf = [sub(N[i], NJ[i], p) for i in range(NV)]
            out.append((nf, Nnf))
        chosen = out
        info['gauge_dim'] = None
    Ds = [D for D, N in chosen]; Ns = [N for D, N in chosen]
    # sanity: lifts satisfy N(F) = A D
    for D, N in zip(Ds, Ns):
        lhs = {}
        for i in range(NV): lhs = add(lhs, mul(N[i], Fi[i], p), p)
        assert lhs == mul(Ap, D, p), 'lift check failed'
    info['nt'] = len(Ds); info['ns'] = len(syz)
    return Ds, Ns, syz, info

def B_poly(D, N, p):
    divN = {}; ND = {}
    for i in range(NV):
        if N[i]:
            divN = add(divN, pderiv(N[i], i, p), p)
            ND = add(ND, mul(N[i], pderiv(D, i, p), p), p)
    return sub(mul(D, divN, p), ND, p)

def build(F, A, q, d, c, p, mode='sector', tsub='all', addJ=True, verbose=True):
    t0 = time.time()
    Ds, Ns, syz, info = colon_and_lifts(F, A, q, d, c, p, tsub)
    nt, ns = len(Ds), len(syz)
    if nt == 0:
        return None, info
    Fp = reduce_mod(F, p); Ap = reduce_mod(A, p); cA = char_of(Ap)
    sector = cadd(cA, cadd(c, c)) if c is not None else None
    mons2d = monos_char(2 * d, sector) if c is not None else monomials(2 * d)
    Fi = [pderiv(Fp, i, p) for i in range(NV)]
    eqs = {}   # equation index -> {key: coeff}
    def addeq(coords, key):
        for r, v in coords.items():
            eqs.setdefault(r, {})
            eqs[r][key] = (eqs[r].get(key, 0) + v) % p
    if mode == 'sector':
        Fm = [shift(Fp, m, 1, p) for m in (monos_char(2 * d - 6, sector) if c is not None else monomials(2 * d - 6))] if 2 * d >= 6 else []
        red = Reducer(Fm, mons2d, p)
        for a in range(nt):
            for b in range(nt):
                key = ('t', min(a, b), 't', max(a, b))
                addeq(red.coords(B_poly(Ds[a], Ns[b], p)), key)
            for k in range(ns):
                addeq(red.coords(B_poly(Ds[a], syz[k], p)), ('t', a, 's', k))
        neqF = len(red.standard)
        info.update({'mode': 'sector', 'neq': neqF, 'red_rank': red.rank})
        if addJ:
            # Jacobian relaxation quadrics as extra (consequence) equations, indexed by ('J', r)
            Jg = [shift(Fi[i], m, 1, p) for i in range(NV) for m in (monos_char(2 * d - 5, cadd(sector, VARCHAR[i])) if c is not None else monomials(2 * d - 5))]
            redJ = Reducer(Jg, mons2d, p)
            for a in range(nt):
                for b in range(nt):
                    key = ('t', min(a, b), 't', max(a, b))
                    cs = redJ.coords(B_poly(Ds[a], Ns[b], p))
                    for r, v in cs.items():
                        eqs.setdefault(('J', r), {})
                        eqs[('J', r)][key] = (eqs[('J', r)].get(key, 0) + v) % p
            info['neqJ'] = len(redJ.standard)
    else:  # relaxJ: pure-t quadrics modulo J(F)_{2d}
        Jg = [shift(Fi[i], m, 1, p) for i in range(NV) for m in (monos_char(2 * d - 5, cadd(sector, VARCHAR[i])) if c is not None else monomials(2 * d - 5))]
        redJ = Reducer(Jg, mons2d, p)
        for a in range(nt):
            for b in range(nt):
                key = ('t', min(a, b), 't', max(a, b))
                addeq(redJ.coords(B_poly(Ds[a], Ns[b], p)), key)
        info.update({'mode': 'relaxJ', 'neq': len(redJ.standard), 'ring_dim_2d': len(redJ.standard)})
        ns = 0
    eqs = {r: e for r, e in eqs.items() if any(v for v in e.values())}
    info.update({'nt': nt, 'ns': ns, 'neq_nonzero': len(eqs), 'build_s': round(time.time() - t0, 1), 'prime': p})
    if verbose: print(f'[build d={d} c={c} mode={mode} tsub={tsub}] {info}', file=sys.stderr, flush=True)
    return {'eqs': list(eqs.values()), 'nt': nt, 'ns': ns, 'Ds': Ds, 'Ns': Ns, 'syz': syz, 'p': p, 'd': d, 'c': c, 'mode': mode}, info

def build_fast(F, A, q, d, c, p, mode='sector', tsub='all', addJ=True, verbose=True, block=4000):
    """Same output as build(), using nmod_mpoly products and matrix normal forms (for large systems)."""
    import w2_fast
    t0 = time.time()
    Ds, Ns, syz, info = colon_and_lifts(F, A, q, d, c, p, tsub)
    nt, ns = len(Ds), len(syz)
    if nt == 0:
        return None, info
    Fp = reduce_mod(F, p); Ap = reduce_mod(A, p); cA = char_of(Ap)
    sector = cadd(cA, cadd(c, c)) if c is not None else None
    mons2d = monos_char(2 * d, sector) if c is not None else monomials(2 * d)
    Fi = [pderiv(Fp, i, p) for i in range(NV)]
    ctx = w2_fast.ctx_for(p)
    if mode == 'relaxJ':
        syz = []; ns = 0
    prods = w2_fast.B_products(ctx, Ds, Ns, syz, p, verbose=verbose)
    eqs = {}
    def addeq(coords, key, tag=None):
        for r, v in coords.items():
            rr = r if tag is None else (tag, r)
            eqs.setdefault(rr, {})
            eqs[rr][key] = (eqs[rr].get(key, 0) + v) % p
    if mode == 'sector':
        Fm = [shift(Fp, m, 1, p) for m in (monos_char(2 * d - 6, sector) if c is not None else monomials(2 * d - 6))] if 2 * d >= 6 else []
        red = w2_fast.MatReducer(Fm, mons2d, p, verbose=verbose)
        for s0 in range(0, len(prods), block):
            blk = prods[s0:s0 + block]
            for (key, _), cs in zip(blk, red.coords_block([P for _, P in blk])):
                addeq(cs, key)
        info.update({'mode': 'sector', 'neq': red.nstd, 'red_rank': red.rank})
    if mode == 'relaxJ' or addJ:
        Jg = [shift(Fi[i], m, 1, p) for i in range(NV) for m in (monos_char(2 * d - 5, cadd(sector, VARCHAR[i])) if c is not None else monomials(2 * d - 5))]
        redJ = w2_fast.MatReducer(Jg, mons2d, p, verbose=verbose)
        tt = [(key, P) for key, P in prods if key[2] == 't']
        for s0 in range(0, len(tt), block):
            blk = tt[s0:s0 + block]
            for (key, _), cs in zip(blk, redJ.coords_block([P for _, P in blk])):
                addeq(cs, key, tag='J')
        if mode == 'relaxJ':
            info.update({'mode': 'relaxJ', 'neq': redJ.nstd, 'ring_dim_2d': redJ.nstd})
        else:
            info['neqJ'] = redJ.nstd
    eqs = {r: e for r, e in eqs.items() if any(v for v in e.values())}
    info.update({'nt': nt, 'ns': ns, 'neq_nonzero': len(eqs), 'build_s': round(time.time() - t0, 1), 'prime': p, 'fast': True})
    if verbose: print(f'[build_fast d={d} c={c} mode={mode} tsub={tsub}] {info}', file=sys.stderr, flush=True)
    return {'eqs': list(eqs.values()), 'nt': nt, 'ns': ns, 'Ds': Ds, 'Ns': Ns, 'syz': syz, 'p': p, 'd': d, 'c': c, 'mode': mode}, info

def chart_polys(sysd, a):
    """Substitute t_b = 0 (b < a), t_a = 1.  Returns (names, list of polys as dict monomial-tuple -> coeff) over the
    remaining variables t_{a+1..}, s_*; monomials encoded as sorted tuples of variable names."""
    nt, ns, p = sysd['nt'], sysd['ns'], sysd['p']
    polys = []
    for eq in sysd['eqs']:
        P = {}
        for key, cf in eq.items():
            k1, i1, k2, i2 = key
            # variables
            vs = []
            dead = False
            for kk, ii in ((k1, i1), (k2, i2)):
                if kk == 't':
                    if ii < a: dead = True; break
                    if ii == a: continue
                    vs.append(f't{ii}')
                else:
                    vs.append(f's{ii}')
            if dead: continue
            mon = tuple(sorted(vs))
            P[mon] = (P.get(mon, 0) + cf) % p
        P = {m: v for m, v in P.items() if v}
        if P: polys.append(P)
    names = [f't{i}' for i in range(a + 1, nt)] + [f's{j}' for j in range(ns)]
    return names, polys

def write_ms(names, polys, p, path):
    used = set()
    for P in polys:
        for m in P: used.update(m)
    names = [n for n in names if n in used] or names[:1]
    lines = [','.join(names), str(p)]
    strs = []
    for P in polys:
        terms = []
        for m, v in P.items():
            terms.append(str(v) + (('*' + '*'.join(m)) if m else ''))
        strs.append('+'.join(terms))
    lines.append(',\n'.join(strs))
    with open(path, 'w') as f: f.write('\n'.join(lines) + '\n')
    return len(names), len(polys)

def run_msolve(path, out, threads=4, timeout=3600):
    t0 = time.time()
    try:
        r = subprocess.run(['msolve', '-f', path, '-o', out, '-t', str(threads), '-v', '0'], capture_output=True, text=True, timeout=timeout)
        res = open(out).read().strip() if os.path.exists(out) else ''
    except subprocess.TimeoutExpired:
        return 'timeout', round(time.time() - t0, 1), ''
    dt = round(time.time() - t0, 1)
    if res.startswith('[-1]'): return 'no-solution', dt, res[:200]
    if res.startswith('[1,'): return 'positive-dimensional', dt, res[:2000]
    if res.startswith('[0,'): return 'zero-dimensional', dt, res[:2000]
    return 'other:' + (r.stderr[-300:] if r.stderr else ''), dt, res[:2000]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--d', type=int, required=True); ap.add_argument('--chars', default='all')
    ap.add_argument('--mode', default='sector'); ap.add_argument('--tsub', default='all')
    ap.add_argument('--prime', type=int, default=PM); ap.add_argument('--charts', default='all')
    ap.add_argument('--threads', type=int, default=4); ap.add_argument('--timeout', type=int, default=3600)
    ap.add_argument('--job', required=True); ap.add_argument('--no-msolve', action='store_true'); ap.add_argument('--addJ', action='store_true')
    ap.add_argument('--outdir', default=None); ap.add_argument('--fast', action='store_true')
    a = ap.parse_args()
    F, A, q = w2lib.load()
    if a.chars == 'all':
        chars = [(x, y) for x in range(6) for y in range(6)]
    elif a.chars == 'mixed':
        chars = [None]
    else:
        chars = [(int(s[0]), int(s[1])) for s in a.chars.split(',')]
    outdir = a.outdir or f'{PK}/bilinear_{a.job}'
    os.makedirs(outdir, exist_ok=True)
    recfile = f'{PK}/receipts_{a.job}.jsonl'
    for c in chars:
        sysd, info = (build_fast if a.fast else build)(F, A, q, a.d, c, a.prime, a.mode, a.tsub, addJ=a.addJ)
        tag = f'd{a.d}_c{"mixed" if c is None else f"{c[0]}{c[1]}"}_{a.mode}_{a.tsub}'
        base = {'exp': 'bilinear', 'job': a.job, 'd': a.d, 'char': c, 'mode': a.mode, 'tsub': a.tsub, 'prime': a.prime, 'family': f'all D in sector ({a.tsub})', **info}
        if sysd is None:
            rec = dict(base, result='empty-t-space', seconds=0); w2lib.log(rec)
            with open(recfile, 'a') as f: f.write(json.dumps(rec, default=str) + '\n')
            print(json.dumps(rec, default=str), flush=True); continue
        nt = sysd['nt']
        if a.charts == 'all': charts = list(range(nt))
        elif a.charts == 'rev': charts = list(range(nt - 1, -1, -1))
        elif a.charts == 'random': charts = ['random']
        else: charts = [int(x) for x in a.charts.split(',')]
        for ch in charts:
            t0 = time.time()
            if ch == 'random':
                # random affine normalization sum r_a t_a = 1 (screen only)
                rng = random.Random(1)
                names = [f't{i}' for i in range(nt)] + [f's{j}' for j in range(sysd['ns'])]
                polys = []
                for eq in sysd['eqs']:
                    P = {}
                    for (k1, i1, k2, i2), cf in eq.items():
                        mon = tuple(sorted([f'{k1}{i1}', f'{k2}{i2}']))
                        P[mon] = (P.get(mon, 0) + cf) % a.prime
                    polys.append({m: v for m, v in P.items() if v})
                polys.append({**{(f't{i}',): rng.randrange(1, a.prime) for i in range(nt)}, (): a.prime - 1})
            else:
                names, polys = chart_polys(sysd, ch)
            path = f'{outdir}/{tag}_chart{ch}.ms'
            nv, npol = write_ms(names, polys, a.prime, path)
            rec = dict(base, chart=ch, nvars=nv, npolys=npol)
            if a.no_msolve:
                rec.update(result='written', seconds=round(time.time() - t0, 1))
            else:
                st, dt, res = run_msolve(path, path.replace('.ms', '.out'), a.threads, a.timeout)
                rec.update(msolve=st, result=('no-solution(modp)' if st == 'no-solution' else st), msolve_out=res[:300], seconds=dt)
                if st not in ('no-solution', 'timeout'):
                    print('!!! POSSIBLE SOLUTIONS on chart', ch, res[:500], flush=True)
            w2lib.log(rec)
            with open(recfile, 'a') as f: f.write(json.dumps(rec, default=str) + '\n')
            print(json.dumps({k: rec[k] for k in ('d', 'char', 'mode', 'tsub', 'nt', 'ns', 'chart', 'nvars', 'npolys', 'result', 'seconds')}, default=str), flush=True)

if __name__ == '__main__':
    main()
