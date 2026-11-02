"""Part (4): structured denominator families (cones, block multiples, bihomogeneous pieces, Hessian minors, radial
multiples) at degrees 6..12, character-restricted.

  python3 -u w2_struct.py --degrees 6,7,8,9,10,11,12 --procs 8 [--restarts 200] [--tag struct]

For every structured linear subspace V of S_d (both parities, every sigma-weight k) we compute the E1-admissible
part V ∩ (J(F):A)_d mod p (zero => the whole family is excluded rigorously as supports, by semicontinuity and the
prescribed-support theorem); for nonzero intersections we lift a small-height rational basis (CRT + LLL), compute
exact E1 lifts and run the sector-restricted Koszul varpro over that family.  Explicit single candidates (Hessian
minors etc.) are tested exactly (E1, then the Koszul correction system mod p).
"""
import os, sys, time, json, math, argparse, pickle, traceback, itertools
os.environ['OMP_NUM_THREADS'] = '1'; os.environ['OPENBLAS_NUM_THREADS'] = '1'
sys.path.insert(0, '/tmp/hodge/results/klein-coupled-character')
import numpy as np
import w2lib as L
from hodge import core, fixed_d
from hodge.core import NV, monomials, shift, mul, add, sub, scale, pderiv
import kcc2
from kcc import W, wt, wmonomials, sigma_weight
from kcc2 import Ctx, PI, dense_from_cols
from flint import nmod_mat
from fractions import Fraction
from multiprocessing import Pool

LANE = L.LANE
F, A, q = L.load()
Fi = [pderiv(F, i) for i in range(NV)]
X = [core.parse(v) for v in core.VARS]
K_block = core.parse('X0**5*X1 + X1**5*X2 + X2**5*X0')
G_block = core.parse('Y0**6 + Y1**6 + Y2**6')
M_block = sub(sub(F, K_block), G_block)
EX = {}
for i in range(3):
    EX = add(EX, mul(X[i], Fi[i]))
EY = {}
for i in range(3, 6):
    EY = add(EY, mul(X[i], Fi[i]))

PREDS = {
    'pureX': lambda e: e[3] == e[4] == e[5] == 0,
    'pureY': lambda e: e[0] == e[1] == e[2] == 0,
    'vars_xyzt': lambda e: e[4] == e[5] == 0,
    'vars_xyzr': lambda e: e[3] == e[5] == 0,
    'vars_xyzs': lambda e: e[3] == e[4] == 0,
    'vars_trsx': lambda e: e[1] == e[2] == 0,
    'vars_trsy': lambda e: e[0] == e[2] == 0,
    'vars_trsz': lambda e: e[0] == e[1] == 0,
    'no_x': lambda e: e[0] == 0, 'no_y': lambda e: e[1] == 0, 'no_z': lambda e: e[2] == 0,
    'no_t': lambda e: e[3] == 0, 'no_r': lambda e: e[4] == 0, 'no_s': lambda e: e[5] == 0,
    'Xdeg_le3': lambda e: e[0] + e[1] + e[2] <= 3,
    'Ydeg_le3': lambda e: e[3] + e[4] + e[5] <= 3,
}
for a_ in range(0, 13):
    PREDS[f'Xdeg_{a_}'] = (lambda a_: (lambda e: e[0] + e[1] + e[2] == a_))(a_)


def multiples_basis(P, d, parD, k):
    """Basis P * m for monomials m with deg P + deg m = d, parity/weight matching."""
    dp = core.degree(P); pP = core.parity_of(P); kP = sigma_weight(P)
    if d - dp < 0 or pP is None or kP is None:
        return []
    return [shift(P, m, 1) for m in wmonomials(d - dp, (parD - pP) % 2, (k - kP) % 7)]


def explicit_families(d, parD, k):
    fams = {}
    for name, P in [('K', K_block), ('G', G_block), ('M', M_block), ('K+G', add(K_block, G_block)), ('EX', EX), ('EY', EY)]:
        b = multiples_basis(P, d, parD, k)
        if b:
            fams[f'{name}-multiples'] = b
    for i in range(NV):
        b = multiples_basis(Fi[i], d, parD, k)
        if b and d >= 10:   # d <= 9 multiples of a single polar are inside J-type searches already (wave 1 / varpro J)
            fams[f'F{i}-multiples'] = b
    return fams


def subspace_kernel(d, parD, k, basis, nprimes=4):
    """E1-admissible part of span(basis) (exact polys of degree d, parity parD, weight k): returns
    dict(dim_V, dim_kernel, reps (exact polys, small height), dims_J (dim of V ∩ J_d))."""
    if not basis:
        return {'dim_V': 0, 'dim_kernel': 0, 'reps': []}
    a_ = 3; pA = 1
    tpar = (parD + pA) % 2; tk = (k + a_) % 7
    tgt = wmonomials(d + 6, tpar, tk)
    rows = {m: i for i, m in enumerate(tgt)}

    def kernel_modp(p):
        ctx = Ctx(F, A, p)
        jk, jp = kcc2.jac_block_cols(ctx, d + 6, tpar, tk)
        MJ, _ = dense_from_cols(jp, None, p, rows)
        Lm, nl = MJ.transpose().nullspace()
        Lt = nmod_mat(nl, len(rows), [int(Lm[r, c]) for c in range(nl) for r in range(len(rows))], p)
        cols = [mul(ctx.A, core.reduce_mod(B, p), p) for B in basis]
        AD, _ = dense_from_cols(cols, None, p, rows)
        Mult = Lt * AD
        Kn, nk = Mult.nullspace()
        return [[int(Kn[r, c]) for r in range(len(basis))] for c in range(nk)]
    vec0 = kernel_modp(L.P0)
    out = {'dim_V': len(basis), 'dim_kernel': len(vec0), 'reps': []}
    if not vec0:
        return out
    rat = L.crt_lift_vectors(kernel_modp, nprimes)
    out['lift_ok'] = rat is not None
    if rat is None:
        return out
    red = L.lll_small_basis([L.clear_denominators(v) for v in rat])
    reps = []
    for r in red:
        D = {}
        for c, B in zip(r, basis):
            if c:
                D = add(D, scale(B, Fraction(int(c))))
        if D:
            g = 0
            for v in D.values():
                g = math.gcd(g, abs(Fraction(v).numerator))
            reps.append(D)
    out['reps'] = reps
    # how much of the kernel lies in J_d (Jacobian type)?
    ctx = Ctx(F, A, L.P0)
    jdk, jdp = kcc2.jac_block_cols(ctx, d, parD, k)
    Dm = wmonomials(d, parD, k)
    idx = {m: i for i, m in enumerate(Dm)}
    if jdp:
        Jm, _ = dense_from_cols(jdp, None, L.P0, idx)
        rJ = Jm.rank()
        allc, _ = dense_from_cols(jdp + [core.reduce_mod(D, L.P0) for D in reps], None, L.P0, idx)
        out['nontrivial_mod_J'] = allc.rank() - rJ
    else:
        out['nontrivial_mod_J'] = len(reps)
    return out


def run_family(args):
    d, parD, k, name, kind, restarts, tag = args
    t0 = time.time()
    try:
        if kind == 'pred':
            basis = [{m: Fraction(1)} for m in wmonomials(d, parD, k) if PREDS[name](m)]
        else:
            basis = explicit_families(d, parD, k).get(name, [])
        if not basis:
            return {'task': (d, parD, k, name), 'dim_V': 0}
        info = subspace_kernel(d, parD, k, basis)
        rec = {'family': f'struct:{name}', 'd': d, 'parD': parD, 'k': k, 'dim_V': info['dim_V'], 'dim_E1': info['dim_kernel'],
               'nontrivial_mod_J': info.get('nontrivial_mod_J'), 'lift_ok': info.get('lift_ok'), 'tag': tag,
               'note': f'structured family {name}: E1-admissible subspace dimension mod p (0 => family excluded as supports)'}
        if info['dim_kernel'] and info['reps']:
            fam = []
            for idx, D in enumerate(info['reps']):
                N = L.rational_e1_lift(F, A, D, nprimes=12, restrict=True)
                if N is None:
                    N, relf = L.float_e1_lift(F, A, D, restrict=True)
                fam.append((D, N, ('S', idx)))
            kv = L.KoszulVarPro(F, A, fam, restrict=True, k=k)
            rec.update({'K': kv.K, 'dimP': kv.dimP, 'rows': kv.rows, 'cols': kv.ncols})
            if kv.dimP == 0:
                rec['note'] += '; E1 part consists of F-multiples only'
            else:
                out = kv.search_lm(n_restarts=restarts, seed=d * 100 + k, max_nfev=300)
                rels = np.array([c['rel'] for c in out])
                rec.update({'restarts': len(out), 'best_rel': float(rels.min()), 'median_rel': float(np.median(rels)),
                            'n_lt_1e-6': int((rels < 1e-6).sum()), 'n_lt_1e-10': int((rels < 1e-10).sum()),
                            'n_spurious': int(sum(1 for c in out if c['spurious']))})
                best = out[0]
                if best['rel'] < 1e-10 and not best['spurious']:
                    from w2_varpro import try_candidate
                    res = try_candidate(kv, best, d, parD, k, f'{tag}:{name}', None)
                    rec['candidate'] = res['status']
                # fixed-D exact tests of the reps themselves (tiny families): each rep is an E1-admissible structured D
                if kv.K <= 6:
                    ex = []
                    for D in info['reps']:
                        st, data = L.confirm_exact(F, A, D, restrict=True)
                        ex.append({'status': st, 'rank': data.get('rank'), 'nullity': data.get('nullity'), 'terms': len(D)})
                    rec['reps_exact'] = ex
        rec['seconds'] = round(time.time() - t0, 1)
        L.wlog(rec)
        return {'task': (d, parD, k, name), 'dim_V': rec['dim_V'], 'dim_E1': rec['dim_E1'], 'best_rel': rec.get('best_rel'), 'seconds': rec['seconds']}
    except Exception as ex:
        L.wlog({'family': f'struct:{name}', 'd': d, 'parD': parD, 'k': k, 'tag': tag, 'error': repr(ex), 'traceback': traceback.format_exc()[-1500:]})
        return {'task': (d, parD, k, name), 'error': repr(ex)}


def hessian_minors():
    """Explicit candidates: 2x2 and principal 3x3 minors of Hess(F), block determinants, their products with monomials."""
    H = [[pderiv(Fi[i], j) for j in range(NV)] for i in range(NV)]
    cands = []
    for i in range(NV):
        for j in range(i + 1, NV):
            m2 = sub(mul(H[i][i], H[j][j]), mul(H[i][j], H[i][j]))
            if m2:
                cands.append((f'minor2({i},{j})', m2))

    def det3(M):
        return add(add(sub(mul(M[0][0], sub(mul(M[1][1], M[2][2]), mul(M[1][2], M[2][1]))),
                           mul(M[0][1], sub(mul(M[1][0], M[2][2]), mul(M[1][2], M[2][0])))),
                       mul(M[0][2], sub(mul(M[1][0], M[2][1]), mul(M[1][1], M[2][0])))), {})
    for trip in itertools.combinations(range(NV), 3):
        m3 = det3([[H[i][j] for j in trip] for i in trip])
        if m3:
            cands.append((f'minor3{trip}', m3))
    return cands


def run_explicit(args):
    name, Dstr, tag = args
    t0 = time.time()
    try:
        D = core.parse(Dstr) if isinstance(Dstr, str) else Dstr
        d = core.degree(D); pD = core.parity_of(D); k = sigma_weight(D)
        rec = {'family': 'struct-explicit', 'name': name, 'd': d, 'parD': pD, 'k': k, 'terms': len(D), 'tag': tag}
        if pD is None:
            rec['note'] = 'mixed Y-parity: split into parity parts'
            parts = [{e: v for e, v in D.items() if (e[3] + e[4] + e[5]) % 2 == s} for s in (0, 1)]
        else:
            parts = [D]
        res = []
        for P in parts:
            if not P:
                continue
            if sigma_weight(P) is None:
                res.append({'parity': core.parity_of(P), 'note': 'not a sigma-eigenvector; unrestricted exact test'})
                st, data = L.confirm_exact(F, A, P, restrict=False)
            else:
                st, data = L.confirm_exact(F, A, P, restrict=True)
            res.append({'parity': core.parity_of(P), 'k': sigma_weight(P), 'status': st, 'rank': data.get('rank'), 'nullity': data.get('nullity'),
                        'rows': data.get('rows'), 'cols': data.get('cols')})
            if st == 'solved':
                wit = {'d': core.degree(P), 'D': core.to_str(P), 'N': [core.to_str(n) for n in data['N']], 'Q': core.to_str(data['Q']), 'tag': tag, 'name': name}
                json.dump(wit, open(f'{LANE}/witness-candidate.json', 'w'), indent=1)
                L.wlog({'family': 'WITNESS', 'd': core.degree(P), 'D': core.to_str(P)[:2000], 'note': f'WITNESS FOUND (explicit structured candidate {name})'})
                print('WITNESS FOUND', flush=True)
        rec['results'] = res
        rec['seconds'] = round(time.time() - t0, 1)
        L.wlog(rec)
        return {'name': name, 'd': d, 'results': [(r.get('status'), r.get('rank')) for r in res], 'seconds': rec['seconds']}
    except Exception as ex:
        L.wlog({'family': 'struct-explicit', 'name': name, 'tag': tag, 'error': repr(ex), 'traceback': traceback.format_exc()[-1500:]})
        return {'name': name, 'error': repr(ex)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--degrees', default='6,7,8,9,10,11,12')
    ap.add_argument('--procs', type=int, default=8)
    ap.add_argument('--restarts', type=int, default=200)
    ap.add_argument('--tag', default='struct')
    ap.add_argument('--explicit', action='store_true', help='also test explicit Hessian minors etc.')
    ap.add_argument('--only', default='', help='comma list of family names to run')
    ap.add_argument('--sectors', default='', help='restrict to sectors, e.g. "0:4-6,1:0-6" (parD:k-range)')
    a = ap.parse_args()
    degs = [int(x) for x in a.degrees.split(',')]
    jobs = []
    only = set(a.only.split(',')) if a.only else None
    secs = None
    if a.sectors:
        secs = set()
        for item in a.sectors.split(','):
            pd, kr = item.split(':')
            lo, hi = (kr.split('-') + [kr])[:2] if '-' in kr else (kr, kr)
            for kk in range(int(lo), int(hi) + 1):
                secs.add((int(pd), kk))
    for d in degs:
        for parD in (0, 1):
            for k in range(7):
                if secs is not None and (parD, k) not in secs:
                    continue
                for name in PREDS:
                    if name.startswith('Xdeg_') and name != 'Xdeg_le3' and int(name.split('_')[1]) > d:
                        continue
                    if only and name not in only:
                        continue
                    jobs.append((d, parD, k, name, 'pred', a.restarts, a.tag))
                for name in explicit_families(d, parD, k):
                    if only and name not in only:
                        continue
                    jobs.append((d, parD, k, name, 'explicit', a.restarts, a.tag))
    t0 = time.time()
    print(f'{len(jobs)} structured-family tasks', flush=True)
    with Pool(a.procs) as pool:
        for r in pool.imap_unordered(run_family, jobs):
            if r.get('dim_E1') or r.get('error'):
                print(json.dumps(r, default=str), f'[{time.time()-t0:.0f}s]', flush=True)
    if a.explicit:
        ex = [(n, P, a.tag) for n, P in hessian_minors()]
        with Pool(a.procs) as pool:
            for r in pool.imap_unordered(run_explicit, ex):
                print(json.dumps(r, default=str), f'[{time.time()-t0:.0f}s]', flush=True)
    print('STRUCT_DONE', flush=True)


if __name__ == '__main__':
    main()
