"""Structured subring families screened exactly against the E1-admissible space K_d + J(F)_d.
Families: monomials of degree d (Y-parity par) in r = 2..5 of the linear forms {X0,X1,X2,Y0,Y1,Y2,L,M}
(cones with vertex a P^{5-r}), plus named extras.  For each family computes dim(family ∩ (K_d+J_d)) and
dim(family ∩ J_d) mod p, and rational representatives (2 primes) of the admissible elements.
Usage: python3 w2_families.py <d> <par> [--primes 2] [--maxr 5] [--kernel path] [--tag T]"""
import sys, os, json, time, argparse, itertools
sys.path.insert(0, '/tmp/hodge/results/short-power-presentation')
from w2lib import *

FORMS = {'X0': 'X0', 'X1': 'X1', 'X2': 'X2', 'Y0': 'Y0', 'Y1': 'Y1', 'Y2': 'Y2',
         'L': 'X0 + 2*X1 + 3*X2', 'M': 'Y0 + 3*Y1 + 5*Y2'}
ORDER = ['X0', 'X1', 'X2', 'Y0', 'Y1', 'Y2', 'L', 'M']

def form_polys():
    return {k: core.parse(v) for k, v in FORMS.items()}

def span_dim(names, FP):
    import numpy as np
    M = np.array([[float(FP[n].get(tuple(1 if t == i else 0 for t in range(NV)), 0)) for i in range(NV)] for n in names])
    return int(np.linalg.matrix_rank(M))

def family_monomials(names, d, par, FP):
    """Monomials prod f_i^{e_i} (deg d, Y-parity par) as (exponent tuple, poly)."""
    pars = [core.parity_of(FP[n]) for n in names]
    r = len(names)
    out = []
    pw = {n: [ {(0,)*NV: Fraction(1)} ] for n in names}
    for n in names:
        for k in range(1, d + 1):
            pw[n].append(mul(pw[n][-1], FP[n]))
    for comb in itertools.combinations(range(d + r - 1), r - 1):
        e = []; prev = -1
        for x in comb:
            e.append(x - prev - 1); prev = x
        e.append(d + r - 2 - prev)
        if sum(ei * pi for ei, pi in zip(e, pars)) % 2 != par:
            continue
        P = {(0,)*NV: Fraction(1)}
        for n, ei in zip(names, e):
            if ei: P = mul(P, pw[n][ei])
        out.append((tuple(e), P))
    return out

def load_kernel_modp(d, par, p):
    path = os.path.join(RES, f'kernel_d{d}_p{par}.json')
    if not os.path.exists(path):
        return None
    J = json.load(open(path))
    if J.get('basis_rational'):
        return [reduce_mod(core.parse(s), p) for s in J['basis_rational']]
    if str(p) in J.get('basis_modp', {}):
        return [modp_poly_from_json(o) for o in J['basis_modp'][str(p)]]
    return None

def screen(F, A, d, par, p, families, Kp, verbose=True):
    jref = JRef(F, d, par, p, verbose=verbose)
    sm = jref.std; sidx = {jref.mons[j]: t for t, j in enumerate(sm)}; ns = len(sm)
    Knf = jref.nf_polys(Kp) if Kp else []
    results = {}
    for key, mons in families.items():
        polys = [P for _, P in mons]
        if not polys:
            continue
        nfs = jref.nf_polys([reduce_mod(P, p) for P in polys])
        n = len(polys); k = len(Knf)
        # matrix (ns x (n + k)): columns NF(fam_i), K_j ; nullspace
        flat = [0] * (ns * (n + k))
        W = n + k
        for i, nf in enumerate(nfs):
            for m, v in nf.items(): flat[sidx[m] * W + i] = v
        for j, kv in enumerate(Knf):
            for m, v in kv.items(): flat[sidx[m] * W + n + j] = v
        Mx = nmod_mat(ns, W, flat, p)
        X, nul = Mx.nullspace()
        vecs = []
        for c in range(nul):
            x = [int(X[i, c]) for i in range(n)]; y = [int(X[n + j, c]) for j in range(k)]
            vecs.append((x, y))
        # canonical form: rref of the x-parts (family coordinates)
        if nul:
            Xm = nmod_mat(nul, n, [v for x, y in vecs for v in x], p)
            Xr, rk = Xm.rref()
            xs = [[int(Xr[i, j]) for j in range(n)] for i in range(rk)]
            Ym = nmod_mat(nul, k, [v for x, y in vecs for v in y], p) if k else None
            rank_y = Ym.rank() if k else 0
        else:
            xs, rk, rank_y = [], 0, 0
        results[key] = {'n': n, 'nullity': nul, 'rank_x': rk, 'rank_y': rank_y, 'xs': xs}
        if verbose and nul:
            print(f'[families] d={d} par={par} {key}: n={n} admissible dim={rk} (kernel part rank {rank_y})', file=sys.stderr, flush=True)
    return results

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('d', type=int); ap.add_argument('par', type=int)
    ap.add_argument('--primes', type=int, default=2); ap.add_argument('--maxr', type=int, default=5)
    ap.add_argument('--minr', type=int, default=2); ap.add_argument('--tag', default='')
    ap.add_argument('--threads', type=int, default=4)
    a = ap.parse_args()
    set_threads(a.threads)
    F, A, q, rec = load()
    d, par = a.d, a.par
    FP = form_polys()
    t0 = time.time()
    families = {}
    seen_spans = set()
    import numpy as np
    for r in range(a.minr, a.maxr + 1):
        for names in itertools.combinations(ORDER, r):
            if span_dim(names, FP) < r:
                continue
            # dedupe by the span of the linear forms (row space signature via rref over Q)
            import sympy as sp
            Mx = sp.Matrix([[FP[n].get(tuple(1 if t == i else 0 for t in range(NV)), 0) for i in range(NV)] for n in names])
            sig = tuple(tuple(Mx.rref()[0].tolist()[i]) for i in range(r))
            if sig in seen_spans:
                continue
            seen_spans.add(sig)
            families['+'.join(names)] = family_monomials(names, d, par, FP)
    nfam = len(families)
    print(f'[families] d={d} par={par}: {nfam} families built ({time.time()-t0:.0f}s)', file=sys.stderr, flush=True)
    primes = (PGB + core.PRIMES)[:a.primes]
    allres = {}
    for p in primes:
        Kp = load_kernel_modp(d, par, p)
        if Kp is None:
            print(f'[families] no kernel basis for p={p}; screening against J_d only', file=sys.stderr)
            Kp = []
        allres[p] = screen(F, A, d, par, p, families, Kp)
    # reconstruct rational representatives where nontrivial and consistent across primes
    out = {}
    p0 = primes[0]
    for key, r0 in allres[p0].items():
        rec_ = {'n': r0['n'], 'admissible_dim': r0['rank_x'], 'kernel_part_rank': r0['rank_y'],
                'consistent_primes': all(allres[p][key]['rank_x'] == r0['rank_x'] for p in primes)}
        if r0['rank_x'] and len(primes) >= 2 and rec_['consistent_primes']:
            mons = families[key]
            lists = {p: [ {mons[j][0]: v for j, v in enumerate(allres[p][key]['xs'][i]) if v} for i in range(r0['rank_x'])] for p in primes}
            rat, info = reconstruct_polys(lists)
            if rat:
                elems = []
                for vec in rat:
                    P = {}
                    for e, c in vec.items():
                        # find the poly for exponent e
                        for ee, PP in mons:
                            if ee == e:
                                P = add(P, scale(PP, c)); break
                    elems.append({'coords': {' '.join(map(str, e)): str(c) for e, c in vec.items()}, 'D': core.to_str(P)})
                rec_['rational'] = elems
            else:
                rec_['rational'] = None; rec_['recon_info'] = info
        out[key] = rec_
    path = os.path.join(RES, f'families_d{d}_p{par}' + (f'_{a.tag}' if a.tag else '') + '.json')
    with open(path, 'w') as f:
        json.dump({'d': d, 'par': par, 'primes': primes, 'families': out, 'forms': FORMS}, f, indent=1)
    nontriv = {k: (v['admissible_dim'], v['kernel_part_rank']) for k, v in out.items() if v['admissible_dim']}
    logrec = {'kind': 'families-screen', 'd': d, 'parity': par, 'n_families': nfam, 'primes': primes,
              'nontrivial': nontriv, 'out': path, 'seconds': round(time.time() - t0, 1), 'tag': a.tag}
    w2log(logrec)
    print(json.dumps(logrec), flush=True)

if __name__ == '__main__':
    main()
