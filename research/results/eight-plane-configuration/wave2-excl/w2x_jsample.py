"""Sampling of individual Jacobian denominators D = w(F) (d = 7, 8, 9) with the Koszul-reduced E2 system mod p.
Candidates: random fields (generic inconsistency), Hamiltonian fields w = (G_j d_i - G_i d_j) and C grad G for named/generated G,
gradient fields w = grad G, arrangement-preserving fields with random coefficients, block fields.
Each candidate is a fixed D: modular inconsistency with rank_p M = ncols would be rigorous; otherwise two primes are reported.
Usage: python3 w2x_jsample.py <d> <nworkers> <tag> [p]
"""
import sys, time, json, random, itertools
from multiprocessing import Pool
sys.path.insert(0, '/tmp/hodge/results/eight-plane-configuration/wave2-excl')
from w2x_lib import *
from w2x_koszul import solve_mod
from w2x_hess import field_basis
from w2x_hess_sub import named_forms, subfamily

def gen_G(deg, par, p):
    """named forms of degree deg and parity par, plus monomial multiples of lower named forms."""
    nf = named_forms(p)
    out = {}
    for name, G in nf.items():
        dg = degree(G)
        if dg == deg and parity_of(G) == par: out[name] = G
        elif dg < deg:
            for m in monomials(deg - dg, (par - parity_of(G)) % 2):
                ms = core.to_str({m: 1})
                out[f'{ms}*{name}'] = shift(G, m, 1, p)
    return out

def candidates(d, p, seed=1):
    random.seed(seed)
    cands = []
    deg = d - 5
    for par in (0, 1):
        B = field_basis(d, par)
        # random fields
        for r in range(6):
            w = [{} for _ in range(NV)]
            for (i, m) in B:
                c = random.randint(-3, 3)
                if c: w[i][m] = c % p
            cands.append((f'random_d{d}_par{par}_r{r}', 'random', par, w))
        # arrangement-preserving fields, random coefficients
        arr = subfamily('arr', d, par, p)
        for r in range(3):
            w = [{} for _ in range(NV)]
            for fld in arr:
                c = random.randint(-2, 2) % p
                if c:
                    for i in range(NV): w[i] = add(w[i], scale(fld[i], c, p), p)
            if any(w): cands.append((f'arr_random_d{d}_par{par}_r{r}', 'arr-random', par, w))
        # Hamiltonian pair fields and gradient fields for generated G of degree deg+1
        Gs = gen_G(deg + 1, 0, p); Gs.update({k + '[odd]': v for k, v in gen_G(deg + 1, 1, p).items()})
        names = sorted(Gs)
        random.shuffle(names)
        used = 0
        for name in names:
            G = Gs[name]; pG = parity_of(G); Gi = [pderiv(G, i, p) for i in range(NV)]
            # gradient field: parity(w_{X_i}) = parity(G_{X_i}) = pG must equal par
            if pG == par and any(Gi):
                cands.append((f'grad({name})_d{d}', 'gradient', par, Gi)); used += 1
            # Hamiltonian pairs
            pairs = []
            for i in range(NV):
                for j in range(i + 1, NV):
                    ok_i = ((pG + (1 if j >= 3 else 0)) % 2) == ((par + (1 if i >= 3 else 0)) % 2)
                    ok_j = ((pG + (1 if i >= 3 else 0)) % 2) == ((par + (1 if j >= 3 else 0)) % 2)
                    if ok_i and ok_j: pairs.append((i, j))
            random.shuffle(pairs)
            for (i, j) in pairs[:2]:
                w = [{} for _ in range(NV)]; w[i] = Gi[j]; w[j] = scale(Gi[i], p - 1, p)
                if any(w): cands.append((f'ham({name};{core.VARS[i]},{core.VARS[j]})_d{d}', 'hamiltonian-pair', par, w)); used += 1
            if pairs:
                # full random antisymmetric C
                w = [{} for _ in range(NV)]
                for (i, j) in pairs:
                    c = random.randint(-2, 2) % p
                    if c:
                        w[i] = add(w[i], scale(Gi[j], c, p), p); w[j] = add(w[j], scale(Gi[i], (p - c) % p, p), p)
                if any(w): cands.append((f'ham({name};randC)_d{d}', 'hamiltonian-randC', par, w)); used += 1
            if used > 60: break
    return cands

def run_one(args):
    name, fam, par, w, p, d = args
    try:
        res = solve_mod(w, p)
        D_ok = res['d'] == d
    except Exception as ex:
        res = {'error': f'{type(ex).__name__}: {ex}'}; D_ok = False
    return {'exp': 'c-jsample-koszulE2', 'D': name, 'family': fam, 'parD': par, **res, 'prime': p, 'degree_ok': D_ok,
            'result': ('inconsistent' if res.get('consistent') is False else 'CONSISTENT' if res.get('consistent') else 'error'),
            'rigorous_rank_full': (res.get('rank') == res.get('cols')) if 'rank' in res else None}

if __name__ == '__main__':
    d = int(sys.argv[1]); nw = int(sys.argv[2]); tag = sys.argv[3]; p = int(sys.argv[4]) if len(sys.argv) > 4 else P_NEW
    cands = candidates(d, p)
    print(f'{len(cands)} candidates at d={d}', flush=True)
    logfn = f'{PK}/wave2-excl-remote-{tag}.jsonl'
    with Pool(nw) as pool:
        for e in pool.imap_unordered(run_one, [(n, f, par, w, p, d) for (n, f, par, w) in cands]):
            e['job'] = tag
            print(json.dumps(e), flush=True)
            wlog(e, fn=logfn)
