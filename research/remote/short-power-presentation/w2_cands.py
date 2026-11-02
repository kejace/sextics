"""Named structured candidates screened exactly against K_d + J_d, plus kernel-structure analysis.
Usage: python3 w2_cands.py <d> <par> [--plist p1,p2] [--tag T]
Candidates (degree d, parity par): all 2x2 / 3x3 minors of Hess(F) (degrees 8 / 12), products of lower-degree kernel
representatives with monomials (S_{d-7} * D7, S_{d-7} * D7'), products D7^a D7'^b, block Hessian determinants,
and E_X(F)/E_Y(F)-type multiples.  Reports which candidates lie in K_d + J_d, and the dimension of the 'new' part
K_d / (S_1 K_{d-1} + J_d) using kernel_d{d-1}_p*.json."""
import sys, os, json, time, argparse, itertools
sys.path.insert(0, '/tmp/hodge/results/short-power-presentation')
from w2lib import *

def load_kernel_rational(d, par):
    path = os.path.join(RES, f'kernel_d{d}_p{par}.json')
    if not os.path.exists(path): return None
    J = json.load(open(path))
    return [core.parse(s) for s in J['basis_rational']] if J.get('basis_rational') else None

def load_kernel_modp(d, par, p):
    path = os.path.join(RES, f'kernel_d{d}_p{par}.json')
    if not os.path.exists(path): return None
    J = json.load(open(path))
    if J.get('basis_rational'):
        return [reduce_mod(core.parse(s), p) for s in J['basis_rational']]
    if str(p) in J.get('basis_modp', {}):
        return [modp_poly_from_json(o) for o in J['basis_modp'][str(p)]]
    return None

def det3(M):
    return (mul(M[0][0], sub(mul(M[1][1], M[2][2]), mul(M[1][2], M[2][1])))
            and add(sub(mul(M[0][0], sub(mul(M[1][1], M[2][2]), mul(M[1][2], M[2][1]))),
                        mul(M[0][1], sub(mul(M[1][0], M[2][2]), mul(M[1][2], M[2][0])))),
                    mul(M[0][2], sub(mul(M[1][0], M[2][1]), mul(M[1][1], M[2][0])))))

def candidates(F, A, d, par):
    Fi = [pderiv(F, i) for i in range(NV)]
    H = [[pderiv(Fi[i], j) for j in range(NV)] for i in range(NV)]
    out = {}
    def addc(name, P):
        if P and core.degree(P) == d and core.parity_of(P) == par:
            out[name] = P
    if d == 8:
        for (i, j) in itertools.combinations(range(NV), 2):
            for (k, l) in itertools.combinations(range(NV), 2):
                addc(f'hessminor[{i}{j}|{k}{l}]', sub(mul(H[i][k], H[j][l]), mul(H[i][l], H[j][k])))
    if d == 12:
        for rows in itertools.combinations(range(NV), 3):
            for cols in itertools.combinations(range(NV), 3):
                M = [[H[r][c] for c in cols] for r in rows]
                addc(f'hessminor3[{"".join(map(str,rows))}|{"".join(map(str,cols))}]', det3(M))
    # products of d=7 kernel cones with monomials / each other
    K7e = load_kernel_rational(7, 0); K7o = load_kernel_rational(7, 1)
    cone_o = None
    if os.path.exists(os.path.join(RES, 'D7odd_cone.json')):
        cone_o = core.parse(json.load(open(os.path.join(RES, 'D7odd_cone.json')))['D_cone'])
    reps = {}
    if K7e: reps['D7'] = K7e[0]
    if cone_o: reps["D7'"] = cone_o
    elif K7o: reps["D7'"] = K7o[0]
    for name, D in reps.items():
        if d - 7 >= 0:
            for m in monomials(d - 7, None):
                addc(f'{name}*{"".join(map(str, m))}', shift(D, m, 1))
    if d == 14 and len(reps) == 2:
        addc("D7*D7'", mul(reps['D7'], reps["D7'"])); addc('D7^2', mul(reps['D7'], reps['D7'])); addc("D7'^2", mul(reps["D7'"], reps["D7'"]))
    # block radial derivatives times forms
    EX = {}; EY = {}
    for i in range(3):
        EX = add(EX, shift(Fi[i], tuple(1 if t == i else 0 for t in range(NV)), 1))
        EY = add(EY, shift(Fi[3 + i], tuple(1 if t == 3 + i else 0 for t in range(NV)), 1))
    for name, D in (('EX(F)', EX), ('EY(F)', EY)):
        for m in monomials(d - 6, None):
            addc(f'{name}*{"".join(map(str, m))}', shift(D, m, 1))
    return out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('d', type=int); ap.add_argument('par', type=int)
    ap.add_argument('--plist', default=''); ap.add_argument('--tag', default='')
    a = ap.parse_args()
    F, A, q, rec = load()
    d, par = a.d, a.par
    primes = [int(x) for x in a.plist.split(',')] if a.plist else PGB[:2]
    t0 = time.time()
    cands = candidates(F, A, d, par)
    print(f'[cands] d={d} par={par}: {len(cands)} candidates', file=sys.stderr, flush=True)
    res = {}
    newdim = {}
    for p in primes:
        jref = JRef(F, d, par, p)
        sm = jref.std; sidx = {jref.mons[j]: t for t, j in enumerate(sm)}; ns = len(sm)
        Kp = load_kernel_modp(d, par, p) or []
        Knf = jref.nf_polys(Kp) if Kp else []
        # kernel structure: S_1 * K_{d-1} (both parities contribute to parity par) inside K_d mod J_d
        prev = []
        for pp in (0, 1):
            Kprev = load_kernel_modp(d - 1, pp, p)
            if Kprev:
                for D in Kprev:
                    for i in range(NV):
                        if (pp + (1 if i >= 3 else 0)) % 2 == par:
                            prev.append(shift(D, tuple(1 if t == i else 0 for t in range(NV)), 1, p))
        if prev and Knf:
            pnf = jref.nf_polys(prev)
            def mat(polys):
                return nmod_mat(len(polys), ns, [v.get(jref.mons[sm[t]], 0) for v in polys for t in range(ns)], p)
            r_prev = mat(pnf).rank() if pnf else 0
            r_all = mat(pnf + Knf).rank()
            newdim[p] = {'dim_K_d': len(Knf), 'dim_S1K_{d-1}': r_prev, 'dim_sum': r_all, 'new': len(Knf) - r_prev}
        # candidates
        names = list(cands)
        nfs = jref.nf_polys([reduce_mod(cands[n], p) for n in names]) if names else []
        for name, nf in zip(names, nfs):
            inJ = (not nf)
            if inJ:
                res.setdefault(name, {})[p] = 'in J_d'
                continue
            if not Knf:
                res.setdefault(name, {})[p] = 'not in J_d (no kernel basis available)'; continue
            flat = [0] * (ns * (len(Knf) + 1))
            W = len(Knf) + 1
            for j, kv in enumerate(Knf):
                for m, v in kv.items(): flat[sidx[m] * W + j] = v
            for m, v in nf.items(): flat[sidx[m] * W + len(Knf)] = v
            Mx = nmod_mat(ns, W, flat, p)
            r1 = Mx.rank()
            r0 = nmod_mat(ns, len(Knf), [flat[t * W + j] for t in range(ns) for j in range(len(Knf))], p).rank()
            res.setdefault(name, {})[p] = 'E1-admissible (in K_d + J_d, not in J_d)' if r1 == r0 else 'fails E1'
    summary = {}
    for name, byp in res.items():
        vals = set(byp.values()); summary[name] = vals.pop() if len(vals) == 1 else str(byp)
    adm = {n: s for n, s in summary.items() if s.startswith('E1') or s.startswith('in J')}
    path = os.path.join(RES, f'cands_d{d}_p{par}' + (f'_{a.tag}' if a.tag else '') + '.json')
    json.dump({'d': d, 'par': par, 'primes': primes, 'results': summary, 'kernel_structure': newdim}, open(path, 'w'), indent=1)
    logrec = {'kind': 'named-candidates', 'd': d, 'parity': par, 'n_candidates': len(cands), 'admissible': adm,
              'n_fail_e1': sum(1 for s in summary.values() if s == 'fails E1'), 'kernel_structure': newdim,
              'seconds': round(time.time() - t0, 1), 'tag': a.tag}
    w2log(logrec); print(json.dumps(logrec), flush=True)

if __name__ == '__main__':
    main()
