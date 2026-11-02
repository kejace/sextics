"""Restrict a built relaxation (base.ms, base_meta.json) to the admissible part of structured subfamilies and solve.
Subfamilies: all subring families of monomials in r forms (r in minr..maxr) of {X0,X1,X2,Y0,Y1,Y2,L,M}, the two-cone family,
and (optionally) explicit polynomial lists.  Coordinates (s, w) of an admissible element E (E = sum s_i K_i + w(F)) are
computed mod p by solving the linear system over the kernel reps and the Jacobian generators.
Usage: python3 w2_subfam.py <base> [--minr 2] [--maxr 5] [--singular-timeout 600] [--only name] [--tag T]"""
import sys, os, json, time, argparse, random, itertools
sys.path.insert(0, '/tmp/hodge/results/short-power-presentation')
from w2lib import *
from msparse import parse_ms
from w2_relax import run_msolve, run_singular, load_kernel, gauge_indices
from w2_families import form_polys, family_monomials, span_dim, ORDER

def main():
    ap = argparse.ArgumentParser(); ap.add_argument('base'); ap.add_argument('--minr', type=int, default=2); ap.add_argument('--maxr', type=int, default=5)
    ap.add_argument('--singular-timeout', type=int, default=600); ap.add_argument('--only', default=''); ap.add_argument('--tag', default='')
    ap.add_argument('--minvars', type=int, default=1)
    a = ap.parse_args()
    meta = json.load(open(a.base + '_meta.json')); labels = meta['labels']; p = meta['p']; d = meta['d']
    par = int(meta['parspec']); vars_, p2, eqs = parse_ms(a.base + '.ms'); assert p2 == p
    F, A, q, rec = load(); FP = form_polys()
    K, src = load_kernel(d, par, p); K = K or []
    N = particular_numerators(F, A, K, p, verbose=False) if K else []
    fam = Family(F, A, d, kernel_reps=K, kernel_nums=N, parities=(par,))
    assert [m['label'] for m in fam.members if m['label'] in set(labels)] and fam.K == len(labels) + len(gauge_indices(fam))
    pos = {lab: i for i, lab in enumerate(labels)}
    # generator matrix for coordinates: columns = kernel reps (mod p) then Jacobian generators G_k (mod p)
    gens = [reduce_mod(m['G'], p) for m in fam.members]
    mons = monomials(d, par); midx = {m: i for i, m in enumerate(mons)}
    nG = len(gens)
    Gm = nmod_mat(len(mons), nG, [0] * (len(mons) * nG), p)
    flat = [0] * (len(mons) * nG)
    for k, g in enumerate(gens):
        for e, v in g.items(): flat[midx[e] * nG + k] = v
    Gm = nmod_mat(len(mons), nG, flat, p); del flat
    jref = JRef(F, d, par, p, verbose=False)
    sm = jref.std; sidx = {jref.mons[j]: t for t, j in enumerate(sm)}; ns = len(sm)
    Knf = jref.nf_polys(K) if K else []
    gidx = set(gauge_indices(fam))
    def coords(E):
        """(s, w) coordinates mod p of E in K_d + J_d (None if not admissible); gauge-fixed vector over labels."""
        # solve Gm x = E  (consistent iff E in K_d + J_d)
        Ev = [0] * len(mons)
        for e, v in E.items(): Ev[midx[e]] = v % p
        aug = nmod_mat(len(mons), nG + 1, [0] * (len(mons) * (nG + 1)), p)
        fl = [0] * (len(mons) * (nG + 1))
        for i in range(len(mons)):
            for k in range(nG):
                v = int(Gm[i, k])
                if v: fl[i * (nG + 1) + k] = v
            if Ev[i]: fl[i * (nG + 1) + nG] = Ev[i]
        Mr, rank = nmod_mat(len(mons), nG + 1, fl, p).rref()
        piv = []; j = 0
        for i in range(rank):
            while int(Mr[i, j]) == 0: j += 1
            piv.append(j); j += 1
        if piv and piv[-1] == nG: return None
        x = [0] * nG
        for i, j in enumerate(piv): x[j] = int(Mr[i, nG])
        # gauge fix: subtract gauge multiples to zero the (5, h Y2) coordinates
        for g in gidx:
            c = x[g]
            if c:
                m5 = fam.members[g]['m']; h = list(m5); h[5] -= 1; h = tuple(h)
                for j in range(NV):
                    key = (j, mon_add(h, tuple(1 if t == j else 0 for t in range(NV))))
                    kk = next(k for k, mm in enumerate(fam.members) if mm['tag'] == 'J' and (mm['j'], mm['m']) == key)
                    x[kk] = (x[kk] - c) % p
        return [x[k] for k in range(nG) if k not in gidx]
    families = {}
    seen = set()
    import sympy as sp
    for r in range(a.minr, a.maxr + 1):
        for names in itertools.combinations(ORDER, r):
            if span_dim(names, FP) < r: continue
            Mx = sp.Matrix([[FP[n].get(tuple(1 if t == i else 0 for t in range(NV)), 0) for i in range(NV)] for n in names])
            sig = tuple(tuple(Mx.rref()[0].tolist()[i]) for i in range(r))
            if sig in seen: continue
            seen.add(sig)
            families['+'.join(names)] = [P for _, P in family_monomials(names, d, par, FP)]
    families['XM+YL'] = [P for _, P in family_monomials(['X0', 'X1', 'X2', 'M'], d, par, FP)] + [P for _, P in family_monomials(['Y0', 'Y1', 'Y2', 'L'], d, par, FP)]
    if a.only:
        families = {k: v for k, v in families.items() if k in a.only.split(',')}
    quad = []
    for eq in eqs:
        e2 = {}
        for key, v in eq.items():
            if len(key) == 2: e2[(key[0], key[1])] = v
        if e2: quad.append(e2)
    Kn = len(labels)
    results = {}
    t0 = time.time()
    for name, polys in families.items():
        # admissible subspace: nullspace of [NF(fam) | K] -> elements
        nfs = jref.nf_polys([reduce_mod(P, p) for P in polys]); n = len(polys); k = len(Knf)
        W = n + k; fl = [0] * (ns * W)
        for i, nf in enumerate(nfs):
            for m, v in nf.items(): fl[sidx[m] * W + i] = v
        for j, kv in enumerate(Knf):
            for m, v in kv.items(): fl[sidx[m] * W + n + j] = v
        X, nul = nmod_mat(ns, W, fl, p).nullspace()
        if nul == 0:
            continue
        vecs = []
        for c in range(nul):
            E = {}
            for i in range(n):
                v = int(X[i, c])
                if v: E = add(E, scale(reduce_mod(polys[i], p), v, p), p)
            if E:
                cv = coords(E)
                if cv is not None: vecs.append(cv)
        if not vecs: continue
        # independent columns
        Bm = nmod_mat(len(vecs), Kn, [v for vec in vecs for v in vec], p)
        Br, rk = Bm.rref()
        Bcols = [[int(Br[i, k]) for k in range(Kn)] for i in range(rk)]
        if rk < a.minvars: continue
        sub_eqs = change_basis_quadratic(quad, Kn, Bcols, p)
        base = a.base + '_fam_' + name.replace('+', '_') + (f'_{a.tag}' if a.tag else '')
        write_singular(base + '.sing', sub_eqs, rk, p)
        rng = random.Random(1)
        chart = list(sub_eqs) + [{(): p - 1, **{(i,): rng.randrange(1, p) for i in range(rk)}}]
        write_msolve(base + '_chart.ms', chart, rk, p)
        r1 = run_msolve(base + '_chart.ms', base + '_chart.msout', 2)
        head = r1.get('out_head', '').strip()[:60]
        r2 = run_singular(base + '.sing', a.singular_timeout)
        so = r2.get('stdout', '') if isinstance(r2.get('stdout'), str) else ''
        dimv = [l for l in so.split('\n') if l.startswith('dim') or l.startswith('vdim')]
        results[name] = {'n_vars': rk, 'n_eqs': len(sub_eqs), 'msolve_chart': head, 'singular': dimv, 'singular_rc': r2.get('rc')}
        print(f'[subfam] {name}: vars {rk}, eqs {len(sub_eqs)}, msolve chart {head!r}, singular {dimv} ({time.time()-t0:.0f}s)', file=sys.stderr, flush=True)
    logrec = {'kind': 'relaxation-subfamilies', 'base': a.base, 'd': d, 'parity': par, 'p': p, 'results': results, 'seconds': round(time.time() - t0, 1), 'tag': a.tag}
    w2log(logrec); print(json.dumps(logrec)[:3000], flush=True)

if __name__ == '__main__':
    main()
