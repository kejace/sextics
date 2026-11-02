"""Subfamily-restricted relaxations, pure or mixed parity.  Usage: python3 w2_subfam2.py <base> [--minr 2 --maxr 5]
[--timeout 300] [--only names] [--extra json] [--tag T]. For a mixed base (parspec m) each subring family contributes its even
and odd admissible parts; also tests the pencil span(K_even, K_odd) and 'fermat' (both parities)."""
import sys, os, json, time, argparse, random, itertools
sys.path.insert(0, '/tmp/hodge/results/short-power-presentation')
from w2lib import *
from msparse import parse_ms
from w2_relax import run_msolve, run_singular, load_kernel, gauge_indices
from w2_families import form_polys, family_monomials, span_dim, ORDER

def admissible_coords(F, A, d, par, p, polys, labels_par):
    """gauge-fixed coordinate vectors (over labels_par) of the admissible part of span(polys) (pure parity par)."""
    K, src = load_kernel(d, par, p); K = K or []
    N = particular_numerators(F, A, K, p, verbose=False) if K else []
    fam = Family(F, A, d, kernel_reps=K, kernel_nums=N, parities=(par,))
    gidx = set(gauge_indices(fam))
    gens = [reduce_mod(m['G'], p) for m in fam.members]
    mons = monomials(d, par); midx = {m: i for i, m in enumerate(mons)}; nG = len(gens)
    jref = JRef(F, d, par, p, verbose=False)
    sm = jref.std; sidx = {jref.mons[j]: t for t, j in enumerate(sm)}; ns = len(sm)
    Knf = jref.nf_polys(K) if K else []
    nfs = jref.nf_polys([reduce_mod(P, p) for P in polys]); n = len(polys); k = len(Knf)
    W = n + k; fl = [0] * (ns * W)
    for i, nf in enumerate(nfs):
        for m, v in nf.items(): fl[sidx[m] * W + i] = v
    for j, kv in enumerate(Knf):
        for m, v in kv.items(): fl[sidx[m] * W + n + j] = v
    X, nul = nmod_mat(ns, W, fl, p).nullspace()
    if nul == 0: return []
    # coordinates: solve [gens] x = E for each admissible E (batch: augmented rref with all E's)
    Es = []
    for c in range(nul):
        E = {}
        for i in range(n):
            v = int(X[i, c])
            if v: E = add(E, scale(reduce_mod(polys[i], p), v, p), p)
        if E: Es.append(E)
    if not Es: return []
    R = len(mons); W2 = nG + len(Es); fl = [0] * (R * W2)
    for kk, g in enumerate(gens):
        for e, v in g.items(): fl[midx[e] * W2 + kk] = v
    for c, E in enumerate(Es):
        for e, v in E.items(): fl[midx[e] * W2 + nG + c] = v
    Mr, rank = nmod_mat(R, W2, fl, p).rref(); del fl
    piv = []; j = 0
    for i in range(rank):
        while int(Mr[i, j]) == 0: j += 1
        piv.append(j); j += 1
    vecs = []
    for c, E in enumerate(Es):
        if (nG + c) in set(piv): continue   # not in the span (should not happen)
        x = [0] * nG
        for i, j in enumerate(piv):
            if j < nG: x[j] = int(Mr[i, nG + c])
        for g in gidx:
            cc = x[g]
            if cc:
                m5 = fam.members[g]['m']; h = list(m5); h[5] -= 1; h = tuple(h)
                for jj in range(NV):
                    key = (jj, mon_add(h, tuple(1 if t == jj else 0 for t in range(NV))))
                    kk = next(k_ for k_, mm in enumerate(fam.members) if mm['tag'] == 'J' and (mm['j'], mm['m']) == key)
                    x[kk] = (x[kk] - cc) % p
        lab = [fam.members[k_]['label'] for k_ in range(nG) if k_ not in gidx]
        vecs.append({l: v for l, v in zip(lab, [x[k_] for k_ in range(nG) if k_ not in gidx]) if v})
    return vecs

def main():
    ap = argparse.ArgumentParser(); ap.add_argument('base'); ap.add_argument('--minr', type=int, default=2); ap.add_argument('--maxr', type=int, default=5)
    ap.add_argument('--timeout', type=int, default=300); ap.add_argument('--only', default=''); ap.add_argument('--tag', default='')
    ap.add_argument('--skip-subrings', action='store_true')
    a = ap.parse_args()
    meta = json.load(open(a.base + '_meta.json')); labels = meta['labels']; p = meta['p']; d = meta['d']
    pars = [0, 1] if meta['parspec'] == 'm' else [int(meta['parspec'])]
    vars_, p2, eqs = parse_ms(a.base + '.ms'); assert p2 == p
    F, A, q, rec = load(); FP = form_polys()
    pos = {lab: i for i, lab in enumerate(labels)}; Kn = len(labels)
    # per-parity label sets (for the mixed case the labels of both parities coexist; kernel labels K0.. are in Family order)
    quad = []
    for eq in eqs:
        e2 = {}
        for key, v in eq.items():
            if len(key) == 2: e2[(key[0], key[1])] = v
        if e2: quad.append(e2)
    # kernel label mapping: in the mixed Family the kernel reps are ordered even then odd (as in w2_relax: reps += K per parity)
    klabels = [l for l in labels if l.startswith('K')]
    def embed(vecs_by_par):
        out = []
        for par, vecs in vecs_by_par:
            Kpar, _ = load_kernel(d, par, p); nk = len(Kpar or [])
            # kernel labels for this parity: for mixed, even reps come first
            offset = 0 if (par == 0 or len(pars) == 1) else (len(load_kernel(d, 0, p)[0] or []) if 0 in pars else 0)
            for vec in vecs:
                full = [0] * Kn
                for lab, v in vec.items():
                    if lab.startswith('K'):
                        lab2 = f'K{int(lab[1:]) + offset}'
                    else:
                        lab2 = lab
                    full[pos[lab2]] = (full[pos[lab2]] + v) % p
                out.append(full)
        return out
    subfams = {}
    if not a.skip_subrings:
        seen = set(); import sympy as sp
        for r in range(a.minr, a.maxr + 1):
            for names in itertools.combinations(ORDER, r):
                if span_dim(names, FP) < r: continue
                Mx = sp.Matrix([[FP[n].get(tuple(1 if t == i else 0 for t in range(NV)), 0) for i in range(NV)] for n in names])
                sig = tuple(tuple(Mx.rref()[0].tolist()[i]) for i in range(r))
                if sig in seen: continue
                seen.add(sig)
                subfams['+'.join(names)] = {par: [P for _, P in family_monomials(names, d, par, FP)] for par in pars}
        subfams['XM+YL'] = {par: [P for _, P in family_monomials(['X0', 'X1', 'X2', 'M'], d, par, FP)] + [P for _, P in family_monomials(['Y0', 'Y1', 'Y2', 'L'], d, par, FP)] for par in pars}
    if a.only:
        subfams = {k: v for k, v in subfams.items() if k in a.only.split(',')}
    # special: kernel pencil and fermat family
    special = {}
    special['kernel'] = [[1 if labels[i].startswith('K') and labels[i] == kl else 0 for i in range(Kn)] for kl in klabels]
    ferm = []
    for lab in labels:
        if lab.startswith('J'):
            j = int(lab[1]); m = tuple(int(c) for c in lab.split(':')[1])
            if m[j] >= 1:
                qm = list(m); qm[j] -= 1
                # w_j = z_j * q: any monomial m with m[j] >= 1 qualifies
                v = [0] * Kn; v[pos[lab]] = 1; ferm.append(v)
    special['fermat+kernel'] = special['kernel'] + ferm
    results = {}; t0 = time.time()
    def solve_sub(name, Bvecs):
        if not Bvecs: return
        Bm = nmod_mat(len(Bvecs), Kn, [v for vec in Bvecs for v in vec], p); Br, rk = Bm.rref()
        Bcols = [[int(Br[i, k]) for k in range(Kn)] for i in range(rk)]
        sub_eqs = change_basis_quadratic(quad, Kn, Bcols, p)
        base = a.base + '_sf_' + name.replace('+', '_') + (f'_{a.tag}' if a.tag else '')
        write_singular(base + '.sing', sub_eqs, rk, p)
        rng = random.Random(1); chart = list(sub_eqs) + [{(): p - 1, **{(i,): rng.randrange(1, p) for i in range(rk)}}]
        write_msolve(base + '_chart.ms', chart, rk, p)
        r1 = run_msolve(base + '_chart.ms', base + '_chart.msout', 2, timeout=a.timeout); head = r1.get('out_head', '').strip()[:60]
        r2 = run_singular(base + '.sing', a.timeout); so = r2.get('stdout', '') if isinstance(r2.get('stdout'), str) else ''
        dimv = [l for l in so.split('\n') if l.startswith('dim') or l.startswith('vdim')]
        results[name] = {'n_vars': rk, 'n_eqs': len(sub_eqs), 'msolve_chart': head, 'singular': dimv, 'singular_rc': r2.get('rc')}
        print(f'[subfam] {name}: vars {rk}, eqs {len(sub_eqs)}, msolve chart {head!r}, singular {dimv} ({time.time()-t0:.0f}s)', file=sys.stderr, flush=True)
    for name, B in special.items():
        solve_sub(name, B)
    for name, bypar in subfams.items():
        vecs_by_par = []
        for par in pars:
            if bypar[par]:
                vecs_by_par.append((par, admissible_coords(F, A, d, par, p, bypar[par], None)))
        solve_sub(name, embed(vecs_by_par))
    logrec = {'kind': 'relaxation-subfamilies', 'base': a.base, 'd': d, 'parspec': meta['parspec'], 'p': p, 'results': results, 'seconds': round(time.time() - t0, 1), 'tag': a.tag}
    w2log(logrec); print(json.dumps(logrec)[:3000], flush=True)

if __name__ == '__main__':
    main()
