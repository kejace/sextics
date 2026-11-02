"""Restrict a built relaxation (relax_<base>.ms + _meta.json) to a linear subfamily of the gauge-fixed a-space and solve.
Subfamilies: 'fermat' = {w_j = z_j q_j} (+ kernel directions), 'kernel' = kernel directions only, or a JSON list of
label->coefficient dicts.  Usage: python3 w2_relax_sub.py <base> <subfamily> [--tag T]"""
import sys, os, json, time, argparse, subprocess, random
sys.path.insert(0, '/tmp/hodge/results/short-power-presentation')
from w2lib import *
from msparse import parse_ms
from w2_relax import run_msolve, run_singular

def main():
    ap = argparse.ArgumentParser(); ap.add_argument('base'); ap.add_argument('sub'); ap.add_argument('--tag', default='')
    a = ap.parse_args()
    meta = json.load(open(a.base + '_meta.json')); labels = meta['labels']; p = meta['p']; par = int(meta['parspec']) if meta['parspec'] != 'm' else None
    vars_, p2, eqs = parse_ms(a.base + '.ms'); assert p2 == p
    K = len(labels); pos = {lab: i for i, lab in enumerate(labels)}
    vecs = []; names = []
    for lab in labels:
        if lab.startswith('K'):
            v = [0] * K; v[pos[lab]] = 1; vecs.append(v); names.append(lab)
    if a.sub == 'fermat':
        d = meta['d']
        qpar = par if par is not None else 0
        for j in range(NV):
            for qm in monomials(d - 6, (par + (0 if j < 3 else 1) + (0 if j < 3 else 1)) % 2 if False else None):
                pass
        # even D: q_j in S_{d-6}^{par}; odd D: q_j in S_{d-6}^{par} as well (w_j = z_j q_j has parity par + par_j... check: parity(w_j) must be par + par_j)
        for j in range(NV):
            pj = 0 if j < 3 else 1
            for qm in monomials(d - 6, (par + pj + pj) % 2):   # z_j q_j parity = pj + par(q) = par + pj  => par(q) = par
                m = mon_add(qm, tuple(1 if t == j else 0 for t in range(NV)))
                lab = f'J{j}:{"".join(map(str, m))}'
                if lab in pos:
                    v = [0] * K; v[pos[lab]] = 1; vecs.append(v); names.append(lab)
    elif a.sub == 'kernel':
        pass
    else:
        for item in json.load(open(a.sub)):
            v = [0] * K
            for lab, c in item.items(): v[pos[lab]] = int(Fraction(c)) % p
            vecs.append(v); names.append(str(item)[:40])
    n = len(vecs)
    # equations as {(k,l): v} with k<=l (parse_ms gives keys as sorted tuples of variable indices)
    eqs2 = []
    for eq in eqs:
        e2 = {}
        for key, v in eq.items():
            if len(key) == 2: e2[(key[0], key[1])] = v
            elif len(key) == 1: e2[(key[0],)] = v
            else: e2[()] = v
        eqs2.append(e2)
    quad = [e for e in eqs2 if all(len(k) == 2 for k in e)]
    cols = [[vecs[c][k] for c in range(n)] for k in range(K)]
    Bcols = [[vecs[c][k] for k in range(K)] for c in range(n)]
    sub_eqs = change_basis_quadratic(quad, K, Bcols, p)
    base = a.base + '_sub_' + a.sub + (f'_{a.tag}' if a.tag else '')
    write_msolve(base + '.ms', sub_eqs, n, p)
    write_singular(base + '.sing', sub_eqs, n, p, extra='if (dim(G) >= 1 and dim(G) <= 4) { LIB "primdec.lib"; list L = minAssGTZ(G); print("ncomp = " + string(size(L))); L; }\n')
    rng = random.Random(1)
    chart = list(sub_eqs) + [{(): p - 1, **{(i,): rng.randrange(1, p) for i in range(n)}}]
    write_msolve(base + '_chart.ms', chart, n, p)
    print(f'[sub] {a.base} sub={a.sub}: {n} variables ({names[:4]}...), {len(sub_eqs)} equations', flush=True)
    r2 = run_singular(base + '.sing', 1800); print('[singular hom]', json.dumps(r2)[:1500], flush=True)
    r1 = run_msolve(base + '_chart.ms', base + '_chart.msout', 2); print('[msolve chart]', r1.get('out_head', '')[:300], r1['seconds'], flush=True)
    w2log({'kind': 'relaxation-subfamily', 'base': a.base, 'sub': a.sub, 'n_vars': n, 'n_equations': len(sub_eqs), 'p': p,
           'singular': r2.get('stdout', '')[:600], 'msolve_chart_head': r1.get('out_head', '')[:200]})

if __name__ == '__main__':
    main()
