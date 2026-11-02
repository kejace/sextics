"""Probe the projective dimension of the relaxation variety V(I) mod p by random linear slices:
V ∩ (random P^k) is nonempty iff dim V >= 34-k... (generic slices).  Usage: excl_slice.py <system.json> <kmin> <kmax> [threads]"""
import sys, json, time, random, subprocess, re, os
sys.path.insert(0, '/tmp/hodge/results/short-power-presentation')
from flint import nmod_mat
path = sys.argv[1]; kmin, kmax = int(sys.argv[2]), int(sys.argv[3])
threads = sys.argv[4] if len(sys.argv) > 4 else '2'
S = json.load(open(path)); p = S['info']['p']; nu = len(S['unknowns'])
random.seed(12345)
# independent quadrics
pairs = {}; rows = []
for eq in S['eqs']:
    r = {}
    for (i, j), v in eq:
        c = pairs.setdefault((i, j), len(pairs)); r[c] = v
    rows.append(r)
R, C = len(rows), len(pairs)
flat = [0] * (R * C)
for i, r in enumerate(rows):
    for c, v in r.items(): flat[i * C + c] = v
Mr, rk = nmod_mat(R, C, flat, p).rref()
inv = {c: k for k, c in pairs.items()}
quads = []
for i in range(rk):
    Q = {}
    for c in range(C):
        v = int(Mr[i, c])
        if v: Q[inv[c]] = v
    quads.append(Q)
print(f'{rk} independent quadrics in {nu} unknowns', flush=True)
base = path.replace('.system.json', '')
for k in range(kmin, kmax + 1):
    m = k + 1   # number of slice variables t_1..t_m ; u_i = sum_a T[i][a] t_a
    T = [[random.randrange(p) for _ in range(m)] for _ in range(nu)]
    # substitute: quadric sum_{i<=j} v u_i u_j -> sum over (a,b)
    lines = []
    for Q in quads:
        coef = {}
        for (i, j), v in Q.items():
            for a in range(m):
                for b in range(m):
                    key = (min(a, b), max(a, b))
                    coef[key] = (coef.get(key, 0) + v * T[i][a] * T[j][b]) % p
        terms = [f'{v}*t{a+1}^2' if a == b else f'{v}*t{a+1}*t{b+1}' for (a, b), v in coef.items() if v]
        if terms: lines.append('+'.join(terms))
    ms = f'{base}.slice{k}.ms'
    with open(ms, 'w') as f:
        f.write(','.join(f't{a+1}' for a in range(m)) + '\n' + f'{p}\n' + ',\n'.join(lines) + '\n')
    t0 = time.time()
    out = ms + '.out'
    subprocess.run(['msolve', '-g', '1', '-t', threads, '-f', ms, '-o', out], capture_output=True, timeout=1800)
    txt = ''.join(l for l in open(out).read().splitlines() if not l.startswith('#')).replace(' ', '')
    leads = [x for x in txt.strip('[]:').split(',') if x]
    degs = [sum(int(e.split('^')[1]) if '^' in e else 1 for e in L.split('*')) for L in leads]
    # projectively empty iff a pure power of every variable is a leading monomial
    purepow = set()
    for L in leads:
        facs = L.split('*')
        if len(facs) == 1:
            purepow.add(facs[0].split('^')[0])
    empty = all(f't{a+1}' in purepow for a in range(m))
    print(f'slice P^{k} ({m} vars): {len(leads)} leading monomials, max lead degree {max(degs) if degs else 0} -> ' +
          ('EMPTY' if empty else 'NONEMPTY (projective dim >= 0 in the slice)') + f'  [{time.time()-t0:.1f}s]', flush=True)
