"""Projective emptiness of the degree-7 L-relaxation zero set over F_p (thread's chart + square method, re-implemented).
Input: relax json/npz prefix (constraint matrix C over F_p, pairs, parameter basis).  Steps:
 1. charts c_0 = ... = c_{k-1} = 0, c_k = 1 for k = 0..K-1: Singular slimgb must give the unit ideal (dim -1);
 2. restrict to c_0 = ... = c_{K-1} = 0, iterate square elimination (c_a^2 in the row span => c_a = 0);
 3. homogeneous GB of the restricted ideal: affine dimension 0 => projective zero set empty.
With equal L-ranks over Q and F_p (flatness, verified separately) this transfers to C (thread's argument).
usage: d7_charts.py <prefix> [K=9]"""
import sys, json, subprocess, time, re
import numpy as np
from flint import nmod_mat
prefix = sys.argv[1]; K = int(sys.argv[2]) if len(sys.argv) > 2 else 9
J = json.load(open(prefix + '.json')); Z = np.load(prefix + '.npz')
C = Z['C'] % J['p']; pairs = [tuple(x) for x in Z['pairs']]; p = J['p']; n = J['params']
M = nmod_mat(C.shape[0], C.shape[1], C.ravel().tolist(), p)
R, rank = M.rref()
rows = [[int(R[i, j]) for j in range(C.shape[1])] for i in range(rank)]
print('constraint rank', rank, 'params', n, flush=True)
def polys(rowlist, keep_pairs):
    out = []
    for row in rowlist:
        terms = [f'{row[j]}*c{a}*c{b}' if a != b else f'{row[j]}*c{a}^2' for j, (a, b) in enumerate(keep_pairs) if row[j]]
        if terms: out.append('+'.join(terms))
    return out
t0 = time.time()
# step 1: charts
vars_ = ','.join(f'c{i}' for i in range(n))
sing = [f'ring r = {p}, ({vars_}), dp;', 'option(redSB);', 'ideal I = ' + ',\n'.join(polys(rows, pairs)) + ';',
        'ideal Z; ideal Jk; ideal G; int k;', f'for (k=1; k<={K}; k++) {{', '  Jk = I + Z + ideal(var(k)-1); G = slimgb(Jk);',
        '  print("CHART=" + string(k-1) + " DIM=" + string(dim(G)) + " SIZE=" + string(size(G)));', '  Z = Z + ideal(var(k)); }', 'quit;']
open(prefix + '-charts.sing', 'w').write('\n'.join(sing) + '\n')
out = subprocess.run(['Singular', '-q', prefix + '-charts.sing'], capture_output=True, text=True, timeout=7200)
charts = re.findall(r'CHART=(\d+) DIM=(-?\d+) SIZE=(\d+)', out.stdout)
print('charts:', charts, round(time.time() - t0, 1), 's', flush=True)
ok_charts = all(d == '-1' for (_, d, _) in charts) and len(charts) == K
# step 2: squares on the restricted system
killed = set(range(K)); steps = []
while True:
    cols = [j for j, (a, b) in enumerate(pairs) if a not in killed and b not in killed]
    sub = np.ascontiguousarray(C[:, cols]) % p
    Ms = nmod_mat(sub.shape[0], sub.shape[1], sub.ravel().tolist(), p); Rs, rk = Ms.rref()
    new = set(); reduced = []
    for i in range(rk):
        row = [int(Rs[i, j]) for j in range(len(cols))]; reduced.append(row)
        nz = [j for j, v in enumerate(row) if v]
        if len(nz) == 1 and pairs[cols[nz[0]]][0] == pairs[cols[nz[0]]][1]: new.add(pairs[cols[nz[0]]][0])
    steps.append({'remaining': n - len(killed), 'rank': rk, 'forced': sorted(new)}); print('squares', steps[-1], flush=True)
    if not new: break
    killed |= new
# step 3: GB of the restricted homogeneous ideal
surv = [i for i in range(n) if i not in killed]
keep_pairs = [pairs[j] for j in cols]
sing2 = [f'ring r = {p}, ({",".join(f"c{i}" for i in surv)}), dp;', 'option(redSB);', 'ideal I = ' + ',\n'.join(polys(reduced, keep_pairs)) + ';',
         'ideal G = slimgb(I);', 'print("DIM=" + string(dim(G)) + " SIZE=" + string(size(G)) + " VDIM=" + string(vdim(G)));', 'quit;']
open(prefix + '-reduced.sing', 'w').write('\n'.join(sing2) + '\n')
out2 = subprocess.run(['Singular', '-q', prefix + '-reduced.sing'], capture_output=True, text=True, timeout=7200)
m = re.search(r'DIM=(-?\d+) SIZE=(\d+) VDIM=(-?\d+)', out2.stdout)
res = {'prefix': prefix, 'p': p, 'constraint_rank': rank, 'charts': charts, 'charts_all_empty': ok_charts, 'square_steps': steps, 'remaining_variables': surv,
       'final': m.groups() if m else out2.stdout[-500:], 'projective_zero_set_empty': bool(ok_charts and m and int(m.group(1)) == 0), 'seconds': round(time.time() - t0, 1)}
json.dump(res, open(prefix + '-charts.json', 'w'), indent=1)
print(json.dumps({k: v for k, v in res.items() if k != 'square_steps'}), flush=True)
