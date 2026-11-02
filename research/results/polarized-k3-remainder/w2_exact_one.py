"""Exact (mod p) reduced-E2 decision for one Jacobian-type candidate from w2-candidates.json (by label substring).
usage: python3 w2_exact_one.py "label substring" [p]"""
import os as _os
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    _os.environ.setdefault(_v, '1')
import sys, json, time
sys.path.insert(0, '/tmp/hodge/results/polarized-k3-remainder')
import numpy as np, w2lib as L
from hodge import core
from hodge.core import mul, add, sub, shift
from flint import nmod_mat
key = sys.argv[1]; p = int(sys.argv[2]) if len(sys.argv) > 2 else L.P0
t0 = time.time()
F, A, q = L.load()
C = json.load(open('w2-candidates.json'))
cands = [c for c in C if key in c[1]]
assert len(cands) == 1, [c[1] for c in cands]
fam, label, expr, frame = cands[0]
D = core.parse(expr); D = L.uv_to_xy(D) if frame == 'uv' else D
d = core.degree(D); pD = core.parity_of(D)
basis = L.jac_basis(d, pD); S_D = L.space(d, pD)
Mb = np.zeros((S_D.n, len(basis)), dtype=np.int64)
for k, (i0, m, Dk) in enumerate(basis):
    for e, c in core.reduce_mod(Dk, p).items(): Mb[S_D.idx[e], k] = c
Dv = S_D.vec(D, dtype=np.int64, p=p)
aug = nmod_mat(S_D.n, len(basis) + 1, np.concatenate([Mb, Dv.reshape(-1, 1)], axis=1).ravel().tolist(), p)
Mr, rk = aug.rref()
a = [0] * len(basis); ok = True
for i in range(rk):
    j = 0
    while j <= len(basis) and int(Mr[i, j]) == 0: j += 1
    if j == len(basis): ok = False; break
    a[j] = int(Mr[i, len(basis)])
rec = {'experiment': 'exact-reduced-e2', 'family': fam, 'label': label, 'd': d, 'parity': pD, 'p': p, 'jacobian_type': ok}
if not ok:
    rec['note'] = 'D not in J(F)_d mod p (E1 fails)'
else:
    avec = [{} for _ in range(6)]
    for k, (i0, m, Dk) in enumerate(basis):
        if a[k]: avec[i0][m] = a[k]
    Dp = core.reduce_mod(D, p); Ap = core.reduce_mod(A, p)
    N0 = [mul(Ap, avec[i], p) if avec[i] else {} for i in range(6)]
    S = L.ReducedE2(d, pD, [(Dp, N0)], verbose=False, build_float=False)
    out = S.exact_consistent(Dp, N0, p=p)
    rec.update({k: out[k] for k in ('consistent', 'rank', 'nullity', 'rows', 'cols', 'seconds')})
rec['total_seconds'] = round(time.time() - t0, 1)
L.log_line('w2-exact-local.jsonl', rec); print(json.dumps(rec))
