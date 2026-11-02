"""Validation of w2fast: (a) non-strict accepted q=1 witnesses (Q != 0) must give f ~ 0 at the witness D;
(b) Koszul-perturbed lifts N' + K(p0) must be absorbed (tests T against Vb, independently coded)."""
import sys, json, time, os
sys.path.insert(0, '/tmp/hodge/results/klein-coupled-character')
import numpy as np
import w2fast as wf
from hodge import core
from hodge.core import NV, pderiv, shift, add, sub, scale, monomials
from fractions import Fraction
ROOT = '/home/kejace/devel/sextics'
cat = json.load(open(f'{ROOT}/data/catalogue.json'))
found = []
for iid in cat['instanceIds']:
    if not os.path.exists(f'{ROOT}/data/{iid}/witness.json'):
        continue
    F, A, q, rec = core.load_instance(iid)
    if q != 1:
        continue
    d, D, N, Q = core.load_witness(iid)
    if d > 7:
        continue
    found.append((iid, d, bool(Q), core.parity_of(D)))
print('q=1 witnesses with d<=7:', found, flush=True)
rng = np.random.default_rng(5)
results = {}
for iid, d, hasQ, parD in found[:12]:
    F, A, q, rec = core.load_instance(iid)
    d, D, N, Q = core.load_witness(iid)
    t0 = time.time()
    Np = wf.normalize_e1(F, A, D, N)
    sec = wf.Sector(F, A, d, parD, None)
    # (b) perturb by a random exact Koszul term p0 (small integers) in the sector columns
    Fi = [pderiv(F, i) for i in range(NV)]
    Npert = [dict(n) for n in Np]
    if sec.nP:
        for (i, j, m) in sec.Pcols[: min(sec.nP, 40)]:
            c = int(rng.integers(-3, 4))
            if c == 0:
                continue
            Npert[i] = add(Npert[i], shift(Fi[j], m, c))
            Npert[j] = sub(Npert[j], shift(Fi[i], m, c))
    Bc, Bn = sec.family_from_polys([(D, Np), (D, Npert)])
    # single-element families
    k1 = wf.KVP(sec, Bc[:, :1], Bn[:, :1], exclude_F_multiples=False)
    f1, b1, xP1, r1, rk1 = k1.value(np.array([1.0]))
    k2 = wf.KVP(sec, Bc[:, 1:], Bn[:, 1:], exclude_F_multiples=False)
    f2, b2, xP2, r2, rk2 = k2.value(np.array([1.0]))
    results[iid] = {'d': d, 'Q_nonzero': hasQ, 'rows': sec.nrows, 'nP': sec.nP, 'nQ': sec.nQ,
                    'f_witness': f1, 'b2_witness': b1, 'f_perturbed': f2, 'b2_perturbed': b2, 'rankT': rk2,
                    'seconds': round(time.time() - t0, 1)}
    print(iid, results[iid], flush=True)
json.dump(results, open('w2_validate.out.json', 'w'), indent=1)
