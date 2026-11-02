#!/usr/bin/env python3
"""Local (light): classify the stored low-residual minima of a job: distinct minima, their D (top terms), Jacobian/non-J split,
sigma2 parity, sparsity in the exact colon basis.   python3 w2_analyze_minima.py <job> [d] [a,b] [group]"""
import os, sys, json, glob
os.environ['OPENBLAS_NUM_THREADS'] = '1'
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import w2lib
from hodge import core
from fractions import Fraction
job = sys.argv[1]
MIR = '/tmp/hodge/remote/off-paired/wave2-search'
recs = []
for fn in glob.glob(f'{MIR}/wave2-{job}.w*.jsonl'):
    for line in open(fn):
        try: recs.append(json.loads(line))
        except Exception: pass
summ = json.load(open(f'{MIR}/wave2-{job}.summary.json')) if os.path.exists(f'{MIR}/wave2-{job}.summary.json') else {}
d = int(sys.argv[2]) if len(sys.argv) > 2 else summ['d']; c = tuple(int(x) for x in sys.argv[3].split(',')) if len(sys.argv) > 3 else tuple(summ['sector'])
grp = sys.argv[4] if len(sys.argv) > 4 else summ.get('group', 'G')
w2lib.set_group(grp)
F, A, q = w2lib.load()
basis, gauge, info = w2lib.build_family(F, A, q, d, c, summ.get('family', 'colon'))
fam = w2lib.SectorVarPro(F, A, q, basis, c, gauge=gauge)
jac = w2lib.exact_span_basis(w2lib.jacobian_generators_exact(F, d, c), w2lib.monos_char(d, c), tag='jac')
Jm = np.zeros((len(jac), fam.nD))
for k, P in enumerate(jac):
    for e, v in P.items(): Jm[k, fam.Dm.index(e)] = float(v)
PJ = np.linalg.pinv(Jm) @ Jm if len(jac) else np.zeros((fam.nD, fam.nD))
rels = sorted(r['relres'] for r in recs if 'relres' in r)
print(f'job {job}: {len(recs)} restarts, best {rels[0]:.3e}, median {rels[len(rels)//2]:.3e}')
groups = {}
for r in recs:
    if 'Dc' not in r: continue
    key = f"{r['relres']:.4e}"
    groups.setdefault(key, []).append(r)
for key in sorted(groups, key=float)[:8]:
    rs = groups[key]; r = rs[0]; Dc = np.array(r['Dc']); Dc = Dc / np.linalg.norm(Dc)
    projJ = PJ @ Dc; fracJ = float(np.linalg.norm(projJ))
    # sigma2 parity: swap X2<->Y2 on coefficients
    Ds = np.zeros_like(Dc)
    for i, m in enumerate(fam.Dm):
        m2 = (m[0], m[1], m[5], m[3], m[4], m[2])
        if m2 in fam.Dm: Ds[fam.Dm.index(m2)] = Dc[i]
    s2 = float(Dc @ Ds)
    top = sorted([(abs(v), core.to_str({fam.Dm[i]: Fraction(round(float(v), 3))})) for i, v in enumerate(Dc)], reverse=True)[:6]
    print(f'  minimum relres {key}: hit by {len(rs)} restarts; |J-part| = {fracJ:.3f}; <D, sigma2 D> = {s2:+.3f}; top terms: ' + ', '.join(t[1] for t in top))
